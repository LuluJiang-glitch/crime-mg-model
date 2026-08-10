# MG model for City of London crime data, 2023-2025.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) > 0) {
  base_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  base_dir <- normalizePath(getwd())
}

# base_dir <- "/Users/JIANL0A/Desktop/Crime/L"
setwd(base_dir)
prepare_dir <- file.path(base_dir, "RData", "Prepare")
spatial_dir <- file.path(base_dir, "RData", "Spatial")
mg_dir <- file.path(base_dir, "RData", "MG")
point_dir <- file.path(base_dir, "RData", "Point")
for (d in c(prepare_dir, spatial_dir, mg_dir, point_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
packages <- c("tidyverse", "rnaturalearth", "rnaturalearthdata", "ggplot2", "sf", "osmdata", "dplyr", "readr", "MetricGraph", "fmesher", "INLA", "tibble", "tidyr", "inlabru", "units", "deldir", "purrr", "sp", "fields", "FNN")
invisible(lapply(packages, library, character.only = TRUE))
source(file.path(base_dir, "metric_graph.R"))
source(file.path(base_dir, "Function_Model.R"))

Region <- "City"
crime_type <- "Theft from the person"
# "Theft from the person", "Robbery", "Drugs", "Bicycle theft"
# "Anti-social behaviour", "Criminal damage and arson", "Violence and sexual offences", "Vehicle crime"
# "Burglary", "Shoplifting"
model_years <- 2023:2025
large_graph_file <- file.path(prepare_dir, "02_graph_City_L100.RData")
step_size <- 0.05
graph_crs <- 4326
metric_crs <- 27700

# create boundary
uk <- ne_states(country = "united kingdom", returnclass = "sf")
london <- uk[uk$region == "Greater London", ]
boundary <- st_transform(
  london[london$name %in% c(Region), ],
  crs = graph_crs
) %>%
  st_union()
boundary_sf <- st_as_sf(boundary) %>%
  st_make_valid()


load(file.path(prepare_dir, "01_data_city.RData"))
data_all_raw <- data_city %>%
  filter(`Crime type` == crime_type) %>%
  group_by(Longitude, Latitude, Year) %>%
  summarise(crime = sum(count), .groups = "drop") %>%
  st_drop_geometry()


# Assign ID to every original coordinate.
# Longitude Latitude original_loc_id
loc_data_all_original <- data_all_raw %>%
  distinct(Longitude, Latitude) %>%
  arrange(Longitude, Latitude) %>%
  mutate(original_loc_id = row_number())

# merge close points
merge_res <- merge_points_within_distance(
  loc_df = loc_data_all_original %>%
    select(Longitude, Latitude),
  dist_m = 30,
  crs_longlat = graph_crs,
  crs_meter = metric_crs
)

message(
  "Coordinates merged: original points = ",
  merge_res$n_original,
  "; merged points = ",
  merge_res$n_merged,
  "; threshold = ",
  30,
  " meters."
)


# Map original coordinates to their merged locations.
# original_loc_id Longitude_original Latitude_original loc_id Longitude Latitude
loc_cluster_map <- loc_data_all_original %>%
  select(
    original_loc_id,
    Longitude_original = Longitude,
    Latitude_original = Latitude
  ) %>%
  left_join(
    merge_res$loc_cluster_map %>%
      select(
        original_loc_id,
        loc_id,
        Longitude,
        Latitude
      ),
    by = "original_loc_id"
  )

# Longitude Latitude loc_id n_original_points
loc_data_all_df <- merge_res$loc_merged_df %>%
  select(
    Longitude,
    Latitude,
    loc_id,
    n_original_points
  ) %>%
  arrange(loc_id)

# Coordinates used to construct the Voronoi partition.
# Longitude Latitude (merged)
loc_data_all <- loc_data_all_df %>% select(Longitude, Latitude)


# Re-aggregate crime counts after merging nearby released locations.
# loc_id Longitude Latitude  Year crime
data_all <- data_all_raw %>%
  left_join(
    loc_data_all_original,
    by = c("Longitude", "Latitude")
  ) %>%
  left_join(
    loc_cluster_map %>%
      select(
        original_loc_id,
        loc_id,
        Longitude_merged = Longitude,
        Latitude_merged = Latitude
      ),
    by = "original_loc_id"
  ) %>%
  group_by(
    loc_id,
    Longitude = Longitude_merged,
    Latitude = Latitude_merged,
    Year
  ) %>%
  summarise(
    crime = sum(crime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(loc_id, Year)

data <- data_all %>% filter(Year %in% model_years)
loc_data <- data %>% select(Longitude, Latitude)
loc_data_df <- as.data.frame(loc_data)
colnames(loc_data_df) <- c("Longitude", "Latitude")




# load graph
load(large_graph_file)
graph$clear_observations()
bb <- st_bbox(boundary_sf)
rw <- c(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])

v <- deldir(
  loc_data_all$Longitude,
  loc_data_all$Latitude,
  rw = rw
)

dirsg <- lapply(seq_len(nrow(v$dirsgs)), function(i) {
  sf::st_linestring(
    rbind(
      c(v$dirsgs$x1[i], v$dirsgs$y1[i]),
      c(v$dirsgs$x2[i], v$dirsgs$y2[i])
    )
  )
})

edge_Voronoi <- do.call(sf::st_sfc, dirsg)
sf::st_crs(edge_Voronoi) <- graph_crs
edge_Voronoi <- sf::st_sf(geometry = edge_Voronoi)
edge_MG <- graph$get_edges(format = "sf")
s2_old <- sf::sf_use_s2()
sf::sf_use_s2(FALSE)

boundary_pts <- suppressWarnings(
  st_intersection(
    st_geometry(edge_MG),
    st_boundary(st_geometry(boundary_sf))
  )
)

boundary_pts <- st_sf(geometry = boundary_pts)
boundary_pts <- suppressWarnings(
  st_collection_extract(boundary_pts, "POINT")
)
boundary_pts <- suppressWarnings(
  st_cast(boundary_pts, "POINT")
)

intsec_pts <- suppressWarnings(
  st_intersection(
    st_geometry(edge_Voronoi),
    st_geometry(edge_MG)
  )
)

intsec_pts <- st_sf(geometry = intsec_pts)
intsec_pts <- suppressWarnings(
  st_collection_extract(intsec_pts, "POINT")
)
intsec_pts <- suppressWarnings(
  st_cast(intsec_pts, "POINT")
)

# Keep only Voronoi intersections inside the original City boundary.
keep_intsec <- lengths(
  st_intersects(intsec_pts, boundary_sf)
) > 0L

intsec_pts <- intsec_pts[
  keep_intsec,
  ,
  drop = FALSE
]



# Combine intsec points
intsec_pts <- dplyr::bind_rows(
  boundary_pts,
  intsec_pts
)

coords <- st_coordinates(intsec_pts)[, c("X", "Y"), drop = FALSE]
point_key <- paste(
  round(coords[, 1], 8),
  round(coords[, 2], 8),
  sep = ":"
)
intsec_pts <- intsec_pts[!duplicated(point_key), , drop = FALSE]
coords <- st_coordinates(intsec_pts)
intsec_pts$x <- coords[, "X"]
intsec_pts$y <- coords[, "Y"]

graph$add_observations(
  data = intsec_pts,
  coord_x = "x",
  coord_y = "y",
  data_coords = "spatial",
  verbose = 0
)

graph$observation_to_vertex(mesh_warning = FALSE)
graph$clear_observations()
edge_graph_all <- graph$get_edges(format = "sf") %>% mutate(edge_id = row_number())

inside_edge <- lengths(
  st_covered_by(
    st_geometry(edge_graph_all),
    st_geometry(boundary_sf)
  )
) > 0L

edge_graph_01 <- edge_graph_all[
  inside_edge,
  ,
  drop = FALSE
]

sf::sf_use_s2(s2_old)

message(
  "Large graph edges used by SPDE: ", nrow(edge_graph_all),
  "; edges used for integration inside original boundary: ",
  nrow(edge_graph_01)
)




# Assign every graph edge to the nearest merged crime location.
mid_points <- t(
  vapply(
    edge_graph_01$geometry,
    FUN = function(coords) colMeans(coords),
    FUN.VALUE = numeric(2L)
  )
)

mid_points_df <- data.frame(
  Longitude = mid_points[, 1],
  Latitude = mid_points[, 2]
)

nn_edges <- FNN::get.knnx(
  data = loc_data_all_df[, c("Longitude", "Latitude")],
  query = as.matrix(mid_points_df[, c("Longitude", "Latitude")]),
  k = 1
)

edge_graph_01$loc_id <- loc_data_all_df$loc_id[nn_edges$nn.index[, 1]]
edge_graph_01$Longitude <- loc_data_all_df$Longitude[nn_edges$nn.index[, 1]]
edge_graph_01$Latitude <- loc_data_all_df$Latitude[nn_edges$nn.index[, 1]]

used_loc_ids <- sort(unique(edge_graph_01$loc_id))
all_loc_ids <- loc_data_all_df$loc_id
unused_loc_ids <- setdiff(all_loc_ids, used_loc_ids)
if (length(unused_loc_ids) == 0) {
  message("All merged coordinate points are already used by graph edges; no reassignment is needed.")
  transfer_map <- integer(0)
} else {
  coords_used <- loc_data_all_df %>%
    filter(loc_id %in% used_loc_ids) %>%
    select(Longitude, Latitude)
  
  coords_unused <- loc_data_all_df %>%
    filter(loc_id %in% unused_loc_ids) %>%
    select(Longitude, Latitude)
  
  nn_transfer <- FNN::get.knnx(
    data = coords_used,
    query = coords_unused,
    k = 1
  )
  
  recipient_loc_ids <- used_loc_ids[nn_transfer$nn.index[, 1]]
  
  transfer_map <- recipient_loc_ids
  names(transfer_map) <- unused_loc_ids
  
  message("Created a mapping from unused merged coordinates to the nearest used merged coordinates.")
}

loc_data_all_map <- loc_data_all_df %>% distinct(Longitude, Latitude, loc_id)



# Add counts from unused locations to nearest location
# loc_id_final  Year crime
data1 <- data %>%
  mutate(
    Year = as.integer(Year),
    loc_id_final = ifelse(
      loc_id %in% unused_loc_ids,
      as.integer(transfer_map[as.character(loc_id)]),
      loc_id
    )
  ) %>%
  group_by(loc_id_final, Year) %>%
  summarise(
    crime = sum(crime, na.rm = TRUE),
    .groups = "drop"
  )

edge_graph_01_aug <- edge_graph_01 %>% mutate(data_after = loc_id)

block_map <- match(
  edge_graph_01_aug$data_after,
  sort(unique(edge_graph_01_aug$data_after))
)

idx_all <- integer(0)
u_all <- numeric(0)
w_all <- numeric(0)
b_all <- integer(0)

for (i in seq_len(nrow(edge_graph_01_aug))) {
  
  rr <- simpson_u_weight(
    edge_graph_01_aug$.edge_lengths[i],
    step = step_size
  )
  
  m <- length(rr$u)
  
  idx_all <- c(
    idx_all,
    rep(edge_graph_01_aug$edge_id[i], m)
  )
  
  u_all <- c(u_all, rr$u)
  w_all <- c(w_all, rr$w)
  b_all <- c(b_all, rep(block_map[i], m))
}

ips <- tibble(
  x = as_MGG(
    cbind(idx_all, u_all),
    graph = graph
  ),
  weight = w_all,
  .block = b_all
)


# plot ips
# ips_plot <- data.frame(
#   edge_number = ips$x$index,
#   distance_on_edge = ips$x$where[, 2],
#   ips_id = seq_len(nrow(ips))
# )
# 
# graph$clear_observations()
# 
# graph$add_observations(
#   data = ips_plot,
#   normalized = TRUE
# )
# 
# graph$plot(
#   data = "ips_id",
#   type = "ggplot",
#   vertex_size = 0,
#   data_size = 0.4
# ) +
#   geom_sf(
#     data = boundary_sf,
#     fill = NA,
#     colour = "red",
#     linewidth = 0.7
#   ) +
#   geom_sf(
#     data = voronoi_in_boundary,
#     colour = "grey50",
#     linewidth = 0.3,
#     linetype = "22"
#   ) +
#   geom_sf(
#     data = loc_points_sf,
#     colour = "blue",
#     size = 0.8
#   ) +
#   ggtitle("Distribution of Simpson integration points") +
#   theme_bw()


loc_points_sf <- st_as_sf(
  loc_data_all_df,
  coords = c("Longitude", "Latitude"),
  crs = graph_crs
)

voronoi_in_boundary_raw <- suppressWarnings(
  st_intersection(
    st_geometry(edge_Voronoi),
    st_geometry(boundary_sf)
  )
)

voronoi_in_boundary <- st_sf(geometry = voronoi_in_boundary_raw)
voronoi_geom_type <- as.character(st_geometry_type(voronoi_in_boundary))
voronoi_in_boundary <- voronoi_in_boundary[
  voronoi_geom_type %in% c("LINESTRING", "MULTILINESTRING"),
  ,
  drop = FALSE
]


# ggplot() +
  # geom_sf(data = edge_graph_all, colour = "grey80", linewidth = 0.25) +
  # geom_sf(data = voronoi_in_boundary, colour = "grey55", linewidth = 0.25, linetype = "22") +
  # geom_sf(data = edge_graph_01, colour = "#0072B2", linewidth = 0.5) +
  # geom_sf(data = boundary_sf, fill = NA, colour = "red", linewidth = 0.8) +
  # geom_sf(data = loc_points_sf, colour = "black", size = 0.7, alpha = 0.65) +
  # coord_sf(expand = FALSE) +
  # theme_bw() +
  # labs(
  #   title = paste0("MG integration domain for ", Region),
  #   subtitle = paste0(
  #     "Grey = enlarged graph used by SPDE; blue = integration-domain edges; ",
  #     "red = original boundary"
  #   ),
  #   x = NULL,
  #   y = NULL,
  #   caption = "Black points are merged released crime locations; dashed grey lines are Voronoi boundaries inside the original boundary."
  # )


# plot ips corresponding to y block
# plot_region_df <- edge_graph_01 %>%
#   rename(loc_id_final = loc_id) %>%
#   left_join(
#     y_df %>% select(Year, loc_id_final, crime),
#     by = "loc_id_final"
#   ) %>%
#   mutate(
#     log_crime = log1p(crime)
#   )
# 
# ggplot() +
#   geom_sf(
#     data = plot_region_df,
#     aes(colour = log_crime),
#     linewidth = 0.8
#   ) +
#   geom_sf(
#     data = boundary_sf,
#     fill = NA,
#     colour = "black",
#     linewidth = 0.5
#   ) +
#   facet_wrap(~ Year) +
#   scale_colour_viridis_c(name = "log1p(crime)") +
#   theme_bw()


block_ids <- sort(unique(ips$.block))
n_block <- length(block_ids)

block_map <- tibble(
  loc_id_final = block_ids,
  block_base = seq_len(n_block)
)

ips <- ips %>%
  mutate(loc_id_final = .block) %>%
  left_join(block_map, by = "loc_id_final") %>%
  mutate(.block = block_base) %>%
  select(-block_base)

pte_ips <- as.matrix(
  cbind(
    ips$x$index,
    ips$x$where[, 2]
  )
)

pte_key <- paste(
  pte_ips[, 1],
  formatC(pte_ips[, 2], digits = 17, format = "fg", flag = "#"),
  sep = ":"
)
unique_pte <- !duplicated(pte_key)
ips$spde_loc_id <- match(pte_key, unique(pte_key))

spde_loc_data <- data.frame(
  spde_loc_id = seq_len(sum(unique_pte)),
  edge_number = as.integer(pte_ips[unique_pte, 1]),
  distance_on_edge = pte_ips[unique_pte, 2]
)




# covariates.
key_list <- c("amenity", "highway", "man_made", "railway", "shop")

value_list <- list(
  amenity = c("bar", "nightclub", "bus_station"), # bank, "police"
  railway = c("subway_entrance"),
  shop = c("supermarket", "convenience")  #"department_store"
)

for (key in key_list) {
  for (value in c(value_list[[key]])) {
    cov_name <- paste(key, value, sep = "_")
    
    if (exists(cov_name, envir = .GlobalEnv)) {
      rm(list = cov_name, envir = .GlobalEnv)
    }
    
    message("Processing: ", Region, ": ", cov_name)
    
    path_cov1 <- file.path(prepare_dir, "cov/loc")
    path_cov2 <- file.path(prepare_dir, "cov/dist/MG", crime_type)
    
    dir.create(path_cov1, recursive = TRUE, showWarnings = FALSE)
    dir.create(path_cov2, recursive = TRUE, showWarnings = FALSE)
    
    path_dist <- paste0(path_cov2, "/", key, "_", value, "_", step_size, "m.RData")
    
    if (file.exists(path_dist)) {
      load(path_dist)
      assign(cov_name, cov_dist_geo)
      
    } else {
      
      message("No existing distance file for: ", key, " : ", value)
      
      path_loc <- paste0(
        path_cov1,
        "/",
        key,
        "_",
        value,
        "_loc.RData"
      )
      
      if (file.exists(path_loc)) {
        
        load(path_loc)
        
      } else {
        
        loc_cov <- extract_covariates(
          boundary = boundary_sf,
          key = key,
          value = value
        )
        
        if (is.null(loc_cov) || nrow(as.data.frame(loc_cov)) == 0) {
          message("No data found for ", key, " = ", value, " ... skipping")
          next
        }
        
        data_cov <- check_unique_loc2pte(graph, loc_cov)
        
        data_cov <- cbind(
          data_cov,
          type = rep(
            paste(key, value, sep = ":"),
            nrow(data_cov)
          )
        )
        
        loc_cov <- cbind(
          data_cov$Longitude,
          data_cov$Latitude
        )
        
        save(loc_cov, file = path_loc)
      }
      
      pte_ips_unique_list <- check_pte_loc_unique(
        pte = pte_ips,
        graph = graph
      )
      
      aa <- compute_geo_matdist(
        type = "geo",
        loc_mesh = pte_ips_unique_list$loc_unique,
        pte_mesh = pte_ips_unique_list$pte_unique,
        loc_data = loc_cov,
        graph = graph
      )
      
      if (is.null(dim(aa))) {
        cov_dist_geo <- aa
        cov_dist_geo1 <- aa
      } else {
        cov_dist_geo1 <- apply(aa, 2, min)
        cov_dist_geo <- cov_dist_geo1[pte_ips_unique_list$index]
      }
      
      save(cov_dist_geo, file = path_dist)
      assign(cov_name, cov_dist_geo)
    }
  }
}



# Load the graph distances and transform d to exp(-d).
covariate_store <- list()

for (key in key_list) {
  for (value in c(value_list[[key]])) {
    cov_name <- paste(key, value, sep = "_")
    
    if (exists(cov_name, envir = .GlobalEnv)) {
      rm(list = cov_name, envir = .GlobalEnv)
    }
    
    message("Loading covariate: ", Region, ": ", cov_name)
    
    path_cov2 <- file.path(prepare_dir, "cov/dist/MG", crime_type)
    
    path_dist <- paste0(path_cov2, "/", key, "_", value, "_", step_size, "m.RData")
    
    if (file.exists(path_dist)) {
      
      load(path_dist)
      
      cov_dist_geo <- exp(-cov_dist_geo)
      assign(cov_name, cov_dist_geo)
      
      covariate_store[[cov_name]] <- cov_dist_geo
      
    } else {
      
      message("No existing covariate file for: ", key, " : ", value)
      next
    }
  }
}

covariate_names <- names(covariate_store)
covariate_df <- as_tibble(covariate_store)
ips_cov <- bind_cols(ips, covariate_df)


# ips_plot <- data.frame(
#   edge_number = ips$x$index,
#   distance_on_edge = ips$x$where[, 2],
#   ips_cov[, covariate_names]
# )
# 
# graph$clear_observations()
# 
# graph$add_observations(
#   data = ips_plot,
#   normalized = TRUE
# )
# 
# for (cov in covariate_names) {
# 
#   path_loc <- file.path(
#     prepare_dir,
#     "cov/loc",
#     paste0(cov, "_loc.RData")
#   )
# 
#   load(path_loc)   # load loc_cov
# 
#   loc_cov_df <- data.frame(
#     Longitude = loc_cov[, 1],
#     Latitude = loc_cov[, 2]
#   )
# 
#   p <- graph$plot(
#     data = cov,
#     type = "ggplot",
#     vertex_size = 0,
#     data_size = 1
#   ) +
#     geom_sf(
#       data = boundary_sf,
#       fill = NA,
#       colour = "red",
#       linewidth = 0.6
#     ) +
#     geom_point(
#       data = loc_cov_df,
#       aes(x = Longitude, y = Latitude),
#       colour = "blue",
#       size = 1
#     ) +
#     ggtitle(cov) +
#     theme_bw()
# 
#   print(p)
# }




# Each year uses the same ips but a distinct block index.
ips_final <- purrr::map_dfr(
  model_years,
  function(yy) {
    r <- match(yy, model_years)
    ips_cov %>%
      mutate(
        Year = yy,
        rep = r,
        .block = .block + n_block * (r - 1L)
      )
  }
)


# Create a year-block response table; unobserved counts are zero.
y_df <- tidyr::expand_grid(
  Year = model_years,
  loc_id_final = block_ids
) %>%
  left_join(data1, by = c("loc_id_final", "Year")) %>%
  left_join(block_map, by = "loc_id_final") %>%
  mutate(
    crime = tidyr::replace_na(crime, 0L),
    rep = match(Year, model_years),
    .block = block_base + n_block * (rep - 1L)
  ) %>%
  arrange(.block)

y <- y_df$crime




# Log-Gaussian priors for the marginal standard deviation and range.
p_sigma <- 0.01
p_range <- 0.01

sigma_median <- 1
range_median <- 0.05

mu_sigma <- log(sigma_median)
sd_sigma <- (log(5) - mu_sigma) / qnorm(1 - p_sigma)

mu_range <- log(range_median)
sd_range <- (mu_range - log(0.01)) / qnorm(1 - p_range)

graph$clear_observations()
graph$add_observations(
  data = spde_loc_data,
  normalized = TRUE
)

# nu = 0.5, alpha = 1.
spde_model <- graph_spde(
  graph_object = graph,
  alpha = 1,
  stationary_endpoints = "all",
  parameterization = "matern",
  start_sigma = sigma_median,
  start_range = range_median,
  prior_sigma = list(
    meanlog = mu_sigma,
    sdlog = sd_sigma
  ),
  prior_range = list(
    meanlog = mu_range,
    sdlog = sd_range
  )
)

spde_model_data <- graph_data_spde(
  graph_spde = spde_model,
  loc_name = "loc",
  drop_na = FALSE,
  drop_all_na = FALSE
)[["data"]]

spde_loc_ids <- spde_model_data$spde_loc_id
spde_loc_matrix <- spde_model_data$loc
spde_match <- match(ips_final$spde_loc_id, spde_loc_ids)
ips_final_data <- as.list(ips_final)
ips_final_data$loc <- spde_loc_matrix[spde_match, , drop = FALSE]

agg <- bru_mapper_aggregate(
  rescale = FALSE,
  n_block = n_block * length(model_years),
  type = "logsumexp"
)

all_vars <- c(covariate_names)

all_vars <- all_vars[
  all_vars %in% names(ips_final)
]

state_str <- paste(
  c(all_vars, "spde"),
  collapse = " + "
)

formula_str <- paste0(
  "y ~ Intercept + ibm_eval(",
  "agg, ",
  "input = list(block = .block, weights = weight), ",
  "state = ",
  state_str,
  ")"
)

formula <- as.formula(formula_str)

lik <- bru_obs(
  formula = formula,
  response_data = data.frame(y = y),
  family = "poisson",
  data = ips_final_data,
  allow_combine = TRUE
)

covariate_terms <- all_vars

if (length(covariate_terms) > 0) {
  cmp_str <- paste(
    "y ~ Intercept(1) +",
    paste(covariate_terms, collapse = " + "),
    "+ spde(loc, model = spde_model, replicate = rep)"
  )
} else {
  cmp_str <- paste(
    "y ~ Intercept(1) +",
    "spde(loc, model = spde_model, replicate = rep)"
  )
}

cmp <- as.formula(cmp_str)

fit <- bru(
  components = cmp,
  lik,
  options = list(
    control.fixed = list(
      mean = 0,
      prec = 1 / 10,
      mean.intercept = 0,
      prec.intercept = 1 / 10
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE,
      config = TRUE
    ),
    control.inla = list(
      int.strategy = "eb"
    )
  )
)

summary(fit)

spde_result <- spde_metric_graph_result(
  fit,
  "spde",
  spde_model
)

summary(spde_result)

save(fit, file = file.path(mg_dir, paste0("fit_", crime_type, ".RData")))
message("MG model saved in: ", mg_dir)
