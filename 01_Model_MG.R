# MG model for City of London crime data, 2023-2025.

code_dir <- "~/Desktop/Crime/code_final"
base_dir <- "~/Desktop/Crime/code_final"
prepare_dir <- file.path(base_dir, "RData", "Prepare")
spatial_dir <- file.path(base_dir, "RData", "Spatial")
mg_dir <- file.path(base_dir, "RData", "MG")
point_dir <- file.path(base_dir, "RData", "Point")
for (d in c(prepare_dir, spatial_dir, mg_dir, point_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

packages <- c("tidyverse", "rnaturalearth", "rnaturalearthdata", "ggplot2", "sf", "osmdata", "dplyr", "readr", "rSPDE", "MetricGraph", "fmesher", "INLA", "tibble", "tidyr", "inlabru", "units", "deldir", "purrr", "sp", "fields", "FNN")
invisible(lapply(packages, library, character.only = TRUE))

source("~/Desktop/Crime/code_final/metric_graph.R")
source(file.path(code_dir, "Function_Model.R"))

Region <- "City"
crime_type <- "Bicycle theft"
# "Theft from the person", "Robbery", "Drugs", "Bicycle theft"
# "Anti-social behaviour", "Criminal damage and arson", "Violence and sexual offences", "Vehicle crime"
# "Burglary", "Shoplifting"


# Load data.
load(file.path(prepare_dir, "01_crime_yearly_data_all.RData"))
uk <- ne_states(country = "united kingdom", returnclass = "sf")
london <- uk[uk$region == "Greater London", ]
boundary <- st_transform(
  london[london$name %in% c(Region), ],
  crs = 4326
) %>%
  st_union()
boundary_sf <- st_as_sf(boundary)
crs <- sf::st_crs(boundary_sf)

data_all_raw <- yearly_data %>%
  filter(`Crime type` == crime_type) %>%
  filter(Year %in% c("2023", "2024", "2025")) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  filter(rowSums(st_within(geometry, boundary, sparse = FALSE)) > 0) %>%
  mutate(
    Longitude = st_coordinates(.)[, 1],
    Latitude = st_coordinates(.)[, 2],
    Year = as.numeric(Year)
  ) %>%
  group_by(Longitude, Latitude, Year) %>%
  summarise(crime = sum(count), .groups = "drop") %>%
  st_drop_geometry()

loc_data_all_original <- data_all_raw %>%
  distinct(Longitude, Latitude) %>%
  arrange(Longitude, Latitude) %>%
  mutate(original_loc_id = row_number())

merge_res <- merge_points_within_distance(
  loc_df = loc_data_all_original %>%
    select(Longitude, Latitude),
  dist_m = 30,
  crs_longlat = 4326,
  crs_meter = 27700
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

loc_data_all_df <- merge_res$loc_merged_df %>%
  select(
    Longitude,
    Latitude,
    loc_id,
    n_original_points
  ) %>%
  arrange(loc_id)
loc_data_all <- loc_data_all_df %>% select(Longitude, Latitude)

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

data <- data_all %>%
  filter(Year %in% 2023:2025)

loc_data <- data %>%
  select(Longitude, Latitude)

loc_data_df <- as.data.frame(loc_data)
colnames(loc_data_df) <- c("Longitude", "Latitude")


load("~/Desktop/Crime/code_final/RData/Prepare/02_graph_City.RData")

bb <- st_bbox(boundary)
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
sf::st_crs(edge_Voronoi) <- crs
edge_Voronoi <- sf::st_sf(geometry = edge_Voronoi)

edge_MG <- sf::st_transform(
  graph$get_edges(format = "sf"),
  crs
)


intsec_file <- file.path(prepare_dir, paste0("03_intsec_", Region, "_", crime_type, ".RData"))

if (file.exists(intsec_file)) {
  load(intsec_file)
  message("Loaded existing intsec_pts: ", intsec_file)
} else {
  intsec_pts <- sf::st_intersection(edge_Voronoi, edge_MG)
  intsec_pts <- intsec_pts[sf::st_geometry_type(intsec_pts) == "POINT", ]
  
  coords <- sf::st_coordinates(intsec_pts)[, c("X", "Y")]
  intsec_pts <- data.frame(
    x = coords[, 1],
    y = coords[, 2]
  )
  
  intsec_pts <- sf::st_as_sf(
    intsec_pts,
    coords = c("x", "y"),
    crs = crs
  )
  
  save(intsec_pts, file = intsec_file)
  message("Saved new intsec_pts: ", intsec_file)
}

graph$add_observations(
  data = intsec_pts,
  coord_x = "x",
  coord_y = "y",
  data_coords = "spatial"
)

graph$observation_to_vertex(mesh_warning = TRUE)


edge_graph_01 <- graph$get_edges(format = "sf")

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

edge_graph_01$loc_id <- nn_edges$nn.index[, 1]
edge_graph_01$Longitude <- loc_data_all_df$Longitude[edge_graph_01$loc_id]
edge_graph_01$Latitude <- loc_data_all_df$Latitude[edge_graph_01$loc_id]

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

loc_data_all_map <- loc_data_all_df %>%
  distinct(Longitude, Latitude, loc_id)


if (!("loc_id" %in% names(data_all))) {
  data_all <- data_all %>%
    left_join(loc_data_all_map, by = c("Longitude", "Latitude"))
}

if (!("loc_id" %in% names(data))) {
  data <- data %>%
    left_join(loc_data_all_map, by = c("Longitude", "Latitude"))
}


data1 <- data %>%
  mutate(
    loc_id_final = ifelse(
      loc_id %in% unused_loc_ids,
      as.integer(transfer_map[as.character(loc_id)]),
      loc_id
    )
  ) %>%
  group_by(loc_id_final) %>%
  summarise(
    crime = sum(crime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    loc_data_all_df %>%
      select(loc_id, Longitude, Latitude) %>%
      rename(loc_id_final = loc_id),
    by = "loc_id_final"
  ) %>%
  arrange(loc_id_final)

edge_graph_01_aug <- edge_graph_01 %>%
  mutate(
    edge_id = row_number(),
    data_after = loc_id
  ) %>%
  left_join(
    data1 %>% select(loc_id_final, crime),
    by = c("data_after" = "loc_id_final")
  ) %>%
  mutate(crime = replace_na(crime, 0L))


step_size <- 0.05

ips <- make_ips(
  edge_graph_01_aug,
  graph,
  step = step_size
)

crime_map <- edge_graph_01_aug %>%
  sf::st_drop_geometry() %>%
  select(data_after, crime) %>%
  group_by(data_after) %>%
  summarise(
    crime = if (n_distinct(crime, na.rm = TRUE) <= 1) {
      dplyr::first(na.omit(crime))
    } else {
      max(crime, na.rm = TRUE)
    },
    .groups = "drop"
  )

ips_one <- ips %>%
  mutate(block_orig = .block) %>%
  left_join(
    crime_map,
    by = c("block_orig" = "data_after")
  ) %>%
  mutate(crime = tidyr::replace_na(crime, 0)) %>%
  select(-block_orig)

block_ids <- sort(unique(ips$.block))
n_block <- length(block_ids)

y_df <- ips_one %>%
  arrange(.block) %>%
  group_by(.block) %>%
  summarise(
    y = first(crime),
    .groups = "drop"
  )

y <- y_df %>%
  arrange(.block) %>%
  pull(y)

pte_ips <- as.matrix(
  cbind(
    ips$x$index,
    ips$x$where[, 2]
  )
)

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
        
        cov_sf <- st_as_sf(
          as.data.frame(loc_cov),
          coords = c(1, 2),
          crs = 4326
        )
        
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

model_years <- 2023:2025

ips_final <- purrr::map_dfr(
  model_years,
  function(yy) {
    ips_cov %>%
      mutate(
        Year = yy,
        rep = as.integer(factor(yy, levels = model_years))
      )
  }
)

p_sigma <- 0.01
p_range <- 0.01

sigma_median <- 1
range_median <- 0.05

mu_sigma <- log(sigma_median)
sd_sigma <- (log(5) - mu_sigma) / qnorm(1 - p_sigma)

mu_range <- log(range_median)
sd_range <- (mu_range - log(0.01)) / qnorm(1 - p_range)

rspde_model <- rspde.metric_graph(
  graph,
  nu = 0.5,
  parameterization = "matern",
  start.lstd.dev = mu_sigma,
  start.lrange = mu_range,
  theta.prior.mean = c(mu_sigma, mu_range),
  theta.prior.prec = diag(c(
    1 / sd_sigma^2,
    1 / sd_range^2
  ))
)

agg <- bru_mapper_aggregate(
  rescale = FALSE,
  n_block = n_block,
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
  data = ips_final,
  allow_combine = TRUE
)

covariate_terms <- all_vars

if (length(covariate_terms) > 0) {
  cmp_str <- paste(
    "y ~ Intercept(1) +",
    paste(covariate_terms, collapse = " + "),
    "+ spde(x, model = rspde_model, mapper = bru_mapper(graph), replicate = rep)"
  )
} else {
  cmp_str <- paste(
    "y ~ Intercept(1) +",
    "spde(x, model = rspde_model, mapper = bru_mapper(graph), replicate = rep)"
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

summary(
  rspde.result(
    fit,
    "spde",
    rspde_model
  )
)

save(fit, file = file.path(mg_dir, paste0("fit_", crime_type, "_M30.RData")))
message("MG model saved in: ", mg_dir)
