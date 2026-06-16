# Spatial model for City of London crime data, 2023-2025.

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
crime_type <- "Theft from the person"
# "Theft from the person", "Robbery", "Drugs", "Bicycle theft"
# "Anti-social behaviour", "Criminal damage and arson", "Violence and sexual offences"
# "Burglary", "Shoplifting", "Vehicle crime"

nsub2 <- 5L

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
    Latitude  = st_coordinates(.)[, 2],
    Year = as.numeric(Year)
  ) %>%
  group_by(Longitude, Latitude, Year) %>%
  summarise(
    crime = sum(count),
    .groups = "drop"
  ) %>%
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

data <- data_all %>% filter(Year %in% 2023:2025)

boundary_m <- st_transform(boundary, 27700)

loc_data_all_sf <- st_as_sf(
  loc_data_all_df,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(27700)

loc_data_all_m <- st_coordinates(loc_data_all_sf)

loc_data_all_map <- loc_data_all_df %>%
  distinct(Longitude, Latitude, loc_id) %>%
  arrange(loc_id)

bb <- st_bbox(boundary_m)
rw <- c(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])

v <- deldir(
  loc_data_all_m[, 1],
  loc_data_all_m[, 2],
  rw = rw
)

tiles <- tile.list(v)

V_polys <- lapply(seq_along(tiles), function(i) {
  xy <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  if (!all(xy[1, ] == xy[nrow(xy), ])) {
    xy <- rbind(xy, xy[1, ])
  }
  st_polygon(list(xy))
})

V_sf <- st_sf(
  loc_id = loc_data_all_map$loc_id,
  geometry = st_sfc(V_polys, crs = 27700)
) %>%
  st_make_valid() %>%
  st_intersection(st_make_valid(boundary_m)) %>%
  st_make_valid()

mesh <- inla.mesh.2d(
  boundary = inla.sp2segment(as_Spatial(boundary_m)),
  max.edge = c(100, 500),
  cutoff   = 80,
  offset   = c(100, 300)
)
# plot(mesh)

fm_crs(mesh) <- st_crs(V_sf)
V_sf <- V_sf %>%
  mutate(
    sampler_id = row_number(),
    weight = 1
  )
# plot(V_sf)

ips_sf <- fmesher::fm_int(
  domain   = mesh,
  samplers = V_sf,
  name     = "x",
  int.args = list(method = "direct", nsub2 = nsub2),
  format   = "sf"
)

triangles <- lapply(seq_len(nrow(mesh$graph$tv)), function(i) {
  idx <- mesh$graph$tv[i, ]
  coords <- mesh$loc[idx, ]
  coords <- rbind(coords, coords[1, ])
  st_polygon(list(coords))
})
mesh_sf <- st_sf(geometry = st_sfc(triangles, crs = st_crs(V_sf)))

block_map <- V_sf %>%
  st_drop_geometry() %>%
  select(sampler_id, loc_id)

ips_sf <- ips_sf %>%
  left_join(
    block_map,
    by = c(".block" = "sampler_id")
  ) %>%
  mutate(
    .block = loc_id
  )

coords_mat <- st_coordinates(ips_sf)

ips <- ips_sf %>%
  st_drop_geometry() %>%
  mutate(
    xcoord = coords_mat[, 1],
    ycoord = coords_mat[, 2]
  )

loc_levels <- sort(unique(ips$.block))
K <- length(loc_levels)

loc_map <- tibble(
  loc_id = loc_levels,
  idx = seq_len(K)
)

ips <- ips %>%
  left_join(
    loc_map,
    by = c(".block" = "loc_id")
  ) %>%
  mutate(
    .block = idx
  ) %>%
  select(-idx)

model_years <- 2023:2025
n_rep <- length(model_years)

counts <- data %>%
  mutate(
    Year = as.numeric(Year),
    rep = as.integer(factor(Year, levels = model_years))
  ) %>%
  group_by(loc_id, rep) %>%
  summarise(
    crime = sum(crime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  inner_join(
    loc_map,
    by = "loc_id"
  ) %>%
  select(idx, rep, crime)

grid <- expand.grid(
  idx = seq_len(K),
  rep = seq_len(n_rep)
) %>%
  as_tibble() %>%
  left_join(
    counts,
    by = c("idx", "rep")
  ) %>%
  mutate(
    crime = replace_na(crime, 0L),
    .block = idx + K * (rep - 1L)
  ) %>%
  arrange(.block)

y <- grid$crime

summary(y)
table(y == 0)

ips_final_base <- purrr::map_dfr(seq_len(n_rep), function(r) {
  ips %>%
    mutate(
      rep = r,
      Year = model_years[r],
      .block = .block + K * (r - 1L)
    )
})

# spde <- inla.spde2.matern(mesh, alpha = 2)

p_sigma <- 0.01
p_range <- 0.01

sigma_median <- 1
range_median <- 50

mu_sigma <- log(sigma_median)
sd_sigma <- (log(5) - mu_sigma) / qnorm(1 - p_sigma)

mu_range <- log(range_median)
sd_range <- (mu_range - log(10)) / qnorm(1 - p_range)

d <- 2
nu <- 1
alpha <- 2

C <- gamma(nu) / (gamma(alpha) * (4 * pi)^(d / 2))

B.kappa <- cbind(log(sqrt(8 * nu)), 0, -1)

B.tau <- cbind(
  0.5 * log(C) - nu * log(sqrt(8 * nu)),
  -1,
  nu
)

spde_model <- INLA::inla.spde2.matern(
  mesh = mesh,
  alpha = alpha,
  B.tau = B.tau,
  B.kappa = B.kappa,
  theta.prior.mean = c(mu_sigma, mu_range),
  theta.prior.prec = diag(c(
    1 / sd_sigma^2,
    1 / sd_range^2
  ))
)


spatial_ips_plot <- ips_sf %>%
  st_transform(4326) %>%
  mutate(
    spatial_w = weight
  )

key_list <- c("amenity", "highway", "man_made", "railway", "shop")

value_list <- list(
  amenity = c("bar", "nightclub", "bus_station"), # bank, "police"
  railway = c("subway_entrance"),
  shop = c("supermarket", "convenience")  #"department_store"
)

loc_ips <- cbind(ips$xcoord, ips$ycoord)
colnames(loc_ips) <- c("Longitude", "Latitude")

for (key in key_list) {
  for (value in c(value_list[[key]])) {
    
    if (exists(value, envir = .GlobalEnv)) {
      rm(list = value, envir = .GlobalEnv)
    }
    
    message("Processing spatial covariate: ", Region, ": ", key, "_", value)
    
    path_cov1 <- file.path(prepare_dir, "cov/loc")
    path_cov2 <- file.path(prepare_dir, "cov/dist/Spatial", crime_type)
    
    if (!dir.exists(path_cov1)) {
      dir.create(path_cov1, recursive = TRUE)
    }
    
    if (!dir.exists(path_cov2)) {
      dir.create(path_cov2, recursive = TRUE)
    }
    path_dist <- paste0(path_cov2, "/", key, "_", value, "_Spatial_", nsub2, ".RData")
    
    if (file.exists(path_dist)) {
      
      load(path_dist)
      assign(value, cov_dist_geo)
      
    } else {
      
      message("No existing path for: ", key, " : ", value)
      
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
        
        save(loc_cov, file = path_loc)
      }
      
      loc_cov <- st_coordinates(
        st_transform(
          st_sfc(
            st_multipoint(as.matrix(loc_cov)),
            crs = 4326
          ),
          27700
        )
      )[, 1:2]
      
      colnames(loc_cov) <- c("Longitude", "Latitude")
      
      aa <- fields::rdist(loc_cov, loc_ips)
      
      if (is.null(dim(aa))) {
        cov_dist_geo <- aa
      } else {
        cov_dist_geo <- apply(aa, 2, min)
      }
      
      save(cov_dist_geo, file = path_dist)
      assign(value, cov_dist_geo)
    }
  }
}

covariate_store <- list()

for (key in key_list) {
  for (value in c(value_list[[key]])) {
    
    if (exists(value, envir = .GlobalEnv)) {
      rm(list = value, envir = .GlobalEnv)
    }
    
    message("Loading spatial covariate: ", Region, ": ", key, "_", value)
    
    path_cov2 <- file.path(prepare_dir, "cov/dist/Spatial", crime_type)
    
    path_dist <- paste0(path_cov2, "/", key, "_", value, "_Spatial_", nsub2, ".RData")
    
    if (file.exists(path_dist)) {
      
      load(path_dist)
      
      cov_dist_geo <- exp(-cov_dist_geo / 1000)
      cov_dist_geo_rep <- rep(cov_dist_geo, n_rep)
      
      assign(value, cov_dist_geo_rep)
      covariate_store[[value]] <- cov_dist_geo_rep
      
    } else {
      
      message("No existing covariate file for: ", key, " : ", value)
      next
    }
  }
}

covariate_names <- names(covariate_store)
covariate_df <- as_tibble(covariate_store)
ips_cov <- bind_cols(
  ips_final_base,
  covariate_df
)
ips_final <- ips_cov

agg <- bru_mapper_aggregate(
  rescale = FALSE,
  n_block = length(unique(ips_final$.block)),
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
formula

lik <- bru_obs(
  formula = formula,
  response_data = data.frame(y = y),
  family = "poisson",
  data = ips_final,
  allow_combine = TRUE
)

fixed_terms <- covariate_names[
  covariate_names %in% names(ips_final)
]

cmp_terms <- c(
  fixed_terms,
  "spde(cbind(xcoord, ycoord), model = spde_model, mapper = bru_mapper(mesh), replicate = rep)"
)

cmp_str <- paste(
  "y ~ Intercept(1) +",
  paste(cmp_terms, collapse = " + ")
)

cmp <- as.formula(cmp_str)
cmp

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
save(fit, file = file.path(spatial_dir, paste0("fit_", crime_type, "_nsub", nsub2, ".RData")))
message("Spatial model saved in: ", spatial_dir)
