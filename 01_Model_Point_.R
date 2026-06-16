# Point model for City of London crime data, 2023-2025.

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

name <- Region <- "City"
crime_type <- "Bicycle theft"
# "Theft from the person", "Robbery", "Drugs", "Bicycle theft"
# "Anti-social behaviour", "Criminal damage and arson", "Violence and sexual offences", "Vehicle crime"
# "Burglary", "Shoplifting"

max_edge <- c(100, 500)
cutoff <- 80
offset <- c(100, 300)

out_dir <- point_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

uk <- ne_states(country = "united kingdom", returnclass = "sf")
london <- uk[uk$region == "Greater London", ]
boundary_ll <- st_transform(london[london$name %in% c(Region), ], crs = 4326) %>%
  st_union()
boundary_ll <- st_as_sf(boundary_ll)
boundary <- st_transform(boundary_ll, 27700)
boundary <- st_make_valid(boundary)

load(file.path(prepare_dir, "01_crime_yearly_data_all.RData"))

data_all_raw <- yearly_data %>%
  filter(`Crime type` == crime_type) %>%
  filter(Year %in% c("2023", "2024", "2025")) %>%
  mutate(Year = as.numeric(Year)) %>%
  filter(Year %in% 2023:2025) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  filter(rowSums(st_within(geometry, boundary_ll, sparse = FALSE)) > 0) %>%
  mutate(
    Longitude = st_coordinates(.)[, 1],
    Latitude = st_coordinates(.)[, 2]
  ) %>%
  group_by(Longitude, Latitude, Year) %>%
  summarise(crime = sum(count, na.rm = TRUE), .groups = "drop") %>%
  st_drop_geometry()

loc_data_all_original <- data_all_raw %>%
  distinct(Longitude, Latitude) %>%
  arrange(Longitude, Latitude) %>%
  mutate(original_loc_id = row_number())

merge_res <- merge_points_within_distance(
  loc_df = loc_data_all_original %>% select(Longitude, Latitude),
  dist_m = 30,
  crs_longlat = 4326,
  crs_meter = 27700
)

message(
  "Coordinates merged: original points = ", merge_res$n_original,
  "; merged points = ", merge_res$n_merged,
  "; threshold = ", 30, " meters."
)

loc_cluster_map <- loc_data_all_original %>%
  select(
    original_loc_id,
    Longitude_original = Longitude,
    Latitude_original = Latitude
  ) %>%
  left_join(
    merge_res$loc_cluster_map %>%
      select(original_loc_id, loc_id, Longitude, Latitude),
    by = "original_loc_id"
  )

loc_unique_ll <- merge_res$loc_merged_df %>%
  select(loc_id, Longitude, Latitude, n_original_points) %>%
  arrange(loc_id)

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
  summarise(crime = sum(crime, na.rm = TRUE), .groups = "drop") %>%
  arrange(loc_id, Year)

loc_unique_sf <- st_as_sf(
  loc_unique_ll,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(27700)

loc_unique_m <- st_coordinates(loc_unique_sf)
loc_unique_ll <- loc_unique_ll %>%
  mutate(
    xcoord = loc_unique_m[, 1],
    ycoord = loc_unique_m[, 2]
  )

K <- nrow(loc_unique_ll)

counts <- data_all %>%
  filter(Year %in% 2023:2025) %>%
  group_by(loc_id, Year) %>%
  summarise(crime = sum(crime, na.rm = TRUE), .groups = "drop")

grid <- tidyr::expand_grid(
  Year = 2023:2025,
  loc_id = loc_unique_ll$loc_id
) %>%
  mutate(rep = match(Year, 2023:2025)) %>%
  left_join(counts, by = c("loc_id", "Year")) %>%
  left_join(
    loc_unique_ll %>%
      select(loc_id, Longitude, Latitude, xcoord, ycoord, n_original_points),
    by = "loc_id"
  ) %>%
  mutate(crime = replace_na(crime, 0L)) %>%
  arrange(Year, loc_id)

y <- grid$crime

mesh <- inla.mesh.2d(
  boundary = inla.sp2segment(as_Spatial(boundary)),
  max.edge = max_edge,
  cutoff = cutoff,
  offset = offset
)
plot(mesh)
fm_crs(mesh) <- st_crs(boundary)

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


key_list <- c("amenity", "highway", "man_made", "railway", "shop")

value_list <- list(
  amenity = c("bar", "nightclub", "bus_station"), # bank, "police"
  railway = c("subway_entrance"),
  shop = c("supermarket", "convenience")  #"department_store"
)

loc_pts <- as.matrix(loc_unique_ll[, c("xcoord", "ycoord")])
colnames(loc_pts) <- c("Longitude", "Latitude")

covariate_store <- list()

for (key in key_list) {
  vals <- value_list[[key]]
  if (is.null(vals)) next
  
  for (value in vals) {
    message("Processing point covariate: ", Region, " : ", key, "_", value)
    
    path_cov1 <- file.path(prepare_dir, "cov/loc")
    path_cov2 <- file.path(prepare_dir, "cov/dist/Point", crime_type)
    
    if (!dir.exists(path_cov1)) dir.create(path_cov1, recursive = TRUE)
    if (!dir.exists(path_cov2)) dir.create(path_cov2, recursive = TRUE)
    
    path_dist <- paste0(path_cov2, "/", key, "_", value, "_Point.RData")
    
    if (file.exists(path_dist)) {
      load(path_dist)
    } else {
      path_loc <- paste0(path_cov1, "/", key, "_", value, "_loc.RData")
      
      if (!file.exists(path_loc)) {
        if (exists("extract_covariates", mode = "function")) {
          message("No location file found; extracting OSM covariate locations for ", key, " = ", value)
          loc_cov_raw <- extract_covariates(
            boundary = boundary_ll,
            key = key,
            value = value
          )
          
          if (is.null(loc_cov_raw) || nrow(as.data.frame(loc_cov_raw)) == 0) {
            message("No OSM data found for ", key, " = ", value, "; skipping.")
            next
          }
          
          loc_cov <- as.matrix(as.data.frame(loc_cov_raw)[, 1:2])
          save(loc_cov, file = path_loc)
        } else {
          message("No location file and extract_covariates() not available for ", key, " = ", value, "; skipping.")
          next
        }
      } else {
        load(path_loc)
      }
      
      loc_cov_m <- st_coordinates(
        st_transform(
          st_sfc(st_multipoint(as.matrix(loc_cov)), crs = 4326),
          27700
        )
      )[, 1:2]
      
      colnames(loc_cov_m) <- c("Longitude", "Latitude")
      
      aa <- fields::rdist(loc_cov_m, loc_pts)
      
      if (is.null(dim(aa))) {
        cov_dist_geo <- aa
      } else {
        cov_dist_geo <- apply(aa, 2, min)
      }
      
      save(cov_dist_geo, file = path_dist)
    }
    
    if (length(cov_dist_geo) != K) {
      stop("Covariate length mismatch for ", value, ": length = ", length(cov_dist_geo), ", K = ", K)
    }
    
    covariate_store[[value]] <- exp(-(cov_dist_geo / 1000))
  }
}

covariate_names <- names(covariate_store)

if (length(covariate_names) > 0) {
  covariate_df <- as_tibble(covariate_store) %>% mutate(loc_id = loc_unique_ll$loc_id)
  grid <- grid %>% left_join(covariate_df, by = "loc_id")
}

Nips_final <- grid

all_vars <- c(covariate_names)
all_vars <- all_vars[all_vars %in% names(Nips_final)]

lik <- bru_obs(
  formula = y ~ .,
  family = "poisson",
  response_data = data.frame(y = y),
  data = Nips_final,
  allow_combine = TRUE
)

fixed_terms <- covariate_names[
  covariate_names %in% names(Nips_final)
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
save(fit, file = file.path(point_dir, paste0("fit_", crime_type, ".RData")))
message("Point model saved in: ", point_dir)
