file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) > 0) {
  base_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  base_dir <- normalizePath(getwd())
}

setwd(base_dir)
prepare_dir <- file.path(base_dir, "RData", "Prepare")



metric_graph_file <- file.path(base_dir, "metric_graph.R")
if (!file.exists(metric_graph_file)) {
  download.file(
    "https://raw.githubusercontent.com/karinalilleborge/LinesMetricGraph/main/metric_graph.R",
    destfile = metric_graph_file,
    mode = "wb"
  )
}
source(metric_graph_file)

# Crime Data 01_data_city.RData
# library(tidyverse)
# base_path <- "" # The path you download the crime data
# folders <- list.dirs(base_path, recursive = FALSE)
# read_month <- function(folder) {
#   files <- list.files(folder, pattern = "street.*\\.csv$", full.names = TRUE)
#   
#   map_dfr(files, read_csv, show_col_types = FALSE) %>%
#     select(Month, Longitude, Latitude, `Crime type`) %>%
#     drop_na() %>%
#     count(Month, Longitude, Latitude, `Crime type`, name = "count")
# }
# monthly_data <- map_dfr(folders, read_month)
# yearly_data <- monthly_data %>%
#   mutate(Year = substr(Month, 1, 4)) %>%
#   group_by(Year, Longitude, Latitude, `Crime type`) %>%
#   summarise(count = sum(count), .groups = "drop")
# data_city <- yearly_data %>%
#   filter(Year %in% as.character(model_years)) %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
#   filter(rowSums(st_within(geometry, boundary, sparse = FALSE)) > 0) %>%
#   mutate(
#     Longitude = st_coordinates(.)[, 1],
#     Latitude = st_coordinates(.)[, 2],
#     Year = as.numeric(Year)
#   )
# save(data_city, file = file.path(prepare_dir, "01_data_city.RData"))




# graph
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(osmdata)
library(dplyr)
library(MetricGraph)

uk <- ne_states(country = "united kingdom", returnclass = "sf")
london <- uk[uk$region == "Greater London", ]
aa <- london[london$name %in% "City", ]
boundary <- st_transform(aa, crs = 4326) %>% st_union()
boundary_sf <- st_as_sf(boundary)

path <- file.path(prepare_dir, "02_graph_City_L100.RData")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

graph_buffer_m <- 100

boundary_graph_sf <- boundary_sf %>%
  st_transform(27700) %>%
  st_buffer(dist = graph_buffer_m) %>%
  st_make_valid()

boundary_graph_osm <- boundary_graph_sf %>% st_transform(4326)

ggplot() +
  geom_sf(
    data = boundary_graph_osm,
    fill = NA,
    color = "red",
    linewidth = 1
  ) +
  geom_sf(
    data = boundary_sf,
    fill = NA,
    color = "blue",
    linewidth = 1
  ) +
  theme_minimal() +
  labs(
    title = "Original and Buffered Boundaries",
    subtitle = paste0("Buffer = ", graph_buffer_m, " m"),
    color = NULL
  )


if (!file.exists(path)) {
  call <- boundary_graph_osm %>%
    st_bbox() %>%
    opq() %>%
    add_osm_feature(key = "highway")
  
  lines_inbound1 <- osmdata_sf(call)$osm_lines
  
  lines_inbound2 <- lines_inbound1 %>%
    st_as_sf() %>%
    st_transform(4326)
  
  lines_inbound <- lines_inbound2 %>%
    dplyr::select(geometry) %>%
    st_intersection(boundary_graph_osm) %>%
    as_Spatial()
  
  graph_pre <- graph_components$new(
    edges = lines_inbound,
    perform_merges = TRUE
  )
  
  graph <- graph_pre$get_largest()
  graph$plot(type = "mapview", vertex_size = 0) + mapview::mapview(
    boundary_sf,
    color = "red",
    col.regions = NA,
    lwd = 3,
    layer.name = "Original boundary"
  )
  
  summary(graph)
  save(graph, file = path)
  
} else {
  load(path)
}
