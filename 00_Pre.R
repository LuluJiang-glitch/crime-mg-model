# Crime Data
# library(tidyverse)
# base_path <- "~/Desktop/Crime/Data/Police_Aichive_all"
# folders <- list.dirs(base_path, recursive = FALSE)
# 
# read_month <- function(folder) {
#   files <- list.files(folder, pattern = "street.*\\.csv$", full.names = TRUE)
#   
#   map_dfr(files, read_csv, show_col_types = FALSE) %>%
#     select(Month, Longitude, Latitude, `Crime type`) %>%
#     drop_na() %>%
#     count(Month, Longitude, Latitude, `Crime type`, name = "count")
# }
# 
# monthly_data <- map_dfr(folders, read_month)
# yearly_data <- monthly_data %>%
#   mutate(Year = substr(Month, 1, 4)) %>%
#   group_by(Year, Longitude, Latitude, `Crime type`) %>%
#   summarise(count = sum(count), .groups = "drop")
# save(yearly_data, file = "~/Desktop/Crime/code_final/RData/Prepare/01_crime_yearly_data_all.RData")
load("~/Desktop/Crime/code_final/RData/Prepare/01_crime_yearly_data_all.RData")



# graph
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(osmdata)
library(dplyr)
library(MetricGraph)

source("~/Desktop/Crime_UK/Code/ExampleCode/LinesMetricGraph-main/metric_graph.R")
uk <- ne_states(country = "united kingdom", returnclass = "sf")
london <- uk[uk$region == "Greater London", ]
aa <- london[london$name %in% "City", ]
boundary <- st_transform(aa, crs = 4326) %>% st_union()
boundary_sf <- st_as_sf(boundary)

path <- "~/Desktop/Crime/code_final/RData/Prepare/02_graph_City.RData"
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(path)) {
  call <- st_bbox(boundary) %>%
    opq() %>%
    add_osm_feature(key = "highway")
  
  lines_inbound <- osmdata_sp(call)$osm_lines %>%
    st_as_sf() %>%
    dplyr::select(geometry) %>%
    st_intersection(boundary_sf) %>%
    as_Spatial()
  
  graph_pre <- graph_components$new(edges = lines_inbound, perform_merges = TRUE)
  graph <- graph_pre$get_largest()
  
  summary(graph)
  save(graph, boundary, boundary_sf, file = path)
} else {
  load(path)
}

# graph$plot(vertex_size = 0)







