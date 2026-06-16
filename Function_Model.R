merge_points_within_distance <- function(loc_df, dist_m = 30, crs_longlat = 4326, crs_meter = 27700) {
  if (!all(c("Longitude", "Latitude") %in% names(loc_df))) {
    stop("loc_df must contain Longitude and Latitude.")
  }

  loc_df <- loc_df %>% mutate(original_loc_id = row_number())
  loc_sf <- st_as_sf(loc_df, coords = c("Longitude", "Latitude"), crs = crs_longlat, remove = FALSE)
  loc_sf_m <- st_transform(loc_sf, crs_meter)

  n <- nrow(loc_sf_m)
  parent <- seq_len(n)

  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }

  union_root <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) parent[rb] <<- ra
  }

  near_list <- st_is_within_distance(loc_sf_m, loc_sf_m, dist = units::set_units(dist_m, "m"))
  for (i in seq_along(near_list)) {
    js <- near_list[[i]]
    js <- js[js > i]
    if (length(js) > 0) for (j in js) union_root(i, j)
  }

  cluster_id <- as.integer(factor(vapply(seq_len(n), find_root, integer(1))))
  coord_m <- st_coordinates(loc_sf_m)

  cluster_centers_m <- tibble(
    original_loc_id = loc_df$original_loc_id,
    cluster_id = cluster_id,
    X = coord_m[, "X"],
    Y = coord_m[, "Y"]
  ) %>%
    group_by(cluster_id) %>%
    summarise(X = mean(X), Y = mean(Y), n_original_points = n(), .groups = "drop")

  cluster_centers_sf_m <- st_as_sf(cluster_centers_m, coords = c("X", "Y"), crs = crs_meter, remove = FALSE)
  cluster_centers_sf_ll <- st_transform(cluster_centers_sf_m, crs_longlat)
  center_ll <- st_coordinates(cluster_centers_sf_ll)

  loc_merged_df <- cluster_centers_sf_ll %>%
    st_drop_geometry() %>%
    mutate(Longitude = center_ll[, "X"], Latitude = center_ll[, "Y"], loc_id = cluster_id) %>%
    select(loc_id, Longitude, Latitude, n_original_points) %>%
    arrange(loc_id)

  loc_cluster_map <- loc_df %>%
    mutate(cluster_id = cluster_id) %>%
    left_join(
      loc_merged_df %>% select(loc_id, Longitude_merged = Longitude, Latitude_merged = Latitude),
      by = c("cluster_id" = "loc_id")
    ) %>%
    transmute(
      original_loc_id,
      Longitude_original = Longitude,
      Latitude_original = Latitude,
      loc_id = cluster_id,
      Longitude = Longitude_merged,
      Latitude = Latitude_merged
    )

  list(
    loc_merged_df = loc_merged_df,
    loc_cluster_map = loc_cluster_map,
    n_original = n,
    n_merged = nrow(loc_merged_df)
  )
}




# cov
extract_covariates <- function(boundary, key, value) {
  query <- opq(bbox = st_bbox(boundary)) %>%
    add_osm_feature(key = key, value = value) %>%
    add_osm_feature(key = "name")

  aa <- osmdata_sf(query)
  if (is.null(aa$osm_points) || nrow(aa$osm_points) == 0) {
    message("Warning: no OSM point data found for ", key, " = ", value)
    return(NULL)
  }

  aa_points <- aa$osm_points %>% st_as_sf() %>% st_transform(crs = st_crs(boundary))
  coords <- unique(as.matrix(st_coordinates(aa_points$geometry)))
  coords_sf <- st_as_sf(as.data.frame(coords), coords = c(1, 2), crs = 4326)
  coords_in_boundary <- suppressWarnings(st_intersection(coords_sf, boundary))

  if (nrow(coords_in_boundary) == 0) return(NULL)
  unique(as.matrix(st_coordinates(coords_in_boundary$geometry)))
}





# functions for helping compute geodist
check_unique_loc2pte <- function(graph, loc_data) {
  pte_data <- graph$coordinates(XY = loc_data)
  if (nrow(loc_data) != nrow(unique(pte_data))) {
    pte_data <- unique(pte_data)
    loc_data <- graph$coordinates(PtE = pte_data)
  }

  tibble(
    Longitude = loc_data[, 1],
    Latitude = loc_data[, 2],
    edge_number = pte_data[, 1],
    distance_on_edge = pte_data[, 2]
  )
}

choose_unique <- function(data) {
  data_unique <- unique(data)
  idx <- match(do.call(paste, as.data.frame(data)), do.call(paste, as.data.frame(data_unique)))
  list(data_unique = data_unique, idx = idx)
}

choose_A_notin_B <- function(A, B) {
  A <- as.data.frame(A)
  B <- as.data.frame(B)
  colnames(A) <- c("A1", "A2")
  colnames(B) <- c("B1", "B2")
  str_A <- apply(A, 1, paste, collapse = "_")
  str_B <- apply(B, 1, paste, collapse = "_")
  A[!(str_A %in% str_B), , drop = FALSE]
}

check_pte_loc_unique <- function(pte = NULL, graph) {
  pte0 <- pte
  pte0u <- unique(pte0)
  index01 <- match(apply(pte0, 1, paste, collapse = "_"), apply(pte0u, 1, paste, collapse = "_"))

  loc0 <- graph$coordinates(PtE = pte0u)
  loc0u <- unique(loc0)
  index02 <- match(apply(loc0, 1, paste, collapse = "_"), apply(loc0u, 1, paste, collapse = "_"))

  if (nrow(pte0u) != nrow(loc0u)) {
    pte1 <- graph$coordinates(XY = loc0u)
    pte1u <- unique(pte1)
    index03 <- match(apply(pte1, 1, paste, collapse = "_"), apply(pte1u, 1, paste, collapse = "_"))
    index <- index03[index02][index01]
    pte_unique <- pte1u
  } else {
    index <- index02[index01]
    pte_unique <- pte0u
  }

  loc_unique <- loc0u
  if (nrow(pte_unique) != nrow(loc_unique)) {
    stop("nrow(pte_unique) != nrow(loc_unique)")
  }

  list(pte_unique = pte_unique, loc_unique = loc_unique, index = index)
}


compute_geo_matdist <- function(type = "geo", loc_mesh = NULL, pte_mesh = NULL, loc_data = NULL, graph) {
  aa <- choose_A_notin_B(loc_data, loc_mesh)

  if (nrow(aa) != 0) {
    loc_data1 <- as.matrix(choose_unique(aa)$data_unique)
    colnames(loc_data1) <- colnames(loc_mesh)
    loc_test <- rbind(loc_data1, loc_mesh)

    pte_data <- graph$coordinates(XY = loc_data1)
    pte_data1 <- as.matrix(choose_unique(choose_A_notin_B(pte_data, pte_mesh))$data_unique)
    pte_test <- rbind(pte_data1, pte_mesh)

    dist_mat <- graph$compute_geodist_PtE(
      PtE = pte_test,
      normalized = TRUE,
      include_vertices = FALSE,
      verbose = 1
    )

    idx1 <- match(do.call(paste, as.data.frame(pte_data)), do.call(paste, as.data.frame(pte_test)))
    idx1_with_mesh <- c(idx1, (nrow(pte_data1) + 1):nrow(pte_test))
    dist_mat1 <- dist_mat[idx1_with_mesh, (nrow(pte_data1) + 1):nrow(pte_test)]

    idx2 <- match(do.call(paste, as.data.frame(loc_data)), do.call(paste, as.data.frame(loc_test)))
    dist_mat1 <- dist_mat1[idx2, ]
  } else {
    dist_mat <- graph$compute_geodist_PtE(
      PtE = pte_mesh,
      normalized = TRUE,
      include_vertices = FALSE,
      verbose = 1
    )
    loc_mesh_df <- as.data.frame(loc_mesh)
    loc_data_df <- as.data.frame(loc_data)
    colnames(loc_mesh_df) <- c("Longitude", "Latitude")
    colnames(loc_data_df) <- c("Longitude", "Latitude")
    str_data <- apply(loc_data_df, 1, paste, collapse = "_")
    str_mesh <- apply(loc_mesh_df, 1, paste, collapse = "_")
    col_idx <- which(str_mesh %in% str_data)
    dist_mat1 <- dist_mat[col_idx, , drop = FALSE]
  }

  loc_mesh_df <- as.data.frame(loc_mesh)
  loc_data_df <- as.data.frame(loc_data)
  colnames(loc_mesh_df) <- c("Longitude", "Latitude")
  colnames(loc_data_df) <- c("Longitude", "Latitude")
  str_data <- apply(loc_data_df, 1, paste, collapse = "_")
  str_mesh <- apply(loc_mesh_df, 1, paste, collapse = "_")
  row_idx <- which(str_data %in% str_mesh)
  col_idx <- which(str_mesh %in% str_data)
  if (length(row_idx) > 0 && length(col_idx) > 0) dist_mat1[row_idx, col_idx] <- 0

  dist_mat1
}


simpson_u_weight <- function(L, step, K = NULL) {
  if (is.null(K)) K <- ceiling(L / step)
  if (K < 2) K <- 2
  if (K %% 2 == 1) K <- K + 1

  u <- seq(0, 1, length.out = K + 1)
  coef <- rep(2, K + 1)
  coef[1] <- 1
  coef[K + 1] <- 1
  coef[seq(2, K, by = 2)] <- 4

  list(u = u, w = (L / K / 3) * coef)
}


make_ips <- function(edge_sf, graph, step = 0.0002, K = NULL) {
  stopifnot(".edge_lengths" %in% names(edge_sf), "data_after" %in% names(edge_sf))

  block_map <- match(edge_sf$data_after, sort(unique(edge_sf$data_after)))
  idx_all <- integer(0)
  u_all <- numeric(0)
  w_all <- numeric(0)
  b_all <- integer(0)

  for (i in seq_len(nrow(edge_sf))) {
    rr <- simpson_u_weight(edge_sf$.edge_lengths[i], step = step, K = K)
    m <- length(rr$u)
    idx_all <- c(idx_all, rep.int(i, m))
    u_all <- c(u_all, rr$u)
    w_all <- c(w_all, rr$w)
    b_all <- c(b_all, rep.int(block_map[i], m))
  }

  tibble(x = as_MGG(cbind(idx_all, u_all), graph = graph), weight = w_all, .block = b_all)
}
