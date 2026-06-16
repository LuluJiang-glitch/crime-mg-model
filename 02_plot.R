library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

base_dir <- "~/Desktop/Crime/code_final"
path_spatial <- file.path(base_dir, "RData", "Spatial")
path_mg <- file.path(base_dir, "RData", "MG")
path_point <- file.path(base_dir, "RData", "Point")
out_fig <- file.path(base_dir, "Fig")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

crime_types <- c(
  "Drugs",
  "Bicycle theft",
  "Robbery",
  "Theft from the person"
)

crime_labels <- c(
  "Drugs" = "Drugs",
  "Bicycle theft" = "Bicycle theft",
  "Robbery" = "Robbery",
  "Theft from the person" = "Person theft"
)

term_labels <- c(
  "bar" = "Bar",
  "nightclub" = "Nightclub",
  "bus_station" = "Bus station",
  "subway_entrance" = "Subway entrance",
  "supermarket" = "Supermarket",
  "convenience" = "Convenience store"
)

term_order <- c(
  "bar",
  "nightclub",
  "bus_station",
  "subway_entrance",
  "supermarket",
  "convenience"
)


get_result <- function(file, model, ct) {
  
  if (!file.exists(file)) {
    warning("File does not exist: ", file)
    return(NULL)
  }
  
  e <- new.env()
  load(file, envir = e)
  
  if (!exists("fit", envir = e)) {
    warning("No object named 'fit' in file: ", file)
    return(NULL)
  }
  
  fit <- e$fit
  
  if (is.null(fit$summary.fixed)) {
    fixed <- NULL
  } else {
    
    fixed <- as.data.frame(fit$summary.fixed) %>%
      mutate(
        Crime_Type = ct,
        Model = model,
        term = rownames(fit$summary.fixed),
        lower = `0.025quant`,
        upper = `0.975quant`,
        sig = ifelse(lower * upper > 0, "CI excludes 0", "CI includes 0")
      ) %>%
      select(
        Crime_Type,
        Model,
        term,
        mean,
        lower,
        upper,
        sig
      )
  }
  
  cpo <- fit$cpo$cpo
  
  metrics <- data.frame(
    Crime_Type = ct,
    Model = model,
    NlogLik = if (!is.null(fit$mlik)) -fit$mlik[1] else NA_real_,
    DIC = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
    WAIC = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_,
    LCPO = if (!is.null(cpo)) sum(-log(cpo), na.rm = TRUE) else NA_real_
  )
  
  list(
    fixed = fixed,
    metrics = metrics
  )
}

all_fixed <- data.frame()
all_metrics <- data.frame()

for (ct in crime_types) {
  
  file_mg <- file.path(
    path_mg,
    paste0("fit_", ct, "_M30.RData")
  )
  
  file_spatial <- file.path(
    path_spatial,
    paste0("fit_", ct, "_nsub5.RData")
  )
  
  file_point <- file.path(
    path_point,
    paste0("fit_", ct, ".RData")
  )
  
  res_mg <- get_result(file_mg, "MG", ct)
  res_sp <- get_result(file_spatial, "Spatial", ct)
  res_pt <- get_result(file_point, "Point", ct)
  
  if (!is.null(res_mg)) {
    if (!is.null(res_mg$fixed)) {
      all_fixed <- bind_rows(all_fixed, res_mg$fixed)
    }
    all_metrics <- bind_rows(all_metrics, res_mg$metrics)
  }
  
  if (!is.null(res_sp)) {
    if (!is.null(res_sp$fixed)) {
      all_fixed <- bind_rows(all_fixed, res_sp$fixed)
    }
    all_metrics <- bind_rows(all_metrics, res_sp$metrics)
  }
  
  if (!is.null(res_pt)) {
    if (!is.null(res_pt$fixed)) {
      all_fixed <- bind_rows(all_fixed, res_pt$fixed)
    }
    all_metrics <- bind_rows(all_metrics, res_pt$metrics)
  }
}

all_fixed <- all_fixed %>%
  mutate(
    Model = factor(Model, levels = c("MG", "Spatial", "Point")),
    Crime_Type = factor(Crime_Type, levels = crime_types),
    
    term = case_when(
      term == "amenity_bar" ~ "bar",
      term == "amenity_nightclub" ~ "nightclub",
      term == "amenity_bus_station" ~ "bus_station",
      term == "railway_subway_entrance" ~ "subway_entrance",
      term == "shop_supermarket" ~ "supermarket",
      term == "shop_convenience" ~ "convenience",
      
      term == "indexcov23" ~ "cov23",
      term == "cov23" ~ "cov23",
      
      TRUE ~ term
    ),
    
    sig = factor(
      sig,
      levels = c("CI includes 0", "CI excludes 0")
    )
  ) %>%
  filter(
    !term %in% c(
      "Intercept",
      "index2023",
      "indexcov23",
      "cov23"
    )
  ) %>%
  filter(
    term %in% term_order
  ) %>%
  mutate(
    term = factor(
      term,
      levels = rev(term_order),
      labels = rev(term_labels[term_order])
    )
  )

all_metrics <- all_metrics %>%
  mutate(
    Model = factor(Model, levels = c("MG", "Spatial", "Point")),
    Crime_Type = factor(Crime_Type, levels = crime_types)
  )


pd_fe <- position_dodge(width = 0.75)

p_fe <- ggplot(
  all_fixed,
  aes(
    x = mean,
    y = term,
    color = Model
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey40"
  ) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper,
      group = Model
    ),
    position = pd_fe,
    height = 0.18,
    linewidth = 0.75,
    linetype = "solid"
  ) +
  geom_point(
    aes(
      shape = sig,
      group = Model
    ),
    position = pd_fe,
    size = 3.4,
    stroke = 1.1
  ) +
  scale_shape_manual(
    values = c(
      "CI includes 0" = 1,
      "CI excludes 0" = 16
    ),
    labels = c(
      "CI includes 0" = "CI includes 0",
      "CI excludes 0" = "CI excludes 0"
    ),
    name = "95% CI"
  ) +
  facet_grid(
    . ~ Crime_Type,
    scales = "free_x",
    space = "fixed",
    labeller = labeller(Crime_Type = crime_labels)
  ) +
  theme_bw(base_size = 20) +
  theme(
    strip.text.x = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 19),
    panel.spacing.x = unit(0.8, "lines"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Model"
  )

p_fe

ggsave(
  filename = file.path(out_fig, "FE_all.png"),
  plot = p_fe,
  width = 16,
  height = 8,
  dpi = 300
)

plot_metrics <- all_metrics %>%
  mutate(
    Crime_Type_short = factor(
      crime_labels[as.character(Crime_Type)],
      levels = crime_labels[crime_types]
    )
  ) %>%
  pivot_longer(
    cols = c(DIC, WAIC, LCPO, NlogLik),
    names_to = "Metric",
    values_to = "Value"
  )

p_metrics <- ggplot(
  plot_metrics,
  aes(
    x = Crime_Type_short,
    y = Value,
    fill = Model
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol = 2
  ) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    legend.position = "right",
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 19),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Value",
    fill = "Model"
  )

p_metrics

ggsave(
  filename = file.path(out_fig, "Model_Comparison_Metrics.png"),
  plot = p_metrics,
  width = 13,
  height = 7,
  dpi = 300
)
