library(tidyverse)
library(Shennong)
library(glue)

# TabulaSapiens2022 -----------------------------------------------------------
alignment_dir <- "/home/sduan/data/TabulaSapiens/data/alignments/10x/"
samples <- list.files(alignment_dir)
metrics <- glue("{alignment_dir}/{samples}/outs/metrics_summary.csv")
names(metrics) <- samples
metrics <- metrics[file.exists(metrics)]
head(metrics)

TabulaSapiens2022_metrics_df <- map_dfr(names(metrics), function(x) {
  print(x)
  # x <- "TSP2_BM_vertebralbody_1_2_5prime_TCR"
  df <- read_csv(metrics[[x]])
  df$sample <- x
  return(df)
})
TabulaSapiens2022_metrics_df$study <- "TabulaSapiens2022"

# StevenBWells2025 -------------------------------------------------------
alignment_dir <- "/home/sduan/data/StevenBWells2025/alignments/count/"
samples <- list.dirs(alignment_dir, full.names = FALSE, recursive = FALSE)
metrics <- glue("{alignment_dir}/{samples}/outs/metrics_summary.csv")
names(metrics) <- samples
metrics <- metrics[file.exists(metrics)]
head(metrics)

StevenBWells2025_metrics_df <- map_dfr(names(metrics), function(x) {
  print(x)
  df <- read_csv(metrics[[x]])
  df$sample <- x
  return(df)
})
StevenBWells2025_metrics_df$study <- "StevenBWells2025"

# ZhenlongLi2024 -----------------------------------------------------
alignment_dir <- "/home/sduan/data/ZhenlongLi2024/alignments/"
samples <- list.dirs(alignment_dir, full.names = FALSE, recursive = FALSE)
metrics <- glue("{alignment_dir}/{samples}/outs/metrics_summary.csv")
names(metrics) <- samples
metrics <- metrics[file.exists(metrics)]
metrics
ZhenlongLi2024_metrics_df <- map_dfr(names(metrics), function(x) {
  print(x)
  df <- read_csv(metrics[[x]])
  df$sample <- x
  return(df)
})
ZhenlongLi2024_metrics_df$study <- "ZhenlongLi2024"

# Combine all metrics -----------------------------------------------------
all_metrics_df <- bind_rows(
  ZhenlongLi2024_metrics_df,
  TabulaSapiens_metrics_df,
  StevenBWells2025_metrics_df
)
all_metrics_df$`Estimated Number of Cells` |> sum()

all_metrics_df |>
  ggplot(mapping = aes(
    x = `Median Genes per Cell`/1000, 
    y = `Estimated Number of Cells`/1000, color = study)) +
  geom_point() +
  catplot::theme_cat(show_panel_grid_marjor = "both",
                     show_panel_grid_minor = "both",
                     aspect_ratio = 1) +
  scale_color_brewer(palette = "Paired") 
