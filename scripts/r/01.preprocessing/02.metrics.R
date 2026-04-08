library(tidyverse)
library(Shennong)
library(glue)

# CDominguezConde2022 -----------------------------------------------------
alignment_dir <- "/mnt/rstudio/CDominguezConde2022/alignments/count/"
samples <- list.files(alignment_dir)
metrics <- glue("{alignment_dir}/{samples}/outs/metrics_summary.csv")
names(metrics) <- samples

CDominguezConde2022_metrics_df <- map_dfr(samples, function(x) {
    print(x)
    df <- read_csv(metrics[[x]])
    df$sample <- x
    return(df)
})
CDominguezConde2022_metrics_df$study <- "CDominguezConde2022"

# NataliaJaeger2024 -------------------------------------------------------
alignment_dir <- "/mnt/public/NataliaJaeger2024/alignments/"
samples <- list.dirs(alignment_dir, full.names = FALSE, recursive = FALSE)
metrics <- glue("{alignment_dir}/{samples}/outs/metrics_summary.csv")
names(metrics) <- samples
metrics <- metrics[file.exists(metrics)]

NataliaJaeger2024_metrics_df <- map_dfr(samples, function(x) {
    print(x)
    df <- read_csv(metrics[[x]])
    df$sample <- x
    return(df)
})
NataliaJaeger2024_metrics_df$study <- "NataliaJaeger2024"

# AndrewDHildreth2021 -----------------------------------------------------
alignment_dir <- "/mnt/public/AndrewDHildreth2021/alignments/"
samples <- list.dirs(alignment_dir, full.names = FALSE, recursive = FALSE)
metrics <- glue("{alignment_dir}/{samples}/outs/metrics_summary.csv")
names(metrics) <- samples
metrics <- metrics[file.exists(metrics)]
metrics
AndrewDHildreth2021_metrics_df <- map_dfr(samples, function(x) {
    print(x)
    df <- read_csv(metrics[[x]])
    df$sample <- x
    return(df)
})
AndrewDHildreth2021_metrics_df$study <- "AndrewDHildreth2021"

# ShuaiHe2020 -------------------------------------------------------------
alignment_dir <- "/mnt/public/ShuaiHe2020/alignments/"
samples <- list.dirs(alignment_dir, full.names = FALSE, recursive = FALSE)
metrics <- glue("{alignment_dir}/{samples}/outs/metrics_summary.csv")
names(metrics) <- samples
metrics <- metrics[file.exists(metrics)]
metrics
ShuaiHe2020_metrics_df <- map_dfr(samples, function(x) {
    print(x)
    df <- read_csv(metrics[[x]])
    df$sample <- x
    return(df)
})
ShuaiHe2020_metrics_df$study <- "ShuaiHe2020"

# MarylineFalquet2023 -----------------------------------------------------
alignment_dir <- "/mnt/public/MarylineFalquet2023/alignments/"
samples <- list.dirs(alignment_dir, full.names = FALSE, recursive = FALSE)
metrics <- glue("{alignment_dir}/{samples}/outs/metrics_summary.csv")
names(metrics) <- samples
metrics <- metrics[file.exists(metrics)]
metrics

MarylineFalquet2023_metrics_df <- map_dfr(samples, function(x) {
    print(x)
    df <- read_csv(metrics[[x]])
    df$sample <- x
    return(df)
})
MarylineFalquet2023_metrics_df$study <- "MarylineFalquet2023"

# Combine all metrics -----------------------------------------------------
all_metrics_df <- bind_rows(
    CDominguezConde2022 = CDominguezConde2022_metrics_df,
    NataliaJaeger2024 = NataliaJaeger2024_metrics_df,
    AndrewDHildreth2021 = AndrewDHildreth2021_metrics_df,
    ShuaiHe2020 = ShuaiHe2020_metrics_df,
    MarylineFalquet2023 = MarylineFalquet2023_metrics_df
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
