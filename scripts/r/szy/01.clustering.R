library(tidyverse)
library(Seurat)
library(glue)
library(Shennong)

# seuart_obj <- sn_read("data/processed/scrup/20251117.integrated.qs")
# counts <- sn_read("data/processed/scrup/bpcells/")
# LayerData(seuart_obj, layer = "counts") <- counts
# seuart_obj <- NormalizeData(seuart_obj)
# Idents(seuart_obj) <- "cell_type_level1"
# sn_write(seuart_obj, "data/processed/scrup/20251117.integrated.normalized.qs")

seuart_obj <- sn_read("data/processed/scrup/20251117.integrated.normalized.qs")
# DimPlot(seuart_obj)



