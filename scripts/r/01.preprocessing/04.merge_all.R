library(tidyverse)
library(Seurat)
library(glue)
library(Shennong)

ZhenlongLi2024 <- sn_read("data/processed/ZhenlongLi2024/merged.qs")
ZhenlongLi2024$study <- "ZhenlongLi2024"

StevenBWells2025 <- sn_read("data/processed/StevenBWells2025/merged_corrected.qs")
StevenBWells2025$study <- "StevenBWells2025"
StevenBWells2025 <- JoinLayers(StevenBWells2025, layers = c(
  "counts","decontaminated_counts"
))
StevenBWells2025@meta.data <- StevenBWells2025@meta.data |> 
  select(-c(HTO_maxID, HTO_secondID, HTO_margin, HTO_classification, HTO_classification.global,
         hash.ID, nFeature_HTO, nCount_HTO))

TabulaSapiens <- sn_read("data/processed/TabulaSapiens/merged_corrected.qs")
TabulaSapiens$study <- "TabulaSapiens"
TabulaSapiens <- JoinLayers(TabulaSapiens, layers = c(
  "counts","decontaminated_counts"
))

# scrup <- sn_read("data/processed/scrup/merged.qs")
# scrup$study <- "scrup"

merged <- merge(
  x = ZhenlongLi2024,
  y = list(StevenBWells2025, TabulaSapiens)
)

merged <- JoinLayers(merged, layers = c(
  "counts","decontaminated_counts"
))
merged <- JoinLayers(merged, assay = "ADT")
merged
table(merged$study)
sn_write(
  merged,
  path = "data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens.merged_corrected.qs"
)

counts <- LayerData(merged,layer = "counts", assay = "RNA")
BPCells::write_matrix_dir(counts, "data/processed/bpcells_corrected/counts")

decontaminated_counts <- LayerData(merged,layer = "decontaminated_counts", assay = "RNA")
BPCells::write_matrix_dir(counts, "data/processed/bpcells_corrected/decontaminated_counts")

adt_counts <- LayerData(merged,layer = "counts", assay = "ADT")
BPCells::write_matrix_dir(adt_counts, "data/processed/bpcells_corrected/adt_counts")

sample_metadata <- merged@meta.data |> 
  select(study, orig.ident) |> 
  distinct()
sn_write(
  sample_metadata,# you may not know it
  path = "data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens.sample_metadata.qs"
)

