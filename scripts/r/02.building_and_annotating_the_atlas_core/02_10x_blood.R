library(tidyverse)
library(Shennong)
library(Seurat)

# Load data ---------------------------------------------------------------
if (FALSE) {
  datasets <- c(
    "pbmc1k","pbmc3k","pbmc4k", "pbmc8k"
  )
  
  seurat_obj_lists <- map(datasets, function(d) {
    # d <- "pbmc1k"
    print(d)
    seurat_obj <- sn_load_data(dataset = d)
    raw_counts <- sn_load_data(dataset = d, matrix_type = "raw")
    seurat_obj <- sn_filter_cells(seurat_obj,
                                  features = c("percent.mt","nFeature_RNA"))
    seurat_obj <- seurat_obj %>% 
      subset(nFeature_RNA > 200)
    seurat_obj <- sn_remove_ambient_contamination(
      seurat_obj, raw = raw_counts, method = "soupx"
    )
    seurat_obj <- sn_find_doublets(seurat_obj, layer = "decontaminated_counts")
    seurat_obj
  })
  seurat_obj <- merge(
    x  = seurat_obj_lists[[1]],
    y = seurat_obj_lists[-1],
    add.cell.ids = datasets
  )
  seurat_obj <- JoinLayers(seurat_obj,
                           layers = c("counts", "decontaminated_counts"))
  sn_write(seurat_obj, "data/processed/pbmc.qs")
  
}
# Merged ------------------------------------------------------------------
pbmc <- sn_read("data/processed/pbmc.qs")

seurat_obj_lists <- list.files("data/processed/MarylineFalquet2023/corrected_soupx/seurat_object/",
           full.names = TRUE) %>% 
  lapply(sn_read)
MarylineFalquet2023 <- merge(x = seurat_obj_lists[[1]],
                             y = seurat_obj_lists[-1],
                             add.cell.ids = c("HD1","HD2","H3"))
MarylineFalquet2023 <- JoinLayers(
  MarylineFalquet2023, layers = c("counts", "decontaminated_counts")
)

seurat_obj_lists <- list.files("data/processed/NataliaJaeger2024/corrected_soupx/seurat_object/",
                               full.names = TRUE)[c(1:2,17)] %>% 
  lapply(sn_read)
NataliaJaeger2024 <- merge(x = seurat_obj_lists[[1]],
                             y = seurat_obj_lists[-1],
                             add.cell.ids = c("Blood1","Blood2","Tonsils_PBMC_NK_cells"))
NataliaJaeger2024 <- JoinLayers(
  NataliaJaeger2024, layers = c("counts", "decontaminated_counts")
)
NataliaJaeger2024@assays$HTO <- NULL
NataliaJaeger2024@assays$ADT <- NULL
NataliaJaeger2024 <- NataliaJaeger2024 %>% 
  subset(!(sample %in% c("TS1101923-TotalSeqB", "TS2101923-TotalSeqB")))
LayerData(NataliaJaeger2024, layer = "data.3") <- NULL
LayerData(NataliaJaeger2024, layer = "scale.data.3") <- NULL
NataliaJaeger2024@meta.data <- NataliaJaeger2024@meta.data %>% 
  select(-c(hash.ID, HTO_classification, HTO_classification.global,
            HTO_maxID, HTO_margin, nCount_HTO, nFeature_HTO,
            HTO_secondID))

seurat_obj_lists <- list.files("data/processed/ZhenlongLi2024/corrected_soupx/seurat_object/",
                               full.names = TRUE)[1:2] %>% 
  lapply(sn_read)
ZhenlongLi2024 <- merge(x = seurat_obj_lists[[1]],
                           y = seurat_obj_lists[-1],
                           add.cell.ids = c("H","H2"))
ZhenlongLi2024 <- JoinLayers(
  ZhenlongLi2024, layers = c("counts", "decontaminated_counts")
)

DanielPCaron2025 <- sn_read("~/FACS.qs2")
DanielPCaron2025$cell_type <- DanielPCaron2025$sample
DanielPCaron2025$sample <- "FACS"
DanielPCaron2025@assays$HTO <- NULL
DanielPCaron2025@meta.data <- DanielPCaron2025@meta.data %>% 
  select(-c(hash.ID, HTO_classification, HTO_classification.global,
           HTO_maxID, HTO_margin, nCount_HTO, nFeature_HTO,
           HTO_secondID))

metadata <- sn_read("data/raw/reference/reference.raw.qs")
metadata <- metadata@meta.data
metadata <- metadata %>% 
  filter(tissue_level1=="Blood") %>% 
  filter(study == "CDominguezConde2022") 
seurat_obj_lists <- list(
  sn_read("data/processed/CDominguezConde2022/corrected_soupx/seurat_object/A35-BLD-0-SC-1.qs"),
  sn_read("data/processed/CDominguezConde2022/corrected_soupx/seurat_object/A36-BLD-0-SC-1.qs")
)
CDominguezConde2022 <- merge(
  x = seurat_obj_lists[[1]],
  y = seurat_obj_lists[-1],
  add.cell.ids = c("A35-BLD-0-SC-1", "A36-BLD-0-SC-1")
)
CDominguezConde2022 <- JoinLayers(
  CDominguezConde2022, layers = c("counts", "decontaminated_counts")
)

merged <- merge(
  x = pbmc,
  y = list(MarylineFalquet2023, NataliaJaeger2024,
           DanielPCaron2025, ZhenlongLi2024, CDominguezConde2022)
)

merged <- JoinLayers(
  merged, layers = c("counts", "decontaminated_counts")
)

LayerData(merged, layer = "data.4") <- NULL
LayerData(merged, layer = "scale.data.4") <- NULL

merged

merged <- merged %>% 
  subset(nFeature_RNA > 200)
merged <- sn_filter_genes(merged, min_cells = 3)

sn_write(merged, "data/processed/pbmc.qs")

# Integrated --------------------------------------------------------------
merged <- sn_read("data/processed/pbmc.qs")

filtered_merged <- merged %>% 
  subset(scDblFinder.class_corrected == "singlet") %>% 
  subset(nFeature_RNA > 500) %>% 
  subset(nCount_RNA > 1000)

merged <- merged %>%
  sn_run_cluster(
    integration_method = "coralysis",
    batch_by = "sample",
    layer = "decontaminated_counts",
    integration_control = list(
      icp_args = list(
        threads = 2
      )
    )
  )

# seurat_obj <- FindClusters(
#   seurat_obj, resolution = 1.9, algorithm = 4, random.seed = 717
# )
sn_write(merged, "data/processed/pbmc.integrated.qs2")
DimPlot(merged, label = TRUE) + 
  NoLegend() + NoAxes()
DimPlot(merged, label = TRUE) + 
  NoLegend() + NoAxes()
FeaturePlot(merged, features = "", order = TRUE)
FeaturePlot(merged, features = "CCR6", order = TRUE)

table(merged$sample)
merged %>% 
  subset(sample == "pbmc8k") %>% 
  # DimPlot(group.by = "cell_type",label = TRUE)
  FeaturePlot(features = "MKI67", order = TRUE)

markers <- FindAllMarkers(
  merged,
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.5
)

top_markers <- markers %>% 
  filter(!str_detect(gene, "^ENSG")) %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 20)
split(top_markers$gene, top_markers$cluster)
