library(Shennong)
library(tidyverse)
library(Seurat)


# 20251117 ----------------------------------------------------------------
if (FALSE) {
  seurat_obj <- sn_read("/home/sduan/projects/immune-atlas/data/processed/scrup/harmonized.qs")
  
  metadata <- read.csv("data/processed/scrup/integrated/20251117_metadata.csv",
                       row.names = 1)
  umap_path <- "data/processed/scrup/integrated/20251117_umap.csv"
  umap <- Matrix::as.matrix(read.csv(file = umap_path, 
                                     row.names = 1))
  seurat_obj <- seurat_obj[, rownames(umap)]
  seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(umap, key = "umap")
  seurat_obj <- Seurat::AddMetaData(object = seurat_obj, metadata = metadata %>% 
                                      select(batch, leiden))
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj
  DimPlot(seurat_obj, group.by = "leiden", label = TRUE)
  
  # FeaturePlot(seurat_obj, features = "KLRF1", order = TRUE)
  
  cell_type_markers <- c(
    "PTPRC","CD2", "CD3D",
    "CD4","CD40LG","CD8A","CD8B",
    "KLRF1",
    "AIF1","CD14", 
    "CXCR2","FCGR3B",
    "CD79A","COL1A1","COL1A2",
    "VWF","CDH5", "EPCAM","KRF19",
    "HBB","HBA1"
  )
  # sn_plot_dot(seurat_obj, features = cell_type_markers, group_by = "leiden")
  
  seurat_obj@meta.data <- seurat_obj@meta.data |> 
    mutate(cell_type_level1 = case_when(
      leiden %in% c(23, 25,26) ~ "B cells",
      leiden %in% c(0,2,6) ~ "CD4+ T cells",
      leiden %in% c(5,7,9) ~ "CD8+ T cells",
      leiden %in% c(12) ~ "Endothelial cells",
      leiden %in% c(31,29,28,21,8,1,17,19,20,24,4) ~ "Epithelial cells",
      leiden %in% c(11) ~ "Innate lymphoid cells",
      leiden %in% c(3,13,14,30,15,10) ~ "Myeloid cells",
      leiden %in% c(27,16,22) ~ "Stromal cells",
      leiden %in% c(18) ~ "Erythrocytes",
      TRUE ~ as.character(leiden)
    ))
  seurat_obj$cell_type_level1 <- factor(
    seurat_obj$cell_type_level1, levels = c(
      "CD4+ T cells","CD8+ T cells",
      "B cells", "Innate lymphoid cells",
      "Myeloid cells","Epithelial cells","Stromal cells",
      "Endothelial cells", "Erythrocytes"
    )
  )
  # sn_plot_dim(seurat_obj, group_by = "cell_type_level1",pt_size = 0.5,label = TRUE)
  # table(seurat_obj$cell_type_level1)
  # DimPlot(seurat_obj, group.by = "cell_type_level1", label = TRUE)
  
  sahmi_table <- sn_read("/home/sduan/data/scrup/object/sahmi.table.qs")
  table(colnames(seurat_obj) %in% sahmi_table$barcode)
  bacteria_barcodes <- sahmi_table |> 
    filter(kingdom == "Bacteria") |> pull(barcode)
  eukaryota_barcodes <- sahmi_table |> 
    filter(kingdom == "Eukaryota") |> pull(barcode)
  viruses_barcodes <- sahmi_table |> 
    filter(kingdom == "Viruses") |> pull(barcode)
  
  seurat_obj$bacteria <- if_else(colnames(seurat_obj) %in% bacteria_barcodes, "Detected", "No detected")
  seurat_obj$eukaryota <- if_else(colnames(seurat_obj) %in% eukaryota_barcodes, "Detected", "No detected")
  seurat_obj$viruses <- if_else(colnames(seurat_obj) %in% viruses_barcodes, "Detected", "No detected")
  
  sn_write(seurat_obj, "/home/sduan/projects/immune-atlas/data/processed/scrup/20251117.integrated.qs")
  
  # sn_plot_feature(seurat_obj, features = "MKI67")
  sn_write(seurat_obj, "/home/sduan/projects/immune-atlas/data/processed/scrup/20251117.integrated.qs")
}

# Endothelial cells -------------------------------------------------------
seurat_obj <- sn_read("/home/sduan/projects/immune-atlas/data/processed/scrup/20251117.integrated.qs")
counts <- LayerData(seurat_obj, layer = "counts")
BPCells::write_matrix_dir(counts, dir = "data/processed/scrup/integrated/bpcells/")
counts <- sn_read("data/processed/scrup/integrated/bpcells/")
LayerData(seurat_obj, layer = "counts") <- counts
seurat_obj <- NormalizeData(seurat_obj)
sn_write(seurat_obj, "/home/sduan/projects/immune-atlas/data/processed/scrup/20251117.integrated.normalized.qs")

sub_seurat_obj <- seurat_obj |> 
  subset(cell_type_level1 == "Endothelial cells")

metadata <- read.csv("data/processed/scrup/integrated/20251117_metadata.csv",
                     row.names = 1)
umap_path <- "data/processed/scrup/integrated/20251117_umap.csv"
umap <- Matrix::as.matrix(read.csv(file = umap_path, 
                                   row.names = 1))
seurat_obj <- seurat_obj[, rownames(umap)]
seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(umap, key = "umap")
seurat_obj <- Seurat::AddMetaData(object = seurat_obj, metadata = metadata %>% 
                                    select(batch, leiden))
