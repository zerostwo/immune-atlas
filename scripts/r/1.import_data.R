library(tidyverse)
library(Seurat)
library(glue)

# (Domínguez Conde C et al., Science, 2022) -------------------------------
samples <- list.files("/mnt/rstudio/CDominguezConde2022/alignments/count")

data_dirs <- glue(
    "/mnt/rstudio/CDominguezConde2022/alignments/count/{samples}/outs/filtered_feature_bc_matrix/"
)
names(data_dirs) <- samples

seurat_obj_lists <- map(names(data_dirs), function(x) {
    print(x)
    counts <- Read10X(data.dir = data_dirs[[x]])
    seurat_obj <- CreateSeuratObject(counts = counts, project = x)
    colnames(seurat_obj) <- paste(x, colnames(seurat_obj), sep = "_")
    return(seurat_obj)
})

seurat_obj <- merge(
    seurat_obj_lists[[1]],
    y = seurat_obj_lists[-1]
)
seurat_obj <- JoinLayers(seurat_obj)

# Mitochondrial genes percentage
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-",
    col.name = "percent.mt"
)
# Ribosomal genes percentage
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^RPS|^RPL",
    col.name = "percent.ribo"
)
# percent.hb
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^HB[^(P)]",
    col.name = "percent.hb"
)
qs::qsave(seurat_obj, "data/raw/object/seurat/CDominguezConde2022.qs")


# (Jaeger N et al., Nature Immunology, 2024) ------------------------------
samples <- list.files("/mnt/public/NataliaJaeger2024/alignments/")
data_dirs <- glue(
    "/mnt/public/NataliaJaeger2024/alignments/{samples}/outs/filtered_feature_bc_matrix/"
)
names(data_dirs) <- str_replace_all(samples, "_", "-")
data_dirs <- data_dirs[dir.exists(data_dirs)]
data_dirs <- data_dirs[-length(data_dirs)]

seurat_obj_lists <- map(names(data_dirs), function(x) {
    print(x)
    counts <- Read10X(data.dir = data_dirs[[x]])
    seurat_obj <- CreateSeuratObject(counts = counts, project = x)
    colnames(seurat_obj) <- paste(x, colnames(seurat_obj), sep = "_")
    return(seurat_obj)
})

seurat_obj <- merge(
    seurat_obj_lists[[1]],
    y = seurat_obj_lists[-1]
)
seurat_obj <- JoinLayers(seurat_obj)

# Mitochondrial genes percentage
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-",
    col.name = "percent.mt"
)
# Ribosomal genes percentage
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^RPS|^RPL",
    col.name = "percent.ribo"
)
# percent.hb
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^HB[^(P)]",
    col.name = "percent.hb"
)
qs::qsave(seurat_obj, "data/raw/object/seurat/NataliaJaeger2024.qs")


# (Hildreth A D et al., Nature Immunology, 2021) --------------------------
samples <- list.files("/mnt/public/AndrewDHildreth2021/alignments")
data_dirs <- glue(
    "/mnt/public/AndrewDHildreth2021/alignments/{samples}/outs/filtered_feature_bc_matrix/"
)
data_dirs <- data_dirs[dir.exists(data_dirs)]

seurat_obj_lists <- map(names(data_dirs), function(x) {
    print(x)
    counts <- Read10X(data.dir = data_dirs[[x]])
    seurat_obj <- CreateSeuratObject(counts = counts, project = x)
    colnames(seurat_obj) <- paste(x, colnames(seurat_obj), sep = "_")
    return(seurat_obj)
})
seurat_obj <- merge(
    seurat_obj_lists[[1]],
    y = seurat_obj_lists[-1]
)
seurat_obj <- JoinLayers(seurat_obj)
# Mitochondrial genes percentage
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-",
    col.name = "percent.mt"
)
# Ribosomal genes percentage
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^RPS|^RPL",
    col.name = "percent.ribo"
)
# percent.hb
seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^HB[^(P)]",
    col.name = "percent.hb"
)
qs::qsave(seurat_obj, "data/raw/object/seurat/AndrewDHildreth2021.qs")

# Merge -------------------------------------------------------------------
s1 <- qs::qread("data/raw/object/seurat/CDominguezConde2022.qs")
s2 <- qs::qread("data/raw/object/seurat/NataliaJaeger2024.qs")
s3 <- qs::qread("data/raw/object/seurat/AndrewDHildreth2021.qs")
seurat_obj <- merge(s1, y = c(s2, s3))
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj
qs::qsave(seurat_obj, "data/raw/object/seurat/merged.qs")

seurat_obj <- qs::qread("data/raw/object/seurat/merged.qs")

# Cells identified as low-quality or germ line in the original studies were excluded, and only cells meeting the following criteria were retained: 500–8,000 genes, 1,000–100,000 gene counts, and less than 20% mitochondrial gene counts.
filtered_seurat_obj <- seurat_obj |> 
    subset(nFeature_RNA > 500 & nFeature_RNA < 8000 & nCount_RNA > 1000 & nCount_RNA < 100000 & percent.mt < 20)
filtered_seurat_obj
# We then excluded samples with fewer than 50 high-quality cells.
keep_samples <- names(table(filtered_seurat_obj$orig.ident)[table(filtered_seurat_obj$orig.ident) >= 50])
filtered_seurat_obj <- subset(filtered_seurat_obj, orig.ident %in% keep_samples)

keep_genes <- rowSums(LayerData(filtered_seurat_obj, layer = "counts")) > 10
table(keep_genes)
filtered_seurat_obj <- filtered_seurat_obj[keep_genes, ]

# Run pca
filtered_seurat_obj <- filtered_seurat_obj |> 
    NormalizeData()

# Find variable features b
    ScaleData(vars.to.regress = c("percent.mt","nFeature_RNA"))

filtered_seurat_obj <- RunPCA(filtered_seurat_obj)
ElbowPlot(filtered_seurat_obj, ndims = 50)
DimPlot(filtered_seurat_obj, group.by = "orig.ident") +
    NoLegend()

filtered_seurat_obj <- RunUMAP(filtered_seurat_obj, dims = 1:30)
filtered_seurat_obj <- FindNeighbors(filtered_seurat_obj, dims = 1:30)
filtered_seurat_obj <- FindClusters(filtered_seurat_obj, resolution = 0.5)

DimPlot(filtered_seurat_obj, reduction = "umap", label = TRUE) +
    NoLegend()

# Mast cells
FeaturePlot(filtered_seurat_obj, c("CD79A"))
DimPlot()
