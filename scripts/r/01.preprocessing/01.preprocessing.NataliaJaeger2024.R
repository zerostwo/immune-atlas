library(tidyverse)
library(Shennong)
library(glue)
library(Seurat)

outdir <- sn_set_path("data/processed/NataliaJaeger2024/corrected")

dirs <- sn_list_10x_paths("/mnt/public/NataliaJaeger2024/alignments/")
for (sample in names(dirs)) {
    print(sample)
    
    so_path <- glue("{outdir}/seurat_object/{sample}.qs")
    if (!file.exists(so_path)) {
        seurat_obj <- sn_initialize_seurat_object(x = dirs[sample])
        seurat_obj <- sn_filter_cells(seurat_obj, features = c("percent.mt", "nFeature_RNA"))
        seurat_obj <- sn_filter_genes(seurat_obj, min_cells = 1)
        seurat_obj <- sn_remove_ambient_contamination(seurat_obj)
        seurat_obj <- sn_find_doublets(seurat_obj, layer = "decontaminated_counts")
        
        sn_write(seurat_obj, path = so_path)
    } else {
        seurat_obj <- sn_read(so_path)
    }
    
    bpcells_so_path <- glue("{outdir}/bpcells/{sample}.qs")
    if (!file.exists(bpcells_so_path)) {
        bpcells_path <- glue("{outdir}/bpcells/{sample}")
        LayerData(seurat_obj, layer = "counts") |>
            sn_write(path = bpcells_path,
                     to = "bpcells",
                     overwrite = TRUE)
        counts <- sn_read(bpcells_path)
        LayerData(seurat_obj, layer = "counts") <- counts
        
        bpcells_path <- glue("{outdir}/bpcells/{sample}_corrected")
        LayerData(seurat_obj, layer = "counts") |>
            sn_write(path = bpcells_path,
                     to = "bpcells",
                     overwrite = TRUE)
        counts <- sn_read(bpcells_path)
        LayerData(seurat_obj, layer = "decontaminated_counts") <- counts
        
        sn_write(seurat_obj, path = bpcells_so_path)
    }
}

if (FALSE) {
    mat <- sn_read(
        "/mnt/public/NataliaJaeger2024/alignments/Tonsils_PBMC_NK_cells/outs/filtered_feature_bc_matrix.h5"
    )
    counts <- mat$`Gene Expression`
    adt <- mat$`Antibody Capture`
    
    hto_ids <- rownames(adt)[str_detect(rownames(adt), "TotalSeq")]
    
    hto <- adt[rownames(adt) %in% hto_ids, ]
    
    seurat_obj <- sn_initialize_seurat_object(x = counts, species = "human")
    seurat_obj[["ADT"]] <- CreateAssay5Object(counts = adt)
    seurat_obj[["HTO"]] <- CreateAssay5Object(counts = hto)
    keep_cells <- names(colSums(hto))[colSums(hto) > 0]
    seurat_obj <- seurat_obj |> subset(cells = keep_cells)
    # seurat_obj <- seurat_obj |>
    #     subset(percent.mt <= 20) |>
    #     subset(nFeature_RNA > 500) |>
    #     subset(nCount_RNA > 1000)
    seurat_obj <- sn_filter_cells(seurat_obj, features = c("percent.mt", "nFeature_RNA"))
    
    raw <- sn_read(
        "/mnt/public/NataliaJaeger2024/alignments/Tonsils_PBMC_NK_cells/outs/raw_feature_bc_matrix/"
    )
    seurat_obj <- sn_remove_ambient_contamination(seurat_obj, raw = raw$`Gene Expression`)
    
    seurat_obj <- NormalizeData(seurat_obj, layer = "decontaminated_counts")
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mean.var.plot")
    seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
    seurat_obj <- NormalizeData(seurat_obj,
                                assay = "HTO",
                                normalization.method = "CLR")
    
    seurat_obj <- HTODemux(seurat_obj,
                           assay = "HTO",
                           positive.quantile = 0.99)
    
    table(seurat_obj$HTO_classification.global)
    Idents(seurat_obj) <- "HTO_classification.global"
    seurat_obj <- subset(seurat_obj, idents = "Singlet")
    seurat_obj$sample <- seurat_obj$HTO_maxID
    print(table(seurat_obj$sample))
    seurat_obj <- sn_find_doublets(seurat_obj, layer = "decontaminated_counts")
    
    sample <- "Tonsils_PBMC_NK_cells"
    sn_write(seurat_obj, path = glue("{outdir}/seurat_object/{sample}.qs"))
    
    bpcells_path <- glue("{outdir}/bpcells/{sample}")
    LayerData(seurat_obj, layer = "counts") |>
        sn_write(path = bpcells_path,
                 to = "bpcells",
                 overwrite = TRUE)
    counts <- sn_read(bpcells_path)
    LayerData(seurat_obj, layer = "counts") <- counts
    seurat_obj$orig.ident <- seurat_obj$sample
    colnames(seurat_obj) <- paste0(seurat_obj$sample, "_", colnames(seurat_obj))
    sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
    
}
# Clustering --------------------------------------------------------------
if (FALSE) {
    list.files("data/processed/NataliaJaeger2024/corrected/seurat_object/")
    # 1. Blood
    blood1 <- sn_read("data/processed/NataliaJaeger2024/corrected/seurat_object/Blood1.qs")
    blood2 <- sn_read("data/processed/NataliaJaeger2024/corrected/seurat_object/Blood2.qs")
    blood <- merge(blood1,
                   y = blood2,
                   add.cell.ids = c("Blood1", "Blood2"))
    blood <- JoinLayers(blood, layers = c("counts", "decontaminated_counts"))
    table(blood$scDblFinder.class_corrected)
    quantile(blood$nFeature_RNA)
    quantile(blood$percent.mt)
    quantile(blood$nFeature_RNA_corrected)
    blood <- sn_filter_genes(blood, gene_class = "coding")
    blood <- blood |>
        subset(nFeature_RNA_corrected > 200) |>
        subset(scDblFinder.class_corrected == "singlet")
    blood <- sn_run_cluster(blood, batch = "sample", layer = "decontaminated_counts")
    sn_plot_dim(blood, label = TRUE, pt_size = 1.5)
    sn_plot_dim(
        blood,
        label = TRUE,
        pt_size = 1.5,
        group_by = "scDblFinder.class_corrected",
        palette = "OkabeIto"
    )
    sn_plot_dim(
        blood,
        label = TRUE,
        pt_size = 1.5,
        group_by = "sample",
        palette = "OkabeIto"
    )
    sn_plot_dim(
        blood,
        label = TRUE,
        pt_size = 1.5,
        group_by = "Phase",
        palette = "OkabeIto"
    )
    sn_plot_feature(blood, features = "IL32", pt_size = 1.5)
    blood <- sn_find_de(blood, analysis = "markers", layer = "decontaminated_counts")
    sn_plot_dot(
        tonsil,
        features = c(
            "PTPRC",
            "IL7R",
            "CD3D",
            "CD3E",
            "CD3G",
            "CD8A",
            "CD8B",
            "CD4",
            "CD40LG",
            "KLRF1",
            "EOMES",
            "PRDM1",
            "GPR34",
            "ENTPD1",
            "GATA3",
            "PTGDR2",
            "HPGDS",
            "KLRG1",
            "IL32",
            "IL10RA",
            "GRIP1",
            "RORC",
            "KIT",
            "PCDH9"
        ),
        group_by = "cell_type_level2"
    )
    markers <- sn_get_de_result(blood, top_n = 20)
    split(markers$gene, markers$cluster)
    blood@meta.data <- blood@meta.data |>
        mutate(
            cell_type_level2 = case_when(
                seurat_clusters %in% c(3, 4) ~ "ILCPs",
                seurat_clusters %in% c(6, 9, 10, 0, 1) ~ "ILC2s",
                TRUE ~ "Undefined"
            )
        )
    sn_plot_dim(
        blood,
        label = TRUE,
        pt_size = 1.5,
        group_by = "cell_type_level2",
        palette = "OkabeIto"
    )
    sn_write(blood,
             "data/processed/NataliaJaeger2024/corrected/annotation/Blood.qs")
    
    # 2. Tonsil
    tonsil1 <- sn_read("data/processed/NataliaJaeger2024/corrected/seurat_object/Tonsil1.qs")
    tonsil2 <- sn_read("data/processed/NataliaJaeger2024/corrected/seurat_object/Tonsil2.qs")
    tonsil3 <- sn_read("data/processed/NataliaJaeger2024/corrected/seurat_object/Tonsil3.qs")
    tonsil4 <- sn_read("data/processed/NataliaJaeger2024/corrected/seurat_object/Tonsil4.qs")
    tonsil5 <- sn_read("data/processed/NataliaJaeger2024/corrected/seurat_object/Tonsil5.qs")
    tonsil_nk <- sn_read(
        "data/processed/NataliaJaeger2024/corrected/seurat_object/Tonsils_PBMC_NK_cells.qs"
    )
    tonsil_nk <- tonsil_nk |>
        subset(sample != "PBMC101923-TotalSeqB")
    # tonsil <- merge(
    #     tonsil1,
    #     y = c(tonsil2, tonsil3, tonsil4, tonsil5, tonsil_nk),
    #     add.cell.ids = c(
    #         "Tonsil1",
    #         "Tonsil2",
    #         "Tonsil3",
    #         "Tonsil4",
    #         "Tonsil5",
    #         "Tonsils_PBMC_NK_cells"
    #     )
    # )
    tonsil <- merge(tonsil4,
                    y = c(tonsil5),
                    add.cell.ids = c("Tonsil4", "Tonsil5"))
    
    tonsil <- JoinLayers(tonsil, layers = c("counts", "decontaminated_counts"))
    tonsil@assays$ADT <- NULL
    tonsil@assays$HTO <- NULL
    LayerData(tonsil, layer = "data.6") <- NULL
    LayerData(tonsil, layer = "scale.data.6") <- NULL
    
    table(tonsil$scDblFinder.class_corrected)
    quantile(tonsil$nFeature_RNA)
    quantile(tonsil$percent.mt)
    quantile(tonsil$nFeature_RNA_corrected)
    tonsil <- sn_filter_genes(tonsil, gene_class = "coding")
    tonsil <- tonsil |>
        subset(nFeature_RNA_corrected > 200) |>
        subset(scDblFinder.class_corrected == "singlet")
    tonsil <- sn_run_cluster(tonsil, batch = "sample", layer = "decontaminated_counts")
    sn_plot_dim(tonsil, label = TRUE, pt_size = 1.5)
    sn_plot_dim(
        tonsil,
        label = TRUE,
        pt_size = 1.5,
        group_by = "scDblFinder.class_corrected",
        palette = "OkabeIto"
    )
    sn_plot_dim(
        tonsil,
        label = TRUE,
        pt_size = 1.5,
        split_by = "sample",
        palette = "OkabeIto"
    )
    sn_plot_dim(
        tonsil,
        label = TRUE,
        pt_size = 1.5,
        group_by = "Phase",
        palette = "OkabeIto"
    )
    sn_plot_feature(tonsil, features = "LYZ", pt_size = 1.5)
    tonsil <- sn_find_de(tonsil, analysis = "markers", layer = "decontaminated_counts")
    markers <- sn_get_de_result(tonsil, top_n = 20)
    split(markers$gene, markers$cluster)
    
    tonsil@meta.data <- tonsil@meta.data |>
        mutate(
            cell_type_level2 = case_when(
                seurat_clusters %in% c(11, 6, 7) ~ "ILC2s",
                seurat_clusters %in% c(0, 1, 3, 9, 2) ~ "ILC3s",
                # seurat_clusters %in% c(4) ~ "NK cells",
                TRUE ~ "Undefined"
            )
        )
    sn_plot_dim(
        tonsil,
        label = TRUE,
        pt_size = 1.5,
        group_by = "cell_type_level2",
        palette = "OkabeIto"
    )
    sn_write(tonsil,
             "data/processed/NataliaJaeger2024/corrected/annotation/Tonsil.qs")
    
    # 3. Lung
    
}
