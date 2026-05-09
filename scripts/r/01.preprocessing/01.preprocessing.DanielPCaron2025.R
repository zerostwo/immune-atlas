library(tidyverse)
library(Shennong)
library(glue)
library(Seurat)

outdir <- sn_set_path("data/processed/DanielPCaron2025/corrected_soupx")

dirs <- sn_list_10x_paths("/home/sduan/data/DanielPCaron2025",
                          path_type = "filtered_h5")
raw <- sn_list_10x_paths("/home/sduan/data/DanielPCaron2025",
                         path_type = "raw_h5")
samples <- names(dirs)
# samples <- list.files(file.path(outdir, "seurat_object"),
#                       pattern = "qs2") %>% str_remove_all(".qs2")
for (sample in samples) {
    print(sample)
    so_path <- glue("{outdir}/seurat_object/{sample}.qs2")
    if (!file.exists(so_path)) {
        mat <- sn_read(dirs[sample])
        counts <- mat$`Gene Expression`
        adt <- mat$`Antibody Capture`

        adt_ids <- rownames(adt)[str_detect(rownames(adt), "TotalSeq")]
        hto_ids <- rownames(adt)[!(rownames(adt) %in% adt_ids)]

        hto <- adt[rownames(adt) %in% hto_ids, ]
        adt <- adt[rownames(adt) %in% adt_ids, ]

        seurat_obj <- sn_initialize_seurat_object(x = counts, species = "human")
        seurat_obj[["ADT"]] <- CreateAssay5Object(counts = adt)

        if (length(hto_ids) > 1) {
            seurat_obj[["HTO"]] <- CreateAssay5Object(counts = hto)
        }

        seurat_obj <- sn_filter_cells(seurat_obj, features = c("nFeature_RNA", "percent.mt"))

        raw_counts <- sn_read(raw[sample])
        raw_counts <- raw_counts$`Gene Expression`
        seurat_obj <- sn_remove_ambient_contamination(seurat_obj,
                                                      method = "soupx",
                                                      raw = raw_counts)
        seurat_obj <- sn_find_doublets(seurat_obj, layer = "decontaminated_counts")

        seurat_obj <- NormalizeData(seurat_obj, layer = "decontaminated_counts")
        seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mean.var.plot")
        seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

        if (length(hto_ids) > 1) {
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
        } else {
            seurat_obj$sample <- hto_ids
        }

        print(table(seurat_obj$scDblFinder.class_corrected))
        sn_write(seurat_obj, path = so_path)
    } else {
        seurat_obj <- sn_read(so_path)
    }

    LayerData(seurat_obj, "data") <- NULL
    LayerData(seurat_obj, "scale.data") <- NULL
    seurat_obj@assays$HTO <- NULL

    bpcells_path <- glue("{outdir}/bpcells/{sample}")
    LayerData(seurat_obj, layer = "counts") |>
        sn_write(path = bpcells_path,
                 to = "bpcells",
                 overwrite = TRUE)
    counts <- sn_read(bpcells_path)
    LayerData(seurat_obj, layer = "counts") <- counts

    bpcells_path <- glue("{outdir}/bpcells/{sample}_corrected")
    LayerData(seurat_obj, layer = "decontaminated_counts") |>
        sn_write(path = bpcells_path,
                 to = "bpcells",
                 overwrite = TRUE)
    counts <- sn_read(bpcells_path)
    LayerData(seurat_obj, layer = "decontaminated_counts") <- counts

    sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs2"))
}
