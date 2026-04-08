library(tidyverse)
library(Shennong)
library(glue)
library(Seurat)

outdir <- sn_set_path("data/processed/DanielPCaron2025/corrected")

dirs <- sn_list_10x_paths("/home/sduan/data/DanielPCaron2025/alignments/count/")

samples <- names(dirs)
for (sample in samples) {
    # sample <- samples[1]
    print(sample)
    so_path <- glue("{outdir}/seurat_object/{sample}.qs")
    if (!file.exists(so_path)) {
        ff <- glue("{dirs[sample]}/filtered_feature_bc_matrix")
        mat <- sn_read(ff)
        counts <- mat$`Gene Expression`
        adt <- mat$`Antibody Capture`
        
        adt_ids <- rownames(adt)[str_detect(rownames(adt), "TotalSeq")]
        hto_ids <- rownames(adt)[!(rownames(adt) %in% adt_ids)]
        
        hto <- adt[rownames(adt) %in% hto_ids,]
        adt <- adt[rownames(adt) %in% adt_ids,]
        
        seurat_obj <- sn_initialize_seurat_object(
            x = counts, species = "human"
        )
        seurat_obj[["ADT"]] <- CreateAssay5Object(counts = adt)
        
        if (length(hto_ids) > 1) {
            seurat_obj[["HTO"]] <- CreateAssay5Object(counts = hto)
        }
        
        seurat_obj <- sn_filter_cells(seurat_obj, features = c("nFeature_RNA","percent.mt"))
        
        raw_counts <- sn_read(str_replace(ff, "filtered", "raw"))
        raw_counts <- raw_counts$`Gene Expression`
        seurat_obj <- sn_remove_ambient_contamination(
            seurat_obj, raw = raw_counts
        )
        seurat_obj <- sn_find_doublets(
            seurat_obj, layer = "decontaminated_counts"
        )
        
        seurat_obj <- NormalizeData(seurat_obj, layer = "decontaminated_counts")
        seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mean.var.plot")
        seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
        
        if (length(hto_ids) > 1) {
            seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")
            
            seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.99)
            
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
        
        LayerData(seurat_obj, "data") <- NULL
        LayerData(seurat_obj, "scale.data") <- NULL
        seurat_obj@assays$HTO <- NULL
        
        bpcells_path <- glue("{outdir}/bpcells/{sample}")
        LayerData(seurat_obj,layer = "counts") |> 
            sn_write(path = bpcells_path, to = "bpcells",overwrite = TRUE)
        counts <- sn_read(bpcells_path)
        LayerData(seurat_obj, layer = "counts") <- counts
        
        bpcells_path <- glue("{outdir}/bpcells/{sample}_corrected")
        LayerData(seurat_obj,layer = "counts") |> 
            sn_write(path = bpcells_path, to = "bpcells",overwrite = TRUE)
        counts <- sn_read(bpcells_path)
        LayerData(seurat_obj, layer = "decontaminated_counts") <- counts
        
        sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
    }
}
# 
# for (sample in names(filtered_data_dirs)) {
#     # sample <- samples[1]
#     print(sample)
#     so_path <- glue("{outdir}/seurat_object/{sample}.qs")
#     if (!file.exists(so_path)) {
#         mat <- sn_read(filtered_data_dirs[sample])
#         counts <- mat$`Gene Expression`
#         adt <- mat$`Antibody Capture`
#         
#         adt_ids <- rownames(adt)[str_detect(rownames(adt), "TotalSeq")]
#         hto_ids <- rownames(adt)[!(rownames(adt) %in% adt_ids)]
#         
#         hto <- adt[rownames(adt) %in% hto_ids,]
#         adt <- adt[rownames(adt) %in% adt_ids,]
#         
#         seurat_obj <- sn_initialize_seurat_object(
#             x = counts, species = "human"
#         )
#         seurat_obj[["ADT"]] <- CreateAssay5Object(counts = adt)
#         
#         if (length(hto_ids) > 1) {
#             seurat_obj[["HTO"]] <- CreateAssay5Object(counts = hto)
#         }
#         
#         seurat_obj <- seurat_obj |>
#             subset(percent.mt <= 20) |>
#             subset(nFeature_RNA > 500) |>
#             subset(nCount_RNA > 1000)
#         # seurat_obj <- seurat_obj |> 
#         #     subset(percent.mt <= 20) |> 
#         #     subset(nFeature_RNA > 200)
#         
#         seurat_obj <- NormalizeData(seurat_obj)
#         seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mean.var.plot")
#         seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
#         
#         if (length(hto_ids) > 1) {
#             seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")
#             
#             seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.99)
#             
#             table(seurat_obj$HTO_classification.global)
#             Idents(seurat_obj) <- "HTO_classification.global"
#             seurat_obj <- subset(seurat_obj, idents = "Singlet")
#             seurat_obj$sample <- seurat_obj$HTO_maxID
#             print(table(seurat_obj$sample))
#             
#             seurat_obj_lists <- map(names(table(seurat_obj$sample)), function(x){
#                 print(x)
#                 temp <- seurat_obj |> subset(sample == x)
#                 temp <- sn_find_doublets(
#                     temp, dbr_sd = 1)    
#                 temp
#             })
#             seurat_obj <- merge(seurat_obj_lists[[1]], seurat_obj_lists[-1])
#             seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
#             seurat_obj <- JoinLayers(seurat_obj, assay = "HTO")
#             seurat_obj <- JoinLayers(seurat_obj, assay = "ADT")
#         } else {
#             seurat_obj$sample <- hto_ids
#             seurat_obj <- sn_find_doublets(
#                 seurat_obj, dbr_sd = 1
#             )
#         }
#         print(table(seurat_obj$scDblFinder.class))
#         sn_write(seurat_obj, path = so_path)
#         
#         LayerData(seurat_obj, "data") <- NULL
#         LayerData(seurat_obj, "scale.data") <- NULL
#         seurat_obj@assays$HTO <- NULL
#         
#         bpcells_path <- glue("{outdir}/bpcells/{sample}")
#         LayerData(seurat_obj,layer = "counts") |> 
#             sn_write(path = bpcells_path, to = "bpcells",overwrite = TRUE)
#         counts <- sn_read(bpcells_path)
#         LayerData(seurat_obj, layer = "counts") <- counts
#         sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
#     }
# }
