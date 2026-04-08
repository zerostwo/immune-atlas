# library(tidyverse)
# library(Shennong)
# library(glue)
# library(Seurat)
#
# outdir <- sn_set_path("data/processed/EMadissoon2020")
#
# samples <- list.files("/mnt/public/EMadissoon2020/alignments/", pattern = "^HCA")
# filtered_data_dirs <- glue("/mnt/public/EMadissoon2020/alignments/{samples}/outs/filtered_feature_bc_matrix")
# names(filtered_data_dirs) <- samples
# filtered_data_dirs <- filtered_data_dirs[dir.exists(filtered_data_dirs)]
# for (sample in names(filtered_data_dirs)) {
#     # sample <- samples[1]
#     print(sample)
#
#     seurat_obj <- sn_initialize_seurat_object(
#         x = filtered_data_dirs[sample], species = "human",
#         min_features = 200
#     )
#     seurat_obj <- seurat_obj |>
#         subset(percent.mt <= 20)
#     seurat_obj <- sn_find_doublets(
#         seurat_obj, dbr_sd = 1
#     )
#     print(table(seurat_obj$scDblFinder.class))
#     sn_write(seurat_obj, path = glue("{outdir}/seurat_object/{sample}.qs"))
#
#     bpcells_path <- glue("{outdir}/bpcells/{sample}")
#     LayerData(seurat_obj,layer = "counts") |>
#         sn_write(path = bpcells_path, to = "bpcells",overwrite = TRUE)
#     counts <- sn_read(bpcells_path)
#     LayerData(seurat_obj, layer = "counts") <- counts
#     sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
# }

library(tidyverse)
library(Shennong)
library(glue)
library(Seurat)

outdir <- sn_set_path("data/processed/EMadissoon2020/corrected")

dirs <- sn_list_10x_paths("/mnt/public/EMadissoon2020/alignments")
samples <- names(dirs)
samples <- samples[!(samples %in% c("HCATisStabAug177276393"))]
for (sample in samples) {
    print(sample)
    so_path <- glue("{outdir}/seurat_object/{sample}.qs")
    if (!file.exists(so_path)) {
        seurat_obj <- sn_initialize_seurat_object(x = dirs[sample])
        seurat_obj <- sn_filter_cells(seurat_obj, features = c("percent.mt", "nFeature_RNA"))
        seurat_obj <- sn_filter_genes(seurat_obj, min_cells = 1)
        seurat_obj <- sn_remove_ambient_contamination(seurat_obj)
        seurat_obj <- sn_find_doublets(seurat_obj, layer = "decontaminated_counts")
        
        sn_write(seurat_obj,
                 path = so_path)
        
        bpcells_path <- glue("{outdir}/bpcells/{sample}")
        LayerData(seurat_obj, layer = "counts") |>
            sn_write(path = bpcells_path,
                     to = "bpcells",
                     overwrite = TRUE)
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

# sample_metadata <- sn_read("/mnt/public/EMadissoon2020/13059_2019_1906_MOESM2_ESM.xlsx",
#                     sheet = 2)
# donor_metadata <- sn_read("/mnt/public/EMadissoon2020/13059_2019_1906_MOESM2_ESM.xlsx",
#                            sheet = 1)
# metadata <- sample_metadata |>
#     left_join(donor_metadata, by = c("patient" = "ID")) |>
#     select(-c(`DBD/DCD`, `Tissues collected`))
