library(tidyverse)
library(Shennong)
library(glue)
library(Seurat)

outdir <- sn_set_path("data/processed/TabulaSapiens/corrected")

dirs <- sn_list_10x_paths("/home/sduan/data/TabulaSapiens/data/alignments/10x")
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

# library(tidyverse)
# library(Shennong)
# library(glue)
# library(Seurat)
#
# outdir <- sn_set_path("data/processed/20260328/TabulaSapiens")
#
# samples <- list.files("/home/sduan/data/TabulaSapiens/data/alignments/10x")
# filtered_data_dirs <- glue("/home/sduan/data/TabulaSapiens/data/alignments/10x/{samples}/outs/filtered_feature_bc_matrix/")
# names(filtered_data_dirs) <- samples
# filtered_data_dirs <- filtered_data_dirs[dir.exists(filtered_data_dirs)]
# length(filtered_data_dirs)
#
# for (sample in names(filtered_data_dirs)) {
#   # sample <- samples[1]
#   print(sample)
#   so_path <- glue("{outdir}/seurat_object/{sample}.qs")
#   if (!file.exists(so_path)) {
#     seurat_obj <- sn_initialize_seurat_object(
#       x = filtered_data_dirs[sample], species = "human",
#       min_features = 100
#     )
#     raw_counts <- sn_read(str_replace(filtered_data_dirs[sample], "filtered", "raw"))
#     seurat_obj <- sn_remove_ambient_contamination(
#       seurat_obj, raw = raw_counts
#     )
#     seurat_obj <- sn_filter_cells(seurat_obj, features = c("nFeature_RNA","percent.mt"))
#     seurat_obj <- sn_find_doublets(
#       seurat_obj, layer = "decontaminated_counts"
#     )
#     print(table(seurat_obj$scDblFinder.class_corrected))
#     sn_write(seurat_obj, path = so_path)
#
#     # seurat_obj <- seurat_obj |>
#     #   subset(percent.mt <= 20)
#     # seurat_obj <- sn_find_doublets(
#     #   seurat_obj, dbr_sd = 1
#     # )
#     # print(table(seurat_obj$scDblFinder.class))
#     # sn_write(seurat_obj, path = so_path)
#
#     # bpcells_path <- glue("{outdir}/bpcells/{sample}")
#     # LayerData(seurat_obj,layer = "counts") |>
#     #   sn_write(path = bpcells_path, to = "bpcells",overwrite = TRUE)
#     # counts <- sn_read(bpcells_path)
#     # LayerData(seurat_obj, layer = "counts") <- counts
#     # sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
#   }
# }
