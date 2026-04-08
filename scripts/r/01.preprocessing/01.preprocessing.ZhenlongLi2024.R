library(tidyverse)
library(Shennong)
library(glue)
library(Seurat)

outdir <- sn_set_path("data/processed/ZhenlongLi2024/corrected")

dirs <- sn_list_10x_paths("/home/sduan/data/ZhenlongLi2024/alignments/")
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
    
    sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
  }
}

files <- list.files(glue("{outdir}/bpcells"),
                    pattern = "qs",
                    full.names = TRUE)
ZhenlongLi2024 <- map(files, sn_read)
ZhenlongLi2024 <- merge(x = ZhenlongLi2024[[1]],
                        y = ZhenlongLi2024[-1],
                        add.cell.ids = str_remove_all(basename(files), ".qs"))

ZhenlongLi2024 <- JoinLayers(ZhenlongLi2024, layers = c("counts", "decontaminated_counts"))
ZhenlongLi2024$study <- "ZhenlongLi2024"

sn_write(
  ZhenlongLi2024,
  "/home/sduan/projects/immune-atlas/data/processed/ZhenlongLi2024/merged.qs"
)

if (FALSE) {
  # Health
  blood1 <- sn_read("data/processed/ZhenlongLi2024/corrected/seurat_object/H.qs")
  blood2 <- sn_read("data/processed/ZhenlongLi2024/corrected/seurat_object/H2.qs")
  blood <- merge(blood1, y = blood2, add.cell.ids = c("H", "H2"))
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
  sn_plot_feature(blood, features = "GZMB", pt_size = 1.5)
  blood <- sn_find_de(blood, analysis = "markers", layer = "decontaminated_counts")
  sn_plot_dot(
    blood,
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
    # group_by = "cell_type_level2"
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
  
  # outdir <- sn_set_path("data/processed/20260328/ZhenlongLi2024")
  #
  # samples <- list.files("/home/sduan/data/ZhenlongLi2024/alignments")
  # filtered_data_dirs <- glue("/home/sduan/data/ZhenlongLi2024/alignments/{samples}/outs/filtered_feature_bc_matrix/")
  # names(filtered_data_dirs) <- samples
  # filtered_data_dirs <- filtered_data_dirs[dir.exists(filtered_data_dirs)]
  # length(filtered_data_dirs)
  
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
  #     # bpcells_path <- glue("{outdir}/bpcells/{sample}")
  #     # LayerData(seurat_obj,layer = "counts") |>
  #     #   sn_write(path = bpcells_path, to = "bpcells",overwrite = TRUE)
  #     # counts <- sn_read(bpcells_path)
  #     # LayerData(seurat_obj, layer = "counts") <- counts
  #     # sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
  #   }
  # }
  
  # files <- list.files(
  #   glue("{outdir}/bpcells"),
  #   pattern = "qs",
  #   full.names = TRUE
  # )
  # ZhenlongLi2024 <- map(files, sn_read)
  # ZhenlongLi2024 <- merge(
  #   x = ZhenlongLi2024[[1]],
  #   y = ZhenlongLi2024[-1]
  # )
  # ZhenlongLi2024 <- JoinLayers(ZhenlongLi2024)
  # ZhenlongLi2024$study <- "ZhenlongLi2024"
  #
  # sn_write(ZhenlongLi2024, "/home/sduan/projects/immune-atlas/data/processed/ZhenlongLi2024/merged.qs")
}