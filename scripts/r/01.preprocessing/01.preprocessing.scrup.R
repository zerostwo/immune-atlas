if (TRUE) {
  library(tidyverse)
  library(Shennong)
  library(glue)
  library(Seurat)
  
  outdir <- sn_set_path("data/processed/scrup/cellranger/corrected")
  dirs <- sn_list_10x_paths("/home/sduan/data/scrup/alignments/cellranger/")
  sn_write(dirs, "/home/sduan/data/scrup/alignments/cellranger.dirs.qs")
  samples <- names(dirs)
  
  excluded_samples <- c("SRX24386716", "SRX14989053", "SRX14989029",
                        "SRX14989027")
  samples <- samples[!(samples %in% excluded_samples)]
  length(samples)
  
  for (sample in samples) {
    print(sample)
    so_path <- glue("{outdir}/seurat_object/{sample}.qs")
    if (!file.exists(so_path)) {
      seurat_obj <- sn_initialize_seurat_object(x = dirs[sample])
      if (ncol(seurat_obj) > 100) {
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
  }
}
# cellranger --------------------------------------------------------------
if (FALSE) {
  library(tidyverse)
  library(Shennong)
  library(glue)
  library(Seurat)
  
  outdir <- sn_set_path("data/processed/20260328/scrup/cellranger")
  
  samples <- list.files("/home/sduan/data/scrup/alignments/cellranger", pattern = "^SRX")
  filtered_data_dirs <- glue(
    "/home/sduan/data/scrup/alignments/cellranger/{samples}/outs/filtered_feature_bc_matrix/"
  )
  names(filtered_data_dirs) <- samples
  filtered_data_dirs <- filtered_data_dirs[dir.exists(filtered_data_dirs)]
  samples <- names(filtered_data_dirs)
  
  excluded_samples <- c("SRX24386716")
  samples <- samples[!(samples %in% excluded_samples)]
  length(samples)
  # sample <- samples[1]
  for (sample in samples) {
    # sample <- "SRX24386716"
    print(sample)
    so_path <- glue("{outdir}/seurat_object/{sample}.qs")
    if (!file.exists(so_path)) {
      seurat_obj <- sn_initialize_seurat_object(x = filtered_data_dirs[sample],
                                                species = "human",
                                                min_features = 100)
      if (ncol(seurat_obj) > 100) {
        raw_counts <- sn_read(str_replace(filtered_data_dirs[sample], "filtered", "raw"))
        seurat_obj <- sn_remove_ambient_contamination(seurat_obj, raw = raw_counts)
        seurat_obj <- sn_filter_cells(seurat_obj, features = c("nFeature_RNA", "percent.mt"))
        seurat_obj <- sn_find_doublets(seurat_obj, layer = "decontaminated_counts")
        print(table(seurat_obj$scDblFinder.class_corrected))
        sn_write(seurat_obj, path = so_path)
        
        # seurat_obj <- seurat_obj |>
        #   subset(percent.mt <= 20)
        # seurat_obj <- sn_find_doublets(seurat_obj, dbr_sd = 1)
        # print(table(seurat_obj$scDblFinder.class))
        # sn_write(seurat_obj, path = so_path)
        
        # bpcells_path <- glue("{outdir}/bpcells/{sample}")
        # LayerData(seurat_obj, layer = "counts") |>
        #   sn_write(path = bpcells_path,
        #            to = "bpcells",
        #            overwrite = TRUE)
        # counts <- sn_read(bpcells_path)
        # LayerData(seurat_obj, layer = "counts") <- counts
        # sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
      }
    }
  }
}

# STAR --------------------------------------------------------------------
if (FALSE) {
  library(tidyverse)
  library(Shennong)
  library(glue)
  library(Seurat)
  
  outdir <- sn_set_path("data/processed/scrup/")
  
  samples <- list.files("/home/sduan/data/scrup/alignments", pattern = "^SRX")
  filtered_data_dirs <- glue("/home/sduan/data/scrup/alignments/{samples}/outs/Gene/filtered/")
  names(filtered_data_dirs) <- samples
  # sample <- samples[1]
  for (sample in samples) {
    print(sample)
    so_path <- glue("{outdir}/seurat_object/{sample}.qs")
    if (!file.exists(so_path)) {
      seurat_obj <- sn_initialize_seurat_object(x = filtered_data_dirs[sample],
                                                species = "human",
                                                min_features = 200)
      seurat_obj <- seurat_obj |>
        subset(percent.mt <= 20)
      seurat_obj <- sn_find_doublets(seurat_obj, dbr_sd = 1)
      print(table(seurat_obj$scDblFinder.class))
      sn_write(seurat_obj, path = so_path)
      
      bpcells_path <- glue("{outdir}/bpcells/{sample}")
      LayerData(seurat_obj, layer = "counts") |>
        sn_write(path = bpcells_path,
                 to = "bpcells",
                 overwrite = TRUE)
      counts <- sn_read(bpcells_path)
      LayerData(seurat_obj, layer = "counts") <- counts
      sn_write(seurat_obj, path = glue("{outdir}/bpcells/{sample}.qs"))
      
    }
    
  }
  
}
