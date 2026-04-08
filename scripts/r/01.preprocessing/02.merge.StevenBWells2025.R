if (TRUE) {
  library(tidyverse)
  library(Seurat)
  library(qs)
  library(glue)
  library(Shennong)
  
  outdir <- sn_set_path(
    "/home/sduan/projects/immune-atlas/data/processed/StevenBWells2025/batches/corrected"
  )
  
  files <- list.files(
    "/home/sduan/projects/immune-atlas/data/processed/StevenBWells2025/corrected/bpcells",
    pattern = "qs",
    full.names = TRUE
  )
  head(files)
  
  batches <- split(files, cut(seq_along(files), 10, labels = FALSE))
  
  walk2(batches, seq_along(batches), function(batch_files, i) {
    qs_file <- glue("{outdir}/batch_{i}.qs")
    if (!file.exists(qs_file)) {
      print(i)
      batch_files <- batches[[i]]
      seurat_obj_lists <- map(batch_files, function(x) {
        # x <- batch_files[1]
        s <- sn_read(x)
        # s@assays$HTO <- NULL
        # s <- JoinLayers(s, assay = "RNA")
        # s <- JoinLayers(s, assay = "ADT")
        # s <- JoinLayers(s, assay = "HTO")
        s
      })
      
      merged <- merge(x = seurat_obj_lists[[1]],
                      y = seurat_obj_lists[-1],
                      add.cell.ids = str_remove_all(basename(batch_files),".qs"))
      LayerData(merged, "data") <- NULL
      LayerData(merged, "scale.data") <- NULL
      
      merged <- JoinLayers(
        object = merged,
        assay = "RNA",
        layers = c("counts", "decontaminated_counts")
      )
      merged <- JoinLayers(object = merged, assay = "ADT")
      
      sn_write(merged, path = qs_file)
      
      message("Saved batch ", i, " with ", length(batch_files), " samples")
    }
  })
}

if (TRUE) {
  # library(tidyverse)
  library(dplyr)
  library(purrr)
  library(Seurat)
  library(glue)
  library(Shennong)
  
  seurat_obj_lists <- list.files(
    "/home/sduan/projects/immune-atlas/data/processed/StevenBWells2025/batches/corrected",
    pattern = ".qs$",
    full.names = TRUE
  ) %>%
    map(sn_read)
  merged <- merge(seurat_obj_lists[[1]], seurat_obj_lists[-1])
  merged <- JoinLayers(merged)
  sn_write(
    merged,
    "/home/sduan/projects/immune-atlas/data/processed/StevenBWells2025/merged_corrected.qs"
  )
}