library(Seurat)
library(tidyverse)
library(glue)

sample <- "582C_CZI-IA9924320"
file <- glue("/home/sduan/ia/alignments/count/{sample}/outs/filtered_feature_bc_matrix.h5")

samples <- list.files("/home/sduan/ia/alignments/count")
files <- glue("/home/sduan/ia/alignments/count/{samples}/outs/filtered_feature_bc_matrix.h5")
names(files) <- samples

files <- files[file.exists(files)]

walk(names(files), function(sample){
  # sample <- "759B_CZI-IA12953809"
  print(sample)
  seurat_obj_path <- glue("/home/sduan/ia/alignments/count/{sample}/outs/seurat_obj.qs")
  
  if (!file.exists(seurat_obj_path)) {
    file <- files[sample]
    data <- Read10X_h5(file)
    counts <- data$`Gene Expression`
    ac <- data$`Antibody Capture`
    adt_ids <- rownames(ac)[str_detect(rownames(ac), "TotalSeq")]
    hto_ids <- rownames(ac)[!(rownames(ac) %in% adt_ids)]
    
    obs <- read.csv(glue("/home/sduan/ia/alignments/count/{sample}/outs/obs.csv.gz"),
                    row.names = 1)
    metadata <- obs %>% 
      select(hash_tag, hash_status)
    
    if (length(hto_ids) > 1) {
      hto <- ac[rownames(ac) %in% hto_ids,]
    } else {
      metadata$hash_tag <- hto_ids
      metadata$hash_status <- "singlet"
    }
    adt <- ac[rownames(ac) %in% adt_ids,]
    
    seurat_obj <- CreateSeuratObject(counts, meta.data = metadata)
    seurat_obj[["ADT"]] <- CreateAssay5Object(counts = adt)
    if (length(hto_ids) > 1) {
      seurat_obj[["HTO"]] <- CreateAssay5Object(counts = hto)
    }
    seurat_obj
    colnames(seurat_obj) <- paste0(seurat_obj$hash_tag, "_", colnames(seurat_obj))
    print(table(seurat_obj$hash_tag))
    qs::qsave(seurat_obj, seurat_obj_path)
    gc()
  } else {
    print("Skip")
  }
 
})

seurat_obj_path <- glue("/home/sduan/ia/alignments/count/{samples}/outs/seurat_obj.qs")

seurat_obj <- qs::qread(seurat_obj_path[[1]])

seurat_obj

seurat_obj_lists <- map(seurat_obj_path, function(x) {
  print(x)
  s <- qs::qread(x)
  s <- s %>% 
    subset(hash_status == "singlet")
  s
})


