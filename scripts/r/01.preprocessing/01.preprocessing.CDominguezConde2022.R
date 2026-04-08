library(tidyverse)
library(Shennong)
library(glue)
library(Seurat)

outdir <- sn_set_path("data/processed/CDominguezConde2022/corrected")

dirs <- sn_list_10x_paths("/mnt/rstudio/CDominguezConde2022/alignments/count")
samples <- names(dirs)
samples <- samples[!(samples %in% c("A29-SKM-0-SC-1", "A36-SKM-0-SC-1", "A37-DUO-0-SC-1"))]
for (sample in samples) {
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
