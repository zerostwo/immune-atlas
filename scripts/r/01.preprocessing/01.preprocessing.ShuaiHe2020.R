library(tidyverse)
library(Shennong)
library(glue)
library(Seurat)

outdir <- sn_set_path("data/processed/ShuaiHe2020/corrected")

dirs <- sn_list_10x_paths("/mnt/public/ShuaiHe2020/alignments")
samples <- names(dirs)

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

# Clustering --------------------------------------------------------------
if (FALSE) {
    seurat_obj <- sn_read("data/processed/ShuaiHe2020/corrected/bpcells/SRR13075720.qs")
    seurat_obj
    seurat_obj$scDblFinder.class_corrected |> table()
    seurat_obj <- seurat_obj |> 
        subset(nFeature_RNA_corrected > 200) |> 
        subset(scDblFinder.class_corrected == "singlet")
    seurat_obj <- sn_filter_genes(seurat_obj, gene_class = "coding")
    
    seurat_obj <- sn_run_cluster(seurat_obj, layer = "decontaminated_counts")
    sn_plot_dim(seurat_obj, label = TRUE, pt_size = 1.5)
    sn_plot_dim(seurat_obj, label = TRUE, pt_size = 1.5,
                group_by = "Phase")
    pp <- sn_calculate_rogue(seurat_obj, cluster = "seurat_clusters")
    sn_plot_feature(seurat_obj, features = "CLDN5", pt_size = 1.5)
}
