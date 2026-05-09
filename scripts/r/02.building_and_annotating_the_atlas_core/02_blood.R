library(Shennong)
library(tidyverse)
library(Seurat)

outdir <- sn_set_path("data/processed/atlas_core/blood")

merged <- sn_read("data/raw/reference/reference.raw.qs")

ss <- merged %>% 
  # subset(sample == "TSP21_Blood_NA_10X_1_1")
  subset(study == "ZhenlongLi2024")

blood <- merged %>%
  subset(dataset_tier == "Core") %>%       # nFeature_RNA_corrected > 1000
  subset(dataset_role == "Reference") %>%  # sample from health state
  subset(dataset_set == "Discovery") %>%   # exclude validation set
  subset(tissue_level1 == "Blood")         # focus on blood

# filter out doublets and samples with too few cells
blood <- blood %>%
  subset(scDblFinder.class_corrected == "singlet")

keep_samples <- blood@meta.data %>%
  count(sample) %>%
  filter(n > 100) %>%
  pull(sample)

blood <- blood %>%
  subset(sample %in% keep_samples)

sn_write(blood, path = file.path(outdir, "blood.qs"))

# sample-level pseudobulk expression PCA ----------------------------------
blood$sample <- paste0(blood$study, "_", blood$sample)


robust_scale <- function(x) {
  s <- mad(x, constant = 1, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - median(x, na.rm = TRUE)) / s
}

sample_qc <- blood@meta.data %>%
  group_by(sample, study, assay, donor) %>%
  summarise(
    n_cells = n(),
    median_nCount = median(nCount_RNA_corrected, na.rm = TRUE),
    median_nFeature = median(nFeature_RNA_corrected, na.rm = TRUE),
    median_mt = median(percent.mt, na.rm = TRUE),
    median_ribo = median(percent.ribo, na.rm = TRUE),
    median_hb = median(percent.hb, na.rm = TRUE),
    median_decontX = median(decontX_contamination, na.rm = TRUE),
    median_dbl = median(scDblFinder.score_corrected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    qc_score =
      robust_scale(log10(n_cells + 1)) +
      robust_scale(median_nFeature) +
      0.5 * robust_scale(median_nCount) -
      robust_scale(median_mt) -
      robust_scale(median_hb) -
      robust_scale(median_decontX) -
      robust_scale(median_dbl)
  )

primary_samples <- sample_qc %>%
  filter(
    n_cells >= 300,
    median_mt <= quantile(median_mt, 0.9, na.rm = TRUE),
    median_decontX <= quantile(median_decontX, 0.9, na.rm = TRUE),
    median_dbl <= quantile(median_dbl, 0.9, na.rm = TRUE)
  ) %>%
  pull(sample)

blood_primary <- subset(blood, sample %in% primary_samples)

DefaultAssay(blood_primary) <- "RNA"

blood_primary <- blood_primary %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(
    features = VariableFeatures(.),
    vars.to.regress = c("percent.mt", "nCount_RNA_corrected")
  ) %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.3)

sample_comp <- blood_primary@meta.data %>%
  count(sample, seurat_clusters, name = "n_cluster_cells") %>%
  group_by(sample) %>%
  mutate(prop = n_cluster_cells / sum(n_cluster_cells)) %>%
  ungroup() %>%
  select(sample, seurat_clusters, prop) %>%
  pivot_wider(
    names_from = seurat_clusters,
    values_from = prop,
    values_fill = 0,
    names_prefix = "comp_cluster_"
  )

# library(edgeR)
# library(limma)

pb_counts <- AggregateExpression(
  blood_primary,
  group.by = "sample",
  assays = "RNA",
  slot = "counts",
  verbose = FALSE
)$RNA

pb_meta <- blood_primary@meta.data %>%
  mutate(sample = str_replace_all(sample, "_", "-")) %>%
  distinct(sample, study, assay, donor, gender, age) %>%
  slice(match(colnames(pb_counts), sample))

keep_genes <- rowSums(pb_counts >= 10) >= max(3, ceiling(0.1 * ncol(pb_counts)))
pb_counts <- pb_counts[keep_genes, , drop = FALSE]

dge <- edgeR::DGEList(pb_counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
pb_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

pb_logcpm_corr <- limma::removeBatchEffect(pb_logcpm, batch = pb_meta$assay)

gene_var <- apply(pb_logcpm_corr, 1, var)
top_pb_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(2000, length(gene_var))]

pb_pca <- prcomp(t(pb_logcpm_corr[top_pb_genes, ]), center = TRUE, scale. = TRUE)

sample_expr <- as_tibble(pb_pca$x[, 1:10], rownames = "sample") %>%
  rename_with( ~ paste0("expr_", .x), -sample)

sample_repr <- sample_qc %>%
  filter(sample %in% primary_samples) %>%
  left_join(sample_comp, by = "sample") %>%
  mutate(sample = str_replace_all(sample, "_", "-")) %>%
  left_join(sample_expr, by = "sample")

pick_medoid <- function(mat) {
  d <- as.matrix(dist(mat))
  rownames(mat)[which.min(rowSums(d))]
}

feature_cols <- c(grep("^expr_PC", names(sample_repr), value = TRUE),
                  grep("^comp_cluster_", names(sample_repr), value = TRUE))

sample_mat <- sample_repr %>%
  select(all_of(feature_cols)) %>%
  as.matrix()

rownames(sample_mat) <- sample_repr$sample
sample_mat <- scale(sample_mat)

global_medoid <- pick_medoid(sample_mat)
global_medoid

# 方法 A：按 sample-level clustering 选每簇 medoid
sample_hc <- hclust(dist(sample_mat), method = "ward.D2")

k <- 6
sample_cluster <- cutree(sample_hc, k = k)

sample_repr <- sample_repr %>%
  mutate(sample_cluster = sample_cluster[sample])

cluster_medoids <- split(sample_repr$sample, sample_repr$sample_cluster) %>%
  purrr::imap_dfr( ~ {
    tibble(sample_cluster = .y, medoid = pick_medoid(sample_mat[.x, , drop = FALSE]))
  })

cluster_medoids
# 方法 B：按 study 再选 medoid
study_medoids <- sample_repr %>%
  group_by(study) %>%
  group_modify( ~ {
    ids <- .x$sample
    tibble(medoid = pick_medoid(sample_mat[ids, , drop = FALSE]))
  }) %>%
  ungroup()

study_medoids

candidate_samples <- c(global_medoid, cluster_medoids$medoid, study_medoids$medoid) %>%
  unique()

candidate_table <- sample_repr %>%
  filter(sample %in% candidate_samples) %>%
  arrange(desc(qc_score)) %>%
  select(
    sample,
    study,
    assay,
    donor,
    n_cells,
    median_nFeature,
    median_mt,
    median_decontX,
    median_dbl,
    qc_score
  )

candidate_table

sample_repr %>%
  mutate(label = if_else(sample %in% candidate_table$sample, sample, "")) %>%
  ggplot(aes(x =  expr_PC1, y = expr_PC2, colour = study)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = label), size = 4)


# Annotation --------------------------------------------------------------
blood@meta.data <- blood@meta.data %>%
  mutate(sample = str_replace_all(sample, "_", "-"))
global_medoid_seurat_obj <- subset(blood, sample == "NataliaJaeger2024-Blood1")
global_medoid_seurat_obj <- sn_filter_genes(global_medoid_seurat_obj, gene_class = "coding")
global_medoid_seurat_obj <- sn_run_cluster(global_medoid_seurat_obj,
                                           vars_to_regress = "percent.mt",
                                           layer = "decontaminated_counts")
DimPlot(global_medoid_seurat_obj, label = TRUE)
sn_plot_feature(global_medoid_seurat_obj,
                features = "IL4I1",
                pt_size = 2)
features <- marker_genes %>% 
  filter(high_hierarchy_cell_type %in% c("T cells", "ILC", "ILC precursor")) %>% 
  pull(human) %>% 
  unique()

sn_plot_dot(global_medoid_seurat_obj, features = features)

global_medoid_seurat_obj <- sn_find_de(global_medoid_seurat_obj, analysis = "markers")
sn_plot_dot(global_medoid_seurat_obj, features = c("PTPRC","IL7R",
                                                   "CD3D","CD3E","CD3G",
                                                   "CD8A","CD8B","CD4","CD40LG"))
top_markers <- sn_get_de_result(global_medoid_seurat_obj,top_n = 10)
split(top_markers$gene, top_markers$cluster)



if (FALSE) {
  # Compare counts and decontaminated_counts
  samples <- unique(merged$sample)
  adata <- map_dfr(samples, function(x) {
    print(x)
    s <- merged %>%
      subset(sample == x)
    tibble(
      sample = x,
      counts = sum(rowSums(LayerData(s, layer = "counts"))),
      decontaminated_counts = sum(rowSums(
        LayerData(s, layer = "decontaminated_counts")
      ))
    )
  })
}

