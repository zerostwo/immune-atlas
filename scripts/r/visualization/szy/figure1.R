library(Seurat)
library(tidyverse)
library(Shennong)
library(glue)
library(patchwork)

features <- c(
  "PTPRC",
  "CD3D","CD3E",
  "CD8A","CD8B",
  "CD79A","MS4A1",
  "GNLY","GZMB",
  "LYZ","AIF1",
  "EPCAM","KRT19",
  "TAGLN","COL1A2",
  "PECAM1","VWF",
  "HBB","HBA1"
)

outdir <- sn_set_path("results/figures/szy/figure1")

seurat_obj <- sn_read("/home/sduan/projects/immune-atlas/data/processed/scrup/20251117.integrated.qs")
seurat_obj@misc$cell_type_level1_colors <- c(
  "#7ec1e4", # CD4+ T cells
  "#335aa1", # CD8+ T cells
  "#c398f2", # B cells
  "#1d6dab", # Innate lymphoid cells
  "#5f358f", # Myeloid cells
  "#e8cf71", # Epithelial cells
  "#dd7c3f", # Stromal cells
  "#fc8f8f", # Endothelial cells
  "#ea3c30"  # Erythroid cells
)

cell_type_level1_colors <- seurat_obj@misc$cell_type_level1_colors
Idents(seurat_obj) <- "cell_type_level1"
names(cell_type_level1_colors) <- levels(seurat_obj)

set.seed(717)
disease_colors <- colorRampPalette(seurat_obj@misc$cell_type_level1_colors)(23) |> 
  sample(23)
names(disease_colors) <- unique(seurat_obj$harmonized_disease)

# seurat_obj <- seurat_obj |> 
#     subset(cell_type_level1 != "Erythroid cells")
p1 <- sn_plot_dim(seurat_obj, pt_size = 0.5, label = TRUE) + NoAxes() +
  labs(title = "Cell type (n = 3,363,615)") +
  scale_color_manual(values = cell_type_level1_colors) + NoLegend()
p1
ggsave(filename = glue("{outdir}/figure1a.pdf"), p1, width = 3, height = 3)

p2 <- sn_plot_dot(seurat_obj, features = features, group_by = "cell_type_level1") +
  scale_color_distiller(palette = "Purples",direction = 1)
p2
ggsave(filename = glue("{outdir}/figure1b.pdf"), p2, width = 5, height = 2.5)

metadata <- seurat_obj@meta.data

round(prop.table(table(metadata$bacteria))*100, 2)
p1 <- sn_plot_dim(
  seurat_obj, group_by = "bacteria", pt_size = 0.5,title = "Bacteria (7.77%)"
) + scale_color_manual(values = c("#6b2578", "lightgrey")) +
  NoAxes()
p1

round(prop.table(table(metadata$viruses))*100, 2)
p2 <- sn_plot_dim(
  seurat_obj, group_by = "viruses", pt_size = 0.5,title = "Viruses (0.02%)"
) + scale_color_manual(values = c("#6b2578", "lightgrey")) +
  NoAxes()

round(prop.table(table(metadata$eukaryota))*100, 2)
p3 <- sn_plot_dim(
  seurat_obj, group_by = "eukaryota", pt_size = 0.5,title = "Eukaryota (1.15%)"
) + scale_color_manual(values = c("#6b2578", "lightgrey")) +
  NoAxes()
p <- p1+p2+p3 + plot_layout(guides = "collect")&theme(
  legend.position = "bottom"
)
p
ggsave(filename = glue("{outdir}/figure1c.pdf"), p, width = 6, height = 3)

p1 <- metadata |> 
  filter(bacteria == "Detected") |>
  count(cell_type_level1) |> 
  mutate(cell_type_level1 = fct_reorder(
    cell_type_level1, n
  )) |>
  ggplot(aes(y = cell_type_level1, x=n/1000,fill = cell_type_level1)) +
  geom_col(width = 0.7) +
  geom_text(aes(x = 0.1,label = cell_type_level1), size = 8/.pt, hjust = 0) +
  scale_fill_manual(values = cell_type_level1_colors) +
  catplot::theme_cat(show_title = "x",show_text = "x") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(
    0,70000/1000
  ), breaks = c(0,10,20,30,40,50,60,70)) +
  NoLegend() +
  labs(x = "Number of cells (thousands)", title = "Bacteria")
p1
p2 <- metadata |> 
  count(cell_type_level1, bacteria) |> 
  pivot_wider(
    names_from = bacteria,
    values_from = n
  ) |> 
  mutate(prop = Detected / (Detected + `No detected`)*100) |> 
  mutate(prop = if_else(
    cell_type_level1 == "Erythrocytes", prop/10, prop
  )) |> 
  mutate(
    cell_type_level1 = fct_reorder(
      cell_type_level1, prop
    )
  ) |> 
  ggplot(aes(y = cell_type_level1, x = prop, fill = cell_type_level1)) +
  geom_col(width = 0.7) +
  geom_text(aes(x = 0.1,label = cell_type_level1), size = 8/.pt, hjust = 0) +
  NoLegend() +
  labs(
    x = "Proportion (%)", title = "Bacteria"
  ) +
  catplot::theme_cat(show_title = "x",show_text = "x") +
  scale_fill_manual(values = cell_type_level1_colors)  +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(
    0,20
  )) + NoLegend()
p2
p3 <- table(metadata$bacteria, metadata$harmonized_disease) |>
  prop.table(margin = 2) |> 
  as.data.frame() |> 
  filter(Var1 == "Detected") |> 
  mutate(Freq = Freq *100,
         Var2 = fct_reorder(Var2, Freq)) |>
  ggplot(aes(x = Freq, y = Var2, fill = Var1)) +
  geom_col(fill= "#6b2578", width = 0.7) +
  catplot::theme_cat(show_title = "x") +
  labs(
    x = "Proportion (%)",
    y = "Disease",
    title = "Bacteria"
  )  +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(
    0,40
  )) 
p3

adata <- metadata |>
  count(cell_type_level1, bacteria, harmonized_disease, sample_type)
p4 <- adata |>
  filter(bacteria == "Detected") |>
  ggplot(aes(y = harmonized_disease, x=n,fill = cell_type_level1)) +
  geom_col(position = "fill",  width = 0.7) +
  scale_fill_manual(values = cell_type_level1_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)))+
  catplot::theme_cat(show_title = "x", legend_position = "top", show_legend_title = FALSE) +
  labs(x = "Proportion (%)") + 
  NoLegend()

p <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
p
ggsave(filename = glue("{outdir}/figure1e.pdf"), p, width = 8, height = 3)


adata_a <- table(seurat_obj$orig.ident, seurat_obj$bacteria) |> 
  as.data.frame() |> 
  pivot_wider(
    names_from = Var2,
    values_from = Freq
  ) |> 
  mutate(
    prop = Detected / (Detected + `No detected`)*100
  ) |> 
  rename(sample = Var1) 
adata_b <- seurat_obj@meta.data |> 
  select(
    orig.ident, assay, harmonized_disease, harmonized_tissue, harmonized_sample_type 
  ) |> distinct()
adata <- adata_a |>
  full_join(
    adata_b, by = c("sample"="orig.ident")
  ) 

x_order <- adata |> 
  group_by(harmonized_disease) |> 
  summarise(prop = median(prop, na.rm = TRUE)) |> 
  arrange(desc(prop)) |> 
  pull(harmonized_disease)
adata$disease <- adata$harmonized_disease
adata$sample_type <- adata$harmonized_sample_type

set.seed(717)
p1 <- adata |> 
  mutate(
    disease = factor(disease, levels = x_order)
  ) |>
  ggplot(aes(x = disease, y = prop, color = disease)) + 
  geom_boxplot(outlier.shape = NA, width = 0.7, fatten = 1, staplewidth = 0.2) +
  geom_jitter(width = 0.3, size = 1) +
  catplot::theme_cat(show_title = "y", x_text_angle = 45) +
  labs(y = "Proportion (%)") +NoLegend() +
  scale_color_manual(values = disease_colors)
p1

p2 <- adata |> 
  filter(disease %in% c("ESCA")) |>
  ggplot(aes(x = sample_type, y = prop, color= sample_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fatten = 1, staplewidth = 0.2) +
  geom_jitter(width = 0.3, size = 1) +
  # facet_wrap(
  #   .~disease
  # ) + 
  ggsignif::geom_signif(comparisons = list(
    c("Normal", "Tumor")
  ),color="black", textsize = 8/.pt, size = 0.5/.pt) +
  catplot::theme_cat(x_text_angle = 45, show_title = "y") +
  NoLegend() +
  scale_color_manual(values = c(
    "lightgrey", "#6b2578"
  )) +
  labs(
    y = "Proportion (%)"
  ) +
  scale_y_continuous(limits = c(
    0, 100
  ))
p2

p3 <- adata |> 
  filter(sample_type %in% c("Normal", "Tumor")) |> 
  group_by(disease,sample_type) |> 
  summarise(
  mean = mean(prop, na.rm = TRUE)
) |> pivot_wider(
  names_from = sample_type,
  values_from = c(mean)
) |> pivot_longer(-disease) |> 
  mutate(value = if_else(is.na(value), 0, value)) |> 
  ungroup() |> 
  ggplot(aes(x = name, y = disease, size= value, color = value)) +
  geom_point() +
  scale_size_continuous(range = c(1, 3.5)) +
  catplot::theme_cat(show_panel_grid_marjor= "both",
                     x_text_angle = 45,show_title = "none") +
  scale_color_distiller(palette = "Purples", direction = 1) 
p3

p <- p1 + p3 + p2 + plot_layout(widths = c(6,1,1))
p
ggsave(filename = glue("{outdir}/figure1f.pdf"), p, width = 8, height = 3)

