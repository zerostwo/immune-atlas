library(Seurat)
library(tidyverse)
library(Shennong)

seurat_obj <- sn_read("/home/duansq/projects/sahmi/data/raw/v1/pan_cancer.qs")
normal <- seurat_obj |> 
    subset(sample_type == "Normal")
normal_mat <- FetchData(normal, c("YTHDF2", "disease","cell_type_level1","study"))
normal_mat <- normal_mat |> 
    mutate(disease = paste0(disease,"_", study)) |>
    filter(cell_type_level1 != "Erythroid cells") |>
    group_by(disease, cell_type_level1) |> 
    summarise(YTHDF2 = mean(YTHDF2)) |> 
    pivot_wider(names_from = disease, values_from = YTHDF2) |> 
    column_to_rownames("cell_type_level1")

tumor <- seurat_obj |> 
    subset(sample_type == "Tumor")
tumor_mat <- FetchData(tumor, c("YTHDF2", "disease","cell_type_level1","study"))
tumor_mat <- tumor_mat |> 
    mutate(disease = paste0(disease,"_", study)) |>
    filter(cell_type_level1 != "Erythroid cells") |>
    group_by(disease, cell_type_level1) |> 
    summarise(YTHDF2 = mean(YTHDF2)) |> 
    pivot_wider(names_from = disease, values_from = YTHDF2) |> 
    column_to_rownames("cell_type_level1")

colnames(tumor_mat) <- paste0(colnames(tumor_mat),"_Tumor")
colnames(normal_mat) <- paste0(colnames(normal_mat),"_Normal")
mat <- cbind(tumor_mat, normal_mat)
pheatmap::pheatmap(
    mat,
    scale = "column",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
    fontsize = 8,
    border_color = "black",
    treeheight_row = 8,
    treeheight_col = 8,
    cellwidth = 10,
    cellheight = 10,
    filename = "./pan-cancer_YTHDF2_heatmap.pdf"
)

a <- table(seurat_obj$cell_type_level1, seurat_obj$sample) |> 
    prop.table(margin = 2) |> 
    as.data.frame()
colnames(a) <- c("cell_type_level1", "sample", "proportion")
b <- FetchData(seurat_obj, c("YTHDF2", "sample","cell_type_level1","disease"))
b <- b |>
    group_by(sample, cell_type_level1,disease) |> 
    summarise(YTHDF2 = mean(YTHDF2)) |> 
    ungroup()

adata <- left_join(a, b, by = c("sample", "cell_type_level1"))
head(adata)
p1 <- adata |> 
    filter(cell_type_level1 == "T cells") |>
    ggplot(aes(x = YTHDF2, y = proportion)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm",color = "#bd2e34") +
    catplot::theme_cat(show_panel_grid_marjor = "both",
                       show_panel_grid_minor = "both",
                       aspect_ratio = 1) +
    labs(x = "YTHDF2",
         y = "Proportion (%)",title = "T cell infiltration")
    # facet_wrap(~cell_type_level1, scales = "free") 
aa <- adata |> 
    filter(cell_type_level1 == "T cells") 
cor.test(aa$YTHDF2, aa$proportion)

tt <- seurat_obj |> 
    subset(cell_type_level1 == "T cells")
a <- FetchData(tt, c("sample","YTHDF2","REL",
                     "CD274","FUT4","TIGIT",
                     "IFNG","GZMB",
                     "LAG3","CTLA4", "PDCD1", "disease"))
a <- a |> 
    group_by(sample, disease) |> 
    summarise(
        YTHDF2 = mean(YTHDF2),
        REL = mean(REL),
        CD274 = mean(CD274),
        FUT4 = mean(FUT4),
        TIGIT = mean(TIGIT),
        LAG3 = mean(LAG3),
        CTLA4 = mean(CTLA4),
        PDCD1 = mean(PDCD1),
        IFNG = mean(IFNG),
        GZMB = mean(GZMB)
    ) |> 
    ungroup()

# calculate correlation of disese and gene YTHDF2
res <- map_dfr(
    c("BRCA","CRC","HNSC","LC","PAAD","RC"),
    function(d) {
        # d <- "LC"
        df <- a |> 
            filter(disease == d)
        temp_res <- map_dfr(
            c("REL",
              "CD274","FUT4","TIGIT",
              "LAG3","CTLA4", "PDCD1",
              "IFNG","GZMB"),
            function(gene) {
                # gene <- "REL"
                cor.test(
                    df[[gene]],
                    df$YTHDF2,
                    # method = "pearson"
                ) |> 
                    broom::tidy() |> 
                    mutate(
                        gene = gene,
                        disease = d
                    )
            }
        )
        temp_res
    }
)
res |> 
    filter(p.value < 0.05) |>
    ggplot(aes(x = disease, y = gene)) +
    geom_point(aes(size = -log10(p.value), color = estimate)) 



res <- map_dfr(
    c("REL",
      "CD274","FUT4","TIGIT",
      "LAG3","CTLA4", "PDCD1",
      "IFNG","GZMB"),
    function(gene) {
        model <- lm(as.formula(glue::glue("{gene} ~ YTHDF2")), data = a)
        summary_model <- summary(model)
        data.frame(
            gene = gene,
            estimate = summary_model$coefficients["YTHDF2", "Estimate"],
            p_value = summary_model$coefficients["YTHDF2", "Pr(>|t|)"],
            r_squared = summary_model$r.squared
        )
    }
)
p2 <- res |> 
    filter(!(gene %in% c("FUT4", "CD274","GZMB"))) |> 
    mutate(gene = fct_reorder(gene, -p_value)) |>
    ggplot(aes(x =-log10(p_value), y = gene, size = estimate)) +
    geom_segment(
        aes(x = 0, xend = -log10(p_value), y = gene, yend = gene),
        linewidth = 0.8,
        color = "lightgrey",
        show.legend = FALSE
    ) +
    geom_point(
        aes(size = estimate),
        color = "black",
        fill = "#bd2e34",
        shape = 21,
        stroke = 0.4
    ) +
    scale_x_continuous(limits = c(0, 5),expand = c(0,0)) +
    catplot::theme_cat(aspect_ratio = 1.5)  +
    labs(
        x = "-log10(p value)",
        y = "Genes",
        size = "Estimate", title = "Correlation with YTHDF2"
    )

p <- p1 + p2
ggsave("./pan-cancer_YTHDF2_Tcell_infiltration_correlation.pdf", p, width = 6, height = 2)
