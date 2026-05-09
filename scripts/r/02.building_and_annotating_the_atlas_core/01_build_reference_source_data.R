library(tidyverse)
library(Seurat)
library(glue)
library(Shennong)

# Part 1 ------------------------------------------------------------------
if (TRUE) {
  files <- list.files(
    "data/processed/EMadissoon2020/corrected_soupx/bpcells",
    pattern = "qs2",
    full.names = TRUE
  )
  # files <- files[-3]
  EMadissoon2020 <- map(files, sn_read)
  EMadissoon2020 <- merge(
    x = EMadissoon2020[[1]],
    y = EMadissoon2020[-1],
    add.cell.ids = str_remove_all(basename(files), ".qs2")
  )
  EMadissoon2020 <- JoinLayers(EMadissoon2020, layers = c("counts", "decontaminated_counts"))
  EMadissoon2020$study <- "EMadissoon2020"
  EMadissoon2020
  
  files <- list.files(
    "data/processed/DanielPCaron2025/corrected_soupx/bpcells",
    pattern = "qs2",
    full.names = TRUE
  )
  # files <- files[-3]
  DanielPCaron2025 <- map(files, sn_read)
  DanielPCaron2025 <- merge(
    x = DanielPCaron2025[[1]],
    y = DanielPCaron2025[-1],
    add.cell.ids = str_remove_all(basename(files), ".qs2")
  )
  DanielPCaron2025 <- JoinLayers(DanielPCaron2025,
                                 layers = c("counts", "decontaminated_counts"))
  DanielPCaron2025 <- JoinLayers(DanielPCaron2025, assay = "ADT")
  DanielPCaron2025$study <- "DanielPCaron2025"
  DanielPCaron2025
  adt_counts <- LayerData(DanielPCaron2025, assay = "ADT", layer = "counts")
  rownames(adt_counts) <- str_replace(rownames(adt_counts), "-TotalSeqA", "")
  # 1. Define a standardization dictionary
  adt_name_map <- c(
    "B7-H4" = "VTCN1",
    "Cadherin-11" = "CDH11",
    "Galectin-9" = "LGALS9",
    "Mac-2" = "LGALS3",
    "Notch-1" = "NOTCH1",
    "Notch-3" = "NOTCH3",
    "TNF-alpha" = "TNF",
    "VEGFR-3" = "FLT4",
    "TIM-4" = "TIMD4",
    "NKp80" = "KLRF1",
    "TSLPR" = "CRLF2",
    "LOX-1" = "OLR1",
    "Podoplanin" = "PDPN",
    "GARP" = "LRRC32",
    "LAP" = "TGFB1-LAP",
    "GP130" = "IL6ST",
    "IgE" = "IGHE",
    "TRA-TRB" = "TCRab",
    "TRG-TRD" = "TCRgd",
    "HLA-ABC" = "HLA-I",
    "HLA-DR-DP-DQ" = "HLA-II",
    "HLA-DR" = "HLA-DR",
    "MICA-MICB" = "MICA-MICB"
  )
  
  # 2. Define features to drop
  drop_pattern <- regex(paste(
    c(
      "isotype",
      "^IgG-Fc$",
      "^IGHG_FC$",
      "^hamster-igg",
      "^igg[0-9a-z-]*",
      "^TotalSeq[A-C]$"
    ),
    collapse = "|"
  ), ignore_case = TRUE)
  
  # 3. Clean and standardize ADT feature names
  adt_features_raw <- rownames(adt_counts)
  
  adt_feature_tbl <- tibble(feature_raw = adt_features_raw) |>
    mutate(
      # Remove TotalSeq suffix if still present
      feature_clean = str_remove(feature_raw, "-TotalSeq[A-C]$"),
      
      # Standardize a few style issues first
      feature_clean = str_replace_all(feature_clean, fixed("–"), "-"),
      feature_clean = str_squish(feature_clean),
      
      # Map non-standard names to more standard names
      feature_std = recode(feature_clean, !!!adt_name_map, .default = feature_clean),
      
      # Mark features to drop
      drop = str_detect(feature_std, drop_pattern)
    )
  
  # View mapping table
  adt_feature_tbl
  
  # 4. Keep only useful features
  adt_counts2 <- adt_counts[!adt_feature_tbl$drop, , drop = FALSE]
  
  # Replace rownames with standardized names
  rownames(adt_counts2) <- adt_feature_tbl |>
    filter(!drop) |>
    pull(feature_std)
  
  # Check results
  cat(rownames(adt_counts2), sep = "\n")
  DanielPCaron2025[["ADT"]] <- CreateAssay5Object(counts = adt_counts2)
  
  files <- list.files(
    "data/processed/ShuaiHe2020/corrected_soupx/bpcells",
    pattern = "qs",
    full.names = TRUE
  )
  ShuaiHe2020 <- map(files, sn_read)
  ShuaiHe2020 <- merge(
    x = ShuaiHe2020[[1]],
    y = ShuaiHe2020[-1],
    add.cell.ids = str_remove_all(basename(files), ".qs")
  )
  ShuaiHe2020 <- JoinLayers(ShuaiHe2020, layers = c("counts", "decontaminated_counts"))
  ShuaiHe2020$study <- "ShuaiHe2020"
  ShuaiHe2020
  
  files <- list.files(
    "data/processed/MarylineFalquet2023/corrected_soupx/bpcells",
    pattern = "qs",
    full.names = TRUE
  )
  MarylineFalquet2023 <- map(files, sn_read)
  MarylineFalquet2023 <- merge(
    x = MarylineFalquet2023[[1]],
    y = MarylineFalquet2023[-1],
    add.cell.ids = str_remove_all(basename(files), ".qs")
  )
  MarylineFalquet2023 <- JoinLayers(MarylineFalquet2023,
                                    layers = c("counts", "decontaminated_counts"))
  MarylineFalquet2023$study <- "MarylineFalquet2023"
  MarylineFalquet2023
  
  files <- list.files(
    "data/processed/AndrewDHildreth2021/corrected_soupx/bpcells",
    pattern = "qs",
    full.names = TRUE
  )
  AndrewDHildreth2021 <- map(files, sn_read)
  AndrewDHildreth2021 <- merge(
    x = AndrewDHildreth2021[[1]],
    y = AndrewDHildreth2021[-1],
    add.cell.ids = str_remove_all(basename(files), ".qs")
  )
  AndrewDHildreth2021 <- JoinLayers(AndrewDHildreth2021,
                                    layers = c("counts", "decontaminated_counts"))
  AndrewDHildreth2021$study <- "AndrewDHildreth2021"
  AndrewDHildreth2021
  
  files <- list.files(
    "data/processed/CDominguezConde2022/corrected_soupx/bpcells",
    pattern = "qs",
    full.names = TRUE
  )
  CDominguezConde2022 <- map(files, sn_read)
  CDominguezConde2022 <- merge(
    x = CDominguezConde2022[[1]],
    y = CDominguezConde2022[-1],
    add.cell.ids = str_remove_all(basename(files), ".qs")
  )
  CDominguezConde2022 <- JoinLayers(CDominguezConde2022,
                                    layers = c("counts", "decontaminated_counts"))
  CDominguezConde2022$study <- "CDominguezConde2022"
  CDominguezConde2022
  
  files <- list.files(
    "data/processed/NataliaJaeger2024/corrected_soupx/bpcells",
    pattern = "qs",
    full.names = TRUE
  )
  NataliaJaeger2024 <- map(files, sn_read)
  NataliaJaeger2024 <- merge(
    x = NataliaJaeger2024[[1]],
    y = NataliaJaeger2024[-1],
    add.cell.ids = str_remove_all(basename(files), ".qs")
  )
  NataliaJaeger2024 <- JoinLayers(NataliaJaeger2024,
                                  layers = c("counts", "decontaminated_counts"))
  NataliaJaeger2024$study <- "NataliaJaeger2024"
  NataliaJaeger2024@assays$ADT <- NULL
  NataliaJaeger2024@assays$HTO <- NULL
  # NataliaJaeger2024@assays$RNA$scale.data <- NULL
  NataliaJaeger2024@assays$RNA$scale.data.17 <- NULL
  NataliaJaeger2024@assays$RNA$data.17 <- NULL
  # NataliaJaeger2024@assays$RNA$data <- NULL
  NataliaJaeger2024
  
  merged <- merge(
    x = CDominguezConde2022,
    y = list(
      AndrewDHildreth2021,
      NataliaJaeger2024,
      MarylineFalquet2023,
      ShuaiHe2020,
      DanielPCaron2025,
      EMadissoon2020
    )
  )
  merged <- JoinLayers(merged, layers = c("counts", "decontaminated_counts"))
  merged <- JoinLayers(merged, assay = "ADT")
  merged
  table(merged$study)
  
  prefix <- "AndrewDHildreth2021_CDominguezConde2022_NataliaJaeger2024_MarylineFalquet2023_ShuaiHe2020_DanielPCaron2025_EMadissoon2020"
  merged_path <- glue("data/processed/{prefix}/merged.qs2")
  
  # LayerData(merged, "data") <- NULL
  # LayerData(merged, "scale.data") <- NULL
  
  counts <- LayerData(merged, layer = "counts", assay = "RNA")
  BPCells::write_matrix_dir(counts,
                            glue("data/processed/{prefix}/bpcells/counts"),
                            overwrite = TRUE)
  counts <- sn_read(glue("data/processed/{prefix}/bpcells/counts"))
  LayerData(merged, layer = "counts", assay = "RNA") <- counts
  counts <- LayerData(merged, layer = "decontaminated_counts", assay = "RNA")
  BPCells::write_matrix_dir(
    counts,
    glue("data/processed/{prefix}/bpcells/decontaminated_counts"),
    overwrite = TRUE
  )
  counts <- sn_read(glue("data/processed/{prefix}/bpcells/decontaminated_counts"))
  LayerData(merged, layer = "decontaminated_counts", assay = "RNA") <- counts
  adt <- LayerData(merged, layer = "counts", assay = "ADT")
  BPCells::write_matrix_dir(adt,
                            glue("data/processed/{prefix}/bpcells/adt"),
                            overwrite = TRUE)
  adt <- sn_read(glue("data/processed/{prefix}/bpcells/adt"))
  LayerData(merged, layer = "counts", assay = "ADT") <- adt
  
  sn_write(merged, path = merged_path)
}
# Merge -------------------------------------------------------------------
if (FALSE) {
  merged1 <- sn_read(
    "data/processed/AndrewDHildreth2021_CDominguezConde2022_NataliaJaeger2024_MarylineFalquet2023_ShuaiHe2020_DanielPCaron2025_EMadissoon2020/merged.qs2"
  )
  merged2 <- sn_read("data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens.merged_corrected.qs2")
  counts <- sn_read(
    "data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens2022/bpcells_corrected/counts/"
  )
  LayerData(merged2, layer = "counts") <- counts
  
  counts <- sn_read(
    "data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens2022/bpcells_corrected/decontaminated_counts/"
  )
  LayerData(merged2, layer = "decontaminated_counts") <- counts
  
  adt_counts <- sn_read(
    "data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens2022/bpcells_corrected/adt_counts/"
  )
  rownames(adt_counts) <- str_replace(rownames(adt_counts), "-TotalSeqC", "")
  
  # 1. Rename dictionary
  adt_name_map <- c(
    "LOX-1"   = "OLR1",
    "HLA-ABC" = "HLA-I",
    "HLA-DR"  = "HLA-DR",
    "TRA-TRB" = "TCRab"
  )
  
  # 2. Define drop pattern for isotype / IgG controls
  drop_pattern <- regex(paste(
    c(
      "isotype",
      "^armenian-hamster-igg",
      "^mouse-igg[0-9a-z-]*",
      "^rat-igg[0-9a-z-]*",
      "^hamster-igg[0-9a-z-]*",
      "^igg[0-9a-z-]*"
    ),
    collapse = "|"
  ), ignore_case = TRUE)
  
  # 3. Build feature table
  feature_tbl <- tibble(feature_raw = rownames(adt_counts)) |>
    mutate(
      feature_clean = str_squish(feature_raw),
      feature_std = recode(feature_clean, !!!adt_name_map, .default = feature_clean),
      drop = str_detect(feature_std, drop_pattern)
    )
  
  # Check what will be dropped
  feature_tbl |>
    filter(drop)
  
  # 4. Apply filtering and renaming
  adt_counts_clean <- adt_counts[!feature_tbl$drop, , drop = FALSE]
  
  rownames(adt_counts_clean) <- feature_tbl |>
    filter(!drop) |>
    pull(feature_std)
  
  # LayerData(merged2, layer = "counts", assay = "ADT") <- adt_counts
  merged2[["ADT"]] <- CreateAssay5Object(counts = adt_counts_clean)
  
  merged <- merge(x = merged1, y = merged2)
  merged <- JoinLayers(merged, layers = c("counts", "decontaminated_counts"))
  merged <- JoinLayers(merged, assay = "ADT")
  sn_write(merged, "data/raw/reference/20260508.reference.raw.qs2")
  
  counts <- LayerData(merged, layer = "counts", assay = "RNA")
  BPCells::write_matrix_dir(counts, "data/raw/reference/bpcells/counts", overwrite = TRUE)
  counts <- sn_read("data/raw/reference/bpcells/counts")
  LayerData(merged, layer = "counts", assay = "RNA") <- counts
  
  counts <- LayerData(merged, layer = "decontaminated_counts", assay = "RNA")
  BPCells::write_matrix_dir(counts,
                            "data/raw/reference/bpcells/decontaminated_counts",
                            overwrite = TRUE)
  counts <- sn_read("data/raw/reference/bpcells/decontaminated_counts")
  LayerData(merged, layer = "decontaminated_counts", assay = "RNA") <- counts
  
  adt <- LayerData(merged, layer = "counts", assay = "ADT")
  BPCells::write_matrix_dir(adt, "data/raw/reference/bpcells/adt", overwrite = TRUE)
  adt <- sn_read("data/raw/reference/bpcells/adt")
  LayerData(merged, layer = "counts", assay = "ADT") <- adt
  
  sn_write(merged, "data/raw/reference/20260508.reference.raw.qs2")
  
  merged <- sn_read("data/raw/reference/reference.raw.qs2")
  
  table(merged$scDblFinder.class_corrected, merged$study) |>
    as.data.frame() |>
    ggplot(aes(x = Freq, y = Var2, fill = Var1)) +
    geom_col(position = "fill") +
    catplot::theme_cat(show_title = "x", aspect_ratio = 1) +
    scale_fill_brewer(palette = "Paired") +
    scale_x_continuous(expand = c(0, 0))
  
}

# QC ----------------------------------------------------------------------
if (FALSE) {
  merged <- sn_read("data/raw/reference/reference.raw.qs2")
  glimpse(merged@meta.data)
  quantile(merged$percent.mt)
  quantile(merged$nFeature_RNA)
  quantile(merged$nCount_RNA)
  table(merged$study)
  table(merged$scDblFinder.class_corrected)
 
  ia_sample_metadata <- sn_read("data/metadata//IA_sample_spreadsheet.xlsx", sheet = 2)
  sample_metadata <- sn_read("data/metadata/20260101.reference.metadata.xlsx", sheet = "sample")
  sample_metadata <- sample_metadata |>
    left_join(
      ia_sample_metadata |> select(
        sample_id,
        fresh_or_frozen,
        sorting,
        stimulation,
        cell_type,
        gex_chem,
        organ
      ),
      by = c("sample" = "sample_id")
    ) |>
    mutate(
      assay = if_else(is.na(assay), gex_chem, assay),
      tissue_level1 = if_else(is.na(tissue_level1), tissue_level2, tissue_level1)
    ) |>
    mutate(assay = case_when(
      assay == "5'v2" ~ "10x5'v2",
      TRUE ~ assay
    ))
  
  sample_metadata |>
    count(tissue_level1, study) |>
    ggplot(aes(x = study, y = tissue_level1, size = n)) +
    geom_point() +
    catplot::theme_cat(
      show_panel_grid_marjor = "both",
      show_title = "none",
      x_text_angle = 45
    ) +
    coord_fixed()
  
  donor_metadata <- sn_read("data/metadata/20260101.reference.metadata.xlsx", sheet = "donor")
  donor_metadata <- donor_metadata |>
    mutate(age = round(as.numeric(age), 1),
           BMI = round(as.numeric(BMI), 2))
  donor_metadata |>
    ggplot(aes(x = BMI)) +
    geom_histogram(bins = 14, color = "white") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 15))
  
  donor_metadata |>
    ggplot(aes(x = age)) +
    geom_histogram(bins = 14, color = "white") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 20))
  
  donor_metadata |>
    count(gender) |>
    ggplot(aes(x = gender, y = n)) +
    geom_col()
  
  donor_metadata |>
    count(ethnicity) |>
    ggplot(aes(x = ethnicity, y = n)) +
    geom_col() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
    catplot::theme_cat(x_text_angle = 45)
  
  metadata <- donor_metadata |>
    left_join(sample_metadata, by = c("study", "donor")) |>
    mutate(assay = case_when(gex_chem == "5'v1" ~ "10x5'v1", TRUE ~ assay)) |>
    select(-c(organ, gex_chem, order, note)) |>
    rename(
      bmi = BMI,
      tobacco = Tobacco,
      alcohol = Alcohol,
      cmv = CMV,
      ebv = EBV,
      diabetes = Diabetes,
      cancer = Cancer,
      hypertension = Hypertension,
      iv_drug_abuse = IV_drug_abuse,
      toxo = TOXO
    )
  
  merged@meta.data <- merged@meta.data |>
    select(
      -c(
        hash.ID,
        # prefix,
        # barcode,
        # donor
        nCount_HTO,
        nFeature_HTO,
        HTO_maxID,
        HTO_secondID,
        HTO_margin,
        HTO_classification,
        HTO_classification.global
      )
    )
  
  merged$dataset_role <- "Reference" # Reference, Query, or Other
  merged$dataset_tier <- "Core" # Core, Extended, or Other
  merged$dataset_set <- "Discovery" # Discovery, Validation, or Other
  merged@meta.data <- merged@meta.data %>% 
    mutate(study = if_else(study == "TabulaSapiens", "TabulaSapiens2022", study))
  merged@meta.data <- merged@meta.data |>
    rownames_to_column("cell_id") |>
    # select(-donor) |>
    left_join(metadata, by = c("study", "sample")) |>
    column_to_rownames("cell_id")
  
  merged@meta.data <- merged@meta.data |>
    mutate(
      donor = case_when(
        sample == "PBMC101923-TotalSeqB" ~ "PBMC101923-TotalSeqB",
        sample %in% c(
          "Central-Memory-CD4-T-cells",
          "Central-Memory-CD8-T-cells",
          "Effector-Memory-CD4-T-cells",
          "Effector-Memory-CD8-T-cells",
          "Monocytes",
          "CD8-TEMRA",
          "CD8-TEMRA",
          "Naive-CD4-T-cells",
          "Naive-CD8-T-cells"
        ) ~ "FACS",
        sample == "TS1101923-TotalSeqB" ~ "TS1101923-TotalSeqB",
        sample == "TS2101923-TotalSeqB" ~ "TS2101923-TotalSeqB",
        sample == "PBMC101923-TotalSeqB" ~ "PBMC101923-TotalSeqB",
        TRUE ~ donor
      ),
      tissue_level1 = case_when(
        donor == "FACS" ~ "Blood",
        sample == "PBMC101923-TotalSeqB" ~ "Blood",
        sample %in% c("TS1101923-TotalSeqB", "TS2101923-TotalSeqB") ~ "Tonsil",
        TRUE ~ tissue_level1
      ),
      assay = case_when(
        study == "DanielPCaron2025" ~ "10x3'v3.1",
        study == "EMadissoon2020" ~ "10x3'v2",
        sample %in% c(
          "TS1101923-TotalSeqB",
          "TS2101923-TotalSeqB",
          "PBMC101923-TotalSeqB"
        ) ~ "10x3'v3",
        TRUE ~ assay
      ),
      dataset_set = case_when(donor == "FACS" ~ "Validation", TRUE ~ dataset_set),
      dataset_role = case_when(
        sample %in% c(
          "CD_IEL2",
          "Lung2_fibrotic_recipient",
          "Lung4_fibrotic_recipient",
          "P",
          "P2",
          "P3",
          "P4"
        ) ~ "Query",
        TRUE ~ dataset_role
      )
    )
  
  merged@meta.data <- merged@meta.data %>% 
    mutate(
      dataset_tier = case_when(
        nFeature_RNA_corrected >= 1000 ~ "Core",
        nFeature_RNA_corrected >= 200 & nFeature_RNA_corrected < 1000 ~ "Extended",
        TRUE ~ "Other"
      )
    )
  table(merged$dataset_tier)
  
  merged$batch <- paste0(merged$study, "_", merged$assay)
  merged$sample <- paste0(merged$study, "_", merged$sample)
  
  # merged@meta.data |>
  #   select(study, donor, sample, tissue_level1, assay, disease) |>
  #   distinct() |> View()
  # 
  table(merged$tissue_level1, merged$study) |>
    prop.table(margin = 2) |>
    as.data.frame() |>
    # mutate(Var1 = if_else(Var1 == "Blood", "Blood", "Others")) %>% 
    ggplot(aes(x = Freq, y = Var2, fill = Var1)) +
    geom_col()
  
  # keep_samples <- merged@meta.data %>% 
  #   filter(dataset_tier == "Core",
  #          dataset_role == "Reference",
  #          dataset_set == "Discovery") %>% 
  #   count(sample) %>% 
  #   arrange(n)
  # merged@meta.data %>% 
  #   filter(dataset_tier == "Core",
  #          dataset_role == "Reference",
  #          dataset_set == "Discovery") %>% 
  #   count(scDblFinder.class_corrected)
  
  sn_write(merged, "data/raw/reference/20260508.reference.raw.qs2")
  
}

# Seed samples ------------------------------------------------------------
if (FALSE) {
  merged <- sn_read("data/raw/reference.qs2")
  NataliaJaeger2024 <- subset(merged, study == "NataliaJaeger2024")
  NataliaJaeger2024_seed <- NataliaJaeger2024 |>
    subset(scDblFinder.class == "singlet") |>
    subset(dataset_tier == "Core") |>
    subset(dataset_role == "Reference")
  LayerData(NataliaJaeger2024_seed, layer = "counts") <- as(LayerData(NataliaJaeger2024_seed, layer = "counts"), "dgCMatrix")
  sn_write(NataliaJaeger2024_seed, "data/raw/reference/seeds/NataliaJaeger2024.h5ad")
}

# Export h5ad -------------------------------------------------------------
if (FALSE) {
  merged <- sn_read("data/raw/reference.qs2")
  # Notably, HVG selection was performed after removing specific genes, including
  # immunoglobulin genes, T cell receptor genes, ribosome protein-coding genes,
  # heat shock proteins-associated genes, and mitochondrial genes.
  block_genes <- sn_get_signatures(
    species = "human",
    category =  c(
      "tcr",
      #T cell receptor genes
      "immunoglobulins",
      # immunoglobulin genes
      "ribo",
      # ribosome protein-coding genes
      "mito",
      # mitochondrial genes
      "heatshock" # heat shock proteins-associated genes
      # "noncoding",
      # "pseudogenes",
      # "g1s",
      # "g2m"
    )
  )
  sn_write(block_genes |> as.data.frame(),
           "data/metadata/block_genes.txt")
  
  g2m <- sn_get_signatures(species = "human", category = "g2m")
  sn_write(g2m |> as.data.frame(), "data/metadata/g2m_genes.txt")
  g1s <- sn_get_signatures(species = "human", category = "g1s")
  sn_write(g1s |> as.data.frame(), "data/metadata/g1s_genes.txt")
  
  gtf <- rtracklayer::import("data/metadata/genes.gtf")
  gene_type_df <- as.data.frame(gtf) |>
    select(gene_id, gene_name, gene_type) |> distinct()
  sn_write(gene_type_df, "data/metadata/gene_types.txt")
  gene_type_df <- sn_read("data/metadata/gene_types.txt")
  
  keep_genes <- gene_type_df |>
    filter(gene_type != "lncRNA") |>
    pull(gene_name)
  merged <- subset(merged, features = keep_genes)
  
  reference <- merged |>
    subset(scDblFinder.class == "singlet") |>
    subset(dataset_role == "Reference") |>
    subset(dataset_tier == "Core") |>
    subset(dataset_set == "Discovery")
  
  if (FALSE) {
    seurat_obj <- merged |>
      subset(sample == "Blood1")
    seurat_obj <- seurat_obj |>
      subset(scDblFinder.class == "singlet") |>
      subset(dataset_role == "Reference") |>
      subset(dataset_tier == "Core") |>
      subset(dataset_set == "Discovery")
    seurat_obj <- sn_remove_ambient_contamination(seurat_obj)
    seurat_obj <- sn_run_cluster(seurat_obj, layer = "decontaminated_counts")
    DimPlot(seurat_obj, label = TRUE, group.by = "seurat_clusters")
    DimPlot(seurat_obj, label = TRUE, group.by = "scDblFinder.class")
    DimPlot(seurat_obj, label = TRUE, group.by = "dataset_tier")
    DimPlot(seurat_obj, label = TRUE, group.by = "Phase")
    FeaturePlot(
      seurat_obj,
      features = "percent.mt",
      order = TRUE,
      max.cutoff = "q95"
    )
    FeaturePlot(
      seurat_obj,
      features = "",
      order = TRUE,
      max.cutoff = "q95"
    )
    seurat_obj <- sn_find_de(
      seurat_obj,
      analysis = "markers",
      min_pct = 0.3,
      logfc_threshold = 0.5
    )
    sn_plot_dot(seurat_obj,
                features = c("CD3D", "CD3E", "CD3G", "KLRF1", "GRIP1", "KIT"))
    sn_plot_feature(seurat_obj, features = "GRIP1")
    
  }
  
  keep_samples <- reference@meta.data |>
    count(sample) |>
    filter(n > 100) |>
    pull(sample)
  reference <- reference |> subset(sample %in% keep_samples)
  # reference <- sn_filter_genes(reference, min_cells = 10)
  
  counts <- LayerData(reference, layer = "counts", assay = "RNA")
  BPCells::write_matrix_anndata_hdf5(counts, "data/raw/reference/counts.h5") # 81.1 GB ->
  adt <- LayerData(reference, layer = "counts", assay = "ADT")
  BPCells::write_matrix_anndata_hdf5(adt, "data/raw/reference/adt.h5")
  
  metadata <- reference@meta.data |>
    rownames_to_column("barcode")
  sn_write(metadata, "data/raw/reference/metadata.csv")
}