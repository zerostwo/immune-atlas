library(tidyverse)
library(Seurat)
library(glue)
library(Shennong)

# files <- list.files(
#     "data/processed/EMadissoon2020/corrected/bpcells/",
#     pattern = "qs",
#     full.names = TRUE
# )
# # files <- files[-3]
# EMadissoon2020 <- map(files, sn_read)
# EMadissoon2020 <- merge(x = EMadissoon2020[[1]],
#                         y = EMadissoon2020[-1],
#                         add.cell.ids = str_remove_all(basename(files), ".qs"))
# EMadissoon2020 <- JoinLayers(EMadissoon2020, layers = c("counts", "decontaminated_counts"))
# EMadissoon2020$study <- "EMadissoon2020"
# EMadissoon2020

files <- list.files(
    "data/processed/DanielPCaron2025/corrected/bpcells",
    pattern = "qs",
    full.names = TRUE
)
# files <- files[-3]
DanielPCaron2025 <- map(files, sn_read)
DanielPCaron2025 <- merge(x = DanielPCaron2025[[1]],
                          y = DanielPCaron2025[-1],
                          add.cell.ids = str_remove_all(basename(files), ".qs"))
DanielPCaron2025 <- JoinLayers(DanielPCaron2025, layers = c("counts", "decontaminated_counts"))
DanielPCaron2025$study <- "DanielPCaron2025"
DanielPCaron2025

files <- list.files("data/processed/ShuaiHe2020/corrected/bpcells",
                    pattern = "qs",
                    full.names = TRUE)
ShuaiHe2020 <- map(files, sn_read)
ShuaiHe2020 <- merge(x = ShuaiHe2020[[1]], y = ShuaiHe2020[-1],
                     add.cell.ids = str_remove_all(basename(files), ".qs"))
ShuaiHe2020 <- JoinLayers(ShuaiHe2020,
                          layers = c("counts", "decontaminated_counts"))
ShuaiHe2020$study <- "ShuaiHe2020"

files <- list.files(
    "data/processed/MarylineFalquet2023/bpcells",
    pattern = "qs",
    full.names = TRUE
)
MarylineFalquet2023 <- map(files, sn_read)
MarylineFalquet2023 <- merge(x = MarylineFalquet2023[[1]], y = MarylineFalquet2023[-1], add.cell.ids = str_remove_all(basename(files), ".qs"))
MarylineFalquet2023 <- JoinLayers(MarylineFalquet2023, layers = c("counts", "decontaminated_counts"))
MarylineFalquet2023$study <- "MarylineFalquet2023"

files <- list.files(
    "data/processed/AndrewDHildreth2021/bpcells",
    pattern = "qs",
    full.names = TRUE
)
AndrewDHildreth2021 <- map(files, sn_read)
AndrewDHildreth2021 <- merge(x = AndrewDHildreth2021[[1]], y = AndrewDHildreth2021[-1], add.cell.ids = str_remove_all(basename(files), ".qs"))
AndrewDHildreth2021 <- JoinLayers(AndrewDHildreth2021, layers = c("counts", "decontaminated_counts"))
AndrewDHildreth2021$study <- "AndrewDHildreth2021"

files <- list.files(
    "data/processed/CDominguezConde2022/bpcells",
    pattern = "qs",
    full.names = TRUE
)
CDominguezConde2022 <- map(files, sn_read)
CDominguezConde2022 <- merge(x = CDominguezConde2022[[1]], y = CDominguezConde2022[-1],
                             add.cell.ids = str_remove_all(basename(files), ".qs"))
CDominguezConde2022 <- JoinLayers(CDominguezConde2022, layers = c("counts", "decontaminated_counts"))
CDominguezConde2022$study <- "CDominguezConde2022"

files <- list.files(
    "data/processed/NataliaJaeger2024/bpcells",
    pattern = "qs",
    full.names = TRUE
)
NataliaJaeger2024 <- map(files, sn_read)
NataliaJaeger2024 <- merge(x = NataliaJaeger2024[[1]], y = NataliaJaeger2024[-1], add.cell.ids = str_remove_all(basename(files), ".qs"))
NataliaJaeger2024 <- JoinLayers(NataliaJaeger2024, layers = c("counts", "decontaminated_counts"))
NataliaJaeger2024$study <- "NataliaJaeger2024"

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
merged
table(merged$study)

prefix <- "AndrewDHildreth2021_CDominguezConde2022_NataliaJaeger2024_MarylineFalquet2023_ShuaiHe2020_DanielPCaron2025_EMadissoon2020"
merged_path <- glue("data/processed/{prefix}.merged.qs")

LayerData(merged, "data") <- NULL
LayerData(merged, "scale.data") <- NULL

merged <- JoinLayers(merged, assay = "ADT")

sn_write(merged, path = merged_path)
sample_metadata <- merged@meta.data |>
    select(study, orig.ident) |>
    distinct()
sample_metadata_path <- glue("data/processed/{prefix}.sample_metadata.qs")
sn_write(sample_metadata, path = sample_metadata_path)

merged <- sn_read(path = merged_path)
merged

# Merge -------------------------------------------------------------------
merged1 <- sn_read(
    "data/processed/AndrewDHildreth2021_CDominguezConde2022_NataliaJaeger2024_MarylineFalquet2023_ShuaiHe2020_DanielPCaron2025.merged.qs"
)
merged2 <- sn_read("data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens.merged.qs")
counts <- sn_read(
    "data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens2022/bpcells/counts"
)
LayerData(merged2, layer = "counts") <- counts
adt_counts <- sn_read(
    "data/processed/ZhenlongLi2024_StevenBWells2025_TabulaSapiens2022/bpcells/adt_counts/"
)
LayerData(merged2, layer = "counts", assay = "ADT") <- adt_counts

merged <- merge(x = merged1, y = merged2)
merged <- JoinLayers(merged)
merged <- JoinLayers(merged, assay = "ADT")
merged <- JoinLayers(merged, assay = "HTO")


sn_write(merged, "data/processed/20260108.reference.raw.qs")


table(merged$scDblFinder.class, merged$study) |>
    as.data.frame() |>
    ggplot(aes(x = Freq, y = Var2, fill = Var1)) +
    geom_col(position = "fill") +
    catplot::theme_cat(show_title = "x", aspect_ratio = 1) +
    scale_fill_brewer(palette = "Paired") +
    scale_x_continuous(expand = c(0, 0))

# QC ----------------------------------------------------------------------
merged <- sn_read("data/processed/20260108.reference.raw.qs")
glimpse(merged@meta.data)
quantile(merged$percent.mt)
quantile(merged$nFeature_RNA)
quantile(merged$nCount_RNA)
table(merged$study)
table(merged$scDblFinder.class)

merged@meta.data <- merged@meta.data |>
    rownames_to_column("cell_id") |>
    mutate(
        prefix  = str_replace(cell_id, "_[ACGTN]+-[0-9]+$", ""),
        barcode = str_extract(cell_id, "[ACGTN]+-[0-9]+$")
    ) |>
    column_to_rownames("cell_id")
# mutate(sample = if_else(is.na(sample), orig.ident, sample))
# study, donor, age, gender, sample, tissue
table(merged$sample)
merged@meta.data <- merged@meta.data |>
    mutate(donor = str_split_fixed(orig.ident, "-", 2)[, 1],
           sample = if_else(is.na(sample), prefix, sample)) |>
    mutate(
        donor = case_when(
            study == "ShuaiHe2020" ~ "ShuaiHe2020",
            donor %in% c("SRR12423009", "SRR12423015") ~ "Lean_donor1",
            donor %in% c("SRR12423010", "SRR12423016") ~ "Lean_donor2",
            donor %in% c("SRR12423011", "SRR12423017") ~ "Lean_donor3",
            donor %in% c("SRR12423012", "SRR12423018") ~ "Obese_donor1",
            donor %in% c("SRR12423013", "SRR12423019") ~ "Obese_donor2",
            donor %in% c("SRR12423014", "SRR12423020") ~ "Obese_donor3",
            donor %in% c("SRR12437066", "SRR12437068") ~ "CD200r_lean",
            donor %in% c("SRR12437067", "SRR12437069") ~ "CD200r_obese",
            TRUE ~ donor
        ),
        study = if_else(study == "TabulaSapiens", "TabulaSapiens2022", study)
    )
# study, donor, age, gender
donor_metadata <- merged@meta.data |>
    select(study, donor) |>
    distinct() |>
    as_tibble()
sn_write(donor_metadata,
         "data/metadata/20251229.reference.donor_metadata.csv")

# study, donor, sample, tissue
sample_metadata <- merged@meta.data |>
    select(study, donor, sample) |>
    distinct() |>
    mutate(tissue = str_split_fixed(sample, "_", 3)[, 2]) |>
    mutate(tissue = if_else(tissue == "", str_split_fixed(sample, "-", 3)[, 2], tissue)) |>
    mutate(
        tissue = case_when(
            tissue %in% c("BoneMarrow", "BMA", "BM") ~ "Bone marrow",
            tissue %in% c("LIV") ~ "Liver",
            tissue %in% c("SPL") ~ "Spleen",
            tissue %in% c("LNG") ~ "Lung",
            tissue %in% c("SKN") ~ "Skin",
            tissue %in% c("BLO", "BLD") ~ "Blood",
            tissue %in% c("LymphNodes", "LymphNode") ~ "Lymph node",
            tissue %in% c("THY") ~ "Thymus",
            tissue %in% c("OME") ~ "Omentum",
            tissue %in% c("SKM") ~ "Skeletal muscle",
            tissue %in% c("MLN") ~ "Mesenteric lymph node",
            tissue %in% c("LLN") ~ "Lung-associated lymph node",
            tissue %in% c("TCL") ~ "Transverse colon",
            tissue %in% c("DUO") ~ "Duodenum",
            tissue %in% c("SCL") ~ "Sigmoid colon",
            tissue %in% c("CAE") ~ "Caecum",
            tissue %in% c("SI") ~ "Small intestine",
            tissue %in% c("FAT") ~ "Fat",
            tissue %in% c("ILN") ~ "Inguinal lymph node",
            tissue %in% c("BAL") ~ "Bronchoalveolar lavage",
            tissue %in% c("SalivaryGland") ~ "Salivary gland",
            tissue %in% c("JEJEPI") ~ "Jejunum epithelium",
            tissue %in% "JEJLP" ~ "Jejunum lamina propria",
            tissue %in% c("COLEPI") ~ "Colon epithelium",
            tissue %in% c("COLLP") ~ "Colon lamina propria",
            tissue %in% c("ILE", "TIL") ~ "Ileum",
            tissue %in% c("TLN") ~ "Thoracic lymph node",
            study %in% "NataliaJaeger2024" ~ "",
            TRUE ~ tissue
        )
    ) |>
    as_tibble()
table(sample_metadata$tissue) |> sort()
sn_write(sample_metadata,
         "data/metadata/20251229.reference.sample_metadata.csv")

# m <- Azimuth:::LoadH5ADobs("/mnt/public/CDominguezConde2022/download/CDominguezConde2022.h5ad")
ia_sample_metadata <- sn_read("data/metadata/IA_sample_spreadsheet.xlsx", sheet = 2)
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
    mutate(
        assay = case_when(
            study == "ZhenlongLi2024"  ~ "10x3'v3.1",
            study == "ShuaiHe2020" ~ "10x5'v1",
            study == "MarylineFalquet2023" ~ "10x3'v3"
            TRUE ~ assay
        )
    )

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
    mutate(age = as.numeric(age), BMI = as.numeric(BMI))
donor_metadata |>
    ggplot(aes(x = BMI)) +
    geom_histogram(bins = 14, color = "white") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 14))

donor_metadata |>
    ggplot(aes(x = age)) +
    geom_histogram(bins = 14, color = "white") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 14))
# scale_x_continuous(expand = c(0,0), limits = c(10, 90))

donor_metadata |>
    count(gender) |>
    ggplot(aes(x = gender, y = n)) +
    geom_col()
# scale_y_continuous(expand = c(0,0), limits = c(0, 14)) +
# scale_x_continuous(expand = c(0,0), limits = c(10, 90))


donor_metadata |>
    count(ethnicity) |>
    ggplot(aes(x = ethnicity, y = n)) +
    geom_col() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +
    catplot::theme_cat(x_text_angle = 45)
