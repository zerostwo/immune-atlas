if (FALSE) {
  library(tidyverse)
  library(Seurat)
  library(qs)
  library(glue)
  library(Shennong)
  
  outdir <- sn_set_path("/home/sduan/projects/immune-atlas/data/processed/scrup/batches")
  
  files <- list.files("/home/sduan/projects/immune-atlas/data/processed/scrup/bpcells", pattern = "qs", full.names = TRUE)
  head(files)
  
  batches <- split(files, cut(seq_along(files), 50, labels = FALSE))
  
  walk2(batches, seq_along(batches), function(batch_files, i) {
    qs_file <- glue("{outdir}/batch_{i}.qs")
    if (!file.exists(qs_file)) {
      print(i)
      # print(batch_files)
      # i <- 48
      batch_files <- batches[[i]]
      seurat_obj_lists <- map(batch_files, function(x){
        s <- sn_read(x)
        # LayerData(s, "spliced") <- NULL
        # LayerData(s, "unspliced") <- NULL
        # LayerData(s, "ambiguous") <- NULL
        s
      })
      
      merged <- merge(x = seurat_obj_lists[[1]],
                      y = seurat_obj_lists[-1])
      # merged <- JoinLayers(object = merged,
      #                      layers = c("counts","spliced", "unspliced", "ambiguous"))
      merged <- JoinLayers(object = merged)
      
      sn_write(merged, path = qs_file)
      
      message("Saved batch ", i, " with ", length(batch_files), " samples")
    }
  })
}

if (FALSE) {
  library(tidyverse)
  library(Seurat)
  library(glue)
  library(Shennong)
  
  # batch_lists <- list.files("/home/sduan/projects/immune-atlas/data/processed/scrup/batches/", full.names = TRUE,
  #                           pattern = "qs")
  # bpcells_root <- sn_set_path("/home/sduan/projects/immune-atlas/data/processed/scrup")
  # 
  # seurat_obj_lists <- map(
  #   set_names(batch_lists, nm = tools::file_path_sans_ext(basename(batch_lists))),
  #   function(fp) {
  #     message("Processing: ", fp)
  #     s <- sn_read(fp)
  #     
  #     # drop_layers <- intersect(c("spliced", "unspliced", "ambiguous"), Layers(s))
  #     # if (length(drop_layers)) {
  #     #   for (ln in drop_layers) LayerData(s, ln) <- NULL
  #     # }
  #     
  #     counts <- LayerData(s, "counts")
  #     bp_path <- glue("{bpcells_root}/{tools::file_path_sans_ext(basename(fp))}")
  #     
  #     if (dir_exists(bp_path)) {
  #       message("Found existing BPCells dir: ", bp_path, " — reusing")
  #     } else {
  #       message("Writing BPCells dir: ", bp_path)
  #       BPCells::write_matrix_dir(counts, bp_path)
  #     }
  #     
  #     LayerData(s, "counts") <- BPCells::open_matrix_dir(bp_path)
  #     qs::qsave(s, glue("{bpcells_root}/{tools::file_path_sans_ext(basename(fp))}.qs"))
  #     s
  #   }
  # )
  # 
  seurat_obj_lists <- list.files("/home/sduan/projects/immune-atlas/data/processed/scrup/batches/",pattern = ".qs$",
             full.names = TRUE) %>% 
    map(sn_read)
  merged <- merge(seurat_obj_lists[[1]],
                      seurat_obj_lists[-1])
  merged <- JoinLayers(merged)
  sn_write(merged, "/home/sduan/projects/immune-atlas/data/processed/scrup/merged.qs")
}

if (FALSE) {
  library(tidyverse)
  library(Seurat)
  library(glue)
  library(Shennong)
  
  merged <- sn_read("/home/sduan/projects/immune-atlas/data/processed/scrup/merged.qs")
  merged <- merged |> 
    subset(scDblFinder.class == "singlet")
  # Cells identified as low-quality or germ line in the original studies were
  # excluded, and only cells meeting the following criteria were retained:
  # 500–8,000 genes, 1,000–100,000 gene counts, and less than 20% mitochondrial
  # gene counts
  filtered_merged <- merged |> 
    subset(nCount_RNA > 1000) |> 
    subset(nCount_RNA < 100000) |> 
    subset(nFeature_RNA < 8000) 
  # We then excluded samples with fewer than 50 high-quality cells.
  keep_samples <- filtered_merged@meta.data |> 
    count(orig.ident) |> 
    filter(n > 100) |> 
    pull(orig.ident)
  filtered_merged <- filtered_merged |> 
    subset(orig.ident %in% keep_samples)
  sn_write(filtered_merged, "/home/sduan/projects/immune-atlas/data/processed/scrup/stringent_qc.merged.qs")
  filtered_merged <- sn_read("/home/sduan/projects/immune-atlas/data/processed/scrup/stringent_qc.merged.qs")
  library(googlesheets4)
  ss <- "1iK1LJDQCukjIOPlhrTbLDEgh_xKBAYxetuOxXGbPaQM"
  sample_metadata <- read_sheet(ss = ss, sheet = "final")
  sample_metadata <- sample_metadata |> 
    filter(sample != "SRX24386715") |> # Mixed
    filter(sample %in% keep_samples)
  
  # Harmonize metadata ------------------------------------------------------
  # 1. Harmonize disease names
  sample_metadata <- sample_metadata |> 
    mutate(harmonized_disease = case_when(
      disease %in% c("Head and Neck Squamous Cell Carcinoma",
                     "Hypopharygeal Squamous Cell Carcinoma",
                     "Oral Squamous Cell Carcinoma",
                     "head and neck squamous cancer",
                     "esophagus squamous cell carcinoma",
                     "HPV+ positive head and neck cancer",
                     "head and neck squamous cell carcinoma",
                     "head and neck squamous carcinoma",
                     "laryngeal squamous cell carcinoma",
                     "Laryngeal Squamous Cell Carcinoma") ~ "HNSCC",
      disease %in% c("Papillary Thyroid Carcinoma",
                     "Anaplastic Thyroid Cancer",
                     "thyroid carcinoma",
                     "anaplastic thyroid cancer") ~ "THCA",
      disease %in% c("Lung Adenocarcinoma",
                     "Lung Squamous Cell Carcinoma",
                     "Multiple Primary Lung Cancer",
                     "Non-small Cell Lung Cancer",
                     "Lung Cancer",
                     "lung cancer",
                     "lung squamous cell carcinoma",
                     "non-small cell lung cancer",
                     "lung adenocarcinoma") ~ "LC",
      disease %in% c("Hepatocellular Cancer",
                     "liver cancer",
                     "hepatocellular carcinoma") ~ "HCC",
      disease %in% c("Gallbladder Carcinoma",
                     "Intrahepatic Cholangiocarcinoma",
                     "Biliary Tract Cancer",
                     "cholangiocarcinoma",
                     "intrahepatic cholangiocarcinoma",
                     "Distal Cholangiocarcinoma") ~ "BTC",
      disease %in% c("Neuroblastoma",
                     "neuroblastoma") ~ "NB",
      disease %in% c("Prostate Cancer",
                     "prostate cancer") ~ "PRAD",
      disease %in% c("") ~ "UCEC",
      disease %in% c("") ~ "MM",
      disease %in% c("Acute Myeloid Leukemia",
                     "acute myeloid leukemia") ~ "AML",
      disease %in% c("") ~ "ALL",
      disease %in% c("chronic lymphocytic leukemia") ~ "CLL",
      disease %in% c("chronic myelomonocytic leukemia") ~ "CML",
      disease %in% c("nasopharyngeal carcinoma") ~ "NPC",
      disease %in% c("Esophageal Squamous Cell Carcinoma",
                     "Esophageal Adenocarcinoma",
                     "esophageal adenocarcinoma",
                     "esophageal squamous cell carcinoma",
                     "esophageal cancer") ~ "ESCA",
      disease %in% c("Breast Cancer",
                     "Invasive Breast Carcinoma",
                     "Triple Negative Breast Cancer",
                     "Invasive Ductal Carcinoma",
                     "breast cancer",
                     "Breast Ductal Carcinoma in Situ") ~ "BRCA",
      disease %in% c("Gastro-esophageal Adenocarcinoma","Gastric Cancer",
                     "gastric cancer",
                     "Stomach Adenocarcinoma") ~ "GC",
      disease %in% c("Pancreatic Ductal Adenocarcinoma",
                     "pancreatic neuroendocrine tumor",
                     "Pancreatic Cancer","pancreatic cancer",
                     "pancreatic ductal adenocarcinoma") ~ "PACA",
      disease %in% c("Clear Cell Renal Cell Carcinoma",
                     "clear renal cell carcinoma",
                     "renal cell carcinoma",
                     "Renal Medullary Carcinoma",
                     "renal clear cell carcinoma",
                     "renal medullary carcinoma",
                     "clear cell renal cell carcinoma",
                     "Papillary Renal Cell Carcinoma") ~ "RC",
      disease %in% c("Colorectal Cancer",
                     "colorectal cancer") ~ "CRC",
      disease %in% c("Ovarian Carcinoma",
                     "ovarian cancer",
                     "ovarian adenocarcinoma",
                     "high-grade serous ovarian cancer",
                     "High-grade Serous Ovarian Carcinoma") ~ "OV",
      disease %in% c("") ~ "FTC",
      disease %in% c("Melanoma","Acral Melanoma",
                     "melanoma","uveal melanoma",
                     "Cutaneous Melanoma") ~ "MELA",
      disease %in% c("Basal Cell Carcinoma",
                     "basal cell carcinoma") ~ "BCC",
      disease %in% c("Cutaneous Squamous Cell Carcinoma",
                     "Squamous Cell Carcinoma in Situ",
                     "cutaneous squamous cell carcinoma",
                     "squamous cell carcinoma",
                     "Squamous Cell Carcinoma") ~ "cSCC",
      disease %in% c("Merkel Cell Carcinoma",
                     "Merkel cell carcinoma") ~ "MCC",
      disease %in% c("Glioblastoma", "Glioma",
                     "glioblastoma") ~ "GL",
      disease %in% c("Bladder Cancer",
                     "Muscle-invasive Urothelial Bladder Cancer") ~ "BLCA",
      disease %in% c("Cervical Cancer", "cervical cancer") ~ "CESC",
      disease %in% c("Medulloblastoma") ~ "MB",
      disease %in% c("Testicular Seminoma") ~ "TGCT",
      # disease %in% c("") ~ "OS",
      # disease %in% c("") ~ "NET",
      # disease %in% c("") ~ "LYM",
      # disease %in% c("") ~ "GCTB",
      # disease %in% c("") ~ "CHC",
      # Special cases
      # experiment_accession %in% c("SRX21311407","SRX21311409") ~ "PF",
      # experiment_accession == "SRX21311402" ~ "CD",
      # bioproject == "PRJNA1003807" & is.na(disease) ~ "Normal",
      # disease == "GSM5788501" ~ "cSCC",
      study == "SRP454227" ~ "Normal",
      TRUE ~ disease
    ))
  table(sample_metadata$harmonized_disease) |> sort()
  # 2. Harmonize tissue / sample_site names
  # csem_metadata$tissue <- csem_metadata$sample_site
  sample_metadata <- sample_metadata |> 
    mutate(tissue = ifelse(is.na(tissue), tissue_raw, tissue)) |>
    mutate(
      harmonized_tissue = case_when(
        tissue %in% c("") ~ "Heart",
        tissue %in% c("") ~ "Vasculature",
        tissue %in% c("Sigmoid",
                      "Sigmoid colon",
                      "Sigmoid colon",
                      "colon",
                      "Descending colon",
                      "Transverse colon",
                      "Ascending colon") ~ "Colon",
        tissue %in% c("bile duct") ~ "Bile duct",
        tissue %in% c("") ~ "Gall bladder",
        tissue %in% c("bone") ~ "Bone",
        tissue %in% c("ascites") ~ "Ascites",
        tissue %in% c("esophagus") ~ "Esophagus",
        tissue %in% c("liver") ~ "Liver",
        tissue %in% c("Oral cavity",
                      "Tongue/Tonsil",
                      "Saliva gland",
                      "Right lower gingiva",
                      "Laryngeal","Left nasopharynx",
                      "Left neck","Right mouth floor",
                      "Parotid gland",
                      "head and neck",
                      "parotid gland",
                      "unknown",
                      "Tonsil",
                      "head","nasopharynx",
                      "Right tongue") ~ "Head and neck",
        tissue %in% c("") ~ "Rectum",
        tissue %in% c("cervix") ~ "Cervix",
        tissue %in% c("kidney") ~ "Kidney",
        tissue %in% c("Ileum LPL","Ileum IEL") ~ "Small intestine",
        tissue %in% c("stomach") ~ "Stomach",
        tissue %in% c("") ~ "Fat",
        tissue %in% c("Pancreas",
                      "pancreas",
                      "Pancreas body/tail",
                      "Pancreas head") ~ "Pancreas",
        tissue %in% c("Iymph node",
                      "Pelvic lymph node",
                      "lymph node",
                      "Renal hilus lymph node") ~ "Lymph node",
        tissue %in% c("Thyroid", "thyroid") ~ "Thymus",
        tissue %in% c("") ~ "Skeletal muscle",
        tissue %in% c("Retina", "eye") ~ "Eye",
        tissue %in% c("breast") ~ "Breast",
        tissue %in% c("ovary") ~ "Ovary",
        tissue %in% c("prostate") ~ "Prostate",
        tissue %in% c("") ~ "Testis",
        tissue %in% c("") ~ "Uterus",
        tissue %in% c("vaginal apex") ~ "Vagina",
        tissue %in% c("lung") ~ "Lung",
        tissue %in% c("") ~ "Trachea",
        tissue %in% c("Right arm",
                      "Right chest",
                      "skin",
                      "Left cheek",
                      "Nose","Left arm") ~ "Skin",
        tissue %in% c("") ~ "Bladder",
        tissue %in% c("") ~ "Ureter",
        tissue %in% c("Brain left cerebellum",
                      "Brain left occipital lobe",
                      "Brain right temporal lobe",
                      "Brain right remporal lobe",
                      "Brain left frontal lobe",
                      "Brain right parietal lobe",
                      "Brain left temporal lobe",
                      "Cerebellum",
                      "Brain right frontal lobe",
                      "Brain cerebellum",
                      "Brain left parietooccipital region",
                      "Cerebellum and spine",
                      "Spine",
                      "brain",
                      "Brain",
                      "Splenial extension into brain left parietal lobe") ~ "CNS",
        tissue %in% c("bone marrow") ~ "Bone marrow",
        tissue %in% c("blood") ~ "Blood",
        tissue %in% c("adrenal gland") ~ "Adrenal gland",
        tissue %in% c("Peritoneum") ~ "Peritoneum",
        # tissue %in% c("mixed sample (Tonsil and blood)") ~ "Mixed",
        # Special cases
        # study_accession == "SRP421224" & tissue == "intestine" ~ "Colon",
        # study_accession %in% c("SRP362712","SRP415592") & tissue == "intestine" ~ "Esophagus",
        TRUE ~ tissue
      )
    )
  table(sample_metadata$harmonized_tissue)
  # 3. Harmonize sample_type names
  sample_metadata <- sample_metadata |> 
    mutate(
      harmonized_sample_type = case_when(
        sample_type %in% c("PBMC") ~ "Blood",
        sample_type %in% c("adjacent normal") ~ "Normal",
        sample_type %in% c("Tumour", "tumor tissue") ~ "Tumor",
        sample_type %in% c("metastatic tumor") ~ "Metastatic",
        sample_type == "other tissue from tumor patient" ~ harmonized_tissue,
        sample_type %in% c("CD", "fibrotic") ~ "Inflammation",
        # Special cases
        # bioproject == "PRJNA1003807" ~ "Normal",
        # experiment_accession %in% c("SRX21311407","SRX21311402",
        #                             "SRX21311409") ~ "Inflammation",
        TRUE ~ sample_type
      ))
  table(sample_metadata$harmonized_sample_type)
  
  # 4. Harmonize assay
  sample_metadata <- sample_metadata %>%
    mutate(
      assay = case_when(
        # 10x 3' chemistries
        strand == "Forward" & cb_length == 14 & umi_length == 10 & cb_whitelist == "737K-august-2016.txt" ~ "10x 3' v1",
        strand == "Forward" & cb_length == 16 & umi_length == 10 & cb_whitelist == "737K-august-2016.txt" ~ "10x 3' v2",
        strand == "Forward" & cb_length == 16 & umi_length == 12 & cb_whitelist %in% c("3M-february-2018.txt", "3M-february-2018_TRU.txt") ~ "10x 3' v3",
        strand == "Forward" & cb_length == 16 & umi_length == 12 & cb_whitelist %in% c("4M-with-alts-february-2019.txt", "4M-with-alts-february-2019_TRU.txt") ~ "10x 3' v4",
        
        # 10x 5' chemistries
        strand == "Reverse" & cb_length == 16 & umi_length == 10 & cb_whitelist == "737K-august-2016.txt" ~ "10x 5' v1",
        strand == "Reverse" & cb_length == 16 & umi_length == 12 & cb_whitelist %in% c("3M-february-2018.txt", "3M-february-2018_TRU.txt") ~ "10x 5' v2",
        strand == "Reverse" & cb_length == 16 & umi_length == 12 & cb_whitelist %in% c("4M-with-alts-february-2019.txt", "4M-with-alts-february-2019_TRU.txt") ~ "10x 5' v3",
        
        # ---- Multiome (ARC) ----
        cb_length == 16 & umi_length == 12 &
          cb_whitelist == "737K-arc-v1.txt" ~ "10x Multiome GEX (ARC v1)",
        cb_length == 16 & (is.na(umi_length) | umi_length == 0) &
          cb_whitelist == "737K-arc-v1.txt" ~ "10x Multiome ATAC (ARC v1)",
        
        # ---- Feature Barcode (ADT/HTO) ----
        # 若有 library_type 列：优先用它判断
        # !is.na(library_type) & library_type %in% c("Antibody Capture","Multiplexing Capture") ~ "10x Feature Barcode (ADT/HTO)",
        # 没有 library_type 时的经验规则（3M/4M 且 UMI=10 常见于 ADT/HTO）
        cb_length == 16 & umi_length == 10 &
          cb_whitelist %in% c("3M-february-2018.txt","3M-february-2018_TRU.txt",
                              "4M-with-alts-february-2019.txt","4M-with-alts-february-2019_TRU.txt") ~ "10x Feature Barcode (ADT/HTO)",
        
        # 其它未知或不在表中的
        TRUE ~ "unknown"
      )
    )
  table(sample_metadata$assay)
  
  sample_metadata |> 
    filter(harmonized_tissue != "Mixed") |> 
    group_by(study,harmonized_disease) |>
    count(harmonized_disease) |> 
    ungroup() |> 
    ggplot(aes(x = study, y = harmonized_disease,
               size = n))+
    geom_point() +
    catplot::theme_cat(show_panel_grid = "both",
                       show_title = "none",
                       x_text_angle = 90) +
    coord_fixed()
  
  sample_metadata |> 
    filter(harmonized_tissue != "Mixed") |> 
    group_by(study,harmonized_sample_type, harmonized_tissue) |>
    count(harmonized_tissue) |> 
    ungroup() |> 
    ggplot(aes(x = study, y = harmonized_tissue,
               size = n)) +
    geom_point() +
    catplot::theme_cat(show_panel_grid = "both",
                       show_title = "none",
                       x_text_angle = 90) +
    coord_fixed()
  
  metadata |> 
    filter(harmonized_tissue != "Mixed") |> 
    group_by(harmonized_sample_type) |>
    count(harmonized_disease) |> 
    ungroup() |> 
    ggplot(aes(x = harmonized_disease, y = harmonized_sample_type,
               size = n))+
    geom_point() +
    catplot::theme_cat(show_panel_grid = "both",
                       show_title = "none",
                       x_text_angle = 90) +
    coord_fixed()
  
  
  write_csv(sample_metadata, "data/processed/scrup/harmonized.sample_metadata.csv")
  sample_metadata <- sn_read("data/processed/scrup/harmonized.sample_metadata.csv")
  sample_metadata <- sample_metadata |> 
    filter(harmonized_sample_type != "Inflammation") |> 
    filter(assay %in% c("10x 3' v2", "10x 3' v3", "10x 5' v1")) |> 
    select(-c(counts, barcode, sahmi_counts)) 
  keep_samples <- sample_metadata |> 
    pull(sample)
  filtered_merged <- filtered_merged |> 
    subset(orig.ident %in% unique(keep_samples))
  
  filtered_merged@meta.data <- filtered_merged@meta.data |> 
    rownames_to_column("barcode") |> 
    left_join(sample_metadata, by = c("orig.ident"="sample")) |> 
    column_to_rownames("barcode")
  filtered_merged
  sn_write(filtered_merged, "/home/sduan/projects/immune-atlas/data/processed/scrup/harmonized.qs")
  filtered_merged <- sn_read("/home/sduan/projects/immune-atlas/data/processed/scrup/harmonized.qs")
  # Notably, HVG selection was performed after removing specific genes,
  # including immunoglobulin genes, T cell receptor genes, ribosome
  # protein-coding genes, heat shock proteins-associated genes, and
  # mitochondrial genes.
  block_genes <- sn_get_signatures(
    species = "human", category = c("tcr","immunoglobulins",
                                    "mito", "ribo","heatshock",
                                    "noncoding", "pseudogenes")
  )
  filtered_merged <- filtered_merged[!(rownames(filtered_merged) %in% block_genes),]
  filtered_merged
  filtered_merged <- sn_filter_genes(filtered_merged, min_cells = 100)
  # filtered_merged <- NormalizeData(filtered_merged)
  # filtered_merged <- sn_score_cell_cycle(filtered_merged, species = "human")
  gtf <- rtracklayer::import("/home/sduan/reference_genomes/10xgenomics/GRCh38/genes/genes.gtf.gz")
  gtf <- as.data.frame(gtf)
  protein_genes <- gtf |> 
    select(gene_type, gene_name) |> 
    filter(gene_type %in% c("protein_coding")) |> 
    distinct() |> 
    pull(gene_name)
  filtered_merged <- filtered_merged[rownames(filtered_merged) %in% protein_genes,]
  LayerData(filtered_merged, "counts") |> 
    BPCells::write_matrix_10x_hdf5("data/processed/scrup/harmonized.counts.h5")
  filtered_merged@meta.data <- filtered_merged@meta.data |> 
    select(-study.x) |> 
    rename(study = study.y) |> 
    mutate(batch = paste0(study,"_",assay,"_",harmonized_tissue))
  filtered_merged@meta.data |> 
    rownames_to_column("barcode") |> 
    write_csv("data/processed/scrup/harmonized.cell_metadata.csv")
  # # seurat_obj <- qread("/home/sduan/scrup/object/merged.qs")
  # sahmi_table <- sn_read("/home/sduan/data/scrup/object/sahmi.table.qs")
  # metadata <- sn_read("/home/sduan/data/scrup/object/metadata.qs")
  # metadata <- metadata %>% 
  #   select(-c(barcode, counts, sahmi_counts))
  # seurat_obj <- seurat_obj[,seurat_obj$orig.ident %in% metadata$sample]
  # seurat_obj@meta.data <- seurat_obj@meta.data %>% 
  #   rownames_to_column("barcode") %>% 
  #   rename(sample=orig.ident) %>% 
  #   right_join(metadata, by = c("sample")) %>% 
  #   column_to_rownames("barcode")
  # seurat_obj <- seurat_obj %>% 
  #   subset(harmonized_disease != "BCC")
  # sort(table(seurat_obj$harmonized_disease))
  # 
  # sahmi_table <- sahmi_table |> 
  #   filter(barcode %in% colnames(seurat_obj))
  # 
  # table(sahmi_table$kingdom)
  # # Bacteria Eukaryota   Viruses 
  # # 1089030     95160      4937 
  # bacteria_barcodes <- sahmi_table |> 
  #   filter(kingdom == "Bacteria") |> 
  #   pull(barcode) |> unique()
  # eukaryota_barcodes <- sahmi_table |> 
  #   filter(kingdom == "Eukaryota") |> 
  #   pull(barcode) |> unique()
  # viruses_barcodes <- sahmi_table |>
  #   filter(kingdom == "Viruses") |> 
  #   pull(barcode) |> unique()
  # 
  # seurat_obj@meta.data <- seurat_obj@meta.data |> 
  #   rownames_to_column("barcode") |>
  #   mutate(
  #     bacteria = if_else(
  #       barcode %in% bacteria_barcodes, "Detected", "Not detected"
  #     ),
  #     eukaryota = if_else(
  #       barcode %in% eukaryota_barcodes, "Detected", "Not detected"
  #     ),
  #     viruses = if_else(
  #       barcode %in% viruses_barcodes, "Detected", "Not detected"
  #     )
  #   ) |> 
  #   column_to_rownames("barcode")
  # 
  # table(seurat_obj$bacteria)
  # table(seurat_obj$viruses)
  # table(seurat_obj$eukaryota)
  # 
  # seurat_obj@misc$sahmi_table <- sahmi_table
  # 
  # gtf <- rtracklayer::import("/home/sduan/scrup/reference/GRCh38/v48/sources/filtered.gtf")
  # gtf <- as.data.frame(gtf)
  # gene_type <- gtf %>% select(seqnames,gene_name, gene_type) %>% distinct() %>% arrange(gene_type)
  # keep_genes <- gene_type %>% 
  #   filter(gene_type != "lncRNA") %>% 
  #   pull(gene_name)
  # # a <- rownames(seurat_obj)[!(rownames(seurat_obj) %in% gene_type$gene_name)] %>% str_remove_all("\\.1") %>% str_remove_all("\\.2") 
  # seurat_obj <- seurat_obj[rownames(seurat_obj) %in% keep_genes,]
  # 
  # m <- seurat_obj@meta.data %>% 
  #   mutate(id = paste0(study,"_", harmonized_disease))
  # # need2remove_ids <- c(
  # #   "SRP439674_LC", "SRP327758_ESCA",
  # #   "SRP257868_PRAD", "SRP417464_PRAD",
  # #   "SRP458159_PACA", "SRP373187_PACA",
  # #   "SRP362712_ESCA", "SRP479490_LC",
  # #   "SRP330210_RC","SRP359815_PACA",
  # #   "SRP425933_PACA", "SRP435639_PACA",
  # #   "SRP443229_LC", "SRP396819_PACA",
  # #   "SRP277925_PACA", "SRP388760_RC",
  # #   "SRP374837_RC", "SRP389538_PRAD",
  # #   "SRP400172_PACA", "SRP348120_GL", "SRP370877_GL",
  # #   "SRP359928_ESCA", "SRP272677_PACA", "SRP376216_PACA",
  # #   "SRP396471_HNSCC", "SRP460816_HNSCC",
  # #   "SRP341560_HNSCC", "SRP355223_HNSCC"
  # # )
  # keep_cells <- m %>% 
  #   filter(!(id %in% need2remove_ids)) %>% 
  #   rownames()
  # 
  # pp <- m %>% 
  #   mutate(id = paste0(study,"_", harmonized_disease)) %>% 
  #   filter(harmonized_sample_type == "Tumor") %>% 
  #   count(id, bacteria) %>% 
  #   pivot_wider(names_from = "bacteria", values_from = "n", values_fill = 0) %>% 
  #   mutate(prop = round(Detected/(Detected+`Not detected`)*100,2))
  # 
  # seurat_obj <- seurat_obj[, colnames(seurat_obj) %in% keep_cells]
  # seurat_obj <- seurat_obj %>% 
  #   subset(scDblFinder.class=="singlet")
  # seurat_obj@misc$sahmi_table <- seurat_obj@misc$sahmi_table[seurat_obj@misc$sahmi_table$barcode %in% colnames(seurat_obj),]
  # seurat_obj$batch <- paste0(seurat_obj$study, "_", seurat_obj$assay)
  # qs::qsave(seurat_obj, "/home/sduan/projects/sahmi/data/20250917.pan-cancer.seurat_obj.qs")
  # 
  # seurat_obj <- qs::qread("/home/sduan/projects/sahmi/data/20250917.pan-cancer.seurat_obj.qs")
  # seurat_obj <- NormalizeData(seurat_obj)
  # 
  # library(SignatuR)
  # data(SignatuR)
  # g1s_genes <- GetSignature(SignatuR$Hs$Programs$cellCycle.G1S)[[1]]
  # write_delim(g1s_genes %>% as.data.frame(), "data/g1s_genes.txt",
  #             col_names = FALSE)
  # g2m_genes <- GetSignature(SignatuR$Hs$Programs$cellCycle.G2M)[[1]]
  # write_delim(g2m_genes %>% as.data.frame(), "data/g2m_genes.txt",
  #             col_names = FALSE)
  # mito_genes <- GetSignature(SignatuR$Hs$Compartments$Mito)[[1]]
  # ribo_genes <- GetSignature(SignatuR$Hs$Compartments$Ribo)[[1]]
  # tcr_genes <- GetSignature(SignatuR$Hs$Compartments$TCR)[[1]]
  # immunoglobulins_genes <- GetSignature(SignatuR$Hs$Compartments$Immunoglobulins)[[1]]
  # heatshock_genes <- GetSignature(SignatuR$Hs$Programs$HeatShock)[[1]]
  # exclude_genes <- c(mito_genes,
  #                    ribo_genes,
  #                    tcr_genes,
  #                    immunoglobulins_genes,
  #                    heatshock_genes)
  # write_delim(exclude_genes %>% as.data.frame(), "data/exclude_genes.txt",
  #             col_names = FALSE)
  # seurat_obj <- CellCycleScoring(seurat_obj, s.features = g1s_genes, 
  #                                g2m.features = g2m_genes)
  # 
  # counts <- LayerData(seurat_obj,"counts")
  # BPCells::write_matrix_10x_hdf5(counts, "/home/sduan/projects/sahmi/data/20250917.pan-cancer.counts.h5")
  # metadata <- seurat_obj@meta.data %>% 
  #   rownames_to_column("barcode")
  # write_csv(metadata,"/home/sduan/projects/sahmi/data/20250917.pan-cancer.metadata.csv")

}


