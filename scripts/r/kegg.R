suppressPackageStartupMessages({
    library(httr)
    library(readr)
    library(stringr)
    library(dplyr)
    library(purrr)
    library(tidyr)
})

#---------------------------
# Helpers
#---------------------------
kegg_get_text <- function(endpoint, query = NULL) {
    base <- "https://rest.kegg.jp"
    url  <- paste0(base, "/", endpoint, if (!is.null(query)) paste0("/", query) else "")
    resp <- httr::GET(url)
    httr::stop_for_status(resp)
    httr::content(resp, as = "text", encoding = "UTF-8")
}

# Parse KEGG "list" / "link" style: each line "key\tvalue"
parse_kegg_kv <- function(txt, col1 = "from", col2 = "to") {
    readr::read_tsv(I(txt), col_names = c(col1, col2), show_col_types = FALSE, progress = FALSE)
}

# Extract classification (level1/2/3) from br08901 (KEGG Pathway hierarchy)
# br08901 structure uses letter codes A/B/C/D etc; D lines contain pathway map IDs like "map00010"
parse_br08901_pathway_hierarchy <- function(txt) {
    lines <- str_split(txt, "\n", simplify = TRUE)
    lines <- lines[nzchar(lines)]
    
    lvl1 <- NA_character_
    lvl2 <- NA_character_
    lvl3 <- NA_character_
    
    out <- vector("list", length(lines))
    j <- 0L
    
    for (i in seq_along(lines)) {
        ln <- lines[[i]]
        # Typical patterns:
        # A <Level1>
        # B  <Level2>
        # C   <Level3>
        # D    map00010  Glycolysis / Gluconeogenesis
        if (str_detect(ln, "^A\\s")) {
            lvl1 <- str_trim(str_remove(ln, "^A\\s+"))
            lvl2 <- NA_character_; lvl3 <- NA_character_
        } else if (str_detect(ln, "^B\\s")) {
            lvl2 <- str_trim(str_remove(ln, "^B\\s+"))
            lvl3 <- NA_character_
        } else if (str_detect(ln, "^C\\s")) {
            lvl3 <- str_trim(str_remove(ln, "^C\\s+"))
        } else if (str_detect(ln, "^D\\s")) {
            j <- j + 1L
            # Extract map id and name
            # Example: "D    map00010  Glycolysis / Gluconeogenesis"
            m <- str_match(ln, "^D\\s+([^\\s]+)\\s+(.+)$")
            map_id <- m[,2]
            map_name <- m[,3]
            out[[j]] <- tibble(
                map_id = map_id,
                map_name = map_name,
                class_level1 = lvl1,
                class_level2 = lvl2,
                class_level3 = lvl3
            )
        }
    }
    bind_rows(out[seq_len(j)])
}

# Batch KEGG get pathway entries (to reduce rate-limits)
kegg_get_pathway_entries <- function(pathway_ids, batch_size = 10, sleep_sec = 0.2) {
    pathway_ids <- unique(pathway_ids)
    batches <- split(pathway_ids, ceiling(seq_along(pathway_ids) / batch_size))
    
    res <- map_dfr(batches, function(b) {
        txt <- kegg_get_text("get", paste(b, collapse = "+"))
        Sys.sleep(sleep_sec)
        tibble(raw = str_split(txt, "\n///\n", simplify = FALSE)[[1]])
    })
    
    res <- res %>%
        mutate(raw = str_trim(raw)) %>%
        filter(nzchar(raw))
    
    res
}

# From a single KEGG pathway entry text, extract genes (hsa:xxxx)
extract_genes_from_entry <- function(entry_txt) {
    # GENE lines look like:
    # GENE        3098  HK1; hexokinase 1 [KO:K00844]
    #             3101  HK2; ...
    # We'll capture: EntrezID and Symbol
    gene_lines <- str_split(entry_txt, "\n", simplify = TRUE)
    gene_lines <- gene_lines[str_detect(gene_lines, "^GENE\\s|^\\s{12}\\d+\\s")]
    if (length(gene_lines) == 0) return(tibble(entrez_id = character(), symbol = character()))
    
    gene_lines <- str_replace(gene_lines, "^GENE\\s+", "")
    gene_lines <- str_trim(gene_lines)
    
    # Extract leading numeric ID and symbol before ';'
    m <- str_match(gene_lines, "^(\\d+)\\s+([^;\\s]+)")
    tibble(
        entrez_id = m[,2],
        symbol = m[,3]
    ) %>% filter(!is.na(entrez_id), !is.na(symbol))
}

#---------------------------
# 1) Get all human pathways (hsa)
#---------------------------
path_list_txt <- kegg_get_text("list", "pathway/hsa")
pathway_tbl <- parse_kegg_kv(path_list_txt, col1 = "pathway_id", col2 = "pathway_name") %>%
    mutate(
        pathway_id = str_replace(pathway_id, "^path:", ""),  # e.g. hsa00010
        pathway_name = str_trim(pathway_name)
    )

#---------------------------
# 2) Build pathway->classification map using br08901 hierarchy
#---------------------------
br_txt <- kegg_get_text("get", "hsa01100")
hier_tbl <- parse_br08901_pathway_hierarchy(br_txt)

# Map "map00010" -> "hsa00010"
hier_tbl <- hier_tbl %>%
    mutate(pathway_id = str_replace(map_id, "^map", "hsa"))

pathway_tbl <- pathway_tbl %>%
    left_join(hier_tbl %>% select(pathway_id, class_level1, class_level2, class_level3),
              by = "pathway_id")

# Note: a small number of human pathways may not be in br08901 or may map differently.
# They will have NA classification.

#---------------------------
# 3) Fetch genes for each pathway (KEGG get)
#---------------------------
entries_tbl <- kegg_get_pathway_entries(pathway_tbl$pathway_id, batch_size = 10, sleep_sec = 0.2)

# Extract ENTRY ID from each raw block to link back
entries_tbl <- entries_tbl %>%
    mutate(
        pathway_id = str_match(raw, "^ENTRY\\s+([^\\s]+)\\s+Pathway")[,2]
    ) %>% filter(!is.na(pathway_id))

pathway_gene_tbl <- entries_tbl %>%
    transmute(pathway_id, genes = map(raw, extract_genes_from_entry)) %>%
    unnest(genes) %>%
    distinct(pathway_id, entrez_id, symbol)

# Join pathway metadata + classification
pathway_gene_tbl <- pathway_gene_tbl %>%
    left_join(pathway_tbl, by = "pathway_id") %>%
    relocate(pathway_id, pathway_name, class_level1, class_level2, class_level3, entrez_id, symbol)

#---------------------------
# 4) Also provide a named list (pathway_id -> symbols)
#---------------------------
pathway_genes_list <- pathway_gene_tbl %>%
    group_by(pathway_id) %>%
    summarise(genes = list(unique(symbol)), .groups = "drop") %>%
    { setNames(.$genes, .$pathway_id) }

#---------------------------
# Outputs
#---------------------------
pathway_tbl
pathway_gene_tbl
pathway_genes_list

# Optionally write out:
# write_csv(pathway_tbl, "kegg_hsa_pathways_with_class.csv")
# write_csv(pathway_gene_tbl, "kegg_hsa_pathway_gene_long.csv")
# saveRDS(pathway_genes_list, "kegg_hsa_pathway_genes_list.rds")
