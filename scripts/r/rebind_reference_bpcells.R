#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BPCells)
  library(SeuratObject)
  library(Shennong)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/r/rebind_reference_bpcells.R <path/to/reference.raw.qs> <path/to/bpcells_dir> [path/to/output.qs]"
    ),
    call. = FALSE
  )
}

reference_qs <- normalizePath(args[[1]], mustWork = TRUE)
bpcells_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_qs <- if (length(args) >= 3) {
  normalizePath(args[[3]], mustWork = FALSE)
} else {
  sub("\\.qs$", ".rebound.qs", reference_qs)
}

counts_dir <- file.path(bpcells_dir, "counts")
decont_dir <- file.path(bpcells_dir, "decontaminated_counts")
adt_dir <- file.path(bpcells_dir, "adt")

required_dirs <- c(counts_dir, decont_dir, adt_dir)
missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_dirs) > 0) {
  stop(
    paste(
      "Missing BPCells directories:",
      paste(missing_dirs, collapse = ", ")
    ),
    call. = FALSE
  )
}

ref <- sn_read(reference_qs)

LayerData(ref, assay = "RNA", layer = "counts") <- BPCells::open_matrix_dir(counts_dir)
LayerData(ref, assay = "RNA", layer = "decontaminated_counts") <- BPCells::open_matrix_dir(decont_dir)
LayerData(ref, assay = "ADT", layer = "counts") <- BPCells::open_matrix_dir(adt_dir)

sn_write(ref, output_qs)

message("Rebound reference saved to: ", output_qs)
