#!/usr/bin/env Rscript
# 
# Download all raw HGMM pools from the 10x website and save one .h5ad file with raw counts per pool

suppressPackageStartupMessages({
  library(DropletUtils)
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
  library(yaml)
})

DATASETS <- list(
  hgmm10k_3v4 = "https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_hgmm_3p_gemx_Multiplex/10k_hgmm_3p_gemx_Multiplex_count_raw_feature_bc_matrix.h5",
  hgmm10k_5v3 = "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/10k_hgmm_5p_gemx_Multiplex/10k_hgmm_5p_gemx_Multiplex_count_raw_feature_bc_matrix.h5",
  hgmm5k_3v3  = "https://cf.10xgenomics.com/samples/cell-exp/8.0.0/5k_hgmm_3p_nextgem_Multiplex/5k_hgmm_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5",
  hgmm5k_5v2  = "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/5k_hgmm_5p_nextgem_Multiplex/5k_hgmm_5p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"
)

# arg parsing
source("src/common/cli.R")
p <- arg_parser("data_download module (raw, discovers its own pools)")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_argument(p, "--mixed_species", flag = TRUE, help = "whether this dataset mixes human and mouse cells")
p <- add_argument(p, "--batch_var", type = "character", help = "batch column name")
p <- add_argument(p, "--sample_var", type = "character", help = "sample column name")
p <- add_argument(p, "--labels_var", type = "character", help = "cell type / ground-truth labels column name")
args <- parse_args(p, argv = strip_output_args())

# handle null values
for (k in c("batch_var", "sample_var", "labels_var")) {
  if (!is.null(args[[k]]) && tolower(args[[k]]) %in% c("none", "null", "")) {
    args[[k]] <- NULL
  }
}

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  for (pool_name in names(DATASETS)) {
    url <- DATASETS[[pool_name]]
    cat(sprintf("  --- pool: %s ---\n", pool_name))

    pool_dir <- file.path(args$output_dir, pool_name)
    dir.create(pool_dir, showWarnings = FALSE, recursive = TRUE)

    h5_path <- file.path(pool_dir, "raw_download.h5")
    cat(sprintf("  downloading %s\n", url))
    download.file(url, h5_path, mode = "wb", quiet = TRUE)
    cat(sprintf("  wrote: %s (%.1f MB)\n", h5_path, file.size(h5_path) / 1e6))

    cat("  reading with read10xCounts ..\n")
    sce <- read10xCounts(h5_path, col.names = TRUE)
    counts(sce) <- as(counts(sce), "dgCMatrix")
    colData(sce)$pool_id <- pool_name
    cat(sprintf("  %s: %d barcodes x %d genes\n", pool_name, ncol(sce), nrow(sce)))
    file.remove(h5_path)

    out_path <- file.path(pool_dir, "raw.h5ad")
    # x_mapping: without it, anndataR leaves X empty and puts the assay in
    # layers["counts"] instead -- invisible to any plain scanpy.read_h5ad() consumer.
    write_h5ad(sce, out_path, x_mapping = "counts")
    cat(sprintf("  wrote: %s\n", out_path))
  }

  message("  writing properties.yaml ..")
  properties_path <- file.path(args$output_dir, paste0(args$name, "_properties.yaml"))
  write_yaml(
    list(
      mixed_species = isTRUE(args$mixed_species),
      batch_var     = args$batch_var,
      sample_var    = args$sample_var,
      labels_var    = args$labels_var
    ),
    properties_path
  )
  cat(sprintf("  wrote: %s\n", properties_path))
}

if (sys.nframe() == 0L) {
  main()
}
