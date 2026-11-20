#!/usr/bin/env Rscript
# Download raw 10x pools for a dataset collection.
#
#`--dataset_name` selects which dataset to download

suppressPackageStartupMessages({
  library(DropletUtils)
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
  library(yaml)
  library(data.table)
  library(magrittr)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("raw data_download module")
p <- add_base_args(p) # --output_dir, --name
p <- add_argument(p, "--dataset_name", type = "character", help = "which dataset collection to download")
p <- add_argument(p, "--batch_var", type = "character", help = "batch column name")
p <- add_argument(p, "--sample_var", type = "character", help = "sample column name")
p <- add_argument(p, "--tenx_assay", type = "character", help = "10x assay type (e.g., polyA or flex)")
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

# Download a 10x raw h5 matrix into pool_dir and read it in as an SCE, with
# the rowData Type ("Gene Expression"/"Antibody Capture") renamed to feature_type. 
download_and_read_10x <- function(h5_url, pool_dir, h5_filename = "raw_download.h5") {
  h5_path <- file.path(pool_dir, h5_filename)
  cat(sprintf("  downloading %s\n", h5_url))
  download.file(h5_url, h5_path, mode = "wb", quiet = TRUE)
  cat(sprintf("  wrote: %s (%.1f MB)\n", h5_path, file.size(h5_path) / 1e6))

  cat("  reading with read10xCounts ..\n")
  sce <- read10xCounts(h5_path, col.names = TRUE)
  counts(sce) <- as(counts(sce), "dgCMatrix")

  colnames(rowData(sce))[colnames(rowData(sce)) == "Type"] <- "feature_type"

  list(sce = sce, h5_path = h5_path)
}

# Download a cellranger multi config CSV and read just its [samples] section
download_samples_section <- function(pool_lu_url, pool_dir) {
  pool_lu_path <- file.path(pool_dir, "pool_lu.csv")
  cat(sprintf("  downloading %s\n", pool_lu_url))
  download.file(pool_lu_url, pool_lu_path, mode = "wb", quiet = TRUE)
  cat(sprintf("  wrote: %s (%.1f MB)\n", pool_lu_path, file.size(pool_lu_path) / 1e6))

  cat("  reading pool to sample mapping ..\n")
  # cellranger multi config file has 4 sections (each its own table with
  # different columns) -- [samples] is the one with the sample<->tag mapping
  lines <- readLines(pool_lu_path)
  start <- which(lines == "[samples]") + 1
  end <- if (any(lines[start:length(lines)] == "")) start + which(lines[start:length(lines)] == "")[1] - 2 else length(lines)

  list(samples = fread(text = lines[start:end]), pool_lu_path = pool_lu_path)
}

# Write the final per-pool h5ad
write_pool_outputs <- function(sce, pool_dir, pool_id, cleanup_paths, pool_lookup = NULL) {
  cat(sprintf("  %s: %d barcodes x %d features\n", pool_id, ncol(sce), nrow(sce)))
  for (path in cleanup_paths) file.remove(path)

  out_path <- file.path(pool_dir, "raw.h5ad")
  if (is.null(pool_lookup)) {
    write_h5ad(sce, out_path, x_mapping = "counts")
  } else {
    metadata(sce)$pool_lookup <- pool_lookup
    write_h5ad(sce, out_path, x_mapping = "counts", uns_mapping = c(pool_lookup = "pool_lookup"))
  }
  cat(sprintf("  wrote: %s\n", out_path))
}

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  write("loading datasets ..", stderr())

  if (args$dataset_name == "hgmm_10x") {

    url_ls <- list(
      hgmm10k_3v4 = "https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_hgmm_3p_gemx_Multiplex/10k_hgmm_3p_gemx_Multiplex_count_raw_feature_bc_matrix.h5",
      hgmm10k_5v3 = "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/10k_hgmm_5p_gemx_Multiplex/10k_hgmm_5p_gemx_Multiplex_count_raw_feature_bc_matrix.h5",
      hgmm5k_3v3  = "https://cf.10xgenomics.com/samples/cell-exp/8.0.0/5k_hgmm_3p_nextgem_Multiplex/5k_hgmm_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5",
      hgmm5k_5v2  = "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/5k_hgmm_5p_nextgem_Multiplex/5k_hgmm_5p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"
    )

    for (pool_name in names(url_ls)) {
      cat(sprintf("  --- pool: %s ---\n", pool_name))
      pool_dir <- file.path(args$output_dir, pool_name)
      dir.create(pool_dir, showWarnings = FALSE, recursive = TRUE)

      dl <- download_and_read_10x(url_ls[[pool_name]], pool_dir)
      sce <- dl$sce
      colData(sce)$pool_id <- pool_name
      colData(sce)$species <- NA_character_ # species not known at this stage for each cell
      colData(sce)$sample_id <- NA_character_ # sample_id not known at this stage

      # add species to rowData based on gene prefix
      genome_prefix <- sub("^([^_]+)_.*", "\\1", rownames(sce))
      species <- rep(NA_character_, length(genome_prefix))
      species[grepl("^(hg|GRCh)", genome_prefix)] <- "human"
      species[grepl("^(mm|GRCm)", genome_prefix)] <- "mouse"
      rowData(sce)$species <- species

      write_pool_outputs(sce, pool_dir, pool_name, cleanup_paths = dl$h5_path)
    }

  } else if (args$dataset_name == "human_pbmcs_hto_10x") {

    h5_url <- "https://cf.10xgenomics.com/samples/cell-exp/9.0.0/60k_Human_TotalSeqB_hashing_3p_gem-x_Multiplex/60k_Human_TotalSeqB_hashing_3p_gem-x_Multiplex_count_raw_feature_bc_matrix.h5"
    pool_lu_url <- "https://cf.10xgenomics.com/samples/cell-exp/9.0.0/60k_Human_TotalSeqB_hashing_3p_gem-x_Multiplex/60k_Human_TotalSeqB_hashing_3p_gem-x_Multiplex_config.csv"
    pool_id <- "totalseqB_hashed_pool"

    pool_dir <- file.path(args$output_dir, pool_id)
    dir.create(pool_dir, showWarnings = FALSE, recursive = TRUE)

    dl <- download_and_read_10x(h5_url, pool_dir)
    sce <- dl$sce
    colData(sce)$pool_id <- pool_id
    colData(sce)$species <- "human"
    colData(sce)$sample_id <- NA_character_ # sample_id not known at this stage

    samples <- download_samples_section(pool_lu_url, pool_dir)
    pool_lu <- samples$samples %>%
      .[, .(hto_id = unlist(strsplit(hashtag_ids, "\\|"))), by = sample_id] %>%
      .[, pool_id := ..pool_id] %>%
      setcolorder(c("pool_id", "sample_id", "hto_id"))

    # check that values in hto_id column match "Antibody Capture" features
    hto_features <- rowData(sce)$feature_type == "Antibody Capture"
    hto_ids <- rowData(sce)$ID[hto_features]
    if (!all(pool_lu$hto_id %in% hto_ids)) {
      stop("Some hto_id values in the pool lookup table do not match 'Antibody Capture' features in the dataset")
    }

    write_pool_outputs(sce, pool_dir, pool_id,
                        cleanup_paths = c(dl$h5_path, samples$pool_lu_path),
                        pool_lookup = pool_lu)

  } else if (args$dataset_name == "mouse_cell_lines_flex_10x") {

    h5_url <- "https://cf.10xgenomics.com/samples/cell-exp/7.1.0/4plex_mouse_cell_lines_multiplex_Multiplex/4plex_mouse_cell_lines_multiplex_Multiplex_count_raw_feature_bc_matrix.h5"
    pool_lu_url <- "https://cf.10xgenomics.com/samples/cell-exp/7.1.0/4plex_mouse_cell_lines_multiplex_Multiplex/4plex_mouse_cell_lines_multiplex_Multiplex_config.csv"
    probe_bcs_lu_url <- "https://cf.10xgenomics.com/supp/cell-exp/probeset/probe-barcodes-fixed-rna-profiling.txt"
    pool_id <- "cell_lines_pool"

    pool_dir <- file.path(args$output_dir, pool_id)
    dir.create(pool_dir, showWarnings = FALSE, recursive = TRUE)

    probe_bcs_lu_path <- file.path(pool_dir, "probe_bcs_lu.txt")
    cat(sprintf("  downloading %s\n", probe_bcs_lu_url))
    download.file(probe_bcs_lu_url, probe_bcs_lu_path, mode = "wb", quiet = TRUE)
    cat(sprintf("  wrote: %s (%.1f MB)\n", probe_bcs_lu_path, file.size(probe_bcs_lu_path) / 1e6))

    # read probe barcodes lookup table (full fixed RNA profiling probe panel)
    cat(sprintf("  reading %s\n", probe_bcs_lu_path))
    probe_bcs_lu <- fread(probe_bcs_lu_path, header = FALSE, col.names = c("barcode_raw", "barcode_ranger", "probe_id"))

    # each probe barcode sequence must map to exactly one probe_id
    n_probe_id_per_barcode <- probe_bcs_lu[, .(n = uniqueN(probe_id)), by = barcode_ranger]
    if (any(n_probe_id_per_barcode$n != 1)) {
      stop("probe_bcs_lu has barcode_ranger value(s) that map to more than one probe_id")
    }

    dl <- download_and_read_10x(h5_url, pool_dir)
    sce <- dl$sce

    # each column in the matrix is a 24 base sequence combining the cell GEM barcode and the sample probe barcode
    # first 16 bp: GEM barcode
    # last 8 bp: sample probe barcode sequence
    # build 'gemID_probeID' cell ids 
    cd <- as.data.table(colData(sce))
    cd[, gem_barcode := substr(Barcode, 1, 16)]
    cd[, probe_barcode := substr(Barcode, 17, 24)]
    cd[probe_bcs_lu, on = .(probe_barcode = barcode_ranger), probe_id := i.probe_id]
    if (any(is.na(cd$probe_id))) {
      stop("some cells have a probe_barcode with no match in probe_bcs_lu")
    }
    cd[, cell_id := paste0(gem_barcode, "_", probe_id)]
    cd[, c("Barcode", "Sample", "gem_barcode", "probe_barcode") := NULL]

    colData(sce) <- DataFrame(cd, row.names = cd$cell_id)
    # colData<- renames only colData; the counts matrix still stores the raw
    # 10x barcodes (set by `counts(sce) <- ...` in download_and_read_10x), and
    # anndataR::write_h5ad checks those against the obs names. Store the
    # counts again so the matrix carries the new cell ids too.
    counts(sce) <- counts(sce)
    colData(sce)$pool_id <- pool_id
    colData(sce)$species <- "mouse"
    colData(sce)$sample_id <- NA_character_ # sample_id not known at this stage

    samples <- download_samples_section(pool_lu_url, pool_dir)
    pool_lu <- samples$samples %>%
      .[, .(sample_id, probe_id = probe_barcode_ids)] %>%
      .[, pool_id := ..pool_id] %>%
      setcolorder(c("pool_id", "sample_id", "probe_id"))

    # keep only cells whose probe was actually used for labelling in this pool
    sce <- sce[, colData(sce)$probe_id %in% pool_lu$probe_id]

    # print the number of cells per probe_id
    cells_per_probe <- table(colData(sce)$probe_id)
    cat("  cells per probe_id:\n")
    for (pid in names(cells_per_probe)) {
      cat(sprintf("    %s: %d\n", pid, cells_per_probe[[pid]]))
    }

    write_pool_outputs(sce, pool_dir, pool_id,
      cleanup_paths = c(dl$h5_path, probe_bcs_lu_path, samples$pool_lu_path),
      pool_lookup = pool_lu
    )
    
  }

  message("  writing properties.yaml ..")
  properties_path <- file.path(args$output_dir, paste0(args$name, "_properties.yaml"))
  write_yaml(
    list(
      batch_var     = args$batch_var,
      sample_var    = args$sample_var,
      labels_var    = args$labels_var,
      tenx_assay    = args$tenx_assay
    ),
    properties_path
  )
  cat(sprintf("  wrote: %s\n", properties_path))
}

if (sys.nframe() == 0L) {
  main()
}
