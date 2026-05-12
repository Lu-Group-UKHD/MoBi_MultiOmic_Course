#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(pbapply)
  library(stringr)
  library(SingleCellExperiment)
  library(scrapper)
  library(SparseArray)
  library(scran)
  library(scater)
  library(sketchR)
})

CITE_GSMS <- c(
  "GSM5613744", "GSM5613745", "GSM5613746",
  "GSM5613747", "GSM5613748", "GSM5613749", "GSM5613750",
  "GSM5613751", "GSM5613752",
  "GSM5613753", "GSM5613754", "GSM5613755",
  "GSM5613756",
  "GSM5613757", "GSM5613758", "GSM5613759",
  "GSM5613760", "GSM5613761",
  "GSM5613762", "GSM5613763",
  "GSM5613764", "GSM5613765",
  "GSM5613766", "GSM5613767",
  # GSM5613768 is scRNA-only, skipped
  "GSM5613769", "GSM5613770", "GSM5613771",
  "GSM5613772", "GSM5613773",
  "GSM5613774", "GSM5613775", "GSM5613776",
  "GSM5613777", "GSM5613778",
  "GSM5613779", "GSM5613780",
  "GSM5613781",
  "GSM5613782",
  "GSM5613783", "GSM5613784", "GSM5613785", "GSM5613786",
  "GSM5613787",
  "GSM5613788", "GSM5613789", "GSM5613790"
)

usage <- function() {
  cat(
    "\nUsage:\n",
    "  Rscript read_gse185381_cite.R <gse_dir> <output>\n\n",
    "Arguments:\n",
    "  <gse_dir>  Directory containing the downloaded GEO RAW files.\n",
    "  <output>   Output .rds file or output directory.\n\n",
    "Examples:\n",
    "  Rscript read_gse185381_cite.R data/GSE185381_RAW output\n",
    "  Rscript read_gse185381_cite.R data/GSE185381_RAW output/GSE185381_cite_sce.rds\n\n",
    sep = ""
  )
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  usage()
}

gse_dir <- args[[1]]
output <- args[[2]]

if (!dir.exists(gse_dir)) {
  stop("GSE directory does not exist: ", gse_dir)
}

if (grepl("\\.rds$", output, ignore.case = TRUE)) {
  output_file <- output
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
} else {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output, "GSE185381_cite_sce.rds")
}

message("Input directory:  ", normalizePath(gse_dir))
message("Output file:      ", output_file)

find_one_file <- function(files, gsm, pattern, ignore.case = TRUE) {
  hits <- grep(pattern, files, value = TRUE, ignore.case = ignore.case)
  hits <- hits[grepl(gsm, hits, fixed = TRUE)]
  
  if (length(hits) == 0) {
    stop("No file found for ", gsm, " matching pattern: ", pattern)
  }
  
  if (length(hits) > 1) {
    stop(
      "Multiple files found for ", gsm, " matching pattern: ", pattern, "\n",
      paste(hits, collapse = "\n")
    )
  }
  
  hits[[1]]
}

read_mtx <- function(path) {
  if (grepl("\\.gz$", path)) {
    Matrix::readMM(gzfile(path))
  } else {
    Matrix::readMM(path)
  }
}

normalize_cell_name <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace("^X", "") |>
    stringr::str_replace_all(":", ".") |>
    stringr::str_replace_all("-", ".")
}

extract_barcode_from_meta_cell <- function(x) {
  normalized <- normalize_cell_name(x)
  parts <- stringr::str_split(normalized, stringr::fixed("."))
  vapply(parts, function(y) y[[length(y)]], character(1))
}

read_one_gsm <- function(gsm, gse_dir) {
  message("Reading ", gsm)
  
  files <- list.files(gse_dir, pattern = gsm, full.names = FALSE)
  
  matrix_file <- find_one_file(files, gsm, "matrix")
  features_file <- find_one_file(files, gsm, "features")
  barcodes_file <- find_one_file(files, gsm, "barcodes")
  meta_file <- find_one_file(files, gsm, "meta")
  adt_file <- find_one_file(files, gsm, "ADT")
  
  matrix_path <- file.path(gse_dir, matrix_file)
  features_path <- file.path(gse_dir, features_file)
  barcodes_path <- file.path(gse_dir, barcodes_file)
  meta_path <- file.path(gse_dir, meta_file)
  adt_path <- file.path(gse_dir, adt_file)
  
  mat <- read_mtx(matrix_path)
  # features <- read.table(
  #   features_path,
  #   row.names = 1,
  #   sep = "\t",
  #   header = FALSE,
  #   stringsAsFactors = FALSE,
  #   check.names = FALSE
  # )
  
  ###
  
  # get correct gene expression
  
  features <- read.table(
    features_path,
    row.names = 1,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  # After row.names = 1, the columns are usually:
  #   V1 = feature symbol/name
  #   V2 + V3 = split feature type, e.g. "Gene" + "Expression"
  # or
  #   V2 + V3 = "Antibody" + "Abundance"
  #
  # We reconstruct the feature type instead of assuming one clean column.
  if (ncol(features) < 2) {
    stop(
      "Feature file for ", gsm,
      " has fewer than two columns after using the first column as row names."
    )
  }
  
  feature_symbol <- features[[1]]
  
  feature_type <- apply(
    features[, -1, drop = FALSE],
    1,
    function(x) {
      paste(na.omit(as.character(x)), collapse = " ")
    }
  )
  
  feature_type <- stringr::str_squish(feature_type)
  
  features$feature_symbol <- feature_symbol
  features$feature_type <- feature_type
  
  gene_rows <- which(features$feature_type == "Gene Expression")
  
  if (length(gene_rows) == 0) {
    stop(
      "No 'Gene Expression' rows found in feature file for ", gsm, ". ",
      "Observed feature types were: ",
      paste(unique(features$feature_type), collapse = ", ")
    )
  }
  
  ###
  
  barcodes <- read.table(
    barcodes_path,
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]
  
  meta <- read.csv(
    meta_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  adt_counts <- read.csv(
    adt_path,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (ncol(mat) != length(barcodes)) {
    stop(
      "Number of RNA matrix columns does not match barcode file for ", gsm,
      ": ncol(mat) = ", ncol(mat),
      ", length(barcodes) = ", length(barcodes)
    )
  }
  
  if (!"cell" %in% colnames(meta)) {
    stop("Metadata for ", gsm, " does not contain a 'cell' column.")
  }
  
  colnames(mat) <- barcodes
  
  rna_barcodes <- stringr::str_replace(colnames(mat), "-1$", "")
  meta$cell <- normalize_cell_name(meta$cell)
  meta$barcode <- extract_barcode_from_meta_cell(meta$cell)
  
  matched_meta <- match(rna_barcodes, meta$barcode)
  keep <- which(!is.na(matched_meta))
  
  if (length(keep) == 0) {
    stop("No RNA barcodes matched metadata cells for ", gsm)
  }
  
  mat <- mat[, keep, drop = FALSE]
  meta <- meta[matched_meta[keep], , drop = FALSE]
  
  rownames(meta) <- meta$cell
  colnames(mat) <- meta$cell
  
  # if ("V3" %in% colnames(features)) {
  #   gene_rows <- which(features$V3 == "Gene")
  #   if (length(gene_rows) == 0) {
  #     stop("Feature file for ", gsm, " has V3 column but no rows with V3 == 'Gene'.")
  #   }
  #   
  #   mat <- mat[gene_rows, , drop = FALSE]
  #   features <- features[gene_rows, , drop = FALSE]
  # } else {
  #   warning("Feature file for ", gsm, " has no V3 column; keeping all rows.")
  # }
  
  ### use new gene rows:
  
  mat <- mat[gene_rows, , drop = FALSE]
  features <- features[gene_rows, , drop = FALSE]
  
  sce <- SingleCellExperiment(
    assays = list(counts = mat),
    rowData = S4Vectors::DataFrame(features),
    colData = S4Vectors::DataFrame(meta)
  )
  
  sce$sample_gsm <- gsm
  
  colnames(adt_counts) <- normalize_cell_name(colnames(adt_counts))
  
  missing_in_adt <- setdiff(colnames(sce), colnames(adt_counts))
  extra_in_adt <- setdiff(colnames(adt_counts), colnames(sce))
  
  if (length(missing_in_adt) > 0) {
    stop(
      "ADT matrix for ", gsm, " is missing ", length(missing_in_adt),
      " cells present in RNA/SCE. First missing cell: ", missing_in_adt[[1]]
    )
  }
  
  if (length(extra_in_adt) > 0) {
    message(
      "ADT matrix for ", gsm, " has ", length(extra_in_adt),
      " extra cells not present in RNA/SCE; dropping them."
    )
  }
  
  adt_counts <- adt_counts[, colnames(sce), drop = FALSE]
  
  adt_sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(adt_counts))
  )
  
  stopifnot(identical(colnames(adt_sce), colnames(sce)))
  
  altExp(sce, "ADT") <- adt_sce
  
  message(
    "  RNA genes: ", nrow(sce),
    " | cells: ", ncol(sce),
    " | ADTs: ", nrow(altExp(sce, "ADT"))
  )
  
  sce
}

sce_list <- pbapply::pblapply(
  CITE_GSMS,
  read_one_gsm,
  gse_dir = gse_dir
)

message("Intersecting RNA feature universe across samples.")

common_genes <- Reduce(intersect, lapply(sce_list, rownames))

if (length(common_genes) == 0) {
  stop("No common RNA genes found across all samples.")
}

sce_list <- lapply(sce_list, function(sce) {
  sce[common_genes, ]
})

message("Combining samples.")

sce <- do.call(cbind, sce_list)

message("Final object:")
message("  RNA genes: ", nrow(sce))
message("  Cells:     ", ncol(sce))
message("  ADTs:      ", nrow(altExp(sce, "ADT")))

adt_sum <- sum(SummarizedExperiment::assay(altExp(sce, "ADT"), "counts"))

if (adt_sum == 0) {
  warning("ADT count matrix sums to zero. This suggests the ADT import still failed.")
} else {
  message("  ADT count sum: ", adt_sum)
}

### geometricsketch

counts(sce) <- counts(sce) |>
  as("CsparseMatrix")

num.threads = 16

sf <- scrapper::centerSizeFactors(colSums(counts(sce)))
sce <- normalizeRnaCounts.se(sce, sf)
sce <- chooseRnaHvgs.se(sce, num.threads=num.threads, block = sce$samples, top = 3000)
sce <- runPca.se(sce, features=rowData(sce)$hvg, num.threads=num.threads)

reducedDim(sce, "PCA") |> class()

idx <- sketchR::geosketch(
  reducedDim(sce, "PCA"),
  N = 30000,
  replace = FALSE,
  seed = 42
)

sub <- sce[, idx]

# library(alabaster.sce)
# library(alabaster.base)

# saveObject(sce, file.path(output, "lasry_all_artifactdb"))

# saveObject(sub, file.path(output, "lasry_sketch_artifactdb"))
library(anndataR)

path <- output

counts(sce) <- counts(sce) |>
  as("CsparseMatrix")
logcounts(sce) <- logcounts(sce) |>
  as("CsparseMatrix")

counts(altExp(sce, "ADT")) <- as.matrix(counts(altExp(sce, "ADT")))


adata <- anndataR::as_AnnData(sce)
anndataR::write_h5ad(adata, file.path(path, "lasry_all_rna.h5ad"))

adata <- anndataR::as_AnnData(altExp(sce, "ADT"))
anndataR::write_h5ad(adata, file.path(path, "lasry_all_adt.h5ad"))

adata <- anndataR::as_AnnData(sub)
anndataR::write_h5ad(adata, file.path(path, "lasry_sketch_rna.h5ad"))

adata <- anndataR::as_AnnData(altExp(sub, "ADT"))
anndataR::write_h5ad(adata, file.path(path, "lasry_sketch_adt.h5ad"))

message("Saved: ", output)