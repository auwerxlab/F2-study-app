#!/usr/bin/env Rscript
# Pre-process RNA-seq hierarchical clustering for the Shiny app.
#
# Computes the liver-transcriptome sample dendrogram (all / males / females, as
# in Scripts/RNA-seq/03-RNAseq_HC.R) and saves the hclust objects plus the
# per-sample metadata / phenotypes. The app draws the dendrogram with strain-
# coloured branches and phenotype/sex/strain colour bars (dendextend) at runtime.
#
# Run once (where liver_dge.RDS is available):
#   Rscript Scripts/preprocess_rnaseq_hc.R [path/to/liver_dge.RDS] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
dge_path <- if (length(args) >= 1) args[1] else "Data/liver_dge.RDS"
out_dir  <- if (length(args) >= 2) args[2] else "Data/RNAseq"

if (!file.exists(dge_path)) {
  stop("liver_dge.RDS not found at '", dge_path, "'. ",
       "Pass the path as the first argument, e.g.\n",
       "  Rscript Scripts/preprocess_rnaseq_hc.R Aim2_Step2_manuscript/Data/liver_dge.RDS")
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Reading ", dge_path, " ...")
liver_dge <- readRDS(dge_path)

cpm     <- liver_dge$cpm                 # genes x samples (log-CPM)
samples <- as.data.frame(liver_dge$samples)
if (is.null(cpm) || is.null(samples)) stop("liver_dge is missing $cpm or $samples")
if (is.null(samples$Sample_ID)) stop("liver_dge$samples has no Sample_ID column")

# Cluster samples (Euclidean distance on log-CPM, complete linkage), and align
# the metadata to the clustering order so the app's colour bars line up.
build_set <- function(cpm_mat, samples_df) {
  x  <- t(cpm_mat)                        # samples x genes
  hc <- stats::hclust(stats::dist(x))     # labels = rownames(x) = Sample_ID
  meta <- samples_df[match(hc$labels, samples_df$Sample_ID), , drop = FALSE]
  rownames(meta) <- NULL
  list(hc = hc, meta = meta)
}

message("Clustering (all samples) ...")
sets <- list(all = build_set(cpm, samples))

sex_col <- samples$Sex
is_m <- sex_col %in% c("m", "M", "male", "Male")
is_f <- sex_col %in% c("f", "F", "female", "Female")
if (sum(is_m) >= 3) { message("Clustering (males) ...");   sets$male   <- build_set(cpm[, is_m, drop = FALSE], samples[is_m, , drop = FALSE]) }
if (sum(is_f) >= 3) { message("Clustering (females) ..."); sets$female <- build_set(cpm[, is_f, drop = FALSE], samples[is_f, , drop = FALSE]) }

# --- Catalog the phenotype options (numeric columns, friendly labels) -------
meta_all <- sets$all$meta
exclude  <- c("Sample_ID", "ID", "lib.size", "norm.factors", "group",
              "orderrank", "files", "Sample_ID.1")
num_cols <- names(meta_all)[vapply(meta_all, is.numeric, logical(1))]
num_cols <- setdiff(num_cols, exclude)

pheno_labels <- c(
  TIMP1_24                         = "Plasma TIMP-1",
  ALAT_24                          = "Plasma ALT",
  GDF15_24                         = "Plasma GDF-15",
  KIM1_24                          = "Plasma KIM-1",
  BW_perGain_24                    = "Body weight gain (%)",
  Fat_perc_22                      = "Fat mass (%)",
  Liver_g_24                       = "Liver weight (g)",
  Liver_g_RelBWSac_24              = "Liver weight (g/g BW)",
  sum_all_vacuoles_percentage_24   = "Steatosis (%)",
  Fibrosis_perc_24                 = "Fibrosis (%)",
  Heart_g_24                       = "Heart weight (g)",
  Kidney_g_24                      = "Kidney weight (g)",
  Spleen_g_24                      = "Spleen weight (g)"
)
tidy_label <- function(col) if (col %in% names(pheno_labels)) unname(pheno_labels[col]) else gsub("_", " ", col)
curated_present <- intersect(names(pheno_labels), num_cols)
rest <- sort(setdiff(num_cols, curated_present))
ordered_cols <- c(curated_present, rest)
pheno_vars <- stats::setNames(ordered_cols, vapply(ordered_cols, tidy_label, character(1)))

# Default colour bars used in the paper figure (those present)
default_bars <- intersect(
  c("TIMP1_24", "ALAT_24", "BW_perGain_24", "Fat_perc_22",
    "Liver_g_RelBWSac_24", "sum_all_vacuoles_percentage_24", "Fibrosis_perc_24"),
  ordered_cols)

out <- list(
  sets          = sets,
  pheno_vars    = pheno_vars,
  default_bars  = default_bars,
  strain_levels = if (is.factor(samples$Strain)) levels(samples$Strain) else sort(unique(as.character(samples$Strain))),
  n_genes       = nrow(cpm),
  built         = as.character(Sys.Date())
)

out_file <- file.path(out_dir, "rnaseq_hc.rds")
saveRDS(out, out_file)

# Tiny UI sidecar: only the dropdown choices, so rnaseq_hc_ui() can build its
# controls without decompressing the full clustering bundle.
saveRDS(list(set_names = names(sets), pheno_vars = pheno_vars,
             default_bars = default_bars, n_genes = out$n_genes, built = out$built),
        file.path(out_dir, "rnaseq_hc_meta.rds"))

message("\nDone. Saved ", out_file)
message("  sample sets : ", paste(names(sets), collapse = ", "))
message("  phenotypes  : ", length(pheno_vars))
message("  default bars: ", paste(default_bars, collapse = ", "))
