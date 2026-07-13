#!/usr/bin/env Rscript
# Pre-process RNA-seq PCA for the Shiny app.
#
# Computes the liver-transcriptome PCA (all samples, males-only, females-only,
# as in Scripts/RNA-seq/02-RNAseq_PCA.R) and saves a small, self-contained RDS
# holding just the principal-component coordinates, the % variance per axis and
# the per-sample metadata / phenotypes. The app then loads this without needing
# edgeR / FactoMineR or the full expression matrix at runtime.
#
# Run once (where liver_dge.RDS is available):
#   Rscript Scripts/preprocess_rnaseq_pca.R [path/to/liver_dge.RDS] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
dge_path <- if (length(args) >= 1) args[1] else "Data/liver_dge.RDS"
out_dir  <- if (length(args) >= 2) args[2] else "Data/RNAseq"

if (!file.exists(dge_path)) {
  stop("liver_dge.RDS not found at '", dge_path, "'. ",
       "Pass the path as the first argument, e.g.\n",
       "  Rscript Scripts/preprocess_rnaseq_pca.R Aim2_Step2_manuscript/Data/liver_dge.RDS")
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Reading ", dge_path, " ...")
liver_dge <- readRDS(dge_path)

cpm     <- liver_dge$cpm                 # genes x samples (log-CPM)
samples <- as.data.frame(liver_dge$samples)
if (is.null(cpm) || is.null(samples)) stop("liver_dge is missing $cpm or $samples")
if (is.null(samples$Sample_ID)) stop("liver_dge$samples has no Sample_ID column")

# --- PCA helper: faithful to the paper (FactoMineR, scaled) with a base fallback ---
compute_pca <- function(cpm_mat, ncp = 10) {
  x <- t(cpm_mat)                          # samples x genes
  keep <- apply(x, 2, function(col) is.finite(stats::var(col)) && stats::var(col) > 0)
  x <- x[, keep, drop = FALSE]
  ncp <- min(ncp, nrow(x) - 1, ncol(x))

  if (requireNamespace("FactoMineR", quietly = TRUE)) {
    p <- FactoMineR::PCA(x, graph = FALSE, ncp = ncp, scale.unit = TRUE)
    coord <- as.data.frame(p$ind$coord)
    eig   <- p$eig[, 2]                     # % variance per dim
  } else {
    message("  FactoMineR not installed; falling back to prcomp(scale. = TRUE)")
    p <- stats::prcomp(x, center = TRUE, scale. = TRUE)
    coord <- as.data.frame(p$x[, seq_len(ncp), drop = FALSE])
    colnames(coord) <- paste0("Dim.", seq_len(ncol(coord)))
    eig <- (p$sdev^2 / sum(p$sdev^2) * 100)[seq_len(ncp)]
  }
  rownames(coord) <- rownames(x)
  list(coord = coord, eig = as.numeric(eig))
}

# Attach metadata/phenotypes to PC coordinates, aligned by Sample_ID
build_set <- function(cpm_mat, samples_df) {
  pca <- compute_pca(cpm_mat)
  meta <- samples_df[match(rownames(pca$coord), samples_df$Sample_ID), , drop = FALSE]
  dim_cols <- grep("^Dim\\.", names(meta), value = TRUE)   # avoid name clashes
  if (length(dim_cols)) meta <- meta[, setdiff(names(meta), dim_cols), drop = FALSE]
  coord <- cbind(pca$coord, meta)
  list(coord = coord, eig = pca$eig)
}

message("Computing PCA (all samples) ...")
set_all <- build_set(cpm, samples)

sex_col <- samples$Sex
is_m <- sex_col %in% c("m", "M", "male", "Male")
is_f <- sex_col %in% c("f", "F", "female", "Female")

sets <- list(all = set_all)
if (sum(is_m) >= 3) {
  message("Computing PCA (males) ...")
  sets$male <- build_set(cpm[, is_m, drop = FALSE], samples[is_m, , drop = FALSE])
}
if (sum(is_f) >= 3) {
  message("Computing PCA (females) ...")
  sets$female <- build_set(cpm[, is_f, drop = FALSE], samples[is_f, , drop = FALSE])
}

# --- Catalog the colour-by options ---------------------------------------
dims <- grep("^Dim\\.", names(set_all$coord), value = TRUE)

# Discrete attributes always offered first
discrete_vars <- c(Strain = "Strain", Sex = "Sex")
discrete_vars <- discrete_vars[discrete_vars %in% names(set_all$coord)]

# Continuous phenotypes = numeric columns that are not PC dims or bookkeeping
exclude <- c(dims, "Sample_ID", "ID", "lib.size", "norm.factors", "group",
             "orderrank", "files", "Sample_ID.1")
num_cols <- names(set_all$coord)[vapply(set_all$coord, is.numeric, logical(1))]
num_cols <- setdiff(num_cols, exclude)

# Friendly labels for the phenotypes used in the paper; fall back to a tidy name
pheno_labels <- c(
  Liver_g_24                       = "Liver weight (g)",
  Liver_g_RelBWSac_24              = "Liver weight (g/g BW)",
  sum_all_vacuoles_percentage_24   = "Steatosis (%)",
  Fibrosis_perc_24                 = "Fibrosis (%)",
  ALAT_24                          = "Plasma ALT",
  TIMP1_24                         = "Plasma TIMP-1",
  GDF15_24                         = "Plasma GDF-15",
  KIM1_24                          = "Plasma KIM-1",
  BW_perGain_24                    = "Body weight gain (%)",
  Fat_perc_22                      = "Fat mass (%)",
  Heart_g_24                       = "Heart weight (g)",
  Kidney_g_24                      = "Kidney weight (g)",
  Spleen_g_24                      = "Spleen weight (g)"
)
tidy_label <- function(col) {
  if (col %in% names(pheno_labels)) return(unname(pheno_labels[col]))
  gsub("_", " ", col)
}
# Order: the curated phenotypes first (those present), then the rest alphabetically
curated_present <- intersect(names(pheno_labels), num_cols)
rest <- sort(setdiff(num_cols, curated_present))
ordered_cols <- c(curated_present, rest)
pheno_vars <- stats::setNames(ordered_cols, vapply(ordered_cols, tidy_label, character(1)))

out <- list(
  sets          = sets,
  dims          = dims,
  discrete_vars = discrete_vars,
  pheno_vars    = pheno_vars,
  strain_levels = if (is.factor(samples$Strain)) levels(samples$Strain) else sort(unique(as.character(samples$Strain))),
  n_genes       = nrow(cpm),
  built         = as.character(Sys.Date())
)

out_file <- file.path(out_dir, "rnaseq_pca.rds")
saveRDS(out, out_file)

# Tiny UI sidecar: only the dropdown choices, so rnaseq_pca_ui() can build its
# controls without decompressing the full PCA bundle.
saveRDS(list(set_names = names(sets), dims = dims,
             discrete_vars = discrete_vars, pheno_vars = pheno_vars,
             n_genes = out$n_genes, built = out$built),
        file.path(out_dir, "rnaseq_pca_meta.rds"))

message("\nDone. Saved ", out_file)
message("  sample sets : ", paste(names(sets), collapse = ", "))
message("  PC axes     : ", length(dims))
message("  phenotypes  : ", length(pheno_vars))
