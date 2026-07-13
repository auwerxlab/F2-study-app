#!/usr/bin/env Rscript
# Pre-process the liver proteomics PCA for the Shiny app.
#
# Mirrors Scripts/Proteomics/02-Proteomics-PCA.R: builds a Gene.Symbol x sample
# abundance matrix from Prot_filt, log2(x + 1), runs a scaled PCA (FactoMineR,
# faithful to the paper, with a prcomp fallback) for all samples / males / females,
# and saves a small self-contained bundle (PC coordinates, % variance per axis,
# per-sample metadata + phenotypes). The app then loads this small RDS instead of
# the 7.4 GB Filtered_prots.Rdata and needs neither FactoMineR nor edgeR at runtime.
#
# Phenotypes (liver weight, TIMP-1, fibrosis, ...) for the colour-by options come
# from liver_dge$samples, joined by sample ID. They are OPTIONAL: if liver_dge is
# unavailable the PCA is still built and can be coloured by Strain / Sex / Generation.
#
# Run once (where the data is available):
#   Rscript Scripts/preprocess_proteomics_pca.R \
#     [Filtered_prots.Rdata] [liver_dge.RDS] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
prot_path <- if (length(args) >= 1) args[1] else "Data/Filtered_prots.Rdata"
dge_path  <- if (length(args) >= 2) args[2] else "Data/liver_dge.RDS"
out_dir   <- if (length(args) >= 3) args[3] else "Data/Proteomics"

if (!file.exists(prot_path)) stop("Filtered_prots.Rdata not found at '", prot_path, "'.")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Loading ", prot_path, " (this is large) ...")
pe <- new.env(); load(prot_path, envir = pe)
Prot_filt <- pe$Prot_filt
metadata  <- as.data.frame(pe$metadata)
if (is.null(Prot_filt) || is.null(metadata)) stop("Filtered_prots.Rdata must contain Prot_filt and metadata")

# --- Gene.Symbol x sample abundance matrix (as in 02-Proteomics-PCA.R) -----------
annot_cols <- intersect(c("protein.groups", "PG.FastaHeaders"), names(Prot_filt))
m <- Prot_filt[, setdiff(names(Prot_filt), annot_cols), drop = FALSE]
m <- m[!is.na(m$Gene.Symbol) & m$Gene.Symbol != "", , drop = FALSE]
data_cols <- setdiff(names(m), "Gene.Symbol")
# For duplicated Gene.Symbols keep the row detected in the most samples (fewest NAs)
m$.na <- rowSums(is.na(m[, data_cols, drop = FALSE]))
m <- m[order(m$Gene.Symbol, m$.na), ]
m <- m[!duplicated(m$Gene.Symbol), ]
rownames(m) <- m$Gene.Symbol
mat     <- as.matrix(m[, data_cols, drop = FALSE])     # genes x samples
mat_log <- log2(mat + 1)
sample_ids <- colnames(mat_log)                        # HDP-xxxx
message("  matrix: ", nrow(mat_log), " proteins x ", ncol(mat_log), " samples")

# --- Per-sample metadata: Strain / Sex / Generation from the proteomics metadata,
# --- enriched with numeric phenotypes from liver_dge$samples (joined by ID). ------
meta <- metadata[match(sample_ids, metadata$ID), , drop = FALSE]
meta$ID        <- sample_ids
meta$Sample_ID <- sample_ids                           # key used by build_set()

phen <- tryCatch({
  if (!file.exists(dge_path)) { message("  liver_dge not found - no phenotype colours"); NULL }
  else {
    dge <- readRDS(dge_path)
    s <- as.data.frame(dge$samples)
    key <- if ("ID" %in% names(s) && any(s$ID %in% sample_ids)) "ID"
           else if ("Sample_ID" %in% names(s) && any(as.character(s$Sample_ID) %in% sample_ids)) "Sample_ID"
           else NA_character_
    if (is.na(key)) { message("  liver_dge$samples has no ID matching the proteomics samples"); NULL }
    else {
      num <- names(s)[vapply(s, is.numeric, logical(1))]
      num <- setdiff(num, names(meta))                 # don't duplicate existing cols
      s2 <- s[, c(key, num), drop = FALSE]
      names(s2)[1] <- "ID"; s2$ID <- as.character(s2$ID)
      message("  joined ", length(num), " phenotype columns from liver_dge")
      s2
    }
  }
}, error = function(e) { message("  (phenotypes skipped: ", conditionMessage(e), ")"); NULL })

if (!is.null(phen)) {
  meta <- merge(meta, phen, by = "ID", all.x = TRUE, sort = FALSE)
  meta <- meta[match(sample_ids, meta$ID), , drop = FALSE]   # restore matrix-column order
}

# --- PCA (scaled), faithful to FactoMineR with a prcomp fallback -----------------
compute_pca <- function(mat_gs, ncp = 10) {
  x <- t(mat_gs)                                        # samples x genes
  # FactoMineR imputes missing values to the column mean; do it explicitly so the
  # prcomp fallback (which cannot handle NA) behaves the same.
  cm <- colMeans(x, na.rm = TRUE)
  na_idx <- which(is.na(x), arr.ind = TRUE)
  if (nrow(na_idx)) x[na_idx] <- cm[na_idx[, "col"]]
  keep <- apply(x, 2, function(col) is.finite(stats::var(col)) && stats::var(col) > 0)
  x <- x[, keep, drop = FALSE]
  ncp <- min(ncp, nrow(x) - 1, ncol(x))
  if (requireNamespace("FactoMineR", quietly = TRUE)) {
    p <- FactoMineR::PCA(x, graph = FALSE, ncp = ncp, scale.unit = TRUE)
    coord <- as.data.frame(p$ind$coord); eig <- p$eig[, 2]
  } else {
    message("  FactoMineR not installed; using prcomp(scale. = TRUE)")
    p <- stats::prcomp(x, center = TRUE, scale. = TRUE)
    coord <- as.data.frame(p$x[, seq_len(ncp), drop = FALSE])
    colnames(coord) <- paste0("Dim.", seq_len(ncol(coord)))
    eig <- (p$sdev^2 / sum(p$sdev^2) * 100)[seq_len(ncp)]
  }
  rownames(coord) <- rownames(x)
  list(coord = coord, eig = as.numeric(eig))
}

build_set <- function(mat_sub, samples_df) {
  pca <- compute_pca(mat_sub)
  md <- samples_df[match(rownames(pca$coord), samples_df$Sample_ID), , drop = FALSE]
  dim_cols <- grep("^Dim\\.", names(md), value = TRUE)
  if (length(dim_cols)) md <- md[, setdiff(names(md), dim_cols), drop = FALSE]
  list(coord = cbind(pca$coord, md), eig = pca$eig)
}

message("Computing PCA (all samples) ...")
set_all <- build_set(mat_log, meta)

is_m <- meta$Sex %in% c("m", "M", "male", "Male")
is_f <- meta$Sex %in% c("f", "F", "female", "Female")
sets <- list(all = set_all)
if (sum(is_m) >= 3) { message("Computing PCA (males) ...");   sets$male   <- build_set(mat_log[, is_m, drop = FALSE], meta[is_m, , drop = FALSE]) }
if (sum(is_f) >= 3) { message("Computing PCA (females) ..."); sets$female <- build_set(mat_log[, is_f, drop = FALSE], meta[is_f, , drop = FALSE]) }

# --- Catalogue the colour-by options (mirrors preprocess_rnaseq_pca.R) -----------
dims <- grep("^Dim\\.", names(set_all$coord), value = TRUE)

discrete_vars <- c(Strain = "Strain", Sex = "Sex", Generation = "Generation")
discrete_vars <- discrete_vars[discrete_vars %in% names(set_all$coord)]

exclude <- c(dims, "Sample_ID", "ID", "lib.size", "norm.factors", "group",
             "orderrank", "files", "Sample_ID.1")
num_cols <- names(set_all$coord)[vapply(set_all$coord, is.numeric, logical(1))]
num_cols <- setdiff(num_cols, exclude)

pheno_labels <- c(
  Liver_g_24                     = "Liver weight (g)",
  Liver_g_RelBWSac_24            = "Liver weight (g/g BW)",
  sum_all_vacuoles_percentage_24 = "Steatosis (%)",
  Fibrosis_perc_24               = "Fibrosis (%)",
  ALAT_24                        = "Plasma ALT",
  TIMP1_24                       = "Plasma TIMP-1",
  GDF15_24                       = "Plasma GDF-15",
  KIM1_24                        = "Plasma KIM-1",
  BW_perGain_24                  = "Body weight gain (%)",
  Fat_perc_22                    = "Fat mass (%)",
  Heart_g_24                     = "Heart weight (g)",
  Kidney_g_24                    = "Kidney weight (g)",
  Spleen_g_24                    = "Spleen weight (g)"
)
tidy_label <- function(col) if (col %in% names(pheno_labels)) unname(pheno_labels[col]) else gsub("_", " ", col)
curated_present <- intersect(names(pheno_labels), num_cols)
ordered_cols <- c(curated_present, sort(setdiff(num_cols, curated_present)))
pheno_vars <- stats::setNames(ordered_cols, vapply(ordered_cols, tidy_label, character(1)))

out <- list(
  sets          = sets,
  dims          = dims,
  discrete_vars = discrete_vars,
  pheno_vars    = pheno_vars,
  strain_levels = if (is.factor(meta$Strain)) levels(meta$Strain) else sort(unique(as.character(meta$Strain))),
  n_proteins    = nrow(mat),
  built         = as.character(Sys.Date())
)

out_file <- file.path(out_dir, "proteomics_pca.rds")
saveRDS(out, out_file)

# Tiny UI sidecar: only the dropdown choices, so proteomics_pca_ui() can build its
# controls without decompressing the full PCA bundle.
saveRDS(list(set_names = names(sets), dims = dims,
             discrete_vars = discrete_vars, pheno_vars = pheno_vars,
             n_proteins = out$n_proteins, built = out$built),
        file.path(out_dir, "proteomics_pca_meta.rds"))

message("\nDone. Saved ", out_file)
message("  sample sets : ", paste(names(sets), collapse = ", "))
message("  PC axes     : ", length(dims))
message("  proteins    : ", out$n_proteins)
message("  phenotypes  : ", length(pheno_vars))

