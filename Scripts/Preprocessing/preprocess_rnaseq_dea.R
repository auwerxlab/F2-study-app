#!/usr/bin/env Rscript
# Pre-process RNA-seq differential expression for the Shiny app.
#
# Fits the limma-voom models from Scripts/RNA-seq/04-RNAseq_DEA.R and saves the
# per-contrast result tables (topTable + gene names). The app then browses them
# (volcano + table) without needing limma/edgeR or the counts at runtime.
#
# Contrasts:
#   - Sex (male vs female), F2 samples
#   - Fat mass (F2 females)
#   - Multivariate liver-disease model (F2, complete cases), one table per term:
#       Fibrosis, Sex, plasma TIMP-1, liver weight (g/g BW), steatosis, fat mass
#
# Run once (where liver_dge.RDS is available; needs limma + edgeR):
#   Rscript Scripts/preprocess_rnaseq_dea.R [path/to/liver_dge.RDS] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
dge_path <- if (length(args) >= 1) args[1] else "Data/liver_dge.RDS"
out_dir  <- if (length(args) >= 2) args[2] else "Data/RNAseq"

if (!file.exists(dge_path)) {
  stop("liver_dge.RDS not found at '", dge_path, "'. Pass its path as the first argument.")
}
for (pkg in c("limma", "edgeR")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Reading ", dge_path, " ...")
liver_dge <- readRDS(dge_path)
samples <- as.data.frame(liver_dge$samples)
counts  <- liver_dge$counts_filtered
gct     <- liver_dge$geneConversionTable
if (is.null(counts) || is.null(samples)) stop("liver_dge missing $counts_filtered or $samples")

gene_cols <- intersect(c("gene_id", "gene_name", "gene_biotype"), names(gct))

# Fit one model (mask + formula) and return one tidy table per requested coef.
fit_model <- function(mask, formula, coefs) {
  sub <- samples[mask, , drop = FALSE]
  design <- model.matrix(formula, data = sub)
  vm  <- limma::voom(counts[, mask, drop = FALSE], design = design)
  fit <- limma::eBayes(limma::lmFit(vm, design = design), robust = TRUE)

  out <- list()
  for (key in names(coefs)) {
    coef <- coefs[[key]]
    if (!coef %in% colnames(design)) {
      message("  ! coef not found, skipping: ", coef)
      next
    }
    tt <- limma::topTable(fit, coef = coef, n = Inf, sort.by = "p")
    tt$gene_id <- rownames(tt)
    tt <- merge(tt, gct[, gene_cols, drop = FALSE], by = "gene_id", all.x = TRUE)
    keep <- c("gene_id", "gene_name", "gene_biotype",
              "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
    tt <- tt[, intersect(keep, names(tt)), drop = FALSE]
    tt <- tt[order(tt$adj.P.Val, tt$P.Value), ]
    rownames(tt) <- NULL
    out[[key]] <- list(table = tt, n_samples = nrow(sub))
  }
  out
}

contrasts <- list()

# --- 1) Sex: male vs female (F2) ---
message("Fitting: Sex (male vs female) ...")
mask_sex <- with(samples, Generation == "F2" & !is.na(Sex))
res <- fit_model(mask_sex, ~Sex, c(sex = "Sexm"))
if (!is.null(res$sex)) contrasts$sex <- c(res$sex,
  list(label = "Sex: male vs female",
       description = "Differential expression between sexes (F2). Positive logFC = higher in males."))

# --- 2) Fat mass (F2 females) ---
message("Fitting: Fat mass (females) ...")
mask_fat <- with(samples, Generation == "F2" & Sex == "f" & !is.na(Fat_perc_22))
res <- fit_model(mask_fat, ~Fat_perc_22, c(fat_female = "Fat_perc_22"))
if (!is.null(res$fat_female)) contrasts$fat_female <- c(res$fat_female,
  list(label = "Fat mass (females)",
       description = "Association with fat mass (% at 22 wk) in F2 females."))

# --- 3) Multivariate liver-disease model (F2, complete cases) ---
message("Fitting: multivariate liver-disease model ...")
mask_fib <- with(samples, Generation == "F2" &
                   !is.na(Fibrosis_perc_24) & !is.na(TIMP1_24) &
                   !is.na(Liver_g_RelBWSac_24) & !is.na(sum_all_vacuoles_percentage_24) &
                   !is.na(Fat_perc_22))
fib_coefs <- c(
  fibrosis  = "log(Fibrosis_perc_24)",
  sex_adj   = "Sexm",
  timp1     = "log(TIMP1_24)",
  liverw    = "Liver_g_RelBWSac_24",
  steatosis = "sum_all_vacuoles_percentage_24",
  fat_adj   = "Fat_perc_22"
)
res <- fit_model(
  mask_fib,
  ~log(Fibrosis_perc_24) + Sex + log(TIMP1_24) + Liver_g_RelBWSac_24 +
    sum_all_vacuoles_percentage_24 + Fat_perc_22,
  fib_coefs)
fib_labels <- c(
  fibrosis  = "Fibrosis (%) [adj.]",
  sex_adj   = "Sex (male) [adj.]",
  timp1     = "Plasma TIMP-1 [adj.]",
  liverw    = "Liver weight (g/g BW) [adj.]",
  steatosis = "Steatosis (%) [adj.]",
  fat_adj   = "Fat mass (%) [adj.]"
)
for (key in names(fib_coefs)) {
  if (is.null(res[[key]])) next
  contrasts[[key]] <- c(res[[key]],
    list(label = unname(fib_labels[key]),
         description = paste0("Association with ", unname(fib_labels[key]),
                              " in the multivariate liver-disease model (F2, complete cases), ",
                              "adjusted for the other terms.")))
}

contrast_choices <- setNames(names(contrasts), vapply(contrasts, function(x) x$label, character(1)))

out <- list(
  contrasts        = contrasts,
  contrast_choices = contrast_choices,
  built            = as.character(Sys.Date())
)
out_file <- file.path(out_dir, "rnaseq_dea.rds")
saveRDS(out, out_file)

# Tiny UI sidecar: only the contrast dropdown choices, so rnaseq_dea_ui() can
# build its controls without decompressing the full DEA result tables.
saveRDS(list(contrast_choices = contrast_choices, built = out$built),
        file.path(out_dir, "rnaseq_dea_meta.rds"))

message("\nDone. Saved ", out_file)
message("  contrasts: ", paste(names(contrasts), collapse = ", "))
