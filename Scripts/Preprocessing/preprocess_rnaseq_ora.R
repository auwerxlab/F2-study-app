#!/usr/bin/env Rscript
# Pre-process set-specific GO over-representation (ORA) for the RNA-seq DEA tab.
#
# Reproduces the UpSet/enrichGO analysis in Scripts/RNA-seq/04-RNAseq_DEA.R:
# using the multivariate liver-disease model, split the significant genes for
# each pair of phenotype coefficients into "A-specific", "B-specific" and
# "shared" sets (separately for up- and down-regulated), and run enrichGO (BP)
# on each. The app then shows the diverging up/down GO dot plot per combination.
#
# Reuses the precomputed DEA tables (Data/RNAseq/rnaseq_dea.rds) -- no limma
# refit. Needs clusterProfiler + org.Mm.eg.db (Ensembl mouse gene IDs).
#
# Run after preprocess_rnaseq_dea.R:
#   Rscript Scripts/preprocess_rnaseq_ora.R [dea_rds] [out_dir] [fdr]

args     <- commandArgs(trailingOnly = TRUE)
dea_path <- if (length(args) >= 1) args[1] else "Data/RNAseq/rnaseq_dea.rds"
out_dir  <- if (length(args) >= 2) args[2] else "Data/RNAseq"
fdr      <- if (length(args) >= 3) as.numeric(args[3]) else 0.05

if (!file.exists(dea_path)) stop("rnaseq_dea.rds not found at '", dea_path, "'. Run preprocess_rnaseq_dea.R first.")
for (pkg in c("clusterProfiler", "org.Mm.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
}
# enrichGO() below is called with `OrgDb = org.Mm.eg.db` (a bare symbol), so the
# package must be ATTACHED, not merely loaded: requireNamespace() alone leaves that
# object unresolved, so every enrichGO() call errors with "object 'org.Mm.eg.db' not
# found" - swallowed by the tryCatch in run_ora() - and EVERY set comes back empty
# (the app then shows "No enriched GO terms for this set" for all combinations).
suppressPackageStartupMessages(library(org.Mm.eg.db))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dea <- readRDS(dea_path)

# Coefficients of the multivariate liver-disease model (same fit -> comparable).
multivar <- c("fibrosis", "sex_adj", "timp1", "liverw", "steatosis", "fat_adj")
coefs <- intersect(multivar, names(dea$contrasts))
if (length(coefs) < 2) stop("Need >= 2 multivariate-model contrasts in the DEA bundle; found: ",
                            paste(coefs, collapse = ", "))
coef_choices <- setNames(coefs, vapply(coefs, function(k) dea$contrasts[[k]]$label, character(1)))

universe <- unique(dea$contrasts[[coefs[1]]]$table$gene_id)

# Up / down significant gene sets per coefficient at the chosen FDR.
sig_sets <- lapply(coefs, function(k) {
  tb <- dea$contrasts[[k]]$table
  sig <- !is.na(tb$adj.P.Val) & tb$adj.P.Val < fdr
  list(up = tb$gene_id[sig & tb$logFC > 0], down = tb$gene_id[sig & tb$logFC < 0])
})
names(sig_sets) <- coefs

empty_ora <- function() data.frame(
  ID = character(), Description = character(), GeneRatio = character(),
  GeneRatioNum = numeric(), Count = integer(), pvalue = numeric(),
  p.adjust = numeric(), geneID = character(), stringsAsFactors = FALSE)

run_ora <- function(genes, universe) {
  genes <- unique(genes[!is.na(genes)])
  if (length(genes) < 5) return(empty_ora())
  do_enrich <- function(g, uni) tryCatch(
    clusterProfiler::enrichGO(gene = g, ont = "BP", universe = uni,
                              OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                              pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
    error = function(e) NULL)
  e <- do_enrich(genes, universe)
  # Ensembl IDs sometimes carry a version suffix (e.g. ENSMUSG...".4") that fails
  # to map to GO and yields zero terms; if nothing came back, retry with the
  # versions stripped (a common cause of "no enriched terms" for large gene sets).
  if (is.null(e) || nrow(e@result) == 0) {
    e <- do_enrich(sub("\\.\\d+$", "", genes), sub("\\.\\d+$", "", universe))
  }
  if (is.null(e) || nrow(e@result) == 0) return(empty_ora())
  r <- e@result
  r$GeneRatioNum <- vapply(strsplit(r$GeneRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  keep <- c("ID", "Description", "GeneRatio", "GeneRatioNum", "Count", "pvalue", "p.adjust", "geneID")
  r <- r[, intersect(keep, names(r)), drop = FALSE]
  r <- r[order(r$p.adjust), , drop = FALSE]
  # Keep the most significant terms even when none clear the old 0.25 cutoff, so
  # the app can still display the top enrichments (its p.adjust slider / auto-tune
  # decides what to show). Top 100 keeps the stored bundle small.
  if (nrow(r) > 100) r <- r[seq_len(100), , drop = FALSE]
  rownames(r) <- NULL
  r
}

pairs <- list()
for (i in seq_len(length(coefs) - 1)) {
  for (j in (i + 1):length(coefs)) {
    a <- coefs[i]; b <- coefs[j]
    key <- paste(sort(c(a, b)), collapse = "|")
    ord <- sort(c(a, b)); a <- ord[1]; b <- ord[2]   # canonical order
    message("ORA: ", a, " vs ", b, " ...")
    seg <- list(); counts <- list(up = c(), down = c())
    for (dir in c("up", "down")) {
      A <- sig_sets[[a]][[dir]]; B <- sig_sets[[b]][[dir]]
      a_only <- setdiff(A, B); b_only <- setdiff(B, A); common <- intersect(A, B)
      counts[[dir]] <- c(a_only = length(a_only), b_only = length(b_only), common = length(common))
      seg$a_only[[dir]] <- run_ora(a_only, universe)
      seg$b_only[[dir]] <- run_ora(b_only, universe)
      seg$common[[dir]] <- run_ora(common, universe)
    }
    pairs[[key]] <- list(a = a, b = b,
                         a_label = unname(coef_choices[coef_choices == a]),
                         b_label = unname(coef_choices[coef_choices == b]),
                         counts = counts, seg = seg)
  }
}

# Diagnostic: total enriched GO rows actually stored across all pairs / segments.
# 0 means enrichment produced nothing for everything (e.g. org.Mm.eg.db not attached,
# or the gene IDs are not Ensembl) - the app would then show "No enriched GO terms for
# this set" for every combination. Fix the cause and regenerate rather than shipping
# an empty bundle.
n_terms <- sum(vapply(pairs, function(p)
  sum(vapply(p$seg, function(s) nrow(s$up) + nrow(s$down), integer(1))), integer(1)))
message("  total enriched GO rows stored: ", n_terms)
if (n_terms == 0)
  message("  WARNING: no GO terms for any set - check that org.Mm.eg.db is attached and gene IDs are Ensembl.")

out <- list(
  coef_choices = coef_choices,
  pairs        = pairs,
  fdr          = fdr,
  ontology     = "BP",
  built        = as.character(Sys.Date())
)
out_file <- file.path(out_dir, "rnaseq_ora.rds")
saveRDS(out, out_file)

# Tiny UI sidecar: only the coefficient dropdown choices, so rnaseq_dea_ui() can
# build its ORA controls without decompressing the full ORA bundle.
saveRDS(list(coef_choices = coef_choices, fdr = fdr, built = out$built),
        file.path(out_dir, "rnaseq_ora_meta.rds"))

message("\nDone. Saved ", out_file)
message("  coefficients: ", paste(coefs, collapse = ", "))
message("  pairs: ", length(pairs), " at FDR < ", fdr)
