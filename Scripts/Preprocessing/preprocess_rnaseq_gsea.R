#!/usr/bin/env Rscript
# Pre-process GO GSEA for the RNA-seq DEA tab.
#
# Reads the precomputed DEA contrast tables (Data/RNAseq/rnaseq_dea.rds) and runs
# GO gene-set enrichment (gseGO, biological process) per contrast, ranking genes
# by sign(logFC) * -log10(adj.P.Val) as in Scripts/RNA-seq/04-RNAseq_DEA.R.
# Saves the enrichment tables; the app draws the dot plot + table at runtime.
#
# No limma refit needed -- it reuses the DEA tables. Needs clusterProfiler +
# org.Mm.eg.db (gene IDs are Ensembl mouse).
#
# Run after preprocess_rnaseq_dea.R:
#   Rscript Scripts/preprocess_rnaseq_gsea.R [dea_rds] [out_dir]

args     <- commandArgs(trailingOnly = TRUE)
dea_path <- if (length(args) >= 1) args[1] else "Data/RNAseq/rnaseq_dea.rds"
out_dir  <- if (length(args) >= 2) args[2] else "Data/RNAseq"

if (!file.exists(dea_path)) {
  stop("rnaseq_dea.rds not found at '", dea_path, "'. Run preprocess_rnaseq_dea.R first.")
}
for (pkg in c("clusterProfiler", "org.Mm.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
}
# gseGO() below is called with `OrgDb = org.Mm.eg.db` (a bare symbol), so the package
# must be ATTACHED, not merely loaded: requireNamespace() alone leaves that object
# unresolved, so gseGO() errors with "object 'org.Mm.eg.db' not found" for every
# contrast and no GSEA results are stored.
suppressPackageStartupMessages(library(org.Mm.eg.db))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dea <- readRDS(dea_path)
contrasts <- dea$contrasts

# Ranked gene list for GSEA: sign(logFC) * -log10(adj.P), Inf capped, sorted desc.
make_ranking <- function(tbl) {
  r <- sign(tbl$logFC) * -log10(tbl$adj.P.Val)
  names(r) <- tbl$gene_id
  r <- r[!is.na(r) & names(r) != "" & !is.na(names(r))]
  fin <- r[is.finite(r)]
  hi <- if (length(fin)) max(fin) else 1
  lo <- if (length(fin)) min(fin) else -1
  r[r ==  Inf] <- hi + abs(hi) * 0.1 + 1
  r[r == -Inf] <- lo - abs(lo) * 0.1 - 1
  r <- r[!duplicated(names(r))]
  sort(r, decreasing = TRUE)
}

gsea <- list()
for (key in names(contrasts)) {
  message("GSEA: ", key, " ...")
  rank <- make_ranking(contrasts[[key]]$table)
  if (length(rank) < 10) { message("  too few genes, skipping"); next }

  gs <- tryCatch(
    clusterProfiler::gseGO(geneList = rank, ont = "BP", OrgDb = org.Mm.eg.db,
                           keyType = "ENSEMBL", pvalueCutoff = 1, verbose = FALSE),
    error = function(e) { message("  gseGO failed: ", conditionMessage(e)); NULL })
  if (is.null(gs) || nrow(gs@result) == 0) { message("  no enrichment"); next }

  res <- gs@result
  res$n_leading <- vapply(strsplit(res$core_enrichment, "/"), length, integer(1))
  res$GeneRatio <- res$n_leading / res$setSize
  keep <- c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust",
            "n_leading", "GeneRatio", "core_enrichment")
  res <- res[, intersect(keep, names(res)), drop = FALSE]
  rownames(res) <- NULL

  gsea[[key]] <- list(label = contrasts[[key]]$label, result = res)
  message("  ", nrow(res), " terms (", sum(res$p.adjust < 0.05, na.rm = TRUE), " at p.adjust<0.05)")
}

out <- list(gsea = gsea, ontology = "BP", built = as.character(Sys.Date()))
out_file <- file.path(out_dir, "rnaseq_gsea.rds")
saveRDS(out, out_file)
message("\nDone. Saved ", out_file)
message("  contrasts with GSEA: ", paste(names(gsea), collapse = ", "))
