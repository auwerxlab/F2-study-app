#!/usr/bin/env Rscript
# Pre-process the protein vs RNA-seq comparison for the Shiny app.
#
# Reproduces the core of Scripts/Protein_RNAseq_comparison/20-Prot_RNAseq_comparison.R:
# joins liver proteomics with RNA-seq (log2 protein vs log2 CPM per mouse/gene),
# and computes RNA-protein correlations per sample, per gene, and per gene within
# each liver cluster. Saves a compact bundle the app browses (distributions,
# per-gene table, gene-level scatter) without the raw omics matrices at runtime.
#
# Run once (where the proteomics + liver_dge data live; needs dplyr):
#   Rscript Scripts/preprocess_prot_rna.R [prots_rdata] [liver_dge] [gct_csv] [cluster_csv] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
prots_rdata <- if (length(args) >= 1) args[1] else "Data/Filtered_prots.Rdata"
dge_path    <- if (length(args) >= 2) args[2] else "Data/liver_dge.RDS"
gct_csv     <- if (length(args) >= 3) args[3] else "Data/geneConversionTable.csv"
cluster_csv <- if (length(args) >= 4) args[4] else "Data/allData_wide_clean_noOutlier_liver_cluster.csv"
out_dir     <- if (length(args) >= 5) args[5] else "Data/Comparison"

for (f in c(prots_rdata, dge_path, gct_csv, cluster_csv)) {
  if (!file.exists(f)) stop("Required input not found: ", f)
}
if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
suppressMessages({ library(dplyr) })
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Filtered_prots.Rdata provides (at least) Prot_filt_m and metadata.
e <- new.env()
load(prots_rdata, envir = e)
if (!exists("Prot_filt_m", envir = e)) stop("Filtered_prots.Rdata must contain 'Prot_filt_m'")
if (!exists("metadata", envir = e))    stop("Filtered_prots.Rdata must contain 'metadata'")
Prot_filt_m <- as.data.frame(e$Prot_filt_m)
metadata    <- as.data.frame(e$metadata)

liver_dge <- readRDS(dge_path)
gct <- read.csv(gct_csv)[, c("gene_id", "gene_name")]
colnames(gct) <- c("gene_id", "Gene.Symbol")
data_wide <- read.csv(cluster_csv)

strain_levels <- c("C57BL/6J", "129S1/SvImJ", "CAST/EiJ", "PWK/PhJ", "B6CASTF1", "129SPWKF1", "B6CAST-129SPWK-F2")
strain_colors <- setNames(c("#000000", "#FF9300", "#107F40", "#0000FE", "#FF2600", "#942192", "#C0C0C0"), strain_levels)

# --- Build long protein/RNA table (per mouse ID, per gene) ----------------
RNA <- as.data.frame(liver_dge$cpm)[-(1:3), ]          # drop first 3 (as in script)
RNA$gene_id <- rownames(RNA)
RNA_m <- reshape2::melt(RNA, id.vars = "gene_id", variable.name = "Sample_ID", value.name = "log2CPM")
RNA_m <- merge(gct, RNA_m)
RNA_m <- RNA_m[RNA_m$Gene.Symbol %in% Prot_filt_m$Gene.Symbol, ]
RNA_m <- merge(metadata[, c("Sample_ID", "ID")], RNA_m)

Prot_RNA <- merge(Prot_filt_m, RNA_m)
Prot_RNA <- Prot_RNA %>%
  dplyr::rename(log2CPM_RNA = log2CPM) %>%
  mutate(log2_Protein = log2(Prot_quantif)) %>%
  left_join(data_wide %>% dplyr::select(ID, liver_cluster), by = "ID")

# Keep one row per mouse x gene (max protein quantification)
Prot_RNA_unique <- Prot_RNA %>%
  group_by(ID, gene_id) %>%
  slice_max(order_by = Prot_quantif, n = 1, with_ties = FALSE) %>%
  ungroup()

# --- Correlations ---------------------------------------------------------
safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}
safe_cortest <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(c(est = NA_real_, p = NA_real_))
  tt <- try(suppressWarnings(cor.test(x[ok], y[ok], method = method)), silent = TRUE)
  if (inherits(tt, "try-error")) return(c(est = NA_real_, p = NA_real_))
  c(est = unname(tt$estimate), p = tt$p.value)
}

corr_per_sample <- Prot_RNA_unique %>%
  group_by(ID) %>%
  summarise(corr_pearson  = safe_cor(log2_Protein, log2CPM_RNA, "pearson"),
            corr_spearman = safe_cor(log2_Protein, log2CPM_RNA, "spearman"),
            n_genes = sum(is.finite(log2_Protein) & is.finite(log2CPM_RNA)),
            .groups = "drop") %>%
  left_join(metadata[, intersect(c("ID", "Sample_ID", "Strain", "Sex"), names(metadata))], by = "ID") %>%
  left_join(data_wide %>% dplyr::select(ID, liver_cluster), by = "ID")
corr_per_sample$Strain <- factor(corr_per_sample$Strain, levels = strain_levels)

corr_per_gene <- Prot_RNA_unique %>%
  group_by(gene_id) %>%
  summarise(corr_pearson  = safe_cor(log2_Protein, log2CPM_RNA, "pearson"),
            corr_spearman = safe_cor(log2_Protein, log2CPM_RNA, "spearman"),
            n_mice = sum(is.finite(log2_Protein) & is.finite(log2CPM_RNA)),
            .groups = "drop") %>%
  left_join(gct, by = "gene_id")

corr_per_gene_cluster <- Prot_RNA_unique %>%
  filter(!is.na(liver_cluster), is.finite(log2_Protein), is.finite(log2CPM_RNA)) %>%
  group_by(liver_cluster, gene_id) %>%
  filter(n() >= 5) %>%
  summarise(
    corr_pearson  = safe_cortest(log2_Protein, log2CPM_RNA, "pearson")["est"],
    p_pearson     = safe_cortest(log2_Protein, log2CPM_RNA, "pearson")["p"],
    corr_spearman = safe_cortest(log2_Protein, log2CPM_RNA, "spearman")["est"],
    p_spearman    = safe_cortest(log2_Protein, log2CPM_RNA, "spearman")["p"],
    n_mice = n(), .groups = "drop") %>%
  group_by(liver_cluster) %>%
  mutate(adj_pearson  = p.adjust(p_pearson, "BH"),
         adj_spearman = p.adjust(p_spearman, "BH")) %>%
  ungroup() %>%
  left_join(gct, by = "gene_id")

# Compact long table for the gene-level scatter
pairs <- Prot_RNA_unique %>%
  dplyr::select(gene_id, ID, log2_Protein, log2CPM_RNA, liver_cluster) %>%
  left_join(gct, by = "gene_id") %>%
  left_join(metadata[, intersect(c("ID", "Strain"), names(metadata))], by = "ID") %>%
  as.data.frame()

gene_syms <- corr_per_gene[!is.na(corr_per_gene$Gene.Symbol) & corr_per_gene$Gene.Symbol != "", ]
gene_choices <- setNames(gene_syms$gene_id, gene_syms$Gene.Symbol)
gene_choices <- gene_choices[order(names(gene_choices))]

cluster_levels <- sort(unique(stats::na.omit(Prot_RNA_unique$liver_cluster)))
cluster_palette <- c("#0696c7", "#FFD700", "#ff8800", "#6e0280", "#11999e", "#e84a5f")
cluster_colors <- setNames(cluster_palette[seq_along(cluster_levels)], as.character(cluster_levels))

out <- list(
  corr_per_sample       = as.data.frame(corr_per_sample),
  corr_per_gene         = as.data.frame(corr_per_gene),
  corr_per_gene_cluster = as.data.frame(corr_per_gene_cluster),
  pairs                 = pairs,
  gene_choices          = gene_choices,
  strain_levels         = strain_levels,
  strain_colors         = strain_colors,
  cluster_levels        = as.character(cluster_levels),
  cluster_colors        = cluster_colors,
  n_genes               = length(unique(Prot_RNA_unique$gene_id)),
  built                 = as.character(Sys.Date())
)
out_file <- file.path(out_dir, "prot_rna.rds")
saveRDS(out, out_file)

# Tiny UI sidecar: just the fields the Data-comparison UI needs to build its
# dropdowns, so compare_ui() can render controls without decompressing the full
# 82 MB bundle. Regenerated alongside prot_rna.rds so it never goes stale.
meta <- list(
  cluster_levels = as.character(cluster_levels),
  n_genes        = out$n_genes,
  built          = out$built
)
saveRDS(meta, file.path(out_dir, "prot_rna_meta.rds"))

message("Done. Saved ", out_file)
message("  samples: ", nrow(corr_per_sample), " | genes: ", nrow(corr_per_gene),
        " | clusters: ", paste(cluster_levels, collapse = ", "))
