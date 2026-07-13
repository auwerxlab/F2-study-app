#!/usr/bin/env Rscript
# Preprocess QTL mapping data: convert raw TSV files to RDS for the Shiny app.
# Usage: Rscript Scripts/preprocess_qtl_data.R

data_dir <- "Data/Figure1_all_strain"

# --- Read chromosome lengths ---
chrom_df <- read.table(file.path(data_dir, "chr_GRCm38.length"),
                       header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chrom_df) <- c("chr", "chr_length")

chrom_df$chr <- as.numeric(chrom_df$chr)
chrom_df <- chrom_df[order(chrom_df$chr), ]
chrom_df$chr <- as.character(chrom_df$chr)

cum_lengths <- cumsum(as.numeric(chrom_df$chr_length))
chrom_df$cum_start <- c(0, head(cum_lengths, -1))
chrom_df$cum_mid   <- chrom_df$cum_start + as.numeric(chrom_df$chr_length) / 2

saveRDS(chrom_df, file.path(data_dir, "chr_GRCm38.rds"))
cat("Saved chr_GRCm38.rds\n")

# --- Helper: read and annotate one QTL file ---
process_qtl <- function(filepath, chrom_df) {
  df <- read.table(filepath, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "pos", "lod", "lod_threshold", "phenotype", "significance")
  df$chr <- as.character(df$chr)

  df <- merge(df, chrom_df[, c("chr", "cum_start")], by = "chr", all.x = TRUE)
  df$global_pos <- as.numeric(df$pos) + df$cum_start

  df$sig_cat <- ifelse(df$significance == "non_significant", "NS",
                ifelse(df$significance == "Significant_005", "0.01<=P<=0.05",
                ifelse(df$significance == "Significant_001", "P<=0.01", NA)))
  df
}

# --- Process male and female ---
pheno_all <- character(0)
for (sex in c("male", "female")) {
  infile  <- file.path(data_dir, paste0(sex, "_005_001"))
  outfile <- file.path(data_dir, paste0(sex, "_005_001.rds"))
  cat("Processing", infile, "...\n")
  df <- process_qtl(infile, chrom_df)
  saveRDS(df, outfile)
  pheno_all <- union(pheno_all, df$phenotype)
  cat("Saved", outfile, "(", nrow(df), "rows )\n")
}

# Tiny UI sidecar: the sorted phenotype names, so qtl_ui() can populate its
# phenotype selector without decompressing a ~4.6 MB mapping RDS just to list them.
saveRDS(sort(unique(pheno_all)), file.path(data_dir, "qtl_pheno_choices.rds"))
cat("Saved qtl_pheno_choices.rds (", length(pheno_all), "phenotypes )\n")

cat("Done.\n")
