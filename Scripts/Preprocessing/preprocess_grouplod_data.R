#!/usr/bin/env Rscript
# Pre-process Group LOD data: convert TSV to RDS, pre-join with phenotype map,
# and split by group for fast per-group loading in the Shiny app.
#
# Run once:
#   Rscript Scripts/preprocess_grouplod_data.R

lod_dir <- "Data/LOD"
out_dir <- file.path(lod_dir, "rds")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Read phenotype map
map_file <- file.path(lod_dir, "fix_interst_phenotype_GB")
mp <- read.table(map_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mp$Raw_name    <- as.character(mp$Raw_name)
mp$Simply_name <- as.character(mp$Simply_name)
mp$Group       <- as.character(mp$Group)

for (sex in c("male", "female")) {
  fname <- if (sex == "male") "add_male.merge005001" else "add_female.merge005001"
  message("Reading ", fname, " ...")
  qtl <- read.table(file.path(lod_dir, fname), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, check.names = FALSE)
  qtl$Chr       <- as.character(qtl$Chr)
  qtl$Pos_bp    <- as.numeric(qtl$Pos_bp)
  qtl$Lod       <- as.numeric(qtl$Lod)
  qtl$x_mb      <- qtl$Pos_bp / 1e6
  qtl$Phenotype <- as.character(qtl$Phenotype)
  qtl$Status    <- as.character(qtl$Status)

  # Join with phenotype map
  df <- merge(qtl, mp, by.x = "Phenotype", by.y = "Raw_name", all.x = FALSE)
  df$display_pheno <- ifelse(is.na(df$Simply_name) | df$Simply_name == "",
                             df$Phenotype, df$Simply_name)
  df <- df[is.finite(df$x_mb) & is.finite(df$Lod) & !is.na(df$Chr) &
           !is.na(df$display_pheno) & df$display_pheno != "", ]

  # Split by group and save each as a separate RDS
  groups <- unique(df$Group)
  for (g in groups) {
    df_g <- df[df$Group == g, ]
    out_file <- file.path(out_dir, paste0(sex, "_", gsub("[^A-Za-z0-9_]", "_", g), ".rds"))
    saveRDS(df_g, out_file)
    message("  Saved: ", out_file, " (", nrow(df_g), " rows)")
  }

  # Also save the full joined dataset as one RDS per sex (fallback)
  full_file <- file.path(out_dir, paste0(sex, "_full.rds"))
  saveRDS(df, full_file)
  message("  Saved full: ", full_file, " (", nrow(df), " rows)")
}

# Save group list
saveRDS(sort(unique(mp$Group)), file.path(out_dir, "group_list.rds"))
# Save phenotype map
saveRDS(mp, file.path(out_dir, "phenotype_map.rds"))

message("\nDone. RDS files saved to: ", out_dir)
