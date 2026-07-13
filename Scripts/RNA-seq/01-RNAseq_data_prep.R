rm(list = ls())
library(tidyverse)
library(cowplot)
library(edgeR)
library(dplyr)
source("./Scripts/00-common.R")

# prepare txdb, useful for later
library(GenomicFeatures)
txdb <- GenomicFeatures::makeTxDbFromGFF("./Data/Mus_musculus.GRCm38.92.gtf")
AnnotationDbi::saveDb(txdb, file = "./Mus_musculus.GRCm38.92.txdb")
system("mv ./Mus_musculus.GRCm38.92.txdb ./Data/Mus_musculus.GRCm38.92.txdb")

gtf <- rtracklayer::import('./Data/Mus_musculus.GRCm38.92.gtf')
gtf_df <- as.data.frame(gtf)

write.csv(gtf_df, file = "./Data/Mus_musculus.GRCm38.92.gtf.df.csv")

geneConversionTable <- gtf_df %>% dplyr::select(seqnames, type, gene_id, gene_name, gene_biotype)
ung <- geneConversionTable$gene_id %>% unique()
geneConversionTable <- geneConversionTable[match(ung, geneConversionTable$gene_id),]
write.csv(geneConversionTable, file = "./Data/geneConversionTable.csv")


###########
count_files <- list.files(path = "./RNAseq_liver_scitas/output/aligned/", recursive = T, pattern = "ReadsPerGene.out.tab", full.names = T)
length(count_files)
samples <- gsub("(.)+\\/", "", gsub("\\.ReadsPerGene.out.tab", "", count_files))
## two samples are mislabled: 1204 --> 2004 ; 5481 --> 5841
samples <- plyr::mapvalues(samples, from = c("1204", "5481"), to = c("2004", "5841"))
## there is a switch between two samples (based on later hierarchical clustering and genotyping): 2001 (129S1, male) <-> 4671 (F2, female)
w1 <- which(samples == "2001")
w2 <- which(samples == "4671")

samples[w1] <- "4671"
samples[w2] <- "2001"

liver_dge <- readDGE(files = count_files, columns = c(1,4),labels = samples)
dim(liver_dge)
liver_dge$samples$Sample_ID <- samples

## add gene conversion table
geneConversionTable <- read.csv(file = "./Data/geneConversionTable.csv")
liver_dge$geneConversionTable <- geneConversionTable

## get the whole phenotype matrix
pheno_melt <- read.csv("./Data/allData_melt_clean_noOutlier.csv")
pheno_melt$Strain <- factor(pheno_melt$Strain, levels = names(strainColors))
pheno_melt$Sample_ID <- substr(pheno_melt$ID, 7, 10)

pheno_wide <- read.csv("./Data/allData_wide_clean_noOutlier.csv")
pheno_wide$Strain <- factor(pheno_wide$Strain, levels = names(strainColors))
pheno_wide$Sample_ID <- substr(pheno_wide$ID, 7, 10)


liver_dge$samples <- left_join(liver_dge$samples, pheno_wide, by = "Sample_ID")

liver_dge <- calcNormFactors(liver_dge)
keep <- edgeR::filterByExpr(liver_dge, design = model.matrix(~Strain + Sex, data = liver_dge$samples))
table(keep)
liver_dge$counts_filtered <- liver_dge$counts[keep,]
liver_dge$cpm <- cpm(liver_dge$counts_filtered, log = T)
liver_dge$y <- voom(liver_dge$counts_filtered, design = model.matrix(~Strain + Sex, data = liver_dge$samples))

saveRDS(liver_dge, file = "./Data/liver_dge.RDS")

