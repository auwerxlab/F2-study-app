rm(list = ls())
library(tidyverse)
library(cowplot)
library(edgeR)
library(FactoMineR)
library(ComplexHeatmap)
library(dendextend)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(plotly)
library(ggrepel)
source("./Scripts/00-common.R")
################################# Data loading ###################################

liver_dge <- readRDS("./Data/liver_dge.RDS")

########################## Z score and color functions ############################

scale2 <- function(x, na.rm = T) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
color.gradient <- function(x, colors=colors, colsteps=colsteps) {
  return( colorRampPalette(colors) (colsteps) 
          [ findInterval(x, seq(-max(abs(x),na.rm = T),max(abs(x), na.rm = T), length.out=colsteps)) ] 
  )
}
scale_and_color <- function(x, na.rm = T, colors = c("blue", "white", "red"), colsteps = 256) {
  x<-color.gradient(x=scale2(x), colors=colors, colsteps = colsteps)
  x[is.na(x)]<-"grey70"
  return(x)
}

################ Hierarchical clustering #####################
### Hierarchical clustering ###
hc <- hclust(dist(liver_dge$cpm %>% t()))
hcm <- hclust(dist(liver_dge$cpm[,liver_dge$samples$Sex == "m"] %>% t()))
hcf <- hclust(dist(liver_dge$cpm[,liver_dge$samples$Sex == "f"] %>% t()))

dend <- as.dendrogram(hc) %>% reorder(wts = liver_dge$samples$TIMP1_24)
dendm <- as.dendrogram(hcm) %>% reorder(wts = liver_dge$samples[liver_dge$samples$Sex == "m",]$TIMP1_24)
dendf <- as.dendrogram(hcf) %>% reorder(wts = liver_dge$samples[liver_dge$samples$Sex == "f",]$TIMP1_24)

#new
dend <- as.dendrogram(hc) %>% reorder(wts = liver_dge$samples$Fibrosis_perc_24)
dendm <- as.dendrogram(hcm) %>% reorder(wts = liver_dge$samples[liver_dge$samples$Sex == "m",]$Fibrosis_perc_24)
dendf <- as.dendrogram(hcf) %>% reorder(wts = liver_dge$samples[liver_dge$samples$Sex == "f",]$Fibrosis_perc_24)


# Export the dendrogram to compare to the kidney
#save(hc,dend,file="./Data/Liver_dendrogram.Rdata")

# selected_phenos <- c(BW_perGain_24,"Rel. Fat gain", "Liver_g_RelBWSac_24", "ALAT_24", "TIMP.1")
liver_phenos <- liver_dge$samples %>%
  mutate('Plasma TIMP-1' = log(TIMP1_24),
         'Plasma GDF-15' = log(GDF15_24),
         'Plasma KIM-1' = log(KIM1_24),
         'Plasma ALT' = log(ALAT_24),
         'BW gain(%)' = BW_perGain_24,
         'Liver Weight(g/g BW)' = Liver_g_RelBWSac_24,
         'Steatosis'= log(sum_all_vacuoles_percentage_24),
         "Fibrosis"=log(Fibrosis_perc_24),
         'Fat mass(%)' = Fat_perc_22) %>%
  dplyr::select('Plasma TIMP-1', 'Plasma ALT', 'BW gain(%)', 'Fat mass(%)', 'Liver Weight(g/g BW)', 'Steatosis',"Fibrosis")


# liver_phenos_scaled <- liver_phenos %>% mutate_all(., scale2)
liver_phenos_colors <- liver_phenos %>% mutate_all(., scale_and_color)
#liver_phenos_test_colors <- liver_phenos_test %>% mutate_all(., scale_and_color)

pdf(file = "./Plots/14-liver_hclust_all.pdf", useDingbats = F, width = 15,height = 5)
par(mar = c(10,6.1,4.1,2.1))
dend %>% color_branches(col = strainColors[liver_dge$samples$Strain[match(hc$labels[order.dendrogram(dend)], liver_dge$samples$Sample_ID)]]) %>% 
  plot(horiz=F, axes=FALSE, leaflab ="none", main = "All liver transcriptomes")
lbls <- cbind(liver_phenos_colors,
              "Sex" = sexColors[liver_dge$samples$Sex],
              "Strain" = strainColors[liver_dge$samples$Strain])
colored_bars(colors = lbls, dend = dend, rowLabels = colnames(lbls), y_shift = -5, cex.rowLabels = 0.8)

dev.off()

#pdf(file = "./Plots/Liver/Giorgia/GB02-liver_hclust_males.pdf", useDingbats = F, width = 15,height = 5)
#par(mar = c(10,6.1,4.1,2.1))
#dendm %>% color_branches(col = strainColors[liver_dge$samples$Strain[match(hcm$labels[order.dendrogram(dendm)], liver_dge$samples$Sample_ID)]]) %>% 
#  plot(horiz=F, axes=FALSE, leaflab ="none", main = "Male liver transcriptomes")
#lbls <- cbind(liver_phenos_colors[match(hcm$labels, liver_dge$samples$Sample_ID),],
#              "Sex" = sexColors[liver_dge$samples$Sex[match(hcm$labels, liver_dge$samples$Sample_ID)]], 
#              "Strain" = strainColors[liver_dge$samples$Strain[match(hcm$labels, liver_dge$samples$Sample_ID)]])
#colored_bars(colors = lbls, dend = dendm, rowLabels =  colnames(lbls), y_shift = -5, cex.rowLabels = 0.8)
#dev.off()

#pdf(file = "./Plots/Liver/02-liver_hclust_females.pdf", useDingbats = F, width = 15,height = 5)
#par(mar = c(10,6.1,4.1,2.1))
#dendf %>% color_branches(col = strainColors[liver_dge$samples$Strain[match(hcf$labels[order.dendrogram(dendf)], liver_dge$samples$Sample_ID)]]) %>% 
#  plot(horiz=F, axes=FALSE, leaflab ="none", main = "Female liver transcriptomes")
#lbls <- cbind(liver_phenos_colors[match(hcf$labels, liver_dge$samples$Sample_ID),],
#              "Sex" = sexColors[liver_dge$samples$Sex[match(hcf$labels, liver_dge$samples$Sample_ID)]], 
#              "Strain" = strainColors[liver_dge$samples$Strain[match(hcf$labels, liver_dge$samples$Sample_ID)]])
#colored_bars(colors = lbls, dend = dendf, rowLabels =  colnames(lbls), y_shift = -5, cex.rowLabels = 0.8)
#dev.off()

### a plot with mitochondrial respiration genes

GOCC <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC")
GOBP <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
genesets<-unique(GOCC[grep("MITOCHON",GOCC$gs_name),]$gs_name)
genesets_oxPhos<-unique(GOBP[grep("GOBP_OXIDATIVE_PHOSPHORYLATION",GOBP$gs_name),]$gs_name)
selected_gene_names<-unique(GOCC[grep("GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX",GOCC$gs_name),]$gene_symbol)
selected_gene_names_oxphos<-unique(GOBP[grep("GOBP_OXIDATIVE_PHOSPHORYLATION",GOBP$gs_name),]$gene_symbol)
# oxphos <- readxl::read_excel("./Data/OxphosGenes.xlsx", skip = 1)
# oxphos <- oxphos %>% mutate(gene_name = gsub("^ND", "mt-ND", `Gene Symbol`))
# tt_oxphos <- tt[match(oxphos$gene_name, tt$gene_name),]
# tt_oxphos %>% filter(!is.na(logFC)) %>% arrange(adj.P.Val)
selected_genes <- liver_dge$geneConversionTable[match(selected_gene_names_oxphos,liver_dge$geneConversionTable$gene_name),]
tt<-as.data.frame(liver_dge$y$E[-(1:3),])
selected_genes<-selected_genes[!is.na(selected_genes$gene_id)&selected_genes$gene_id %in% row.names(tt),]

selected_genes_exp <- tt[selected_genes$gene_id ,] %>% t() %>% as.data.frame()
#selected_genes_exp <- liver_dge$y$E[selected_genes$gene_id %>% as.character(),] %>% t() %>% as.data.frame()
colnames(selected_genes_exp) <- paste0("Liver tx: ", selected_genes$gene_name)

selected_genes_exp <- selected_genes_exp %>% mutate_all(., scale_and_color)

pdf(file = "./Plots/14-liver_hclust_oxphos.pdf", useDingbats = F, width = 15,height = 170)
par(mar = c(800,6.1,4.1,2.1))
dend  %>% color_branches(col = strainColors[liver_dge$samples$Strain[match(hc$labels[order.dendrogram(dend)], liver_dge$samples$Sample_ID)]]) %>%
  plot(horiz=F, axes=FALSE, leaflab ="none", main = "All liver transcriptomes")
lbls <- cbind(selected_genes_exp,
              liver_phenos_colors,
              "Sex" = sexColors[liver_dge$samples$Sex],
              "Strain" = strainColors[liver_dge$samples$Strain])
colored_bars(colors = lbls, dend = dend, rowLabels = colnames(lbls), y_shift = 0, cex.rowLabels = 0.8)
dev.off()

### a plot with some "key" liver gene expression included

selected_gene_names <-  c("Uqcrc1","Atp5d","Cox5a","Ndufb9", "Nfkb1", "Il1b", "Col3a1", "Col4a1","Timp1","Gdf15")
selected_genes <- liver_dge$geneConversionTable[match(selected_gene_names,liver_dge$geneConversionTable$gene_name),]
selected_genes<-selected_genes[!is.na(selected_genes$gene_id)&selected_genes$gene_id %in% row.names(tt),]
selected_genes_exp <- liver_dge$y$E[selected_genes$gene_id %>% as.character(),] %>% t() %>% as.data.frame()
colnames(selected_genes_exp) <- paste0("Liver tx: ", selected_genes$gene_name)

liver_dge$geneConversionTable[grep("mt-",liver_dge$geneConversionTable$gene_name),]$gene_name

selected_genes_exp <- selected_genes_exp %>% mutate_all(., scale_and_color)
pdf(file = "./Plots/14-liver_hclust_all_withexpression_mito.pdf", useDingbats = F, width = 8,height = 4)
par(mar = c(12,6.1,4.1,2.1))
dend  %>% color_branches(col = strainColors[liver_dge$samples$Strain[match(hc$labels[order.dendrogram(dend)], liver_dge$samples$Sample_ID)]]) %>%
  plot(horiz=F, axes=FALSE, leaflab ="none", main = "All liver transcriptomes")
lbls <- cbind(selected_genes_exp,
              liver_phenos_colors,
              "Sex" = sexColors[liver_dge$samples$Sex],
              "Strain" = strainColors[liver_dge$samples$Strain])
colored_bars(colors = lbls, dend = dend, rowLabels = colnames(lbls), y_shift = 0, cex.rowLabels = 0.8)
dev.off()

### a plot with some signature genes from recent papers (PMID:31398325, PMID:37075704, PMID:37037945,PMID:40310463, PMID:38478630)

selected_gene_names <-  c("Clstn1","Pdgfra", "Vim","Vwf" ,"Jag1",
                          "Trem2","Gpnmb",
                          #"Cd9",
                          "Ccl2","Col1a1", "Timp1",
                          "Adamtsl2",
                          #"Akr1b10",
                          #"Cfhr4",not found
                          #"Postn",
                          #"Ms4a7", 
                          "Lgals3", 
                          #"Hsd17b13",
                          #"Pnpla3", 
                          "Ctsd")
selected_genes <- liver_dge$geneConversionTable[match(selected_gene_names,liver_dge$geneConversionTable$gene_name),]
selected_genes<-selected_genes[!is.na(selected_genes$gene_id)&selected_genes$gene_id %in% row.names(tt),]
selected_genes_exp <- liver_dge$y$E[selected_genes$gene_id %>% as.character(),] %>% t() %>% as.data.frame()
colnames(selected_genes_exp) <- paste0("Liver tx: ", selected_genes$gene_name)

liver_dge$geneConversionTable[grep("mt-",liver_dge$geneConversionTable$gene_name),]$gene_name

selected_genes_exp <- selected_genes_exp %>% mutate_all(., scale_and_color)
pdf(file = "./Plots/14-liver_hclust_all_withexpression_b.pdf", useDingbats = F, width = 8,height = 12)
par(mar = c(12,6.1,4.1,2.1))
dend  %>% color_branches(col = strainColors[liver_dge$samples$Strain[match(hc$labels[order.dendrogram(dend)], liver_dge$samples$Sample_ID)]]) %>%
  plot(horiz=F, axes=FALSE, leaflab ="none", main = "All liver transcriptomes")
lbls <- cbind(selected_genes_exp,
              liver_phenos_colors,
              "Sex" = sexColors[liver_dge$samples$Sex],
              "Strain" = strainColors[liver_dge$samples$Strain])
colored_bars(colors = lbls, dend = dend, rowLabels = colnames(lbls), y_shift = 0, cex.rowLabels = 0.8)
dev.off()

# liver_dge$samples[order.dendrogram(dend),] %>% pull(Strain) %>% head()
# tmp <- liver_dge$samples[order.dendrogram(dend),] %>% select(Strain, Sex, Sample_ID)

