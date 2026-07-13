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

mds <- plotMDS(liver_dge$counts_filtered, gene.selection = "common", top = 500)


############## PCA ################
pca <- FactoMineR::PCA(liver_dge$cpm %>% t(), graph = F)
pca$ind$coord <- cbind(pca$ind$coord %>% as.data.frame(), liver_dge$samples)
pca$ind$coord$orderrank <- rank(pca$ind$coord$TIMP1_24,ties.method="first")
pca$ind$coord <- pca$ind$coord[order(pca$ind$coord$TIMP1_24, decreasing=F), ]

ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = Strain, shape = Sex)) + 
  geom_point(aes(shape=Sex
                 #,size = TIMP1_24
  ))+
  scale_color_manual(values = strainColors) +
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(filename = "./Plots/13-liver_dge_PCA_1_2.pdf", width =4.3, height = 2.5, units = c("in"), useDingbats = F)

summary(pca$ind$coord$BW_perGain_24)

bw <- ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = BW_perGain_24
                                #, shape = Sex
                                #,order=orderrank
)) + 
  geom_point(size = 1)+
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
  )+
  labs(color = NULL)+
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  theme_cowplot()+
  ggtitle("Body weight(%gain)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
bw

summary(pca$ind$coord$Liver_g_24)
lg <- ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = log1p(Liver_g_24)
                                #, shape = Sex
                                #,order=orderrank
)) + 
  geom_point(size = 1)+
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
  )+
  labs(color = NULL)+
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  theme_cowplot()+
  ggtitle("Liver weight(g)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
lg

summary(pca$ind$coord$Liver_g_RelBWSac_24)
lr <- ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = Liver_g_RelBWSac_24
                                #, shape = Sex
                                #,order=orderrank
)) + 
  geom_point(size = 1)+
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
  )+
  labs(color = NULL)+
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  theme_cowplot()+
  ggtitle("Liver weight(g/g BW)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
lr

summary(pca$ind$coord$TIMP1_24)
tmp <- ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = log1p(TIMP1_24)
                                 #, shape = Sex
                                 #,order=orderrank
)) + 
  geom_point(size = 1)+
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
  )+
  labs(color = NULL)+
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  theme_cowplot()+
  ggtitle("Plasma TIMP1")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
tmp

summary(pca$ind$coord$Fibrosis_perc_24)
fp <- ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = log1p(Fibrosis_perc_24)
                                #, shape = Sex
                                #,order=orderrank
)) + 
  geom_point(size = 1)+
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
  )+
  labs(color = NULL)+
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  theme_cowplot()+
  ggtitle("Fibrosis")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
fp

lims <- quantile(pca$ind$coord$ALAT_24, c(0.01, 0.99), na.rm=T) 
mid_val <- median(range(pca$ind$coord$ALAT_24, na.rm=T))

summary(pca$ind$coord$ALAT_24)
alt <- ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = log1p(ALAT_24)
                                 #, shape = Sex
                                 #,order=orderrank
)) + 
  geom_point(size = 1)+
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
  )+
  labs(color = NULL)+
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  theme_cowplot()+
  ggtitle("ALT")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
alt

summary(pca$ind$coord$sum_all_vacuoles_percentage_24)

lf <- ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = log1p(sum_all_vacuoles_percentage_24)
                                #, shape = Sex
                                #,order=orderrank
)) + 
  geom_point(size = 1)+
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
  )+
  labs(color = NULL)+
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  theme_cowplot()+
  ggtitle("Steatosis")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
lf 
ggsave(filename = "./Plots/13-liver_dge_PCA_1_2_v2.pdf", width =6.5, height = 4, units = c("in"), useDingbats = F)

all <- plot_grid(bw, lg,tmp,alt, lf,fp, ncol=2)
ggsave(all,file="./Plots/13-transcriptome_PCA_pheno.pdf",
       width = 12, height = 13, units = c("cm"))  

ggplot(pca$ind$coord, aes(x = Dim.3, y = Dim.4, color = Strain, shape = Sex)) + 
  geom_point()+
  scale_color_manual(values = strainColors) +
  xlab(paste0("PC3 (", round(pca$eig[3,2],2),"%)")) +
  ylab(paste0("PC4 (", round(pca$eig[4,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
ggsave(filename = "./Plots/13-liver_dge_PCA_3_4.pdf", width =4.3, height = 2.5, units = c("in"), useDingbats = F)

ggplot(pca$ind$coord, aes(x = Dim.3, y = Dim.4, color = TIMP1_24, shape = Sex)) + 
  geom_point(aes(shape=Sex
                 ,size = Liver_g_RelBWSac_24
  ))+
  #scale_color_manual(values = strainColors) +
  scale_color_gradient(low = "yellow", high = "darkblue")+
  xlab(paste0("PC3 (", round(pca$eig[3,2],2),"%)")) +
  ylab(paste0("PC4 (", round(pca$eig[4,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(filename = "./Plots/13-liver_dge_PCA_3_4.pdf", width =4.3, height = 2.5, units = c("in"), useDingbats = F)



ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = Sex)) + 
  geom_point()+
  scale_color_manual(values = sexColors) +
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(filename = "./Plots/13-liver_dge_PCA_m_f_1_2.pdf", width =3.3, height = 2.5, units = c("in"), useDingbats = F)

ggplot(pca$ind$coord, aes(x = Dim.1, y = Dim.2, color = Sex
                          #, shape = Sex
)) + 
  geom_point(aes(
    #shape=Sex,
    size = sum_all_vacuoles_percentage_24
  ))+
  scale_color_manual(values = sexColors) +
  #scale_color_gradient(low = "yellow", high = "darkblue")+
  #scale_color_gradient(low = "darkblue", high = "orange")+
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(filename = "./Plots/13-liver_dge_PCA_1_2_sex_v2.pdf", width =6.5, height = 4, units = c("in"), useDingbats = F)


ggplot(pca$ind$coord, aes(x = Dim.3, y = Dim.4, color = Sex)) + 
  geom_point()+
  scale_color_manual(values = sexColors) +
  xlab(paste0("PC3 (", round(pca$eig[3,2],2),"%)")) +
  ylab(paste0("PC4 (", round(pca$eig[4,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(filename = "./Plots/13-liver_dge_PCA_m_f_3_4.pdf", width =3.3, height = 2.5, units = c("in"), useDingbats = F)

library("factoextra")
contrib <- fviz_contrib(pca, choice = "var", axes = 1, top = 300,
                        #palette = phenoColors
)
contrib

############## PCA males and females separately################
cpm <- liver_dge[["cpm"]]
samples <- liver_dge[["samples"]]
samples_m <- samples[grepl("m", samples$Sex),]
sample_ID_m <- samples_m$Sample_ID
cpm_m <- cpm[, colnames(cpm) %in% sample_ID_m]
pca_m <- FactoMineR::PCA(cpm_m %>% t(), graph = F)

pca_m$ind$coord <- cbind(pca_m$ind$coord %>% as.data.frame(), samples_m)

samples_f <- samples[grepl("f", samples$Sex),]
sample_ID_f <- samples_f$Sample_ID
cpm_f <- cpm[, colnames(cpm) %in% sample_ID_f]
pca_f <- FactoMineR::PCA(cpm_f %>% t(), graph = F)

pca_f$ind$coord <- cbind(pca_f$ind$coord %>% as.data.frame(), samples_f)

ggplot(pca_m$ind$coord, aes(x = Dim.1, y = Dim.2, color = Strain)) + 
  geom_point()+
  scale_color_manual(values = strainColors) +
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(filename = "./Plots/11-liver_dge_PCA_males_1_2.pdf", width =4.3, height = 2.5, units = c("in"), useDingbats = F)

ggplot(pca_m$ind$coord, aes(x = Dim.1, y = Dim.2, color = TIMP1_24)) + 
  geom_point(aes(size = sum_all_vacuoles_percentage_24))+
  scale_color_gradient(low = "yellow", high = "darkblue")+
  #scale_color_manual(values = strainColors) +
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size = 12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(filename = "./Plots/13-liver_dge_PCA_males_1_2_v2.pdf", width =7, height = 4, units = c("in"), useDingbats = F)


p1 <- fviz_contrib(pca_m, choice = "var", axes = 1, top = 40,
                   #palette = phenoColors
)
p1
contrib <- pca_m[["var"]][["contrib"]]
contrib <- as.data.frame(contrib)
contrib <- contrib[-c(1:3),]
contrib <- contrib[contrib$Dim.1>0.01,]
genes <- rownames(contrib)
list_genes <- list(genes)
names(list_genes) <- c("males")

go <- clusterProfiler::compareCluster(
  geneClusters = list_genes,
  fun = "enrichGO",
  # gene = x[,gene_id],
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = 'ENSEMBL',
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable=TRUE,
  minGSSize = 5,
  pvalueCutoff = 1,
) 

go_res <- go@compareClusterResult

ggplot(pca_f$ind$coord, aes(x = Dim.1, y = Dim.2, color = TIMP1_24)) + 
  geom_point(aes(size = sum_all_vacuoles_percentage_24))+
  scale_color_gradient(low = "yellow", high = "darkblue")+
  #scale_color_manual(values = strainColors) +
  xlab(paste0("PC1 (", round(pca$eig[1,2],2),"%)")) +
  ylab(paste0("PC2 (", round(pca$eig[2,2],2),"%)")) +
  # facet_grid(~Sex) +
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12),
        #axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size = 12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(filename = "./Plots/13-liver_dge_PCA_females_1_2_v2.pdf", width =7, height = 4, units = c("in"), useDingbats = F)



