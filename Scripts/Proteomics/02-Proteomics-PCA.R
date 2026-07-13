library(dplyr)
source("./Scripts/Proteomics/00-common.R")
load(file="./Aim2_Step2_manuscript/Data/Filtered_prots.Rdata")
liver_dge <- readRDS("./Aim2_Step2_manuscript/Data/liver_dge.RDS")



prot_matrix <- Prot_filt[,-c(2,3)]
# Remove rows with missing/empty Gene.Symbol
prot_matrix <- prot_matrix[
  !is.na(prot_matrix$Gene.Symbol) &
    prot_matrix$Gene.Symbol != "",
]

# Count NAs across data columns, excluding Gene.Symbol
data_cols <- setdiff(names(prot_matrix), "Gene.Symbol")
prot_matrix$na_count <- rowSums(is.na(prot_matrix[, data_cols]))

# For duplicated Gene.Symbols, keep row with fewer NAs
prot_matrix <- prot_matrix[
  order(prot_matrix$Gene.Symbol, prot_matrix$na_count),
]

prot_matrix <- prot_matrix[!duplicated(prot_matrix$Gene.Symbol), ]
#prot_matrix$Gene.Symbol[duplicated(prot_matrix$Gene.Symbol)]
# Remove helper column
prot_matrix$na_count <- NULL

# Set rownames and convert to matrix
rownames(prot_matrix) <- prot_matrix$Gene.Symbol

prot_matrix_mat <- as.matrix(
  prot_matrix[, !names(prot_matrix) %in% "Gene.Symbol"]
)
prot_matrix_mat_log <- log2(prot_matrix_mat + 1)
#range(prot_matrix_mat_log, na.rm = TRUE)
#sum(prot_matrix_mat_log == 0, na.rm = TRUE)
pca <- FactoMineR::PCA(prot_matrix_mat_log %>% t(), graph = F)

pca_coords <- cbind(
  ID = rownames(pca$ind$coord),
  as.data.frame(pca$ind$coord)
)

samples_pca <- samples %>%
  inner_join(pca_coords, by = "ID")



samples_pca$orderrank <- rank(samples_pca$TIMP1_24,ties.method="first")
samples_pca <- samples_pca[order(samples_pca$TIMP1_24, decreasing=F), ]


ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = Strain, shape = Sex)) + 
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
#ggsave(filename = "./Plots/22-liver_proteomics_PCA_1_2.pdf", width =4.3, height = 2.5, units = c("in"), useDingbats = F)

ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = Sex)) + 
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
#ggsave(filename = "./Plots/22-liver_proteomics_sex_PCA_1_2.pdf", width =3.3, height = 2.5, units = c("in"), useDingbats = F)


summary(samples_pca$BW_perGain_24)

bw <- ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = BW_perGain_24
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

summary(samples_pca$Liver_g_24)
lg <- ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = log1p(Liver_g_24)
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

summary(samples_pca$Liver_g_RelBWSac_24)
lr <- ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = Liver_g_RelBWSac_24
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

summary(samples_pca$TIMP1_24)
tmp <- ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = log1p(TIMP1_24)
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

summary(samples_pca$Fibrosis_perc_24)
fp <- ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = log1p(Fibrosis_perc_24)
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

lims <- quantile(samples_pca$ALAT_24, c(0.01, 0.99), na.rm=T) 
mid_val <- median(range(samples_pca$ALAT_24, na.rm=T))

summary(samples_pca$ALAT_24)
alt <- ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = log1p(ALAT_24)
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

summary(samples_pca$sum_all_vacuoles_percentage_24)

lf <- ggplot(samples_pca, aes(x = Dim.1, y = Dim.2, color = log1p(sum_all_vacuoles_percentage_24)
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
#ggsave(filename = "./Plots/22-liver_proteomics_PCA_1_2_v2.pdf", width =5, height = 4, units = c("in"), useDingbats = F)

all <- plot_grid(bw, lg,tmp,alt, lf,fp, ncol=2)
all
#ggsave(all,file="./Plots/22-proteomics_PCA_pheno.pdf",
#       width = 12, height = 13, units = c("cm"))
