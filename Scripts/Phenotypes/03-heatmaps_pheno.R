rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(readxl)
library(ggplot2)
library(janitor)
library(cowplot)
library(reshape2)
library(openxlsx)
library(lubridate)
library(DescTools)
library(multcompView)
library(ggcorrplot)
library(grid)
library(ComplexHeatmap)
library(tidyHeatmap)
source("./Scripts/Phenotypes/00-common.R")

data_melt <- read.csv(file = './Aim2_Step2_manuscript/Data/allData_melt_clean_noOutlier.csv')
data_wide <- read.csv(file = './Aim2_Step2_manuscript/Data/allData_wide_clean_noOutlier.csv')
dfplot_outlier <- read.csv(file = './Aim2_Step2_manuscript/Data/dfplot_outlier.csv')

col<- colorRampPalette(c("red", "white", "blue"))(256)
colours <- list('Strain' = strainColors,
                'Generation' = c('F0' = '#FFFFFF', 'F1' = '#818589',
                                 'F2' = "#000000"),
                "Sex" = sexColors)

########heatmap real values##########
########Fat%#########
dfPlot <- data_wide[,c(1,2,4, 28,32,36)]
dfPlot <- dfPlot[!is.na(dfPlot$Fat_perc_22),]
hm.identity <- as.data.frame(dfPlot[,c(1,2,3)])
rownames(hm.identity) <- hm.identity[,1]
hm.identity <- hm.identity[,-1]
hm.identity$Strain <- factor(hm.identity$Strain,levels = c("C57BL/6J","129S1/SvImJ", "CAST/EiJ","PWK/PhJ", 
                                                                           "B6CASTF1","129SPWKF1","B6CAST-129SPWK-F2"))
phenotype_names <- colnames(dfPlot)
#pheno_names <- read.csv(file = './Aim2_Step2_manuscript/Data/phenotype_names.csv')
#pheno_names <- pheno_names[-c(1:3),]

dfPlot$Strain_Sex <- paste0(dfPlot$Strain,"_",dfPlot$Sex)
dfPlot <- reshape2::melt(dfPlot, id.vars = c("ID", "Strain", "Sex", "Strain_Sex"))
dfPlot <- dfPlot[, c(1,5,6)]
dfPlot <- dfPlot %>% spread(ID, value)
colnames(dfPlot)[colnames(dfPlot)=="variable"] <- "Phenotype"
#dfPlot$Pheno_new <- pheno_names$Phenotype_new[match(dfPlot$Phenotype_old,pheno_names$Phenotype_old)]
df <- dfPlot[,-c(1)]
rownames(df) <- dfPlot$Phenotype
#df <- t(df)
df.matrix <- as.matrix(df)


#heatmap(df.matrix, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)
col_ann <- HeatmapAnnotation(df=hm.identity, 
                             which ="col", col = colours, annotation_name_gp = gpar(fontface="bold"))


fat_heatmap <- Heatmap(df.matrix, 
                       name = "Fat %", #title of legend
                       cluster_columns = T, cluster_rows = F,
                       #row_title = "Phenotypes", 
                       row_names_gp = gpar(fontsize = 10),show_column_names = F, 
                       top_annotation = col_ann, na_col = "black", 
                       #column_split = 4, 
                       #right_annotation = ha,
                       column_gap = unit(2, "mm"), row_gap = unit(2, "mm"), border = F) # Text size for row names

fat_heatmap
save_pdf(
  fat_heatmap,
  "./Plots/03_fat_heatmap.pdf",
  width = 25,
  height = 5,
  units = "cm")





########heatmap z scores##########

########tissue weight##########
########grams##########
dfPlot <- data_wide[,c(1,2,4, 40:43)]
dfPlot <- dfPlot[!is.na(dfPlot$Heart_g_24),]
hm.identity <- as.data.frame(dfPlot[,c(1,2,3)])
rownames(hm.identity) <- hm.identity[,1]
hm.identity <- hm.identity[,-1]
hm.identity$Strain <- factor(hm.identity$Strain,levels = c("C57BL/6J","129S1/SvImJ", "CAST/EiJ","PWK/PhJ", 
                                                           "B6CASTF1","129SPWKF1","B6CAST-129SPWK-F2"))
phenotype_names <- colnames(dfPlot)
#pheno_names <- read.csv(file = './Aim2_Step2_manuscript/Data/phenotype_names.csv')
#pheno_names <- pheno_names[-c(1:3),]

dfPlot$Strain_Sex <- paste0(dfPlot$Strain,"_",dfPlot$Sex)
dfPlot <- reshape2::melt(dfPlot, id.vars = c("ID", "Strain", "Sex", "Strain_Sex"))
dfPlot <- dfPlot %>% group_by(variable) %>% dplyr::mutate(zscore = (value - mean(value, na.rm = T))/sd(value, na.rm = T) ) %>%
  group_by(variable, ID)  %>% dplyr::summarise(mean_zscore = mean(zscore, na.rm = T))
dfPlot <- dfPlot %>% spread(ID, mean_zscore)
colnames(dfPlot)[colnames(dfPlot)=="variable"] <- "Phenotype"


#dfPlot$Pheno_new <- pheno_names$Phenotype_new[match(dfPlot$Phenotype_old,pheno_names$Phenotype_old)]
df <- dfPlot[,-c(1)]
rownames(df) <- dfPlot$Phenotype
#df <- t(df)
df.matrix <- as.matrix(df)


#heatmap(df.matrix, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)
col_ann <- HeatmapAnnotation(df=hm.identity, 
                             which ="col", col = colours, annotation_name_gp = gpar(fontface="bold"))


tissue_heatmap <- Heatmap(df.matrix, 
                       name = "Fat %", #title of legend
                       cluster_columns = T, cluster_rows = T,
                       row_title = "Phenotypes", 
                       row_names_gp = gpar(fontsize = 10),show_column_names = F, 
                       top_annotation = col_ann, na_col = "black", 
                       #column_split = 5, 
                       #right_annotation = ha,
                       column_gap = unit(2, "mm"), row_gap = unit(2, "mm"), border = F) # Text size for row names

tissue_heatmap
save_pdf(
  tissue_heatmap,
  "./Plots/03_tissue_g_zscore_heatmap.pdf",
  width = 25,
  height = 6,
  units = "cm")

########relative to BW##########
data_wide <- read.csv(file = './Aim2_Step2_manuscript/Data/allData_wide_clean_noOutlier.csv')
dfPlot <- data_wide[,c(1,2,4, 44:47)]
dfPlot <- dfPlot[!is.na(dfPlot$Heart_g_RelBWSac_24),]
dfPlot <- dfPlot[!is.na(dfPlot$Liver_g_RelBWSac_24),]
dfPlot <- dfPlot[!is.na(dfPlot$Kidney_g_RelBWSac_24),]
dfPlot <- dfPlot[!is.na(dfPlot$Spleen_g_RelBWSac_24),]
hm.identity <- as.data.frame(dfPlot[,c(1,2,3)])
rownames(hm.identity) <- hm.identity[,1]
hm.identity <- hm.identity[,-1]
hm.identity$Strain <- factor(hm.identity$Strain,levels = c("C57BL/6J","129S1/SvImJ", "CAST/EiJ","PWK/PhJ", 
                                                           "B6CASTF1","129SPWKF1","B6CAST-129SPWK-F2"))
phenotype_names <- colnames(dfPlot)
#pheno_names <- read.csv(file = './Aim2_Step2_manuscript/Data/phenotype_names.csv')
#pheno_names <- pheno_names[-c(1:3),]

dfPlot$Strain_Sex <- paste0(dfPlot$Strain,"_",dfPlot$Sex)
dfPlot <- reshape2::melt(dfPlot, id.vars = c("ID", "Strain", "Sex", "Strain_Sex"))
dfPlot <- dfPlot %>% group_by(variable) %>% dplyr::mutate(zscore = (value - mean(value, na.rm = T))/sd(value, na.rm = T) ) %>%
  group_by(variable, ID)  %>% dplyr::summarise(mean_zscore = mean(zscore, na.rm = T))
dfPlot <- dfPlot %>% spread(ID, mean_zscore)
colnames(dfPlot)[colnames(dfPlot)=="variable"] <- "Phenotype"


#dfPlot$Pheno_new <- pheno_names$Phenotype_new[match(dfPlot$Phenotype_old,pheno_names$Phenotype_old)]
df <- dfPlot[,-c(1)]
rownames(df) <- dfPlot$Phenotype
#df <- t(df)
df.matrix <- as.matrix(df)


#heatmap(df.matrix, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)
col_ann <- HeatmapAnnotation(df=hm.identity, 
                             which ="col", col = colours, annotation_name_gp = gpar(fontface="bold"))


tissue_rel_heatmap <- Heatmap(df.matrix, 
                          name = "zscore", #title of legend
                          cluster_columns = T, cluster_rows = T,
                          row_title = "Phenotypes", 
                          row_names_gp = gpar(fontsize = 10),show_column_names = F, 
                          top_annotation = col_ann, na_col = "black", 
                          #column_split = 3, 
                          #right_annotation = ha,
                          column_gap = unit(2, "mm"), row_gap = unit(2, "mm"), border = F) # Text size for row names

tissue_rel_heatmap
save_pdf(
  tissue_rel_heatmap,
  "./Plots/03_tissue_rel_zscore_heatmap_noOutlier_noNA.pdf",
  width = 25,
  height = 6,
  units = "cm")


#######liver heatmap############
#allData_wide <- allData_wide[,-c(4:13)]
#temporary script to include the fibrosis
#data_wide <- read.csv(file = './Aim2_Step2_manuscript/Data/allData_wide_clean_noOutlier.csv')


dfPlot <- data_wide %>% dplyr::select(ID, Strain, Sex, 
                                      BW_perGain_22,
                                      Fat_perc_22, 
                                      Liver_g_RelBWSac_24, 
                                      Liver_g_24,
                                      ALAT_24, 
                                      #ASAT_24,
                                      TIMP1_24, 
                                      #GDF15_24, 
                                      #Cholesterol_24,
                                      #large_vacoules_area_percentage_24, 
                                      #small_vacuoles_area_percentage_24, 
                                      sum_all_vacuoles_percentage_24,
                                      Fibrosis_perc_24
                                      #Perc_Coll_IV_Positive_Tissue_24
                                      )
#dfPlot <- dfPlot[grepl("f",dfPlot$Sex), ]
dfPlot <- dfPlot[!is.na(dfPlot$Liver_g_RelBWSac_24),]
dfPlot <- dfPlot[!is.na(dfPlot$Liver_g_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$BW_perGain_22),]
dfPlot <- dfPlot[!is.na(dfPlot$Fat_perc_22),]
#dfPlot <- dfPlot[!is.na(dfPlot$ASAT),]
dfPlot <- dfPlot[!is.na(dfPlot$ALAT_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$ASAT_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$GDF15_24),]
dfPlot <- dfPlot[!is.na(dfPlot$TIMP1_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$Cholesterol_24),]
dfPlot <- dfPlot[!is.na(dfPlot$sum_all_vacuoles_percentage_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$large_vacoules_area_percentage_24),]
dfPlot <- dfPlot[!is.na(dfPlot$Fibrosis_perc_24),]
hm.identity <- as.data.frame(dfPlot[,c(1,2,3)])
rownames(hm.identity) <- hm.identity[,1]
hm.identity <- hm.identity[,-1]
hm.identity$Strain <- factor(hm.identity$Strain,levels = c("C57BL/6J","129S1/SvImJ", "CAST/EiJ","PWK/PhJ", 
                                                           "B6CASTF1","129SPWKF1","B6CAST-129SPWK-F2"))
phenotype_names <- colnames(dfPlot)



dfPlot$Strain_Sex <- paste0(dfPlot$Strain,"_",dfPlot$Sex)
dfPlot <- reshape2::melt(dfPlot, id.vars = c("ID", "Strain", "Sex", "Strain_Sex"))
dfPlot <- dfPlot %>% group_by(variable) %>% dplyr::mutate(zscore = (value - mean(value, na.rm = T))/sd(value, na.rm = T) ) %>%
  group_by(variable, ID)  %>% dplyr::summarise(mean_zscore = mean(zscore, na.rm = T))
dfPlot <- dfPlot %>% spread(ID, mean_zscore)
colnames(dfPlot)[colnames(dfPlot)=="variable"] <- "Phenotype"
#dfPlot$Pheno_new <- pheno_names$Phenotype_new[match(dfPlot$Phenotype_old,pheno_names$Phenotype_old)]
df <- dfPlot[,-c(1)]
rownames(df) <- dfPlot$Phenotype
#df <- t(df)
df.matrix <- as.matrix(df)


#heatmap(df.matrix, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)

col_ann <- HeatmapAnnotation(df=hm.identity, 
                             which ="col", col = colours, annotation_name_gp = gpar(fontface="bold"))
#col_order <- order(df.matrix[6, ])

liver_heatmap <- Heatmap(df.matrix, 
                       name = "z-score", #title of legend
                       cluster_columns = T, cluster_rows = T,
                       #row_title = "Phenotypes", 
                       row_names_gp = gpar(fontsize = 10),show_column_names = F, 
                       top_annotation = col_ann, na_col = "black",
                       #column_order = col_order,
                       column_km = 4,
                       #column_split = 2, 
                       #right_annotation = ha,
                       column_gap = unit(2, "mm"), row_gap = unit(2, "mm"), border = F) # Text size for row names

liver_heatmap
save_pdf(
  liver_heatmap,
  "./Plots/03_liver_heatmap.pdf",
  width = 27,
  height = 5.2,
  units = "cm")

ht2 = draw(liver_heatmap)
col <- column_order(ht2)


column_names_list <- lapply(col, function(pos) colnames(df.matrix)[pos])


# Flatten the list and create a corresponding identifier for each element
df <- data.frame(
  Column_Names = unlist(column_names_list),
  List_Name = rep(names(column_names_list), lengths(column_names_list))
)

# Print the result
print(df)
colnames(df) <- c("ID", "liver_cluster")
unique(df$ID)

data_melt <- merge(data_melt, df, by = "ID", all.x = TRUE)
data_wide <- merge(data_wide, df, by = "ID", all.x = TRUE)

counts <- ddply(data_wide, .(data_wide$liver_cluster
                          , data_wide$Sex
), nrow)

names(counts) <- c("Cluster"
                   , "Sex"
                   , "Freq")

counts_clust <- ddply(data_wide, .(data_wide$liver_cluster), nrow)
names(counts_clust) <- c("Cluster", "Freq")
counts_clust$Sex <- "both"
counts_clust <- counts_clust[, c(1,3,2)]
counts <- rbind(counts, counts_clust)
counts <- counts[!is.na(counts$Cluster),]
#library(viridis)
b1 <- ggplot(counts, aes(fill=Cluster, 
                         x=Sex, 
                         y=Freq)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(name = "Cluster",values = c("#0696c7","#FFD700","#ff8800","#6e0280"))+
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        axis.ticks.x = element_blank(),
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
b1
ggsave(filename ="03-cluster_freq.pdf",plot=b1,path ="./Plots",width =2.5, height = 2, units = c("in"))

b2 <- ggplot(counts[!grepl("both", counts$Sex),], aes(fill=Sex, x=Cluster, y=Freq)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(name = "Sex",values = sexColors)+
  theme_cowplot()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8),
        #axis.line = element_blank(),
        axis.ticks.x = element_blank(),
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
b2
ggsave(filename ="03-cluster_freq_m_f.pdf",plot=b2,path ="./Plots",width =2.5, height = 2, units = c("in"))

# Remove rows where liver_cluster column has NA
df <- data_melt[!is.na(data_melt$liver_cluster), ]

g1 <- ggplot(df[grepl("Fat_perc_22",df$unique_variable) 
                       |grepl("Liver_g_RelBWSac_24",df$unique_variable)
                       |grepl("sum_all_vacuoles_percentage_24",df$unique_variable)
                       |grepl("TIMP1_24",df$unique_variable)
                       |grepl("ALAT_24",df$unique_variable)
                |grepl("Fibrosis_perc_24",df$unique_variable)
                |grepl("Cholesterol_24",df$unique_variable)
                & !grepl("HDL-Cholesterol_24",df$unique_variable)
                & !grepl("LDL-Cholesterol_24",df$unique_variable)
                |grepl("Liver_g_24",df$unique_variable)
                |grepl("BW_perGain_24",df$unique_variable)
                |grepl("GDF15_24",df$unique_variable)
                |grepl("ACR_22",df$unique_variable)
                |grepl("TGL_24",df$unique_variable),],
             aes(y = Value, 
                 x = liver_cluster)) +
  geom_violin(aes(color = liver_cluster, fill=liver_cluster))+
  geom_boxplot(width=0.1, outlier.shape = NA,position=position_dodge(.9))+
  #geom_point(mapping = aes(group = ID),
  #           size = 0.05, 
  #           #shape = 21,
  #           position = position_dodge(0.1)) +
  scale_color_manual(name = "Cluster",values = c("#0696c7","#FFD700","#ff8800","#6e0280"))+
  scale_fill_manual(name = "Cluster",values = c("#0696c7","#FFD700","#ff8800","#6e0280"))+
  facet_wrap(~ unique_variable , scales = "free")+
  #scale_alpha_manual(values = c(0.2,1))+
  ggtitle("Liver Clusters") +
  ylab("")+
  theme_bw()+
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
        strip.text = element_text(size=8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        #legend.position = "bottom",
        panel.background = element_blank(),
        #strip.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
g1
ggsave(filename ="03-cluster_phenotypes.pdf",plot=g1,path ="./Plots",
       width =6.8, height = 4, units = c("in"))
g2 <- ggplot(df,
             aes(y = Value, 
                 x = liver_cluster)) +
  geom_violin(aes(color = liver_cluster, fill=liver_cluster))+
  geom_boxplot(width=0.1, outlier.shape = NA,position=position_dodge(.9))+
  #geom_point(mapping = aes(group = ID),
  #           size = 0.05, 
  #           #shape = 21,
  #           position = position_dodge(0.1)) +
  scale_color_manual(name = "Cluster",values = c("#0696c7","#FFD700","#ff8800","#6e0280"))+
  scale_fill_manual(name = "Cluster",values = c("#0696c7","#FFD700","#ff8800","#6e0280"))+
  facet_wrap(~ unique_variable , scales = "free")+
  #scale_alpha_manual(values = c(0.2,1))+
  ggtitle("Liver Clusters") +
  ylab("")+
  theme_bw()+
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
g2

#######liver heatmap --- F2 only############
#allData_wide <- allData_wide[,-c(4:13)]
dfPlot <- data_wide %>% dplyr::select(ID, Strain, Generation, Sex, 
                                      #BW_perGain_22,
                                      Fat_perc_22, 
                                      Liver_g_RelBWSac_24, 
                                      ALAT_24, 
                                      TIMP1_24, 
                                      #GDF15_24, 
                                      #Cholesterol_24,
                                      #large_vacoules_area_percentage_24, 
                                      #small_vacuoles_area_percentage_24, 
                                      sum_all_vacuoles_percentage_24,
                                      Fibrosis_perc_24
                                      #, 
                                      #Perc_Coll_IV_Positive_Tissue_24
)
dfPlot <- dfPlot[grepl("F2", dfPlot$Generation),]
dfPlot <- dfPlot[,-3]
#dfPlot <- dfPlot[grepl("f",dfPlot$Sex), ]
dfPlot <- dfPlot[!is.na(dfPlot$Liver_g_RelBWSac_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$BW_perGain_22),]
dfPlot <- dfPlot[!is.na(dfPlot$Fat_perc_22),]
#dfPlot <- dfPlot[!is.na(dfPlot$ASAT),]
dfPlot <- dfPlot[!is.na(dfPlot$ALAT_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$GDF15_24),]
dfPlot <- dfPlot[!is.na(dfPlot$TIMP1_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$Cholesterol_24),]
dfPlot <- dfPlot[!is.na(dfPlot$sum_all_vacuoles_percentage_24),]
dfPlot <- dfPlot[!is.na(dfPlot$Fibrosis_perc_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$large_vacoules_area_percentage_24),]
#dfPlot <- dfPlot[!is.na(dfPlot$Perc_Coll_IV_Positive_Tissue_24),]
hm.identity <- as.data.frame(dfPlot[,c(1,2,3)])
rownames(hm.identity) <- hm.identity[,1]
hm.identity <- hm.identity[,-1]
hm.identity$Strain <- factor(hm.identity$Strain,levels = c("B6CAST-129SPWK-F2"))
phenotype_names <- colnames(dfPlot)



dfPlot$Strain_Sex <- paste0(dfPlot$Strain,"_",dfPlot$Sex)
dfPlot <- reshape2::melt(dfPlot, id.vars = c("ID", "Strain", "Sex", "Strain_Sex"))
dfPlot <- dfPlot %>% group_by(variable) %>% dplyr::mutate(zscore = (value - mean(value, na.rm = T))/sd(value, na.rm = T) ) %>%
  group_by(variable, ID)  %>% dplyr::summarise(mean_zscore = mean(zscore, na.rm = T))
dfPlot <- dfPlot %>% spread(ID, mean_zscore)
colnames(dfPlot)[colnames(dfPlot)=="variable"] <- "Phenotype"
#dfPlot$Pheno_new <- pheno_names$Phenotype_new[match(dfPlot$Phenotype_old,pheno_names$Phenotype_old)]
df <- dfPlot[,-c(1)]
rownames(df) <- dfPlot$Phenotype
#df <- t(df)
df.matrix <- as.matrix(df)


#heatmap(df.matrix, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)

col_ann <- HeatmapAnnotation(df=hm.identity, 
                             which ="col", col = colours, annotation_name_gp = gpar(fontface="bold"))

liver_heatmap <- Heatmap(df.matrix, 
                         name = "z-score", #title of legend
                         cluster_columns = T, cluster_rows = F,
                         #row_title = "Phenotypes", 
                         row_names_gp = gpar(fontsize = 10),show_column_names = F, 
                         top_annotation = col_ann, na_col = "black",
                         column_km = 4,
                         #cluster_column_slices = FALSE, 
                         #right_annotation = ha,
                         column_gap = unit(2, "mm"), row_gap = unit(2, "mm"), border = F) # Text size for row names

liver_heatmap

save_pdf(
  liver_heatmap,
  "./Plots/03-liver_heatmap_F2only.pdf",
  width = 27,
  height = 5.2,
  units = "cm")
