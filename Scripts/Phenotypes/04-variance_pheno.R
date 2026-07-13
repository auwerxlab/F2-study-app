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
library(stringr)
#library(rawswap)
library(ggpubr)
library(ggh4x)
library(ggbeeswarm)
library(scales)
library(EnvStats)
library(data.table)
library(naniar)
library(tigerstats)

source("./Scripts/Phenotypes/00-common.R")
sexColors2 <- c("both"="#A9A9A9", "m" = "#00A2FF", "f" = "#FF644E")

Phenotypes<- read.csv(file = './Aim2_Step2_manuscript/Data/allData_melt_clean_noOutlier.csv')

############################## data formatting  ####################################
unique(Phenotypes$name)
Phenotypes$name<-gsub("_g_RelBWSac"," weight",Phenotypes$Variable)
Phenotypes$name<-gsub("Lean_g","Lean mass",Phenotypes$name)
Phenotypes$name<-gsub("Lean_perc","Lean mass",Phenotypes$name)
Phenotypes$name<-gsub("Fat_g","Fat mass",Phenotypes$name)
Phenotypes$name<-gsub("Fat_perc","Fat mass",Phenotypes$name)
Phenotypes$name[Phenotypes$name == "BW"] <- "Body Weight" 
Phenotypes$name[Phenotypes$name == "BW_Sac"] <- "Body Weight" 
Phenotypes$name[Phenotypes$name == "BW_perGain"] <- "Body Weight" 
#Phenotypes$name<-gsub("_perc"," weight",Phenotypes$name)
Phenotypes$name<-gsub("_DAY","",Phenotypes$name)
Phenotypes$name<-gsub("_NIGHT","",Phenotypes$name)
Phenotypes <- Phenotypes[!grepl("delta", Phenotypes$name),]
unique(Phenotypes$name)
Phenotypes$name<-gsub("_g"," weight",Phenotypes$name)

Phenotypes$name[Phenotypes$name == "Ammonia"] <- "Ammonia (plasma)"
Phenotypes$name[Phenotypes$name == "Glycemia"] <- "Glycemia (plasma)"
Phenotypes$name[Phenotypes$name == "Urine_vol"] <- "Volume (urine)"
Phenotypes$name[Phenotypes$name == "Urine_albumin"] <- "Albumin (urine)"
Phenotypes$name[Phenotypes$name == "Urine_creatinin"] <- "Creatinine (urine)"
Phenotypes$name[Phenotypes$name == "Urine_urea"] <- "Urea (urine)"
Phenotypes$name[Phenotypes$name == "ASAT"] <- "ASAT (plasma)"
Phenotypes$name[Phenotypes$name == "ALAT"] <- "ALAT (plasma)"
Phenotypes$name[Phenotypes$name == "Cholesterol"] <- "Chol (plasma)"
Phenotypes$name[Phenotypes$name == "HDL-Cholesterol"] <- "HDL-Chol (plasma)"
Phenotypes$name[Phenotypes$name == "LDL-Cholesterol"] <- "LDL-Chol (plasma)"
Phenotypes$name[Phenotypes$name == "TGL"] <- "TGL (plasma)"
Phenotypes$name[Phenotypes$name == "Crea"] <- "Creatinine (plasma)"
Phenotypes$name[Phenotypes$name == "Urea"] <- "Urea (plasma)"
Phenotypes$name[Phenotypes$name == "Free Fatty Acids"] <- "FFA (plasma)"
Phenotypes$name[Phenotypes$name == "GDF15"] <- "GDF-15 (plasma)"
Phenotypes$name[Phenotypes$name == "TIMP1"] <- "TIMP-1 (plasma)"
Phenotypes$name[Phenotypes$name == "XAMB"] <- "Activity"
Phenotypes$name[Phenotypes$name == "RER"] <- "Respiratory Exchange Ratio"
Phenotypes$name[Phenotypes$name == "EE"] <- "Energy Expenditure"
Phenotypes$name[Phenotypes$name == "Food_intake_CLAMS"] <- "Food Intake"
Phenotypes$name[Phenotypes$name == "Fibrosis_perc"] <- "Histology (liver)"
Phenotypes$name[Phenotypes$name == "large_vacoules_area_percentage"] <- "Histology (liver)"
Phenotypes$name[Phenotypes$name == "small_vacuoles_area_percentage"] <- "Histology (liver)"
Phenotypes$name[Phenotypes$name == "sum_all_vacuoles_percentage"] <- "Histology (liver)"
Phenotypes$name[Phenotypes$name == "Perc_Coll_IV_Positive_Tissue"] <- "Histology (liver)"
unique(Phenotypes$name)


unique(Phenotypes$name)
Phenotypes <- Phenotypes[!grepl("Col IV raw values", Phenotypes$name),]
Phenotypes <- Phenotypes[!grepl("KIM1", Phenotypes$name),]
Phenotypes <- Phenotypes[!grepl("VO2", Phenotypes$name),]
Phenotypes <- Phenotypes[!grepl("VCO2", Phenotypes$name),]
Phenotypes <- Phenotypes[!grepl("BW_loss_perc_CLAMS", Phenotypes$name),]
#Phenotypes$name<-gsub( "% Coll IV positive tissue normalized", "% Col IV pos",Phenotypes$name,fixed=T)

Phenotypes$name_week<-paste0(Phenotypes$name," (",Phenotypes$Week,"w)")

############# Which pheno are available at what timepoint? #####################

Counts_condition<-as.data.frame(table(Phenotypes[,c("Week","name")]))
Counts_condition<-Counts_condition[Counts_condition$Freq>0,]

df = data.frame(Week = "24", name="RNAseq (liver)", Freq = "569")
Counts_condition = rbind(Counts_condition,df)

df2 = data.frame(Week = "24", name="Proteomics (liver)", Freq = "569")
Counts_condition = rbind(Counts_condition,df2)

Counts_condition$name <- factor(Counts_condition$name,levels = c("Body Weight","Fat mass",
                                                                 "Lean mass","Liver weight",
                                                                 "Spleen weight","Kidney weight",
                                                                 "Heart weight","Activity","Food Intake","Energy Expenditure",
                                                                 "Respiratory Exchange Ratio",
                                                                 "Glycemia (plasma)", "Ammonia (plasma)",
                                                                 "ALAT (plasma)","ASAT (plasma)","Chol (plasma)",
                                                                 "HDL-Chol (plasma)","LDL-Chol (plasma)",
                                                                 "TGL (plasma)","FFA (plasma)","Creatinine (plasma)", 
                                                                 "Urea (plasma)", "TIMP-1 (plasma)", "GDF-15 (plasma)","ACR",
                                                                 "Albumin (urine)","Creatinine (urine)","Urea (urine)","Volume (urine)", "Histology (liver)",                                                                                
                                                                 "RNAseq (liver)","Proteomics (liver)"))

#Counts_condition$Week <- as.numeric(as.character(Counts_condition$Week))       

plot_pheno_availability<-ggplot(Counts_condition,aes(x = Week,y=name))+
  geom_point(color="darkgrey")+
  scale_y_discrete(limits=rev)+
  theme_minimal()+
  theme(
    #axis.line = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.line = element_line(colour = 'black', size = 0.3),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size=8),
    axis.title.y = element_blank(),
    legend.text=element_text(size=8),
    legend.title = element_text(size=8),
    legend.position = "bottom",
    panel.background = element_blank(),
    #strip.text = element_text(size = 6),
    legend.key.size = unit(0.5, "cm"))
plot_pheno_availability
ggsave(plot_pheno_availability,file="./Plots/04-Phenotype_availability_timepoint.pdf",
       width = 8.5, height = 9.5, units = c("cm"))
############################# plot steatosis variance ##################
steatosis <- ggplot(Phenotypes[grepl("sum_all_vacuoles_percentage_24",Phenotypes$unique_variable),], 
                    aes(x = reorder(ID, Value), y = Value, color = Sex, fill = Sex)) +
  geom_bar(stat="identity")+
  scale_color_manual(name = "Sex",values = sexColors)+
  scale_fill_manual(name = "Sex",values = sexColors)+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 30))+
  #scale_alpha_manual(values = c(0.3,1))+
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title=element_text(size=8),
        legend.title = element_text(size=7),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Steatosis (%)")+
  xlab("Strain")
#theme(axis.text.x = element_blank(),
#      axis.text=element_text(size=6),
#      axis.title=element_text(size=6.5),
#      legend.text=element_text(size=4),
#      legend.title = element_text(size=5),
#      axis.line = element_line(colour = 'black', size = 0.4),
#      axis.ticks.y = element_line(size=0.4),
#      axis.ticks.length.y = unit(.07, "cm"),
#     axis.ticks.x = element_blank(),
#      legend.key.size = unit(0.3, "cm"))+
steatosis
ggsave(useDingbats = FALSE, steatosis,file="./Plots/04-steatosis_variation_histo.pdf",
       width = 7, height = 2.5, units = c("in"))

fibrosis <- ggplot(Phenotypes[grepl("Fibrosis_perc_24",Phenotypes$unique_variable),], 
                   aes(x = reorder(ID, Value), y = Value, color = Sex, fill = Sex)) +
  geom_bar(stat="identity")+
  scale_color_manual(name = "Sex",values = sexColors)+
  scale_fill_manual(name = "Sex",values = sexColors)+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 30))+
  #scale_alpha_manual(values = c(0.3,1))+
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title=element_text(size=8),
        legend.title = element_text(size=7),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Steatosis (%)")+
  xlab("Strain")
#theme(axis.text.x = element_blank(),
#      axis.text=element_text(size=6),
#      axis.title=element_text(size=6.5),
#      legend.text=element_text(size=4),
#      legend.title = element_text(size=5),
#      axis.line = element_line(colour = 'black', size = 0.4),
#      axis.ticks.y = element_line(size=0.4),
#      axis.ticks.length.y = unit(.07, "cm"),
#     axis.ticks.x = element_blank(),
#      legend.key.size = unit(0.3, "cm"))+
fibrosis
#ggsave(useDingbats = FALSE, steatosis,file="./Plots/04-steatosis_variation_histo.pdf",
#       width = 7, height = 2.5, units = c("in"))

steatosis2 <- ggplot(Phenotypes[grepl("sum_all_vacuoles_percentage_24",Phenotypes$unique_variable),], 
                     aes(x = reorder(ID, Value), y = Value, color = Strain, fill = Strain)) +
  geom_bar(stat="identity")+
  scale_color_manual(name = "Strain",values = strainColors)+
  scale_fill_manual(name = "Strain",values = strainColors)+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 30))+
  #scale_alpha_manual(values = c(0.3,1))+
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title=element_text(size=8),
        legend.title = element_text(size=7),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Steatosis (%)")+
  xlab("Strain")
#theme(axis.text.x = element_blank(),
#      axis.text=element_text(size=6),
#      axis.title=element_text(size=6.5),
#      legend.text=element_text(size=4),
#      legend.title = element_text(size=5),
#      axis.line = element_line(colour = 'black', size = 0.4),
#      axis.ticks.y = element_line(size=0.4),
#      axis.ticks.length.y = unit(.07, "cm"),
#     axis.ticks.x = element_blank(),
#      legend.key.size = unit(0.3, "cm"))+
steatosis2
ggsave(useDingbats = FALSE, steatosis2,file="./Plots/04-steatosis_variation_histo.pdf",
       width = 7, height = 2.5, units = c("in"))



############################## calculate phenotype stats  ####################################
data_wide <- read.csv(file = './Aim2_Step2_manuscript/Data/allData_wide_clean_noOutlier.csv')
data_wide_m <- data_wide[grepl("m", data_wide$Sex),]
data_wide_f <- data_wide[grepl("f", data_wide$Sex),]


#####both sexes
df <- data_wide[,c(5:120)]
#favstats(~BW_6,data=data_wide)
stats <- lapply(df, favstats)
stats <- rbindlist(stats,use.names=TRUE, idcol = TRUE)

variance <- lapply(df, var, na.rm = TRUE)
variance <- unlist(variance)
variance <- as.data.frame(variance)

variance$Phenotypes <- rownames(variance)

stats <- merge(stats, variance, by.x=".id", by.y="Phenotypes")
stats$cv <- stats$sd/stats$mean
stats$Sex <- "both"


#####males
df_m <- data_wide_m[,c(5:115)]
#favstats(~BW_6,data=data_wide)
stats_m <- lapply(df_m, favstats)
stats_m <- rbindlist(stats_m,use.names=TRUE, idcol = TRUE)

variance_m <- lapply(df_m, var, na.rm = TRUE)
variance_m <- unlist(variance_m)
variance_m <- as.data.frame(variance_m)
colnames(variance_m)<- c("variance")

variance_m$Phenotypes <- rownames(variance_m)

stats_m <- merge(stats_m, variance_m, by.x=".id", by.y="Phenotypes")
stats_m$cv <- stats_m$sd/stats_m$mean
stats_m$Sex <- "m"


#####females
df_f <- data_wide_f[,c(5:115)]
#favstats(~BW_6,data=data_wide)
stats_f <- lapply(df_f, favstats)
stats_f <- rbindlist(stats_f,use.names=TRUE, idcol = TRUE)

variance_f <- lapply(df_f, var, na.rm = TRUE)
variance_f <- unlist(variance_f)
variance_f <- as.data.frame(variance_f)
colnames(variance_f)<- c("variance")

variance_f$Phenotypes <- rownames(variance_f)

stats_f <- merge(stats_f, variance_f, by.x=".id", by.y="Phenotypes")
stats_f$cv <- stats_f$sd/stats_f$mean
stats_f$Sex <- "f"

stats <- rbind(stats,stats_m, stats_f)

bw <- stats[grepl("BW", stats$.id) 
            & !grepl("perGain", stats$.id)
            & !grepl("Sac", stats$.id)
            & !grepl("loss", stats$.id)
            ,]

bw$.id <- factor(bw$.id,levels = c("BW_6","BW_7","BW_8",
                                   "BW_10","BW_12","BW_14",
                                   "BW_16","BW_18","BW_20",
                                   "BW_22","BW_24"))


v1 <- ggplot(data=bw[!grepl("both", bw$Sex),], aes(x=.id, y=variance, group=Sex, fill=Sex)) +
  geom_bar(position="dodge", stat="identity")+
  xlab("")+
  ylab("variance")+
  scale_fill_manual(values = sexColors2)+
  theme_minimal()+
  theme(#axis.text.x = element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.text=element_text(size=8),
    #legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_line(size=0.3),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    #legend.position = "none"
    legend.key.size = unit(0.05, "cm")
  )
v1
ggsave(filename ="04-bw_variance.pdf",plot=v1,
       path ="./Plots",width =3, height = 1.7, units = c("in"))


fat <- stats[grepl("Fat", stats$.id) 
             & !grepl("perc", stats$.id)
             & !grepl("Acids", stats$.id)
             ,]

fat$.id <- factor(fat$.id,levels = c("Fat_g_6","Fat_g_14","Fat_g_22"))


v2 <- ggplot(data=fat[!grepl("both", fat$Sex),], aes(x=.id, y=variance, group=Sex, fill=Sex)) +
  geom_bar(position="dodge", stat="identity")+
  xlab("")+
  ylab("Variance")+
  scale_fill_manual(values = sexColors2)+
  theme_minimal()+
  theme(#axis.text.x = element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.text=element_text(size=8),
    #legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_line(size=0.3),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    #legend.position = "none"
    legend.key.size = unit(0.05, "cm")
  )
v2
ggsave(filename ="04-fat_g_variance.pdf",plot=v2,
       path ="./Plots",width =2, height = 1.7, units = c("in"))


lean <- stats[grepl("Lean", stats$.id) 
              & !grepl("perc", stats$.id),]

lean$.id <- factor(lean$.id,levels = c("Lean_g_6","Lean_g_14","Lean_g_22"))


v3 <- ggplot(data=lean[!grepl("both", lean$Sex),], aes(x=.id, y=variance, group=Sex, fill=Sex)) +
  geom_bar(position="dodge", stat="identity")+
  xlab("")+
  ylab("Variance")+
  scale_fill_manual(values = sexColors2)+
  theme_minimal()+
  theme(#axis.text.x = element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.text=element_text(size=8),
    #legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_line(size=0.3),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    #legend.position = "none"
    legend.key.size = unit(0.05, "cm")
  )
v3
ggsave(filename ="04-lean_g_variance.pdf",plot=v3,
       path ="./Plots",width =2, height = 1.7, units = c("in"))


fatlean <- stats[grepl("_g", stats$.id) 
                 & !grepl("24", stats$.id)
                 #& !grepl("Acids", stats$.id)
                 ,]

fatlean$.id <- factor(fatlean$.id,levels = c("Fat_g_6","Fat_g_14","Fat_g_22",
                                             "Lean_g_6","Lean_g_14","Lean_g_22"))


v4 <- ggplot(data=fatlean[!grepl("both", fatlean$Sex),], aes(x=.id, y=variance, group=Sex, fill=Sex)) +
  geom_bar(position="dodge", stat="identity")+
  xlab("")+
  ylab("Variance")+
  scale_fill_manual(values = sexColors2)+
  theme_minimal()+
  theme(#axis.text.x = element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.text=element_text(size=8),
    #legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_line(size=0.3),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    #legend.position = "none"
    legend.key.size = unit(0.05, "cm")
  )
v4
ggsave(filename ="04-fatlean_g_variance.pdf",plot=v4,
       path ="./Plots",width =2.5, height = 1.7, units = c("in"))

vn <- plot_grid(v1,v4,nrow = 1)
ggsave("GB03_bw_fatlean_g_variance.pdf",plot=vn, path="./Plots/Giorgia",width = 7, height = 1.7, units = c("in"))

tissue <- stats[grepl("Liver", stats$.id) 
                | grepl("Heart", stats$.id)
                | grepl("Kidney", stats$.id)
                | grepl("Spleen", stats$.id)
                #& !grepl("RelBW", stats$.id)
                ,]

tissue <- tissue[!grepl("Sac",tissue$.id)
                 ,]

t <- ggplot(data=tissue, aes(x=.id, y=cv, group=Sex, fill=Sex)) +
  geom_bar(position="dodge", stat="identity")+
  xlab("")+
  ylab("CV")+
  scale_fill_manual(values = sexColors2)+
  theme_minimal()+
  theme(#axis.text.x = element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.text=element_text(size=8),
    #legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_line(size=0.3),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    #legend.position = "none"
    legend.key.size = unit(0.05, "cm")
  )

t


############coefficient of variation################
pheno_sel <- stats[!grepl("Urine",stats$.id)
                   & !grepl("ACR", stats$.id)
                   & !grepl("delta", stats$.id)
                   & !grepl("ColIV", stats$.id)
                   & !grepl("KIM1", stats$.id)
                   ,]
pheno_sel <- pheno_sel[!is.na(pheno_sel$cv),]

b <- ggplot(data=pheno_sel[grepl("both", pheno_sel$Sex),], 
            aes(x=cv, y=reorder(.id,cv))) +
  geom_point(size=2, shape=18, aes(color=Sex))+
  #ylim(0,1.5)+
  xlab("Coefficient of Variation")+
  ylab("")+
  #facet_wrap(~Sex, nrow = 3)+
  scale_color_manual(values = sexColors2)+
  theme_minimal()+
  theme(#axis.text.x = element_blank(),
    axis.text=element_text(size=7),
    axis.title=element_text(size=7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.text=element_text(size=8),
    #legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_line(size=0.3),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
    #legend.key.size = unit(0.1, "cm")
  )
b
ggsave(filename ="GB03_cv_phenotypes.pdf",plot=b,
       path ="./Plots/Giorgia",width = 2.8, height = 10, units = c("in"))

m <- ggplot(data=pheno_sel[grepl("m", pheno_sel$Sex),], 
            aes(x=cv, y=reorder(.id,cv))) +
  geom_point(size=2, shape=18, aes(color=Sex))+
  #ylim(0,1.5)+
  xlab("Coefficient of Variation")+
  ylab("")+
  #facet_wrap(~Sex, nrow = 3)+
  scale_color_manual(values = sexColors2)+
  theme_minimal()+
  theme(#axis.text.x = element_blank(),
    axis.text=element_text(size=7),
    axis.title=element_text(size=7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.text=element_text(size=8),
    #legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_line(size=0.3),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
    #legend.key.size = unit(0.1, "cm")
  )
m
ggsave(filename ="GB03_cv_phenotypes_m.pdf",plot=m,
       path ="./Plots/Giorgia",width = 2.8, height = 10, units = c("in"))

f <- ggplot(data=pheno_sel[grepl("f", pheno_sel$Sex),], 
            aes(x=cv, y=reorder(.id,cv))) +
  geom_point(size=2, shape=18, aes(color=Sex))+
  #ylim(0,1.5)+
  xlab("Coefficient of Variation")+
  ylab("")+
  #facet_wrap(~Sex, nrow = 3)+
  scale_color_manual(values = sexColors2)+
  theme_minimal()+
  theme(#axis.text.x = element_blank(),
    axis.text=element_text(size=7),
    axis.title=element_text(size=7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.text=element_text(size=8),
    #legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_line(size=0.3),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
    #legend.key.size = unit(0.1, "cm")
  )
f
ggsave(filename ="GB03_cv_phenotypes_f.pdf",plot=f,
       path ="./Plots/Giorgia",width = 2.8, height = 10, units = c("in"))

cv <- plot_grid(b,m,f, nrow = 1)

ggsave(filename ="GB03_cv_phenotypes_merged.pdf",plot=cv,
       path ="./Plots/Giorgia",width = 7.2, height = 9, units = c("in"))

