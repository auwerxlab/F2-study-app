rm(list = ls())
library(tidyverse)
library(cowplot)
library(openxlsx)
theme_set(theme_cowplot())
library("data.table")
library("ggplot2")
library("lubridate")
library(readxl)
library(reshape2)
library(ggplot2)
library(plyr)
library(RColorBrewer)

source("./Scripts/Phenotypes/00-common.R")

bw <- read_excel("Aim2_Step2_manuscript/phenotype_data/A_BW/Aim2Step2_BodyWeight.xlsx")
###### PWK male mice are in reality 129SPWKF1 ######
bw$Strain[which((bw$Strain == "PWK/PhJ") & (bw$Sex == "m"))] <- "129SPWKF1"
###take PWK from step 1###
bw_aim2a <- read_excel("Aim2_Step2_manuscript/phenotype_data/A_BW/Aim2a_BodyWeight_pwk.xlsx")
bw <- rbind(bw, bw_aim2a)
colnames(bw)[colnames(bw)=="Animal ID"] <- "Pyrat_ID"
colnames(bw)[colnames(bw)=="Date of Birth"] <- "DOB"

  
bw_time_course <- bw[bw$Parameter=="BW", ]
bw_time_course <- bw_time_course[!is.na(bw_time_course$Pyrat_ID),]
#bw_time_course <- bw_time_course[!is.na(bw_time_course$Strain),]
#bw_time_course <- bw_time_course[!is.na(bw_time_course$`Week 24`),]
#remove the vet cases and comments
bw_time_course <- bw_time_course[,!grepl("Vet",colnames(bw_time_course)) 
                                 & !grepl("Comments",colnames(bw_time_course))]

bw_time_course <- data.frame(bw_time_course)
bw_time_course[,11:21] <- sapply(bw_time_course[,11:21],as.numeric)

setDT(bw_time_course) 
bw_time_course[Strain=="*C57BL/6J", Strain:="C57BL/6J"]

bw_time_course$Generation <- "F0"
bw_time_course[which(bw_time_course$Strain == "B6CAST-129SPWK-F2"),"Generation"] <- 'F2'
bw_time_course[which(bw_time_course$Strain == "B6CASTF1" | bw_time_course$Strain == "129SPWKF1"),"Generation"] <- 'F1'
#bw_time_course$Strain[which((bw_time_course$Strain == "PWK/PhJ") & 
#                              (bw_time_course$Sex == "m"))] <- "129SPWKF1"
#get id_vars except week
id_vars <- colnames(bw_time_course)[!grepl("Week",colnames(bw_time_course))]
#get only the week
meas_vars <- colnames(bw_time_course)[grepl("Week",colnames(bw_time_course))]

plyr::count(bw_time_course, vars = c("Strain", "Sex"))

bw_time_course <- bw_time_course[!is.na(bw_time_course$Week.22),]

plyr::count(bw_time_course, vars = c("Strain", "Sex"))

#re-shape the data
bw_time_course_melt <- melt(data = bw_time_course,id.vars = id_vars,
                            measure.vars = meas_vars,variable.name = "Week",
                            value.name = "BW_value")

#remove "week" from the variable name 
bw_time_course_melt$Week <- gsub("Week.","",bw_time_course_melt$Week)
bw_time_course_melt$Week <- as.numeric(bw_time_course_melt$Week)
bw_time_course_melt <- bw_time_course_melt[!is.na(as.numeric(bw_time_course_melt$BW_value)),]

#create a new column for the initial body weight next to bw
bw_initial <- ddply(bw_time_course_melt,c("Pyrat_ID"),function(x)x$BW_value[x$Week==6])

bw_time_course_melt$bw_initial <- 
  bw_initial$V1[match(bw_time_course_melt$Pyrat_ID,bw_initial$Pyrat_ID)]

#create a new column for the final body weight next to bw
bw_final <- dlply(bw_time_course_melt,c("Pyrat_ID"),function(y)y$BW_value[y$Week==24])
bw_final <- ldply(bw_final, data.frame)
colnames(bw_final)[2] <- "V1"
bw_time_course_melt$bw_final <- bw_final$V1[match(bw_time_course_melt$Pyrat_ID,bw_final$Pyrat_ID)]

#calculate percent BW gain
bw_time_course_melt$BW_NormBW6_per <- 100 * bw_time_course_melt$BW_value / bw_time_course_melt$bw_initial
bw_time_course_melt$BW_perGain <- 100*(bw_time_course_melt$BW_value-bw_time_course_melt$bw_initial)/bw_time_course_melt$bw_initial

bw_time_course_melt$Sex <- factor(bw_time_course_melt$Sex,levels = c("m","f"))
bw_time_course_melt$Generation <- factor(bw_time_course_melt$Generation,levels = c("F0","F1","F2"))
bw_time_course_melt$Strain <- factor(bw_time_course_melt$Strain,levels = c("C57BL/6J","129S1/SvImJ", "CAST/EiJ","PWK/PhJ", 
                                                                           "B6CASTF1","129SPWKF1","B6CAST-129SPWK-F2"))


#####Body weight graphs######
######initial body weight histo######
bw1 <- ggplot(bw_time_course_melt[grepl("6",bw_time_course_melt$Week) 
                           & grepl("m",bw_time_course_melt$Sex) & !grepl("16",bw_time_course_melt$Week),], 
       aes(x = reorder(Pyrat_ID, bw_initial), y = bw_initial, color = Strain, fill = Strain)) +
  geom_bar(stat="identity")+
  scale_color_manual(name = "Strain",values = strainColors)+
  scale_fill_manual(name = "Strain",values = strainColors)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30))+
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
  ylab("Initial Body Weight (g)")+
  xlab("Strain")+
  ggtitle("Males")

bw1
ggsave(useDingbats = FALSE, bw1,file="./Plots/02_bw_initial_males.pdf",
       width = 7, height = 2, units = c("in"))


bw2 <- ggplot(bw_time_course_melt[grepl("6",bw_time_course_melt$Week) 
                           & grepl("f",bw_time_course_melt$Sex) & !grepl("16",bw_time_course_melt$Week),], 
       aes(x = reorder(Pyrat_ID, bw_initial), y = bw_initial, color = Strain, fill = Strain)) +
  geom_bar(stat="identity")+
  scale_color_manual(name = "Strain",values = strainColors)+
  scale_fill_manual(name = "Strain",values = strainColors)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30))+
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
  ylab("Initial Body Weight (g)")+
  xlab("Strain")+
  ggtitle("Females")

bw2
ggsave(useDingbats = FALSE, bw2,file="./Plots/02_bw_initial_females.pdf",
       width = 7, height = 2, units = c("in"))

bw3 <- ggplot(bw_time_course_melt[grepl("6",bw_time_course_melt$Week) & !grepl("16",bw_time_course_melt$Week),], 
       aes(x = reorder(Pyrat_ID, bw_initial), y = bw_initial, color = Sex, fill = Sex)) +
  geom_bar(stat="identity")+
  scale_color_manual(name = "Sex",values = sexColors)+
  scale_fill_manual(name = "Sex",values = sexColors)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30))+
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
  ylab("Initial Body Weight (g)")+
  xlab("Strain")+
  ggtitle("All")

bw3
ggsave(useDingbats = FALSE, bw3,file="./Plots/02_bw_initial.pdf",
       width = 7, height = 2, units = c("in"))


######final body weight histo######
bw4 <- ggplot(bw_time_course_melt[grepl("24",bw_time_course_melt$Week) 
                             & grepl("m",bw_time_course_melt$Sex),], 
         aes(x = reorder(Pyrat_ID, bw_final), y = bw_final, color = Strain, fill = Strain)) +
    geom_bar(stat="identity")+
    scale_color_manual(name = "Strain",values = strainColors)+
    scale_fill_manual(name = "Strain",values = strainColors)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 70))+
    #scale_alpha_manual(values = c(0.3,1))+
    theme(strip.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylab("Final Body Weight (g)")+
    xlab("Strain")+
    ggtitle("Males")

bw4
ggsave(useDingbats = FALSE, bw4,file="./Plots/02_bw_final_males.pdf",
       width = 7, height = 2, units = c("in"))

  
bw5 <- ggplot(bw_time_course_melt[grepl("24",bw_time_course_melt$Week) 
                             & grepl("f",bw_time_course_melt$Sex),], 
         aes(x = reorder(Pyrat_ID, bw_final), y = bw_final, color = Strain, fill = Strain)) +
    geom_bar(stat="identity")+
    scale_color_manual(name = "Strain",values = strainColors)+
    scale_fill_manual(name = "Strain",values = strainColors)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 70))+
    #scale_alpha_manual(values = c(0.3,1))+
    theme(strip.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylab("Final Body Weight (g)")+
    xlab("Strain")+
    ggtitle("Females")

bw5
ggsave(useDingbats = FALSE, bw5,file="./Plots/02_bw_final_females.pdf",
       width = 7, height = 2, units = c("in"))

ggplot(bw_time_course_melt[grepl("24",bw_time_course_melt$Week),], 
         aes(x = reorder(Pyrat_ID, bw_final), y = bw_final, color = Strain, fill = Strain)) +
    geom_bar(stat="identity")+
    scale_color_manual(name = "Strain",values = strainColors)+
    scale_fill_manual(name = "Strain",values = strainColors)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 70))+
    #scale_alpha_manual(values = c(0.3,1))+
    theme(strip.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylab("Final Body Weight (g)")+
    xlab("Strain")+
    ggtitle("Final Body Weight All")


  
  
    
ggplot(bw_time_course_melt[grepl("24",bw_time_course_melt$Week),], 
         aes(x = reorder(Pyrat_ID, bw_final), y = bw_final, color = Sex, fill = Sex)) +
    geom_bar(stat="identity")+
    scale_color_manual(name = "Sex",values = sexColors)+
    scale_fill_manual(name = "Sex",values = sexColors)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 70))+
    #scale_alpha_manual(values = c(0.3,1))+
    theme(strip.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylab("Final Body Weight (g)")+
    xlab("Strain")+
    ggtitle("All")
  
  

bw_time_course_melt$Generation_Sex <- paste0(bw_time_course_melt$Generation,"_",bw_time_course_melt$Sex)
bw_time_course_melt$Generation_Sex <- factor(bw_time_course_melt$Generation_Sex,levels = c("F0_m","F0_f","F1_m","F1_f","F2_m","F2_f"))
bw_time_course_melt$Strain_Sex <- paste0(bw_time_course_melt$Strain,"_",bw_time_course_melt$Sex)
bw_time_course_melt$Strain_Sex <- factor(bw_time_course_melt$Strain_Sex,
                                     levels = c("C57BL/6J_m","C57BL/6J_f",
                                                "129S1/SvImJ_m","129S1/SvImJ_f", 
                                                "CAST/EiJ_m","CAST/EiJ_f",
                                                "PWK/PhJ_m", "PWK/PhJ_f",
                                                "B6CASTF1_m","B6CASTF1_f",
                                                "129SPWKF1_m","129SPWKF1_f",
                                                "B6CAST-129SPWK-F2_m","B6CAST-129SPWK-F2_f"))


#COLORCODE <- c("#FF9300","#107F40","#0000FE","#FF2600", "#942192","#929292")
#strainColors <- COLORCODE
#######boxplots#######
ggplot(bw_time_course_melt[grepl("6",bw_time_course_melt$Week) 
                           & grepl("m",bw_time_course_melt$Sex) & !grepl("16",bw_time_course_melt$Week),], 
       aes(x = Generation, y = bw_initial, group = Generation)) +
  geom_boxplot(mapping = aes(group = Generation), color = c("grey","grey","grey"),outlier.shape = NA)+
  #geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.75), 
  #dotsize = 0.85, aes(col= Strain,fill = Strain, group = Strain_Diet, alpha=Diet))+
  #scale_color_manual(values= c("grey","grey","grey"))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8,dodge.width = 0.8), 
              size=2.5, aes(colour = Strain, group = Generation))+
  scale_color_manual(name = "Strain",values = strainColors)+
  #scale_alpha_manual(values = c(0.4,1))+
  #theme(strip.background = element_blank(),panel.background = element_blank(),axis.line = element_line())+
  ylab("Initial Body Weight (g)")+
  xlab("Generation")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.background = element_blank(),panel.background = element_blank(),axis.line = element_line())


ggplot(bw_time_course_melt[grepl("6",bw_time_course_melt$Week) 
                           & grepl("f",bw_time_course_melt$Sex) & !grepl("16",bw_time_course_melt$Week),], 
       aes(x = Generation, y = bw_initial, group = Generation)) +
  geom_boxplot(mapping = aes(group = Generation), color = c("grey","grey","grey"),outlier.shape = NA)+
  #geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.75), 
  #dotsize = 0.85, aes(col= Strain,fill = Strain, group = Strain_Diet, alpha=Diet))+
  #scale_color_manual(values= c("grey","grey","grey"))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8,dodge.width = 0.8), 
              size=2.5, aes(colour = Strain, group = Generation))+
  scale_color_manual(name = "Strain",values = strainColors)+
  #scale_alpha_manual(values = c(0.4,1))+
  #theme(strip.background = element_blank(),panel.background = element_blank(),axis.line = element_line())+
  ylab("Initial Body Weight (g)")+
  xlab("Generation")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.background = element_blank(),panel.background = element_blank(),axis.line = element_line())


aggregate(bw_initial ~ Generation_Sex, data = bw_time_course_melt, max) #30.40
aggregate(bw_initial ~ Generation_Sex, data = bw_time_course_melt, min) #10.30
aggregate(bw_final ~ Generation_Sex, data = bw_time_course_melt, max) #66.55
aggregate(bw_final ~ Generation_Sex, data = bw_time_course_melt, min) #13.12


#######BW timecourse#######
#COLORCODE <- c("#FF9300","#107F40","#0000FE","#FF2600", "#942192","#929292")
#strainColors <- COLORCODE
#single values
ggplot(bw_time_course_melt[!is.na(bw_time_course_melt$bw_final),], 
       aes(x = Week, y = BW_perGain, color = Strain)) +
  geom_point(aes(alpha = Sex), shape = 1, size = 0.5, stroke =0.35) +
  geom_path(aes(color = Strain, alpha = Sex, group = Pyrat_ID), size =0.8) +
  #stat_summary(geom = "line", fun.y = "median", mapping = aes(alpha = Diet), lwd=0.5, size = 0.5)+
  scale_color_manual(name = "Strain",values = strainColors) +
  facet_wrap(~Strain, nrow = 1, scales = "fixed") +
  theme(strip.background = element_blank(),panel.background = element_blank(),axis.line = element_line())+
  scale_alpha_manual(values = c(1,0.3))+
  scale_x_continuous(breaks = seq(6, 24, by = 6))+
  ylab("Body Weight Gain [% Start]")+
  xlab("Age [Weeks]") 




ggplot(bw_time_course_melt, 
       aes(x = Week, y = BW_perGain, color = Strain)) +
  geom_point(shape = 1, size = 0.5, stroke =0.35) +
  geom_path(aes(color = Strain, group = Pyrat_ID), size =0.8) +
  #stat_summary(geom = "line", fun.y = "median", mapping = aes(alpha = Diet), lwd=0.5, size = 0.5)+
  scale_color_manual(name = "Strain",values = strainColors) +
  facet_grid(Sex ~ Strain, scales = "fixed")+
  #facet_wrap(~Sex+Strain, nrow = 2, scales = "fixed") +
  theme(strip.background = element_blank(),panel.background = element_blank(),axis.line = element_line())+
  #scale_alpha_manual(values = c(1,0.3))+
  scale_x_continuous(breaks = seq(6, 24, by = 6))+
  ylab("Body Weight Gain [% Start]")+
  xlab("Age [Weeks]")+ 
  guides(color = "none")


###################with the PWK mice from Step 1, because in step 2 the PWK were in reality 129-PWK cross##################
  data_wide <- read.csv(file = 'Aim2_Step2_manuscript/Data/allData_wide_clean_noOutlier_pwk.csv')
  df <- data_wide[,c(1:4,16:26)]
  id_vars <- colnames(df)[!grepl("BW",colnames(df))]
  #get only the week
  meas_vars <- colnames(df)[grepl("BW",colnames(df))]
  
  
  #re-shape the data
  df_melt <- melt(data = df,id.vars = id_vars,
                  measure.vars = meas_vars,variable.name = "Variable",
                  value.name = "Value")
  df_melt$Week <- gsub('[^_]*_', '',df_melt$Variable)
  
  df_melt$Week <- as.numeric(as.character(df_melt$Week))
  df_melt$Sex <- factor(df_melt$Sex,levels = c("m","f"))
  df_melt$Strain <- factor(df_melt$Strain,levels = c("C57BL/6J","129S1/SvImJ", "CAST/EiJ","PWK/PhJ", 
                                                     "B6CASTF1","129SPWKF1","B6CAST-129SPWK-F2"))
  ggplot(df_melt, 
         aes(x = Week, y = Value, color = Strain)) +
    geom_point(aes(alpha = Sex), shape = 1, size = 0.5, stroke =0.35) +
    geom_path(aes(color = Strain, alpha = Sex, group = ID), size =0.8) +
    #stat_summary(geom = "line", fun.y = "median", mapping = aes(alpha = Diet), lwd=0.5, size = 0.5)+
    scale_color_manual(name = "Strain",values = strainColors) +
    facet_wrap(~Strain, nrow = 1, scales = "fixed") +
    scale_x_continuous(breaks = seq(6, 24, by = 6))+
    theme(strip.background = element_blank(),panel.background = element_blank(),axis.line = element_line())+
    scale_alpha_manual(values = c(1,0.3))+
    ylab("Body Weight Gain [% Start]")+
    xlab("Age [Weeks]") 
  
  bw <- ggplot(df_melt, 
               aes(x = Week, y = Value, color = Strain)) +
    #geom_point(shape = 1, size = 0.5, stroke =0.35) +
    geom_path(aes(color = Strain, group = ID), size =0.8) +
    #stat_summary(geom = "line", fun.y = "median", mapping = aes(alpha = Diet), lwd=0.5, size = 0.5)+
    scale_color_manual(name = "Strain",values = strainColors) +
    facet_grid(Sex ~ Strain, scales = "fixed")+
    scale_x_continuous(breaks = seq(6, 24, by = 6))+
    theme(strip.background = element_blank(),panel.background = element_blank(),axis.line = element_line())+
    scale_alpha_manual(values = c(1,0.3))+
    ylab("Body Weight Gain [% Start]")+
    xlab("Age [Weeks]")+
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

