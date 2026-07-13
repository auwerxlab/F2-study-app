rm(list=ls())
source(file = "./Scripts/Proteomics/00-common.R")

################## Libraries ##################
library(dplyr)
library(readxl)
################ data loading #################
Prot<-read.csv(file="./Aim2_Step2_manuscript/Data/proteomics/Auwerx_Liver_4databasesearch_PG_FINAL(in).csv")
metadata<-read.csv(file = "./Aim2_Step2_manuscript/Data/Sample_metadata.csv")

#####check how many samples we sent and how many we have data from#########
file <- "./Aim2_Step2_manuscript/Data/proteomics/F2_MICE_LIST_CoonLab_randomized_weight_20mg_v2.xlsx"
sheets <- excel_sheets(file)
samples_list <- lapply(sheets, function(s) {
  read_excel(file, sheet = s, col_types = "text")
}) %>% 
  bind_rows()

samples_missing <- samples_list %>%
  filter(`Weight (20mg)` == "no tube")

samples_sent <- samples_list %>%
  filter(`Weight (20mg)` != "no tube")

samples_received <- colnames(Prot)[-c(1:3)]
samples_received <- gsub("^HDP_", "HDP-", samples_received)
sum(samples_sent$Id %in% samples_received)
#536
missing_data <- samples_sent$Id[!samples_sent$Id %in% samples_received]
#14 mice are missing from the data:  "HDP-015291", "HDP-012309", "HDP-015829", "HDP-014078", "HDP-014787", 
#"HDP-014934", "HDP-015099" ,"HDP-014904","HDP-016188", "HDP-014777" ,"HDP-013231", "HDP-014102", "HDP-013614" ,"HDP-013229"


###################### QC #####################
colnames(Prot)
count_missing_row<- cbind(Prot[,grep("HDP",colnames(Prot),invert = T)],
                          N_missing=rowSums(is.na(Prot[,grep("HDP",colnames(Prot))])))

hist(count_missing_row$N_missing)

count_missing_col<- data.frame(ID=colnames(Prot[,grep("HDP",colnames(Prot))]),
                               N_missing=colSums(is.na(Prot[,grep("HDP",colnames(Prot))])))

hist(count_missing_col$N_missing)

test<-count_missing_row[count_missing_row$N_missing<(539/4),]
test$Gene.Symbol
# FILTER: min half genes detected per sample (only one is bad), and min half samples detected per gene
Prot_filt<-Prot[Prot$protein.groups %in% count_missing_row[count_missing_row$N_missing<(539/2),]$protein.groups,
                colnames(Prot)%in%c(colnames(Prot)[grep("HDP",colnames(Prot),invert = T)],
                                    count_missing_col[count_missing_col$N_missing<6000,]$ID)]


count_missing_row_f<- cbind(Prot_filt[,grep("HDP",colnames(Prot_filt),invert = T)],
                            N_missing=rowSums(is.na(Prot_filt[,grep("HDP",colnames(Prot_filt))])))

hist(count_missing_row_f$N_missing)

count_missing_col_f<- data.frame(ID=colnames(Prot_filt[,grep("HDP",colnames(Prot_filt))]),
                                 N_missing=colSums(is.na(Prot_filt[,grep("HDP",colnames(Prot_filt))])))

hist(count_missing_col_f$N_missing)
length(unique(Prot_filt$Gene.Symbol))
colnames(Prot_filt)<-gsub("_","-", colnames(Prot_filt))
metadata<-metadata[metadata$ID%in%colnames(Prot_filt),]

Prot_filt_m<-reshape2::melt(Prot_filt,id.vars=c("Gene.Symbol","protein.groups","PG.FastaHeaders"),variable.name="ID",value.name="Prot_quantif")
Prot_filt_m<-merge(metadata,Prot_filt_m)

save(Prot_filt_m,Prot_filt,metadata,file="./Data/Filtered_prots.Rdata")




