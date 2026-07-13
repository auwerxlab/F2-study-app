####################################### Libraries ###########################################
library(reshape2)
library(plyr) # For reformatting data
library(ggplot2) # For plotting
library(RColorBrewer) # For color palettes
library(scales) # For "oob =squish" Makes out-of-bounds values not considered as NA in the color scale
library(cowplot)
library(data.table)
library(reshape2)

####################################### Palettes ###########################################

LEVELS_STRAINS <- c("C57BL/6J","129S1/SvImJ", "CAST/EiJ","PWK/PhJ", "B6CASTF1","129SPWKF1","B6CAST-129SPWK-F2")
#PALETTE <- "Dark2"
# C57BL/6J		000000
# DBA	AA7942
# AJ	942192
# 129S1 FF9300
# WSB	FF2600
# CAST107F40
# PWK	0000FE
COLORCODE <- c("#000000","#FF9300","#107F40","#0000FE","#FF2600", "#942192","#C0C0C0")
strainColors <- COLORCODE
names(strainColors) <- LEVELS_STRAINS

sexColors <- c("m" = "#00A2FF", "f" = "#FF644E")

Heatmap_palette <- c(rev(brewer.pal(7,"Blues")),"white","white","white",brewer.pal(7,"Reds"))

####################################### Graphical Themes ###########################################

theme_graphs <-theme_bw()+
  theme(axis.text=element_text(size=12,color="black"),
        plot.title=element_text(face = "bold",size=14,hjust = 0.5),
        plot.subtitle = element_text(size=13,hjust = 0.5),
        axis.title=element_text(size=12,hjust = 0.5),
        legend.title = element_blank(),
        legend.text=element_text(size=13),
        legend.position = "bottom",
        strip.text = element_text(face = "bold",size=12),
        strip.background = element_blank(),
        #strip.background = element_rect(fill = "gray88" ,colour="black", size = 0.3),
        element_line(color = "gray25", linewidth = 0.5),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing = unit(0, "lines"))


theme_heatmaps <-theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),#Center the title
    axis.text.y = element_text(size=9,face="bold"),
    axis.text.x = element_text(size=9.5,angle = 45,hjust = 1, vjust = 1,face="bold"),
    legend.title = element_text(size=10,face="bold"),
    legend.text = element_text(size=10,hjust = 0.5),
    strip.text = element_text(face = "bold",size=12),
    strip.background = element_blank(),
    element_line(color = "gray25", size = 0.5),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
    panel.spacing = unit(0.1, "lines"),
    panel.grid = element_blank()
  )

