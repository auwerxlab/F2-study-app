rm(list=ls())
source(file = "./Scripts/P00-Proteomics-Colorcode_themes_common_functions.R")
################## Libraries ##################
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tidyHeatmap)
library(tidyverse)
library(patchwork) 
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(ggpubr)
################## data loading ###############
load(file="./Data/Filtered_prots.Rdata")
liver_dge <- readRDS("./Data/liver_dge.RDS")
ENS_to_gene_name<-read.csv(file="./Data/geneConversionTable.csv")
ENS_to_gene_name<-ENS_to_gene_name[,c("gene_id","gene_name")]
colnames(ENS_to_gene_name)<-c("gene_id","Gene.Symbol")
###############################################

RNA<-as.data.frame(liver_dge$cpm)[-(1:3),]
RNA$gene_id<-row.names(RNA)
RNA_m<-reshape2::melt(RNA,id.vars=c("gene_id"), variable.name="Sample_ID",value.name = "log2CPM")
RNA_m<-merge(ENS_to_gene_name,RNA_m)
RNA_m<-RNA_m[RNA_m$Gene.Symbol%in%Prot_filt_m$Gene.Symbol,]
RNA_m<-merge(metadata[,c("Sample_ID","ID")],RNA_m)


Prot_RNA<-merge(Prot_filt_m,RNA_m)
Prot_RNA <- Prot_RNA %>%
  dplyr::rename(log2CPM_RNA = log2CPM)
Prot_RNA <- Prot_RNA %>%
  mutate(log2_Protein = log2(Prot_quantif))

data_wide <- read.csv(file = './Data/allData_wide_clean_noOutlier_liver_cluster.csv')
Prot_RNA <- Prot_RNA %>%
  left_join(data_wide %>% dplyr::select(ID, liver_cluster), by = "ID")

#proteins_per_sample <- Prot_RNA %>%
#  filter(
#    !is.na(log2_Protein),
#    !is.nan(log2_Protein)
#  ) %>%
#  group_by(ID) %>%
#  summarise(
#    n_proteins = n_distinct(gene_id),
#    .groups = "drop"
#  )

mice_per_cluster <- Prot_RNA %>%
  group_by(liver_cluster) %>%
  summarise(
    n_mice = n_distinct(ID),
    .groups = "drop"
  )

#proteins_per_sample <- proteins_per_sample %>%
#  left_join(data_wide %>% select(ID, liver_cluster), by = "ID")


num_rna_prot <- length(unique(Prot_RNA$gene_id))
num_rna_prot
#7921
Prot_RNA_unique <- Prot_RNA %>%
  group_by(ID, gene_id) %>%   # group by mouse and gene
  slice_max(order_by = Prot_quantif, n = 1) %>%  # keep row with max Prot_quantif
  ungroup()


corr_per_sample <- Prot_RNA_unique %>%
  group_by(ID) %>%
  summarise(
    RNA_Protein_Corr = cor(
      log2_Protein,
      log2CPM_RNA,
      method = "pearson",
      use = "complete.obs"
    )
  )
corr_per_sample<-merge(metadata,corr_per_sample)
corr_per_sample <- corr_per_sample %>%
  left_join(data_wide %>% dplyr::select(ID, liver_cluster), by = "ID")
corr_per_sample$RNA_Protein_Corr<-as.numeric(corr_per_sample$RNA_Protein_Corr)
corr_per_sample$Strain<-factor(corr_per_sample$Strain,LEVELS_STRAINS)


ggplot(corr_per_sample, aes(x = as.factor(RNA_Protein_Corr), y = RNA_Protein_Corr, fill=Strain)) +
  geom_col() +
  theme_minimal() +
  scale_fill_manual(values = COLORCODE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "RNA–Protein Correlation per Sample",
    x = "Sample",
    y = "Correlation Coefficient"
  )

range(corr_per_sample$RNA_Protein_Corr)

prot_rna_per_sample_dist <- ggplot(corr_per_sample, aes(x = RNA_Protein_Corr)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black", linewidth = 0.5) +
  theme_classic() +
  labs(
    x = "RNA–Protein Correlation(r)",
    y = "Count"
    #,
    #title = "Distribution of RNA–Protein Correlation per Sample"
  ) +
  theme(
    text = element_text(size = 8),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    plot.title=element_text(size=8),
    legend.title = element_text(size=8),
    legend.text = element_text(size=8),
    line = element_line(linewidth = 0.5)
  )
prot_rna_per_sample_dist
ggsave(useDingbats = FALSE, prot_rna_per_sample_dist,file="./Plots/20-prot_RNA_corr_per_sample_dist_pearson.pdf",
       width = 2, height = 1.8, units = c("in"))

prot_rna_barplot <- corr_per_sample %>%
  arrange(RNA_Protein_Corr) %>%
  ggplot(aes(x = reorder(ID, RNA_Protein_Corr),
             y = RNA_Protein_Corr, fill=Strain)) +
  geom_col() +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title=element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values = COLORCODE)+
  labs(
    title = "RNA–Protein Correlation per Sample (r)",
    x = "Sample",
    y = "Correlation Coefficient"
  )
prot_rna_barplot
ggsave(useDingbats = FALSE, prot_rna_barplot,file="./Plots/20-prot_RNA_corr_barplot_pearson.pdf",
       width = 7, height = 2, units = c("in"))

range(corr_per_sample$RNA_Protein_Corr, na.rm = TRUE)

a <- ggplot(corr_per_sample,
            aes(x = RNA_Protein_Corr, color = factor(liver_cluster))) +
  geom_density(linewidth = 0.5) +
  scale_color_manual(
    name = "Cluster",
    values = c("#0696c7", "#FFD700", "#ff8800", "#6e0280")
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),       # smaller legend keys
    legend.spacing.y = unit(0.1, "cm"),      # reduce spacing between legend items
    legend.box.spacing = unit(0.2, "cm")     # reduce spacing around legend box
  ) +
  labs(
    #title = "RNA–Protein Correlation per Gene",
    x = "Spearman correlation coefficient",
    y = "Density"
  )

a
aggregate(RNA_Protein_Corr ~ liver_cluster,
          data = corr_per_sample,
          FUN = median,
          na.rm = TRUE)

corr_per_gene <- Prot_RNA_unique %>%
  group_by(gene_id) %>%
  summarise(
    RNA_Protein_Corr = cor(log2_Protein, log2CPM_RNA, method="pearson",use = "complete.obs")
  )
corr_per_gene <- corr_per_gene %>%
  left_join(ENS_to_gene_name %>% dplyr::select(gene_id, Gene.Symbol), by = "gene_id")
anyDuplicated(corr_per_gene$Gene.Symbol)

p <- ggplot(corr_per_gene, aes(x = RNA_Protein_Corr)) +
  geom_histogram(
    breaks = seq(-0.7, 0.9, by = 0.1),
    #binwidth = 0.1, 
    #boundary = -0.9,        # ensures the first bin starts at -0.7
    fill = "grey70", color = "black",
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    plot.title=element_text(size=8),
    legend.title = element_text(size=8),
    legend.text = element_text(size=8),
    line = element_line(linewidth = 0.5)
  )+
  labs(
    #title = "Distribution of RNA–Protein Correlation per Gene",
    x = "RNA-protein correlation(r)",
    y = "Genes(n)"
  )
p
ggsave(
  filename = "./Plots/20-RNA_Protein_correlation_per_gene_histo_pearson.pdf",
  plot = p,
  width = 2.2,
  height = 1.8,
  units = "in"
)


gene_counts_cluster <- Prot_RNA_unique %>%
  group_by(liver_cluster, gene_id) %>%
  summarise(n_mice = sum(!is.na(log2_Protein) & !is.na(log2CPM_RNA)), .groups = "drop")

corr_per_gene_cluster <- Prot_RNA_unique %>%
  filter(!is.na(liver_cluster)) %>%                  # remove samples with NA cluster
  group_by(liver_cluster, gene_id) %>%
  filter(!is.na(log2_Protein) & !is.na(log2CPM_RNA)) %>%
  filter(n() >= 5) %>%                               # keep genes with >=10 mice per cluster
  summarise(
    RNA_Protein_Corr = cor(log2_Protein, log2CPM_RNA,
                           method = "spearman",
                           use = "complete.obs"),
    n_mice = n(),
    .groups = "drop"
  )

corr_per_gene_cluster_spear <- Prot_RNA_unique %>%
  filter(!is.na(liver_cluster)) %>%                  # remove samples with NA cluster
  group_by(liver_cluster, gene_id) %>%
  filter(!is.na(log2_Protein) & !is.na(log2CPM_RNA)) %>%
  filter(n() >= 5) %>%                               # keep genes with >=10 mice per cluster
  summarise(
    RNA_Protein_Corr = cor(log2_Protein, log2CPM_RNA,
                           method = "spearman",
                           use = "complete.obs"),
    n_mice = n(),
    .groups = "drop"
  )


corr_per_gene_cluster_spear <- Prot_RNA_unique %>%
  filter(!is.na(liver_cluster)) %>%                  
  group_by(liver_cluster, gene_id) %>%
  filter(!is.na(log2_Protein) & !is.na(log2CPM_RNA)) %>%
  filter(n() >= 5) %>%                               # keep genes with >=5 mice per cluster
  summarise(
    # Run cor.test safely even if data are constant
    test = list(try(cor.test(log2_Protein, log2CPM_RNA,
                             method = "spearman",
                             use = "complete.obs"), silent = TRUE)),
    n_mice = n(),
    .groups = "drop"
  ) %>% 
  mutate(RNA_Protein_Corr = map_dbl(test, ~ if(inherits(.x, "try-error")) NA_real_ else .x$estimate), 
         p_value = map_dbl(test, ~ if(inherits(.x, "try-error")) NA_real_ else .x$p.value)) %>%
  dplyr::select(-test) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))    # FDR correction

corr_per_gene_cluster_pears <- Prot_RNA_unique %>%
  filter(!is.na(liver_cluster)) %>%                  
  group_by(liver_cluster, gene_id) %>%
  filter(!is.na(log2_Protein) & !is.na(log2CPM_RNA)) %>%
  filter(n() >= 5) %>%                               # keep genes with >=5 mice per cluster
  summarise(
    # Run cor.test safely even if data are constant
    test = list(try(cor.test(log2_Protein, log2CPM_RNA,
                             method = "pearson",
                             use = "complete.obs"), silent = TRUE)),
    n_mice = n(),
    .groups = "drop"
  ) %>% 
  mutate(RNA_Protein_Corr = map_dbl(test, ~ if(inherits(.x, "try-error")) NA_real_ else .x$estimate), 
         p_value = map_dbl(test, ~ if(inherits(.x, "try-error")) NA_real_ else .x$p.value)) %>%
  dplyr::select(-test) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

df <- aggregate(RNA_Protein_Corr ~ liver_cluster,
          data = corr_per_gene_cluster_pears,
          FUN = median,
          na.rm = TRUE)

b <- ggplot(df, aes(x = RNA_Protein_Corr,
              y = factor(liver_cluster), fill = factor(liver_cluster))) +
  geom_col() +
  scale_fill_manual(
    name = "Cluster",
    values = c("#0696c7", "#FFD700", "#ff8800", "#6e0280")
  ) +
  labs(
    x = "Median RNA–Protein corr.",
    y = "Cluster",
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),       # smaller legend keys
    legend.spacing.y = unit(0.1, "cm"),      # reduce spacing between legend items
    legend.box.spacing = unit(0.2, "cm")     # reduce spacing around legend box
  )
b
ggsave(
  filename = "./Plots/20-RNA_Protein_correlation_median_per_cluster_pearson.pdf",
  plot = b,
  width = 1.9,
  height = 1,
  units = "in"
)
summary_counts_pears <- corr_per_gene_cluster_pears %>%
  mutate(
    direction = case_when(
      adj_p_value < 0.05 & RNA_Protein_Corr > 0 ~ "positive",
      adj_p_value < 0.05 & RNA_Protein_Corr < 0 ~ "negative",
      TRUE ~ "not_significant"
    )
  ) %>%
  count(liver_cluster, direction, name = "n_genes") %>%
  pivot_wider(
    names_from  = direction,
    values_from = n_genes,
    values_fill = 0
  )



summary_counts_spear <- corr_per_gene_cluster_spear %>%
  mutate(
    direction = case_when(
      adj_p_value < 0.05 & RNA_Protein_Corr > 0 ~ "positive",
      adj_p_value < 0.05 & RNA_Protein_Corr < 0 ~ "negative",
      TRUE ~ "not_significant"
    )
  ) %>%
  count(liver_cluster, direction, name = "n_genes") %>%
  pivot_wider(
    names_from  = direction,
    values_from = n_genes,
    values_fill = 0
  )

genes_per_cluster <- corr_per_gene_cluster %>%
  group_by(liver_cluster) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    .groups = "drop"
  )

genes_per_cluster


p <- ggplot(corr_per_gene_cluster_pears,
            aes(x = RNA_Protein_Corr, color = factor(liver_cluster))) +
  geom_density(linewidth = 0.5) +
  scale_color_manual(
    name = "Cluster",
    values = c("#0696c7", "#FFD700", "#ff8800", "#6e0280")
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),       # smaller legend keys
    legend.spacing.y = unit(0.1, "cm"),      # reduce spacing between legend items
    legend.box.spacing = unit(0.2, "cm")     # reduce spacing around legend box
  ) +
  labs(
    #title = "RNA–Protein Correlation per Gene",
    x = "Pearson correlation coefficient",
    y = "Density"
  )

p


ggsave(
  filename = "./Plots/20-RNA_Protein_correlation_per_gene_density_by_cluster_pearson.pdf",
  plot = p,
  width = 2.3,
  height = 1.5,
  units = "in"
)

q <- ggplot(corr_per_gene_cluster_spear,
            aes(x = RNA_Protein_Corr, color = factor(liver_cluster))) +
  geom_density(linewidth = 0.5) +
  scale_color_manual(
    name = "Cluster",
    values = c("#0696c7", "#FFD700", "#ff8800", "#6e0280")
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),       # smaller legend keys
    legend.spacing.y = unit(0.1, "cm"),      # reduce spacing between legend items
    legend.box.spacing = unit(0.2, "cm")     # reduce spacing around legend box
  ) +
  labs(
    #title = "RNA–Protein Correlation per Gene",
    x = "Spearman correlation coefficient",
    y = "Density"
  )

q


ggsave(
  filename = "./Plots/20-RNA_Protein_correlation_per_gene_density_by_cluster_spearman.pdf",
  plot = q,
  width = 2.3,
  height = 1.5,
  units = "in"
)



ggplot(corr_per_gene_cluster, aes(x = RNA_Protein_Corr)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 0.1,
    boundary = -0.7,
    fill = "grey80",
    color = "black"
  ) +
  geom_density(color = "red", linewidth = 1) +
  facet_wrap(~ liver_cluster, ncol = 2) +
  theme_minimal() +
  labs(
    title = "RNA–Protein Correlation Distribution per Gene by Liver Cluster",
    x = "Correlation coefficient",
    y = "Density"
  )

ggplot(corr_per_gene_cluster,
       aes(x = RNA_Protein_Corr, 
           #fill = factor(liver_cluster), 
           color = factor(liver_cluster))) +
  scale_color_manual(name = "Cluster",values = c("#0696c7","#FFD700","#ff8800","#6e0280"))+
  #scale_fill_manual(name = "Cluster",values = c("#0696c7","#FFD700","#ff8800","#6e0280"))+
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 0.1,
    boundary = -0.7,
    alpha = 0.1,
    position = "identity"
  ) +
  geom_density(linewidth = 1.1) +
  theme_minimal() +
  labs(
    title = "RNA–Protein Correlation Distribution by Liver Cluster",
    x = "Correlation coefficient",
    y = "Density",
    #fill = "Cluster",
    color = "Cluster"
  )

corr_wide <- corr_per_gene_cluster_spear %>%
  dplyr::select(gene_id, liver_cluster, RNA_Protein_Corr) %>%
  pivot_wider(
    names_from = liver_cluster,
    values_from = RNA_Protein_Corr
  )
anyDuplicated(corr_wide$gene_id)

pairs <- tibble(
  x = c("1", "1", "1"
        #, "2", "2","3"
        ),
  y = c("2", "3", "4"
        #, "3", "4","4"
        )
)

#plot_df <- pairs %>% 
#  mutate(pair = paste(x, "vs", y)) %>%
#  left_join(
#    corr_wide %>% mutate(row = row_number()),
#    by = character()
#  )
plot_df <- map2_dfr(pairs$x, pairs$y, ~ {
  tibble(
    gene_id = corr_wide$gene_id,       # keep gene_id
    x_val   = corr_wide[[.x]],
    y_val   = corr_wide[[.y]],
    pair    = paste(.x, "vs", .y)
  )
})
plot_df <- plot_df %>%
  mutate(category = case_when(
    abs(x_val) <= 0.3 & abs(y_val) <= 0.3 ~ "both_low",
    x_val > 0.3 & y_val > 0.3 ~ "both_positive",
    x_val < -0.3 & y_val < -0.3 ~ "both_negative",
    x_val > 0.3 & y_val >= -0.3 & y_val <= 0.3 ~ "x_positive",
    y_val > 0.3 & x_val >= -0.3 & x_val <= 0.3 ~ "y_positive",
    x_val < -0.3 & y_val >= -0.3 & y_val <= 0.3 ~ "x_negative",
    y_val < -0.3 & x_val >= -0.3 & x_val <= 0.3 ~ "y_negative",
    x_val < -0.3 & y_val > 0.3 ~ "x_negative_y_positive",
    x_val > 0.3 & y_val < -0.3 ~ "x_positive_y_negative",
    TRUE ~ "other"   # fallback if any point doesn't fit above
  ))

plot_df <- plot_df %>%
  filter(!is.na(y_val))

category_counts <- plot_df %>%
  group_by(pair, category) %>%
  summarise(n_points = n(), .groups = "drop")
# Optional: order categories
category_levels <- c("both_low",
  "both_positive", "x_positive", "y_positive", "x_positive_y_negative",
  "both_negative", "x_negative", "y_negative", "x_negative_y_positive"
)
category_counts$category <- factor(category_counts$category, levels = category_levels)
b <- ggplot(category_counts, aes(x = pair, y = n_points, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(
    "both_low"              ="#D3D3D3",
    "both_positive"         = "#D73027",
    "x_positive"            = "#FDAE61",
    "y_positive"            = "#E78AC3" ,
    "x_positive_y_negative" = "#6e0280",
    "both_negative"         = "#4575B4",
    "x_negative"            = "#9ECAE1",
    "y_negative"            = "#66C2A5",
    "x_negative_y_positive" = "#FFD700"
  )) +
  labs(x = "Pair", y = "Number of RNA-protein pairs", fill = "Category") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )
b
# Save to PDF
ggsave("./Plots/20-RNA_prot_corr_cluster_barplot.pdf", plot = b, width = 3.3, height = 2.1)

# Save to PNG
ggsave("./Plots/20-RNA_prot_corr_cluster_barplot.png", plot = b, width = 3.5, height = 2.5, dpi = 300)

# Filter out low_corr points
plot_df_filtered <- plot_df %>%
  filter(category != "both_low", !is.na(y_val))

ordered_categories <- c(
  "both_positive",
  "x_positive",
  "y_positive",
  "x_positive_y_negative",
  "both_negative",
  "x_negative",
  "y_negative"
)
plot_df_filtered <- plot_df_filtered %>%
  mutate(category = factor(category, levels = ordered_categories))
category_colors <- c(
  #"both_low"              = "white",
  "both_positive"         = "#D73027",  # red
  "both_negative"         = "#4575B4",  # blue
  "x_positive"            = "#FDAE61",  # orange
  "y_positive"            = "#E78AC3" ,  # magenta
  "x_negative"            = "#9ECAE1",  # light blue / purple
  "y_negative"            = "#66C2A5",  # green
  #"x_negative_y_positive" = "#FFD700",  # turquoise
  "x_positive_y_negative" = "#6e0280"   # magenta
)


# Build the plot
p <- ggplot(plot_df_filtered, aes(x_val, y_val, color = category)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray50") +
  facet_wrap(~ pair, scales = "free") +
  scale_color_manual(values = category_colors) +
  labs(x = "", y = "") +
  theme_bw(base_size = 8) +   # base text size 8
  theme(
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
# Save to PDF
ggsave("./Plots/20-RNA_prot_corr_cluster_scatter_plot.pdf", plot = p, width = 8, height = 2.2)

# Save to PNG
ggsave("./Plots/20-RNA_prot_corr_cluster_scatter_plot.png", plot = p, width = 8, height = 2.2, dpi = 300)

# Scatter plot with density contours
a <- ggplot(plot_df_filtered, aes(x_val, y_val, color = category)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_density_2d(color = "black", alpha = 0.5) +
  #geom_density_2d_filled(alpha = 0.2)+# add 2D density contours
  geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray50") +
  facet_wrap(~ pair, scales = "free") +
  scale_color_manual(values = c(
    "both_positive"         = "#D73027",
    "both_negative"         = "#4575B4",
    "x_positive"            = "#FDAE61",
    "y_positive"            = "#E78AC3" ,  # magenta
    "x_negative"            = "#9ECAE1",
    "y_negative"            = "#66C2A5",  # green
    "x_negative_y_positive" = "#FFD700",
    "x_positive_y_negative" = "#6e0280"
  )) +
  labs(x = "", y = "") +
  theme_bw(base_size = 8) +
  theme(
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
a
# Save to PDF
ggsave("./Plots/20-RNA_prot_corr_cluster_scatter_plot_density.pdf", plot = a, width = 8, height = 2.2)

# Save to PNG
ggsave("./Plots/20-RNA_prot_corr_cluster_scatter_plot_density.png", plot = a, width = 8, height = 2.2, dpi = 300)

# Scatter plot with density contours all points
b <- ggplot(plot_df, aes(x_val, y_val, color = category)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_density_2d(color = "black", alpha = 0.5) +
  #geom_density_2d_filled(alpha = 0.2)+# add 2D density contours
  geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray50") +
  facet_wrap(~ pair, scales = "free") +
  scale_color_manual(values = c(
    "both_low"              ="#D3D3D3",
    "both_positive"         = "#D73027",
    "both_negative"         = "#4575B4",
    "x_positive"            = "#FDAE61",
    "y_positive"            = "#E78AC3" ,  # magenta
    "x_negative"            = "#9ECAE1",
    "y_negative"            = "#66C2A5",  # green
    "x_negative_y_positive" = "#FFD700",
    "x_positive_y_negative" = "#6e0280"
  )) +
  labs(x = "", y = "") +
  theme_bw(base_size = 8) +
  theme(
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
b
# Save to PDF
ggsave("./Plots/20-RNA_prot_corr_cluster_scatter_plot_density_all_points.pdf", plot = b, width = 8, height = 2.2)

# Save to PNG
ggsave("./Plots/20-RNA_prot_corr_cluster_scatter_plot_density_all_points.png", plot = b, width = 8, height = 2.2, dpi = 300)


# Define the pair and categories of interest
pair_of_interest <- "1 vs 4"
categories_of_interest <- c("x_positive", "y_positive", "y_negative")

# Filter the dataframe
plot_df_filtered <- plot_df %>%
  filter(pair == pair_of_interest,
         category %in% categories_of_interest)
plot_df_filtered <- plot_df_filtered %>%
  mutate(corr_value = case_when(
    category == "x_positive" ~ x_val,
    category %in% c("y_positive", "y_negative") ~ y_val
  ))
plot_df_filtered <- plot_df_filtered %>%
  left_join(ENS_to_gene_name %>% dplyr::select(gene_id, Gene.Symbol), by = "gene_id")


# Split by category
genes_by_category <- plot_df_filtered %>%
  group_by(category) %>%
  summarise(
    gene_list = list(tibble(
      gene_id = Gene.Symbol,
      correlation = corr_value
    )),
    .groups = "drop"
  )
gsea_results <- map(genes_by_category$gene_list, function(df_cat){
  # create ranked vector
  ranked_genes <- setNames(df_cat$correlation, df_cat$gene_id)
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # run pre-ranked GSEA
  GSEA(
    ranked_genes,
    TERM2GENE = msigdbr(species = "Mus musculus") %>%
      dplyr::select(gs_name, gene_symbol),
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
})
summary(plot_df_filtered$corr_value)

genes_by_category <- plot_df_filtered %>%
  group_by(category) %>%
  summarise(
    gene_list = list(unique(gene_id)),
    .groups = "drop"
  )
genes_by_category <- plot_df_filtered %>%
  group_by(category) %>%
  summarise(
    genes = list(unique(dplyr::select(cur_data(), gene_id, Gene.Symbol))),
    .groups = "drop"
  )

mouse_gene_sets <- msigdbr(species = "Mus musculus") %>%
  dplyr::select(gs_name, ensembl_gene)

mouse_gene_sets_filtered <- mouse_gene_sets %>%
  filter(
    grepl("GOBP|HALLMARK|KEGG", gs_name, ignore.case = TRUE)
  )

names_vec <- genes_by_category$category

ora_results <- lapply(genes_by_category$genes, function(df) {
  genes <- unique(df$gene_id)
  
  enricher(
    gene          = genes,
    TERM2GENE     = mouse_gene_sets_filtered,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
})


#names(ora_results) <- genes_by_category$category
names(ora_results) <- names_vec
names(ora_results) <- c("cluster1_positive", "cluster4_negative", "cluster4_positive")

library(clusterProfiler)
library(org.Mm.eg.db)

ora_results_symbol <- lapply(ora_results, function(x){
  if (is.null(x)) return(NULL)
  setReadable(x, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
})

# how many top pathways per category to show
top_n <- 20

# extract tables from enrichment results
df_all <- map2_df(
  ora_results_symbol,
  names(ora_results_symbol),
  ~{
    if (is.null(.x) || nrow(.x@result) == 0) return(NULL)
    
    .x@result %>%
      mutate(
        category = .y,
        neg_log10_padj = -log10(p.adjust)
      )
  }
)

df_all <- df_all %>%
  mutate(
    Description = Description %>%
      str_replace_all("NUCLEAR_TRANSCRIBED_", "") %>%
      str_replace_all("CONTAINING_COMPOUND_", "") %>%
      str_replace_all("NON_MEMBRANE_BOUNDED_", "") %>%
      str_replace_all("REGULATION_OF_", "") %>%
      str_replace_all("ESTABLISHMENT_OF_", "") %>%
      str_replace_all("ENERGY_DERIVATION_BY_", "") %>%
      str_replace_all("GENERATION_OF_", "") %>%
      str_replace_all("_BIOGENESIS", "") %>%
      str_replace_all("GOBP_", "") %>%
      str_replace_all("KEGG_", "") %>%
      str_replace_all("HALLMARK_", "") %>%
      str_replace_all("METABOLIC_PROCESS", "METABOLISM") %>%
      str_replace_all("CATABOLIC_PROCESS", "CATABOLISM") %>%
      str_replace_all("BIOSYNTHETIC_PROCESS", "BIOSYNTHESIS") %>%
      str_replace_all("SIGNAL_TRANSDUCTION", "SIGNALING") %>%
      str_replace_all("_", " ")
  )

df_all <- df_all %>%
  group_by(category, Description) %>%
  slice_max(neg_log10_padj, n = 1, with_ties = FALSE) %>%
  ungroup()

# pick top pathways in each category by adjusted p-value
df_top <- df_all %>%
  group_by(category) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  ungroup()


# reorder descriptions within categories for nicer plotting
df_top <- df_top %>%
  mutate(
    Description = factor(Description, levels = unique(Description[order(category, neg_log10_padj)]))
  )

library(tidytext)
df_top <- df_top %>%
  mutate(
    Description_reordered = reorder_within(Description, neg_log10_padj, category)
  )
df_top <- df_top %>%
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]))
  )


g<- ggplot(df_top, aes(
  x = neg_log10_padj,
  y = Description_reordered,
  size = GeneRatio_num
)) +
  geom_point() +
  facet_wrap(~ category, scales = "free") +
  scale_y_reordered() +
  labs(
   # title = "Top enriched pathways per category",
    x = "-log10 adjusted p-value",
    y = "Pathway",
    size = "GeneRatio"
    ) +
  theme_bw()+
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )+
  scale_size_continuous(range = c(1, 4))
g
ggsave("./Plots/20-Top_enriched_pathways_per_category.pdf", plot = g, width = 11.5, height = 3.6)

genes_per_pathway <- df_top %>%
  filter(p.adjust < 0.05) %>%                     # keep significant pathways
  dplyr::select(category, Description, geneID) %>%
  separate_rows(geneID, sep = "/") %>%            # split into one gene per row
  group_by(category, Description) %>%
  slice_head(n = 20) %>%                            # << number of genes per pathway
  ungroup()

genes_per_pathway <- df_top %>%
  filter(p.adjust < 0.05) %>%                     # keep only significant pathways
  dplyr::select(category, Description, geneID) %>%
  separate_rows(geneID, sep = "/") %>%            # split /-separated genes into rows
  distinct(category, Description, geneID)  

genes_per_category <- genes_per_pathway %>%
  dplyr::select(category, geneID) %>%   # keep only these columns
  distinct()  

selected_genes <- genes_per_category$geneID

corr_changes <- corr_wide %>%
  mutate(
    corr_range = apply(across(-gene_id), 1, function(x)
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  ) %>%
  arrange(desc(corr_range))

corr_changes <- corr_changes %>%
  left_join(
    ENS_to_gene_name %>% dplyr::select(gene_id, Gene.Symbol),
    by = "gene_id"
  )
anyDuplicated(corr_changes$Gene.Symbol)

top_genes <- corr_changes %>%
  arrange(desc(corr_range)) %>%
  slice_head(n = 10) %>%
  pull(Gene.Symbol)



goi <- c("Ndufa13", "Fabp2", "Crot", "Col1a1", "Col5a1", "Ndufv3", "Kynu","Eci3")
goi <- c ("Phka2", "Col1a1", "Lpl", "Ndufa5")
goi <- c("Ndufa2","Sdha","Ndufa6","Mrps7","Cox7b")
goi <- c ("Fabp2", "Pdha1", "Acsl4")
goi <- c ( "Acsl4","Fabp2","Col1a1", "Ndufa2")

#New analysis
goi <- c ( "Plin4","Gpnmb","Fabp5", "Gstm1","Crat")
goi <- c ( "Cela1","Scd1","Hnrnpf", "Thbs1")

plot_df <- Prot_RNA_unique %>%
  filter(
    Gene.Symbol %in% goi,
    !is.na(liver_cluster),
    !is.na(log2_Protein),
    !is.na(log2CPM_RNA)
  )
#plot_df$Gene.Symbol <- factor(
#  plot_df$Gene.Symbol,
#  levels = c("Acsl4", "Fabp2", "Col1a1", "Ndufa2")
#)
p <- ggplot(
  plot_df,
  aes(x = log2CPM_RNA, y = log2_Protein, color = factor(liver_cluster))
) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Gene.Symbol+ liver_cluster , scales = "free", ncol = 4) +
  stat_cor(
    aes(label = ..r.label..),
    method = "pearson",
    label.x = -Inf,
    label.y = Inf,
    hjust = -0.1,
    vjust = 1.1,
    size = 3
  )+
  scale_color_manual(
    name = "Cluster",
    values = c("#0696c7", "#FFD700", "#ff8800", "#6e0280")
  ) +
  theme_classic2(base_size = 8) +
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  ) +
  labs(
    #title = "RNA–Protein Relationship for Top 6 Genes with Largest Correlation Changes",
    x = "Normalized mRNA expression",
    y = "Normalized protein expression"
  )
p

ggsave(
  filename = "./Plots/20-RNA_Protein_correlation_by_cluster_4gene_examples.pdf",
  plot = p,
  width = 6.6,
  height = 6.4,
  units = "in"
)


ggsave(
  filename = "./Plots/20-RNA_Protein_correlation_by_cluster_gene_examples_1.pdf",
  plot = p,
  width = 8,
  height = 5.5,
  units = "in"
)


goi <- c ("Col1a1", "Col1a2", "Col5a1")
plot_df <- Prot_RNA_unique %>%
  filter(
    Gene.Symbol %in% goi,
    !is.na(liver_cluster),
    !is.na(log2_Protein),
    !is.na(log2CPM_RNA)
  )

q <- ggplot(
  plot_df,
  aes(x = log2CPM_RNA, y = log2_Protein, color = factor(liver_cluster))
) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Gene.Symbol + liver_cluster, scales = "free", ncol = 4) +
  stat_cor(
    aes(label = ..r.label..),   # show only r
    method = "pearson",
    label.x.npc = "left",      # top-left corner
    label.y.npc = "top",
    size = 3
  ) +
  scale_color_manual(
    name = "Cluster",
    values = c("#0696c7", "#FFD700", "#ff8800", "#6e0280")
  ) +
  theme_minimal(base_size = 8) +
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  ) +
  labs(
    #title = "RNA–Protein Relationship for Top 6 Genes with Largest Correlation Changes",
    x = "RNA expression (log2CPM)",
    y = "Protein abundance (log2)"
  )
q

ggsave(
  filename = "./Plots/20-RNA_Protein_correlation_by_cluster_gene_examples_2.pdf",
  plot = q,
  width = 8,
  height = 5.5,
  units = "in"
)

goi <- c ("Ndufa2", "Atp6v1f", "Mrps7")
plot_df <- Prot_RNA_unique %>%
  filter(
    Gene.Symbol %in% goi,
    !is.na(liver_cluster),
    !is.na(log2_Protein),
    !is.na(log2CPM_RNA)
  )

r <- ggplot(
  plot_df,
  aes(x = log2CPM_RNA, y = log2_Protein, color = factor(liver_cluster))
) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Gene.Symbol + liver_cluster, scales = "free", ncol = 4) +
  stat_cor(
    aes(label = ..r.label..),   # show only r
    method = "pearson",
    label.x.npc = "left",      # top-left corner
    label.y.npc = "top",
    size = 3
  ) +
  scale_color_manual(
    name = "Cluster",
    values = c("#0696c7", "#FFD700", "#ff8800", "#6e0280")
  ) +
  theme_minimal(base_size = 8) +
  theme(
    strip.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) +
  labs(
    #title = "RNA–Protein Relationship for Top 6 Genes with Largest Correlation Changes",
    x = "RNA expression (log2CPM)",
    y = "Protein abundance (log2)"
  )
r

ggsave(
  filename = "./Plots/20-RNA_Protein_correlation_by_cluster_gene_examples_3.pdf",
  plot = r,
  width = 8,
  height = 5.5,
  units = "in"
)

corr_changes_filtered <- corr_changes %>%
  filter(Gene.Symbol %in% selected_genes)
corr_changes_filtered <- corr_changes_filtered %>%
  filter(corr_range > 0.5)

dfPlot2 <- corr_changes_filtered %>% dplyr::select(Gene.Symbol, "1", "2", 
                                                 "3",
                                                 "4")

hm.identity <- as.data.frame(dfPlot2[,1])
rownames(hm.identity) <- hm.identity[,1]
hm.identity <- hm.identity[,-1]

phenotype_names <- colnames(dfPlot2)
df2 <- dfPlot2[,-c(1)]
rownames(df2) <- dfPlot2$Gene.Symbol
#df <- t(df)
df2.matrix <- as.matrix(df2)
col<- colorRampPalette(c("red", "white", "blue"))(256)
col_fun <- colorRamp2(
  c(min(df2.matrix, na.rm = TRUE), 0, max(df2.matrix, na.rm = TRUE)),
  c("blue", "white", "red")
)
#col_ann <- HeatmapAnnotation(which ="col")
gene_heatmap <- Heatmap(df2.matrix, 
                       name = "RNA-prot correlation", #title of legend
                       cluster_columns = T, cluster_rows = T,
                       row_title = "Cluster", 
                       column_title = "Gene ID",
                       row_names_gp = gpar(fontsize = 6),
                       column_names_gp = gpar(fontsize = 6),
                       show_column_names = T, 
                       show_row_names = T,
                       #top_annotation = col_ann, 
                       #na_col = "black", 
                       #column_split = 9, 
                       #right_annotation = ha,
                       column_gap = unit(2, "mm"), row_gap = unit(2, "mm"), border = F,col = col_fun)
gene_heatmap
save_pdf(
  gene_heatmap,
  "./Plots/20-gene_heatmap.pdf",
  width = 25,
  height = 5,
  units = "cm")

mat <- df2.matrix
mat_scaled <- t(scale(t(mat)))
hc <- hclust(dist(mat_scaled), method = "ward.D2")
plot(hc)
clusters <- cutree(hc, k = 6) 

set.seed(123)
km <- kmeans(mat_scaled, centers = 6)   # try k from 4–10
clusters <- km$cluster

Prot_RNA$Strain<-factor(Prot_RNA$Strain,LEVELS_STRAINS)

mouse_id <- "HDP-012206"

Prot_RNA %>%
  filter(ID == mouse_id) %>%
  ggplot(aes(x = log2CPM, y = log2(Prot_quantif))) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    title = paste("RNA vs Protein for", mouse_id),
    x = "RNA Quantification",
    y = "Protein Quantification"
  )


colnames(corr_changes) <- c("gene_id", "cluster_1","cluster_2", "cluster_3", "cluster_4", "corr_range", "Gene.Symbol")
df_long <- corr_changes %>%
  pivot_longer(
    cols = starts_with("cluster_"),
    names_to = "cluster",
    values_to = "correlation"
  ) %>%
  mutate(
    cluster = gsub("cluster_", "", cluster),   # keep 1,2,3,4
    cluster = factor(cluster, levels = c("1","2","3","4"))
  )
p <- ggplot(df_long, aes(x = cluster, y = correlation, group = gene_id)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2) +
  labs(
    x = "Cluster",
    y = "RNA–protein correlation coefficient",
    title = "Correlation across clusters"
  ) +
  theme_minimal()


#Step1: compute direction changes per gene
threshold <- 0.2
df_pattern <- df_long %>%
  arrange(Gene.Symbol, cluster) %>%
  group_by(Gene.Symbol) %>%
  summarise(
    c1 = correlation[cluster == "1"],
    c2 = correlation[cluster == "2"],
    c3 = correlation[cluster == "3"],
    c4 = correlation[cluster == "4"]
  ) %>%
  mutate(
    d12 = case_when(
      c2 - c1 >  threshold ~  1,
      c2 - c1 < -threshold ~ -1,
      TRUE ~ 0
    ),
    d23 = case_when(
      c3 - c2 >  threshold ~  1,
      c3 - c2 < -threshold ~ -1,
      TRUE ~ 0
    ),
    d34 = case_when(
      c4 - c3 >  threshold ~  1,
      c4 - c3 < -threshold ~ -1,
      TRUE ~ 0
    )
  )



#Step 2 — assign pattern labels
df_pattern <- df_pattern %>%
  mutate(
    d12_lab = case_when(
      d12 == 1 ~ "U",
      d12 == -1 ~ "D",
      TRUE ~ "F"
    ),
    d23_lab = case_when(
      d23 == 1 ~ "U",
      d23 == -1 ~ "D",
      TRUE ~ "F"
    ),
    d34_lab = case_when(
      d34 == 1 ~ "U",
      d34 == -1 ~ "D",
      TRUE ~ "F"
    ),
    pattern = paste(d12_lab, d23_lab, d34_lab, sep = "-")
  )

#Step 3 — count genes per pattern
pattern_counts <- df_pattern %>%
  # clean pattern strings
  mutate(pattern = str_trim(as.character(pattern))) %>%
  # ensure one row per gene
  distinct(Gene.Symbol, pattern) %>%
  # count genes per pattern
  group_by(pattern) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  arrange(desc(n_genes))

c <- ggplot(pattern_counts, aes(x = reorder(pattern, -n_genes), y = n_genes)) +
  geom_col(fill = "#0696c7") +
  theme_minimal() +
  labs(x = "Pattern", y = "Number of genes", title = "Genes per correlation trajectory pattern") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
c
ggsave(useDingbats = FALSE, c,file="./Plots/20-genes_per_correlation_trajectory_pattern.pdf",
       width = 6, height = 2.8, units = c("in"))

#Step 4- merge back to long format
df_long <- df_long %>%
  left_join(df_pattern %>% dplyr::select(Gene.Symbol, pattern), by = "Gene.Symbol")

p <- ggplot(df_long, aes(x = cluster, y = correlation, group = gene_id, color = gene_id)) +
  geom_line(linewidth = 0.4, alpha = 0.4) +
  geom_point(size = 1) +
  facet_wrap(~ pattern, scales = "free_y") +   # split by pattern
  labs(
    x = "Cluster",
    y = "RNA–protein correlation coefficient",
    title = "Correlation across clusters by trajectory pattern"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "none"
  )
p
ggsave(useDingbats = FALSE, p,file="./Plots/20-Correlation_across_clusters_by_trajectory_pattern.pdf",
       width = 11, height = 7, units = c("in"))


genes_by_pattern <- df_long %>%
  group_by(pattern) %>%
  summarise(gene_list = list(gene_id), .groups = "drop") %>%
  deframe()   # convert to named list: names = pattern, values = vector of genes

ora_results_by_pattern <- lapply(genes_by_pattern, function(genes) {
  enricher(
    gene         = genes,
    TERM2GENE    = mouse_gene_sets_filtered,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
})

ora_results_by_pattern_symbol <- lapply(ora_results_by_pattern, function(x){
  if (is.null(x)) return(NULL)
  setReadable(x, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
})

#combine results for visualization

ora_pattern_df <- map2_df(
  ora_results_by_pattern_symbol,
  names(ora_results_by_pattern_symbol),
  ~{
    if (is.null(.x) || nrow(.x@result) == 0) return(NULL)
    .x@result %>%
      mutate(pattern = .y)
  }
)

#plot top pathways per pattern

top_n <- 5

df_top <- ora_pattern_df %>%
  group_by(pattern) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  ungroup()

ggplot(
  df_top,
  aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))
) +
  geom_point(aes(size = GeneRatio), alpha = 0.8) +
  facet_wrap(~ pattern, scales = "free") +
  labs(
    x = "-log10 adjusted p-value",
    y = "Pathway",
    title = "Top enriched pathways per pattern"
  ) +
  theme_minimal(base_size = 7) +                   # globally small
  theme(
    plot.title = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    strip.text = element_text(size = 7),            # facet labels
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing = unit(0.4, "lines"),             # tighter facets
    plot.margin = margin(2, 2, 2, 2, "pt")
  )  

##############transformation to accunt for different n of mice###################
##############linear model###################
library(dplyr)
library(broom)
library(purrr)

# Make sure cluster is factor and cluster 1 is reference
Prot_RNA_unique <- Prot_RNA_unique %>%
  mutate(liver_cluster = factor(liver_cluster)) %>%
  mutate(liver_cluster = relevel(liver_cluster, ref = "1"))

# Function to fit model per gene
fit_gene_model <- function(df) {
  
  # Skip genes with too few samples
  if(nrow(df) < 20) return(NULL)
  
  fit <- lm(log2_Protein ~ log2CPM_RNA * liver_cluster, data = df)
  
  tidy(fit)
}

# Fit models for all genes
model_results <- Prot_RNA_unique %>%
  group_by(gene_id) %>%
  group_modify(~ fit_gene_model(.x)) %>%
  ungroup()

#Reshape results
results_wide <- model_results %>%
  filter(term %in% c(
    "log2CPM_RNA",
    "log2CPM_RNA:liver_cluster2",
    "log2CPM_RNA:liver_cluster3",
    "log2CPM_RNA:liver_cluster4"
  )) %>%
  dplyr::select(gene_id, term, estimate, p.value) %>%
  tidyr::pivot_wider(
    names_from = term,
    values_from = c(estimate, p.value)
  )

#FDR correction
colnames(results_wide) <- gsub(":", "_", colnames(results_wide))
results_wide <- results_wide %>%
  mutate(
    FDR_2vs1 = p.adjust(p.value_log2CPM_RNA_liver_cluster2, method="BH"),
    FDR_3vs1 = p.adjust(p.value_log2CPM_RNA_liver_cluster3, method="BH"),
    FDR_4vs1 = p.adjust(p.value_log2CPM_RNA_liver_cluster4, method="BH")
  )
sig_cluster2 <- results_wide %>%
  filter(FDR_2vs1 < 0.05)

sig_cluster3 <- results_wide %>%
  filter(FDR_3vs1 < 0.05)

sig_cluster4 <- results_wide %>%
  filter(FDR_4vs1 < 0.05)

#make a volcano plot
library(dplyr)

volcano_lm <- bind_rows(
  
  results_wide %>%
    transmute(
      gene_id,
      comparison = "Cluster 2 vs 1",
      delta_slope = estimate_log2CPM_RNA_liver_cluster2,
      FDR = FDR_2vs1
    ),
  
  results_wide %>%
    transmute(
      gene_id,
      comparison = "Cluster 3 vs 1",
      delta_slope = estimate_log2CPM_RNA_liver_cluster3,
      FDR = FDR_3vs1
    ),
  
  results_wide %>%
    transmute(
      gene_id,
      comparison = "Cluster 4 vs 1",
      delta_slope = estimate_log2CPM_RNA_liver_cluster4,
      FDR = FDR_4vs1
    )
)

volcano_lm <- volcano_lm %>%
  mutate(
    neglog10FDR = -log10(FDR),
    significance = case_when(
      FDR < 0.05  
      #& abs(delta_slope) > 0.3
      ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

volcano_lm <- volcano_lm %>%
  left_join(ENS_to_gene_name %>% dplyr::select(gene_id, Gene.Symbol), by = "gene_id")

ggplot(volcano_lm, aes(x = delta_slope, y = neglog10FDR)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Significant" = "red",
    "Not significant" = "grey70"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_text_repel(
    data = volcano_lm %>% filter(significance != "Not significant"),
    aes(label = Gene.Symbol),
    size = 3,
    max.overlaps = 40,      # adjust to avoid clutter
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  facet_wrap(~ comparison, nrow = 1,scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Differential RNA–Protein Coupling (Linear Model)",
    x = expression(Delta~Slope~(Cluster~vs~1)),
    y = expression(-log[10](FDR))
  )





##############Fisher z transformation###################
#compute correlations per gene per cluster
library(dplyr)

cor_results <- Prot_RNA_unique %>%
  group_by(gene_id, liver_cluster) %>%
  summarise(
    n = n(),
    r = cor(log2CPM_RNA, log2_Protein, 
            use = "pairwise.complete.obs", 
            method = "pearson"),
    .groups = "drop"
  )

#Fisher transformation
cor_results <- cor_results %>%
  mutate(
    fisher_z = 0.5 * log((1 + r) / (1 - r))
  )

#reshape to wide format
library(tidyr)

cor_wide <- cor_results %>%
  pivot_wider(
    names_from = liver_cluster,
    values_from = c(r, fisher_z, n),
    names_glue = "{.value}_cl{liver_cluster}"
  )

#differential correlation test
diff_test <- function(zA, zB, nA, nB) {
  se <- sqrt(1/(nA - 3) + 1/(nB - 3))
  z <- (zA - zB) / se
  p <- 2 * pnorm(-abs(z))
  return(list(z = z, p = p))
}

cor_wide <- cor_wide %>%
  rowwise() %>%
  mutate(
    z_2vs1 = (fisher_z_cl2 - fisher_z_cl1) /
      sqrt(1/(n_cl2 - 3) + 1/(n_cl1 - 3)),
    p_2vs1 = 2 * pnorm(-abs(z_2vs1)),
    
    z_3vs1 = (fisher_z_cl3 - fisher_z_cl1) /
      sqrt(1/(n_cl3 - 3) + 1/(n_cl1 - 3)),
    p_3vs1 = 2 * pnorm(-abs(z_3vs1)),
    
    z_4vs1 = (fisher_z_cl4 - fisher_z_cl1) /
      sqrt(1/(n_cl4 - 3) + 1/(n_cl1 - 3)),
    p_4vs1 = 2 * pnorm(-abs(z_4vs1))
  ) %>%
  ungroup()

#FDR correction
cor_wide <- cor_wide %>%
  mutate(
    FDR_2vs1 = p.adjust(p_2vs1, method="BH"),
    FDR_3vs1 = p.adjust(p_3vs1, method="BH"),
    FDR_4vs1 = p.adjust(p_4vs1, method="BH")
  )

sig_cluster2 <- cor_wide %>%
  filter(FDR_2vs1 < 0.05)

sig_cluster3 <- cor_wide %>%
  filter(FDR_3vs1 < 0.05)

sig_cluster4 <- cor_wide %>%
  filter(FDR_4vs1 < 0.05)


#make a volcano plot
volc2 <- cor_wide %>%
  mutate(
    delta_r = r_cl2 - r_cl1,
    neglog10FDR = -log10(FDR_2vs1),
    significance = case_when(
      FDR_2vs1 < 0.05 & abs(delta_r) > 0.3 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )


volcano_df <- bind_rows(
  
  cor_wide %>%
    transmute(
      gene_id,
      comparison = "Cluster 2 vs 1",
      delta_r = r_cl2 - r_cl1,
      FDR = FDR_2vs1
    ),
  
  cor_wide %>%
    transmute(
      gene_id,
      comparison = "Cluster 3 vs 1",
      delta_r = r_cl3 - r_cl1,
      FDR = FDR_3vs1
    ),
  
  cor_wide %>%
    transmute(
      gene_id,
      comparison = "Cluster 4 vs 1",
      delta_r = r_cl4 - r_cl1,
      FDR = FDR_4vs1
    )
)

volcano_df <- volcano_df %>%
  mutate(
    neglog10FDR = -log10(FDR),
    significance = case_when(
      FDR < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

volcano_df <- volcano_df %>%
  mutate(direction = case_when(
    FDR < 0.05 & delta_r > 0  ~ "Increased",
    FDR < 0.05 & delta_r < 0 ~ "Decreased",
    TRUE ~ "Not significant"
  ))



ggplot(volc2, aes(x = delta_r, y = neglog10FDR)) +
  geom_point(aes(color = significance), alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red",
                                "Not significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Differential RNA–Protein Correlation",
    subtitle = "Cluster 2 vs Cluster 1",
    x = expression(Delta~Correlation~(r[2] - r[1])),
    y = expression(-log[10](FDR))
  )

volcano_df <- volcano_df %>%
  left_join(ENS_to_gene_name %>% dplyr::select(gene_id, Gene.Symbol), by = "gene_id")

summary_counts <- volcano_df %>%
  filter(FDR < 0.05) %>%             # keep only significant genes
  mutate(direction = ifelse(delta_r > 0, "Positive", "Negative")) %>%
  group_by(comparison, direction) %>%
  summarise(n_genes = n(), .groups = "drop")
summary_counts

ggplot(volcano_df, aes(x = delta_r, y = neglog10FDR)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Significant" = "red",
    "Not significant" = "grey70"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  facet_wrap(~ comparison, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Differential RNA–Protein Correlation",
    x = expression(Delta~Correlation~(r[cluster] - r[1])),
    y = expression(-log[10](FDR))
  )

volcano_df <- volcano_df %>%
  mutate(Label = case_when(
    comparison == "Cluster 2 vs 1" & Gene.Symbol %in% genes_df_top_pos2 ~ "yes",
    comparison == "Cluster 3 vs 1" & Gene.Symbol %in% genes_df_top_pos ~ "yes",
    TRUE ~ "no"
  ))


# Volcano plot with custom colors
library(ggrepel)
v <- ggplot(volcano_df, aes(x = delta_r, y = neglog10FDR)) +
  geom_point(aes(color = direction), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Increased" = "#ff0000",   # red
    "Decreased" = "#0000ff",   # blue
    "Not significant" = "grey70"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_text_repel(
    data = volcano_df %>% filter(Label == "yes"),
    aes(label = Gene.Symbol),
    size = 1,
    max.overlaps = Inf,      # ensure all "yes" genes are attempted
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  facet_wrap(~ comparison, nrow = 1) +
  theme_classic(base_size = 8) +
  theme(
    strip.text = element_text(size = 8),        # facet labels
    axis.title = element_text(size = 8),        # x/y axis titles
    axis.text = element_text(size = 8),         # x/y axis tick labels
    legend.text = element_text(size = 8),       # legend labels
    legend.title = element_text(size = 8),      # legend title
    plot.title = element_text(size = 8)         # plot title
  ) +
  labs(
    title = "Differential RNA–Protein Correlation",
    x = expression(Delta~Correlation~(r[cluster] - r[1])),
    y = expression(-log[10](FDR))
  )
v 
ggsave(
  filename = "./Plots/20-RNA_Protein_delta_correlation_by_cluster_volcano_big.pdf",
  plot = v,
  width = 20,
  height = 8,
  units = "in"
)

ggsave(
  filename = "./Plots/20-RNA_Protein_delta_correlation_by_cluster_volcano_small.pdf",
  plot = v,
  width = 5,
  height = 2,
  units = "in"
)

#do ORA on significant genes per comparison
# Filter significant genes per comparison
sig_genes_list <- volcano_df %>%
  filter(FDR < 0.05) %>%  # same thresholds as volcano
  group_by(comparison) %>%
  summarise(genes = list(gene_id)) %>%
  deframe()  # produces a named list: comparison -> vector of genes

sig_genes_list_pos <- volcano_df %>%
  filter(FDR < 0.05, direction == "Increased") %>%  # same thresholds as volcano
  group_by(comparison) %>%
  summarise(genes = list(gene_id)) %>%
  deframe()
overlap_1_2 <- intersect(sig_genes_list_pos[[1]], sig_genes_list_pos[[2]])

sig_genes_list_neg <- volcano_df %>%
  filter(FDR < 0.05, direction == "Decreased") %>%  # same thresholds as volcano
  group_by(comparison) %>%
  summarise(genes = list(gene_id)) %>%
  deframe()
overlap_1_2 <- intersect(sig_genes_list_neg[[1]], sig_genes_list_neg[[2]])
####################

genes_by_comparison <- volcano_df %>%
  filter(FDR < 0.05) %>%   # adjust threshold if needed
  group_by(comparison) %>%
  summarise(genes = list(unique(gene_id)), .groups = "drop")
gene_list <- setNames(genes_by_comparison$genes,
                      genes_by_comparison$comparison)

ora_results <- lapply(gene_list, function(genes) {
  
  enricher(
    gene          = genes,
    TERM2GENE     = mouse_gene_sets_filtered,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
})

ora_results_symbol <- lapply(ora_results, function(x){
  if (is.null(x)) return(NULL)
  setReadable(x, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
})

genes_by_comparison_pos <- volcano_df %>%
  filter(FDR < 0.05, direction == "Increased") %>%   # adjust threshold if needed
  group_by(comparison) %>%
  summarise(genes = list(unique(gene_id)), .groups = "drop")
gene_list_pos <- setNames(genes_by_comparison_pos$genes,
                      genes_by_comparison_pos$comparison)

ora_results_pos <- lapply(gene_list_pos, function(genes) {
  
  enricher(
    gene          = genes,
    TERM2GENE     = mouse_gene_sets_filtered,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
})

ora_results_symbol_pos <- lapply(ora_results_pos, function(x){
  if (is.null(x)) return(NULL)
  setReadable(x, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
})

genes_by_comparison_neg <- volcano_df %>%
  filter(FDR < 0.05, direction == "Decreased") %>%   # adjust threshold if needed
  group_by(comparison) %>%
  summarise(genes = list(unique(gene_id)), .groups = "drop")
gene_list_neg <- setNames(genes_by_comparison_neg$genes,
                          genes_by_comparison_neg$comparison)

ora_results_neg <- lapply(gene_list_neg, function(genes) {
  
  enricher(
    gene          = genes,
    TERM2GENE     = mouse_gene_sets_filtered,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
})

ora_results_symbol_neg <- lapply(ora_results_neg, function(x){
  if (is.null(x)) return(NULL)
  setReadable(x, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
})

# how many top pathways per category to show
top_n <- 20

# extract tables from enrichment results
df_all <- map2_df(
  ora_results_symbol,
  names(ora_results_symbol),
  ~{
    if (is.null(.x) || nrow(.x@result) == 0) return(NULL)
    
    .x@result %>%
      filter(p.adjust < 0.05) %>%          # keep only significant pathways
      mutate(
        category = .y,
        neg_log10_padj = -log10(p.adjust)
      )
  }
)
df_all <- df_all %>%
  filter(category != "Cluster 4 vs 1")
df_all <- df_all %>%
  mutate(
    Description = Description %>%
      str_replace_all("NUCLEAR_TRANSCRIBED_", "") %>%
      str_replace_all("CONTAINING_COMPOUND_", "") %>%
      str_replace_all("NON_MEMBRANE_BOUNDED_", "") %>%
      str_replace_all("REGULATION_OF_", "") %>%
      str_replace_all("ESTABLISHMENT_OF_", "") %>%
      str_replace_all("ENERGY_DERIVATION_BY_", "") %>%
      str_replace_all("GENERATION_OF_", "") %>%
      str_replace_all("_BIOGENESIS", "") %>%
      str_replace_all("_PATHWAY", "") %>%
      str_replace_all("GOBP_", "") %>%
      str_replace_all("KEGG_", "") %>%
      str_replace_all("HALLMARK_", "") %>%
      str_replace_all("METABOLIC_PROCESS", "METABOLISM") %>%
      str_replace_all("CATABOLIC_PROCESS", "CATABOLISM") %>%
      str_replace_all("BIOSYNTHETIC_PROCESS", "BIOSYNTHESIS") %>%
      str_replace_all("SIGNAL_TRANSDUCTION", "SIGNALING") %>%
      str_replace_all("_", " ")
  )

df_all <- df_all %>%
  group_by(category, Description) %>%
  slice_max(neg_log10_padj, n = 1, with_ties = FALSE) %>%
  ungroup()

# pick top pathways in each category by adjusted p-value
df_top <- df_all %>%
  group_by(category) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  ungroup()


# reorder descriptions within categories for nicer plotting
df_top <- df_top %>%
  mutate(
    Description = factor(Description, levels = unique(Description[order(category, neg_log10_padj)]))
  )

library(tidytext)
df_top <- df_top %>%
  mutate(
    Description_reordered = reorder_within(Description, neg_log10_padj, category)
  )
df_top <- df_top %>%
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]))
  )


g<- ggplot(df_top, aes(
  x = neg_log10_padj,
  y = Description_reordered,
  size = GeneRatio_num
)) +
  geom_point() +
  facet_wrap(~ category, scales = "free") +
  scale_y_reordered() +
  labs(
    # title = "Top enriched pathways per category",
    x = "-log10 adjusted p-value",
    y = "Pathway",
    size = "GeneRatio"
  ) +
  theme_bw()+
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )+
  scale_size_continuous(range = c(1, 4))
g
ggsave("./Plots/20-Top_enriched_pathways_per_category.pdf", plot = g, width = 11.5, height = 3.6)

################ extract tables from enrichment results - positive corr only#############
df_all_pos <- map2_df(
  ora_results_symbol_pos,
  names(ora_results_symbol_pos),
  ~{
    if (is.null(.x) || nrow(.x@result) == 0) return(NULL)
    
    .x@result %>%
      filter(p.adjust < 0.05) %>%          # keep only significant pathways
      mutate(
        category = .y,
        neg_log10_padj = -log10(p.adjust)
      )
  }
)
df_all_pos <- df_all_pos %>%
  filter(category != "Cluster 4 vs 1")
df_all_pos <- df_all_pos %>%
  mutate(
    Description = Description %>%
      str_replace_all("NUCLEAR_TRANSCRIBED_", "") %>%
      str_replace_all("CONTAINING_COMPOUND_", "") %>%
      str_replace_all("NON_MEMBRANE_BOUNDED_", "") %>%
      str_replace_all("REGULATION_OF_", "") %>%
      str_replace_all("ESTABLISHMENT_OF_", "") %>%
      str_replace_all("ENERGY_DERIVATION_BY_", "") %>%
      str_replace_all("GENERATION_OF_", "") %>%
      str_replace_all("_BIOGENESIS", "") %>%
      str_replace_all("_PATHWAY", "") %>%
      str_replace_all("GOBP_", "") %>%
      str_replace_all("KEGG_", "") %>%
      str_replace_all("HALLMARK_", "") %>%
      str_replace_all("METABOLIC_PROCESS", "METABOLISM") %>%
      str_replace_all("CATABOLIC_PROCESS", "CATABOLISM") %>%
      str_replace_all("BIOSYNTHETIC_PROCESS", "BIOSYNTHESIS") %>%
      str_replace_all("SIGNAL_TRANSDUCTION", "SIGNALING") %>%
      str_replace_all("_", " ")
  )

df_all_pos <- df_all_pos %>%
  group_by(category, Description) %>%
  slice_max(neg_log10_padj, n = 1, with_ties = FALSE) %>%
  ungroup()

# pick top pathways in each category by adjusted p-value
df_top_pos <- df_all_pos %>%
  group_by(category) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  ungroup()


# reorder descriptions within categories for nicer plotting
df_top_pos <- df_top_pos %>%
  mutate(
    Description = factor(Description, levels = unique(Description[order(category, neg_log10_padj)]))
  )

library(tidytext)
df_top_pos <- df_top_pos %>%
  mutate(
    Description_reordered = reorder_within(Description, neg_log10_padj, category)
  )
df_top_pos <- df_top_pos %>%
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]))
  )

filtered <- df_top_pos %>%
  filter(category == "Cluster 3 vs 1")
filtered2 <- df_top_pos %>%
  filter(category == "Cluster 2 vs 1")
genes_df_top_pos <- filtered %>%
  pull(geneID) %>%
  strsplit("/") %>%
  unlist() %>%
  unique()
genes_df_top_pos2 <- filtered2 %>%
  pull(geneID) %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

h<- ggplot(df_top_pos, aes(
  x = neg_log10_padj,
  y = Description_reordered,
  size = GeneRatio_num
)) +
  geom_point() +
  facet_wrap(~ category, scales = "free") +
  scale_y_reordered() +
  labs(
    # title = "Top enriched pathways per category",
    x = "-log10 adjusted p-value",
    y = "Pathway",
    size = "GeneRatio"
  ) +
  theme_bw()+
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )+
  scale_size_continuous(range = c(1, 4))
h
ggsave("./Plots/20-Top_enriched_pathways_per_category_pos.pdf", plot = h, width = 9.8, height = 3.5)

################ extract tables from enrichment results - negative corr only#############
df_all_neg <- map2_df(
  ora_results_symbol_neg,
  names(ora_results_symbol_neg),
  ~{
    if (is.null(.x) || nrow(.x@result) == 0) return(NULL)
    
    .x@result %>%
      filter(p.adjust < 0.05) %>%          # keep only significant pathways
      mutate(
        category = .y,
        neg_log10_padj = -log10(p.adjust)
      )
  }
)
df_all_neg <- df_all_neg %>%
  filter(category != "Cluster 4 vs 1")
df_all_neg <- df_all_neg %>%
  mutate(
    Description = Description %>%
      str_replace_all("NUCLEAR_TRANSCRIBED_", "") %>%
      str_replace_all("CONTAINING_COMPOUND_", "") %>%
      str_replace_all("NON_MEMBRANE_BOUNDED_", "") %>%
      str_replace_all("REGULATION_OF_", "") %>%
      str_replace_all("ESTABLISHMENT_OF_", "") %>%
      str_replace_all("ENERGY_DERIVATION_BY_", "") %>%
      str_replace_all("GENERATION_OF_", "") %>%
      str_replace_all("_BIOGENESIS", "") %>%
      str_replace_all("_PATHWAY", "") %>%
      str_replace_all("GOBP_", "") %>%
      str_replace_all("KEGG_", "") %>%
      str_replace_all("HALLMARK_", "") %>%
      str_replace_all("METABOLIC_PROCESS", "METABOLISM") %>%
      str_replace_all("CATABOLIC_PROCESS", "CATABOLISM") %>%
      str_replace_all("BIOSYNTHETIC_PROCESS", "BIOSYNTHESIS") %>%
      str_replace_all("SIGNAL_TRANSDUCTION", "SIGNALING") %>%
      str_replace_all("_", " ")
  )

df_all_neg <- df_all_neg %>%
  group_by(category, Description) %>%
  slice_max(neg_log10_padj, n = 1, with_ties = FALSE) %>%
  ungroup()

# pick top pathways in each category by adjusted p-value
df_top_neg <- df_all_neg %>%
  group_by(category) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  ungroup()


# reorder descriptions within categories for nicer plotting
df_top_neg <- df_top_neg %>%
  mutate(
    Description = factor(Description, levels = unique(Description[order(category, neg_log10_padj)]))
  )

library(tidytext)
df_top_neg <- df_top_neg %>%
  mutate(
    Description_reordered = reorder_within(Description, neg_log10_padj, category)
  )
df_top_neg <- df_top_neg %>%
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]))
  )


i<- ggplot(df_top_neg, aes(
  x = neg_log10_padj,
  y = Description_reordered,
  size = GeneRatio_num
)) +
  geom_point() +
  facet_wrap(~ category, scales = "free") +
  scale_y_reordered() +
  labs(
    # title = "Top enriched pathways per category",
    x = "-log10 adjusted p-value",
    y = "Pathway",
    size = "GeneRatio"
  ) +
  theme_bw()+
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )+
  scale_size_continuous(range = c(1, 4))
i
ggsave("./Plots/20-Top_enriched_pathways_per_category_neg.pdf", plot = i, width = 9.8, height = 3.5)
####################
ora_go_list <- lapply(sig_genes_list, function(genes) {
  enrichGO(
    gene         = genes,
    OrgDb        = org.Mm.eg.db,
    keyType      = "ENSEMBL",
    ont          = "BP",
    pAdjustMethod= "BH",
    qvalueCutoff = 0.05,
    readable     = TRUE
  )
})

ora_kegg_list <- lapply(sig_genes_list, function(genes) {
  enrichKEGG(
    gene         = genes,
    organism     = "mmu",
    keyType      = "kegg",
    pAdjustMethod= "BH",
    qvalueCutoff = 0.05
  )
})

ora_go_df <- bind_rows(lapply(names(ora_go_list), function(comp) {
  df <- as.data.frame(ora_go_list[[comp]])
  df$comparison <- comp
  df
}))

top_go <- ora_go_df %>%
  group_by(comparison) %>%
  slice_max(order_by = -p.adjust, n = 20)



# Prepare plotting data
plot_df <- ora_go_df %>%
  filter(comparison != "Cluster 4 vs 1") %>%  # remove sparse comparison
  filter(p.adjust < 0.05)  # keep only significant pathways

# Keep top 20 pathways per comparison based on significance
plot_df <- plot_df %>%
  group_by(comparison) %>%
  slice_max(order_by = -log10(p.adjust), n = 20) %>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust))) %>%
  ungroup()

# Horizontal dotplot
ggplot(plot_df, aes(x = -log10(p.adjust), y = Description, size = Count, color = Count)) +
  geom_point() +
  facet_wrap(~ comparison, scales = "free_y") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 20 GO BP pathways for significant differential correlation genes",
    x = expression(-log[10](FDR)),
    y = "GO term",
    color = "Gene count",
    size = "Gene count"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

#do GSEA with r value
library(fgsea)
library(msigdbr)

pathways <- msig %>%
  dplyr::filter(!is.na(ensembl_gene)) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

library(msigdbr)
library(fgsea)
library(dplyr)

# Hallmark
hallmark <- msigdbr(species = "Mus musculus", category = "H")

# Reactome
reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")

# KEGG
kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")

# GO BP
gobp <- msigdbr(
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:BP"
)


all_pathways <- rbind(hallmark, reactome, kegg, gobp)
pathways <- all_pathways%>%
  dplyr::filter(!is.na(ensembl_gene)) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

#create ranked vector per cluster
make_rank_vector <- function(df, r_column) {
  
  tmp <- df %>%
    dplyr::select(gene_id, !!sym(r_column)) %>%
    filter(!is.na(.data[[r_column]])) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  ranks <- tmp[[r_column]]
  names(ranks) <- tmp$gene_id
  
  sort(ranks, decreasing = TRUE)
}

rank_cl1 <- make_rank_vector(cor_wide, "r_cl1")
rank_cl2 <- make_rank_vector(cor_wide, "r_cl2")
rank_cl3 <- make_rank_vector(cor_wide, "r_cl3")
rank_cl4 <- make_rank_vector(cor_wide, "r_cl4")


gsea_cl1 <- fgsea(
  pathways = pathways,
  stats = rank_cl1,
  nperm = 10000
)

gsea_cl2 <- fgsea(
  pathways = pathways,
  stats = rank_cl2,
  nperm = 10000
)

gsea_cl3 <- fgsea(
  pathways = pathways,
  stats = rank_cl3,
  nperm = 10000
)

gsea_cl4 <- fgsea(
  pathways = pathways,
  stats = rank_cl4,
  nperm = 10000
)







