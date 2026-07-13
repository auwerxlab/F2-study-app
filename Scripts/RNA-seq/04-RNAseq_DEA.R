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

################## DEA ###################
# male vs female
dsx <- model.matrix(~Sex, data = liver_dge$samples %>% filter(Generation == "F2") )
vm_sex <- voom(liver_dge$counts_filtered[,(liver_dge$samples$Generation=="F2")], design = dsx)
fit_sex <- lmFit(vm_sex, design = dsx)
fit_sex <- eBayes(fit_sex, robust=TRUE)
summary(decideTests(fit_sex))
tt_male_female <- topTable(fit_sex, n = Inf) %>% rownames_to_column(var = "gene_id") %>% left_join(., liver_dge$geneConversionTable, by = "gene_id") %>% 
  mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% dplyr::arrange(-signed_log_p)
## do gsea
gs <- gseGO(geneList = setNames(tt_male_female$signed_log_p, tt_male_female$gene_id), ont = "BP", OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",pvalueCutoff = 1, nPermSimple = 100000)
p <- gs@result %>% filter(p.adjust < 0.7) %>% ggplot(aes(x = NES, y = -log10(p.adjust), desc = Description)) + geom_point()
ggplotly(p)
gs_toplot <- gs@result %>% split(., sign(gs@result$NES)) %>% lapply(., head, n = 10) %>% do.call(rbind, .) %>% arrange(NES)
gs_toplot$nleading <- sapply(strsplit(gs_toplot$core_enrichment, "/"), length)
gs_toplot$`Gene Ratio` <- gs_toplot$nleading/gs_toplot$setSize
gs_toplot$Description <- gsub('(.{1,30})(\\s|$)', '\\1\n', gs_toplot$Description)
gs_toplot$Description <- gsub('\n$)', '', gs_toplot$Description)
gs_toplot$Description <- factor(gs_toplot$Description, levels = gs_toplot$Description)
gsea_male_vs_female <- gs_toplot %>%  ggplot(aes(x = Description, y = NES, color = -log10(p.adjust), size=`Gene Ratio`)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  coord_flip() +
  scale_color_viridis_c(limits = c(0,NA), direction = -1, option = "magma") +
  scale_size_continuous(range = c(2,10)) +
  theme_cowplot() +
  ylab("Normalized Enrichment Score (NES)") +
  xlab("GO category") +
  ggtitle("Top enriched biological processes\nin sexually dimorphic genes")
gsea_male_vs_female
ggsave(useDingbats = F, filename = "./Plots/17-liver_dge_male_vs_female_gsea.pdf", width = 7, height = 7)
v_male_vs_female <- ggplot(tt_male_female, aes(x = logFC, y = -log10(adj.P.Val), color = logFC, label=gene_name)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  scale_color_gradient2(low = sexColors["f"], high = sexColors["m"], midpoint = 0, mid = "white") +
  theme_cowplot() +
  ggtitle("Differential expression by sex\n(male - female)")
v_male_vs_female
ggplotly(v_male_vs_female)
plot_grid(v_male_vs_female, gsea_male_vs_female, nrow = 1)
ggsave(useDingbats = F, filename = "./Plots/17-liver_dge_male_vs_female_vp_gsea.pdf", width = 15, height = 7)


#final fat percentage
dft <- model.matrix(~Fat_perc_22, data = liver_dge$samples %>% filter(Generation == "F2") %>% filter(Sex == "f"))
vm_fat <- voom(liver_dge$counts_filtered[,(liver_dge$samples$Generation=="F2") & (liver_dge$samples$Sex=="f")], design = dft)
fit_fat <- lmFit(vm_fat, design = dft)
fit_fat <- eBayes(fit_fat, robust=TRUE)
summary(decideTests(fit_fat))
topTable(fit_fat, coef = 2, n = 30)

########### fibrosis ##############
dfib <- model.matrix(~log(Fibrosis_perc_24) + Sex + log(TIMP1_24) + Liver_g_RelBWSac_24 + sum_all_vacuoles_percentage_24 + Fat_perc_22, data = liver_dge$samples %>% 
                       filter(Generation == "F2") %>% 
                       filter(!is.na(Fibrosis_perc_24)) %>% 
                       filter(!is.na(TIMP1_24)) %>% 
                       filter(!is.na(Liver_g_RelBWSac_24))  %>% 
                       filter(!is.na(sum_all_vacuoles_percentage_24))  %>%
                       filter(!is.na(Fat_perc_22)) )
# dtimp_voom <- model.matrix(~ Sex, data = liver_dge$samples %>% filter(Generation == "F2") %>% filter(!is.na(TIMP.1)) %>% filter(!is.na(Liver_g_RelBWSac_24)))
vm_fib <- voom(liver_dge$counts_filtered[,(liver_dge$samples$Generation=="F2") & 
                                           (!is.na(liver_dge$samples$Fibrosis_perc_24))& 
                                           (!is.na(liver_dge$samples$TIMP1_24)) & 
                                           (!is.na(liver_dge$samples$Liver_g_RelBWSac_24)) & 
                                           (!is.na(liver_dge$samples$sum_all_vacuoles_percentage_24)) & 
                                           (!is.na(liver_dge$samples$Fat_perc_22)) ], design = dfib)
fit_fib <- lmFit(vm_fib, design = dfib)
fit_fib <- eBayes(fit_fib, robust=TRUE)
summary(decideTests(fit_fib))

dt <- summary(decideTests(fit_fib))
dt <- as.data.table(dt)
colnames(dt) <- c("Direction", "Phenotype", "Count")
dt_filtered <- dt %>%
  filter(Direction %in% c("Down", "Up"))
library(data.table)
dt_filtered[,count_direction:=ifelse(Direction=="Up", Count,-Count)]
dt_filtered <- dt_filtered %>%
  filter(Phenotype != "(Intercept)")

dt_filtered <- dt_filtered %>%
  mutate(
    Phenotype = recode(
      Phenotype,
      "log(Fibrosis_perc_24)" = "Fibrosis(%)",
      "Sexm" = "Sex(Male)",
      "log(TIMP1_24)" = "plasma TIMP1",
      "Liver_g_RelBWSac_24"= "Liver(g/g BW)",
      "sum_all_vacuoles_percentage_24"="Steatosis AV(%)",
      "Fat_perc_22"= "Fat 22wk(%)"
    )
  )

dt_filtered$Phenotype <- factor(dt_filtered$Phenotype, levels = c(
  "Sex(Male)",
  "Fat 22wk(%)",
  "Liver(g/g BW)",
  "plasma TIMP1",
  "Steatosis AV(%)",
  "Fibrosis(%)"
))
gp1 <- ggplot(dt_filtered, aes(x = Phenotype, y = count_direction, fill = Direction))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Count, group = Direction),
            position = position_stack(vjust = 1.05),
            size = 2)+
  geom_hline(yintercept = 0,  size = 0.3)+
  #facet_wrap(~DurationName, scales = "free_y", nrow = 1)+
  scale_fill_manual(values = c("Down"="blue","Up"="red"))+
  xlab("")+
  ylab("DEGs(n)")+
  theme(#
    #axis.text.x = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    legend.text=element_text(size=8),
    legend.title = element_text(size=8),
    axis.line = element_line(colour = 'black', size = 0.3),
    axis.ticks.y = element_line(size=0.3),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.1, "cm"))
ggsave("17-DEA_count.pdf", plot= gp1, path = "./Plots", width = 2.7, height = 2.1)




tt <- topTable(fit_fib, coef = "log(Fibrosis_perc_24)", n = Inf) %>% rownames_to_column(var = "gene_id") %>% left_join(., liver_dge$geneConversionTable, by = "gene_id") %>% 
  mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% dplyr::arrange(-signed_log_p)
head(tt, n = 20)
# tt <- tt %>% mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% dplyr::arrange(-signed_log_p)
ttsex <- topTable(fit_fib, coef = "Sexm", n = Inf) %>% rownames_to_column(var = "gene_id") %>% left_join(., liver_dge$geneConversionTable, by = "gene_id")  %>% 
  mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% dplyr::arrange(-signed_log_p)
ttliverw <- topTable(fit_fib, coef = "Liver_g_RelBWSac_24", n = Inf) %>% rownames_to_column(var = "gene_id") %>% left_join(., liver_dge$geneConversionTable, by = "gene_id")  %>% 
  mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% dplyr::arrange(-signed_log_p)
ttfat <- topTable(fit_fib, coef = "Fat_perc_22", n = Inf) %>% rownames_to_column(var = "gene_id") %>% left_join(., liver_dge$geneConversionTable, by = "gene_id")  %>% 
  mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% dplyr::arrange(-signed_log_p)
tttimp <- topTable(fit_fib, coef = "log(TIMP1_24)", n = Inf) %>% rownames_to_column(var = "gene_id") %>% left_join(., liver_dge$geneConversionTable, by = "gene_id")  %>% 
  mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% dplyr::arrange(-signed_log_p)
ttlf <- topTable(fit_fib, coef = "sum_all_vacuoles_percentage_24", n = Inf) %>% rownames_to_column(var = "gene_id") %>% left_join(., liver_dge$geneConversionTable, by = "gene_id")  %>% 
  mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% dplyr::arrange(-signed_log_p)


selected_genes <- ttlf[ttlf$gene_name %in% colnames(cor_mat),]
#############try correlation between genes and phenotypes###############
expr <- liver_dge$y$E
pheno <- liver_dge$samples

phenos_of_interest <- c("Liver_g_24","Liver_g_RelBWSac_24","Fibrosis_perc_24",
                        "GDF15_24", "TIMP1_24", "sum_all_vacuoles_percentage_24",
                        "ALAT_24")
genes_of_interest <- ttlf[1:30, 1]
expr_sel <- expr[genes_of_interest, , drop = FALSE]
pheno_sel <- pheno[, phenos_of_interest, drop = FALSE]

# Log-transform specific phenotypes
to_log <- c("Fibrosis_perc_24",
            "GDF15_24", "TIMP1_24","ALAT_24")
pheno_sel[to_log] <- lapply(pheno_sel[to_log], function(x) log1p(x))

# Scale all phenotypes
pheno_scaled <- as.data.frame(scale(pheno_sel))


# Correlation matrix
cor_mat <- sapply(phenos_of_interest, function(p) {
  sapply(genes_of_interest, function(g) {
    cor(expr_sel[g, ], pheno_sel[[p]],
        use = "pairwise.complete.obs",
        method = "spearman")
  })
})

# P-value matrix
pval_mat <- sapply(phenos_of_interest, function(p) {
  sapply(genes_of_interest, function(g) {
    cor.test(expr_sel[g, ], pheno_sel[[p]],
             use = "pairwise.complete.obs",
             method = "spearman")$p.value
  })
})

# Transpose to match pheatmap structure
cor_mat  <- t(cor_mat)
pval_mat <- t(pval_mat)

pval_stars <- ifelse(pval_mat < 0.001, "***",
                     ifelse(pval_mat < 0.01,  "**",
                            ifelse(pval_mat < 0.05,  "*", "")))



colnames(cor_mat) <- liver_dge$geneConversionTable$gene_name[
  match(colnames(cor_mat),
        liver_dge$geneConversionTable$gene_id)
]

colnames(pval_stars) <- colnames(cor_mat)
rownames(pval_stars) <- rownames(cor_mat)

cor_mat[pval_mat >= 0.05] <- NA

library(pheatmap)

pheatmap(
  cor_mat,
  cluster_rows = F,
  cluster_cols = T,
  display_numbers = TRUE,
  main = "Correlation: Gene Expression vs Phenotypes",
  color = colorRampPalette(c("white", "pink","red"))(50)
)

num_mat <- matrix(
  paste0(round(cor_mat, 2), pval_stars),
  nrow = nrow(cor_mat),
  dimnames = dimnames(cor_mat)
)


pheatmap(
  cor_mat,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  display_numbers = num_mat,
  number_color = "black",
  main = "Correlation: Gene Expression vs Phenotypes",
  color = colorRampPalette(c("white", "pink", "red"))(50)
)



cor_mat_steatosis <- cor_mat[6,]

cor_mat_steatosis_df <- data.frame(
  gene_name = names(cor_mat_steatosis),
  correlation = as.numeric(cor_mat_steatosis),
  row.names = NULL
)
ttlf2 <- merge(ttlf, cor_mat_steatosis_df, by = "gene_name", all.x = TRUE)
ttlf2 <- ttlf2[ttlf2$gene_name %in% names(cor_mat_steatosis),]
ttlf2$gene_name <- factor(ttlf2$gene_name,
                          levels = ttlf2$gene_name[order(ttlf2$correlation, decreasing = TRUE)])

top_steatosis_genes <- ggplot(ttlf2, aes(
  x = gene_name,
  y = correlation,
  size = logFC,
  fill = signed_log_p
)) +
  geom_point(
    shape = 21, colour = "black", stroke = 0.5
  ) +
  scale_size_continuous(name = "LogFC") +
  scale_fill_gradientn(
    colors = c( "#FFE6EC", "#FF8AAE", "#D60000"),
    name = "Signed logP"
  )+ 
  
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=8),
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.7, "lines"),   # smaller symbol boxes
    legend.spacing.y = unit(0.1, "lines")
    #plot.title = element_text(size = 8)
    #panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Gene",
    y = "Steatosis correlation coefficient"
    #,
    #title = "Gene Expression with Signed p and Correlation"
  )
ggsave(top_steatosis_genes,file = "./Plots/17-top_steatosis_genes.pdf", width =5.8, height = 2, units = c("in"), useDingbats = F)

#sig_labels <- ifelse(pval_mat < 0.001, "***",
#                     ifelse(pval_mat < 0.01, "**",
#                            ifelse(pval_mat < 0.05, "*", "")))

#pheatmap(
#  cor_mat,
#  cluster_rows = FALSE,
#  cluster_cols = FALSE,
#  display_numbers = sig_labels,
#  main = "Correlation Heatmap (with significance)",
#  color = colorRampPalette(c("blue", "white", "red"))(50)
#)


dtest<- decideTests(fit_fib, lfc = 0)[,c("log(TIMP1_24)",
                                         #"Liver_g_RelBWSac_24",
                                         "sum_all_vacuoles_percentage_24")]
vennDiagram(dtest, include=c("up","down") )

up <- lapply(colnames(dtest), function(x){
  dtest[,x] %>% as.data.frame() %>% filter(.[[1]] == 1) %>% rownames()
})
names(up) <- paste0(colnames(dtest), " - up")
down <- lapply(colnames(dtest), function(x){
  dtest[,x] %>% as.data.frame() %>% filter(.[[1]] == -1) %>% rownames()
})
names(down) <- paste0(colnames(dtest), " - down")

up_comb <- ComplexHeatmap::make_comb_mat(up, mode = "distinct")
up_upset <- ComplexHeatmap::UpSet(up_comb, comb_col = "red", comb_order = order(comb_size(up_comb), decreasing = T), set_order = names(up))
up_grob <- grid.grabExpr(draw(up_upset)) 
down_comb <- ComplexHeatmap::make_comb_mat(down, mode = "distinct")
down_upset <- ComplexHeatmap::UpSet(down_comb, comb_col = "blue", comb_order = order(comb_size(down_comb), decreasing = T), set_order = names(down))
down_grob <- grid.grabExpr(draw(down_upset)) 
plot_grid(up_grob, down_grob, nrow = 2)

## perform enrichment on the TIMP-1 specific terms above
timp1_specific_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "10"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
timp1_specific_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "10"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
dotplot(timp1_specific_up)
dotplot(timp1_specific_down)

## perform enrichment on the fibrosis specific terms above
#fib_specific_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "1000"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#fib_specific_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "1000"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#dotplot(fib_specific_up)
#dotplot(fib_specific_down)

## perform enrichment on the steatosis specific terms above
lf_specific_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "01"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
lf_specific_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "01"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
dotplot(lf_specific_up)
dotplot(lf_specific_down)

## perform enrichment on common terms above
common_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "11"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
common_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "11"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
dotplot(common_up)
dotplot(common_down)


################## perform enrichment on the TIMP-1 specific terms above ##################

#timp1_specific_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "100"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#timp1_specific_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "100"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
timp1_specific_up_dt <- timp1_specific_up@result
timp1_specific_up_dt <- timp1_specific_up_dt[timp1_specific_up_dt$p.adjust<0.05,]
timp1_specific_up_dt$setSize <- as.numeric(gsub(".*/", "", timp1_specific_up_dt$GeneRatio))
timp1_specific_up_dt$GeneRatio <- timp1_specific_up_dt$Count/timp1_specific_up_dt$setSize
#representative upregulated timp1-specific GOBP: "GO:0022613", "GO:0034470", "GO:0042254", "GO:0006397", "GO:0006364", "GO:0000398", "GO:0006613", "GO:0007030", "GO:0006413", "GO:0045047"	
timp1_toplot_up <- timp1_specific_up_dt[timp1_specific_up_dt$ID %in% 
                                          c("GO:0022613", "GO:0034470", "GO:0042254", 
                                            "GO:0006397", "GO:0006364", "GO:0000398", 
                                            "GO:0006613", "GO:0007030", "GO:0006413", "GO:0045047"),] 
#timp1_toplot_up <- timp1_specific_up_dt %>% slice_min(p.adjust, n = 10)
timp1_toplot_up$logP <- -log10(timp1_toplot_up$p.adjust)
timp1_toplot_up$Direction <- "up"

timp1_specific_down_dt <- timp1_specific_down@result
timp1_specific_down_dt <- timp1_specific_down_dt[timp1_specific_down_dt$p.adjust<0.05,]
timp1_specific_down_dt$setSize <- as.numeric(gsub(".*/", "", timp1_specific_down_dt$GeneRatio))
timp1_specific_down_dt$GeneRatio <- timp1_specific_down_dt$Count/timp1_specific_down_dt$setSize
#representative downregulated timp1-specific GOBP: "GO:0009060", "GO:0046395", "GO:0006119", "GO:0009063", "GO:0006754", "GO:0009150", "GO:0008202", "GO:0022904", "GO:0008203", "GO:0006699"
timp1_toplot_down <- timp1_specific_down_dt[timp1_specific_down_dt$ID %in% 
                                              c("GO:0009060", "GO:0046395", "GO:0006119", "GO:0009063", "GO:0006754", 
                                                "GO:0009150", "GO:0008202", "GO:0022904", "GO:0008203", "GO:0006699"),] 
#timp1_toplot_down <- timp1_specific_down_dt %>% slice_min(p.adjust, n = 10)
timp1_toplot_down$logP <- log10(timp1_toplot_down$p.adjust)
timp1_toplot_down$Direction <- "down"

timp1_toplot <- rbind(timp1_toplot_up, timp1_toplot_down)

#gs_toplot$nleading <- sapply(strsplit(gs_toplot$core_enrichment, "/"), length)
#gs_toplot$`Gene Ratio` <- gs_toplot$nleading/gs_toplot$setSize
timp1_toplot$Description <- gsub('(.{1,40})(\\s|$)', '\\1\n', timp1_toplot$Description)
timp1_toplot$Description <- gsub('\n$)', '', timp1_toplot$Description)
timp1_toplot$Description <- factor(timp1_toplot$Description, levels = timp1_toplot$Description)
ora_timp1_specific <- timp1_toplot %>%  ggplot(aes(x = reorder(Description,logP), y = logP, color = Direction, size=GeneRatio)) +
  #geom_hline(yintercept = 0) +
  geom_point() +
  coord_flip() +
  scale_color_manual(values = c("blue","red")) +
  #scale_size_continuous(range = c(2,10)) +
  theme_cowplot() +
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
        legend.key.size = unit(0.5, "cm"))+
  ylab("Signed log10(p.adjust)") +
  xlab("GO category") +
  ggtitle("Timp1 specific GOBP")
ora_timp1_specific



#enrichplot::dotplot(timp1_specific_up)
#enrichplot::dotplot(timp1_specific_down)

#lw_specific_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "010"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#lw_specific_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "010"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#lw_specific_up_dt <- lw_specific_up@result
#lw_specific_up_dt <- lw_specific_up_dt[lw_specific_up_dt$p.adjust<0.05,]
#lw_specific_up_dt$setSize <- as.numeric(gsub(".*/", "", lw_specific_up_dt$GeneRatio))
#lw_specific_up_dt$GeneRatio <- lw_specific_up_dt$Count/lw_specific_up_dt$setSize
#lw_toplot_up <- lw_specific_up_dt[lw_specific_up_dt$ID %in% 
#                                              c("GO:0019395", "GO:0006635", "GO:0006637", "GO:0006734", "GO:0019674", 
#                                                "GO:0006163", "GO:0019362", "GO:0046320", "GO:0046486", "GO:0006695"),] 

#lw_toplot_up <- lw_specific_up_dt %>% slice_min(p.adjust, n = 10)
#lw_toplot_up$logP <- -log10(lw_toplot_up$p.adjust)
#lw_toplot_up$Direction <- "up"

#lw_specific_down_dt <- lw_specific_down@result
#lw_specific_down_dt <- lw_specific_down_dt[lw_specific_down_dt$p.adjust<0.05,]
#lw_specific_down_dt$setSize <- as.numeric(gsub(".*/", "", lw_specific_down_dt$GeneRatio))
#lw_specific_down_dt$GeneRatio <- lw_specific_down_dt$Count/lw_specific_down_dt$setSize
#lw_toplot_down <- lw_specific_down_dt %>% slice_min(p.adjust, n = 10)
#lw_toplot_down$logP <- log10(lw_toplot_down$p.adjust)
#lw_toplot_down$Direction <- "down"

#lw_toplot <- rbind(lw_toplot_up, lw_toplot_down)

#gs_toplot$nleading <- sapply(strsplit(gs_toplot$core_enrichment, "/"), length)
#gs_toplot$`Gene Ratio` <- gs_toplot$nleading/gs_toplot$setSize
#lw_toplot$Description <- gsub('(.{1,40})(\\s|$)', '\\1\n', lw_toplot$Description)
#lw_toplot$Description <- gsub('\n$)', '', lw_toplot$Description)
#lw_toplot$Description <- factor(lw_toplot$Description, levels = lw_toplot$Description)
#ora_lw_specific <- lw_toplot %>%  ggplot(aes(x = reorder(Description,logP), y = logP, color = Direction, size=GeneRatio)) +
#geom_hline(yintercept = 0) +
#  geom_point() +
#  coord_flip() +
#  scale_color_manual(values = c("blue","red")) +
#scale_size_continuous(range = c(2,10)) +
#  theme_cowplot() +
#  theme(panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        plot.title = element_text(size = 8),
#axis.line = element_blank(),
#axis.ticks.x = element_blank(),
#        axis.line = element_line(colour = 'black', size = 0.3),
#        axis.text.x = element_text(size = 8),
#        axis.text.y = element_text(size = 8),
#        axis.title.y = element_text(size=8),
#        axis.title.x = element_text(size = 8),
#        legend.text=element_text(size=8),
#        legend.title = element_text(size=8),
#legend.position = "bottom",
#        panel.background = element_blank(),
#strip.text = element_text(size = 6),
#        legend.key.size = unit(0.5, "cm"))+
#  ylab("Signed log10(p.adjust)") +
#  xlab("GO category") +
# ggtitle("Liver weight specific GOBP")
#ora_lw_specific


#enrichplot::dotplot(lw_specific_up)
#enrichplot::dotplot(lw_specific_down)

#steatosis_specific_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "001"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#steatosis_specific_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "001"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
steatosis_specific_up_dt <- lf_specific_up@result
steatosis_specific_up_dt <- steatosis_specific_up_dt[steatosis_specific_up_dt$p.adjust<0.05,]
steatosis_specific_up_dt$setSize <- as.numeric(gsub(".*/", "", steatosis_specific_up_dt$GeneRatio))
steatosis_specific_up_dt$GeneRatio <- steatosis_specific_up_dt$Count/steatosis_specific_up_dt$setSize
steatosis_toplot_up <- steatosis_specific_up_dt %>% slice_min(p.adjust, n = 10)
steatosis_toplot_up$logP <- -log10(steatosis_toplot_up$p.adjust)
steatosis_toplot_up$Direction <- "up"

steatosis_specific_down_dt <- lf_specific_down@result
steatosis_specific_down_dt <- steatosis_specific_down_dt[steatosis_specific_down_dt$p.adjust<0.05,]
steatosis_specific_down_dt$setSize <- as.numeric(gsub(".*/", "", steatosis_specific_down_dt$GeneRatio))
steatosis_specific_down_dt$GeneRatio <- steatosis_specific_down_dt$Count/steatosis_specific_down_dt$setSize
steatosis_toplot_down <- steatosis_specific_down_dt[steatosis_specific_down_dt$ID %in% 
                                                      c("GO:0006635", "GO:0016042", "GO:0033865", "GO:0046320", "GO:0006639", 
                                                        "GO:0006637", "GO:0000038", "GO:0006641", "GO:0046503", "GO:0034389"),]
#steatosis_toplot_down <- steatosis_specific_down_dt %>% slice_min(p.adjust, n = 10)
steatosis_toplot_down$logP <- log10(steatosis_toplot_down$p.adjust)
steatosis_toplot_down$Direction <- "down"

steatosis_toplot <- rbind(steatosis_toplot_up, steatosis_toplot_down)

#gs_toplot$nleading <- sapply(strsplit(gs_toplot$core_enrichment, "/"), length)
#gs_toplot$`Gene Ratio` <- gs_toplot$nleading/gs_toplot$setSize
steatosis_toplot$Description <- gsub('(.{1,40})(\\s|$)', '\\1\n', steatosis_toplot$Description)
steatosis_toplot$Description <- gsub('\n$)', '', steatosis_toplot$Description)
steatosis_toplot$Description <- factor(steatosis_toplot$Description, levels = steatosis_toplot$Description)
ora_steatosis_specific <- steatosis_toplot %>%  ggplot(aes(x = reorder(Description,logP), y = logP, color = Direction, size=GeneRatio)) +
  #geom_hline(yintercept = 0) +
  geom_point() +
  coord_flip() +
  scale_color_manual(values = c("blue","red")) +
  #scale_size_continuous(range = c(2,10)) +
  theme_cowplot() +
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
        legend.key.size = unit(0.5, "cm"))+
  ylab("Signed log10(p.adjust)") +
  xlab("GO category") +
  ggtitle("Steatosis specific GOBP")
ora_steatosis_specific

#enrichplot::dotplot(steatosis_specific_up)
#enrichplot::dotplot(steatosis_specific_down)

#steatosis_lw_specific_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "011"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#steatosis_lw_specific_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "011"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#steatosis_lw_specific_up_dt <- steatosis_lw_specific_up@result
#steatosis_lw_specific_up_dt <- steatosis_lw_specific_up_dt[steatosis_lw_specific_up_dt$p.adjust<0.05,]
#steatosis_lw_specific_up_dt$setSize <- as.numeric(gsub(".*/", "", steatosis_lw_specific_up_dt$GeneRatio))
#steatosis_lw_specific_up_dt$GeneRatio <- steatosis_lw_specific_up_dt$Count/steatosis_lw_specific_up_dt$setSize
#steatosis_lw_toplot_up <- steatosis_lw_specific_up_dt %>% slice_min(p.adjust, n = 10)
#steatosis_lw_toplot_up$logP <- -log10(steatosis_lw_toplot_up$p.adjust)
#steatosis_lw_toplot_up$Direction <- "up"

#steatosis_lw_specific_down_dt <- steatosis_lw_specific_down@result
#steatosis_lw_specific_down_dt <- steatosis_lw_specific_down_dt[steatosis_lw_specific_down_dt$p.adjust<0.05,]
#steatosis_lw_specific_down_dt$setSize <- as.numeric(gsub(".*/", "", steatosis_lw_specific_down_dt$GeneRatio))
#steatosis_lw_specific_down_dt$GeneRatio <- steatosis_lw_specific_down_dt$Count/steatosis_lw_specific_down_dt$setSize
#steatosis_lw_toplot_down <- steatosis_lw_specific_down_dt %>% slice_min(p.adjust, n = 10)
#steatosis_lw_toplot_down$logP <- log10(steatosis_lw_toplot_down$p.adjust)
#steatosis_lw_toplot_down$Direction <- "down"

#steatosis_lw_toplot <- rbind(steatosis_lw_toplot_up, steatosis_lw_toplot_down)

#gs_toplot$nleading <- sapply(strsplit(gs_toplot$core_enrichment, "/"), length)
#gs_toplot$`Gene Ratio` <- gs_toplot$nleading/gs_toplot$setSize
#steatosis_lw_toplot$Description <- gsub('(.{1,40})(\\s|$)', '\\1\n', steatosis_lw_toplot$Description)
#steatosis_lw_toplot$Description <- gsub('\n$)', '', steatosis_lw_toplot$Description)
#steatosis_lw_toplot$Description <- factor(steatosis_lw_toplot$Description, levels = steatosis_lw_toplot$Description)
#ora_steatosis_lw_specific <- steatosis_lw_toplot %>%  ggplot(aes(x = reorder(Description,logP), y = logP, color = Direction, size=GeneRatio)) +
#geom_hline(yintercept = 0) +
#  geom_point() +
#  coord_flip() +
#  scale_color_manual(values = c("blue","red")) +
#scale_size_continuous(range = c(2,10)) +
#  theme_cowplot() +
#  theme(panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        plot.title = element_text(size = 8),
#axis.line = element_blank(),
#axis.ticks.x = element_blank(),
#        axis.line = element_line(colour = 'black', size = 0.3),
#        axis.text.x = element_text(size = 8),
#        axis.text.y = element_text(size = 8),
#        axis.title.y = element_text(size=8),
#        axis.title.x = element_text(size = 8),
#        legend.text=element_text(size=8),
#        legend.title = element_text(size=8),
#legend.position = "bottom",
#        panel.background = element_blank(),
#strip.text = element_text(size = 6),
#        legend.key.size = unit(0.5, "cm"))+
#  ylab("Signed log10(p.adjust)") +
#  xlab("GO category") +
#  ggtitle("Liver weight&steatosis specific GOBP")
#ora_steatosis_lw_specific
#enrichplot::dotplot(steatosis_lw_specific_up)
#enrichplot::dotplot(steatosis_lw_specific_down)

#steatosis_timp1_specific_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "101"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#steatosis_timp1_specific_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "101"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#steatosis_timp1_specific_up_dt <- steatosis_timp1_specific_up@result
#steatosis_timp1_specific_up_dt <- steatosis_timp1_specific_up_dt[steatosis_timp1_specific_up_dt$p.adjust<0.05,]
#steatosis_timp1_specific_up_dt$setSize <- as.numeric(gsub(".*/", "", steatosis_timp1_specific_up_dt$GeneRatio))
#steatosis_timp1_specific_up_dt$GeneRatio <- steatosis_timp1_specific_up_dt$Count/steatosis_timp1_specific_up_dt$setSize
#steatosis_timp1_toplot_up <- steatosis_timp1_specific_up_dt %>% slice_min(p.adjust, n = 10)
#steatosis_timp1_toplot_up$logP <- -log10(steatosis_timp1_toplot_up$p.adjust)
#steatosis_timp1_toplot_up$Direction <- "up"

#steatosis_timp1_specific_down_dt <- steatosis_timp1_specific_down@result
#steatosis_timp1_specific_down_dt <- steatosis_timp1_specific_down_dt[steatosis_timp1_specific_down_dt$p.adjust<0.05,]
#steatosis_timp1_specific_down_dt$setSize <- as.numeric(gsub(".*/", "", steatosis_timp1_specific_down_dt$GeneRatio))
#steatosis_timp1_specific_down_dt$GeneRatio <- steatosis_timp1_specific_down_dt$Count/steatosis_timp1_specific_down_dt$setSize
#steatosis_timp1_toplot_down <- steatosis_timp1_specific_down_dt %>% slice_min(p.adjust, n = 10)
#steatosis_timp1_toplot_down$logP <- log10(steatosis_timp1_toplot_down$p.adjust)
#steatosis_timp1_toplot_down$Direction <- "down"

#steatosis_timp1_toplot <- rbind(steatosis_timp1_toplot_up, steatosis_timp1_toplot_down)

#gs_toplot$nleading <- sapply(strsplit(gs_toplot$core_enrichment, "/"), length)
#gs_toplot$`Gene Ratio` <- gs_toplot$nleading/gs_toplot$setSize
#steatosis_timp1_toplot$Description <- gsub('(.{1,40})(\\s|$)', '\\1\n', steatosis_timp1_toplot$Description)
#steatosis_timp1_toplot$Description <- gsub('\n$)', '', steatosis_timp1_toplot$Description)
#steatosis_timp1_toplot$Description <- factor(steatosis_timp1_toplot$Description, levels = steatosis_timp1_toplot$Description)
#ora_steatosis_timp1_specific <- steatosis_timp1_toplot %>%  ggplot(aes(x = reorder(Description,logP), y = logP, color = Direction, size=GeneRatio)) +
#geom_hline(yintercept = 0) +
#  geom_point() +
#  coord_flip() +
#  scale_color_manual(values = c("red")) +
#scale_size_continuous(range = c(2,10)) +
#  theme_cowplot() +
#  theme(panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        plot.title = element_text(size = 8),
#axis.line = element_blank(),
#axis.ticks.x = element_blank(),
#       axis.line = element_line(colour = 'black', size = 0.3),
#        axis.text.x = element_text(size = 8),
#        axis.text.y = element_text(size = 8),
#        axis.title.y = element_text(size=8),
#       axis.title.x = element_text(size = 8),
#        legend.text=element_text(size=8),
#        legend.title = element_text(size=8),
#legend.position = "bottom",
#        panel.background = element_blank(),
#strip.text = element_text(size = 6),
#        legend.key.size = unit(0.5, "cm"))+
#  ylab("Signed log10(p.adjust)") +
#  xlab("GO category") +
#  ggtitle("Steatosis&Timp1 specific GOBP")
#ora_steatosis_timp1_specific
#enrichplot::dotplot(steatosis_lw_specific_up)
#enrichplot::dotplot(steatosis_lw_specific_down)



#common_up <- clusterProfiler::enrichGO(extract_comb(up_comb, "111"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
#common_down <- clusterProfiler::enrichGO(extract_comb(down_comb, "111"),  ont = "BP", universe = rownames(liver_dge$cpm), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")

common_up_dt <- common_up@result
common_up_dt <- common_up_dt[common_up_dt$p.adjust<0.05,]
common_up_dt$setSize <- as.numeric(gsub(".*/", "", common_up_dt$GeneRatio))
common_up_dt$GeneRatio <- common_up_dt$Count/common_up_dt$setSize
common_toplot_up <- common_up_dt[common_up_dt$ID %in% 
                                   c("GO:0050900", "GO:0007015", "GO:0001819", "GO:0032943", "GO:0030198", 
                                     "GO:0032635", "GO:0050863", "GO:0042060", "GO:0044770", "GO:0019221",
                                     "GO:0043299","GO:0007093", "GO:1905517", "GO:0050729","GO:0032963",
                                     "GO:0045089","GO:0006260","GO:0034341","GO:0038065","GO:0072593"),]
#common_toplot_up <- common_up_dt %>% slice_min(p.adjust, n = 10)
common_toplot_up$logP <- -log10(common_toplot_up$p.adjust)
common_toplot_up$Direction <- "up"

#common_down_dt <- common_down@result
#common_down_dt <- common_down_dt[common_down_dt$p.adjust<0.05,]
#common_down_dt$setSize <- as.numeric(gsub(".*/", "", common_down_dt$GeneRatio))
#common_down_dt$GeneRatio <- common_down_dt$Count/common_down_dt$setSize
#common_toplot_down <- common_down_dt %>% slice_min(p.adjust, n = 10)
c#ommon_toplot_down$logP <- log10(common_toplot_down$p.adjust)
#common_toplot_down$Direction <- "down"

common_toplot <- common_toplot_up

#gs_toplot$nleading <- sapply(strsplit(gs_toplot$core_enrichment, "/"), length)
#gs_toplot$`Gene Ratio` <- gs_toplot$nleading/gs_toplot$setSize
common_toplot$Description <- gsub('(.{1,40})(\\s|$)', '\\1\n', common_toplot$Description)
common_toplot$Description <- gsub('\n$)', '', common_toplot$Description)
common_toplot$Description <- factor(common_toplot$Description, levels = common_toplot$Description)
ora_common <- common_toplot %>%  ggplot(aes(x = reorder(Description,logP), y = logP, color = Direction, size=GeneRatio)) +
  #geom_hline(yintercept = 0) +
  geom_point() +
  coord_flip() +
  scale_color_manual(values = c("red")) +
  #scale_size_continuous(range = c(2,10)) +
  theme_cowplot() +
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
        legend.key.size = unit(0.5, "cm"))+
  ylab("Signed log10(p.adjust)") +
  xlab("GO category") +
  ggtitle("Common GOBP")
ora_common

pg <- plot_grid(ora_timp1_specific,ora_steatosis_specific,ora_common, nrow=1)
pg
ggsave(filename ="17-ORA_plots.pdf",plot=pg,path ="./Plots",width = 18, height = 7, units = c("in"))
ggsave(filename ="17-timp1_specific_ORA.pdf",plot=ora_timp1_specific,path ="./Plots",width = 4, height = 3.8, units = c("in"))
#ggsave(filename ="17-lw_specific_ORA.pdf",plot=ora_lw_specific,path ="./Plots",width = 3.9, height = 3, units = c("in"))
ggsave(filename ="17-steatosis_specific_ORA.pdf",plot=ora_steatosis_specific,path ="./Plots",width = 3.8, height = 2.5, units = c("in"))
ggsave(filename ="17-common_ORA.pdf",plot=ora_common,path ="./Plots",width = 4.2, height = 3.8, units = c("in"))




################# do gsea for Timp1########################
gst <- gseGO(geneList = setNames(tttimp$signed_log_p, tttimp$gene_id), ont = "BP", OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",pvalueCutoff = 1)
pt <- gst@result %>% filter(p.adjust < 0.05) %>% ggplot(aes(x = -log10(p.adjust), y = NES, desc = Description)) + geom_point()
ggplotly(pt)
gst@result %>% 
  filter(p.adjust < 0.05) %>% 
  filter(NES < 0) %>% 
  select(-core_enrichment, - leading_edge)

gst@result %>% 
  filter(p.adjust < 0.05) %>% 
  filter(NES < 0) %>% 
  
  head(gst@result %>% select(Description, NES, p.adjust), n = 50)

enrichplot::dotplot(gst, showCategory = 10)

############# plot top 10 up and down terms #############
gst_toplot <- gst@result %>% split(., sign(gst@result$NES)) %>% lapply(., head, n = 10) %>% do.call(rbind, .) %>% arrange(NES)

gst_toplot$nleading <- sapply(strsplit(gst_toplot$core_enrichment, "/"), length)
gst_toplot$`Gene Ratio` <- gst_toplot$nleading/gst_toplot$setSize
gst_toplot$Description <- gsub('(.{1,30})(\\s|$)', '\\1\n', gst_toplot$Description)
gst_toplot$Description <- gsub('\n$)', '', gst_toplot$Description)
gst_toplot$Description <- factor(gst_toplot$Description, levels = gst_toplot$Description)

gst_toplot %>%  ggplot(aes(x = Description, y = NES, color = -log10(p.adjust), size=`Gene Ratio`)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  coord_flip() +
  scale_color_viridis_c(limits = c(3,NA), direction = -1, option = "magma") +
  scale_size_continuous(range = c(2,10)) +
  theme_cowplot() +
  ylab("Normalized Enrichment Score (NES)") +
  xlab("GO category") +
  ggtitle("Top enriched biological processes\nin genes associated with plasma TIMP-1 levels")
ggsave(useDingbats = F, filename = "./Plots/17-liver_dge_timp1_top_gsea.pdf", width = 7.5, height = 9)


################## now do sex-specific timp1 DE #################
timp_sex_specific_de <- lapply(c("m","f"), function(sx){
  dtimp <- model.matrix(~log(TIMP1_24) + Liver_g_RelBWSac_24 + sum_all_vacuoles_percentage_24, 
                        data = liver_dge$samples %>% filter(Sex == sx) %>% 
                          filter(Generation == "F2") %>% filter(!is.na(TIMP1_24)) %>% 
                          filter(!is.na(Liver_g_RelBWSac_24))  %>% filter(!is.na(sum_all_vacuoles_percentage_24)) )
  # dtimp_voom <- model.matrix(~ Sex, data = liver_dge$samples %>% filter(Generation == "F2") %>% filter(!is.na(TIMP.1)) %>% filter(!is.na(Liver_g_RelBWSac_24)))
  vm_timp <- voom(liver_dge$counts_filtered[,(liver_dge$samples$Sex ==sx) & 
                                              (liver_dge$samples$Generation=="F2") & 
                                              (!is.na(liver_dge$samples$TIMP1_24)) & 
                                              (!is.na(liver_dge$samples$Liver_g_RelBWSac_24)) & 
                                              (!is.na(liver_dge$samples$sum_all_vacuoles_percentage_24)) ], design = dtimp)
  fit_timp <- lmFit(vm_timp, design = dtimp)
  fit_timp <- eBayes(fit_timp, robust=TRUE)
  tt <- topTable(fit_timp, coef = "log(TIMP1_24)", n = Inf) %>% 
    rownames_to_column(var = "gene_id") %>% 
    left_join(., liver_dge$geneConversionTable, by = "gene_id") %>% 
    mutate(signed_log_p = sign(logFC) * -log10(adj.P.Val)) %>% 
    dplyr::arrange(-signed_log_p)
  return(list(fit=fit_timp, tt = tt))
})
names(timp_sex_specific_de) <- c("m","f")

timp_sex_specific_tt <- left_join(timp_sex_specific_de$m$tt, timp_sex_specific_de$f$tt, by = c("gene_id", "gene_name","seqnames", "gene_biotype"), suffix = c(".m", ".f"))
timp_sex_specific_tt <- timp_sex_specific_tt %>% 
  mutate(sigclass = paste0(adj.P.Val.m < 0.05, adj.P.Val.f < 0.05) %>% 
           plyr::mapvalues(., from = c("TRUETRUE", "TRUEFALSE", "FALSETRUE", "FALSEFALSE"), to = c("Both", "Only Males", "Only Females", "Not. sig")))
timp_sex_specific_tt$sigclass <- factor(timp_sex_specific_tt$sigclass, levels = c("Both", "Only Males", "Only Females", "Not. sig"), ordered = T)

sigclass_colors <- setNames(c("orange", unname(sexColors), "black"), c("Both", "Only Males", "Only Females", "Not. sig"))

############ label some top genes in each group ############
headtail <- function(x,n=10){c(head(x,n), tail(x,n))}
tolabel <- sapply(split(timp_sex_specific_tt, timp_sex_specific_tt$sigclass), function(x){
  if(x$sigclass[1] == "Both"){
    out<- x %>% mutate(averageP = (logFC.m + logFC.f)/2) %>% dplyr::arrange(averageP) %>% pull(gene_name) %>% headtail(n=20)
  }
  if(x$sigclass[1] == "Only Males"){
    out <- x %>%  dplyr::arrange(logFC.m) %>% pull(gene_name) %>% headtail(n=20) %>% return()
  }
  if(x$sigclass[1] == "Only Females"){
    out<- x %>%  dplyr::arrange(logFC.f) %>% pull(gene_name) %>% headtail(n=20) %>% return()
  }
  if(x$sigclass[1] ==  "Not. sig"){
    out<-NULL
  }
  return(out)
})

tolabel <- tolabel %>% unlist() %>% unique()
timp_sex_specific_tt$label <- timp_sex_specific_tt$gene_name
timp_sex_specific_tt$label[!timp_sex_specific_tt$label %in% tolabel] = ""

vmale <- ggplot(timp_sex_specific_tt, aes(x = logFC.m, y = -log10(adj.P.Val.m), color = sigclass
                                          , label = label
)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_label_repel(data = timp_sex_specific_tt  %>% filter(sigclass == "Only Males"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() +
  scale_color_manual(values = sigclass_colors, name = "Significant in:")
vmale
vfemale <- ggplot(timp_sex_specific_tt, aes(x = logFC.f, y = -log10(adj.P.Val.f), color = sigclass, label = label)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_label_repel(data = timp_sex_specific_tt  %>% filter(sigclass == "Only Females"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() +
  scale_color_manual(values = sigclass_colors, name = "Significant in:")
vfemale

logFC_compared <- ggplot(timp_sex_specific_tt, aes(x = logFC.m, y = logFC.f, color = sigclass, label = label)) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_label_repel(data = timp_sex_specific_tt  %>% filter(sigclass == "Both"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() +
  scale_color_manual(values = sigclass_colors, name = "Significant in:")
logFC_compared

vgrid <- plot_grid(get_legend(vmale + theme(legend.position='bottom')),
                   plot_grid(plot_grid(vmale + theme(legend.position = "none"),
                                       vfemale + theme(legend.position = "none"),
                                       align = 'hv', nrow = 1),
                             logFC_compared  + theme(legend.position = "none"),
                             nrow = 2, axis  = "tblr"),
                   nrow = 2, rel_heights = c(1,15))

ggsave(useDingbats = F, plot = vgrid, filename = "./Plots/11-volcano_timp1_male_female_compare.pdf", width = 15, height = 10)

##################volacano plot for Johan, genes############################
x <- timp_sex_specific_tt
x <- x %>% 
  mutate(Direction_m = paste0(adj.P.Val.m < 0.05, logFC.m > 0) %>% 
           plyr::mapvalues(., from = c("TRUETRUE", "TRUEFALSE", "FALSETRUE", "FALSEFALSE"), to = c("Sig_up", "Sig_down", "Not. sig", "Not. sig")))
x$Direction_m <- factor(x$Direction_m, levels = c("Sig_up", "Sig_down", "Not. sig"), ordered = T)


vmale <- ggplot(x, aes(x = logFC.m, y = -log10(adj.P.Val.m), color = Direction_m
                       #, label = gene_name
)) +
  geom_hline(yintercept = -log10(0.05), color="white",linetype="dashed") +
  geom_vline(xintercept = 0, color="white") +
  geom_point() +
  #geom_label_repel(data = timp_sex_specific_tt  %>% filter(sigclass == "Only Males"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() +
  theme(#
    #axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background =  element_rect(fill="black", colour="black"),
    panel.background = element_rect(fill="black", colour="black"),
    axis.text=element_text(colour = 'white', size=10),
    axis.title=element_text(colour = 'white', size=10),
    legend.text=element_text(colour = 'white', size=10),
    legend.title = element_text(colour = 'white', size=10),
    axis.line = element_line(colour = 'white', size = 0.3),
    axis.ticks = element_line(colour = 'white', size=0.3),
    #axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    legend.key.size = unit(0.1, "cm"))+
  scale_color_manual(values = c("red","#29A7DF","grey"))
vmale
ggsave(plot = vmale, filename = "./Plots/Liver/Giorgia/GB02-liver_dge_timp1_volcano_genes_males.pdf", width = 5.7, height = 3.5)

write.csv(timp_sex_specific_tt, "./Data/sex_specific_toptable.csv", row.names=FALSE)
########### gsea on the sex-specific timp-1 association ###############

timp_sex_specific_de_gsea <- lapply(timp_sex_specific_de, function(x){
  gl <- x$tt %>% arrange(-signed_log_p)
  gl <- setNames(gl$signed_log_p, gl$gene_id)
  ggo <- gseGO(geneList = gl, ont = "BP", OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", pvalueCutoff = 1)
  ggo@result <- ggo@result %>% mutate(signed_log_p = -log10(p.adjust) * sign(NES))
  ggo
})

volcano_gsea_dt <- do.call(rbind, lapply(c("m","f"), function(x)timp_sex_specific_de_gsea[[x]]@result %>% mutate(Sex = x)))
volcano_gsea_sexcompare_dt <- volcano_gsea_dt %>% 
  pivot_wider(id_cols = c("Description"), names_from = c("Sex"), values_from = c("NES", "p.adjust", "signed_log_p")) %>%
  mutate(sigclass = paste0(p.adjust_m < 0.05, p.adjust_f < 0.05) %>% 
           plyr::mapvalues(., from = c("TRUETRUE", "TRUEFALSE", "FALSETRUE", "FALSEFALSE"), to = c("Both", "Only Males", "Only Females", "Not. sig")))
volcano_gsea_sexcompare_dt$sigclass <- factor(volcano_gsea_sexcompare_dt$sigclass, levels = c("Both", "Only Males", "Only Females", "Not. sig"), ordered = T)
tolabel_gsea <- sapply(split(volcano_gsea_sexcompare_dt, volcano_gsea_sexcompare_dt$sigclass), function(x){
  if(x$sigclass[1] == "Both"){
    out<- x %>% mutate(averageP = (NES_m + NES_f)/2) %>% dplyr::arrange(averageP) %>% pull(Description) %>% headtail(n=10)
  }
  if(x$sigclass[1] == "Only Males"){
    out <- x %>%  dplyr::arrange(NES_m) %>% pull(Description) %>% headtail(n=10) %>% return()
  }
  if(x$sigclass[1] == "Only Females"){
    out<- x %>%  dplyr::arrange(NES_f) %>% pull(Description) %>% headtail(n=10) %>% return()
  }
  if(x$sigclass[1] ==  "Not. sig"){
    out<-NULL
  }
  return(out)
}) %>% unlist()
volcano_gsea_sexcompare_dt$label <- volcano_gsea_sexcompare_dt$Description
volcano_gsea_sexcompare_dt$label[!volcano_gsea_sexcompare_dt$label %in% tolabel_gsea] <- ""
volcano_gsea_sexcompare_dt$label <- gsub('(.{1,30})(\\s|$)', '\\1\n', volcano_gsea_sexcompare_dt$label)
volcano_gsea_sexcompare_dt$label <- gsub('\n$', '', volcano_gsea_sexcompare_dt$label)


write.csv(volcano_gsea_sexcompare_dt, "./Data/volcano_gsea_sexcompare_dt.csv", row.names=FALSE)

vmale <- ggplot(volcano_gsea_sexcompare_dt, aes(x = -log10(p.adjust_m), y = NES_m, color = sigclass, label = label)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -log10(0.05)) +
  geom_point() +
  geom_label_repel(data = volcano_gsea_sexcompare_dt  %>% filter(sigclass == "Only Males"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() +
  scale_color_manual(values = sigclass_colors, name = "Significant in:")
vmale
vfemale <- ggplot(volcano_gsea_sexcompare_dt, aes(x = -log10(p.adjust_f), y = NES_f, color = sigclass, label = label)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -log10(0.05)) +
  geom_point() +
  geom_label_repel(data = volcano_gsea_sexcompare_dt  %>% filter(sigclass == "Only Females"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() +
  scale_color_manual(values = sigclass_colors, name = "Significant in:")
vfemale

logFC_compared <- ggplot(volcano_gsea_sexcompare_dt, aes(x = NES_m, y = NES_f, color = sigclass, label = label)) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_label_repel(data = volcano_gsea_sexcompare_dt  %>% filter(sigclass == "Both"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() +
  scale_color_manual(values = sigclass_colors, name = "Significant in:")
logFC_compared

vgrid <- plot_grid(get_legend(vmale + theme(legend.position='bottom')),
                   plot_grid(plot_grid(vmale + theme(legend.position = "none"),
                                       vfemale + theme(legend.position = "none"),
                                       align = 'hv', nrow = 1),
                             logFC_compared  + theme(legend.position = "none"),
                             nrow = 2, axis  = "tblr"),
                   nrow = 2, rel_heights = c(1,15))

ggsave(useDingbats = F, plot = vgrid, filename = "./Plots/11-volcano_timp1_gsea_male_female_compare.pdf", width = 15, height = 10)

volcano_gsea_sexcompare_dt %>% arrange(NES_m) %>% group_by(sigclass, keep = T) %>% group_map(~ rbind(head(.x, 5L), tail(.x, 5L)))
volcano_gsea_sexcompare_dt %>% arrange(NES_m) %>% group_by(sigclass) %>% group_split()

################# below some plots for Johan #####################
library(readxl)
library(plyr)
library(ggalt)
library(tidyr)
volcano_gsea_sexcompare_dt <- read_excel("./Data/volcano_gsea_sexcompare_dt_label2.xlsx")
volcano_gsea_sexcompare_dt <- read.csv("./Data/volcano_gsea_sexcompare_dt_label2.csv")
volcano_gsea_sexcompare_dt$label2<- sub("^$", "no label", volcano_gsea_sexcompare_dt$label2)
volcano_gsea_sexcompare_dt$label2 <- volcano_gsea_sexcompare_dt$label2 %>% replace_na('no label')
volcano_gsea_sexcompare_dt$label2 <- factor(volcano_gsea_sexcompare_dt$label2 ,levels = c("no label","Immune","Fibrosis","Sphingolipid Metabolism",
                                                                                          "UPR/ISR/ER stress", "Bile Acid Metabolism", "Mitochondria", "NAD",
                                                                                          "Cholesterol"))


volcano_gsea_sexcompare_dt$order <- ifelse(volcano_gsea_sexcompare_dt$label2=="NA", 1, 2)
gsColors <- c("no label" = "gray70", 
              "Fibrosis" = "#E31A1C",
              #"Sphingolipid Metabolism" = "yellow3",
              "Sphingolipid Metabolism" = "yellow3",
              #"Organ_weight" = "skyblue2", 
              "UPR/ISR/ER stress" = "steelblue4", 
              "Immune" = "gold1",
              "Bile Acid Metabolism" = "#FF7F00",
              #"BW_perGain" = "black",
              "Mitochondria" = "#B2DF8A",
              #"Mitochondria" = "green4",
              "AA Metabolism" = "#6A3D9A", 
              "NAD" = "orchid1", 
              #"Plasma_circ_factors" = "maroon",
              #"Activity" = "#FB9A99",
              "Cholesterol" = "darkturquoise")

vmale <- ggplot(volcano_gsea_sexcompare_dt, aes(x = -log10(p.adjust_m), y = NES_m, 
                                                #color = label2, 
                                                size=-log10(p.adjust_m))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -log10(0.05)) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("no label",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "gray70"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Immune",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "#FB9A99"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Fibrosis",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "#E31A1C"
  ) +
  
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("UPR/ISR/ER stress",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "gold1"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Sphingolipid Metabolism",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "yellow4"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Bile Acid Metabolism",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "#FF7F00"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Mitochondria",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "#B2DF8A"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("NAD",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "orchid1"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Cholesterol",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "darkturquoise"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("AA Metabolism" ,volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "darkturquoise"
  ) +
  #sca
  #scale_color_manual(values = gsColors)+
  #geom_point(aes(color = label2), subset = .(label2 == 'Immune'))+
  #geom_label_repel(data = volcano_gsea_sexcompare_dt  %>% filter(sigclass == "Only Males"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() 
#+
#scale_color_manual(values = sigclass_colors, name = "Significant in:")
vmale
vfemale <- ggplot(volcano_gsea_sexcompare_dt, aes(x = -log10(p.adjust_f), y = NES_f, 
                                                  #color = label2, 
                                                  size=-log10(p.adjust_f))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -log10(0.05)) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("no label",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "gray70"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Immune",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "#FB9A99"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Fibrosis",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "#E31A1C"
  ) +
  
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("UPR/ISR/ER stress",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "gold1"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Sphingolipid Metabolism",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "yellow4"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Bile Acid Metabolism",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "#FF7F00"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Mitochondria",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "#B2DF8A"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("NAD",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "orchid1"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("Cholesterol",volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "darkturquoise"
  ) +
  geom_point(data=volcano_gsea_sexcompare_dt[grepl("AA Metabolism" ,volcano_gsea_sexcompare_dt$label2),], 
             #aes(x = -log10(p.adjust_m), y = NES_m), 
             color = "darkturquoise"
  ) +
  #sca
  #scale_color_manual(values = gsColors)+
  #geom_point(aes(color = label2), subset = .(label2 == 'Immune'))+
  #geom_label_repel(data = volcano_gsea_sexcompare_dt  %>% filter(sigclass == "Only Males"), max.overlaps = 10000, show.legend = F) +
  theme_cowplot() 
vfemale

vboth <- plot_grid(vmale,vfemale, nrow = 1)
vboth
ggsave(plot = vboth, filename = "./Plots/Liver/Giorgia/GB02-volcano_timp1_gsea_male_female_compare.pdf", width = 13, height = 4)

x <- read.csv("./Data/volcano_gsea_sexcompare_dt_sel_m.csv")
y <- read.csv("./Data/volcano_gsea_sexcompare_dt_sel_f.csv")
x$Direction <- "up"
x[which(x$NES_m <0),"Direction"] <- 'down'
y$Direction <- "up"
y[which(y$NES_f <0),"Direction"] <- 'down'

up_m <- ggplot(
  x[grepl("up", x$Direction),],
  aes(y = reorder(Description, -log10(p.adjust_m)), x = -log10(p.adjust_m), col = -log10(p.adjust_m), size = -log10(p.adjust_m))
)+
  geom_point(stat = "identity", show.legend = FALSE)+
  scale_color_gradient(low = "lightgrey", high = "red")+
  #scale_color_manual(values = "red")+
  #facet_grid(~ Direction, scales = "free")+
  ylab("")+
  theme(#
    #axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background =  element_rect(fill="black", colour="black"),
    panel.background = element_rect(fill="black", colour="black"),
    axis.text=element_text(colour = 'white', size=10),
    axis.title=element_text(colour = 'white', size=10),
    legend.text=element_text(colour = 'white', size=10),
    legend.title = element_text(colour = 'white', size=10),
    axis.line = element_line(colour = 'white', size = 0.3),
    axis.ticks = element_line(colour = 'white', size=0.3),
    #axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    legend.key.size = unit(0.1, "cm"))
ggsave(plot = up_m, filename = "./Plots/Liver/Giorgia/GB02-timp1_gsea_up_male_forJohan.pdf", width = 5.5, height = 3)

down_m <- ggplot(
  x[grepl("down", x$Direction),],
  aes(y = reorder(Description, -(-log10(p.adjust_m))), x = -log10(p.adjust_m), col = -log10(p.adjust_m), size = -log10(p.adjust_m))
)+
  geom_point(stat = "identity", show.legend = FALSE)+
  scale_color_gradient(low = "lightgrey", high = "#29A7DF")+
  #scale_color_manual(values = "red")+
  #facet_grid(~ Direction, scales = "free")+
  ylab("")+
  theme(#
    #axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background =  element_rect(fill="black", colour="black"),
    panel.background = element_rect(fill="black", colour="black"),
    axis.text=element_text(colour = 'white', size=10),
    axis.title=element_text(colour = 'white', size=10),
    legend.text=element_text(colour = 'white', size=10),
    legend.title = element_text(colour = 'white', size=10),
    axis.line = element_line(colour = 'white', size = 0.3),
    axis.ticks = element_line(colour = 'white', size=0.3),
    #axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    legend.key.size = unit(0.1, "cm"))
ggsave(plot = down_m, filename = "./Plots/Liver/Giorgia/GB02-timp1_gsea_down_male_forJohan.pdf", width = 6.4, height = 3)

up_f <- ggplot(
  y[grepl("up", x$Direction),],
  aes(y = reorder(Description, -log10(p.adjust_f)), x = -log10(p.adjust_f), col = -log10(p.adjust_f), size = -log10(p.adjust_f))
)+
  geom_point(stat = "identity", show.legend = FALSE)+
  scale_color_gradient(low = "lightgrey", high = "red")+
  #scale_color_manual(values = "red")+
  #facet_grid(~ Direction, scales = "free")+
  ylab("")+
  theme(#
    #axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background =  element_rect(fill="black", colour="black"),
    panel.background = element_rect(fill="black", colour="black"),
    axis.text=element_text(colour = 'white', size=10),
    axis.title=element_text(colour = 'white', size=10),
    legend.text=element_text(colour = 'white', size=10),
    legend.title = element_text(colour = 'white', size=10),
    axis.line = element_line(colour = 'white', size = 0.3),
    axis.ticks = element_line(colour = 'white', size=0.3),
    #axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    legend.key.size = unit(0.1, "cm"))
ggsave(plot = up_f, filename = "./Plots/Liver/Giorgia/GB02-timp1_gsea_up_female_forJohan.pdf", width = 5.5, height = 3)

down_f <- ggplot(
  y[grepl("down", x$Direction),],
  aes(y = reorder(Description, -(-log10(p.adjust_f))), x = -log10(p.adjust_f), col = -log10(p.adjust_f), size = -log10(p.adjust_f))
)+
  geom_point(stat = "identity", show.legend = FALSE)+
  scale_color_gradient(low = "lightgrey", high = "cyan")+
  #scale_color_manual(values = "red")+
  #facet_grid(~ Direction, scales = "free")+
  ylab("")+
  theme(#
    #axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background =  element_rect(fill="black", colour="black"),
    panel.background = element_rect(fill="black", colour="black"),
    axis.text=element_text(colour = 'white', size=10),
    axis.title=element_text(colour = 'white', size=10),
    legend.text=element_text(colour = 'white', size=10),
    legend.title = element_text(colour = 'white', size=10),
    axis.line = element_line(colour = 'white', size = 0.3),
    axis.ticks = element_line(colour = 'white', size=0.3),
    #axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(.07, "cm"),
    axis.ticks.length.x = unit(.07, "cm"),
    #axis.ticks.x = element_blank(),
    legend.key.size = unit(0.1, "cm"))
ggsave(plot = down_f, filename = "./Plots/Liver/Giorgia/GB02-timp1_gsea_down_female_forJohan.pdf", width = 6, height = 3)

gsp <- plot_grid(up_m, up_f, down_m,down_f, nrow=2)
gsp_m <- plot_grid(up_m, down_m, nrow=2)
############### overlap with human liver datasets #######################
source("./Scripts/mouse_to_human_orthodb.R")

all_human_toptables <- readRDS("~/common/Users/vonalven/aim2a_step1_giacomo/Data/final_results_mouse_human_comparison/DEA_results.RDS")
#all_human_toptables <- all_human_toptables[c("GSE135251", "GSE162694")]

#out <- all_human_toptables$GSE135251$tt %>% rownames_to_column(var = "gene_id_human")
#cnv <- convert_human_id_mouse_orthodb(gns = out$gene_id_human)

all_human_toptables_vs_mouse <- setNames(lapply(names(all_human_toptables), function(x){
  out <- all_human_toptables[[x]]$tt %>% rownames_to_column(var = "gene_id_human")
  # convert to mouse IDs
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://asia.ensembl.org/")
  orth <- getBM(
    attributes = c(
      "ensembl_gene_id", 
      "mmusculus_homolog_ensembl_gene"
    ),
    filters = "ensembl_gene_id",
    values = out$gene_id_human,   # your vector
    mart = human
  )
  colnames(orth) <- c("gene_id_human","gene_id")
  gct <- read.csv(file = "./Data/geneConversionTable.csv")
  orth$gene_name <- gct$gene_name[match(orth$gene_id, gct$gene_id)]
  lj <- left_join(out, orth, by = "gene_id_human")
  lj
  #left_join(ttlf, lj, by = "gene_name") %>% mutate(dataset = x)
}),
nm = names(all_human_toptables)  # preserves original names
)

all_human_toptables_vs_mouse <- do.call(rbind, all_human_toptables_vs_mouse)
all_human_toptables_vs_mouse <- all_human_toptables_vs_mouse %>% filter(!is.na(gene_name))
all_human_toptables_vs_mouse_plot_m <- all_human_toptables_vs_mouse %>% 
  filter(adj.P.Val < 0.05) %>%
  ggplot(aes(x = logFC, color = sigclass)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(aes(y = logFC.m, gene_name = gene_name)) +
  facet_grid(~dataset) +
  scale_color_manual(values = sigclass_colors) +
  theme_cowplot() +
  xlab("logFC in human NAS>=4 vs NAS<4") +
  ylab("TIMP-1 logFC in male F2")
ggplotly(all_human_toptables_vs_mouse_plot_m) %>% toWebGL()

all_human_toptables_vs_mouse_plot_f <- all_human_toptables_vs_mouse %>% 
  filter(adj.P.Val < 0.05) %>%
  ggplot(aes(x = logFC, color = sigclass)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(aes(y = logFC.f, gene_name = gene_name)) +
  facet_grid(~dataset) +
  scale_color_manual(values = sigclass_colors) +
  theme_cowplot()+
  xlab("logFC in human NAS>=4 vs NAS<4") +
  ylab("TIMP-1 logFC in female F2")
ggplotly(all_human_toptables_vs_mouse_plot_f)


plot_grid(all_human_toptables_vs_mouse_plot_m + theme(legend.position = "none"),
          all_human_toptables_vs_mouse_plot_f + theme(legend.position = "none"), nrow = 2, align = "hv")
ggsave(filename = "./Plots/Liver/Giorgia/GB02-liver_dge_timp1_vs_human_NAS.pdf", width = 8, height = 8)

# # genes only specifically related to timp-1
# vennDiagram(decideTests(fit_timp))
# timp_spec <- decideTests(fit_timp) %>% as.data.frame() %>% dplyr::filter(`log(TIMP1)` != 0) %>% dplyr::filter(Sexm == 0 & Liver_g_RelBWSac_24 ==0 & Fat_perc_22 == 0 & `(Intercept)` == 0)
# dim(timp_spec)
# timp_spec$gene_id <- rownames(timp_spec)
# timp_spec <- left_join(timp_spec, liver_dge$geneConversionTable) %>% arrange(-`log(TIMP1)`)
# timp_spec_genes_up <- timp_spec %>% filter(`log(TIMP1)` > 0) %>% pull(gene_id)
# ego <- clusterProfiler::enrichGO(gene = timp_spec_genes_up, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", universe = rownames(liver_dge$cpm))
# dotplot(ego)
# ego@result
# 
# liver_dge$geneConversionTable %>% filter(gene_name == "Col3a1") %>% pull(gene_id)
# par(mfrow = c(1,2))
# plot(liver_dge$y$E[liver_dge$geneConversionTable %>% filter(gene_name == "Col3a1") %>% pull(gene_id),], liver_dge$y$E[liver_dge$geneConversionTable %>% filter(gene_name == "Timp1") %>% pull(gene_id),], xlab = "Col3a1 expression", ylab = "Timp1 expression")
# # plot(liver_dge$y$E[liver_dge$geneConversionTable %>% filter(gene_name == "Col3a1") %>% pull(gene_id),], liver_dge$y$E[liver_dge$geneConversionTable %>% filter(gene_name == "Timp1") %>% pull(gene_id),])
# plot(liver_dge$y$E[liver_dge$geneConversionTable %>% filter(gene_name == "Timp1") %>% pull(gene_id),], log10(liver_dge$samples$TIMP1), xlab= "Timp1 expression", ylab="Plasma Timp-1")
# dev.off()


################### plot genes ##########################
selected_gene_names <-  c("Clstn1","Pdgfra", "Vim","Vwf" ,"Jag1","Trem2","Gpnmb","Cd9","Ccl2","Col1a1", "Timp1","Postn")
selected_phenotypes <- c("BW_22","BW_perGain_22","Fat_g_22","Fat_perc_22","Liver_g_24","Liver_g_RelBWSac_24",
                         "ASAT_24","ALAT_24","Cholesterol_24","TIMP1_24","GDF15_24")

liver_phenos <- liver_dge$samples %>%
  dplyr::select(all_of(selected_phenotypes))

dge = liver_dge
gid <- dge$geneConversionTable[match(selected_gene_names,  dge$geneConversionTable$gene_name),"gene_id"]
dtplot <- cbind(liver_phenos, t(liver_dge$y$E[gid,]))
data<- as.matrix(dtplot)
#library(Hmisc)
#corr_data <- Hmisc::rcorr(data, type="spearman")
cormat <- round(cor(data,use = "complete.obs"),2)

library(corrplot)
corrplot(cormat, type = "upper"
         #, order = "hclust", 
         #tl.col = "black", 
         #tl.srt = 45
)

library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri

melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

plotGeneScatter <- function(dge = liver_dge, gname = "Timp1", facet = c("Sex"),
                            x.axis = "TIMP1_24", color = "Strain", xscale = "none", 
                            topt = timp_sex_specific_tt, fit = fit_timp){
  gid <- dge$geneConversionTable[match(gname,  dge$geneConversionTable$gene_name),"gene_id"]
  dtplot <- cbind(dge$samples, "E" = dge$y$E[gid,])
  curtt <- topt %>% dplyr::filter(gene_name == gname)
  
  g <- ggplot(dtplot, aes_string(x = x.axis, y = "E", color = color)) +
    geom_point() +
    stat_cor(method = "pearson", size=2)+
    # geom_abline(slope = lfc, intercept = interc, color = "black") +
    ylab("log2(CPM)") +
    ggtitle(gname, subtitle = paste0("Female logFC = ", round(curtt$logFC.f, 4), "; p.adjust = ", signif(curtt$adj.P.Val.f),"\n", "Male logFC = ", round(curtt$logFC.m, 4), "; p.adjust = ", signif(curtt$adj.P.Val.m)) )
  
  if(color == "Strain"){
    g <- g + scale_color_manual(values = strainColors)
  }
  if(facet!="none"){
    frm <- as.formula(paste0(".~", facet))
    g <- g + facet_grid(frm)
  }
  if(xscale == "log"){
    g <- g + scale_x_log10() + xlab(paste0(x.axis, " (log)"))
  }
  g + theme_cowplot()
}
plotGeneScatter(xscale = "none") + geom_smooth(aes(group = 1), method = "lm")

plist <- lapply(rev(selected_gene_names), 
                function(x){plotGeneScatter(gname = x, xscale = "log")+ 
                    geom_smooth(aes(group = 1), method = "lm") + 
                    theme(legend.position = "none")})

pg <- plot_grid(plotlist = plist)
pgl <- plot_grid(pg, get_legend(plotGeneScatter(xscale = "log") + 
                                  theme(legend.position = "bottom", 
                                        legend.direction = "horizontal")),
                 nrow = 2, rel_heights = c(10,1))
pgl
ggsave(plot = pgl, useDingbats = F, filename = "./Plots/Liver/GB02-selected_genes_vs_TIMP1_24.pdf", width = 18, height = 12)

####################### now plot some of the top up and down genes based on DE ###########################
topgenes <- ttlf %>% split(., sign(ttlf$logFC)) %>% lapply(., function(x){head(x,n = 5) %>% pull(gene_name)}) %>% unlist() %>% factor(., levels = .)

plist <- lapply(topgenes, function(x){plotGeneScatter(gname = x, xscale = "log")+ geom_smooth(aes(group = 1), method = "lm") + theme(legend.position = "none")})
pg <- plot_grid(plotlist = plist, nrow = 2)
pgl <- plot_grid(pg, get_legend(plotGeneScatter(xscale = "log") + theme(legend.position = "bottom", legend.direction = "horizontal")),
                 nrow = 2, rel_heights = c(10,1))
pgl



############## some genes to please johan #####################
## TM6SF2, GCKR, MBOAT7, HSD17B13, IL28B, MERTK –
selected_gene_names_johan = c("Acmsd", "Cox7a2l", "Acaca", "Acacb", "Pnpla3", "Tm6sf2", "Gckr", "Mboat7", "Hsd17b13", "Mertk")

plist <- lapply(selected_gene_names_johan, function(x){
  curtt <- tt%>%filter(gene_name == x)
  plotGeneScatter(gname = x, xscale = "log") + 
    geom_smooth(aes(group = 1), method = "lm") +
    theme(legend.position = "none")
})

pg <- plot_grid(plotlist = plist, nrow = 4)
pgl <- plot_grid(pg, get_legend(plotGeneScatter(xscale = "log") + theme(legend.position = "bottom", legend.direction = "horizontal")),
                 nrow = 2, rel_heights = c(10,1))
pgl
ggsave(useDingbats = F, filename = "./Plots/Liver/02-selected_genes_johan_vs_TIMP1_24.pdf", width = 12, height = 12)

## some volcano plots





## now plot some sex interaction genes (timp-1 by sex)

dplyr::left_join(tt %>% filter(adj.P.Val < 0.05), tt_sexint %>% filter(adj.P.Val < 0.05), by = "gene_id") %>% filter(!is.na(logFC.y)) %>% 
  
  sexIntGenes <- tt_sexint %>% split(., sign(tt_sexint$logFC)) %>% lapply(., function(x){head(x,n = 5) %>% pull(gene_name)}) %>% unlist() %>% factor(., levels = .)

plist <- lapply(sexIntGenes, function(x){plotGeneScatter(gname = x, xscale = "log")+ geom_smooth(aes(group = 1), method = "lm") + theme(legend.position = "none")})
pg <- plot_grid(plotlist = plist, nrow = 2)
pgl <- plot_grid(pg, get_legend(plotGeneScatter(xscale = "log") + theme(legend.position = "bottom", legend.direction = "horizontal")),
                 nrow = 2, rel_heights = c(10,1))
pgl

####################################################################################################
########################### below just some drafts and tests #######################################
####################################################################################################


## biological age prediction, as in https://github.com/emcgl/trap/blob/master/backend/FORMULA-GENERAL-PREDICTOR-GENE_ID.txt ; paper: https://www.nature.com/articles/ncomms9570
gene_predictor <- read.table("./Data/external_data/FORMULA-GENERAL-PREDICTOR-GENE_ID.txt", sep = "]")
gene_predictor <- gene_predictor %>% t()
gene_predictor <- strsplit(gene_predictor[,1], split = "*selection[,", fixed = T)
gene_predictor <- do.call(rbind, gene_predictor)

gene_predictor %>% head()

colnames(gene_predictor) <- c("coef", "gene_symbol_human")
gene_predictor <- as.data.frame(gene_predictor, stringsAsFactors= F)
gene_predictor <- gene_predictor %>% filter(!is.na(gene_symbol_human))
source("./R/mouse_to_human_orthodb.R")
conversion <- convert_human_mouse_orthodb(gene_predictor$gene_symbol_human)
head(conversion)

gene_predictor <- left_join(gene_predictor, conversion, by = "gene_symbol_human")
gene_predictor$coef <- as.numeric(gene_predictor$coef)
gene_predictor$gene_name <- sapply(gene_predictor$gene_name, function(x){if(length(x) == 0){return(NA)}else{x[1]}}) %>% unlist()
gene_predictor$gene_id <- liver_dge$geneConversionTable$gene_id[match(gene_predictor$gene_name,  liver_dge$geneConversionTable$gene_name)]

cpm_df_selected_with_coef <- left_join(gene_predictor, liver_dge$cpm %>% as.data.frame() %>% rownames_to_column(var = "gene_id"))
biological_age_cpm <- sapply(colnames(liver_dge$cpm), function(x){
  (cpm_df_selected_with_coef$coef *cpm_df_selected_with_coef[,x]) %>% sum(na.rm = T)
})
liver_dge$samples$BiologicalAge_cpm <- biological_age_cpm[liver_dge$samples$Sample_ID]
cor.test(liver_dge$samples$TIMP1_24, liver_dge$samples$BiologicalAge, use = "pairwise.complete")
lm(BiologicalAge ~ TIMP1_24 + Sex + TIMP1_24:Sex, data = liver_dge$samples %>% filter(Generation == "F2")) %>% summary()
lm(BiologicalAge ~ TIMP1_24, data = liver_dge$samples %>% filter(Generation == "F2") %>% filter(Sex == "m")) %>% summary()
lm(BiologicalAge ~ TIMP1_24, data = liver_dge$samples %>% filter(Generation == "F2") %>% filter(Sex == "f")) %>% summary()

lm(BiologicalAge ~ TIMP1_24, data = liver_dge$samples %>% filter(Sex == "m")) %>% summary()
lm(BiologicalAge ~ TIMP1_24, data = liver_dge$samples %>% filter(Sex == "f")) %>% summary()

liver_dge$samples %>% ggplot(aes(x = TIMP1_24, y = BiologicalAge_cpm, color = Strain)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1)) +
  facet_grid(~Sex) +
  scale_color_manual(values = strainColors) +
  scale_x_log10() +
  theme_cowplot()

ggsave(useDingbats = FALSE, filename = "./Plots/Liver/02-liver_biological_age_vs_timp1.pdf", width = 12, height = 5)


liver_dge$samples %>% ggplot(aes(x = GDF15_24, y = BiologicalAge_cpm, color = Strain)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1)) +
  facet_grid(~Sex) +
  scale_color_manual(values = strainColors) +
  scale_x_log10() +
  theme_cowplot()

ggsave(useDingbats = FALSE, filename = "./Plots/Liver/02-liver_biological_age_vs_GDF15.pdf", width = 12, height = 5)


## do the same with Timp1 expression in liver

liver_dge$samples %>% mutate(Timp1 = liver_dge$y$E["ENSMUSG00000001131",]) %>% ggplot(aes(x = Timp1, y = BiologicalAge_cpm, color = Strain)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = 1)) +
  facet_grid(~Sex) +
  scale_color_manual(values = strainColors) +
  theme_cowplot()
ggsave(useDingbats = FALSE, filename = "./Plots/Liver/02-liver_biological_age_vs_timp1_liver_expression.pdf", width = 12, height = 5)
# # now apply to rpkm matrix
# txdb <- loadDb("./Data/Mus_musculus.GRCm38.92.txdb")
# 
# ebg <- exonsBy(txdb, by="gene")
# exonicParts(txdb)
# library(parallel)
# 
# gene_lengths <- mclapply(X=ebg, mc.cores = 16, FUN =  function(x){
#   width(reduce(x)) %>% sum()
# })
# gene_lengths <- unlist(gene_lengths)
# rpkm_matrix <- rpkm(liver_dge, gene.length = gene_lengths[match(rownames(liver_dge$counts), names(gene_lengths))], log = T)
# rpkm_df <- rpkm_matrix %>% as.data.frame()
# rpkm_df_selected <- rpkm_df[match(gene_predictor$gene_id, rownames(rpkm_df)),]
# rpkm_df_selected[1:5,1:5]
# 
# rpkm_df_selected_with_coef <- left_join(gene_predictor, rpkm_df_selected %>% rownames_to_column(var = "gene_id"))
# biological_age_rpkm <- sapply(colnames(rpkm_df_selected), function(x){
#   (rpkm_df_selected_with_coef$coef *rpkm_df_selected_with_coef[,x]) %>% sum(na.rm = T)
# })
# liver_dge$samples$BiologicalAge_rpkm <- biological_age_rpkm[liver_dge$samples$Sample_ID]
# liver_dge$samples %>% ggplot(aes(x = log(TIMP1_24), y = BiologicalAge_rpkm, color = Strain)) +
#   geom_point() +
#   geom_smooth(method = "lm", aes(group = 1)) +
#   facet_grid(~Sex) +
#   scale_color_manual(values = strainColors) +
#   theme_cowplot()


library(WGCNA)
cpv <- bicorAndPvalue(liver_dge$samples %>% mutate(sex_male = Sex %>% factor(., levels = c("f", "m")) %>% as.numeric() -1 ))

cpv <- do.call(cbind, lapply(cpv, reshape::melt))[,c(1:3, seq(6,length(cpv)*3,3))] #cpv$bicor %>% reshape::melt() %>% cbind(., reshape::melt(cpv$p))
colnames(cpv) <- gsub("\\.value", "", gsub("bicor.", "", colnames(cpv)))
cpv <- cpv %>% filter(!is.na(value)) %>% filter(X1 != X2) %>% filter(X1 != "Sample_ID") %>% filter(X2 != "Sample_ID")                 

cpv <- cpv %>% filter(! (grepl("BW_", X1)&(grepl("BW_", X2)) ) )
cpv %>% ggplot(aes(x = X1, y=X2, color = value, size = -log10(p))) + geom_point() +
  scale_color_gradient2() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))

cpv %>% filter(X1 %in% c("TIMP1_24", "GDF15_24", "KIM1_24", "ALAT_24", "ASAT","TGL", "Free_Fatty_Acids")) %>% ggplot(aes(x = X1, y=X2, color = value, size = -log10(p))) + geom_point() +
  scale_color_gradient2() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))

cpv %>% filter(X1 %in% c("Kidney_g_24", "Kidney_g_RelBWSac_24")) %>% filter(!grepl("BW", X2)) %>% ggplot(aes(x = X1, y=X2, color = value, size = -log10(p))) + geom_point() +
  scale_color_gradient2() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))



pheno_melt <- read.csv("~/rcp_storage/janssen/aim2_step2/aim2_step2_phenotypes/Data/allData_melt_clean_noOutlier.csv")

rbind(liver_dge$samples %>% mutate(week = 6, acr = ACR_6),
      liver_dge$samples %>% mutate(week = 14, acr = ACR_14),
      liver_dge$samples %>% mutate(week = 22, acr = ACR_22)) %>% 
  filter(week != 6) %>%
  ggplot(aes(x = Kidney_g_RelBWSac_24, y = log(acr+1), color = Strain, group = 1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(Sex~week) +
  scale_color_manual(values = strainColors) +
  theme_cowplot()


rbind(liver_dge$samples %>% mutate(week = 6, acr = ACR_6),
      liver_dge$samples %>% mutate(week = 14, acr = ACR_14),
      liver_dge$samples %>% mutate(week = 22, acr = ACR_22)) %>% 
  # filter(week != 6) %>%
  ggplot(aes(x = Kidney_g_RelBWSac_24, y = log(acr+1), color = Sex, lty=factor(week) )) +
  # geom_point() +
  geom_smooth(method = "lm") +
  # facet_wrap() +
  scale_color_manual(values = sexColors) +
  theme_cowplot()

rbind(liver_dge$samples %>% mutate(week = 6, acr = ACR_6),
      liver_dge$samples %>% mutate(week = 14, acr = ACR_14),
      liver_dge$samples %>% mutate(week = 22, acr = ACR_22)) %>% 
  filter(week != 6) %>%
  ggplot(aes(x = Kidney_g_RelBWSac_24, y = log(acr+1), color = factor(week), lty=factor(week) )) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Sex) +
  # scale_color_manual(values = sexColors) +
  # scale_linetype_manual(values = c(1,5,10))
  theme_cowplot()



rbind(liver_dge$samples %>% mutate(week = 6, acr = ACR_6),
      liver_dge$samples %>% mutate(week = 14, acr = ACR_14),
      liver_dge$samples %>% mutate(week = 22, acr = ACR_22)) %>% 
  filter(week != 6) %>%
  ggplot(aes(x = week, y = acr, color = Strain, group = ID)) +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "lm") +
  facet_wrap(Sex~.) +
  scale_color_manual(values = strainColors) +
  theme_cowplot()


rbind(liver_dge$samples %>% mutate(week = 6, acr = create.post()),
      liver_dge$samples %>% mutate(week = 14, acr = ACR_14),
      liver_dge$samples %>% mutate(week = 22, acr = ACR_22)) %>% 
  filter(week != 6) %>%
  ggplot(aes(x = week, y = acr, color = Strain, group = ID)) +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "lm") +
  facet_wrap(Sex~.) +
  scale_color_manual(values = strainColors) +
  theme_cowplot()


liver_dge$samples %>% mutate(FEUrea = (Urine_urea_22*Crea_24)/(Urine_creatinin_22*Urea_24)) %>%
  filter(FEUrea < 1000) %>%
  ggplot(aes(x = Kidney_g_RelBWSac_24, y = FEUrea, color = Strain)) + 
  scale_color_manual(values = strainColors) +
  geom_point() +
  facet_wrap(Sex~.) +
  theme_cowplot()



# # cohort info
# 
# bw <- readxl::read_xlsx("~/JanssenArchive2/aim2_step2/phenotype_data/A_BW/Aim2Step2_BodyWeight.xlsx")
# mouse_metadata <- bw %>% as.data.frame() %>% 
#   select(Cohort, Strain, Diet, Cage, `Animal ID`, Sex, `Date of Birth`, Comments, `Vet case nr`, `Vet case cause`) %>% 
#   unique() %>% 
#   filter(!is.na(Strain)) %>% 
#   mutate(Sample_ID = substr(`Animal ID`, 7, 10))
# 
# pheno_batch <- left_join(liver_dge$samples, mouse_metadata %>% select(-Sex, -Strain), by = "Sample_ID")
# pheno_batch <- pheno_batch %>%  filter(!is.na(Glycemia_22))
# 
# pheno_batch_f2 <- pheno_batch %>% filter(Generation == "F2") %>% filter(!is.na(Glycemia_22))
# library(lme4)
# 
# linear_mod <- lm(Glycemia_22 ~ as.factor(Cohort), data = pheno_batch )
# mlm <-  lmer(Glycemia_22 ~ (1|Cohort), data = pheno_batch )
# 
# pheno_batch$Glycemia_22_residual <- residuals(linear_mod)
# pheno_batch$Glycemia_22_residualmlm <- residuals(mlm)
# 
# 
# a <- ggplot(pheno_batch, aes(x = factor(Cohort), y = Glycemia_22, color = Strain, group = factor(Cohort))) +
#   geom_violin() +
#   geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
#   scale_color_manual(values = strainColors)
# b <- ggplot(pheno_batch, aes(x = factor(Cohort), y = Glycemia_22_residual, color = Strain, group = factor(Cohort))) +
#   geom_violin() +
#   geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
#   scale_color_manual(values = strainColors)
# c <- ggplot(pheno_batch, aes(x = factor(Cohort), y = Glycemia_22_residualmlm, color = Strain, group = factor(Cohort))) +
#   geom_violin() +
#   geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
#   scale_color_manual(values = strainColors)
# plot_grid(a,b,c, ncol = 1)
# 
# 
# ggplot(pheno_batch_f2, aes(x = Glycemia_22, y = Glycemia_22_residual, fill = Sex)) +
#   geom_point()
# 
# ggplot(pheno_batch_f2, aes(x = Glycemia_22, y = Glycemia_22_residualmlm, fill = Sex)) +
#   geom_point()
# 
# ggplot(pheno_batch_f2, aes(x = Glycemia_22_residual, y = Glycemia_22_residualmlm, fill = Sex)) +
#   geom_point()

