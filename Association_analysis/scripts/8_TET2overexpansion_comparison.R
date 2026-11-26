library(ggplot2)
library(ggrepel)

tet2_growth_program <- readxl::read_excel('../RNA_analysis/data/genesets/download/blood_bld-2023-019620-mmc2.xlsx', )
tet2_growth_program <- data.frame(tet2_growth_program)
colnames(tet2_growth_program) <- c('gene', 'logFC', 'pvalue', 'adjP')
qplot(tet2_growth_program$logFC, tet2_growth_program$pvalue, log = 'y')
tet2_growth_program$geneH <- toupper(tet2_growth_program$gene)
rownames(tet2_growth_program) <- tet2_growth_program$geneH
gns.sign <- tet2_growth_program$gene[tet2_growth_program$adjP < 0.2]
str(gns.sign)
tet2_growth_program$sign <- tet2_growth_program$gene %in% gns.sign
tet2_growth_program$gene_show <- ifelse(tet2_growth_program$sign, tet2_growth_program$gene[tet2_growth_program$sign], NA)
pl <- ggplot(tet2_growth_program, aes(logFC, pvalue, col = sign, label = gene_show)) + geom_point(size = 2) + scale_y_log10() +
  geom_text_repel() + theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = c(1-0.03, 0.03), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
pl
pdf('figures/default/TET2_overexpansion_scatterplot.pdf', width = 7, height = 7)
pl
dev.off()


#########################. compare with  gene modules of clone phenotypes ################
ass.env <- new.env()
load('data/processed/default/association_summary.RData', envir = ass.env)

tbl <- do.call(rbind, lapply(1:length(ass.env$test.basic.list), function(i){
  tbl <- data.frame(ass.env$test.basic.list[[i]])
  tbl$phenotype <- names(ass.env$test.basic.list)[i]
  tbl$gene <- rownames(tbl)
  tbl
}))
tbl$score <- tbl$pheno
tbl$value <- ifelse(tbl$gene %in% gns.sign, 1, 0)
tmp1 <- tbl

# function to prepare data structure for roc curve visualization
prepare_auc_params <- function(tmp1){
  require(dplyr)
  df <- tmp1 %>%
    filter(!is.na(value)) %>%
    mutate(value = factor(value))
  table(df$value, df$phenotype)
  require(pROC)
  roc_list <- df %>%
    group_by(phenotype) %>%
    group_map(~{
      roc_obj <- roc(.x$value, .x$score, quiet = TRUE)
      list(
        phenotype = .y$phenotype,   # ✅ correct way
        roc = roc_obj,
        auc = auc(roc_obj)
      )
    })
  auc_df <- do.call(rbind, lapply(roc_list, function(x) {
    data.frame(phenotype = x$phenotype, AUC = as.numeric(x$auc))
  }))
  #plot
  roc_plot_df <- do.call(rbind, lapply(roc_list, function(x) {
    data.frame(
      specificity = x$roc$specificities,
      sensitivity = x$roc$sensitivities,
      phenotype = x$phenotype
    )
  }))
  #
  auc_labels <- paste0(auc_df$phenotype, " (AUC=", round(auc_df$AUC, 2), ")")
  return(list(roc_plot_df = roc_plot_df, auc_labels = auc_labels))
}


roc_stat <- prepare_auc_params(tmp1)
pl <- ggplot(roc_stat$roc_plot_df, aes(x = 1 - specificity, y = sensitivity, color = phenotype)) +
  #geom_point(size = 0.5) +
  geom_path(size = 3) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  scale_color_discrete(labels = roc_stat$auc_labels) +
  labs(title = paste('TET2 KO overexpansion program'), color = "HSC gene sets") +
  theme(text = element_text(size = 20),
        legend.position = c(0.02, 1-0.02), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
pl


pdf('figures/default/roc_curve_TET2_overexpansion_phenotypes.pdf', width = 7, height = 7)
pl
dev.off()

de.env <- new.env()
load('../RNA_analysis/data/processed//mutant_de.RData', envir = de.env)
tbl <- do.call(rbind, lapply(names(de.env$de.mut2wt.shrink[[1]]), function(KO){
  tbl <- de.env$de.mut2wt.shrink[[1]][[KO]]
  tbl$gene <- rownames(tbl)
  tbl$phenotype <- KO
  tbl
}))
tbl <- data.frame(tbl)
head(tbl)
tbl$score <- tbl$log2FoldChange
tbl$value <- ifelse(tbl$gene %in% gns.shared, 1, 0)
tmp1 <- tbl
#
roc_stat <- prepare_auc_params(tmp1)
pl <- ggplot(roc_stat$roc_plot_df, aes(x = 1 - specificity, y = sensitivity, color = phenotype)) +
  #geom_point(size = 0.5) +
  geom_path(size = 3) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  scale_color_discrete(labels = roc_stat$auc_labels) +
  labs(title = paste('TET2 KO overexpansion program'), color = "HSC gene sets") +
  theme(text = element_text(size = 20),
        legend.position = c(0.02, 1-0.02), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
pl

pdf('figures/default/roc_curve_TET2_overexpansion_mutants.pdf', width = 7, height = 7)
pl
dev.off()


