######## quantify consensus signaure of TF overexpression phenotype
######## perturbation association analysis
library(ggplot2)
library(cowplot)

tf.pheno.env <- new.env()
load('data/processed/TF_trends.RData', envir = tf.pheno.env)
# upload clone association data
clone_data.env <- new.env()
load(file = '../Association_analysis/data/processed/use/association_summary.RData', envir = clone_data.env)
clone_data_list <- clone_data.env$test.basic.list
#
mutant_data.env <- new.env()
load(file = '../RNA_analysis/data/processed/mutant_de.RData', envir = mutant_data.env)

# select TFs for growth
growth.table <- tf.pheno.env$tf.trends.no25.v1m #tf.pheno.env$tf.trends.v1m
head(growth.table)
tf.grow <- growth.table$TF[!is.na(growth.table$p.t) & p.adjust(growth.table$p.t, 'BH') < 0.1 & growth.table$t > 0]
tf.decrease <- growth.table$TF[!is.na(growth.table$p.t) & p.adjust(growth.table$p.t, 'BH') < 0.1 & growth.table$t < 0]
str(tf.grow)
str(tf.decrease)
#
tf.features$growth <- 'stable'
tf.features$growth[tf.features$TF %in% tf.grow] <- 'grow'
tf.features$growth[tf.features$TF %in% tf.decrease] <- 'decrease'
tf.features$growth <- factor(tf.features$growth, level = c('stable', 'grow', 'decrease'))
table(tf.features$growth)

# select TFs for lineage bias
bias.table <- tf.pheno.env$biastot
head(bias.table)
tf.my <- bias.table$TF[bias.table$trend == 'myeloid']
tf.lymp <- bias.table$TF[bias.table$trend == 'B']
#
tf.features$bias <- 'none'
tf.features$bias[tf.features$TF %in% tf.my] <- 'myeloid'
tf.features$bias[tf.features$TF %in% tf.lymp] <- 'lymphoid'
tf.features$bias <- factor(tf.features$bias, level = c('none', 'myeloid'))
table(tf.features$bias)

### test growth DE genes
# cutoffs as in 6_tf_invitro.de.R
retain.detected <- colSums(tf.clust.mat > 0) > 100
table(retain.detected)
retain.gns <- rowSums(cd.tfvitro > 0) / ncol(cd.tfvitro) > 0.05
table(retain.gns)
dds_tf_growth <- DESeqDataSetFromMatrix(countData = tf.clust.mat[retain.gns, retain.detected],
                                 colData = tf.features[retain.detected,],
                                 design =  ~ growth + cluster)
dds_tf_growth <- DESeq(dds_tf_growth)
resultsNames(dds_tf_growth)
res_growth <- results(dds_tf_growth, name = "growth_grow_vs_stable")
res_decrease <- results(dds_tf_growth, name = "growth_decrease_vs_stable")
head(data.frame(res_growth[order(res_growth$pvalue),]), 20)
head(data.frame(res_decrease[order(res_decrease$pvalue),]), 20)
## result: no significant hits

### test consensus myeloid perturbation signature
dds_tf_myeloid <- DESeqDataSetFromMatrix(countData = tf.clust.mat[retain.gns, retain.detected],
                                        colData = tf.features[retain.detected,],
                                        design =  ~ bias + cluster)
dds_tf_myeloid <- DESeq(dds_tf_myeloid)
resultsNames(dds_tf_myeloid)
res_myeloid <- results(dds_tf_myeloid, name = "bias_myeloid_vs_none")
head(data.frame(res_myeloid[order(res_myeloid$pvalue),]), 10)

saveRDS(list(dds_tf_myeloid = dds_tf_myeloid, res_myeloid = res_myeloid,
             dds_tf_growth = dds_tf_growth, res_growth = res_growth, res_decrease = res_decrease), file = 'data/processed/TF_consensus_signatures')

plot_tf_association <- function(res, pheno_name, clone_data_list = NULL, mut_response = NULL, fdr.cutoff = 0.2, pv.cutoff = NULL, recalculate.fdr = FALSE){
  # volcano DE plot
  if (recalculate.fdr) res$padj <- p.adjust(res$pvalue, 'BH')
  res$gene <- rownames(res)
  if (is.null(pv.cutoff)){
    res$sign <- !is.na(res$padj) & res$padj < fdr.cutoff
  }else{
    res$sign <- !is.na(res$pvalue) & res$pvalue < pv.cutoff
  }
  res$gene[!res$sign] <- ''
  pl1 <- ggplot(res, aes(x = log2FoldChange, y = pvalue, label = gene)) +
    geom_point(color = ifelse(res$sign, 'lightblue', "grey60"), alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = 2) + scale_y_log10() + theme_bw() +
    geom_text_repel(max.overlaps = 5, box.padding = 0.5,  point.padding = 0.3,
                    segment.color = 'grey50', size = 3) + labs(x = "in vitro log2 FC", y = "p-vlaue", title = paste("overexpressed TFs with in vivo", pheno_name)) +
    theme(axis.title = element_text(size = 16), plot.title = element_text(size = 16, hjust = 0.5))
  #pl1

  # show correlation with clone phenotype score
  clone_data <- clone_data_list[[pheno_name]]
  gns.shared <- intersect(rownames(res), rownames(clone_data))
  str(gns.shared)
  if (is.null(pv.cutoff)){
    gns.sign <- gns.shared[!is.na(res[gns.shared,]$padj) & res[gns.shared,]$padj < fdr.cutoff]
  }else{
    gns.sign <- gns.shared[!is.na(res[gns.shared,]$pvalue) & res[gns.shared,]$pvalue < pv.cutoff]
  }
  str(gns.sign)
  res_comparison <- data.frame(cbind(log2FoldChange = res[gns.sign,]$log2FoldChange, pheno = clone_data[gns.sign,]$pheno), gene = gns.sign)
  rownames(res_comparison) <- gns.sign
  if (length(gns.sign) > 2){
  cr <- cor.test(res_comparison$log2FoldChange, res_comparison$pheno)
  stats_label <- paste0("r = ", round(cr$estimate, 2),
                        "\np = ", format(cr$p.value, scientific = TRUE, digits = 2))
  pl2 <- ggplot(res_comparison, aes(x = log2FoldChange, y = pheno, label = gene)) +
    geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) +
    geom_text_repel(max.overlaps = 5, box.padding = 0.5,  point.padding = 0.3, segment.color = 'grey50', size = 3) +
    annotate("text", x = -Inf, y = Inf, label = stats_label,
             hjust = -0.2, vjust = 1.5, size = 5, fontface = "italic") +
    labs(x = paste('in vivo ', pheno_name, ' TFs, in vitro logFC', sep = ''), y = paste("clones with ", pheno_name, ', GE shift', sep = '') ) +
    theme(axis.title = element_text(size = 16), plot.title = element_text(size = 16, hjust = 0.5))
  }else {pl2 <- NULL}

  # show correlations with TET2 an DNMT3A HSC alterations
  gns.shared <- intersect(rownames(res), rownames(mut_response))
  str(gns.shared)
  #gns.sign <- gns.shared[res[gns.shared,]$padj < fdr.cutoff & !is.na(res[gns.shared,]$padj)]
  if (is.null(pv.cutoff)){
    gns.sign <- gns.shared[!is.na(res[gns.shared,]$padj) & res[gns.shared,]$padj < fdr.cutoff]
  }else{
    gns.sign <- gns.shared[!is.na(res[gns.shared,]$pvalue) & res[gns.shared,]$pvalue < pv.cutoff]
  }
  str(gns.sign)
  res_comparison <- data.frame(cbind(log2FoldChange = res[gns.sign,]$log2FoldChange,
                                     DNMT3A = mut_response[gns.sign, 'DNMT3A'], TET2 = mut_response[gns.sign, 'TET2']), gene = gns.sign)
  rownames(res_comparison) <- gns.sign
  if (length(gns.sign) > 2){
  # plot DNMT3A
  cr <- cor.test(res_comparison$log2FoldChange, res_comparison$DNMT3A)
  stats_label <- paste0("r = ", round(cr$estimate, 2), "\np = ", format(cr$p.value, scientific = TRUE, digits = 2))
  pl3 <- ggplot(res_comparison, aes(x = log2FoldChange, y = DNMT3A, label = gene)) +
    geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) +
    geom_text_repel(max.overlaps = 5, box.padding = 0.5,  point.padding = 0.3, segment.color = 'grey50', size = 3) +
    annotate("text", x = -Inf, y = Inf, label = stats_label,
             hjust = -0.2, vjust = 1.5, size = 5, fontface = "italic") +
    labs(x = paste('in vivo ', pheno_name, ' TFs, in vitro logFC', sep = ''), y = "DNMT3A, log FC" ) +
    theme(axis.title = element_text(size = 16), plot.title = element_text(size = 16, hjust = 0.5))

  # plot TET2
  cr <- cor.test(res_comparison$log2FoldChange, res_comparison$TET2)
  stats_label <- paste0("r = ", round(cr$estimate, 2), "\np = ", format(cr$p.value, scientific = TRUE, digits = 2))
  pl4 <- ggplot(res_comparison, aes(x = log2FoldChange, y = TET2, label = gene)) +
    geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) +
    geom_text_repel(max.overlaps = 5, box.padding = 0.5,  point.padding = 0.3, segment.color = 'grey50', size = 3) +
    annotate("text", x = -Inf, y = Inf, label = stats_label,
             hjust = -0.2, vjust = 1.5, size = 5, fontface = "italic") +
    labs(x = paste('in vivo ', pheno_name, ' TFs, in vitro logFC', sep = ''), y = "TET2, log FC" ) +
    theme(axis.title = element_text(size = 16), plot.title = element_text(size = 16, hjust = 0.5))
  }else{
    pl3 <- NULL; pl4 <- NULL
  }

  return(list(pl1 = pl1, pl2 = pl2, pl3 = pl3, pl4 = pl4))
}


## plot myeloid phenotype
res <- res_myeloid
pheno_name <- 'myeloid bias'
plot_res <- plot_tf_association(res, pheno_name, clone_data_list, mut_response)

pdf(paste('figures/revision/Oe_association_', pheno_name, '.pdf', sep = ''),
    width = 10, height = 10)
plot_grid(plot_res$pl1, plot_res$pl2, plot_res$pl3, plot_res$pl4,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, scale = 0.9)
dev.off()

## plot growth phenotype
res <- res_growth
pheno_name <- 'growth rate'
fdr.cutoff <- 0.2
pv.cutoff <- 0.01
plot_res <- plot_tf_association(res, pheno_name, clone_data_list, mut_response,
                                fdr.cutoff = fdr.cutoff, pv.cutoff = pv.cutoff)

pdf(paste('figures/revision/Oe_association_', pheno_name, '.pdf', sep = ''),
    width = 10, height = 10)
plot_grid(plot_res$pl1, plot_res$pl2, plot_res$pl3, plot_res$pl4,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, scale = 0.9)
dev.off()


### pathway-centric analysis
# how to choose a representative set of pathways/gene sets to visuzliae?

tbl_dnmt3a <- mutant_data.env$de.mut2wt.shrink[[1]]$DNMT3A
tbl_tet2 <- mutant_data.env$de.mut2wt.shrink[[1]]$TET2
mutant_de_list <- list(
  DNMT3A_up = rownames(tbl_dnmt3a)[tbl_dnmt3a$padj < 0.2 & tbl_dnmt3a$log2FoldChange > 0.1 & !is.na(tbl_dnmt3a$padj)],
  DNMT3A_down = rownames(tbl_dnmt3a)[tbl_dnmt3a$padj < 0.2 & tbl_dnmt3a$log2FoldChange < -0.1 & !is.na(tbl_dnmt3a$padj)],
  TET2_up = rownames(tbl_tet2)[tbl_tet2$padj < 0.2 & tbl_tet2$log2FoldChange > 0.1 & !is.na(tbl_tet2$padj)],
  TET2_down = rownames(tbl_tet2)[tbl_tet2$padj < 0.2 & tbl_tet2$log2FoldChange < -0.1 & !is.na(tbl_tet2$padj)]
)
str(mutant_de_list)
fgsea_tf_mut <- apply(logFC.tf, 2, function(x){
  sumr <- fgsea(pathways = mutant_de_list, stats  = x, minSize = 5, maxSize = 1000)
})

plot_gene_set_volcano <- function(gsea.dat, gene_set_select, tf.grow, tf.decrease, tf.my, fdr.cutoff = 0.05){
  tf_enrich_scores <- do.call(rbind, lapply(1:length(gsea.dat), function(i){
    gsea.dat[[i]][grepl(gene_set_select, gsea.dat[[1]]$pathway), ]
  }))
  tf_enrich_scores <- data.frame(tf_enrich_scores)
  rownames(tf_enrich_scores) <- names(gsea.dat)
  tf_enrich_scores$TF <- rownames(tf_enrich_scores)
  #
  tf_enrich_scores$growth <- NA
  tf_enrich_scores$growth[tf_enrich_scores$TF %in% tf.grow] <- 'up'
  tf_enrich_scores$growth[tf_enrich_scores$TF %in% tf.decrease] <- 'down'
  #
  tf_enrich_scores$myeloid <- NA
  tf_enrich_scores$myeloid[tf_enrich_scores$TF %in% tf.my] <- 'myeloid'
  #
  #
  tf_enrich_scores$sign <- tf_enrich_scores$padj < fdr.cutoff

  pl <- ggplot(tf_enrich_scores, aes(x = NES, y = pval, label = TF)) +
    geom_point(aes(shape = growth, fill = myeloid), color = "black",
               alpha = ifelse(tf_enrich_scores$sign, 1, 0.2), size = 5) +
    #scale_shape_manual(values = c("up" = 17, "down" = 18), na.value = 19) +
    scale_shape_manual(values = c("up" = 24, "down" = 23, "NA" = 21),
                       na.value = 21) +
    scale_fill_manual(values = c("myeloid" = "orange", "NA" = "grey60"), name = "Bias") +
    geom_vline(xintercept = 0, linetype = 2) + scale_y_log10() + theme_bw() +
    geom_text_repel(max.overlaps = 10, box.padding = 0.5,  point.padding = 0.3, segment.color = 'grey50', size = 4) +
    labs(x = "Normalized Enrichment Score ", y = "p-vlaue", title = gene_set_select) +
    theme(axis.title = element_text(size = 16), plot.title = element_text(size = 16, hjust = 0.5))

  return(pl)
}

gsea.dat <- fgsea_tf_mut
pl <- list()
for(i in 1:length(gsea.dat[[1]]$pathway)){
  gene_set_select <- fgsea_tf_mut[[1]]$pathway[i]
  pl.temp <- plot_gene_set_volcano(gsea.dat, gene_set_select, tf.grow, tf.decrease, tf.my, fdr.cutoff = 0.05)
  pl <- append(pl, pl.temp)
}
#
pdf(paste('figures/revision/Mutant_cor_tfs.pdf', sep = ''), width = 10, height = 8)
plot_grid(pl[[1]], pl[[2]], pl[[3]], pl[[4]],
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, scale = 0.9)
dev.off()


gsea.dat <- gsea.tfover.compose
pl <- list()
path.show <- c('MYC_TARGETS_V1', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'Signaling by Receptor Tyrosine Kinases')
for(gene_set_select in path.show){
  pl.temp <- plot_gene_set_volcano(gsea.dat, gene_set_select, tf.grow, tf.decrease, tf.my, fdr.cutoff = 0.05)
  pl <- append(pl, pl.temp)
}
#
pdf(paste('figures/revision/Pathways_cor_tfs.pdf', sep = ''), width = 10, height = 8)
plot_grid(pl[[1]], pl[[2]], pl[[3]], pl[[4]],
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, scale = 0.9)
dev.off()
