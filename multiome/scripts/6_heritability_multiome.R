clone_multiome_env <- new.env()
load(file = './../RNA_clone_analysis/data/processed/clone_multiome_info.RData', envir = clone_multiome_env)
#
#hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
#de.env <- new.env()
#load(file = '../RNA_analysis/data/processed/mutant_de.RData', envir = de.env)
#env.dna <- new.env()
#load('../DNA_analysis/data/processed/DNA_clones.RData', env = env.dna)
obj.combined <- readRDS(file = 'data/processed/multiome_seurat.rds')


### create cell <- clone -> phenotype matrix
cell2clone <- clone.states2[clone_multiome_env$clone.info.multi[clone_multiome_env$cell.clone.multi$cloneid,]$larry,]
cell2clone$genotype <- clone.info[clone_multiome_env$clone.info.multi[clone_multiome_env$cell.clone.multi$cloneid,]$larry,]$genotype
rownames(cell2clone) <- rownames(clone_multiome_env$cell.clone.multi)
str(cell2clone)
#
projection_clusters <- obj.combined$projection
names(projection_clusters) <- unlist(lapply(strsplit(names(projection_clusters), '-'), function(x) x[1]))
#
#cells.in <- unlist(lapply(strsplit(rownames(obj.combined@reductions$pca@cell.embeddings), '-'), function(x)x[1]))
#cells.shared <- intersect(rownames(cell2clone), cells.in)
cells.shared <- rownames(cell2clone)
cells.shared <- cells.shared[!is.na(cell2clone[cells.shared,]$genotype) & !is.na(cell2clone[cells.shared,9])]
str(cells.shared)


r2.cv <- function(ipheno = 9, cluster_state = 1,genotype = 'WT', reduc = reduc, n.samples = 100, kfrac = 0.8){
  r2vec <- unlist(lapply(1:n.samples, function(j){
    cells.all <- cells.shared[cell2clone[cells.shared,]$genotype == genotype]
    cells.all <- cells.all[projection_clusters[cells.all] %in% cluster_state]
    cells.train <- sample(cells.all, size = length(cells.all) * kfrac)
    cells.test <- setdiff(cells.all, cells.train)
    train.df <- data.frame(
      y = cell2clone[cells.train, ipheno], reduc[cells.train, ]#, clust = projection_clusters[cells.train]
    )
    ft <- lm(y ~ ., data = train.df)
    summary(ft)
    #
    test.df <- data.frame(reduc[cells.test, ])
    y.pred <- predict(ft, newdata = test.df)
    y.true <- cell2clone[cells.test, ipheno]
    R2 <- 1 - sum((y.true - y.pred)^2) / sum((y.true - mean(y.true))^2)
    R2 <- cor(y.true, y.pred, method = 'spearman')
    return(R2)
  }))
  return(r2vec)
}


modalities <- c('RNA', 'ATAC')
genotypes <- c('WT', 'TET2', 'DNMT3A')
pheno.clusts <- c(4, 8, 9, 10)
heritability.spearman.cv <- do.call(rbind,lapply(modalities, function(modality){
  print(modality)
  if (modality == 'RNA'){
    reduc <- obj.combined@reductions$pca@cell.embeddings
  }else{
    reduc <- obj.combined@reductions$lsi@cell.embeddings
  }
  nms <- unlist(lapply(strsplit(rownames(reduc), '-'), function(x) x[1]))
  rownames(reduc) <- nms
  #
  do.call(rbind, lapply(genotypes, function(genotype){
    print(genotype)
    do.call(rbind, lapply(pheno.clusts, function(pheno.clust){
      r2.cv.sample <- r2.cv(ipheno = pheno.clust, cluster_state = c(1), genotype = genotype, reduc = reduc, kfrac = 0.9,  n.samples = 300)
      data.frame(m = mean(r2.cv.sample), sd = sd(r2.cv.sample), pheno = pheno.clust, genotype = genotype, modality = modality)
    }))
  }))
}))

###  TET2 vs WT
tbl <- heritability.spearman.cv[heritability.spearman.cv$genotype %in% c('WT', 'TET2'), ]
require(tidyverse)
require(ggrepel)
tbl_wide <- tbl %>%
  pivot_wider(
    names_from = genotype,
    values_from = c(m, sd),
    names_glue = "{genotype}_{.value}"
)
#
text_size <- 30
new_labels <- c("4" = "Lymph", "8" = "GMP", "9" = 'MkP', "10" = "Ery")
pl <- ggplot(tbl_wide, aes(x = WT_m, y = TET2_m, color = modality)) +
  geom_errorbar(aes(ymin = TET2_m - TET2_sd, ymax = TET2_m + TET2_sd), width = 0) +
  geom_errorbarh(aes(xmin = WT_m - WT_sd, xmax = WT_m + WT_sd), height = 0) +
  geom_point(size = 6) +
  # Optional: Label the points by phenotype
  #geom_text(aes(label = pheno), vjust = -1, show.legend = FALSE) +
  geom_text_repel(aes(label = new_labels[as.character(pheno)]), size = 6,
                  box.padding = 0.5,      # Distance between label and point
                  point.padding = 0.3,    # Distance from the point itself
                  segment.color = 'grey50'
                    ) + #, vjust = -1) +
  # Optional: Add a diagonal line to see which genotype is higher
  geom_abline(slope = 1, intercept = 0, size = 3, linetype = "dashed", color = "grey40") +
  lims(x = c(0, 0.6), y = c(0, 0.6)) +
  theme_bw() +
  labs(
    x = "WT HSC to Phenotype, Spearman",
    y =  "TET2 HSC to Phenotype, Spearman",
    title = "HSC to Phenotype heritability"
  ) + theme(text = element_text(size = text_size), axis.title = element_text(size = text_size), legend.text = element_text(size = text_size), legend.title = element_text(size = text_size),
            plot.title = element_text(size = text_size,  hjust = 0.5),
            legend.position = c(0.05, 0.95), legend.justification = c("left", "top"),
            legend.background = element_rect(fill = "white", color = "black")) # Optional: adds a box)
pl
#
pdf('figures/heritability/scatterplot_WT_vs_TET2.pdf', width = 8, height = 8)
pl
dev.off()


###  DNMT3A vs WT
tbl <- heritability.spearman.cv[heritability.spearman.cv$genotype %in% c('WT', 'DNMT3A'), ]
require(tidyverse)
tbl_wide <- tbl %>%
  pivot_wider(
    names_from = genotype,
    values_from = c(m, sd),
    names_glue = "{genotype}_{.value}"
  )
#
text_size <- 30
new_labels <- c("4" = "Lymph", "8" = "GMP", "9" = 'MkP', "10" = "Ery")
pl <- ggplot(tbl_wide, aes(x = WT_m, y = DNMT3A_m, color = modality)) +
  geom_errorbar(aes(ymin = DNMT3A_m - DNMT3A_sd, ymax = DNMT3A_m + DNMT3A_sd), width = 0) +
  geom_errorbarh(aes(xmin = WT_m - WT_sd, xmax = WT_m + WT_sd), height = 0) +
  geom_point(size = 6) +
  # Optional: Label the points by phenotype
  #geom_text(aes(label = pheno), vjust = -1, show.legend = FALSE) +
  geom_text_repel(aes(label = new_labels[as.character(pheno)]), size = 6,
                  box.padding = 0.5,      # Distance between label and point
                  point.padding = 0.3,    # Distance from the point itself
                  segment.color = 'grey50'
  ) + #, vjust = -1) +
  # Optional: Add a diagonal line to see which genotype is higher
  geom_abline(slope = 1, intercept = 0, size = 3, linetype = "dashed", color = "grey40") +
  lims(x = c(0, 0.6), y = c(0, 0.6)) +
  theme_bw() +
  labs(
    x = "WT HSC to Phenotype, Spearman",
    y =  "DNMT3A HSC to Phenotype, Spearman",
    title = "HSC to Phenotype heritability"
  ) + theme(text = element_text(size = text_size), axis.title = element_text(size = text_size), legend.text = element_text(size = text_size), legend.title = element_text(size = text_size),
            plot.title = element_text(size = text_size,  hjust = 0.5),
            legend.position = c(0.05, 0.95), legend.justification = c("left", "top"),
            legend.background = element_rect(fill = "white", color = "black")) # Optional: adds a box)
pl
#
pdf('figures/heritability/scatterplot_WT_vs_DNMT3A.pdf', width = 8, height = 8)
pl
dev.off()
