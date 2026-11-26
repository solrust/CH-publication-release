run <- 'RNA'

projection_clusters <- clust.final
cell2clone <- cell.clone
cell2clone$genotype <- clone.info[cell2clone$cloneid, ]$genotype
cell2clone <- data.frame(cell2clone, clone.pheno[cell2clone$cloneid,])
#
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
modalities <- 'scRNA'
genotypes <- c('WT', 'TET2', 'DNMT3A')
pheno.clusts <- c(4, 8, 9, 10)
if (modalities == 'scRNA'){
  pheno.clusts <- which(colnames(cell2clone) %in% c('growth.rate', 'lymphoid.bias', 'myeloid.bias', 'erythroid.bias'))
}
heritability.spearman.cv <- do.call(rbind,lapply(modalities, function(modality){
  print(modality)
  if (modality == 'RNA'){
    reduc <- obj.combined@reductions$pca@cell.embeddings
  }else if (modality == 'ATAC'){
    reduc <- obj.combined@reductions$lsi@cell.embeddings
  }else{
    reduc <- p2pca
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

###  CH vs WT
tbl <- heritability.spearman.cv[heritability.spearman.cv$genotype %in% c('WT', 'DNMT3A', 'TET2'), ]
require(tidyverse)
tbl_wide <- tbl %>%
  pivot_wider(
    names_from = genotype,
    values_from = c(m, sd),
    names_glue = "{genotype}_{.value}"
  )
#
text_size <- 30
#new_labels <- c("4" = "Lymph", "8" = "GMP", "9" = 'MkP', "10" = "Ery")
if (modality == 'scRNA'){
  new_labels <- c("11" = "Lymph", "12" = "Dn", "13" = 'GMP', "14" = "Ery")
  tbl_wide <- tbl_wide[!tbl_wide$pheno %in% '12',]
}

 pl <- ggplot(data = tbl_wide, aes(x = WT_m)) +
  geom_point(aes(y = DNMT3A_m,  color = "DNMT3A"), size = 6) +
  geom_point(aes(y = TET2_m,  color = "TET2"), size = 6) +
  geom_errorbar(aes(y = DNMT3A_m, xmin = WT_m - WT_sd, xmax = WT_m + WT_sd, color = "DNMT3A"), height = 0) +
  geom_errorbar(aes(y = TET2_m, xmin = WT_m - WT_sd, xmax = WT_m + WT_sd, color = "TET2"), height = 0) +
  geom_errorbar(aes(x = WT_m, ymin = TET2_m - TET2_sd, ymax = TET2_m + TET2_sd, color = 'TET2'), width = 0) +
  geom_errorbar(aes(x = WT_m, ymin = DNMT3A_m - DNMT3A_sd, ymax = DNMT3A_m + DNMT3A_sd, color = 'DNMT3A'), width = 0) +
  geom_abline(slope = 1, intercept = 0, size = 3, linetype = "dashed", color = "grey40") +
  geom_text_repel(aes(y = TET2_m, label = new_labels[as.character(pheno)],, color = "TET2"),
                  size = 6, box.padding = 0.5, point.padding = 0.3, segment.color = 'grey50') +
  geom_text_repel(aes(y = DNMT3A_m, label = new_labels[as.character(pheno)],, color = "DNMT3A"),
                  size = 6, box.padding = 0.5, point.padding = 0.3, segment.color = 'grey50') +
  lims(x = c(0, 0.5), y = c(0, 0.5)) +
  theme_bw() +
  labs(x = "WT HSC to Phenotype, Spearman",
      y =  "CH HSC to Phenotype, Spearman",
      title = "HSC to Phenotype heritability",
  ) +
   scale_color_manual(
     name = "genotype",
     values = c("TET2" = 'darkgreen', "DNMT3A" = 'darkred')
   ) +
   theme(text = element_text(size = text_size), axis.title = element_text(size = text_size), legend.text = element_text(size = text_size), legend.title = element_text(size = text_size),
        plot.title = element_text(size = text_size,  hjust = 0.5),
        legend.position = c(0.01, 0.99), legend.justification = c("left", "top"),
        legend.background = element_rect(fill = "white", color = "black")) # Optional: adds a box)
pl

pdf('figures/heritability/scatterplot_WT_vs_CH.pdf', width = 8, height = 8)
pl
dev.off()
