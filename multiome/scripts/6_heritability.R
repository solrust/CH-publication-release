library(parallel)
load('../DNA_analysis/data/processed/DNA_clones.RData')

### choose HSC cells
cls.hsc <- names(obj.combined$projection)[obj.combined$projection == 1 & !is.na(obj.combined$projection)]
str(cls.hsc)

### specify PC space
FeaturePlot(obj.combined, features = 'PC_5', min.cutoff = 0,
            cols = c('grey90', "darkred"), reduction = 'wnn.umap', ncol = 2, pt.size = I(0.5), combine = FALSE)
#
tPCA <- obj.combined@reductions$pca@cell.embeddings[cls.hsc, ]
#tPCA <- obj.combined@reductions$lsi@cell.embeddings[cls.hsc, ]
rownames(tPCA) <- unlist(lapply(strsplit(rownames(tPCA), '-'), function(x) x[1]))

### choose clones
clone.size <- 50
clones <- rownames(clone.info)[clone.info[, 'size'] >= clone.size]
str(clones)
prclone <- prcomp(clone.states2[clones,], scale. = T)
clone.pheno <- cbind(clone.states2[clones,], prclone$x[clones,], rate = rate[clones])


pheno_score <- clone.pheno[,'PC1']; names(pheno_score) <- rownames(clone.pheno)
#  number of PCs to estimate heritability
x <- mclapply(2:50, function(npcs){
  y <- lapply(1:50, function(j) pheno_heritability(tPCA, pheno_score, cell.clone.multi, n.tcomp = npcs,
                                                   train.percent = 0.6)
  )
  Reduce('+', y) / 50
}, mc.cores = 50)
r2.train <- do.call(rbind, lapply(x, function(y) y[,1]))
r2.test <- do.call(rbind, lapply(x, function(y) y[,2]))
ggplot() + geom_point(aes(2:50, rowMeans(r2.train), col = 'train')) +
  geom_point(aes(2:50, rowMeans(r2.test), col = 'test'))


### estimate phenotype-specific R2
npcs <- 20
n <- 100

pheno_state_list <- lapply(colnames(clone.pheno), function(nm){
  print(nm)
  #nm <- 'PC1'
  #nm <- '4'
  #nm <- 'rate'
  pheno_score <- clone.pheno[,nm]
  names(pheno_score) <- rownames(clone.pheno)
  pheno_stat <- mclapply(1:n, function(j)
    pheno_heritability(tPCA, pheno_score, cell.clone.multi, train.percent = 0.6,
                       n.tcomp = npcs)
    , mc.cores = 20)
  pheno_stat
})
names(pheno_state_list) <- colnames(clone.pheno)
#saveRDS(pheno_state_list, file = 'data/processed/heritability/heritability_rna.rds')
#pheno_state_list_rna <- pheno_state_list

pdf('figures/heritability/boxplots_atac_test.pdf', width = 3.5, height = 4)
for (i in 1:length(pheno_state_list)){
  pheno_stat <- pheno_state_list[[i]]
  pl <- plot_boxplot(pheno_stat, title = names(pheno_state_list)[i], do.test = F)
  print(pl)
}
dev.off()

pv_pheno <- do.call(rbind, lapply(pheno_state_list, function(pheno_stat){
  res <- do.call(rbind, lapply(pheno_stat, function(x)
    data.frame(geno = rownames(x), scores = x[, ifelse(do.test, 2, 1)]) ))
  res$geno[res$geno %in% c('WT1', 'WT2')] <- 'WT'
  res$geno <- factor(res$geno, levels = c('DNMT3A', 'WT', 'TET2'))
  est_pvalue(res)
}))
#saveRDS(pv_pheno, file = 'data/processed/heritability/heritability_pvalues.rds')



#### routines
plot_boxplot <- function(pheno_stat, title = NULL, do.test = TRUE){
  res <- do.call(rbind, lapply(pheno_stat, function(x)
    data.frame(geno = rownames(x), scores = x[, ifelse(do.test, 2, 1)]) ))
  res$geno[res$geno %in% c('WT1', 'WT2')] <- 'WT'
  res$geno <- factor(res$geno, levels = c('DNMT3A', 'WT', 'TET2'))
  pl <- ggplot(res, aes(x = geno, y = scores, col = geno)) + theme_bw() +
    labs(y = expression(paste('HSC state to phenotype, ', R^2)), x = ' genotype',
         title = title) +
    geom_boxplot(width = 0.5, size = 1.5, outlier.shape=NA, alpha = 1) +
    #scale_color_manual(values = col.genotype) +
    #scale_color_manual(values = c('DNMT3A' = 'red', 'TET2' = 'blue', 'WT' = 'grey')) +
    geom_point(position=position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5),
               size = 0.5, alpha = 0.25) + #ylim(0.0, max(res$scores)) +
    theme(text = element_text(size = 20), legend.position = "none")
  pl
  return(pl)
}

pheno_heritability <- function(tPCA = NULL, pheno_score = NULL, cell.clone = NULL, n.tcomp = 10, train.percent = NULL){
  r2 <- do.call(rbind, lapply(
    #unique(cell.clone$id),
    c('WT1', 'WT2', 'TET2', 'DNMT3A'),
    function(sid){
      print(sid)
      if (sid == 'WT1'){
        sid <- 'WT'
        expid <- c('123', '127')
      }else if (sid == 'WT2'){
        sid <- 'WT'
        expid <- c('127', '123')
      }
      if (sid == 'TET2') expid <- '127'
      if (sid == 'DNMT3A') expid <- '123'
      # retain HSCs with clone information
      cls <- rownames(tPCA)
      #
      cloneid <- cell.clone[cls,]$cloneid
      #cloneid <- cell.clone[cls,]$larry
      names(cloneid) <- cls
      #
      cls <- cls[!is.na(cloneid) &
                   cloneid %in% names(pheno_score)[!is.na(pheno_score)] &
                   #cell.clone[cls,]$id %in% sid
                   hash.match[cell.clone[cls,]$id,]$genotype %in% sid &
                   hash.match[cell.clone[cls,]$id,]$exp %in% expid #& cell.clone[cls,]$sample != 'C2'
      ]
      #
      pcomp <- tPCA[cls, 1:n.tcomp]
      pheno <- pheno_score[cloneid[cls]]

      if (is.null(train.percent)){
        smpl.train <- smpl.test <- 1:length(pheno)
      }else{
        smpl.train <- sample(1:length(pheno), length(pheno) * train.percent)
        smpl.test <- setdiff(1:length(pheno), smpl.train)
      }

      # fit clone phenotype to tPCs
      ft <- lm(pheno[smpl.train] ~ pcomp[smpl.train, ])
      coefs <- summary(ft)$coefficients[,1]
      x <- coefs[1] + pcomp[smpl.test, ] %*% coefs[-1]
      r2.train <- summary(ft)$r.squared
      r2.test <- 1 - sum((pheno[smpl.test] - x)^2) / sum((pheno[smpl.test] - mean(pheno[smpl.test]))^2)
      #r2.test <- cor(pheno[smpl.test], x)
      return(c(r2.train, r2.test))
    }))
  rownames(r2) <- c('WT1', 'WT2', 'TET2', 'DNMT3A')#c('WT', 'TET2', 'DNMT3A')#unique(cell.clone$id)
  #rownames(r2) <- unique(cell.clone$id)
  colnames(r2) <- c('train', 'test')
  return(r2)
}

est_pvalue <- function(res){
  x1 <- res$scores[res$geno == 'DNMT3A']
  x2 <- res$scores[res$geno == 'WT']
  x3 <- res$scores[res$geno == 'TET2']
  p1 <- mean(unlist(lapply(x1, function(x) mean(x > x2))))
  p2 <- mean(unlist(lapply(x3, function(x) mean(x < x2))))
  c(DNMT3A = min(p1, 1-p1), TET2 = min(p2, 1-p2))
}
