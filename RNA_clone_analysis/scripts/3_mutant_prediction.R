library(ggplot2)
library(ggrepel)
library(gridExtra)
library(parallel)
library(data.table)
library(DESeq2)
library(MASS)
library(pROC)
source('R/mutant_prediction_routines.R')
load('../RNA_analysis/data/processed/scrna.RData')
load('data/processed/clone_info.RData')
hash.annot <- readRDS(file = "../RNA_analysis/data/processed/hashtag_info.rds")
hash.match <- readRDS(file = "../RNA_analysis/data/processed/hashtag_meta.rds")
expids <- c('exp123', 'exp127')
smpls <- c('RNA1', 'RNA2', 'RNA3', 'RNA4')


redsp <- do.call(rbind, lapply(1:10, function(cl){
  print(cl)
  cls <- names(clust.final)[clust.final == cl]
  cls <- intersect(cls, rownames(cell.clone))
  #
  #cls <- cls[hash.match[sampleid[cls],]$exp == '127']
  #
  geno_score <- as.numeric(as.factor(hash.match[cell.clone[cls,]$id,]$genotype))
  exp_score <- as.numeric(as.factor(hash.match[cell.clone[cls,]$id,]$exp))
  smpl <- paste(cell.clone[cls,]$expid, cell.clone[cls,]$sample, sep = '_')
  #
  ft <- lm(p2pca[cls,] ~ geno_score + smpl)
  reduc <- ft$residuals + as.matrix(geno_score) %*% t(as.matrix(ft$coefficients[2,]))
  reduc <- p2pca[cls,]
  return(reduc)
}))


# learn cluster-specific scores
p.train <- 0.5
for (gen in c('DNMT3A', 'TET2')){
  if (gen == 'TET2'){
    expid <- 127
  }else{
    expid <- 123
  }
# fix cell phenotype
pheno <- as.numeric(hash.match[sampleid,]$genotype == gen)
names(pheno) <- names(sampleid)
pheno[hash.match[sampleid,]$exp != expid] <- NA

mut.sig <- lapply(1:10, function(cl){
  cls <- names(clust.final)[clust.final == cl]
  cls <- intersect(cls, rownames(cell.clone))
  cls1 <- cls[!is.na(pheno[cls])]
  #
  cls.train <- sample(cls1, length(cls1) * p.train)
  #cls.train <- cls1[! hash.match[sampleid[cls1],]$RNA %in% c('1', '3')  ]
  cls.test <- setdiff(cls1, cls.train)
  ft <- ppr(cbind(1, redsp[cls.train,]), pheno[cls.train], nterms = 1, max.terms = 10)

  pred <- predict(ft, cbind(1, redsp[cls1,]))

  data.frame(score = pred, pheno = pheno[cls1], use = ifelse(cls1 %in% cls.train, 'train', 'test'))
})
mut.sig.qq <- lapply(mut.sig, function(tbl){
  #tbl$score[tbl$pheno == 1]
  x <- unlist(lapply(tbl$score, function(sc)
    mean(tbl$score[tbl$use == 'train' & tbl$pheno == 1] >= sc)
  ))
  y <- unlist(lapply(tbl$score, function(sc)
    mean(tbl$score[tbl$use == 'train' & tbl$pheno == 0] >= sc)
  ))
  tbl$qmut <- x
  tbl$qwt <- y
  return(tbl)
})

hash.match1 <- hash.match
hash.match1$RNA <- factor(paste(hash.match$exp, hash.match$RNA, sep = '_'))
#
mut.sig.de <- lapply(1:10, function(cl){
  print(paste('cluster', cl))
  cls <- names(clust.final)[clust.final == cl]
  cls <- intersect(cls, rownames(cell.clone))
  cls1 <- cls[!is.na(pheno[cls])]
  #
  cls.train <- sample(cls1, length(cls1) * p.train)
  #cls.train <- cls1[! hash.match[sampleid[cls1],]$RNA %in% c('1', '3')  ]
  cls.test <- setdiff(cls1, cls.train)

  mat <- do.call(cbind, lapply(unique(sampleid[cls.train]), function(smpl){
    #print(smpl)
    ind <- sampleid[cls.train] %in% smpl
    rowSums(cd[, cls.train][, ind])
  }))
  colnames(mat) <- unique(sampleid[cls.train])

  print('perform DE..')
  # estimate DE genes
  dds_mut <- DESeqDataSetFromMatrix(countData = mat,
                                    colData = hash.match1[colnames(mat), , drop = F],
                                    design = ~ genotype)

  dds_mut <- DESeq(dds_mut)
  res_mut <- lfcShrink(dds_mut, coef=paste("genotype_",gen, "_vs_WT", sep = ''), type="normal")

  # measure cell-specific mutation score
  #gg <- rownames(res_mut)[!is.na(res_mut$padj) & res_mut$pvalue < 1e-1]
  gg <- rownames(res_mut)[order(res_mut$pvalue)][1:100]
  str(gg)
  gg <- gg[ rowSums(cd[gg, cls1] > 0) > length(cls1) * 0.1 ]
  lfc <- res_mut[gg,]$log2FoldChange
  size1 <- rowSums(cd[gg, cls1])

  print('estimate cell-specific scores')
  #cls.est <- sample(cls1, 100)
  cls.est <- cls1
  cell.score <- do.call(rbind, mclapply(cls.est, function(cll){
    ft1 <- glm.nb(cd[gg, cll] ~ lfc + offset(log(size1[gg]+0)),
                  start = c(1, 0))
    summary(ft1)$coefficients[2,c(1,4)]
  }, mc.cores = 30))
  rownames(cell.score) <- cls.est
  #genov <- hash.match[sampleid[cls.est],]$genotype
  #qplot(cell.score[,1], col = genov, geom = 'density')

  run.roc <- FALSE
  if (run.roc){
    ind <- cls1 %in% cls.train
    roc_train <- roc(pheno[ind], cell.score[ind,1])
    roc_test <- roc(pheno[!ind], cell.score[!ind,1])
    plot(roc_test, print.auc = TRUE, main = cl)
    plot.roc(roc_train, add = TRUE, col = "red")
  }

  data.frame(score = cell.score[,1], pheno = pheno[cls1], use = ifelse(cls1 %in% cls.train, 'train', 'test'))
})
#
mut.sig.de.qq <- lapply(mut.sig.de, function(tbl){
  #tbl$score[tbl$pheno == 1]
  x <- unlist(lapply(tbl$score, function(sc)
    mean(tbl$score[tbl$use == 'train' & tbl$pheno == 1] >= sc)
  ))
  y <- unlist(lapply(tbl$score, function(sc)
    mean(tbl$score[tbl$use == 'train' & tbl$pheno == 0] >= sc)
  ))
  tbl$qmut <- x
  tbl$qwt <- y
  return(tbl)
})

xsave <- list(mut.sig = mut.sig, mut.sig.qq = mut.sig.qq, mut.sig.de = mut.sig.de, mut.sig.de.qq = mut.sig.de.qq)
saveRDS(xsave, paste('data/processed/', gen, '_scores.rds', sep =''))
}


### show two ROC (DNMT3A and TET2) curves together
gen <- 'DNMT3A'
mut1 <- readRDS(paste('data/processed/', gen, '_scores.rds', sep =''))
mut1.cell.de <- do.call(rbind, mut1$mut.sig.de.qq)
mut1.clone.stat.de <- summary_mut_clone(mut1.cell.de, clone.info)
#
gen <- 'TET2'
mut2 <- readRDS(paste('data/processed/', gen, '_scores.rds', sep =''))
mut2.cell.de <- do.call(rbind, mut2$mut.sig.de.qq)
mut2.clone.stat.de <- summary_mut_clone(mut2.cell.de, clone.info)

pdf(paste('figures/','mutant_roc_summary.pdf', sep = ''), width = 6, height = 3)
par(mfrow=c(1,2))
for (i in 1:2){
  if (i == 1){
    tbl <- mut1.cell.de
    ind <- tbl$use == 'train'
    roc_train <- roc(tbl$pheno[ind], tbl$qwt[ind])
    roc_test <- roc(tbl$pheno[!ind], tbl$qwt[!ind])
  }else{
    tbl <- mut1.clone.stat.de
    ind <- clone.info$size >= 10
    roc_train <- roc(!clone.info$genotype[ind] == 'WT', tbl[ind,1])
    roc_test <- roc(!clone.info$genotype[ind] == 'WT', tbl[ind,2])
  }
  plot(roc_test, print.auc = FALSE, cex.lab = 2, cex.main = 2, font.main = 1, lwd = 5, col = col.genotype['DNMT3A'], lty = 1,
       main = ifelse(i == 1, 'cell level', 'clone level'))
  plot.roc(roc_train, add = TRUE, lwd = 4, col = col.genotype['DNMT3A'], lty = 2)
  #
  if (i == 1){
    tbl <- mut2.cell.de
    ind <- tbl$use == 'train'
    roc_train <- roc(tbl$pheno[ind], tbl$qwt[ind])
    roc_test <- roc(tbl$pheno[!ind], tbl$qwt[!ind])
  }else{
    tbl <- mut2.clone.stat.de
    ind <- clone.info$size >= 10
    roc_train <- roc(!clone.info$genotype[ind] == 'WT', tbl[ind,1])
    roc_test <- roc(!clone.info$genotype[ind] == 'WT', tbl[ind,2])
  }
  plot.roc(roc_test, add = TRUE, lwd = 4, col = col.genotype['TET2'], lty = 1)
  plot.roc(roc_train, add = TRUE, lwd = 4, col = col.genotype['TET2'], lty = 2)
  legend('bottomright', c('DNMT3A', 'TET2'), col = col.genotype[c('DNMT3A', 'TET2')],
         lwd = 4, bty = 'n', cex = 1)
}
dev.off()


#
ind <- clone.info$size > 50
#
tbl <- mut1.clone.stat.de; gen <- 'DNMT3A'
pl1 <- qplot(tbl[ind,1], tbl[ind,2], main = paste(gen, 'clone score'), xlab = 'train', ylab = 'test',
      #col = clone.info$sample[ind],
      col = clone.info$genotype[ind],
      #xlim = c(0,1), ylim = c(0,1),
      size = clone.info$size[ind], geom = c('point', 'smooth'), method = 'lm') + theme_bw() +
  scale_color_manual(values = col.genotype[c('WT', gen)]) +
  scale_size(guide = 'none') +
  theme(text = element_text(size = 15), legend.title = element_blank(), legend.position = c(0.75, 0.2))
pl1

pdf(paste('figures/mutant_', gen, '_clone_scores.pdf', sep = ''), width = 3.5, height = 3.5)
pl1
dev.off()

#
tbl <- mut2.clone.stat.de; gen <- 'TET2'
pl2 <- qplot(tbl[ind,1], tbl[ind,2], main = paste(gen, 'clone score'), xlab = 'train', ylab = 'test',
             #col = clone.info$sample[ind],
             col = clone.info$genotype[ind],
             #xlim = c(0,1), ylim = c(0,1),
             size = clone.info$size[ind], geom = c('point', 'smooth'), method = 'lm') + theme_bw() +
  scale_color_manual(values = col.genotype[c('WT', gen)]) +
  scale_size(guide = 'none') +
  theme(text = element_text(size = 15), legend.title = element_blank(), legend.position = c(0.7, 0.2))
pl2

pdf(paste('figures/mutant_', gen, '_clone_scores.pdf', sep = ''), width = 3.5, height = 3.5)
pl2
dev.off()


