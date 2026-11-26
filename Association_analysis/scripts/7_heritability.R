load('../DNA_analysis/data/processed/DNA_clones.RData')

### choose HSC cells
cls.hsc <- names(clust.final)[clust.final == 1]
str(cls.hsc)

### specify PC space
tPCA <- p2pca[cls.hsc, ]

### choose clones
clone.size <- 50
clones <- rownames(clone.info)[clone.info[, 'size'] >= clone.size]
str(clones)
prclone <- prcomp(clone.states2[clones,], scale. = T)
clone.pheno <- cbind(clone.states2[clones,], prclone$x[clones,], rate = rate[clones])

pheno_score <- clone.pheno[,'PC1']; names(pheno_score) <- rownames(clone.pheno)
#  number of PCs to estimate heritability
x <- mclapply(2:30, function(npcs){
  y <- lapply(1:50, function(j) pheno_heritability(tPCA, pheno_score, cell.clone, n.tcomp = npcs,
                     train.percent = 0.8)
  )
  Reduce('+', y) / 50
}, mc.cores = 25)
r2.train <- do.call(rbind, lapply(x, function(y) y[,1]))
r2.test <- do.call(rbind, lapply(x, function(y) y[,2]))
ggplot() + geom_point(aes(2:30, rowMeans(r2.train), col = 'train')) +
           geom_point(aes(2:30, rowMeans(r2.test), col = 'test'))


### estimate phenotype-specific R2
npcs <- 15
n <- 100
pheno_state_list <- mclapply(colnames(clone.pheno)[c(1:10, 21)], function(nm){
  print(nm)
  pheno_score <- rank(clone.pheno[,nm])
  names(pheno_score) <- rownames(clone.pheno)
  pheno_stat <- mclapply(1:n, function(j)
    pheno_heritability(tPCA, pheno_score, cell.clone, do.permute = F, do.bootstrap = T, #train.percent = 0.6,
                       n.tcomp = npcs)
  , mc.cores = 1)
  pheno_stat
}, mc.cores = 11)
names(pheno_state_list) <- colnames(clone.pheno)[c(1:10, 21)]
#saveRDS(pheno_state_list, file = 'data/processed/heritability/heritability_estimates.rds')

#
plist <- pheno_state_list
pdf('figures/heritability/boxplots_final.pdf', width = 3.5, height = 4)
for (i in 1:length(plist)){
  pheno_stat <- plist[[i]]
  pl <- plot_boxplot(pheno_stat, title = names(plist)[i], do.test = F)
  print(pl)
}
dev.off()

do.test <- F
pv_pheno <- do.call(rbind, lapply(pheno_state_list, function(pheno_stat){
  res <- do.call(rbind, lapply(pheno_stat, function(x)
    data.frame(geno = rownames(x), scores = x[, ifelse(do.test, 2, 1)]) ))
  res$geno[res$geno %in% c('WT1', 'WT2')] <- 'WT'
  res$geno <- factor(res$geno, levels = c('DNMT3A', 'WT', 'TET2'))
  est_pvalue(res)
}))
#saveRDS(pv_pheno, file = 'data/processed/heritability/heritability_pvalues.rds')



############## routines

plot_boxplot <- function(pheno_stat, title = NULL, do.test = TRUE){
  res <- do.call(rbind, lapply(pheno_stat, function(x)
    data.frame(geno = rownames(x), scores = x[, ifelse(do.test, 2, 1)]) ))
  res$geno[res$geno %in% c('WT1', 'WT2')] <- 'WT'
  res$geno <- factor(res$geno, levels = c('DNMT3A', 'WT', 'TET2'))
  pl <- ggplot(res, aes(x = geno, y = scores, col = geno)) + theme_bw() +
    labs(y = expression(paste('HSC state to phenotype, ', R^2)), x = ' genotype',
         title = title) +
    geom_boxplot(width = 0.5, size = 1.5, outlier.shape=NA, alpha = 1) +
    scale_color_manual(values = col.genotype) +
    #scale_color_manual(values = c('DNMT3A' = 'red', 'TET2' = 'blue', 'WT' = 'grey')) +
    geom_point(position=position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5),
               size = 0.5, alpha = 0.1) + ylim(0.0, max(res$scores)) +
    theme(text = element_text(size = 20), legend.position = "none")
  pl
  return(pl)
}


pheno_heritability <- function(tPCA = NULL, pheno_score = NULL, cell.clone = NULL, n.tcomp = 10, train.percent = NULL,
                               do.permute = FALSE, do.bootstrap = FALSE){
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
    t1 <- Sys.time()
    cls <- rownames(tPCA)
    cell.clone1 <- cell.clone[cls,]
    hash.match1 <- hash.match; hash.match1$sample <- rownames(hash.match)
    hash.match1 <- hash.match1[cell.clone1$id,]; rownames(hash.match1) <- cls
    cls <- cls[!is.na(cell.clone1$cloneid) &
                 cell.clone1$cloneid %in% names(pheno_score)[!is.na(pheno_score)] &
                 #cell.clone[cls,]$id %in% sid
                 hash.match1$genotype %in% sid &
                 hash.match1$exp %in% expid
               ]
    #
    t2 <- Sys.time()
    print(t2-t1)
    pcomp <- tPCA[cls, 1:n.tcomp]
    pheno <- pheno_score[cell.clone[cls,]$cloneid]; names(pheno) <- cls

    if (do.permute){
      do.sample.control <- FALSE
      if (do.sample.control){
      cls.permute <- unlist(lapply(unique(hash.match1[cls,]$sample), function(smpl){
        ind <- hash.match1[cls,]$sample == smpl
        cls1 <- cls[ind]
        cls2 <- sample(cls1)
        names(cls2) <- cls1
        return(cls2)
      }))
      }else{
        cls.permute <- sample(cls)
        names(cls.permute) <- cls
      }
      cls.permute <- cls.permute[cls]
    }else{
      cls.permute <- T
    }

    if (is.null(train.percent)){
      smpl.train <- smpl.test <- 1:length(pheno)
      if (do.bootstrap){
        smpl.train <- smpl.test <- sample(1:length(pheno), replace = T)
      }
    }else{
      smpl.train <- sample(1:length(pheno), length(pheno) * train.percent)
      smpl.test <- setdiff(1:length(pheno), smpl.train)
    }

    # fit clone phenotype to tPCs
    run.ridge <- TRUE
    if (run.ridge){
      cvfit <- cv.glmnet(pcomp[smpl.train, ], pheno[smpl.train])
      best_ridge <- glmnet(pcomp[smpl.train, ], pheno[smpl.train], alpha = 0, lambda = cvfit$lambda.min)

      pred.train <- predict(best_ridge, s = cvfit$lambda.min, newx = pcomp[smpl.train, ])
      r2.train <- 1 - sum((pheno[smpl.train] - pred.train)^2) / sum((pheno[smpl.train] - mean(pheno[smpl.train]))^2)
      #
      pred.test <- predict(best_ridge, s = cvfit$lambda.min, newx = pcomp[smpl.test, ])
      r2.test <- 1 - sum((pheno[smpl.test] - pred.test)^2) / sum((pheno[smpl.test] - mean(pheno[smpl.test]))^2)
      c(r2.train, r2.test)
    }else{
      if (do.permute){
        ft <- lm(pheno[smpl.train][ cls.permute[names(pheno[smpl.train])] ] ~ pcomp[smpl.train, ])
      }else{
        ft <- lm(pheno[smpl.train] ~ pcomp[smpl.train, ])
      }
      coefs <- summary(ft)$coefficients[,1]
      x <- coefs[1] + pcomp[smpl.test, ] %*% coefs[-1]
      r2.train <- summary(ft)$r.squared
      r2.test <- 1 - sum((pheno[smpl.test] - x)^2) / sum((pheno[smpl.test] - mean(pheno[smpl.test]))^2)
    }
    #r2.test <- cor(pheno[smpl.test], x)
    return(c(r2.train, r2.test))
  }))
  rownames(r2) <- c('WT1', 'WT2', 'TET2', 'DNMT3A')#c('WT', 'TET2', 'DNMT3A')#unique(cell.clone$id)
  #rownames(r2) <- unique(cell.clone$id)
  colnames(r2) <- c('train', 'test')
  return(r2)
}

#ft <- lm(pheno[smpl.train] ~ pcomp[smpl.train, ])

library(glmnet)
lambda_seq <- 10^seq(4, -4, by = -.2)
# Using glmnet function to build the ridge regression in r
fit <- glmnet(pcomp[smpl.train, ], pheno[smpl.train], alpha = 0, lambda  = lambda_seq)
summary(fit)
plot(fit)

xpheno <- pheno[smpl.train]
xpheno <- pheno[smpl.train][ cls.permute[names(pheno[smpl.train])] ]
cvfit <- cv.glmnet(pcomp[smpl.train, ], xpheno)
cvfit <- cv.glmnet(pcomp[smpl.train, ], xpheno,
                   lambda = lambda_seq)
plot(cvfit)
cvfit$lambda.min
best_ridge <- glmnet(pcomp[smpl.train, ], xpheno, alpha = 0, lambda = cvfit$lambda.min)
coef(best_ridge)
pred <- predict(best_ridge, s = cvfit$lambda.min, newx = pcomp[smpl.train, ])
coefs <- summary(ft)$coefficients[,1]
x <- coefs[1] + pcomp[smpl.test, ] %*% coefs[-1]
r2.test <- 1 - sum((xpheno - pred)^2) / sum((xpheno - mean(xpheno))^2)
r2.test

est_pvalue <- function(res){
  x1 <- res$scores[res$geno == 'DNMT3A']
  x2 <- res$scores[res$geno == 'WT']
  x3 <- res$scores[res$geno == 'TET2']
  p1 <- mean(unlist(lapply(x1, function(x) mean(x > x2))))
  p2 <- mean(unlist(lapply(x3, function(x) mean(x < x2))))
  c(DNMT3A = min(p1, 1-p1), TET2 = min(p2, 1-p2))
}
