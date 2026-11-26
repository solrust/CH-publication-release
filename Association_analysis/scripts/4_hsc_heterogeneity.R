load('../RNA_analysis/data/processed/scrna.RData')
load('../RNA_analysis/data/processed/mutant_de.RData')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
load('../TF_overexpression/data/processed/TF_trends.RData')
source('R/association_routines.R')
library(Matrix)
library(RColorBrewer)
library(parallel)
library(Matrix)
library(MASS)
library(glmmTMB)
library(ggplot2)
library(data.table)
library(gridExtra)
hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
tfs.list <- readRDS('../RNA_analysis/data/genesets/processed/tfs.rds')
tfs.all <- unlist(tfs.list)
save_folder <- 'default'
data_folder <- paste('data/processed/', save_folder, '/', sep = '')
load(paste(data_folder, 'association_summary.RData', sep = ''))



cls <- names(clust.final)[clust.final == 1]
cutf <- 4
n.cmp <- 10
### generate PCs for HSC cluster
gns <- rownames(gsea.clone$gscore.clone)[apply(abs(gsea.clone$gscore.clone), 1, max) > cutf]
prr <- prcomp(p2counts[cls, gns], scale = T)
prr$x[,2] <- -prr$x[,2]; prr$rotation[,2] <- -prr$rotation[,2]

### project phenotypic signatures
pheno.projs <- do.call(cbind, lapply(1:7, function(i){
  print(i)
  #gns <- rownames(gsea.clone$gscore.clone)[apply(abs(gsea.clone$gscore.clone), 1, max) > cutf]
  gns <- rownames(gsea.clone$gscore.clone)[abs(gsea.clone$gscore.clone[,i])> cutf]
  #prr <- prcomp(p2counts[cls, gns], scale = T)
  y <- (gsea.clone$gscore.clone[gns, i])
  ft <- lm(y ~ prr$rotation[gns, 1:n.cmp])
  y <- prr$x[,1:n.cmp] %*% ft$coefficients[-1]
  return(y)
  #cor(gsea.clone$gscore.clone[gns,i], ft$fitted.values)
}))
colnames(pheno.projs) <- colnames(gsea.clone$gscore.clone)

### explore projections of signatures:
png('figures/hsc/hsc_phenotype_emb_projections.png', width = 3 * 150, height = 3 * 150)
par(mfrow=c(3,3), mar = rep(2, 4))
for (i in 1:7){
  cols <- colorRampPalette(c("blue","grey90", "red" ))(n = 10)[as.numeric(cut(rank(pheno.projs[,i]), 10))]
  plot(emb[cls,], pch = 19, cex = 0.5, col = adjustcolor(cols, 0.5), xaxt = 'n', yaxt = 'n', xlim = c(-2, 5), ylim = c(2, 10),
       main = colnames(pheno.projs)[i], cex.main = 2, font.main = 1)
}
dev.off()

png('figures/hsc/hsc_phenotype_new.png', width = 3 * 150, height = 2 * 150)
#pheno.projs <- pheno.projs1
y <- cor(prr$x[,1:3], pheno.projs)
par(mfrow=c(2,3), mar = c(2, 2, 0, 2))
col.types <- c('RdBu', 'PiYG', 'PRGn')
nms.types <- c('growth', 'output', 'unknown')
reord.pnames <- c("lymphoid bias", "clone size", "growth rate", "erythroid bias", "myeloid bias", "HSC bias", "dendritic bias")
reord.pnames <- rev(reord.pnames)
for (i in 1:3){
  #cols <-  colorRampPalette(c("blue","grey90", "red" ))(n = 10)[as.numeric(cut(rank(prr$x[,i]), 10))]
  cols <-  brewer.pal(10, col.types[i])[as.numeric(cut(rank(-prr$x[,i]), 10))]
  plot(emb[cls,], pch = 19, cex = 0.5, col = adjustcolor(cols, 0.5), bty = 'n', xaxt = 'n', yaxt = 'n',
       #main = nms.types[i],
       font.main=1, cex.main = 2, xlim = c(-2, 5), ylim = c(3, 9))
}
for (i in 1:3){
  #cols <-  brewer.pal(10, col.types[i])[as.numeric(cut(-y[i,], 10))]
  cols <-  brewer.pal(10, col.types[i])[as.numeric(cut(-y[i,], breaks = seq(-1, 1, length.out = 10)))]
  barplot(y[i,reord.pnames], horiz = T, col =  cols[match(reord.pnames, colnames(y))], cex.axis = 2, xlim = c(-1,1), names.arg = '')
}
dev.off()


### test for mutation-specific effects along hsc PCs:
cd.hsc <- cd[,cls]
meta.hsc <- data.frame(hash.match[sampleid[colnames(cd.hsc)], c('genotype', 'exp')],
                       sample = sampleid[colnames(cd.hsc)])
meta.hsc$genotype <- factor(meta.hsc$genotype, levels = c('WT', 'DNMT3A', 'TET2'))
rownames(meta.hsc) <- cls
fpm.hsc <- t(t(cd.hsc) / colSums(cd.hsc))
gns <- rownames(test.basic.list[[1]])
n.cores <- 25

test.pc.mutant <- lapply(1:3, function(ipc){
  print(paste('PC:', ipc))
  print('run test.pc')
  test.pc <- fit.nb(t(cd.hsc[gns,]), meta.hsc, prr$x[cls, ipc],
                    do.batch = TRUE, do.interactions = FALSE,
                    nc.cores = n.cores)
  rownames(test.pc) <- gns

  print('run DNMT3A-specific PC associations')
  ind <- meta.hsc$exp == '123'
  test.pc.dnmt3a <- fit.nb(t(cd.hsc[gns,ind]), meta.hsc[ind,], prr$x[cls, ipc][ind],
                           do.batch = FALSE, do.interactions = TRUE,
                           nc.cores = n.cores)
  rownames(test.pc.dnmt3a) <- gns

  print('run TET2-specific PC associations')
  ind <- meta.hsc$exp == '127'
  test.pc.tet2 <- fit.nb(t(cd.hsc[gns, ind]), meta.hsc[ind,], prr$x[cls, ipc][ind],
                         do.batch = FALSE, do.interactions = TRUE,
                         nc.cores = n.cores)
  rownames(test.pc.tet2) <- gns

  return(list(test.pc = test.pc, test.pc.dnmt3a = test.pc.dnmt3a, test.pc.tet2 = test.pc.tet2))
})
names(test.pc.mutant) <- paste('pc', 1:3, sep = '')

### generate tables
table.pc.mutant <- lapply(1:length(test.pc.mutant), function(ipc){
  test.list <- test.pc.mutant[[ipc]]
  #test.list <- test.pc.mutant[[1]]
  test.tbl <- cbind(test.list$test.pc[,c(1, 5)],
                    test.list$test.pc.dnmt3a[,c(3, 6)],
                    test.list$test.pc.tet2[,c(3, 6)])
  nm <- paste('pc', ipc, sep = '')
  nm <- 'pc'
  colnames(test.tbl) <- c(nm, paste(nm, '.p', sep = ''), paste(nm, 'dnmt3a', sep = ''),  paste(nm, 'dnmt3a.p', sep = ''), paste(nm, 'tet2', sep = ''),  paste(nm, 'tet2.p', sep = ''))
  dpc <- max(prr$x[cls, ipc]) - min(prr$x[cls, ipc])
  test.tbl[,c(1,3,5)] <- test.tbl[,c(1,3,5)] * dpc / log(2)
  return(data.frame(test.tbl))
})
names(table.pc.mutant) <- names(test.pc.mutant)
save(pheno.projs, prr, cd.hsc, meta.hsc, test.pc.mutant, table.pc.mutant, file = 'data/processed/hsc/hsc_axes.RData')
