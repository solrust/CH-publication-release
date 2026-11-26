source('R/larry_calling_routines.R')
library(parallel)
library(ggplot2)
hash.match <- readRDS(file = "../RNA_analysis/data/processed/hashtag_meta.rds")
load('data/processed/clone_info.RData')

multiome.smpls <- c('C1', 'C2', 'T', 'D')
multiome.match <- list(C1 = c('123_A0', '123_A1', '127_B4'),
                       D = c('123_B0', '123_B1', '123_B2', '123_B3'),
                       T = c('127_A3', '127_A4', '127_B0', '127_B1'),
                       C2 = c('127_B2', '127_B3', '123_A2', '127_A4')
                       )

### call and save multiome-larry barcodes
for (smpl in multiome.smpls){
  print(smpl)
  fpath <- paste('../../larry_atac/trex_', smpl, '/reads_clean.txt', sep = '')
  trios <- upload_trios(fpath, pref = paste(smpl, '_', sep = ''))
  larry <- identify_barcodes(trios = trios, do.plot = FALSE)
  saveRDS(larry, paste("data/processed/multiome_larry_",smpl,".rds", sep = ''))
}


### save diagnostic plots
for (smpl in multiome.smpls){
  print(smpl)
  larry <-readRDS(paste("data/processed/multiome_larry_",smpl,".rds", sep = ''))
  #png(paste("figures/multiome_larry_", smpl, ".png", sep = ''), width = 200 * 2, height = 200 * 3)
  par(mfrow=c(3,2), mgp = c(2.5,1,0), mar = c(4, 4, 1, 1))
  visualize_clone_processing(larry, cexlab = 1.5, do.par = F)
  #dev.off()
}

### identify cell-clone relationships
cell.clone.multi <- do.call(rbind, lapply(multiome.smpls, function(smpl){
  larry <-readRDS(paste("data/processed/multiome_larry_", smpl, ".rds", sep = ''))
  cloneid <- larry$cloneid
  cls <- names(cloneid)
  cell.clone <- as.data.frame(cbind(smpl, cloneid[cls]))
  colnames(cell.clone) <- c('sample', 'clone')
  return(cell.clone)
}))
cell.clone.multi <- cell.clone.multi[!is.na(cell.clone.multi$clone),]
cell.clone.multi$cloneid <- paste(cell.clone.multi$sample, '_', cell.clone.multi$clone, sep = '')
saveRDS(cell.clone.multi, file = "data/processed/cell_clone_multi.rds")
dim(cell.clone.multi)
# 40,000 down to 17000 with LARRY barcodes. This is across all samples

### store clone infor
clone.info.multi <- do.call(rbind, lapply(unique(cell.clone.multi$cloneid), function(cl){
  cls <- rownames(cell.clone.multi)[!is.na(cell.clone.multi$clone) & cell.clone.multi$cloneid == cl]
  str(cls)
  xinf <- data.frame(table(cell.clone.multi[cls, c('sample', 'clone')]))
  xinf
  freq <- max(xinf$Freq) / sum(xinf$Freq)
  cbind(xinf[which.max(freq), ], freq)
}))
rownames(clone.info.multi) <- unique(cell.clone.multi$cloneid)
colnames(clone.info.multi)[colnames(clone.info.multi) == 'Freq'] <- 'size'

clone.info.multi$larry <- NA

for (smpl in multiome.smpls){
  #smpl <- multiome.smpls[2]

  clone.stat <- clone.info.multi[clone.info.multi$sample == smpl, ]$size
  names(clone.stat) <- clone.info.multi[clone.info.multi$sample == smpl, ]$clone

   ##
  ind <- clone.info$id %in% multiome.match[[smpl]]
  head(clone.info[ind,])

  clone.dist <- do.call(rbind, lapply(names(clone.stat), function(clname){
    adist(clname, clone.info[ind,]$clone)
  }))
  rownames(clone.dist) <- names(clone.stat); colnames(clone.dist) <- clone.info[ind,]$clone
  clone.dist.summary <- data.frame(dist = apply(clone.dist, 1, min),
                                   index = apply(clone.dist, 1, which.min),
                                   doublet = apply(clone.dist, 1, function(x) sum(x == min(x))-1 ))

  clone.info.multi[clone.info.multi$sample == smpl, ][clone.dist.summary$dist < 2,]$larry <-
    rownames(clone.info[ind,][clone.dist.summary$index[clone.dist.summary$dist < 2],])

  x1 <- clone.info.multi[clone.info.multi$sample == smpl, ]$size
  names(x1) <- clone.info.multi[clone.info.multi$sample == smpl,]$clone
  #
  x2 <- clone.info[ind, ]$size
  names(x2) <-  clone.info[ind, ]$clone

  nms <- unique(c(names(x1), names(x2)))
  plot(x1[nms], x2[nms], log = 'xy')
  #
}

# plot clone sizes between multiome and RNA
nms <- as.character(unique(c(clone.info.multi$clone, clone.info$clone)))
str(nms)
#
x1 <- rep(0, length(nms)); names(x1) <- nms
x1[as.character(clone.info.multi$clone)] <- clone.info.multi$size
#
x2 <- rep(0, length(nms)); names(x2) <- nms
x2[as.character(clone.info$clone)] <- clone.info$size
#
musample <- rep(NA, length(nms)); names(musample) <- nms
musample[as.character(clone.info.multi$clone)] <- as.character(clone.info.multi$sample)
#
pl <- qplot(x1[nms], x2[nms], log = 'xy', col = musample, size = I(1),
      #geom = c('point', 'smooth'), method = 'lm', linewidth = I(0.1),
      xlab = 'multiome', ylab = 'RNA') +
  theme_bw() + geom_abline(linetype = 2) +
  theme(text = element_text(size = 20), legend.title = element_blank(), legend.position = c(0.8, 0.35))
pl
cor(x1[nms], x2[nms], use = 'na.or.complete')
cor(log(x1[nms]+1), log(x2[nms]+1), use = 'na.or.complete')
cor(x1[nms], x2[nms], use = 'na.or.complete', method = 'spearman')
pdf('figures/clones_RNA_multiome.pdf', width = 4, height = 4)
pl
dev.off()

save(cell.clone.multi, clone.info.multi, file = 'data/processed/clone_multiome_info.RData')




