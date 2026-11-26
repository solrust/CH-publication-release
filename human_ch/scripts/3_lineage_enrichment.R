hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
old.scrna <- new.env()
load(file = '../RNA_analysis/data/processed/scrna.RData', envir = old.scrna)

obj.human <- readRDS('data/human_ch.rds')
obj.human <- readRDS('data/human_ch_seurat.rds')

do.joint_emb <- TRUE
if (do.joint_emb){
  metadata_integrated <- readRDS(file = 'data/integrated_seurat_metadata.rds')
  emb.joint <- metadata_integrated[, c('umaprpca_1', 'umaprpca_2')]
}

do.CH <- FALSE
col.genotype <- old.scrna$col.genotype
if (do.CH){
  # CH data
  pca.reduction <- old.scrna$p2pca
  emb <- old.scrna$emb
  if (do.joint_emb) emb <- emb.joint[rownames(emb.joint) %in% rownames(emb),]
  sampleid <- old.scrna$sampleid
  genoid <- hash.match[old.scrna$sampleid, ]$genotype; names(genoid) <- names(sampleid)
  expid = hash.match[old.scrna$sampleid, ]$exp; names(expid) <- names(sampleid)
  metadat <- unique(data.frame(sampleid, genoid)); rownames(metadat) <- metadat$sampleid
  clustid <- old.scrna$clust.final[names(sampleid)]
  clust.hsc2gmp <- c(8,7,6,1)
  clust.hsc <- c(1)
  clust.cols <- old.scrna$clust.final.cols
}else{
  # Vyas data
  pca.reduction <- obj.human@reductions$pca@cell.embeddings
  emb <- obj.human@reductions$umap.pca@cell.embeddings
  if (do.joint_emb) emb <- emb.joint[rownames(emb.joint) %in% rownames(emb),]
  sampleid <- paste(obj.human$Sample, obj.human$Genotype, sep = '_'); names(sampleid) <- names(obj.human$Sample)
  sampleid <- sampleid[obj.human$Genotype %in% c('DNMT3A', 'TET2', 'WT-CH', 'WT-control')]
  genoid <- obj.human$Genotype[names(sampleid)]
  expid <- NA
  metadat <- unique(data.frame(sampleid, genoid)); rownames(metadat) <- metadat$sampleid
  clustid <- obj.human$Cluster[names(sampleid)]
  clust.hsc2gmp <- c('HSC/MPP', 'LMPP', 'Early GMP', 'Late GMP')
  clust.hsc <- c('HSC/MPP')
  clust.cols <- 1:length(clust.hsc2gmp); names(clust.cols) <- clust.hsc2gmp
  #
  names(col.genotype) <-  gsub("^WT$", "WT-CH", names(col.genotype))
}

pr <- pcaMethods::pca(pca.reduction[names(clustid)[clustid %in% clust.hsc2gmp], ], nPcs = 3)
plot(pr@scores[,1], pr@scores[,2])
score <- pr@scores[,1]
score0 <- score
cells.show <- intersect(rownames(emb), names(score))
qplot(emb[cells.show,1], emb[cells.show,2], size = I(0.03),
      col = rank(score[cells.show])) + scale_color_gradient2(low = "blue",  mid = "white",  high = "red",
                                                         midpoint = mean(rank(score[cells.show])))

datx0 <- data.frame(score = rank(score), sid = sampleid[names(score)],
                   genotype = genoid[names(score)],
                   exp = expid[names(score)])
head(datx0)

# select samples
do.tet2 <- FALSE
if (do.CH) {
  if (do.tet2){
    datx <- datx0[datx0$exp == '127',]
  }else{
    datx <- datx0[datx0$exp == '123',]
  }
}else{
  # select TET-CH samples
  gn <- 'DNMT3A'; gn.remove <- 'TET2';
  if (do.tet2) {gn <- 'TET2'; gn.remove <- 'DNMT3A'}
  sid.tet2 <- names(table(obj.human$Sample[which(obj.human$Genotype == gn)]))
  datx <- datx0
  datx <- datx[!datx$genotype %in% gn.remove, ]
}


#TODO: use only WT for quantiles
do.intra <- FALSE
if (do.CH) {
  ref <- 'WT'; step <- 0.02
}else {
  if (do.intra){
    ref <- c('WT-CH'); step <- 0.1
    sel <- rep(FALSE, nrow(datx))
    for (i in 1:length(sid.tet2)) sel <- sel + grepl(sid.tet2[i], datx$sid)
    datx <- datx[sel==1, ]
  }else{
    ref <- c('WT-control'); step <- 0.1

    col.genotype <- c(col.genotype, 'WT-control' = 'lightblue')
    #sel <- rep(FALSE, nrow(datx))
    #for (i in 1:length(sid.tet2)) sel <- sel + grepl(sid.tet2[i], datx$sid)
    #datx <- datx[sel==1, ]
  }
}
score_quant <- as.numeric(quantile(datx$score[datx$genotype %in% ref],
                                   probs = seq(0, 1, step)))
pt.density <- do.call(rbind, lapply(1:(length(score_quant)-1), function(i){
  ind <- datx$score >= score_quant[i] & datx$score < score_quant[i+1]
  clust.dist <- table(clustid[rownames(datx)[ind]])
  cl <- names(clust.dist)[which.max(clust.dist)]
  sid.stat <- table(datx[ind,]$sid)
  data.frame(bin = i, sid = names(sid.stat), clust = cl, stat = as.numeric(sid.stat))
}))
pt.density$genotype <- as.character(metadat[pt.density$sid, ]$genoid)
head(pt.density)
pt.density.genotype <- aggregate(stat ~ genotype + bin, data = pt.density, FUN = sum)
head(pt.density.genotype)
# genotype-level density

# draw by normalizing to HSC density
cells.hsc <- names(clustid)[which(clustid %in% clust.hsc)]
#cells.hsc <- names(clustid)[which(clustid %in% clust.hsc2gmp)]
sid.hsc <- table(sampleid[cells.hsc])
pt.density$stat_norm <- pt.density$stat / as.numeric(sid.hsc[pt.density$sid])

# genotype-level normalized HSC density
sid.names <- unique(pt.density$sid)
hsc.geno.stat <- data.frame(genotype = metadat[sid.names,]$genoid,#hash.match[sid.names, 'genotype'],
                            geno.hsc = as.numeric(sid.hsc[sid.names]))
genoid.hsc <- aggregate(geno.hsc ~ genotype, data = hsc.geno.stat, FUN = sum);
genoid.hsc1 <- genoid.hsc[,2]; names(genoid.hsc1) <- genoid.hsc[,1]
pt.density.genotype$stat_norm <- pt.density.genotype$stat / genoid.hsc1[as.character(pt.density.genotype$genotype)]

ggplot(pt.density) + geom_point(aes(bin, stat_norm,  group = sid, col = genotype), alpha = I(0.3))

span.param <- 1
if (do.CH) span.param <- 0.3
#span.param <- 0.3 # 1 for human CH #0.3 for xeno CH
pl <- ggplot(pt.density) +  geom_point(aes(bin, stat_norm,  group = sid, col = genotype), alpha = I(0.3)) +
  #geom_point(aes(bin, rep(0, nrow(pt.density))), col = clust.cols[pt.density$clust], size = I(3), alpha = I(0.2)) +
  geom_line(aes(bin, stat_norm, group = sid, col = genotype), alpha = I(0.3)) +
  geom_smooth(aes(bin, stat_norm,  group = sid, col = genotype),span = 0.3, linewidth = 0.5,
              se = FALSE) + xlab('pseudotime') + ylab('density to HSC') +
  geom_smooth(data = pt.density.genotype, aes(bin, stat_norm, col = genotype),
              span = span.param, linewidth = 3, fill = 'grey80', alpha = (0.5)) + theme_bw() +
  scale_color_manual(values = col.genotype) +
  theme(axis.title = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 20),
        #legend.position = c(0.05, 0.05), legend.justification = c(-0.1, -0.1))
        legend.position = c(0.1, 0.65), legend.justification = c(0, 0))
pl

pt <- rank(score)
cells.show <- intersect(rownames(emb), names(score))
str(cells.show)
pl_emb <- ggplot() + geom_point(aes(emb[1:10000,1], emb[1:10000,2]), col = 'grey20', alpha = I(0.1), size = I(0.05)) +
                    geom_point(aes(emb[cells.show,1], emb[cells.show,2], col = pt[cells.show]), size = I(0.05)) +
  scale_color_gradient2(low = "blue",  mid = "white",  high = "red", midpoint = mean(pt[cells.show])) +
  xlab('') + ylab('') + theme_bw() + theme(legend.position = 'none')
pl_emb
#
require(gridExtra)
grid.arrange(pl_emb, pl, nrow = 1,  ncol = 2)
#

do.plot <- TRUE
if (do.plot){
  pref <- 'CH'; if (!do.CH) pref <- 'human'
  gn <- 'TET2'; if (!do.tet2) gn <- 'DNMT3A'
  comparison <- ''
  if (do.intra == FALSE) {comparison <- 'inter'}
  flname <- paste('figures/', pref, '_', gn,'_', comparison,'_pt_bias_new.pdf', sep = '')
  pdf(flname, width = 8, height = 4)
  grid.arrange(pl_emb, pl, nrow = 1,  ncol = 2)
  dev.off()
}
