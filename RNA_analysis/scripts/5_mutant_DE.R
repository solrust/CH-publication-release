source('R/project_gset.R'); source('R/draw_go.R')
library(ggrepel)
library(parallel)
library(Matrix)
library(RColorBrewer)
library(DESeq2)
library(fgsea)
library(gridExtra)

#
load(file = 'data/processed/scrna.RData')
#
hash.annot <- readRDS(file = "data/processed/hashtag_info.rds")
hash.match <- readRDS(file = "data/processed/hashtag_meta.rds")
gene_sets <- readRDS("data/genesets/processed/gene_sets.rds")
#

# estimate cluster-specific count matrices
p2clust1 <- clust.final
clust.states.list <-  lapply(1:length(unique(p2clust1)), function(cl){
  print(cl)
  cls <- names(p2clust1)[p2clust1 == cl]
  str(cls)

  #mat <- aggregate(t(as.matrix(cd[, cls])), by = list(sampleid0[cls]), sum)
  #rownames(mat) <- mat[,1]
  #mat <- mat[,-1]
  #mat <- t(mat)

  mat <- do.call(cbind, lapply(unique(sampleid[cls]), function(smpl){
    print(smpl)
    ind <- sampleid[cls] %in% smpl
    rowSums(cd[, cls][, ind])
  }))
  colnames(mat) <- unique(sampleid[cls])

  return(mat)

})#, mc.cores = length(unique(p2clust1)))
length(clust.states.list)

# estimate cluster-specific DE analysis
de.mut2wt.list <- lapply(1:length(unique(p2clust1)), function(cl){
  print(cl)
  mat <- clust.states.list[[cl]]
  dds_mut <- DESeqDataSetFromMatrix(countData = mat,
                                    colData = hash.match[colnames(mat), , drop = F],
                                    design = ~ exp + genotype)

  dds_mut <- DESeq(dds_mut)
  resultsNames(dds_mut)
  return(dds_mut)
})

de.mut2wt.shrink <- lapply(1:length(unique(p2clust1)), function(cl){
  print(cl)
  dds_mut <- de.mut2wt.list[[cl]]
  res_mut.dnmt3a <- lfcShrink(dds_mut, coef="genotype_DNMT3A_vs_WT",
                              #type="apeglm"
                              type="normal"
  )
  res_mut.tet2 <- lfcShrink(dds_mut, coef="genotype_TET2_vs_WT",
                            #type="apeglm"
                            type="normal"
  )
  list(DNMT3A = res_mut.dnmt3a, TET2 = res_mut.tet2)
})

## run GSEA across cell states
fgsea.dnmt3a.clust <- mclapply(de.mut2wt.shrink, function(deinfo){
  #gscore <- sign(deinfo$DNMT3A$log2FoldChange) * -log(deinfo$DNMT3A$pvalue)
  gscore <- deinfo$DNMT3A$log2FoldChange
  names(gscore) <- rownames(deinfo$DNMT3)
  gscore <- gscore[!is.na(gscore)]
  fgsea(pathways = gene_sets, stats = gscore, minSize  = 15, maxSize  = 500)
}, mc.cores = length(de.mut2wt.shrink))
#
fgsea.tet2.clust <- mclapply(de.mut2wt.shrink, function(deinfo){
  #gscore <- sign(deinfo$TET2$log2FoldChange) * -log(deinfo$TET2$pvalue)
  gscore <- deinfo$TET2$log2FoldChange
  names(gscore) <- rownames(deinfo$TET2)
  gscore <- gscore[!is.na(gscore)]
  fgsea(pathways = gene_sets, stats = gscore, minSize  = 15, maxSize  = 500)
}, mc.cores = length(de.mut2wt.shrink))
#
fgsea.dnmt3a.score <- lapply(fgsea.dnmt3a.clust, function(x){
  y <- -log(x$pval) * sign(x$NES)
  names(y) <- x$pathway
  y
})
gg <- unique(names(unlist(fgsea.dnmt3a.score)))
y <- rep(NA, length(gg)); names(y) <- gg
fgsea.dnmt3a.score <- Reduce(cbind, lapply(fgsea.dnmt3a.score, function(x){
  y[names(x)] <- x
  y
}))
colnames(fgsea.dnmt3a.score) <- 1:ncol(fgsea.dnmt3a.score)
#
fgsea.tet2.score <- lapply(fgsea.tet2.clust, function(x){
  y <- -log(x$pval) * sign(x$NES)
  names(y) <- x$pathway
  y
})
gg <- unique(names(unlist(fgsea.tet2.score)))
y <- rep(NA, length(gg)); names(y) <- gg
fgsea.tet2.score <- Reduce(cbind, lapply(fgsea.tet2.score, function(x){
  y[names(x)] <- x
  y
}))
colnames(fgsea.tet2.score) <- 1:ncol(fgsea.tet2.score)
#
fgsea.all.score <- cbind(fgsea.dnmt3a.score, fgsea.tet2.score)
colnames(fgsea.all.score) <- paste( c(rep('DNMT3A', ncol(fgsea.dnmt3a.score)),
                                      rep('TET2', ncol(fgsea.dnmt3a.score))),
                                    colnames(fgsea.all.score))
fgsea.all.score[is.na(fgsea.all.score)] <- 0

# choose non-redundant pathways:
mainPathways.DNMT3A <- mclapply(1:length(de.mut2wt.shrink), function(cl){
  res <- data.frame(de.mut2wt.shrink[[cl]]$DNMT3A)
  gsc <- res$log2FoldChange; names(gsc) <- rownames(res)
  fgseaRes <- fgsea.dnmt3a.clust[[cl]]
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
                                        gene_sets[nchar(names(gene_sets)) < 60],  gsc[!is.na(gsc)],
                                        pval.threshold = 0.001)
  mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
  return(mainPathways)
}, mc.cores = 10)
#
mainPathways.TET2 <- mclapply(1:length(de.mut2wt.shrink), function(cl){
  res <- data.frame(de.mut2wt.shrink[[cl]]$TET2)
  gsc <- res$log2FoldChange; names(gsc) <- rownames(res)
  fgseaRes <- fgsea.tet2.clust[[cl]]
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
                                        gene_sets[ nchar(names(gene_sets)) < 60],  gsc[!is.na(gsc)],
                                        pval.threshold = 0.001)
  mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
  return(mainPathways)
}, mc.cores = 10)


save(clust.states.list, de.mut2wt.list, de.mut2wt.shrink, fgsea.dnmt3a.clust, fgsea.tet2.clust,
     fgsea.dnmt3a.score, fgsea.tet2.score, fgsea.all.score, mainPathways.DNMT3A, mainPathways.TET2,
     file = 'data/processed/mutant_de.RData')

gosign <- apply(fgsea.all.score, 2, function(x) p.adjust(exp(-abs(x)), 'BH') )
gosign <- apply(gosign, 1, function(x) sum(x < 1e-5, na.rm=T) > 0 )
table(gosign)
#igen <- 2
pcut <- 1e-3
for (igen in 1:2){
fl <- paste('figures/summary_', ifelse(igen == 1, 'DNMT3A', 'TET2'),'.pdf', sep = '')
pdf(fl, width = 12, height = 6)
for(cl in 1:length(de.mut2wt.shrink)){
  print(cl)

  res <- data.frame(de.mut2wt.shrink[[cl]][[igen]])
  lbl <- rownames(res); lbl[res$pvalue > pcut] <- NA
  pl1 <- ggplot(res, aes(log2FoldChange, pvalue,  label = lbl,
                         alpha = I(ifelse(pvalue > pcut, 0.1, 1)),
                         #alpha = I(ifelse(res$padj > pcut, 0.1, 1)),
  )) +
    geom_point(size = I(2)) + labs(x = paste(ifelse(igen == 1, 'DNMT3A', 'TET2'), '/ WT, log'), y = 'p-value') +
    geom_text_repel(size = 4, max.overlaps = 10) +
    geom_point(aes(x = max(res$log2FoldChange,na.rm=T)*0.6, y = 0.4), size = I(20), pch = 21, alpha = I(0.1), fill = as.character(clust.final.cols[cl])) +
    annotate(geom = "text", x = max(res$log2FoldChange,na.rm=T)*0.6, y = 0.4, label = paste(cl), size = 10 ) +
    scale_y_log10(limits = c(min(res$pvalue), 1)) + theme_bw() + theme(text = element_text(size = 20))
  #pl1

  ### show top non-redundant pathways
  #gsc <- res$log2FoldChange; names(gsc) <- rownames(res)
  #collapsedPathways <- collapsePathways(fgsea.dnmt3a.clust[[cl]][order(pval)][padj < 0.01],
  #                                      gene_sets[ nchar(names(gene_sets)) < 60],  gsc[!is.na(gsc)],
  #                                      pval.threshold = 0.001)
  #str(collapsedPathways[[1]])
  #mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
  #nchar(mainPathways)
  #plotGseaTable(gene_sets[mainPathways], gsc[!is.na(gsc)], fgseaRes,
  #              gseaParam = 0.5)

  if (igen == 1){
    mainPathways <- mainPathways.DNMT3A[[cl]]
    fgseaRes <- fgsea.dnmt3a.clust[[cl]]
  }else{
    mainPathways <- mainPathways.TET2[[cl]]
    fgseaRes <- fgsea.tet2.clust[[cl]]
  }

  ind <- fgseaRes$pathway %in% mainPathways
  pl2 <- ggplot(fgseaRes[ind,],
                aes(x = factor(pathway, level = fgseaRes[ind,][order(-sign(NES) * log(pval))]$pathway),
                    y = -log10(pval) * sign(NES), fill = as.factor(sign(NES)) )
  ) +
    ylab('p-value, log') + geom_bar(position = position_dodge(), stat = "identity", col = adjustcolor('black', 1)) +
    scale_fill_manual(values = c('-1' = adjustcolor('blue', 0.2), '1' = adjustcolor('red', 0.2))) +
    coord_flip() + theme_bw() +
    theme(text = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), legend.position = 'none') +
    scale_x_discrete(name = "", position = "top")
  #pl2

  #pdf('/home/solrust/CH/draft/figures/mutant/try.pdf', width = 12, height = 6)
  grid.arrange(arrangeGrob(pl1), arrangeGrob(pl2), widths = c(1,  1))#, top = paste('cluster', cl))
  #dev.off()

}
dev.off()
}


### show some GO enrichment plots
# plot MYC depletion
cl <- 1
gsc <- de.mut2wt.shrink[[cl]]$TET2$log2FoldChange; names(gsc) <- rownames(de.mut2wt.shrink[[cl]]$TET2)
path <- 'M5926: HALLMARK_MYC_TARGETS_V1'; nm <- 'MYC targets'
pl <- plotEnrichment_new(gene_sets[[path]], gsc[!is.na(gsc)], ticksSize = 0.2, spreadES = 0.4, lsize = 2) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) + labs(title = nm) #+ geom_line(aes(size = I(3)))
pl
#
fl <- paste('figures/enrichmnt_hsc_tet2_myc.pdf', sep = '')
pdf(fl, width = 4, height = 4)
pl
dev.off()
# plot translation depletion
cl <- 1
gsc <- de.mut2wt.shrink[[cl]]$TET2$log2FoldChange; names(gsc) <- rownames(de.mut2wt.shrink[[cl]]$TET2)
path <- 'R-HSA-72766: Translation'; nm <- 'Translation'
pl <- plotEnrichment_new(gene_sets[[path]], gsc[!is.na(gsc)], ticksSize = 0.2, spreadES = 0.4, lsize = 2) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) + labs(title = nm) #+ geom_line(aes(size = I(3)))
pl
#
fl <- paste('figures/enrichment_hsc_tet2_translation.pdf', sep = '')
pdf(fl, width = 4, height = 4)
pl
dev.off()
# plot NFkB enrichment
cl <- 1
gsc <- de.mut2wt.shrink[[cl]]$DNMT3A$log2FoldChange; names(gsc) <- rownames(de.mut2wt.shrink[[cl]]$DNMT3A)
path <- "M5890: HALLMARK_TNFA_SIGNALING_VIA_NFKB"; nm <- 'TNFA-NFKB'
pl <- plotEnrichment_new(gene_sets[[path]], gsc[!is.na(gsc)], ticksSize = 0.2, spreadES = 0.4, lsize = 2) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) + labs(title = nm) #+ geom_line(aes(size = I(3)))
pl
#
fl <- paste('figures/enrichment_hsc_dnmt3a_nfkb.pdf', sep = '')
pdf(fl, width = 4, height = 4)
pl
dev.off()


### show heatmap across states and genotypes
p2clust1 <- clust.final
demut <- de.mut2wt.shrink
l2f <- cbind(Reduce('cbind', lapply(demut, function(x) x[[1]]$log2FoldChange)), Reduce('cbind', lapply(demut, function(x) x[[2]]$log2FoldChange)))
rownames(l2f) <- rownames(demut[[1]][[1]])
colnames(l2f) <- 1:(2*length(unique(p2clust1)))
#
l2f.annot <- data.frame(genotype = c(rep('DNMT3A', length(unique(p2clust1))), rep('TET2', length(unique(p2clust1)))),
                        state = factor(rep(1:length(unique(p2clust1)), 2)))
rownames(l2f.annot) <- colnames(l2f)
#
pv <- cbind(Reduce('cbind', lapply(demut, function(x) x[[1]]$padj)), Reduce('cbind', lapply(demut, function(x) x[[2]]$padj)))
pv.min <- apply(pv, 1, function(x) min(x, na.rm = T))

ind <- pv.min < 1e-3 & apply(abs(l2f), 1, function(x) max(x, na.rm = T)) > 0.3
sum(ind)
pheat <- pheatmap::pheatmap(l2f[ind,], breaks = seq(-0.6, 0.6, length.out = 100), show_rownames = T, show_colnames = T,
                            annotation_col = l2f.annot, border_color = NA,#'black',
                            treeheight_row = 50, treeheight_col = 4,
                            annotation_colors = list('genotype' = col.genotype[!names(col.genotype) %in% 'WT'], 'state' = clust.final.cols))

fl <- paste('figures/de_heatmap_names.pdf', sep = '')
pdf(fl, width = 5, height = 50)
pheat
dev.off()
#
pheat <- pheatmap::pheatmap(l2f[ind,], breaks = seq(-0.6, 0.6, length.out = 100), show_rownames = F, show_colnames = F,
                            annotation_col = l2f.annot, border_color = NA,#'black',
                            treeheight_row = 50, treeheight_col = 4,
                            annotation_colors = list('genotype' = col.genotype[!names(col.genotype) %in% 'WT'], 'state' = clust.final.cols))

fl <- paste('figures/de_heatmap.pdf', sep = '')
pdf(fl, width = 5, height = 5)
pheat
dev.off()


### plot clsuter-specific DE scatterplot
detot <- do.call(rbind, lapply(1:10, function(cl){
  data.frame(cbind(de.mut2wt.shrink[[cl]][[1]], de.mut2wt.shrink[[cl]][[2]]), cl)
}))
detot$state <- NA
cutp <- 1e-1
detot$state[detot$padj < cutp & detot$padj.1 < cutp] <- 'both'
detot$state[detot$padj < cutp & detot$padj.1 > cutp] <- 'DNMT3A'
detot$state[detot$padj > cutp & detot$padj.1 < cutp] <- 'TET2'
#
ind <-  (detot$padj < 1e-2 & abs(detot$log2FoldChange) > 0.4) |
  (detot$padj < 1e-2 & abs(detot$log2FoldChange.1) > 0.2)
head(ind)
ind <- ind & !is.na(ind)
head(ind)
length(unique(rownames(detot)[ind]))
detot$genes <- NA
detot$genes[ind] <- rownames(detot)[ind]
sum(!is.na(detot$genes))
detot$genes[is.na(detot$state)] <- NA
sum(!is.na(detot$genes))
#
yy <- data.frame(do.call(rbind, lapply(1:10, function(cl){
  cr <- cor(detot$log2FoldChange[detot$cl == cl], detot$log2FoldChange.1[detot$cl == cl], use = 'na.or.complete')
  cr <- round(cr, 2)
  c(cl = cl, r = cr, x = min(detot$log2FoldChange[detot$cl == cl], na.rm = T) + 0.2,
    y = min(detot$log2FoldChange.1[detot$cl == cl], na.rm = T) + 0.05)
})))
#
detot1 <- detot#[sample(1:nrow(detot), 100000),]
pl <- ggplot() + geom_point(data = detot1, aes(log2FoldChange,log2FoldChange.1, col = state),
                            size = I(1), alpha = I(0.75)) +
  geom_text(data = yy, aes(x = x, y = y,  label = r), size = 10) +
  geom_text_repel(data = detot1, aes(x = log2FoldChange, y = log2FoldChange.1, label = genes)) +  # If you want repelling labels
  labs(x = "DNMT3A, logFC", y = "TET2, logFC", color = 'DE') +
  scale_color_manual(values = c(col.genotype, 'both' = 'black') ) +
  facet_wrap(~ cl, nrow = 2, scales = 'free') + theme_bw()  +
  theme(text = element_text(size = 20), legend.text = element_text(size = 20), axis.title = element_text(size = 30))
pl
#
scl <- 4
pdf('figures/de_scatterplots.pdf', width = scl * 5, height = scl * 2)
pl
dev.off()


dstgo1 <- as.dist(1-cor((fgsea.all.score[gosign,])))
dstgo <- dist((fgsea.all.score[gosign, ]))
hc <- hclust(dstgo)
hclgo <- as.factor(cutree(hc, 6))
table(hclgo)
hclgo <- data.frame("GO clusters" = hclgo)
#require(RColorBrewer)
hclgo.cols <- brewer.pal(length(unique(hclgo$GO.clusters)), "Set2")
names(hclgo.cols) <- unique(hclgo$GO.clusters)
# TODO:
# + change cluster colors
# + remove column dendogram?
# + change cluster names: 'GO term clusters', 'GO clusters'
# play with saved heatmap size
x <- data.frame(clust = rep(factor(1:(ncol(fgsea.all.score)/2)),2),
                geno = unlist(lapply(strsplit(colnames(fgsea.all.score), ' '), function(x)x[1])))
x <- data.frame(do.call(rbind, strsplit(colnames(fgsea.all.score), ' ')))
colnames(x) <- c('genotype', 'clust')
rownames(x) <- colnames(fgsea.all.score)
pheat <- pheatmap::pheatmap( t(fgsea.all.score[gosign,]),
                    clustering_distance_rows = dstgo1,
                    clustering_distance_col = dstgo,
                    breaks = seq(-40, 40, length.out = 100),
                    show_colnames  = F,
                    treeheight_col = 0,
                    annotation_col = hclgo,
                    labels_row = x[,2],
                    annotation_row =  x[,1, drop = F],
                    annotation_names_row = F,
                    annotation_colors = list('GO.clusters' = hclgo.cols,
                                             genotype = col.genotype[-1], clust = clust.final.cols),
                    #annotation_legend = F,
)
pheat

pdf('figures/mutant_states_GO_new.pdf', width = 8, height = 3)
pheat
dev.off()

### plot representative cluster-specific GO terms
pl <- draw_go(matcl = fgsea.all.score[gosign,],
              feature.var0 = rownames(fgsea.all.score)[gosign],
              setclust.go = hclgo,
              cols.go = hclgo.cols,
              n.terms = 15, only.go = FALSE, term.length = 70, cut.pval = 60)
pl

pdf('figures/mutant_GOrepresentatives.pdf',
    width = 12, height = 18)
pl
dev.off()


### download DNMT3A R882 Landau's data
dnmt3a_r882 <- readxl::read_excel("data/genesets/download/Landau_DNMT3A_R882.xlsx", skip = 3, sheet = 3)
dnmt3a_r882 <- as.data.frame(dnmt3a_r882)
head(dnmt3a_r882)
#
dnmt3a_r882 <- dnmt3a_r882[,c(7:10)] # HSPC
#dnmt3a_r882 <- dnmt3a_r882[,c(1:4)] # IMP_ME
#dnmt3a_r882 <- dnmt3a_r882[,c(13:16)] # MEP
#
colnames(dnmt3a_r882) <- c('gene', 'p_perm', 'p_lmm', 'logFC')
dnmt3a_r882 <- dnmt3a_r882[!is.na(dnmt3a_r882$gene),]
rownames(dnmt3a_r882) <- make.unique(dnmt3a_r882$gene)

qplot(dnmt3a_r882$logFC, dnmt3a_r882$p_lmm, log = 'y')
qplot(dnmt3a_r882$p_lmm, bins = 200)
head(dnmt3a_r882[order(dnmt3a_r882$p_lmm),])

dmt3a_this <- de.mut2wt.shrink[[1]]$DNMT3A
gg <- intersect(rownames(dnmt3a_r882), rownames(dmt3a_this))
str(gg)
ind.sign <- p.adjust(dmt3a_this$pvalue, 'BH') < 0.5
table(ind.sign)
dmt3a_this$state <- 'NA'
dmt3a_this[!is.na(ind.sign) & ind.sign & dmt3a_this$log2FoldChange > 0,]$state <- 'up'
dmt3a_this[!is.na(ind.sign) & ind.sign & dmt3a_this$log2FoldChange < 0,]$state <- 'down'
#
gg1 <- gg[dmt3a_this[gg,]$state != 'NA']
y <- dnmt3a_r882[gg1,]$logFC
#y <- sign(dnmt3a_r882[gg1,]$logFC) * -log(dnmt3a_r882[gg1,]$p_perm)
p.show <- wilcox.test(dnmt3a_r882[gg,]$logFC[dmt3a_this[gg,]$state == 'up'], dnmt3a_r882[gg,]$logFC[dmt3a_this[gg,]$state == 'down'])$p.value
pl <- ggplot(,aes(x = as.factor(dmt3a_this[gg1,]$state), y = y,#dnmt3a_r882[gg1,]$logFC,
                  color = as.factor(dmt3a_this[gg1,]$state))) +
  labs(y = 'DNMT3A/WT, HSPC log FC (Nam 2022)', x = ' DNMT3A/WT, HSC DE genes') + theme_bw() +
  #scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_boxplot(width = 0.8, outlier.shape=NA, alpha = 1) +
  geom_point(position=position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8), size = 2, alpha = 0.5)  +
  #geom_text(aes(0.5, 0.8),label = as.character(p.show)) +
  theme(text = element_text(size = 20), legend.position = "none")
pl
p.show
#
pdf('figures/DNMT3A_vs_r882.pdf', width = 5, height = 6)
pl
dev.off()
