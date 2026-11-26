library(DESeq2)
library(openxlsx)
source('R/routines.R')
#
obj.combined <- readRDS('data/processed/multiome_seurat.rds')

obj.combined <- FindClusters(obj.combined, graph.name = "wsnn",
                             cluster.name = 'main_clusters',
                             algorithm = 3, resolution = 1, verbose = FALSE)
table(obj.combined$main_clusters)

Idents(obj.combined) <- obj.combined$projection
DimPlot(obj.combined, reduction = "wnn.umap", shuffle = T,
        #group.by = 'projection'
)

### analyse fate-specific DA peaks
clust.deriv <- c(1, 4, 8, 5, 9, 10, 3)
clust.prog <- c(1, 2, 6)
da_commitment <- mclapply(clust.deriv, function(cl){
  #cl <- 4
  mat1 <- clust.atac.list[[cl]]
  #mat2 <- clust.atac.list[[10]]
  mat2 <- Reduce('+', clust.atac.list[clust.prog])
  mat <- cbind(mat1, mat2)
  inf <- cbind(hash.match[colnames(mat), , drop = F],
               cond = c(rep('clust', ncol(mat1)), rep('back', ncol(mat2))) )
  colnames(mat) <- rownames(inf)
  dds_mut <- DESeqDataSetFromMatrix(countData = mat,
                                    colData = inf,
                                    design = ~ cond + exp)

  dds_mut <- DESeq(dds_mut)
  resultsNames(dds_mut)
  res_mut <- lfcShrink(dds_mut, coef="cond_clust_vs_back",
                       #type="apeglm"
                       type="normal"
  )
  return(res_mut)
}, mc.cores = length(clust.deriv))
names(da_commitment) <- clust.deriv
str(da_commitment)
da_commitment_max <- apply(do.call(cbind, lapply(da_commitment, function(x) x$log2FoldChange)), 1, max)
da_commitment_peaks <- lapply(da_commitment, function(x){
  rownames(x)[!is.na(x$padj) & x$padj < 1e-1 &  x$log2FoldChange >= 1 &
                x$log2FoldChange >= da_commitment_max]
})
len <- unlist(lapply(da_commitment_peaks, length))
len
da_commitment_peaks <- da_commitment_peaks[len > 0]
#
commitment.counts <- do.call(rbind, lapply(da_commitment_peaks, function(x){
  x <- intersect(x, rownames(obj.combined@assays$ATAC$counts))
  colSums(obj.combined@assays$ATAC$counts[x,, drop = F])
}))



### calculate motif-specific matrices
motif.matrix <- obj.combined@assays$ATAC@motifs@data
motif.names <- obj.combined@assays$ATAC@motifs@motif.names
motif.counts <- t(motif.matrix) %*% obj.combined@assays$ATAC$counts
#
xx <- unlist(lapply(strsplit(rownames(obj.combined@assays$ATAC$counts), '-'), function(x)x[1]))
ind <- !xx %in% c(paste('chr', c('M', 'X', 'Y'), sep = ''))
motif.counts <- t(motif.matrix[ind,]) %*% obj.combined@assays$ATAC$counts[ind,]
#
meta.feature <- GetAssayData(obj.combined, assay = "ATAC", slot = "meta.features")
motif.gc <- t(motif.matrix[ind,]) %*% meta.feature$GC.percent[ind]
motif.gc <- motif.gc[,1] / colSums(motif.matrix)
#
open.peaks <- AccessiblePeaks(obj.combined, idents = c(1, 2, 6)) #idents = c(0, 1, 5))
# generate GC-controlled peak permutation
gc.groups <- as.numeric(cut(meta.feature$GC.percent, 30))
qplot(gc.groups, meta.feature$GC.percent)
gcl <- lapply(1:length(unique(gc.groups)), function(i) which(gc.groups == i))
gc.groups.control <- vector()
for (i in 1:length(gc.groups)){
  if (i %% 1000 == 0) print(i)
  j <- gc.groups[i]
  k <- sample(gcl[[j]], 1)
  gc.groups.control[i] <- k
  gcl[[j]] <- setdiff(gcl[[j]], k)
}
# estimate GC-controled background motif-peak matrix
motif.matrix.gc <- motif.matrix[gc.groups.control,]
motif.matrix.gc <- Matrix::Matrix(motif.matrix.gc, sparse = TRUE)
colnames(motif.matrix.gc) <- colnames(motif.matrix)
rownames(motif.matrix.gc) <- rownames(motif.matrix)
#
motif.gc.control <- t(motif.matrix.gc) %*% meta.feature$GC.percent
motif.gc.control <- motif.gc.control[,1] / colSums(motif.matrix.gc)
# estimate GC-controlled ATAC counts
motif.counts.gc <- t(motif.matrix.gc) %*% obj.combined@assays$ATAC$counts
#
qplot(motif.gc.control, motif.gc)+geom_abline()
i <- 1
ggplot() + geom_density(aes(meta.feature$GC.percent[motif.matrix[,i] == 1], col = 'real')) +
  geom_density(aes(meta.feature$GC.percent[motif.matrix.gc[,i] == 1], col = 'control'))


################## run association analysis
geno.names <- c('DNMT3A', 'TET2')
pheno.names <- c('1', '4', '5', '8', '10', '9', 'rate', 'size')
pheno.names.print <- c('HSC bias', 'lymphoid bias', 'dendritic bias', 'myeloid bias', 'erythroid bias', 'mk bias', 'growth rate', 'clone size')
names(pheno.names.print) <- pheno.names
# names of commitment peak modules
fate.names <- c('HSC', 'lymphoid', 'myeloid', 'dendritic', 'Mk', 'erythroid', 'cluster3')
names(fate.names) <- c('1', "4", "8",  "5",  "9",  "10", "3")
#


######################## test peak association with clone phenotypes ############################
size.cut <- 50
use.hspc <- T
clusts <- 1
#if (use.hspc) clusts <- c(0, 1, 5)
if (use.hspc) clusts <- c(1, 2, 6)
#
test.geno <- F
if (test.geno){
  names <- geno.names
  batch_name <- 'id'
}else {
  names <- pheno.names
  batch_name <- 'id'
}

test.modules.pheno <- mclapply(1:length(names), function(i){
  pheno.name <- names[i]
  print(pheno.name)

  cell_conditions <- cells_partition(clone.info.multi, pheno.name = pheno.name, size.cut = size.cut,
                                     batch_name = 'id', test.geno = test.geno)

  da_analysis <- perform_da(obj.combined,
                            count.mat = commitment.counts,
                            do.deseq2 = F,
                            cell_conditions = cell_conditions,
                            cell.info = cell.clone.multi,
                            batch_name = batch_name,
                            #seurat_class = 'main_clusters',
                            seurat_class = 'projection',
                            clusts = clusts
  )
  mat_commitment <- da_analysis$mat
  da_analysis <- perform_da(obj.combined,
                            modality = 'ATAC',
                            do.deseq2 = F,
                            cell_conditions = cell_conditions,
                            cell.info = cell.clone.multi,
                            batch_name = batch_name,
                            seurat_class = 'projection',
                            clusts = clusts
  )
  mat_counts <- da_analysis$mat

  list(mat_commitment = mat_commitment, mat_counts = mat_counts, hash.match1 = da_analysis$hash.match1)
}, mc.cores = length(names))
names(test.modules.pheno) <- names
#
require(glmmTMB)
glm.modules.pheno <- mclapply(1:length(test.modules.pheno), function(i){
  print(i)
  info <- test.modules.pheno[[i]]
  res <- do.call(rbind, lapply(1:nrow(info$mat_commitment), function(j){
    y <- as.numeric(info$mat_commitment[j,])
    size <- colSums(info$mat_counts)
    # perform glmmTMB
    cl <- info$hash.match1$clone
    st <- info$hash.match1$clust
    #ft <- glmmTMB::glmmTMB(y ~ cl + st + offset(log(size)), family = nbinom2(link = 'log'))
    ft <- glmmTMB(y ~ cl + st + offset(log(size)), family = nbinom2(link = 'log'))
    summary(ft)$coefficients$cond[2,]
  }))
  rownames(res) <- rownames(info$mat_commitment)
  res
}, mc.cores = length(test.modules.pheno))

fc <- do.call(rbind, lapply(glm.modules.pheno, function(x) x[,1]))
pv <- do.call(rbind, lapply(glm.modules.pheno, function(x) x[,4]))
rownames(fc) <- rownames(pv) <- names(test.modules.pheno)
rownames(fc) <- rownames(pv)  <- pheno.names.print[rownames(fc)]
colnames(fc) <- colnames(pv)  <- fate.names[colnames(fc)]
fc <- t(fc); pv <- t(pv)
#
pheatmap::pheatmap(fc)
pv

tbl <- prepare_table(fc, pv)
pl <-ggplot(tbl) + theme_bw() + geom_point(aes(module, phenotype, size = abs((value)), stroke = p < 0.05, fill = state, alpha = -log10(p)),
                                     col = 'black', shape = 21)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.title = element_blank(), text = element_text(size = 15)) +
  labs(x = 'peak module', y = 'HSC phenotype') + guides(size = FALSE, alpha = FALSE) +
  scale_fill_manual(values = c('up' = adjustcolor('red', 1), 'down' = adjustcolor('blue', 1) ))
pl
#
pdf('figures/peak_modules_phenotypes_full.pdf', width = 3.5, height = 3.5)
pl
dev.off()

tbl <- prepare_table(fc, pv,
                     row.subset = setdiff(rownames(fc), c('cluster3', 'HSC')),
                     #col.subset = setdiff(colnames(fc), c('mk bias', 'clone size', 'growth rate'))
                     )
pl <-ggplot(tbl) + theme_bw() + geom_point(aes(module, phenotype, size = abs((value))^1, stroke = p < 0.05, fill = state, alpha = -log10(p)),
                                           col = 'black', shape = 21)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.title = element_blank(), text = element_text(size = 15)) +
  labs(x = 'peak module', y = 'HSC phenotype') + guides(size = FALSE, alpha = FALSE) +
  scale_fill_manual(values = c('up' = adjustcolor('red', 1), 'down' = adjustcolor('blue', 1) ))
pl
pdf('figures/peak_modules_phenotypes.pdf', width = 3.0, height = 3.0)
pl
dev.off()



####################### test peak association with mutant genotypes ############################
size.cut <- 50
use.hspc <- F
clusts <- 1
#if (use.hspc) clusts <- c(0, 1, 5)
if (use.hspc) clusts <- c(1, 2, 6)
#
test.geno <- TRUE
if (test.geno){
  names <- geno.names
  batch_name <- 'id'
}else {
  names <- pheno.names
  batch_name <- 'id'
}

test.modules.geno <- mclapply(1:length(names), function(i){
  pheno.name <- names[i]
  print(pheno.name)

  cell_conditions <- cells_partition(clone.info.multi, pheno.name = pheno.name, size.cut = size.cut,
                                     batch_name = batch_name, test.geno = test.geno)

  da_analysis <- perform_da(obj.combined,
                            #modality = 'ATAC',
                            count.mat = commitment.counts,
                            do.deseq2 = F,
                            cell_conditions = cell_conditions,
                            cell.info = cell.clone.multi,
                            batch_name = batch_name,
                            #seurat_class = 'main_clusters',
                            seurat_class = 'projection',
                            clusts = clusts
  )
  mat_commitment <- da_analysis$mat
  da_analysis <- perform_da(obj.combined,
                            modality = 'ATAC',
                            do.deseq2 = F,
                            cell_conditions = cell_conditions,
                            cell.info = cell.clone.multi,
                            batch_name = batch_name,
                            #seurat_class = 'main_clusters',
                            seurat_class = 'projection',
                            clusts = clusts
  )
  mat_counts <- da_analysis$mat

  list(mat_commitment = mat_commitment, mat_counts = mat_counts, hash.match1 = da_analysis$hash.match1)
}, mc.cores = length(names))
names(test.modules.geno) <- names
#
glm.modules.geno <- mclapply(1:length(test.modules.geno), function(i){
  print(i)
  info <- test.modules.geno[[i]]
  res <- do.call(rbind, lapply(1:nrow(info$mat_commitment), function(j){
    y <- as.numeric(info$mat_commitment[j,])
    size <- colSums(info$mat_counts)
    #size <- colSums(info$mat_commitment)
    # perform glmmTMB
    cl <- info$hash.match1$clone
    st <- info$hash.match1$clust
    ft <- glmmTMB(y ~ cl + #st +
                             offset(log(size)), family=nbinom2(link = 'log'))
    summary(ft)$coefficients$cond[2,]
  }))
  rownames(res) <- rownames(info$mat_commitment)
  res
}, mc.cores = length(test.modules.geno))
fc <- do.call(rbind, lapply(glm.modules.geno, function(x) x[,1]))
pv <- do.call(rbind, lapply(glm.modules.geno, function(x) x[,4]))
rownames(fc) <- rownames(pv) <- names(test.modules.geno)
colnames(fc) <- colnames(pv)  <- fate.names[colnames(fc)]
fc <- t(fc); pv <- t(pv)
#
pheatmap::pheatmap(fc)
fc
pv

tbl <- prepare_table(fc, pv)
pl <- ggplot(tbl) + theme_bw() + geom_point(aes(module, phenotype, size = abs((value)), stroke = p < 0.05, fill = state, alpha = -log10(p)),
                                           col = 'black', shape = 21)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.title = element_blank(), text = element_text(size = 15)) +
  labs(x = 'peak module', y = '') + guides(size = FALSE, alpha = FALSE) +
  scale_size_continuous(range = c(0, 10))  +
  scale_fill_manual(values = c('up' = adjustcolor('red', 1), 'down' = adjustcolor('blue', 1) ))
pl
#
pdf('figures/peak_modules_genotypes_norm_all.pdf', width = 6, height = 3)
pl
dev.off()


###################### test motif peaks association with clone phenotypes #######################
size.cut <- 50
use.hspc <- T
clusts <- 1
#if (use.hspc) clusts <- c(0, 1, 5)
if (use.hspc) clusts <- c(1, 2, 6)
#
test.geno <- F
if (test.geno){
  names <- geno.names
  batch_name <- 'id'
}else {
  names <- pheno.names
  batch_name <- 'id'
}

test.motifs.pheno <- mclapply(1:length(names), function(i){
  pheno.name <- names[i]
  print(pheno.name)

  cell_conditions <- cells_partition(clone.info.multi, pheno.name = pheno.name, size.cut = size.cut,
                                     batch_name = batch_name, test.geno = test.geno)

  da_analysis <- perform_da(obj.combined,
                            #modality = 'ATAC',
                            count.mat = motif.counts,
                            #do.deseq2 = T,
                            cell_conditions = cell_conditions,
                            cell.info = cell.clone.multi,
                            batch_name = batch_name,
                            #seurat_class = 'main_clusters',
                            seurat_class = 'projection',
                            clusts = clusts
  )

  list(cell_conditions = cell_conditions, da_analysis = da_analysis)
}, mc.cores = length(names))
names(test.motifs.pheno) <- names
#

pdf('figures/motif_phenotpyes_volcano.pdf')
for (i in which(names(test.motifs.pheno) %in% names(pheno.names.print))){
  print(i)
  res <- test.motifs.pheno[[i]]$da_analysis$res
  rownames(res) <- as.factor(as.character(unlist(motif.names)))
  pl <- qplot(res$log2FoldChange, res$pvalue, #col = motif.gc,
              log = 'y', label = rownames(res), na.rm=T) + geom_text_repel() +
    theme_bw() + labs(title  = paste(pheno.names.print[names(test.motifs.pheno)[i]]),
                      x = 'logFC', y = 'p-value') + theme(text = element_text(size = 20))
    #geom_hline(yintercept = min(res$pvalue[ res$padj > 0.2]), linetype = 2)#+
    #scale_color_gradientn(colors = c("blue", "darkred"), values = scales::rescale(c(min(motif.gc), mean(motif.gc), max(motif.gc))))
  print(pl)
}
dev.off()


fc <- do.call(rbind, lapply(test.motifs.pheno, function(x) x$da_analysis$res$log2FoldChange))
pv <- do.call(rbind, lapply(test.motifs.pheno, function(x)  x$da_analysis$res$pvalue))
padj <- do.call(rbind, lapply(test.motifs.pheno, function(x)  x$da_analysis$res$padj))
rownames(fc) <- rownames(pv)<- rownames(padj) <- names(test.motifs.pheno)
colnames(fc) <- colnames(pv) <- colnames(padj)  <- unlist(motif.names)[rownames(motif.counts)]
fc <- t(fc); pv <- t(pv); padj <- t(padj);
#
nms <- intersect(names(pheno.names.print), colnames(fc))
fc <- fc[, nms]; pv <- pv[,nms]; padj <- padj[,nms]
colnames(fc) <- colnames(pv) <- colnames(padj) <- pheno.names.print[nms]
#
pheatmap::pheatmap(fc[apply(padj,1,min) < 3e-1, apply(padj,2,min) < 3e-1])


ind <- apply(padj,1,min) < 2e-1
nms <- setdiff(colnames(fc), 'mk bias')
tbl <- prepare_table(fc[ind, nms], pv[ind, nms])
pl <- ggplot(tbl) + theme_bw() + geom_point(aes(module, phenotype, size = abs((value)),
      stroke = p < 0.01, fill = state, alpha = -log10(p)), col = 'black', shape = 21)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
        legend.title = element_blank(), text = element_text(size = 20)) +
  labs(x = 'TF motif', y = 'HSC phenotype' ) + guides(size = FALSE, alpha = FALSE) +
  scale_size_continuous(range = c(0, 8))  +
  scale_fill_manual(values = c('up' = adjustcolor('red', 1), 'down' = adjustcolor('blue', 1) ))
pl
#
pdf('figures/motif_phenotypes_bubble.pdf', width = 10, height = 5)
print(pl)
dev.off()

# save data tables
nms <- intersect(names(pheno.names.print), colnames(fc))
ll <- lapply(1:length(test.motifs.pheno), function(i){
  ord <- order(test.motifs.pheno[[i]]$da_analysis$res$pvalue)
  res <- test.motifs.pheno[[i]]$da_analysis$res[ord,]
  rownames(res) <- motif.names[rownames(res)]
  res <- data.frame(res)
  return(res)
})
names(ll) <- pheno.names.print[names(test.motifs.pheno)]
write.xlsx(ll, '../resource/motif_associations.xlsx', row.names = TRUE, col.names = TRUE)





###################### test motif peaks association with genotypes ##################
size.cut <- 50
use.hspc <- T#TRUE #T
clusts <- 1
#if (use.hspc) clusts <- c(0, 1, 5)
if (use.hspc) clusts <- c(1, 2, 6)
#
test.geno <- T
if (test.geno){
  names <- geno.names
  batch_name <- 'id'
}else {
  names <- pheno.names
  batch_name <- 'id'
}


myc.ct <- colSums(obj.combined@assays$ATAC$counts[names(myc_bound_peaks), ])
motif.counts.ext <- rbind(motif.counts,
                          'MYC' = myc.ct
                          )

motif.counts.ext1 <- t(t(motif.counts.ext) / colSums(motif.counts.ext))
xx <- cor(as.matrix(t(motif.counts.ext1)), as.matrix(t(motif.counts.ext1)))
pheatmap::pheatmap(xx)

test.motifs.geno <- mclapply(1:length(names), function(i){
  pheno.name <- names[i]
  print(pheno.name)

  cell_conditions <- cells_partition(clone.info.multi, pheno.name = pheno.name, size.cut = size.cut,
                                     batch_name = batch_name, test.geno = test.geno)

  da_analysis <- perform_da(obj.combined,
                            #modality = 'ATAC',
                            count.mat = motif.counts.ext,
                            #do.deseq2 = T,
                            cell_conditions = cell_conditions,
                            cell.info = cell.clone.multi,
                            batch_name = batch_name,
                            #seurat_class = 'main_clusters', #'projection',
                            seurat_class = 'projection',
                            clusts = clusts
  )

  run_control <- FALSE
  if (run_control){
  da_analysis_gc <- perform_da(obj.combined,
                            #modality = 'ATAC',
                            count.mat = motif.counts.gc,
                            #do.deseq2 = T,
                            cell_conditions = cell_conditions,
                            cell.info = cell.clone.multi,
                            batch_name = batch_name,
                            #seurat_class = 'main_clusters', #'projection',
                            seurat_class = 'projection',
                            clusts = clusts
  )
  }else{
    da_analysis_gc <- NULL
  }

  list(cell_conditions = cell_conditions, da_analysis = da_analysis, da_analysis_gc = da_analysis_gc)
}, mc.cores = length(names))
names(test.motifs.geno) <- names
#

pdf('figures/motif_genotypes_volcano.pdf')
for (i in 1:2){
  print(i)
  res <- test.motifs.geno[[i]]$da_analysis$res
  rownames(res) <- c(as.factor(as.character(unlist(motif.names))), 'MYC_Chip')
  pl <- qplot(res$log2FoldChange, res$pvalue, #col = motif.gc,
              log = 'y', label = rownames(res), na.rm=T) + geom_text_repel() +
    labs(x = 'logFC', y = 'p-value', title  = names(test.motifs.geno)[i]) + theme_bw() +
    geom_hline(yintercept = ifelse(sum(res$padj <= 0.2)>1, min(res$pvalue[ res$padj > 0.2]), NA), linetype = 2) +
  scale_color_gradientn(colors = c("blue",  'grey70', "darkred"), values = scales::rescale(c(min(motif.gc), mean(motif.gc), max(motif.gc))))
  print(pl)
}
dev.off()

lfc1 <- test.motifs.geno$DNMT3A$da_analysis$res$log2FoldChange
lfc2 <- test.motifs.geno$TET2$da_analysis$res$log2FoldChange
nms <- unlist(motif.names)[rownames(test.motifs.geno$DNMT3A$da_analysis$res)]
pl <- ggplot(,aes(x = lfc1, y = lfc2, label = nms,  col = motif.gc)) + geom_point() +
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + geom_text_repel() +
  labs(x = 'logFC, DNMT3A', y = 'logFC, TET2', title  = 'DNMT3A vs TET2') + theme_bw() +
  scale_color_gradientn(colors = c("blue", 'grey70' ,"darkred"), values = scales::rescale(c(min(motif.gc), mean(motif.gc), max(motif.gc))))
pl
#
pdf('figures/motif_genotypes_comparison.pdf')
pl
dev.off()

# save data tables
ll <- lapply(1:length(test.motifs.geno), function(i){
  ord <- order(test.motifs.geno[[i]]$da_analysis$res$pvalue)
  res <- test.motifs.geno[[i]]$da_analysis$res[ord,]
  rownames(res) <- motif.names[rownames(res)]
  res <- data.frame(res)
  return(res)
})
names(ll) <- names(test.motifs.geno)
write.xlsx(ll, '../resource/motif_associations_genotypes.xlsx', row.names = TRUE, col.names = TRUE)

x <- read.xlsx('../resource/motif_associations_genotypes.xlsx', sheet = 1)
head(x)
y <- read.xlsx('../resource/motif_associations_genotypes.xlsx', sheet = 2)
head(y)
#
rownames(x) <- x[[1]]
rownames(y) <- y[[1]]
tfs.show <- c('HLF', 'ID1', 'ID3', 'ID4', 'TEF', 'BARHL2', 'PAX7')
tfs.show <- tfs.show[tfs.show %in% rownames(x)]
qplot(x[tfs.show,]$log2FoldChange, y[tfs.show,]$log2FoldChange)
x[tfs.show, ]
y[tfs.show, ]

#
glm.modules.geno <- lapply(1:length(test.motifs.geno), function(i){
  print(i)
  info <- test.motifs.geno[[i]]
  res <- do.call(rbind, mclapply(1:nrow(info$da_analysis$mat), function(j){
    if(j %% 100 == 0) print(j)
    y <- as.numeric(info$da_analysis$mat[j,])
    size <- as.numeric(info$da_analysis_gc$mat[j,])
    # perform glmmTMB
    cl <- info$da_analysis$hash.match1$clone
    st <- info$da_analysis$hash.match1$clust
    ft <- glmmTMB::glmmTMB(y ~ cl + st + offset(log(size)), family=nbinom2(link = 'log'))
    summary(ft)$coefficients$cond[2,]
  }, mc.cores = 20))
  rownames(res) <- rownames(info$da_analysis$mat)
  res
})
names(glm.modules.geno) <- names(test.motifs.geno)
fc <- do.call(rbind, lapply(glm.modules.geno, function(x) x[,1]))
pv <- do.call(rbind, lapply(glm.modules.geno, function(x) x[,4]))
rownames(fc) <- rownames(pv) <- names(test.motifs.geno)
colnames(fc) <- colnames(pv)  <- unlist(motif.names)[rownames(glm.modules.geno[[1]])]
fc <- t(fc); pv <- t(pv)
#
ind <- p.adjust(pv[,1], 'BH') < 0.1
qplot(fc[,1], pv[,1], log = 'y', col = ind, label = rownames(fc)) + geom_text_repel()
ind <- p.adjust(pv[,2], 'BH') < 0.2
qplot(fc[,2], pv[,2], log = 'y', col = ind, label = rownames(fc)) + geom_text_repel()
qplot(fc[,1], fc[,2], label = rownames(fc)) + geom_text_repel()
#
pheatmap::pheatmap(fc[apply(pv,1,min) < 1e-2, ])
pv

# save data tables
ll <- lapply(1:length(glm.modules.geno), function(i){
  ord <- order(test.motifs.geno[[i]]$da_analysis$res$pvalue)
  res <- test.motifs.geno[[i]]$da_analysis$res[ord,]
  rownames(res) <- motif.names[rownames(res)]
  res <- data.frame(res)
  return(res)
})
names(ll) <- names(test.motifs.geno)
write.xlsx(ll, '../resource/motif_associations_genotypes.xlsx', row.names = TRUE, col.names = TRUE)



pdf('figures/motif_genotypes_gcontrol_volcano1.pdf')
for (i in 1:2){
  print(i)
  res <- glm.modules.geno[[i]]
  rownames(res) <- as.factor(as.character(unlist(motif.names)))
  pl <- qplot(res[,1], res[,4], col = motif.gc,
              log = 'y', label = rownames(res), na.rm=T) + geom_text_repel() +
    labs(x = 'logFC', y = 'p-value', title  = names(glm.modules.geno)[i]) + theme_bw() +
    scale_color_gradientn(colors = c("blue",  'grey70', "darkred"), values = scales::rescale(c(min(motif.gc), mean(motif.gc), max(motif.gc))))
  print(pl)
}
dev.off()
#
lfc1 <- glm.modules.geno[[1]][,1]
lfc2 <- glm.modules.geno[[2]][,1]
nms <- unlist(motif.names)[rownames(glm.modules.geno[[1]])]
pl <- ggplot(,aes(x = lfc1, y = lfc2, label = nms,  col = motif.gc)) + geom_point() +
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + geom_text_repel() +
  labs(x = 'logFC, DNMT3A', y = 'logFC, TET2', title  = 'DNMT3A vs TET2') + theme_bw() +
  scale_color_gradientn(colors = c("blue", 'grey70' ,"darkred"), values = scales::rescale(c(min(motif.gc), mean(motif.gc), max(motif.gc))))
pl
#
pdf('figures/motif_genotypes_gcontrol_comparison1.pdf')
pl
dev.off()





