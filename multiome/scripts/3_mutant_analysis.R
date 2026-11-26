library(DESeq2)
library(ggrepel)

hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
load('../RNA_clone_analysis/data/processed/clone_multiome_info.RData')
de.env <- new.env()
load(file = '../RNA_analysis/data/processed/mutant_de.RData', envir = de.env)
env.dna <- new.env()
load('../DNA_analysis/data/processed/DNA_clones.RData', env = env.dna)
obj.combined <- readRDS(file = 'data/processed/multiome_seurat.rds')

rownames(cell.clone.multi) <- paste(rownames(cell.clone.multi), '1', sep = '-')
# retain rate information for clones detected in multiome
x <- intersect(names(env.dna$rate), clone.info.multi$larry)
rate <- env.dna$rate[x]
clons <- clone.info.multi$larry
summary.size <- aggregate(clone.info[clons,]$size, by = list(clone.info[clons,]$id), sum)
rownames(summary.size) <- summary.size[,1]; summary.size <- summary.size[,-1,drop = F]
#
clone.info.multi <- cbind(clone.info.multi, clone.info[clons, c('id', 'genotype')],
                          'rate' = rate[clons],
                          'size' = log(clone.info[clons,]$size / summary.size[clone.info[clons,]$id,]),
                          clone.states2[clons, c('1', '4', '5', '8', '9', '10')])
cell.clone.multi <- cbind(cell.clone.multi, clone.info.multi[cell.clone.multi$cloneid,] )
# add expgeno
ids <- cell.clone.multi$id
cell.clone.multi$expgeno <- paste(hash.match[ids,]$exp, hash.match[ids,]$genotype, sep = '_')
cell.clone.multi$expgeno[grepl('NA', cell.clone.multi$expgeno)] <- NA
head(cell.clone.multi)

### construct and save a cell metadata file
x <- cell.clone.multi[, c('sample', 'cloneid', 'larry', 'id', 'genotype')]
x <- x[intersect(rownames(x), colnames(obj.combined@assays$ATAC$counts)), ]
cells <- colnames(obj.combined@assays$ATAC$counts)
cells.add <- cells[!cells%in%rownames(x)]
new_rows <- as.data.frame(matrix(NA, ncol = ncol(x), nrow = length(cells.add)))
colnames(new_rows) <- colnames(x)
rownames(new_rows) <- cells.add
new_rows$sample <- unlist(lapply(strsplit(cells.add, '_'), function(x) x[1]))
new_rows$genotype[new_rows$sample == 'C1'] <- 'WT'
new_rows$genotype[new_rows$sample == 'D'] <- 'DNMT3A'
new_rows$genotype[new_rows$sample == 'T'] <- 'TET2'
x <- rbind(x, new_rows)
x <- x[cells,]
x <- cbind(x, projection = obj.combined$projection[rownames(x)], seurat_clusters = obj.combined$seurat_clusters[rownames(x)])
x <- cbind(x, obj.combined@reductions$wnn.umap@cell.embeddings[rownames(x),])
write.table(x, file = "data/processed/metainfo.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = '\t')


### partition clusters in WT/Mut:
geno <- as.character(cell.clone.multi$genotype); names(geno) <- rownames(cell.clone.multi)
geno[!is.na(geno) & geno == 'WT'] <- paste(cell.clone.multi[!is.na(geno) & geno == 'WT', ]$genotype,
                                           cell.clone.multi[!is.na(geno) & geno == 'WT', ]$sample, sep = '_')
str(geno)
obj.combined$geno <- geno
obj.combined$projmut <- paste(obj.combined$geno, obj.combined$projection, sep = '_')
obj.combined$projmut[grepl('NA', obj.combined$projmut)] <- NA


# estimate cluster-specific count matrices
p2clust1 <- obj.combined$projection
#p2clust1 <- obj.combined$main_clusters
n.clust <- length(unique(p2clust1))-1
clust.ids <- 1:n.clust# - 1
count.mat <- cbind(obj.combined@assays$RNA$counts.1, obj.combined@assays$RNA$counts.2, obj.combined@assays$RNA$counts.3, obj.combined@assays$RNA$counts.4)
clust.rna.list <-  lapply(clust.ids, function(cl){
  print(cl)
  cls <- names(p2clust1)[!is.na(p2clust1) & p2clust1 == cl]
  print(paste('cells:',length(cls)))
  #table(cell.clone.multi[cls,]$id)

  smpls <- unique(cell.clone.multi[cls,]$id)
  smpls <- smpls[!is.na(smpls)]

  #smpls <- unique(cell.clone.multi[cls,]$sample)
  #smpls <- smpls[!is.na(smpls)]

  mat <- do.call(cbind, lapply(smpls, function(smpl){
    #print(smpl)
    ind <- cell.clone.multi[cls,]$id %in% smpl
    #ind <- cell.clone.multi[cls,]$sample %in% smpl
    if (sum(ind) > 0){
      rowSums(count.mat[, cls][, ind, drop = F])
    }else{
      rep(0, nrow(count.mat))
    }
  }))
  colnames(mat) <- smpls

  return(mat)

})#, mc.cores = length(unique(p2clust1)))
#
count.mat <- obj.combined@assays$ATAC$counts
#count.mat <- cbind(obj.combined@assays$RNA$counts.1, obj.combined@assays$RNA$counts.2, obj.combined@assays$RNA$counts.3, obj.combined@assays$RNA$counts.4)
clust.atac.list <-  lapply(clust.ids, function(cl){
  print(cl)
  cls <- names(p2clust1)[!is.na(p2clust1) & p2clust1 == cl]
  print(paste('cells:',length(cls)))
  #table(cell.clone.multi[cls,]$id)

  smpls <- unique(cell.clone.multi[cls,]$id)
  smpls <- smpls[!is.na(smpls)]

  #smpls <- unique(cell.clone.multi[cls,]$sample)
  #smpls <- smpls[!is.na(smpls)]

  mat <- do.call(cbind, lapply(smpls, function(smpl){
    #print(smpl)
    ind <- cell.clone.multi[cls,]$id %in% smpl
    #ind <- cell.clone.multi[cls,]$sample %in% smpl
    if (sum(ind) > 0){
      rowSums(count.mat[, cls][, ind, drop = F])
    }else{
      rep(0, nrow(count.mat))
    }
  }))
  colnames(mat) <- smpls

  return(mat)

})#, mc.cores = length(unique(p2clust1)))
length(clust.atac.list)



### run joint DE analysis across progenitor clusters
mat <- do.call(cbind, clust.rna.list)
clid <- unlist(lapply(1:length(clust.rna.list), function(i) rep(i, ncol(clust.rna.list[[i]])) ))
hash.match1 <- hash.match[colnames(mat), , drop = F]
hash.match1$clid <- clid
colnames(mat) <- rownames(hash.match1)
# DE for dnmt3a
ind1 <- hash.match1$exp == '123' #& hash.match1$clid %in% c(1, 2)#, 6, 7)
dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                  design = ~ genotype + clid)
dds_mut <- DESeq(dds_mut)
resultsNames(dds_mut)
res_rna_dnmt3a <- lfcShrink(dds_mut, coef="genotype_DNMT3A_vs_WT", type="normal")
sum(res_rna_dnmt3a$padj < 0.1, na.rm = T)
#
# DE for tet2
ind1 <- hash.match1$exp == '127' #& hash.match1$clid %in% c(1, 2)#, 6, 7)
dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                  design = ~ genotype + clid)
dds_mut <- DESeq(dds_mut)
resultsNames(dds_mut)
res_rna_tet2 <- lfcShrink(dds_mut, coef="genotype_TET2_vs_WT", type="normal")
sum(res_rna_tet2$padj < 0.1, na.rm = T)
#
# comparison of DNMT3A and TET2 DE trends
ind1 <- T
table(ind1)
pl <- qplot(res_rna_dnmt3a$log2FoldChange[ind1], res_rna_tet2$log2FoldChange[ind1], alpha = I(0.2)) + theme_bw() +
  labs(x = 'DNMT3A DE', y = 'TET2 DE') + theme(text = element_text(size = 20))
pl
png('figures/DE_dnmt3a_tet2.png')
pl
dev.off()
# testing significance of anti-correlation
cr.test <- cor.test(res_rna_dnmt3a$log2FoldChange[ind1], res_rna_tet2$log2FoldChange[ind1],
                method = "pearson")
cr.test
# comparison with previous results
cl <- 1
gg <- intersect(rownames(res_rna_dnmt3a), rownames(de.env$de.mut2wt.shrink[[cl]][['DNMT3A']]))
str(gg)
res <- cbind(res_rna_dnmt3a[gg,]$log2FoldChange, res_rna_tet2[gg,]$log2FoldChange,
             de.env$de.mut2wt.shrink[[cl]][['DNMT3A']][gg,]$log2FoldChange, de.env$de.mut2wt.shrink[[cl]][['TET2']][gg,]$log2FoldChange)
colnames(res) <- c('DNMT3A-multi', 'TET2-multi', 'DNMT3A-rna', 'TET2-rna')
cr <- cor(res, use = 'na.or.complete')
diag(cr) <- NA
pheat <- pheatmap::pheatmap(cr, breaks = seq(-0.2, 0.2, length.out = 100))
png('figures/DE_sample_heatmap.png')
pheat
dev.off()

### run joint DA peak analysis across progenitor clusters
mat <- do.call(cbind, clust.atac.list)
clid <- unlist(lapply(1:length(clust.atac.list), function(i) rep(i, ncol(clust.atac.list[[i]])) ))
hash.match1 <- hash.match[colnames(mat), , drop = F]
hash.match1$clid <- clid
colnames(mat) <- rownames(hash.match1)
# DE for dnmt3a
ind1 <- hash.match1$exp == '123' & hash.match1$clid %in% c(1,2,6)#c(0, 1, 5)#c(1, 2)#, 6, 7)
dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                  design = ~ genotype + clid)
dds_mut <- DESeq(dds_mut)
resultsNames(dds_mut)
res_atac_dnmt3a <- lfcShrink(dds_mut, coef="genotype_DNMT3A_vs_WT", type="normal")
sum(res_atac_dnmt3a$padj < 0.1, na.rm = T)
#
# DE for tet2
ind1 <- hash.match1$exp == '127' & hash.match1$clid  %in% c(1,2,6) #c(0, 1, 5) #c(1, 2) #, 6, 7)
dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                  design = ~ genotype + clid)
dds_mut <- DESeq(dds_mut)
resultsNames(dds_mut)
res_atac_tet2 <- lfcShrink(dds_mut, coef="genotype_TET2_vs_WT", type="normal")
sum(res_atac_tet2$padj < 0.1, na.rm = T)
#
# comparison of DNMT3A and TET2 DE trends
ind1 <- T
table(ind1)
pl <- qplot(res_atac_dnmt3a$log2FoldChange[ind1], res_atac_tet2$log2FoldChange[ind1], alpha = I(0.2)) + theme_bw() +
  labs(x = 'DNMT3A DE', y = 'TET2 DE') + theme(text = element_text(size = 20))
pl
png('figures/DA_dnmt3a_tet2.png')
pl
dev.off()
# testing significance of anti-correlation
cr.test <- cor.test(res_atac_dnmt3a$log2FoldChange[ind1], res_atac_tet2$log2FoldChange[ind1],
                    method = "pearson")
cr.test

# create clusters:
# unite clusters 1,2,6 across mutants
x <- paste(cell.clone.multi$sample, cell.clone.multi$genotype, sep = '_')
x[is.na(cell.clone.multi$genotype)] <- NA
names(x) <- rownames(cell.clone.multi)
head(x)
table(x)
obj.combined$progeno <- x
Idents(obj.combined) <- obj.combined$progeno

pert <- 'TET2'
#pert <- 'DNMT3A'
if (pert == 'DNMT3A'){
  res_mut <- res_atac_dnmt3a
  fname <- 'figures/DNMT3A_coverage_plot.pdf'
}else if (pert == 'TET2'){
  res_mut <- res_atac_tet2
  fname <- 'figures/TET2_coverage_plot.pdf'

}
pdf(fname)
for (i in 1:30){
  print(i)
  x <- res_mut[order(res_mut$padj),]
  #head(x[!is.na(x$log2FoldChange) & x$log2FoldChange > 0, ], 30)

  reg <- rownames(x[!is.na(x$log2FoldChange) & x$log2FoldChange > 0, ])[i]
  #i <- 1
  #reg <- rownames(res_mut[order(res_mut$padj),])[i]
  pl <- CoveragePlot(
    object = obj.combined,
    region = reg, #'chr2-20305087-20305984',
    #features = "MS4A1",
    #expression.assay = "SCT",
    #idents = idents.plot,
    idents = unique(obj.combined$progeno)[!is.na(unique(obj.combined$progeno))],
    extend.upstream = 10000,
    extend.downstream = 10000
  )
  print(pl)

}
dev.off()



# save data tables
ll <- list(data.frame(res_atac_dnmt3a[order(res_atac_dnmt3a$pvalue),]),
           data.frame(res_atac_tet2[order(res_atac_tet2$pvalue),]))
names(ll) <- c('DNMT3A', 'TET2')
write.xlsx(ll, '../resource/DApeaks_genotypes.xlsx', row.names = TRUE, col.names = TRUE)


x <- readxl::read_xlsx('../resource/DApeaks_genotypes.xlsx', sheet = 2)
sum(x$padj < 1e-1, na.rm=T)





