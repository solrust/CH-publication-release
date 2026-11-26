library(Seurat)
library(ggplot2)
library(DESeq2)

### load human CH data
ct <- read.csv('data/count.txt', sep = '\t', row.names = 1)
dim(ct)
ct1 <- Matrix::as.matrix(ct)
#
data_genotyping <- read.csv('data/TARGET-seq_genotyping_data.tsv', sep = '\t')
metadata <- read.csv('data/TARGET-seq_metadata.tsv', sep = '\t')

### load cell cycle data
gprogs <- readRDS("../RNA_analysis/data/genesets/processed/fate_markers.rds")
genes.s <- gprogs[['S phase']]; genes.g2m <- gprogs[['G2M']]
table.cycle <- read.table("../RNA_analysis/data/genesets/download/CellCycleGeneList_1134.txt",sep="\t",header=TRUE);table.cycle = table.cycle[,c("PHASE","NAME")]
cycle_genes <- unlist(lapply(table.cycle$NAME,function(name) {strsplit(as.character(name)," ")[[1]][1]}))
cycle_genes <- unique(c(cycle_genes, genes.s, genes.g2m))
#
gs.atlas <- readxl::read_excel("../RNA_analysis/data/genesets/download/Level-3M-MarkerGenes-matrix.xlsx",
                          col_names = T, sheet = 2)
gs.atlas <- as.data.frame(gs.atlas)
rownames(gs.atlas) <- gs.atlas[,1]; gs.atlas <- gs.atlas[,-1]

### process primary human data
obj.human <- CreateSeuratObject(counts = ct1, meta.data = metadata)
obj.human[["RNA"]] <- split(obj.human[["RNA"]], f = obj.human$Sample)
obj.human <- NormalizeData(obj.human)
obj.human <- FindVariableFeatures(obj.human)
obj.human <- ScaleData(obj.human)
obj.human <- RunPCA(obj.human)
obj.human <- FindNeighbors(obj.human, dims = 1:30, reduction = "pca")
obj.human <- FindClusters(obj.human, resolution = 2, cluster.name = "clusters")
obj.human <- RunUMAP(obj.human, reduction = "pca", dims = 1:30, reduction.name = "umap.pca")
obj.human
obj.human$mutation <- data_genotyping$mutation
obj.human$genotype <- data_genotyping$genotype
DimPlot(obj.human, reduction = "umap.pca")
DimPlot(obj.human, reduction = "umap.pca",
        group.by = c('Sample', 'Genotype', 'Cluster', 'Lin_CD34_CD38'), ncol = 2)

# add Neu and Mono scores
obj.human <- AddModuleScore(
  object = obj.human,
  features = list(rownames(gs.atlas)[gs.atlas[,'preNeu'] > 0.2],
                  rownames(gs.atlas)[gs.atlas[,'Mono-1'] > 0.2]),
  name = c("preNeu", "Mono"), nbin = 24)
FeaturePlot(obj.human, reduction = "umap.pca",
            features = c('preNeu1', 'Mono2'),
            cols = c("grey95", "red"), pt.size = 0.5, order = TRUE)
#saveRDS(obj.human, file = 'data/human_ch.rds')

### analysis of cycling cells
obj.human <- CellCycleScoring(obj.human, s.features = genes.s, g2m.features = genes.g2m, set.ident = TRUE)
obj.human <- RunPCA(obj.human, features = c(genes.s, genes.g2m),
                       reduction.name = 'pca_cc')
DimPlot(obj.human, reduction = "pca_cc", group.by = c('Phase'),
        size = I(0.02), alpha = I(0.5))




### perform DE reanalysis
head(metadata[, c('Sample', 'Sample_type', 'Cluster', 'Genotype')])
ids.keep <- rownames(metadata)[metadata$Genotype %in% c('DNMT3A', 'TET2', 'WT-CH', 'WT-control')]
metadata.summary <- unique(metadata[ids.keep, c('Sample', 'Cluster', 'Genotype')])
# annotate control, tet2 or dnmt3 samples
tbl.s2g <- table(metadata.summary[,c('Sample', 'Genotype')])[,c('DNMT3A', 'TET2')]
metadata.summary$sample_state <- 'control'
metadata.summary$sample_state[which(metadata.summary$Sample %in% rownames(tbl.s2g)[tbl.s2g[,'DNMT3A'] > 0])] <- 'DNMT3A'
metadata.summary$sample_state[which(metadata.summary$Sample %in% rownames(tbl.s2g)[tbl.s2g[,'TET2'] > 0])] <- 'TET2'
metadata.summary$sample_state[which(metadata.summary$Sample %in% rownames(tbl.s2g)[ tbl.s2g[,'DNMT3A'] > 0 & tbl.s2g[,'TET2'] > 0])] <- 'DNMT3A_TET2'
head(metadata.summary)

## create cluster-specific pseudobulk count matrices
count.pseudobulk <- do.call(cbind, lapply(1:nrow(metadata.summary), function(i){
  print(i)
  cls1 <- metadata$Sample %in% metadata.summary$Sample[i] &
    metadata$Cluster %in% metadata.summary$Cluster[i] &
    metadata$Genotype %in% metadata.summary$Genotype[i] & rownames(metadata) %in% ids.keep
  cells <- metadata$Cell[cls1]
  rowSums(ct1[, gsub(pattern = "\\.", replacement = "-", x = colnames(ct1)) %in% cells, drop = F])
}))

## make simple visualization of genotype-specific GE
#head(res_tet2_hsc$res_mut[order(-res_tet2_hsc$res_mut$log2FoldChange),],20)
gn <- 'CTC-425O23.2' # 'ARHGEF26-AS1'
mut <- 'TET2'
indx <- metadata.summary$Cluster %in% subpops & #c('HSC/MPP') &
  metadata.summary$sample_state %in% mut
x <- count.pseudobulk[gn, indx] / colMeans(count.pseudobulk[, indx])
smpl <- metadata.summary[indx,]$Sample
ggplot(metadata.summary[indx, ], aes(x = Genotype, y = x, col = smpl, group = smpl)) +
  geom_point(size = 3) + geom_line() + facet_wrap(~ Cluster, scales = "free_y") + theme_bw()
if (mut == 'TET2') res_tet2$res_mut[gn,] else res_dnmt3a$res_mut[gn,]


## focus on specific cell states
subpops <- c('HSC/MPP', 'Early GMP', 'Late GMP', 'EBMP', 'EMPP', 'LMPP', 'LMPP cycling', 'MLP',
             'EryP cycling', 'EryP non-cycling', 'Erythroblasts')
subpops <- 'HSC/MPP'

## test intra-sample DNMT3A vs WT cells
inds <- metadata.summary$sample_state %in% c('DNMT3A', 'DNMT3A_TET2')
inds <- inds & metadata.summary$Genotype %in% c('WT-CH', 'DNMT3A')
inds <- inds & metadata.summary$Cluster %in% subpops
metadata.summary[inds, , drop = F]
table(inds)
dsq2_DNMT3A_intra <- DESeqDataSetFromMatrix(countData = count.pseudobulk[, inds],
                                            colData = metadata.summary[inds, , drop = F],
                                            design = ~ Genotype + Sample)
dsq2_DNMT3A_intra <- DESeq(dsq2_DNMT3A_intra)
res_DNMT3A_intra <- results(dsq2_DNMT3A_intra, contrast = c("Genotype", "DNMT3A", "WT-CH"))
res_DNMT3A_intra_shrink <- lfcShrink(dds = dsq2_DNMT3A_intra,
                                     res = res_DNMT3A_intra,
                                     contrast = c("Genotype", "DNMT3A", "WT-CH"),
                                     type = "normal")
head(res_DNMT3A_intra_shrink)

## test intra-sample TET2 vs WT cells
inds <- metadata.summary$sample_state %in% c('TET2', 'DNMT3A_TET2')
inds <- inds & metadata.summary$Genotype %in% c('WT-CH', 'TET2')
inds <- inds & metadata.summary$Cluster %in% subpops
metadata.summary[inds, , drop = F]
table(inds)
dsq2_TET2_intra <- DESeqDataSetFromMatrix(countData = count.pseudobulk[, inds],
                                            colData = metadata.summary[inds, , drop = F],
                                            design = ~ Genotype + Sample)

dsq2_TET2_intra <- DESeq(dsq2_TET2_intra)
res_TET2_intra <- results(dsq2_TET2_intra, contrast = c("Genotype", "TET2", "WT-CH"))
res_TET2_intra_shrink <- lfcShrink(dds = dsq2_TET2_intra,
                                     res = res_TET2_intra,
                                     contrast = c("Genotype", "TET2", "WT-CH"),
                                     type = "normal")
head(res_TET2_intra_shrink)


## test inter-sample DNMT3A_WT cells vs control samples
inds <- metadata.summary$sample_state %in% c('DNMT3A', 'control')
inds <- metadata.summary$sample_state %in% c('control') |
  (metadata.summary$sample_state %in% c('DNMT3A') & metadata.summary$Genotype %in% 'WT-CH')
inds <- inds & metadata.summary$Cluster %in% subpops
metadata.summary[inds, , drop = F]
table(inds)
dsq2_DNMT3A_wt_inter <- DESeqDataSetFromMatrix(countData = count.pseudobulk[, inds],
                                            colData = metadata.summary[inds, , drop = F],
                                            design = ~ sample_state)
dsq2_DNMT3A_wt_inter <- DESeq(dsq2_DNMT3A_wt_inter)
res_DNMT3A_wt_inter <- results(dsq2_DNMT3A_wt_inter, contrast = c("sample_state", "DNMT3A", "control"))
res_DNMT3A_wt_inter_shrink <- lfcShrink(dds = dsq2_DNMT3A_wt_inter,
                                     res = res_DNMT3A_wt_inter,
                                     contrast = c("sample_state", "DNMT3A", "control"),
                                     type = "normal")
head(res_DNMT3A_wt_inter_shrink)

## test inter-sample DNMT3A_Mut cells vs control samples
inds <- metadata.summary$sample_state %in% c('DNMT3A', 'control')
inds <- metadata.summary$sample_state %in% c('control') |
  (metadata.summary$sample_state %in% c('DNMT3A') & metadata.summary$Genotype %in% 'DNMT3A')
inds <- inds & metadata.summary$Cluster %in% subpops
metadata.summary[inds, , drop = F]
table(inds)
dsq2_DNMT3A_mut_inter <- DESeqDataSetFromMatrix(countData = count.pseudobulk[, inds],
                                               colData = metadata.summary[inds, , drop = F],
                                               design = ~ sample_state)
dsq2_DNMT3A_mut_inter <- DESeq(dsq2_DNMT3A_mut_inter)
res_DNMT3A_mut_inter <- results(dsq2_DNMT3A_mut_inter, contrast = c("sample_state", "DNMT3A", "control"))
res_DNMT3A_mut_inter_shrink <- lfcShrink(dds = dsq2_DNMT3A_mut_inter,
                                        res = res_DNMT3A_mut_inter,
                                        contrast = c("sample_state", "DNMT3A", "control"),
                                        type = "normal")
head(res_DNMT3A_mut_inter_shrink)

## test inter-sample TET2_WT cells vs control samples
inds <- metadata.summary$sample_state %in% c('TET2', 'control')
inds <- metadata.summary$sample_state %in% c('control') |
  (metadata.summary$sample_state %in% c('TET2') & metadata.summary$Genotype %in% 'WT-CH')
inds <- inds & metadata.summary$Cluster %in% subpops
metadata.summary[inds, , drop = F]
table(inds)
dsq2_TET2_wt_inter <- DESeqDataSetFromMatrix(countData = count.pseudobulk[, inds],
                                               colData = metadata.summary[inds, , drop = F],
                                               design = ~ sample_state)
dsq2_TET2_wt_inter <- DESeq(dsq2_TET2_wt_inter)
res_TET2_wt_inter <- results(dsq2_TET2_wt_inter, contrast = c("sample_state", "TET2", "control"))
res_TET2_wt_inter_shrink <- lfcShrink(dds = dsq2_TET2_wt_inter,
                                        res = res_TET2_wt_inter,
                                        contrast = c("sample_state", "TET2", "control"),
                                        type = "normal")
head(res_TET2_wt_inter_shrink)

## test inter-sample TET2_Mut cells vs control samples
inds <- metadata.summary$sample_state %in% c('TET2', 'control')
inds <- metadata.summary$sample_state %in% c('control') |
  (metadata.summary$sample_state %in% c('TET2') & metadata.summary$Genotype %in% 'TET2')
inds <- inds & metadata.summary$Cluster %in% subpops
metadata.summary[inds, , drop = F]
table(inds)
dsq2_TET2_mut_inter <- DESeqDataSetFromMatrix(countData = count.pseudobulk[, inds],
                                             colData = metadata.summary[inds, , drop = F],
                                             design = ~ sample_state)
dsq2_TET2_mut_inter <- DESeq(dsq2_TET2_mut_inter)
res_TET2_mut_inter <- results(dsq2_TET2_mut_inter, contrast = c("sample_state", "TET2", "control"))
res_TET2_mut_inter_shrink <- lfcShrink(dds = dsq2_TET2_mut_inter,
                                      res = res_TET2_mut_inter,
                                      contrast = c("sample_state", "TET2", "control"),
                                      type = "normal")
head(res_TET2_mut_inter_shrink)

saveRDS(file = 'data/human_DE_analysis.rds', list(DNMT3A_intra = res_DNMT3A_intra_shrink, TET2_intra = res_TET2_intra_shrink,
             DNMT3A_wt_inter = res_DNMT3A_wt_inter_shrink, TET2_wt_inter = res_TET2_wt_inter_shrink,
             DNMT3A_mut_inter = res_DNMT3A_mut_inter_shrink, TET2_mut_inter = res_TET2_mut_inter_shrink
             ))


de_human <- readRDS('data/human_DE_analysis.rds')

# compare with paper's DE results
names(de_human)
names(human_ch_del)
nm1 <- "DNMT3A_intra"; nm2 <- "DNMT3A_mut"
nm1 <- "TET2_intra"; nm2 <- "TET2_mut"
nm1 <- "DNMT3A_wt_inter"; nm2 <- "DNMT3A_env"
nm1 <- "DNMT3A_mut_inter"; nm2 <- "DNMT3A_env"

gns.shared <- intersect(rownames(de_human[[nm1]]), rownames(human_ch_del[[nm2]]))
qplot(de_human[[nm1]][gns.shared,]$log2FoldChange, human_ch_del[[nm2]][gns.shared,]$log2FC,
      main = paste(nm1, 'vs', nm2))
cor(de_human[[nm1]][gns.shared,]$log2FoldChange, human_ch_del[[nm2]][gns.shared,]$log2FC,
    use = 'na.or.complete')


