library(GenomicRanges)

# upload gene sets
gene_sets <- readRDS('../RNA_analysis/data/genesets/processed/gene_sets.rds')
go_myc <- gene_sets[grepl('MYC', names(gene_sets))][[1]]

# 2. Load the public MYC ChIP-seq peaks
myc_chip_gr <- import("data/processed/ENCFF462OYX.bigBed")
peaks.gr <- granges(obj.combined[["ATAC"]])

peaks.gr <- granges(obj.combined[["ATAC"]])
names(peaks.gr) <- rownames(obj.combined[["ATAC"]])

myc_bound_peaks <- subsetByOverlaps(x = peaks.gr, ranges = myc_chip_gr)
peakpadj <- apply(do.call(cbind, lapply(da_commitment, function(x) x[names(myc_bound_peaks),]$padj)), 1, min)
myc_bound_peaks1 <- myc_bound_peaks#[!is.na(peakpadj) & peakpadj < 0.2, ]

# identify genes that are close to MYC-bound peaks
peak_gene_mapping <- ClosestFeature(
  object = obj.combined,
  regions = myc_bound_peaks1
)
peak_gene_mapping$peak_id <- names(myc_bound_peaks1)
head(peak_gene_mapping)
# identify confident MYC target gene list
custom_myc_targets <- peak_gene_mapping %>%
  filter(abs(distance) < 0.5e+3) %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(gene_name) %>%
  unique()
length(custom_myc_targets)

# clone phenotypes
i <- 1
tbl <- clone.ass$test.basic.list[[i]]
gscores <- tbl$pheno
# mutant DE
tbl <- de.env$de.mut2wt.shrink[[1]]$TET2
gscores <- tbl$log2FoldChange
#
names(gscores) <- rownames(tbl)
str(gscores)
#
ind1 <- names(gscores) %in% custom_myc_targets
ind1 <- names(gscores) %in% go_myc
ind1 <- names(gscores) %in% intersect(go_myc, custom_myc_targets)
#ind1 <- names(gscores) %in% setdiff(go_myc, custom_myc_targets)
qplot(gscores, col = ind1, geom = 'density')
#quantile(gscores, na.rm=T); quantile(gscores[ind1], na.rm=T)
mean(gscores[!ind1] < 0, na.rm=T)
mean(gscores[ind1] < 0, na.rm=T)
wilcox.test(gscores[ind1], gscores[!ind1])
wilcox.test(gscores[names(gscores) %in% setdiff(go_myc, custom_myc_targets)], gscores[names(gscores) %in% intersect(go_myc, custom_myc_targets)])


### run MYC footprinting analysis
# create three tiers of MYC binding peaks
myc_peak_tier1 <- names(myc_bound_peaks)
str(myc_peak_tier1)
# MYC motifs that are close to DE genes
tbl <- de.env$de.mut2wt.shrink[[1]]$TET2
genes_de <- rownames(tbl)[!is.na(tbl$padj) & tbl$padj < 0.2]
str(genes_de)
ind <- peak_gene_mapping$gene_name %in% genes_de & peak_gene_mapping$distance < 2e+3
table(ind)
myc_peak_tier2 <- peak_gene_mapping$query_region[ind]
str(myc_peak_tier2)
# close to downregulated MYC targets
ind1 <- ind & peak_gene_mapping$gene_name %in% go_myc
table(ind1)
myc_peak_tier3 <- peak_gene_mapping$query_region[ind1]
str(myc_peak_tier3)

library(Signac)
library(motifmatchr)
library(JASPAR2020) # Or the version you are using
library(TFBSTools)

str(which(unlist(motif.names) %in% 'GATA2'))
unlist(motif.names)[which(unlist(motif.names) %in% 'GATA2')]

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
# add motif information
obj.combined <- AddMotifs(obj.combined, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
obj.combined <- Footprint(
  object = obj.combined,
  motif.name = 'MYC', #c("MA0036.3"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
)
pl <- PlotFootprint(obj.combined, features = c("MYC"))
pl + patchwork::plot_layout(ncol = 1)

# 1. Convert your character vector to a GRanges object
myc_peak_tier3_gr <- StringToGRanges(myc_peak_tier1)
#
#rand.peaks <- sample(names(peaks.gr), 30000)
#rand.peaks <- names(peaks.gr)
#rand.peaks <- names(which(motif.matrix[, names(unlist(motif.names))[which(unlist(motif.names) %in% 'GATA2')]]))
#myc_peak_tier3_gr <- StringToGRanges(rand.peaks)
#
tf <- 'GATA2'
tf.name <- names(unlist(motif.names))[which(unlist(motif.names) %in% tf)]
#pfm <- getMatrixSet(x = JASPAR2020, opts = list(name = tf, species = 9606, all_versions = FALSE))
# new way
pfm <- obj.combined@assays$ATAC@motifs@pwm[tf.name]
motif.pos <- obj.combined@assays$ATAC@motifs@positions[[tf.name]]
#
myc_positions_restricted <- matchMotifs(
  pwms = pfm,
  subject = myc_peak_tier3_gr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  out = "positions"
)
#
myc_gr <- myc_positions_restricted[[1]]
#
motif.pos <- obj.combined@assays$ATAC@motifs@positions[[tf.name]]
# subset of MYC chip-seq peaks
myc_peak_tier_gr <- StringToGRanges(myc_peak_tier1)
motif.pos.restricted <- subsetByOverlaps(
  x = motif.pos, ranges = myc_peak_tier_gr,
  type = "within" # Ensures the whole motif is inside the peak
)
#
myc_gr_centered <- resize(motif.pos.restricted, #myc_gr,
                          width = 400, fix = "center")
#myc_gr_centered <- resize(myc_gr, width = 400, fix = "center")
#
#DefaultAssay(obj.combined) <- "ATAC"
#Key(obj.combined[["ATAC"]]) <- "atac_"
#
obj.combined <- Footprint(
  object = obj.combined,
  regions = myc_gr_centered,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "ATAC",
  label = "MYC_Restricted_Aligned",
  key = 'footprint_GATA2',
  #group_by = "dataset"
)
pl <- PlotFootprint(obj.combined, features = c("footprint_GATA2"),
                    group.by = 'projection')
pl + patchwork::plot_layout(ncol = 1) & ylim(-0.5, 2)


obj.combined <- Footprint(
  object = obj.combined,
  regions = myc_gr_centered,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  key = 'fooptrint_MYC',
  assay = "ATAC",
  label = "MYC_Restricted"
)

myc_footprint_data <- Footprint(
  object = obj.combined,
  regions = myc_gr_centered,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  label = "MYC_Restricted",
  key = "footprint_"  # Match your successful run
)
enrichment_data <- GetFootprintData(myc_footprint_data, features = "MYC_Restricted")



###################################### previous code: #################
myc.binding.counts <- t(motif.matrix[names(myc_bound_peaks), ]) %*% obj.combined@assays$ATAC$counts[names(myc_bound_peaks), ]

myc.ct <- colSums(motif.matrix[names(myc_bound_peaks), ])
back.ct <- colSums(motif.matrix[, ])

obj.combined@assays$ATAC$counts[names(myc_bound_peaks), ]

cl <- 10
ind.myc <- rownames(da_commitment[[cl]]) %in% names(myc_bound_peaks)
myc.nm <- names(unlist(motif.names))[which(unlist(motif.names) %in% 'MYC')]
ind.myc.atac <- rownames(da_commitment[[cl]]) %in% rownames(motif.matrix)[which(motif.matrix[,myc.nm])]

ggplot() + geom_point(aes(x = da_commitment[[cl]]$log2FoldChange, y = da_commitment[[cl]]$pvalue,
                          col = ind.myc)) +
  geom_point(aes(x = da_commitment[[cl]]$log2FoldChange, y = da_commitment[[cl]]$pvalue,
                   col = ind.myc.atac)) + scale_y_log10()

cl <- 1
ind1 <- da_commitment[[cl]]$padj < 0.2
ggplot() + #geom_density(aes(x = da_commitment[[cl]]$log2FoldChange, col = 'background')) +
  geom_density(aes(x = da_commitment[[cl]]$log2FoldChange[ind1], col = 'sign')) +
  #geom_density(aes(x = da_commitment[[cl]]$log2FoldChange[ind.myc], col = 'Myc chip')) +
  geom_density(aes(x = da_commitment[[cl]]$log2FoldChange[ind.myc&ind1], col = 'Myc chip sign')) +
  #geom_density(aes(x = da_commitment[[cl]]$log2FoldChange[ind.myc.atac], col = 'MYC atac')) +
  geom_density(aes(x = da_commitment[[cl]]$log2FoldChange[ind.myc.atac&ind1], col = 'MYC motif sign')) +
  geom_density(aes(x = da_commitment[[cl]]$log2FoldChange[ind.myc&ind.myc.atac&ind1], col = 'MYC chip+motif sign'))
names(da_commitment)[cl]

## mutant DA
# save data tables
ll <- list(data.frame(res_atac_dnmt3a[order(res_atac_dnmt3a$pvalue),]),
           data.frame(res_atac_tet2[order(res_atac_tet2$pvalue),]))
names(ll) <- c('DNMT3A', 'TET2')


da_mutants <- readxl::read_excel('../resource/DApeaks_genotypes.xlsx',sheet = 2)
head(da_mutants)

# Accessibility of MYC binding regions:
# lymp: down, myeloid: up, dentritic: down?, ery/gran/mg: up, mlp (clust 3): no change
cl <- 7
tbl.deseq2 <- da_mutants#
#tbl.deseq2 <- da_commitment[[cl]]
#ind1 <- tbl.deseq2$padj < 0.2
ind1 <- tbl.deseq2$pvalue < 0.05
table(ind1)
ggplot() + #geom_density(aes(x = tbl.deseq2$log2FoldChange, col = 'background')) +
  geom_density(aes(x = tbl.deseq2$log2FoldChange[ind1], col = 'sign')) +
  #geom_density(aes(x = tbl.deseq2$log2FoldChange[ind.myc], col = 'Myc chip')) +
  geom_density(aes(x = tbl.deseq2$log2FoldChange[ind.myc&ind1], col = 'Myc chip sign')) +
  #geom_density(aes(x = tbl.deseq2$log2FoldChange[ind.myc.atac], col = 'MYC atac')) +
  geom_density(aes(x = tbl.deseq2$log2FoldChange[ind.myc.atac&ind1], col = 'MYC motif sign')) +
  geom_density(aes(x = tbl.deseq2$log2FoldChange[ind.myc&ind.myc.atac&ind1], col = 'MYC chip+motif sign'))
x1 <- table(tbl.deseq2$log2FoldChange[ind1] > 0)
x1[2]/x1[1]
x2 <- table(tbl.deseq2$log2FoldChange[ind.myc&ind1] > 0)
x2[2]/x2[1]
x3 <- table(tbl.deseq2$log2FoldChange[ind.myc.atac&ind1] > 0)
x3[2]/x3[1]
x4 <- table(tbl.deseq2$log2FoldChange[ind.myc&ind.myc.atac&ind1] > 0)
x4[2]/x4[1]
names(da_commitment)[cl]


motif.matrix.ext <- cbind(motif.matrix,
                          MYC = rownames(motif.matrix) %in% names(myc_bound_peaks)
)
#
obj.combined <- RunChromVAR(
  object = obj.combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix = motif.matrix.ext[, c('MYC', colnames(motif.matrix.ext)[1:3])],
  new.assay.name = "chromvar_myc"
)

###### analysis of motif activities across cells
motif.object <- obj.combined@assays$ATAC@motifs
chvar <- obj.combined@assays$chromvar_myc@data
chvar[is.na(chvar)] <- 0
rownames(chvar) <- as.character(unlist(motif.object@motif.names)[rownames(chvar)])
#
gg <- intersect(rownames(chvar), rownames(obj.combined@assays$SCT$scale.data))
str(gg)
head(sort(apply(chvar[gg,,drop=F], 1, sd), decreasing = T), 20)
gn <- 'MYC'#'GATA2' #'GATA2' #'CEBPD'
#gn <- names(sort(diag(xx), decreasing = T))[17]
#gn <- names(sort(apply(chvar[gg,], 1, sd), decreasing = T))[8]
#pdf('figures/motif_expression_patterns.pdf', width = 8, height = 4)
for (gn in names(sort(apply(chvar[gg,], 1, sd), decreasing = T))){
  print(gn)
  gn <- 'MYC'
  motif.name <- ConvertMotifID(obj.combined, name = gn)
  gene_plot <- FeaturePlot(obj.combined, features = paste("sct", gn, sep = '_'),
                           reduction = 'wnn.umap')#umap.naive.rna 'umap.harmony')
  motif_plot <- FeaturePlot(obj.combined,
                            #features = 'MYC',
                            features = 'chromvarmyc__MYC',
                            min.cutoff = 0, cols = c("lightgrey", "darkred"),
                            reduction = 'wnn.umap')
  print(gene_plot | motif_plot)
}


# 1. Pull the MYC ChromVAR scores manually
# Assumes 'MYC' is the row name in that specific assay
myc_scores <- obj.combined@assays$chromvar_myc@data["MYC", ]
# 2. Add it to the metadata as a new column
obj.combined$MYC_motif_score <- myc_scores
# 3. Plot the metadata column directly
FeaturePlot(
  object = obj.combined,
  features = 'MYC_motif_score',
  reduction = 'wnn.umap',
  cols = c("lightgrey", "darkred")
)

# Plot the metadata column with the centered gradient
FeaturePlot(
  object = obj.combined,
  features = 'MYC_motif_score',
  reduction = 'wnn.umap'
) + scale_colour_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0
)

clone.ass <- new.env()
load('../Association_analysis/data/processed/use/association_summary.RData', envir = clone.ass)

head(clone.ass$gsea.clone$gsea.clone[[1]])
i <- 2
xx <- do.call(rbind, lapply(clone.ass$gsea.clone$gsea.clone, function(x)
 x[grepl('HALLMARK_MYC_TARGETS_V1', x$pathway),]
))
xx <- data.frame(xx)
rownames(xx) <- colnames(clone.ass$gsea.clone$gscore.clone)


#

