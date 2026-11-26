library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(parallel)
library(glmGamPoi)
library(GenomicRanges)
#
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)


old.scrna <- new.env()
load(file = '../RNA_analysis/data/processed/scrna.RData', envir = old.scrna)
#
load(file = 'data/processed/multiome_input.RData')
sample.names <- names(inputdatal)
annotation.clusters <- readRDS(file = "data/processed/clusters_annotation.rds")

### joint integration analysis
seuratl <- mclapply(sample.names, function(sample.name){
  rna_counts <- inputdatal[[sample.name]]$`Gene Expression`
  atac_counts <- inputdatal[[sample.name]]$Peaks
  atac_counts <- atac_counts[grepl('chr', rownames(atac_counts)), ]
  hemo <- CreateSeuratObject(counts = rna_counts)


  chrom_assay <- CreateChromatinAssay(counts = atac_counts,
                                      sep = c(":", "-"), genome = 'hg38',
                                      fragments = fragsl[[sample.name]],
                                      min.cells = 10)
  hemo[["ATAC"]] <- chrom_assay

  DefaultAssay(hemo) <- "RNA"
  hemo <- SCTransform(hemo, verbose = FALSE) %>% RunPCA() %>%
    RunUMAP(dims = 1:50, reduction.name = 'umap.naive.rna', reduction.key = 'naive_rnaUMAP_')

  DefaultAssay(hemo) <- "ATAC"
  hemo <- RunTFIDF(hemo)
  hemo <- FindTopFeatures(hemo, min.cutoff = 'q0')
  hemo <- RunSVD(hemo)
  hemo <- RunUMAP(hemo, reduction = 'lsi', dims = 2:50,
                  reduction.name = "umap.naive.atac", reduction.key = "naiveatacUMAP_")

  # calculate a WNN graph
  hemo <- FindMultiModalNeighbors(hemo, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  hemo <- RunUMAP(hemo, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  hemo <- FindClusters(hemo, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

  hemo$dataset <- sample.name
  hemo
}, mc.cores = length(sample.names))
names(seuratl) <- sample.names

### combine datasets
obj.combined <- merge(
  x = seuratl[[1]],
  y = seuratl[-1],
  add.cell.ids = sample.names
)
obj.combined
table(obj.combined$dataset)
obj.combined$projection <- annotation.clusters$projection.clusters
### add ATAC meta information
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
seqlevelsStyle(annotation) <- 'UCSC' # "NCBI"
#
DefaultAssay(obj.combined) <- "ATAC"
Annotation(obj.combined) <- annotation
obj.combined <- NucleosomeSignal(obj.combined)
obj.combined <- TSSEnrichment(obj.combined)
# estimage GC content
obj.combined <- RegionStats(obj.combined, genome = BSgenome.Hsapiens.UCSC.hg38)

### filter out 'poor' cells
obj.combined <- subset(
  x = obj.combined,
  subset = nCount_ATAC > 500 & nCount_RNA > 500 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
obj.combined

### RNA pre-processing
DefaultAssay(obj.combined) <- "RNA"
obj.combined <- SCTransform(obj.combined, verbose = FALSE, vst.flavor = "v2") %>%
  RunPCA() %>% RunUMAP(dims = 1:50,
                       reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
### ATAC pre-processing
DefaultAssay(obj.combined) <- "ATAC"
obj.combined <- RunTFIDF(obj.combined)
obj.combined <- FindTopFeatures(obj.combined, min.cutoff = 'q0')
obj.combined <- RunSVD(obj.combined)
obj.combined <- RunUMAP(obj.combined, reduction = 'lsi', dims = 2:50,
                        reduction.name = "umap.atac", reduction.key = "atacUMAP_")
#

### Estimate gene activities using accessibility of promoters
gene.activities <- GeneActivity(obj.combined)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
obj.combined[['Promoter']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(obj.combined) <- 'Promoter'
obj.combined <- SCTransform(obj.combined,
                            assay = 'Promoter', new.assay.name = 'SCT_prom',
                            verbose = FALSE, vst.flavor = "v2")
obj.combined <- RunPCA(obj.combined ,reduction.name = 'pca.prom', assay = "SCT_prom")
obj.combined <- RunUMAP(obj.combined, dims = 1:50, reduction = 'pca.prom',
                       reduction.name = 'umap.promoter', reduction.key = 'promUMAP_')

### calculate a joint ATAC-RNA WNN graph
obj.combined <- FindMultiModalNeighbors(obj.combined, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
obj.combined <- RunUMAP(obj.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
obj.combined <- FindClusters(obj.combined, graph.name = "wsnn",
                             algorithm = 3, resolution = 1, verbose = FALSE)


### motif accessibility analysis
# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(obj.combined) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(obj.combined), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.matrix.score <- CreateMotifMatrix(features = granges(obj.combined), pwm = pwm_set, genome = 'hg38', score = TRUE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
obj.combined <- SetAssayData(obj.combined, assay = 'ATAC', slot = 'motifs',
                             new.data = motif.object)

# Note that this step can take 30-60 minutes
obj.combined <- RunChromVAR(
  object = obj.combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)
#
chvar <- obj.combined@assays$chromvar@data
chvar[is.na(chvar)] <- 0
rownames(chvar) <- as.character(unlist(motif.object@motif.names)[rownames(chvar)])
qplot(rowMeans(chvar), apply(chvar, 1, sd))


saveRDS(obj.combined, file = 'data/processed/multiome_seurat.rds')

