library(matrixStats)
library(ArchR)
addArchRThreads(threads = 1)

# Create Arrow Files (Point to your original .fragments.tsv.gz files)
my_input_files <- c(
  "C1" = "data/raw/atac_fragments_C1.tsv.gz",
  "C2" = "data/raw/atac_fragments_C2.tsv.gz",
  "T"  = "data/raw/atac_fragments_T.tsv.gz",
  "D"  = "data/raw/atac_fragments_D.tsv.gz"
)
ArrowFiles <- createArrowFiles(
  inputFiles = my_input_files,
  sampleNames = names(my_input_files),
  minTSS = 4,
  minFrags = 1000,
  force = TRUE
)

# Extract metadata from seurat object for ArchR
metadata <- obj.combined@meta.data
metadata$cell_id <- rownames(metadata)
metadata$cell_id_archr <- gsub("_", "#", rownames(metadata))
# match Archr and seurat cells
cell_indices <- match(proj$cellNames, metadata$cell_id_archr)
proj <- addCellColData(
  ArchRProj = proj,
  data = paste0(metadata$seurat_clusters)[cell_indices],
  name = "seurat_clusters",
  cells = proj$cellNames,
  force = TRUE
)
head(proj$seurat_clusters)
# add annotation of original clusters
proj <- addCellColData(
  ArchRProj = proj,
  data = paste0(metadata$projection)[cell_indices],
  name = "main_clusters",
  cells = proj$cellNames,
  force = TRUE
)
head(proj$main_clusters)
head(getCellColData(proj))
# QC check
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p
p1 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges",
  baseSize = 10
)
## 1
p1
## Picking joint bandwidth of 0.918
p1 <- plotFragmentSizes(ArchRProj = proj)
p1
p2 <- plotTSSEnrichment(ArchRProj = proj)
p2

# dim reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  force = TRUE
)
# clustering
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = T
)
table(proj$Clusters)
# UMAP
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                    name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                    name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                    name = "seurat_clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, p3, type = "h")

markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA2", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A"#TCells
)
p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
p
p$GATA2

# add peaks
peaks.gr <- granges(obj.combined[["ATAC"]])
head(peaks.gr)
peaks.gr1 <- peaks.gr[!grepl('chrM', peaks.gr@seqnames),]
proj <- addPeakSet(proj, peakSet = peaks.gr1, force = T)
proj <- addPeakMatrix(proj, force = T)

# identify cluster marker peaks
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerPeaks
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
markerList$C1

# Add Motif Annotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2020", force = T)
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# motif enrichment analysis
proj <- addBgdPeaks(proj, method = 'ArchR')
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotVarDev
#
motifs <- c("GATA2", "CEBPA", "EBF1", "IRF4", "MYC")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% c("SREBF1_22", "MYCN_302")]
markerMotifs
#
p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "MotifMatrix",
  name = sort(markerMotifs),
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)
p2 <- lapply(p, function(x){
  x + guides(color = "none", fill = "none") +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

# footprinting
motifPositions <- getPositions(proj)
motifPositions
#
motifs <- c("GATA2", "MYC")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% c("SREBF1_22", "MYCN_302")]
markerMotifs
# background model
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters", force = TRUE)
# get footprints
seFoot <- getFootprints(
  ArchRProj = proj,
  positions = motifPositions[markerMotifs],
  groupBy = "Clusters"
)
# visualize
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj,
  normMethod = "Subtract", #"Subtract", # "Divide" # "None"
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)


############# tailored footprinting analysis ######################
head(cell.clone.multi)
cell_id_archr <- gsub("_", "#", rownames(cell.clone.multi))
cell_indices <- match(proj$cellNames, cell_id_archr)
# add genotype information
proj <- addCellColData(
  ArchRProj = proj,
  data = paste0(cell.clone.multi$genotype)[cell_indices],
  name = "genotype",
  cells = proj$cellNames,
  force = TRUE
)
proj$genotype[proj$genotype == 'NA'] <- NA
# annotate HSPC genotypes
proj$genotype_hspc <- NA
ind.hspc <- !is.na(proj$main_clusters) & !is.na(proj$genotype) & proj$main_clusters %in% c(1, 6, 2)
proj$genotype_hspc[ind.hspc] <- proj$genotype[ind.hspc]
table(proj$genotype_hspc)

# add clone growth phenotype
breaks <- quantile(cell.clone.multi$'9', probs = seq(0, 1, 0.25), na.rm = TRUE)
rate_quartile <- cut(cell.clone.multi$'9',  breaks = breaks, include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4"))
proj <- addCellColData(
  ArchRProj = proj,
  data = paste0(rate_quartile)[cell_indices],
  name = "growth",
  cells = proj$cellNames,
  force = TRUE
)
proj$growth[proj$growth == 'NA'] <- NA
# annotate HSPC genotypes
proj$growth_hspc <- NA
ind.hspc <- !is.na(proj$main_clusters) & !is.na(proj$genotype) & proj$main_clusters %in% c(1, 6, 2)
proj$growth_hspc[ind.hspc] <- proj$growth[ind.hspc]
table(proj$growth_hspc)

#
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "growth_hspc", force = TRUE)
# select markers
motifs <- c("IRF4", "GATA2", 'ZBTB18', 'MYC')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% c("SREBF1_22", "MYCN_302")]
markerMotifs
#
seFoot_hspc <- getFootprints(
  ArchRProj = proj,
  positions = motifPositions[markerMotifs],
  groupBy = "growth_hspc"
)
#
plotFootprints(
  seFoot = seFoot_hspc,
  ArchRProj = proj,
  normMethod = "Subtract", #"Subtract", # "Divide" # "None"
  plotName = "HSPCs_ery_Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

### look into MYC Chip-seq peaks
seMYC_de <- getFootprints(
  ArchRProj = proj,
  positions = GRangesList(MYCpeak = myc_gr_centered),
  groupBy = "ery_hspc"
)
#
plotFootprints(
  seFoot = seMYC_de,
  ArchRProj = proj,
  normMethod = "Subtract", #"Subtract", # "Divide" # "None"
  plotName = "HSPCs_growth_chip_myc_de_Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

