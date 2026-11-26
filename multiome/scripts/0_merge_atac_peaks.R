library(Seurat)
library(Signac)
library(Azimuth)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(glmGamPoi)
library(parallel)
#
library(GenomicRanges)
#library(future)

# the 10x hdf5 file contains both data types.
sample.names <- c('C1', 'C2', 'D', 'T')
inputdatasets <- lapply(sample.names, function(sample.name){
  Read10X_h5(paste('../../larry_atac/multiome_larry_', sample.name, '/outs/filtered_feature_bc_matrix.h5', sep = ''))
})


### merge peaks across samples
# read ATAC peaks
peaksl <- lapply(sample.names, function(sample.name){
  peaks <- read.table(
    file <- paste('../../larry_atac/multiome_larry_', sample.name, '/outs/atac_peaks.bed', sep = ''),
    col.names = c("chr", "start", "end")
  )
})
names(peaksl) <- sample.names

# convert to genomic ranges
gr <- lapply(peaksl, function(peaks) makeGRangesFromDataFrame(peaks))
combined.peaks <- GenomicRanges::reduce(c(gr[[1]], gr[[2]], gr[[3]], gr[[4]]))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
mdl <- lapply(sample.names, function(sample.name){
  read.table(
    file = paste('../../larry_atac/multiome_larry_', sample.name, '/outs/per_barcode_metrics.csv', sep = ''),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )#[-1, ] # remove the first row
})

# perform an initial filtering of low count cells
mdl1 <- lapply(mdl, function(md){
  ind <- md$atac_fragments > 100 | md$gex_mapped_reads > 500
  md[ind,]
})
names(mdl1) <- sample.names

# create fragment objects
for (sample.name in sample.names){
  path <- paste('../../larry_atac/multiome_larry_', sample.name, '/outs/atac_fragments.tsv.gz', sep = '')
  file.copy(from = path, to = paste("data/raw/atac_fragments_", sample.name,".tsv.gz", sep = ''))
}

fragsl <- lapply(sample.names, function(sample.name){
  x <- CreateFragmentObject(
    path = paste('../../larry_atac/multiome_larry_', sample.name, '/outs/atac_fragments.tsv.gz', sep = ''),
    #path = paste("data/raw/atac_fragments_", sample.name,".tsv.gz", sep = ''),
    cells = rownames(mdl1[[sample.name]])
  )
})
names(fragsl) <- sample.names

# quantify peaks in each dataset
atac.countsl <- mclapply(sample.names, function(sample.name){
  FeatureMatrix(
    fragments = fragsl[[sample.name]],
    features = combined.peaks,
    cells = rownames(mdl1[[sample.name]])
  )
}, mc.cores = length(sample.names))
names(atac.countsl) <- sample.names


# upload expression count matrices
inputdatal <- lapply(sample.names, function(sample.name){
  inputdata.10x <- Read10X_h5(paste('../../larry_atac/multiome_larry_', sample.name, '/outs/filtered_feature_bc_matrix.h5', sep = ''))
})
names(inputdatal) <- sample.names

for (i in 1:4){
  v1 <- colnames(inputdatal[[i]]$Peaks)
  v2 <- colnames(atac.countsl[[i]])
  cls <- intersect(v1, v2)
  inputdatal[[i]]$`Gene Expression` <- inputdatal[[i]]$`Gene Expression`[,cls]
  inputdatal[[i]]$Peaks <- atac.countsl[[i]][, cls]
}

save(atac.countsl, inputdatal, fragsl, file = 'data/processed/multiome_input.RData')

