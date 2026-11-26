library(Matrix)
library(pagoda2)
library(ggplot2)
library(parallel)
library(dplyr)
library(RColorBrewer)



#
hash.annot <- readRDS(file = "data/processed/hashtag_info.rds")
hash.match <- readRDS(file = "data/processed/hashtag_meta.rds")

# read in scRNA-seq datasets
cdl <- lapply(1:8, function(i){
  if (i < 5){
    cd <- read.10x.matrices(paste('data/raw/exp123/RNA/RNA', i, '/filtered_feature_bc_matrix',sep=''))
  }else{
    cd <- read.10x.matrices(paste('data/raw/exp127/RNA/RNA', i-4 , '/filtered_feature_bc_matrix',sep=''))
  }
  return(cd)
})
unlist(lapply(cdl, function(cd) ncol(cd)))

# rename cells
for (i in 1:length(hash.annot)){
  if (i < 5){
    label <-  paste('123', i, sep ='_')
  }else{
    label <-  paste('127', i-4, sep ='_')
  }
  colnames(cdl[[i]]) <- unlist(lapply(colnames(cdl[[i]]), function(nm) unlist(strsplit(unlist(strsplit(nm, '_'))[2], '-'))[1] ))
  colnames(cdl[[i]]) <- paste(label, colnames(cdl[[i]]), sep = '_')
}

# combine samples
cd <- Reduce(cbind, cdl)
sampleid0 <- sampleid <- unlist(lapply(1:length(hash.annot), function(i){
  expid <- '123'
  if (i >= 5) expid <- '127'
  x <- paste(expid, hash.annot[[i]], sep = '_')
  names(x) <- names(hash.annot[[i]])
  x[hash.annot[[i]] == 'empty'] <- 'empty'
  x[hash.annot[[i]] == 'doublets'] <- 'doublets'
  return(x)
}))

# remove empty cells
cells <- intersect(colnames(cd), names(sampleid0))
str(cells)
cells <- cells[! sampleid0[cells] %in% c('empty')]
length(cells)
cd <- cd[, cells]
rownames(cd) <- make.unique(rownames(cd))
sampleid <- sampleid0[cells]


# exclude cell cycle-correlated genes
p2init <- Pagoda2$new(cd, log.scale = TRUE, n.cores = 20, batch = sampleid)
cycle_genes <- readRDS("data/genesets/processed/cycle_genes.rds")
cycle_genes <- intersect(cycle_genes, colnames(p2init$counts))
cycle_genes <- cycle_genes[apply(p2init$counts[,cycle_genes], 2, sd) > 0]
str(cycle_genes)
pr.cycle <- prcomp(p2init$counts[, cycle_genes], center = TRUE, scale. = TRUE)
par(mfrow=c(1,1))
plot(pr.cycle$x[,1:2], cex = 0.1, pch = 19, col = adjustcolor('grey60', 0.2))
ccor <- cor(pr.cycle$x[,1], as.matrix(p2init$counts))
ccor <- ccor[1,]
qplot(ccor)
head(sort(ccor, decreasing = T))
#
cycle.cutoff <- 0.2
sum(ccor > cycle.cutoff, na.rm = T); sum(ccor < -cycle.cutoff, na.rm = T)
#
gexcl <- unique(c(cycle_genes, names(ccor)[!is.na(ccor) & abs(ccor) >= cycle.cutoff]))
str(gexcl)
rm(p2init)

# run Pagoda2 without mito cells
p2 <- Pagoda2$new(cd[!rownames(cd) %in% gexcl, ], log.scale=TRUE, n.cores=20, batch = sampleid)
p2$adjustVariance(plot=T, gam.k=10)
p2$calculatePcaReduction(nPcs=100, n.odgenes=5e3, maxit = 300)
p2$makeKnnGraph(k=30, type='PCA', center=T, distance='cosine');
#p2$getKnnClusters(method = igraph::multilevel.community, type = "PCA", name = "multilevel")
p2$getEmbedding(type='PCA', embeddingType = 'UMAP', name = 'UMAP1', distance = 'cosine', n_neighbors = 40)
#
p2$getKnnClusters(method = igraph::cluster_leiden, type = "PCA", resolution_parameter = 4, objective_function = 'modularity', name = "leiden")
length(unique(p2$clusters$PCA$leiden))
p2$plotEmbedding(type = "PCA", embeddingType = 'UMAP1', clusterType = 'leiden')
p2$embeddings$PCA$UMAP <- p2$embeddings$PCA$UMAP1
#saveRDS(p2, file = 'data/p2objects/p2.rds')

# perform downstream cleanup
p2 <- readRDS(file = 'data/p2objects/p2.rds')

# choose mito-active cells
genes.mt <- rownames(cd)[grepl("MT-",rownames(cd))]
str(genes.mt)
pr.mito0 <- prcomp(p2$counts[,genes.mt], center = TRUE, scale. = TRUE)
mitcor <- cor(pr.mito0$x[,1], as.matrix(p2$counts))
mitcor <- mitcor[1,]
qplot(abs(mitcor))
genes.mt1 <- unique(c(names(mitcor)[!is.na(mitcor) & abs(mitcor) > 0.4], genes.mt))
str(genes.mt1)
pr.mito <- prcomp(p2$counts[,genes.mt1], center = TRUE, scale. = TRUE)
#pr.mito$x[,1] <- -pr.mito$x[,1]
plot(pr.mito$x[,1:2], pch = 19, cex = 0.05, col = adjustcolor('grey60', 0.5))
mito.clust <- unlist(lapply(1:length(unique(p2$clusters$PCA$leiden)), function(cl){
  median(pr.mito$x[,1][p2$clusters$PCA$leiden == cl])
}))
qplot(1:length(mito.clust), mito.clust)

## annotate mito cells
mito.cutoff <- 11 #11
cells.mito <- c(rownames(pr.mito$x)[pr.mito$x[,1] >= mito.cutoff],
                names(p2$clusters$PCA$leiden)[p2$clusters$PCA$leiden %in% which(mito.clust > mito.cutoff)] )
cells.mito <- unique(cells.mito)
str(cells.mito)
# set up doublet.cutoff
dsize <- table(as.factor(p2$clusters$PCA$leiden)[sampleid0[rownames(p2$embeddings$PCA$UMAP)] == 'doublets'])
clsize <- table(as.factor(p2$clusters$PCA$leiden))
doublet.cutoff <- 0.3
cells.doublet <- names(p2$clusters$PCA$leiden)[p2$clusters$PCA$leiden %in% which(dsize/clsize >= doublet.cutoff)]
str(cells.doublet)
cells.doublet <- unique(c(cells.doublet, names(sampleid)[sampleid == 'doublets'] ))
str(cells.doublet)

png("figures/cleanup_mito_doublets.png", width = 250 * 2, height = 250 * 2)
par(mfrow=c(2,2), mar = c(5,5,2,2), mgp = c(3,1,0))
cex <- 2
### plot mito
cols <- colorRampPalette(c("grey90","grey60", "red" ))(n = 10)[as.numeric(cut(pr.mito$x[,1], 10))]
plot(p2$embeddings$PCA$UMAP, col = cols, pch = 19, cex = 0.05, main = 'mito-score',
     cex.main = cex, font.main = 1, cex.lab = cex, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')#, bty = 'n')
legend('bottomleft', c('high', 'low'), pch = 19, cex = cex, col = c('red', 'grey90'), bty = 'n')
hist(pr.mito$x[,1], breaks = 100, xlab = 'mito-score', ylab = '', main = '',
     cex.lab = cex, cex.axis = 1)
abline(v = mito.cutoff, col = adjustcolor('red', 0.5), lty = 1, lwd = 3)
### plot doublets
cols <- ifelse(sampleid[rownames(p2$embeddings$PCA$UMAP)] == 'doublets', 'blue', 'grey60')
plot(p2$embeddings$PCA$UMAP, col = cols, pch = 19, cex = 0.05, main = 'hashtag doublets',
     cex.main = cex, font.main = 1, cex.lab = cex, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')#, bty = 'n')
plot(as.numeric(clsize), as.numeric(dsize) / as.numeric(clsize), log = 'x',
     ylim = c(0, 0.5), pch = 21, cex = 1.5, cex.lab = cex, bty = 'n', cex.main = cex, font.main = 1, xlab = 'cluster size', ylab = 'fraction doublets',
     col = 'black', bg = adjustcolor(ifelse(dsize/clsize > 0.3, 'blue', 'grey60'), 0.5))
abline(h = doublet.cutoff, col = adjustcolor('blue', 0.5), lty = 1, lwd = 3)
dev.off()

png("figures/remove_mito_doublets.png", width = 300, height = 300)
par(mfrow=c(1,1), mar = c(2,2,2,2))
cols <- rep('grey60', nrow(p2$embeddings$PCA$UMAP))
cols[rownames(p2$embeddings$PCA$UMAP) %in% cells.mito] <- 'red'
cols[rownames(p2$embeddings$PCA$UMAP) %in% cells.doublet] <- 'blue'
plot(p2$embeddings$PCA$UMAP, col = cols, pch = 19, cex = 0.02, #main = 'hashtag doublets',
         cex.main = cex, font.main = 1, cex.lab = cex, #xaxt = 'n', yaxt = 'n',
         xlab = '', ylab = '')#, bty = 'n')
points(p2$embeddings$PCA$UMAP[cols == 'red'], col = cols[cols == 'red'], pch = 19, cex = 0.02, #main = 'hashtag doublets',
           cex.main = cex, font.main = 1, cex.lab = cex, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')#, bty = 'n')
points(p2$embeddings$PCA$UMAP[cols == 'blue'], col = cols[cols == 'blue'], pch = 19, cex = 0.02, #main = 'hashtag doublets',
           cex.main = cex, font.main = 1, cex.lab = cex, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')#, bty = 'n')
legend('bottomleft', c('mito', 'doublet'), pch = 19, cex = cex, col = c('red', 'blue'), bty = 'n')
dev.off()


# clean up cells and re-run pagoda2
rm(p2)
cells.remove <- c(cells.mito, cells.doublet)
str(cells.remove)
cls <- setdiff(colnames(cd), cells.remove)
str(cls)
p2filter <- Pagoda2$new(cd[, cls], log.scale=TRUE, n.cores=20, batch = sampleid[cls])
p2filter$adjustVariance(plot=T, gam.k=10)
p2filter$calculatePcaReduction(nPcs=100, maxit = 300)
str(p2filter$misc$odgenes)
odgns <- setdiff(p2filter$misc$odgenes, gexcl)
str(odgns)
p2filter$calculatePcaReduction(nPcs=30, maxit = 300, n.odgenes = 2000, name = 'PCAncc',
                               odgenes = odgns)
# make PCA-based analysis
p2filter$makeKnnGraph(k=30, type='PCAncc', center=T, distance='cosine')
p2filter$getEmbedding(type='PCAncc', embeddingType = 'UMAP', name = 'UMAP',
                      #distance = 'cosine',
                      distance = 'euclidean',
                      #spread = 2, n_neighbors = 15)
                      spread = 2, n_neighbors = 30)
#
p2filter$makeGeneKnnGraph(nPcs = 30)
p2filter$getKnnClusters(method = igraph::cluster_leiden, type = "PCAncc",
                        resolution_parameter = 1, objective_function = 'modularity', name = "leiden")
p2filter$getKnnClusters(method = igraph::cluster_leiden, type = "PCAncc",
                        resolution_parameter = 0.6, objective_function = 'modularity', name = "leiden_middle") # fine: 0.6
#length(unique(p2filter$clusters$PCAncc$leiden))
p2filter$plotEmbedding(type = "PCAncc", embeddingType = 'UMAP', clusterType = 'leiden_middle')
emb <- p2filter$embeddings$PCAncc$UMAP
p2filter$plotEmbedding(type = "PCAncc", embeddingType = 'UMAP', clusterType = 'leiden')

### combine clusters into final clustering scheme:
clust <- p2filter$clusters$PCAncc$leiden
clust.final <- as.numeric(clust); names(clust.final) <- names(clust)
clust.final[clust %in% c(12, 17)] <- NA # remove pre-B cells
clust.final[clust %in% c(8,9)] <- 8 # aggregate CMP
#clust.final[clust %in% c(5,6,7, 2)] <- 5 # aggregate lymphoid prog
clust.final[clust %in% c(5,6,7)] <- 5 # aggregate lymphoid prog
#clust.final[clust %in% c(10, 11)] <- 11 # aggregate myeloid
clust.final[clust %in% c(4,13,15)] <- 4 # aggregate dendritic
#clust.final[clust %in% c(1, 2)] <- 1 # aggregate MLP
table(clust.final)
# re-number clusters
clust.final <- do.call(
  recode,
  c(list(clust.final), setNames(rank(unique(clust.final)), unique(clust.final)))
)
names(clust.final) <- names(clust)
clust.final <- clust.final[!is.na(clust.final)]
#
clust.final1 <- clust.final
clust.final1[clust.final == 3] <- 1
clust.final1[clust.final == 1] <- 2
clust.final1[clust.final == 2] <- 3
clust.final1[clust.final == 5] <- 4
clust.final1[clust.final == 4] <- 5
clust.final <- clust.final1
#
p2filter$clusters$PCAncc$final <- clust.final
p2filter$getDifferentialGenes(type='PCAncc', verbose=T, clusterType='final')
#
#saveRDS(p2filter, file = 'data/p2objects/p2filter.rds')


# fix colors of clusters
clust.final.cols <- brewer.pal(length(unique(clust.final)), "Paired")
names(clust.final.cols) <- 1:length(unique(clust.final))
clust.final.names <- c('HSC', 'MLP', NA, 'L', 'D', 'CMP', 'M', 'M', 'Eo/Bas', 'Mg/E')
# record italic gene names
genes.ital <- unlist(lapply(rownames(cd), function(gn){
  gn1 <- paste(substr(gn, 1, 1), tolower(substr(gn, 2, nchar(gn))), sep = '')
  bquote(italic(.(gn1)))
}))
names(genes.ital) <- rownames(cd)
# annotate genotype colors
col.genotype <- adjustcolor(c('black', 'red', 'green'), 0.5)
names(col.genotype) <- c('WT', 'DNMT3A', 'TET2')
# fix embedding
emb <- p2filter$embeddings$PCAncc$UMAP
p2counts <- p2filter$counts
p2pca <- p2filter$reductions$PCAncc
diffGenes <- p2filter$diffgenes
p2clusts <- p2filter$clusters
#
save(cd, p2counts, p2pca, diffGenes, p2clusts, sampleid, emb, clust.final, col.genotype, clust.final.cols, genes.ital,
     file = 'data/processed/scrna.RData')



