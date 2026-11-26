hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
load('../RNA_clone_analysis/data/processed/clone_multiome_info.RData')
de.env <- new.env()
load(file = '../RNA_analysis/data/processed/mutant_de.RData', envir = de.env)
env.dna <- new.env()
load('../DNA_analysis/data/processed/DNA_clones.RData', env = env.dna)
obj.combined <- readRDS(file = 'data/processed/multiome_seurat.rds')

x <- read.table( file = "data/processed/metainfo.txt")
load('../RNA_clone_analysis/data/processed/clone_info.RData')
titl <- c('myeloid biased', 'lymphoid biased', 'HSC biased')
indxs <- c(16, 213, 396)
cols <- c('#0000FF', '#8000FF', '#FFA500')
emb <- obj.combined@reductions$wnn.umap@cell.embeddings
png("figures/clones_multiome_examples.png", width = 750, height = 300)
par(mfrow = c(1,3), mar = c(1,1,4,1))
for (i in 1:length(indxs)){
  cln <- rownames(clone.info)[indxs[i]]
  y <- x[x$larry == cln,]
  cls <- rownames(y)[!is.na(y$larry)]
  sz <- length(cls)
  plot(emb, pch = 19, cex = 0.1, col = 'grey', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
       main = titl[i], font.main = 1, cex.main = 3)
  points(emb[cls,], xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
         #   pch = 21, cex = 1, lwd = 1.5, col = cols[i], bg = 'white')
         pch = 21, cex = 1, lwd = 1.5, col = cols[i], bg = 'white')

  legend('topleft', paste(#clone.info$genotype[cln],
    'clone', substr(cln, 1, 5),
    #clone.info$clone[cln],
    #'\n',
    '\nsize =', sz#clone.info$size[cln]
    #'id:', clone.info$id[cln],
  ), bty = 'n', cex = 2)
}
dev.off()
