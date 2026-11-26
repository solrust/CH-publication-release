library(Matrix)
#library(ggplot2)
#library(parallel)
#library(pagoda2)

###
hash.list <- lapply(1:8, function(i){
  if (i < 5){
    exper <- '123'
  }else{
    exper <- '127'
    i <- i - 4
  }
  matrix_dir = paste("data/raw/exp",exper,"/hashtags/RNA_",i,"/Results/umi_count/",sep="")

  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

  hashtags <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  #colnames(mhashtag) = paste("one_", barcode.names$V1, "-1", sep="")
  colnames(hashtags) = barcode.names$V1
  rownames(hashtags) = feature.names$V1
  return(hashtags)
})

# for exp123
hash.id <- list(c(1, 2), c(3, 3), c(1, 3), c(1, 2))
hash.id <- c(hash.id, list(c(1, 2), c(3, 1), c(3, 3), c(1, 2)))
hash.match <- as.data.frame(
  rbind( c(1, 1, "A0", '123', 'WT'),
         c(1, 2, "A1", '123', 'WT'),
         c(2, 3, "A2", '123', 'WT'),
         c(3, 1, "B0", '123', 'DNMT3A'),
         c(3, 3, "B1", '123', 'DNMT3A'),
         c(4, 1, "B2", '123', 'DNMT3A'),
         c(4, 2, "B3", '123', 'DNMT3A'),
         ###exp127:
         c(1, 1, "A3", '127', 'TET2'),
         c(1, 2, "A4", '127', 'TET2'),
         c(2, 3, "B0", '127', 'TET2'),
         c(2, 1, "B1", '127', 'TET2'),
         c(3, 3, "B2", '127', 'WT'),
         c(4, 1, "B3", '127', 'WT'),
         c(4, 2, "B4", '127', 'WT')))
colnames(hash.match) <- c("RNA", "hashtag", "mouse", 'exp', 'genotype')
hash.match$RNA <- as.numeric(hash.match$RNA)
hash.match$hashtag <- as.numeric(hash.match$hashtag)
rownames(hash.match) <- paste(hash.match$exp, hash.match$mouse, sep = '_')
hash.match$genotype <- factor(hash.match$genotype, levels = c('WT', 'DNMT3A', 'TET2'))

cutoffs <- c(150, 90)
png("figures/hashtag_decomposition.png", width = 150 * 3, height = 150 * 2)
#pdf("/home/solrust/CH/draft/figures/hashtag_decomposition.pdf", width = 2 * 3, height = 2 * 2)
cex <- 1.5 #1.5
show.legend <- FALSE
par(mfrow=c(2, 3), mar = c(5,5,1.5,1.5), mgp = c(3,1,0))
for (i in #1:length(hash.list)
     setdiff(1:length(hash.list), c(2, 7)) # do not show samples with hashtags
){
  if (i < 5){
    cutoff <- cutoffs[1]
  }else{
    cutoff <- cutoffs[2]
  }
  plot(hash.list[[i]][hash.id[[i]][1],], hash.list[[i]][hash.id[[i]][2],],
       xlim = c(1, 1e+4), ylim = c(1, 1e+4),
       pch = 19, cex = 0.1, log = 'xy', col = adjustcolor('black', 0.1),
       cex.lab = cex,
       xaxt='n', yaxt='n', #cex.axis = 4,#1.5,
       xlab = paste("hashtag",hash.id[[i]][1]), ylab = paste("hashtag",hash.id[[i]][2]) )
  axis(side=1, at=c(1, 1e+2, 1e+4), cex.axis = cex)
  axis(side=2, at=c(1, 1e+2, 1e+4), cex.axis = cex)
  abline(h=cutoff, v=cutoff, col = "red")
  title(paste("RNA_",i,sep=""))
  x12 <- sum(hash.list[[i]][hash.id[[i]][1],] > cutoff & hash.list[[i]][hash.id[[i]][2],] > cutoff)
  x1 <- sum(hash.list[[i]][hash.id[[i]][1],] > cutoff & hash.list[[i]][hash.id[[i]][2],] < cutoff)
  x2 <- sum(hash.list[[i]][hash.id[[i]][1],] < cutoff & hash.list[[i]][hash.id[[i]][2],] > cutoff)
  if (i != 2 & show.legend){
    legend("topright", as.character(c(paste("h1=",x1), paste("h2=",x2), paste("double=",x12))), bty="n")
  }
}
dev.off()


hash.annot <- lapply(1:length(hash.list), function(i){
  print(i)
  if (i < 5){
    cutoff <- cutoffs[1]
    exper <- '123'
    j <- i
  }else{
    cutoff <- cutoffs[2]
    exper <- '127'
    j <- i - 4
  }
  tannot <- rep('empty', length = ncol(hash.list[[i]]))
  names(tannot) <- colnames(hash.list[[i]])
  ind <-  hash.list[[i]][hash.id[[i]][1],] > cutoff & hash.list[[i]][hash.id[[i]][2],] > cutoff
  tannot[ind] <- 'doublets'
  ind <- hash.list[[i]][hash.id[[i]][1],] > cutoff & hash.list[[i]][hash.id[[i]][2],] < cutoff
  tannot[ind] <- hash.match$mouse[hash.match$exp == exper & hash.match$RNA == j & hash.match$hashtag == hash.id[[i]][1]]
  ind <- hash.list[[i]][hash.id[[i]][2],] > cutoff & hash.list[[i]][hash.id[[i]][1],] < cutoff
  tannot[ind] <- hash.match$mouse[hash.match$exp == exper & hash.match$RNA == j & hash.match$hashtag == hash.id[[i]][2]]
  return(tannot)
})
hash.annot[[2]][hash.annot[[2]] == "doublets"] <- 'A2'
hash.annot[[7]][hash.annot[[7]] == "doublets"] <- 'B2'
table(hash.annot[[2]])
table(hash.annot[[7]])
for (i in 1:length(hash.annot)){
  if (i < 5){
    label <-  paste('123', i, sep ='_')
  }else{
    label <-  paste('127', i-4, sep ='_')
  }
  names(hash.annot[[i]]) <- paste(label, names(hash.annot[[i]]), sep = '_')
}
saveRDS(hash.annot, file = "data/processed/hashtag_info.rds")
saveRDS(hash.match, file = "data/processed/hashtag_meta.rds")

# TODO
# save hash.match too? I need it later, right?
# use a single file to store hash info?

