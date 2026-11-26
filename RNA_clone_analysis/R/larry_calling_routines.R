identify_barcodes <- function(fpath = NULL, trios = NULL,
                              filter.barcode.cov = 5, cutoff.hclust = 3, cutoff.multi = 0.8, do.plot = TRUE){
  # parameters:
  #do.plot <- TRUE
  #filter.barcode.cov <- 5
  #cutoff.hclust <- 3
  #cutoff.multi <- 0.8

  # larry object stores intermediate data on barcode calling
  larry <- list()

  # read in RNA barcodes
  if (is.null(trios)){
    larry$trios <- upload_trios(fpath)
  }else{
    larry$trios <- trios
  }

  # filter out barodes with less overall small number of reads in the dataset
  #length(larry$trios$clone_id)
  #sum(table(larry$trios$clone_id) < filter.barcode.cov)
  nms <- names(table(larry$trios$clone_id))[table(larry$trios$clone_id) >= filter.barcode.cov]
  larry$trios <- larry$trios[larry$trios$clone_id %in% nms, ]
  dim(larry$trios)

  # create a cell-barcode matrix
  larry$mat <- create_clone_matrix(larry)

  # estimate hamming distance between barcode sequences
  larry$bardist <- hamming_barcodes(larry)
  dim(larry$bardist)

  # cluster barcodes with similar sequences
  hclust.barcode <- hclust(as.dist(larry$bardist), method = "single")
  if (do.plot){
    par(mfrow=c(1,1), mar = c(4,4,1,1))
    plot(hclust.barcode, labels =  rep("", length(hclust.barcode$order)), cex = 0.1)
    abline(h=cutoff.hclust, col = "red")
  }
  clust <- cutree(hclust.barcode, h = cutoff.hclust)
  length(unique(clust))
  larry$hclust.barcode <- hclust.barcode
  larry$clust <- clust
  larry$cutoff.hclust <- cutoff.hclust

  # build tSNE and explore cluster quality
  tsn <- Rtsne::Rtsne(as.dist(t(larry$bardist)), perplexity = 20)
  larry$tsn <- tsn
  if (do.plot){
    cols <- rainbow(length(unique(clust)))
    plot(tsn$Y, col = sample(cols)[clust], pch = 19, cex = 0.2)
    # double check that clusters are OK
    par(mfrow=c(4,4), mar = c(0.5, 0.5, 0.5, 0.5))
    for (i in 1:16){
      plot(tsn$Y, col = as.factor(clust == i), pch = 19, cex = 0.1,
           xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')

    }
  }

  # choose representative cluster barcodes and check their similarity
  larry$barcode.counts <- rowSums(larry$mat > 0)
  larry$cluster.barcodes <- unlist(lapply(1:length(unique(clust)), function(cl){
    names(which.max(larry$barcode.counts[clust == cl]))
  }))
  larry$cluster.abund <- unlist(lapply(1:length(unique(clust)), function(cl){
    sum(larry$barcode.counts[clust == cl])
  }))
  larry$clustdist <- adist(larry$cluster.barcodes, larry$cluster.barcodes)

  # collapse cluster barcodes:
  larry$mat.collapsed <- filter_out_barcodes(larry, #cutoff = cutoff, cutdist = cutdist,
                                             sum.barcodes = T)
  if (do.plot){
    par(mfrow=c(1,1), mar=c(4,4,1,1))
    plot(apply(larry$mat.collapsed, 2, sum), apply(larry$mat.collapsed, 2, max))
    abline(a=0,b=1,col='red')
  }

  # merge barcodes correlated across cells (likely were co-injected in the original cell)
  #barcor <- cor(t(larry$mat.collapsed > 1))
  barcor <- cor(t(larry$mat.collapsed))
  barcor[is.na(barcor)] <- 0
  hcl <- hclust(as.dist(1 - barcor), method = "single")
  cls <- cutree(hcl, h = cutoff.multi)
  larry$cutoff.multi <- cutoff.multi
  larry$cls <- cls
  larry$hcl = hcl
  if (do.plot){
    par(mar=c(4,4,1,1))
    plot(hcl, labels =  rep("", length(hcl$order)), cex = 0.5)
    abline(h=cutoff.multi, col = "red")
  }

  # collapse co-occuring barcodes in final cell-clone matrix
  abund <- rowSums(larry$mat.collapsed > 0)
  cls.barcodes <- unlist(lapply(1:length(unique(cls)), function(cl)
    names(which.max(abund[cls == cl]))))
  mat.clones <- do.call(rbind, mclapply(1:length(unique(cls)), function(cl){
    colSums(larry$mat.collapsed[cls == cl, , drop=F])
  }, mc.cores = 20))
  rownames(mat.clones) <- cls.barcodes
  larry$mat.clones <- mat.clones

  if (do.plot){
    par(mfrow=c(2,2), mgp = c(2,1,0), mar = c(4,4,0.1, 0.1))
    hist(colSums(mat.clones > 0), xlab = "distinct barcodes per cell", main = '')
    hist(colSums(mat.clones > 1), xlab = "distinct barcodes (>1 UMI) per cell", main = '')
    hist(colSums(mat.clones)+1e-1, plot = T, xlab = "barcode UMIs per cell", main = '')
    plot(apply(mat.clones, 2, function(x) max(x)), colSums(mat.clones) - apply(mat.clones, 2, function(x) max(x)),
         xlab = "major barcode UMIs", ylab = "non-major barcode UMIs", col = adjustcolor("black", 0.1), pch = 21, cex = 0.5)
    abline(a=0,b=1,col="blue")
    abline(v = 5, col = "green"); abline(a=0, b = 0.05, col = "red")
  }

  # call clones in individual cells
  cloneid <- rownames(mat.clones)[apply(mat.clones, 2, which.max)]
  names(cloneid) <- colnames(mat.clones)
  clone.call <- apply(mat.clones, 2, max) >= 5 &
    apply(mat.clones, 2, max)/colSums(mat.clones) > 0.85
  cloneid[!clone.call] <- NA

  larry$cloneid <- cloneid

  return(larry)
}


visualize_clone_processing <- function(larry, cexlab = 2, do.par = F){
  if (do.par){
    par(mfrow=c(3,2), mgp = c(2,1,0), mar = c(3, 3, 1, 1))
  }
  # show barcode sequence similarity
  #plot(larry$hclust.barcode, xlab = '', labels = FALSE, cex = 0.1, main = "barcode sequences",
  #     cex.main = cexlab, cex.lab = cexlab, font.main = 1, lwd = 0.2)

  hcd <- as.dendrogram(larry$hclust.barcode)
  plot(hcd, xlab = '', ylab = 'height', leaflab = "none", cex = 0.1, main = "barcode sequences",
       cex.main = cexlab, cex.lab = cexlab, font.main = 1, edgePar = list(lwd = 0.3))
  abline(h=larry$cutoff.hclust, col = "red")

  # show barcode correlation across cells
  #plot(larry$hcl, labels =  rep("", length(larry$hcl$order)), xlab = "", cex = 0.1,  main = "barcode co-expression",
  #     cex.main = cexlab, cex.lab = cexlab, font.main = 1, lwd = 0.3)
  #abline(h=larry$cutoff.multi, col = "red")

  hcd <- as.dendrogram(larry$hcl)
  plot(hcd, leaflab = "none", xlab = "", cex = 0.1,  main = "barcode co-expression",
       cex.main = cexlab, cex.lab = cexlab, font.main = 1, edgePar = list(lwd = 0.3))
  abline(h=larry$cutoff.multi, col = "red")


  # clonal summary:
  hist(colSums(larry$mat.clones > 1), breaks = 40, dborder = 'black', border.lwd = 3, xlim = c(0, 8),
       xlab = "clones (> 1 UMI) per cell", main = '',
       cex.lab = cexlab)
  hist(colSums(larry$mat.clones), plot = T, breaks = 500, log = 'x', xlab = "total barcode UMIs per cell", main = '',
       cex.lab = cexlab)
  plot(apply(larry$mat.clones, 2, function(x) max(x)), colSums(larry$mat.clones) - apply(larry$mat.clones, 2, function(x) max(x)),
       xlab = "major clone, UMIs", ylab = "non-major clone, UMIs", col = adjustcolor("black", 0.05), pch = 21, cex = 0.5,
       cex.lab = cexlab)
  abline(v = 5, col = "green", lwd = 0.5); abline(a=0, b = 0.15, col = "red", lwd = 0.5)
  # called clones
  hist(table(larry$cloneid), plot = T, xlab = "cells per clone", main = '',
       cex.lab = cexlab)
}

filter_out_barcodes <- function(larry, #cutoff = 4, cutdist = 2,
                                sum.barcodes = TRUE){
  #ind.filter <- larry$dist2center >= cutoff | dist2second - dist2center < cutdist
  #barcodes.retain1 <- larry$barcodes.retain
  #barcodes.retain1[larry$barcodes.retain][ind.filter] <- FALSE
  #clust1 <- larry$clust[!ind.filter]
  #mat.retained <- larry$mat[barcodes.retain1, ]

  clust1 <- larry$clust
  mat.retained <- larry$mat
  mat.collapsed <- do.call(rbind, mclapply(1:length(unique(clust1)), function(cl){
    if (sum.barcodes){
      colSums(mat.retained[clust1 == cl, ,drop = F])
    }else{
      apply(mat.retained[clust1 == cl, ,drop = F], 2, max)
    }
  }, mc.cores = 20))
  rownames(mat.collapsed) <- larry$cluster.barcodes
  return(mat.collapsed)
}

hamming_barcodes <- function(larry, n.cores = 20){
  mat <- larry$mat
  bardist <- do.call(cbind, mclapply(1:nrow(mat), function(i){
    as.numeric(adist(rownames(mat)[i],rownames(mat)))
  }, mc.cores = n.cores))
  rownames(bardist) <- colnames(bardist) <- rownames(mat)
  return(bardist)
}

create_clone_matrix <- function(larry){
  trios <- larry$trios
  cls <- unique(trios$X.cell_id)
  brcd <- unique(trios$clone_id)
  mat <- do.call(cbind, mclapply(cls, function(cl){
    x <- table(trios[trios$X.cell_id == cl, 2])
    x1 <- rep(0, length(brcd)); names(x1) <- brcd
    x1[names(x)] <- x
    return(x1)
  }, mc.cores = 20))
  rownames(mat) <- brcd
  colnames(mat) <- paste('one_', cls, '-1', sep = '')
  colnames(mat) <- cls
  return(mat)
}

upload_trios <- function(fpath, pref = ''){
  trios <- read.csv(fpath, sep = '\t')
  trios <- unique(trios)
  trios <- trios[table(trios$umi) == 1, ]
  trios <- trios[,c(1,3)]
  trios[,1] <- paste(pref, trios[,1], sep = '')
  #trios[,2] <- paste(pref, trios[,2], sep = '')
  return(trios)
}



