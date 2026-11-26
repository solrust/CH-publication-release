make_plots <- function(motif_stat){
  p1 <- qplot(motif_stat$fgseal[[1]]$NES, motif_stat$fgseal[[1]]$padj, log = 'y', label = rownames(motif_stat$fgseal[[1]])) + theme_bw() +
    labs(x = 'NES', y = 'p-adj', title = 'GSEA, all peaks') + geom_text_repel() # + theme(text = element_text(size = 20))
  p2 <- qplot(motif_stat$fgseal[[2]]$NES, motif_stat$fgseal[[2]]$padj, log = 'y', label = rownames(motif_stat$fgseal[[2]])) + theme_bw() +
    labs(x = 'NES', y = 'p-adj', title = 'GSEA, peaks p<=0.1') + geom_text_repel() # + theme(text = element_text(size = 20))
  gg <- intersect(rownames(motif_stat$fgseal[[1]]), rownames(motif_stat$fgseal[[2]]))
  p3 <- qplot(motif_stat$fgseal[[1]][match(gg, rownames(motif_stat$fgseal[[1]])),]$NES, motif_stat$fgseal[[2]][match(gg, rownames(motif_stat$fgseal[[2]])),]$NES, label = gg) + theme_bw() +
    labs(x = 'NES, all', y = 'NES, p <= 0.1') + geom_text_repel() # + theme(text = element_text(size = 20))
  #
  p4 <-  qplot(motif_stat$chisql[[1]][,1], p.adjust(motif_stat$chisql[[1]][,2], 'BH'), log = 'y', label = rownames(motif_stat$chisql[[1]])) + theme_bw() +
    labs(x = 'enrichment', y = 'p-adj', title = 'peaks all') + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() # + theme(text = element_text(size = 20))
  p5 <-  qplot(motif_stat$chisql[[2]][,1], p.adjust(motif_stat$chisql[[2]][,2], 'BH'), log = 'y', label = rownames(motif_stat$chisql[[2]])) + theme_bw() +
    labs(x = 'enrichment', y = 'p-adj', title = 'peaks p<=0.1') + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() # + theme(text = element_text(size = 20))
  p6 <- qplot(motif_stat$chisql[[1]][,1], motif_stat$chisql[[2]][,1], label = rownames(motif_stat$chisql[[1]])) + theme_bw() +
    labs(x = 'enrichment, all', y = 'entichment, p<=0.1') + geom_text_repel() # + theme(text = element_text(size = 20))
  #
  gg <- intersect(rownames(motif_stat$fgseal[[1]]), rownames(motif_stat$chisql[[1]]))
  p7 <- qplot(motif_stat$fgseal[[1]][match(gg, rownames(motif_stat$fgseal[[1]])),]$NES, motif_stat$chisql[[1]][gg,1], label = gg) + theme_bw() +
    labs(x = 'GSEA', y = 'entichment', title = 'all peaks') + geom_text_repel() # + theme(text = element_text(size = 20))
  gg <- intersect(rownames(motif_stat$fgseal[[2]]), rownames(motif_stat$chisql[[2]]))
  p8 <- qplot(motif_stat$fgseal[[2]][match(gg, rownames(motif_stat$fgseal[[2]])),]$NES, motif_stat$chisql[[2]][gg,1], label = gg) + theme_bw() +
    labs(x = 'GSEA', y = 'entichment', title = 'peaks p<=0.1') + geom_text_repel() # + theme(text = element_text(size = 20))
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, p7 = p7, p8 = p8))
}


cells_partition <- function(clone.info.multi,
                            pheno.name = NULL, size.cut = 50, batch_name = 'id', test.geno = TRUE){
  # alternative batch_name = 'expgeno'

  # identify high/low clones
  if (!test.geno){
    pheno.score <- clone.info.multi[, pheno.name]
    #pheno.score <- clone.info.multi[, '10'] - clone.info.multi[, '8']
    med <- median(pheno.score[clone.info.multi$size >= size.cut], na.rm = T)
    clones.high <- clone.info.multi$larry[which(pheno.score >= med & clone.info.multi$size >= size.cut)]
    clones.low <- clone.info.multi$larry[which(pheno.score < med & clone.info.multi$size >= size.cut)]
    cells.high <- rownames(cell.clone.multi)[cell.clone.multi$larry %in% clones.high]
    cells.low <- rownames(cell.clone.multi)[cell.clone.multi$larry %in% clones.low]
  }else{
    if (pheno.name == 'DNMT3A'){
      expid <- '123'
      geno <- 'DNMT3A'
    }else{
      expid <- '127'
      geno <- 'TET2'
    }
    cells.high <- rownames(cell.clone.multi)[cell.clone.multi$genotype %in% geno & grepl(expid, cell.clone.multi$id)]
    cells.low <- rownames(cell.clone.multi)[cell.clone.multi$genotype %in% 'WT' & grepl(expid, cell.clone.multi$id)]
  }
  cell_conditions <- list(cells.low, cells.high)
  return(cell_conditions)
}




perform_da <- function(obj.combined, cell_conditions = NULL,
                       modality = NULL,  count.mat = NULL,
                       cell.info = NULL, batch_name = NULL,
                       do.deseq2 = TRUE, shrinkage = FALSE,
                       seurat_class = 'seurat_clusters', clusts = NULL, old.dataset = FALSE){
  if (is.null(clusts)) stop('provide clusters')
  if (!old.dataset){
    if (!is.null(modality)){
      if (modality == 'ATAC'){
        count.mat <- obj.combined@assays$ATAC$counts
      }else if (modality == 'RNA'){
        count.mat <- obj.combined@assays$RNA$counts
      }
    }
    Idents(obj.combined) <- seurat_class
    p2clust1 <- Idents(obj.combined)
  }else{
    p2clust1 <- old.scrna$clust.final
    count.mat <- old.scrna$cd
  }

  # calculate pseudo-bulk condition-specific matrices
  print('calculate pseudo-bulk matrices')
  clust.clone.list <- lapply(1:2, function(i){
    print(i)
    cells <- cell_conditions[[i]]
    matl <- lapply(clusts, function(clust){
      print(paste('clust', clust))
      print('p2clust1')
      print(str(p2clust1))
      cls <- names(p2clust1)[!is.na(p2clust1) & p2clust1 %in% clust]
      print('cls')
      print(str(cls))
      print('cells')
      print(str(cells))
      cls <- intersect(cls, cells)
      smpls <- unique(cell.info[cls, batch_name])
      print('cell.info')
      print(str(cell.info))
      print('done')
      print('cls')
      print(str(cls))
      print('done')
      print('batch_name')
      print(str(batch_name))
      smpls <- unique(as.character(cell.info[cls, batch_name]))
      if (old.dataset) smpls <- unique(cell.info[cls,]$id)
      smpls <- smpls[!is.na(smpls)]
      print('here')
      print(smpls)
      mat <- do.call(cbind, lapply(smpls, function(smpl){
        print(paste('smpl', smpl))
        if (old.dataset){
          ind <- cell.info[cls,]$id %in% smpl
        }else{
          ind <- cell.info[cls, batch_name] %in% smpl
        }
        if (sum(ind) > 0){
          rowSums(count.mat[, cls][, ind, drop = F])
        }else{
          rep(0, nrow(count.mat))
        }
      }))
      print(str(mat))
      colnames(mat) <- smpls

      return(mat)
    })
    names(matl) <- clusts
    return(matl)
  })

  # combine in a single matrix and record metadata
  mat <- do.call(cbind, do.call(cbind, clust.clone.list))
  #
  hash.match1 <- do.call(rbind, lapply(1:length(clust.clone.list), function(i){
    do.call(rbind, lapply(1:length(clust.clone.list[[i]]), function(j){
      cbind(colnames(clust.clone.list[[i]][[j]]), i, names(clust.clone.list[[i]])[j])
    }))
  }))
  hash.match1 <- data.frame(hash.match1)
  colnames(hash.match1) <- c('batch', 'clone', 'clust')
  hash.match1$clone <- as.factor(hash.match1$clone)
  hash.match1$clust <- as.factor(hash.match1$clust)
  colnames(mat) <- rownames(hash.match1)

  ind <- colSums(mat) > 0
  mat <- mat[, ind]
  hash.match1 <- hash.match1[ind,]

  # perform DESeq2 analysis
  if (do.deseq2){
    print('perform DESeq2 analysis')
    ind1 <- TRUE
    #print(hash.match1[ind1,])
    #print(dim(mat[,ind1]))
    #print(length(length(unique(hash.match1$clust))))
    #print(hash.match1$clone)
    if (length(unique(hash.match1$clust)) > 1){
      dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                        design = ~ clone + clust)
    }else{
      dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                        design = ~ clone)
    }
    dds_mut <- DESeq(dds_mut)
    print(resultsNames(dds_mut))
    if (shrinkage){
      res <- lfcShrink(dds_mut, coef="clone_2_vs_1", type="normal")
    }else{
      res <-  results(dds_mut, name="clone_2_vs_1")
    }
  }else{
    dds_mut <- NULL; res <- NULL
  }

  return(list(mat = mat, hash.match1 = hash.match1, dds = dds_mut, res = res))
}


estimate_motif_enrichment <- function(obj.combined, da_analysis,  motif.matrix = NULL, motif.names = NULL, pcuts = c(1, 1e-1),
                                      permute = FALSE, run_peaks = TRUE,
                                      features.match = c('GC.percent', "percentile"), run_fgsea = FALSE){
  res <- da_analysis$res
  # run fgsea for motif enrichment
  if (run_fgsea){
    print('perform fgsea for motif enrichment')
    fgseal <- lapply(pcuts, function(pcut){
      print(pcut)
      pa.scores <- res$log2FoldChange
      if (permute) pa.scores <- sample(pa.scores)
      names(pa.scores) <- rownames(res)
      pa.scores <- pa.scores[!is.na(pa.scores) & rank(res$pvalue) <= pcut * sum(!is.na(res$pvalue)) ]
      pks.sets <- lapply(1:ncol(motif.matrix), function(i){
        pks <- names(which(motif.matrix[, i]))
        pks <- intersect(names(pa.scores), pks)
        pks
      })
      names(pks.sets) <- colnames(motif.matrix)
      fgseaRes <- fgsea::fgsea(pathways = pks.sets, pa.scores,
                               minSize=20, maxSize = 30000)#, nperm=1000)
      rownames(fgseaRes) <- unlist(motif.names)[fgseaRes$pathway]

      return(fgseaRes)
    })
  }else{
    fgseal <- NULL
  }

  # run chisq test for motif enrichment
  print('perform fisher test for motif enrichment')
  meta.feature <- GetAssayData(obj.combined, assay = "ATAC", slot = "meta.features")
  open.peaks <- AccessiblePeaks(obj.combined, idents = clusts)
  # TODO: use all 20 cores?
  chisql <- lapply(pcuts, function(pcut){
    scores <- do.call(rbind, mclapply(1:ncol(motif.matrix), function(i){
      if (i %% 100 == 0) print(i)
      pa.scores <- res$log2FoldChange
      if (permute) pa.scores <- sample(pa.scores)
      names(pa.scores) <- rownames(res)
      pa.scores <- pa.scores[!is.na(pa.scores) &  rank(res$pvalue) <= pcut * sum(!is.na(res$pvalue)) ]#res$pvalue <= pcut]
      #
      pks <- names(which(motif.matrix[, i]))
      pks <- intersect(names(pa.scores), pks)
      peaks.matched <- setdiff(names(pa.scores), pks)
      if (run_peaks){
        peaks.matched <- MatchRegionStats(
          meta.feature = meta.feature[open.peaks,],
          query.feature = meta.feature[pks, ],
          features.match = features.match, n = 50000)
      }
      polars.tf <- table(as.factor(pa.scores > 0)[pks])
      polars <- table(pa.scores[peaks.matched] > 0)
      tst <- fisher.test(cbind(polars.tf, polars))
      scores <- (polars.tf[2] / polars.tf[1]) / (polars[2] / polars[1])
      c(scores, tst$p.value, polars.tf, polars)
    }, mc.cores = ifelse(run_peaks, 20, 1)))
    rownames(scores) <- unlist(motif.names)[colnames(motif.matrix)]
    return(scores)
  })


  return(list(fgseal = fgseal, chisql = chisql))
}



prepare_table <- function(fc, pv, row.subset = NULL, col.subset = NULL){
  if (!is.null(row.subset)) {
    fc <- fc[row.subset,]
    pv <- pv[row.subset,]
  }
  if (!is.null(col.subset)) {
    fc <- fc[, col.subset]
    pv <- pv[, col.subset]
  }
  fc.hl <- hclust(dist(fc))
  fc.hl1 <- hclust(dist(t(fc)))
  print(fc.hl$order)
  print(fc.hl1$order)
  #
  tbl <- data.table::melt(fc, na.rm = F)
  tbl1 <- data.table::melt(pv, na.rm = F)
  tbl <- cbind(tbl, p = tbl1$value)
  tbl$state <-  ifelse(tbl$value > 0, 'up', 'down')
  colnames(tbl) <- c('module', 'phenotype', 'value', 'p', 'state')
  tbl$module <- factor(tbl$module, levels = rownames(fc)[fc.hl$order])
  tbl$phenotype <- factor(tbl$phenotype, levels = colnames(fc)[fc.hl1$order])
  return(tbl)
}


