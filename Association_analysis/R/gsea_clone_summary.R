gsea_clone_summary <- function(test.basic.list, gene_sets, genes = NULL, do.DE = FALSE, do.permute = FALSE, minSize = 15, maxSize = 1000){
  print('gscore.clone..')
  gscore.clone <- do.call(cbind, lapply(test.basic.list, function(test.basic.nm){
    if (do.DE){
      gscore <- test.basic.nm[,'pheno']
    }else{
      gscore <- -log(test.basic.nm[,'p.pheno']) * sign(test.basic.nm[,'pheno'])
    }
    if (do.permute){
      gscore <- sample(gscore)
    }
    names(gscore) <- rownames(test.basic.nm)
    if (!is.null(genes)){
      gscore <- gscore[intersect(genes, names(gscore))]
    }
    gscore[is.na(gscore)] <- 0
    return(gscore)
  }))
  #
  print("gsea.clone")
  gsea.clone <- mclapply(1:length(test.basic.list), function(i){
    test.basic.nm <- test.basic.list[[i]]
    gscore <- gscore.clone[,i]
    #gscore <- -log(test.basic.nm[,'p.pheno']) * sign(test.basic.nm[,'pheno'])
    #names(gscore) <- rownames(test.basic.nm)
    #gscore[is.na(gscore)] <- 0
    fgseaRes <- fgsea(pathways = gene_sets, stats  = gscore, minSize  = minSize, maxSize  = maxSize)
    return(fgseaRes)
  }, mc.cores = length(test.basic.list))

  #
  print("gsea.score")
  gsea.score <- do.call(cbind, lapply(gsea.clone, function(fgseaRes){
    x <- -log(unlist(fgseaRes[,2])) * sign(unlist(fgseaRes[,5]))
    names(x) <- as.character(fgseaRes$pathway)
    x
  }))
  return(list(gscore.clone = gscore.clone, gsea.clone = gsea.clone, gsea.score = gsea.score))
}
