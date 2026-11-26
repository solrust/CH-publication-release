project_gset <- function(gset, mat, pr = NULL){
  gset <- intersect(gset, colnames(mat))
  gset <- gset[colSums(mat[,gset] > 0) > 0]
  if (length(gset) > 0){
    av <- rowMeans(mat[, gset])
    if (is.null(pr)){
      if (!is.null(mat)){
        pr <- prcomp(mat[, gset], center = TRUE, scale. = TRUE)
        gscore <- pr$rotation[,1]
        pr <- pr$x
      }else{
        stop('requires pr or mat')
      }
    }else{
      gscore <- NULL
    }
    ft <- lm(av ~ pr)
    av.fit <- ft$fitted.values

    return(list(av = av, av.fit = av.fit, gscore = gscore))
  }else{
    NULL
  }
}
