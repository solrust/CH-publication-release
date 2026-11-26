summary_mut_clone <- function(tbl, clone.info){
  mut.clone.stat.de <- do.call(rbind, lapply(rownames(clone.info), function(cln){
    print(cln)
    #cln <- rownames(clone.info)[nrow(clone.info) - 30]
    cls <- rownames(cell.clone)[cell.clone$cloneid == cln]
    cls1 <- intersect(cls, rownames(tbl))
    #str(cls1)
    print(length(cls1))
    ind <- tbl[cls1,]$use == 'train'
    #table(ind)
    if (length(cls1) >= 1){
      #predx1 <- mean(tbl[cls1[ind], ]$score)
      #predx2 <- mean(tbl[cls1[!ind], ]$score)
      predx1 <-  mean(tbl[cls1[ind], ]$qmut)
      predx2 <-  mean(tbl[cls1[!ind], ]$qmut)
      if (sum(!ind) > 1 & sum(ind) > 1){
        tst <- ks.test(tbl[cls1[ind], ]$qmut, "punif", 0, 1)
        pv <- tst$p.value

        x <- runif(length(cls1[ind]), 0, 1)
        tst1 <- wilcox.test(tbl[cls1[ind], ]$qmut, x)
        pv1 <- tst1$p.value

        tst2 <- ks.test(tbl[cls1[!ind], ]$qmut, "punif", 0, 1)
        pv2 <- tst2$p.value
      }else{
        pv <- NA
        pv1 <- NA
        pv2 <- NA
      }

      #predx2 <-  mean(tbl[cls1[ind], ]$qwt)
      return(c(predx1, predx2, pv, pv1, pv2))

    }else{
      rep(c(NA, NA))
    }
  }))
  rownames(mut.clone.stat.de) <- rownames(clone.info)
  return(mut.clone.stat.de)
}
