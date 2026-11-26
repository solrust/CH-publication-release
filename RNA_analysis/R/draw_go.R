draw_go <- function(matcl, feature.var0, setclust.go, cols.go,
                    n.terms = 15, only.go = FALSE, term.length = 70, cut.pval = 40){
  pvmin <- apply(abs(matcl[feature.var0,]), 1, max)
  xx <- do.call(rbind, lapply(1:length(unique(setclust.go[,1])), function(cl){
    print(cl)
    cl1 <- unique(setclust.go[,1])[cl]
    nms <- feature.var0[setclust.go[feature.var0, ] == cl1]
    if (only.go) nms <- nms[grepl('GO:', nms)]
    nms <- nms[nchar(nms) < term.length]
    xx <- as.data.frame(sort(pvmin[nms], decreasing = T)[1:min(n.terms, length(nms))])
    xx <- cbind(xx, cl)
    colnames(xx) <- c('pvmin', 'cl')
    if (cl < length(unique(setclust.go[,1]))) xx <- rbind(xx, 'blank' = NA)
    xx
  }))
  xx[!is.na(xx[,1]) & xx[,1] > cut.pval,1] <- cut.pval
  yy <- rownames(xx)
  pos.blank <- grepl('blank', yy)
  yy[pos.blank] <- strrep(" ", 1:sum(pos.blank))
  yy <- factor(yy, levels = rev(yy))
  pl <- ggplot(xx, aes(yy, pvmin) ) +
    geom_bar(position = position_dodge(), stat = "identity", bg = cols.go[xx$cl], col = adjustcolor('black', 0.2)) +
    #geom_hline(aes(yintercept = yy), size = 1.5, color="red", group = 1) +
    geom_vline(aes(xintercept = yy), size = 1, color=ifelse(pos.blank, 'black', NA), group = 1) +
    coord_flip() + theme_bw() +  theme(text = element_text(size = 20)) +
    scale_x_discrete(name = "", position = "top") +
    scale_y_continuous(name = "p-value, -log", #limits = c(0, cut.pval),
                       breaks = seq(0, 60, by = 10), labels = seq(0,  60, by =  10))
  return(pl)
}

