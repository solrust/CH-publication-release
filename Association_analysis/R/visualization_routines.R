plot_associations <- function(de, pv, pheno.name = NULL, cex.text = 2){
  if (is.null(pheno.name)){
    pheno.name <- 'phenotype'
  }
  par(mfrow=c(1,3), mar=c(5,5,2,2), mgp = c(3,1,0))
  plot(de, pv, xlim = c(max(-10, min(de, na.rm = T)), min(10, max(de, na.rm = T))),
       log = 'y', main = pheno.name, font.main = 1, xlab = 'DE, log', ylab = 'p-value',
       pch = 19, cex = 1, col = adjustcolor('darkgrey', 0.5),
       cex.lab = cex.text, cex.axis = cex.text,  cex.main = cex.text)
  abline(h=0.05,v=0, lty = 2)
  if (exists('tf.list')){
    for (i in 1:length(tfs.list)){
      gns.int <- intersect(tfs.list[[i]], rownames(test.basic))
      points(de[gns.int], pv[gns.int],
             pch = 19, cex = 1, col = tfs.cols[i])
    }
    legend('bottomright', names(tfs.list), col = tfs.cols, pch = 19, bty = 'n')
  }
  #
  hist(pv, breaks = 100, col = 'grey', cex.lab = cex.text, cex.axis = cex.text, xlab = 'p-value', ylab = '', main = '')
  fdr <- p.adjust(pv, 'BH')
  sfdr <- sort(fdr); scum <- (cumsum(sfdr))
  plot(c(0, sfdr), c(0, 1:length(sfdr)), xlim = c(-0.01, 0.5), ylim = c(0, length(sfdr[sfdr <= 0.5])+0), #ylim = c(0, max(scum[sfdr < 0.5])),
       pch = 19, cex = 0.1, type = 's', lwd = 2,
       xlab = 'FDR', ylab = 'genes', cex.lab = cex.text, cex.axis = cex.text)
  abline(v= seq(0, 0.5, length.out = 4), lty = 2, lwd = 0.5)
  abline(h= seq(0, length(sfdr[sfdr <= 0.5]), length.out = 4), lty = 2, lwd = 0.5)
  #abline(h=c(max(scum[sfdr < 0.1]), max(scum[sfdr < 0.25])), lty = 2)
}


show_gene <- function(mat = NULL, gn, clone.resid, clone.states, sign.clones = T, title = gn, xname = 'clone score', yname =  'mRNA abundance',
                      text.size = 15){
  if (is.null(mat)){
    mat <- expr.hsc[,]
  }
  DE.clone <- ifelse(sign.clones, 1, 0.01)
  Geno <- factor(clone.states$genotype); Exp <- clone.states$exp
  ft <- lm(mat[sign.clones, gn] ~ clone.resid[sign.clones])
  pl <- qplot(clone.resid, mat[, gn], xlab = xname, ylab = yname, #log = 'y',
              ylim = c(min(0, min(mat[sign.clones, gn], na.rm=T)), max(mat[sign.clones, gn], na.rm=T)),
              #alpha = DE.clone,
              #colour = Exp,
              #fill = I('white'),
              shape = Exp,
              #shape = as.factor(Genotype),
              size = I(3),
              stroke = I(0.5),
              main = title
  )+ theme_bw() + geom_point(aes(shape = Exp, colour = Geno), #colour = col.genotype,
                             size = 2.2, alpha = I(0.8)) +
    geom_smooth(aes(colour = Geno), method = 'lm', size = I(1), alpha = I(0.1))+
    scale_color_manual(values = adjustcolor(col.genotype, 1)) +
    theme(axis.text = element_text(size = text.size),
          axis.title = element_text(size = text.size),
          legend.text = element_text(size = text.size),
          plot.title = element_text(size = text.size)#, legend.title = element_blank(),
    ) # + scale_color_manual(values = adjustcolor(col.genotype, 1))
  #pl <- pl + ggrepel::geom_text_repel(aes(clone.resid, mat[, gn], label = 'er'), size = 2, max.overlaps = 100)
  return(pl)

}
