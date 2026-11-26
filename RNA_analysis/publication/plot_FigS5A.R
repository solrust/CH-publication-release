# Plot Supplementary Fig. 5A from exported count pairs (no raw counts needed).
stopifnot(file.exists('RNA_analysis/publication/source_data/FigS5A_hashtag_count_pairs.tsv.gz'))
x <- read.delim(gzfile('RNA_analysis/publication/source_data/FigS5A_hashtag_count_pairs.tsv.gz'))
outdir <- 'RNA_analysis/publication/figures'
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
pdf(file.path(outdir, 'FigS5A_reproduced.pdf'), width=9, height=6)
par(mfrow=c(2,3), mar=c(4.3,4.3,2.1,1), mgp=c(2.6, 0.8, 0))
for (j in sort(unique(x$panel_order))) {
  z <- x[x$panel_order==j, ]
  # Original base-R log axes omit nonpositive values; keep them in Source Data.
  keep <- is.finite(z$x_count) & is.finite(z$y_count) & z$x_count>0 & z$y_count>0
  plot(z$x_count[keep], z$y_count[keep], log='xy', xlim=c(1,10000),
       ylim=c(1,10000), pch=19, cex=.15, col=adjustcolor('black', alpha.f=.12),
       xlab=paste('hashtag', z$hashtag_x[1]),
       ylab=paste('hashtag', z$hashtag_y[1]), main=paste('RNA_', z$RNA_run[1]),
       xaxt='n', yaxt='n')
  axis(1, at=c(1,100,10000)); axis(2, at=c(1,100,10000))
  abline(h=z$threshold[1], v=z$threshold[1], col='red')
}
dev.off()
cat('Created ', file.path(outdir, 'FigS5A_reproduced.pdf'), '\n')
