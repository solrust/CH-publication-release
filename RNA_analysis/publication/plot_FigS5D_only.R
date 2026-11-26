# S5D-only renderer. Original S5D code extracted unchanged from plot_FigS5CD.R.
# S5C is explicitly deferred; no S5C table is read or rendered here.
outdir <- 'RNA_analysis/publication/figures'
sourcedir <- 'RNA_analysis/publication/source_data'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
files <- list.files(sourcedir,'^FigS5D_program_[0-9]+_scores\\.tsv\\.gz$',full.names=TRUE)
if (!length(files)) stop('No FigS5D program scores found.')
files <- sort(files)
pdf(file.path(outdir,'FigS5D_reproduced_review.pdf'),width=12,height=12)
par(mfrow=c(4,4),mar=c(.6,.6,2,.6))
for (f in files) {
  x <- read.delim(gzfile(f))
  pal <- colorRampPalette(c('grey90','grey80','red'))(20)
  col <- pal[as.integer(cut(x$display_score, breaks=20))]
  plot(x$umap_1,x$umap_2,col=col,pch=19,cex=.15,
       xaxt='n',yaxt='n',xlab='',ylab='',main=x$program[1],asp=1)
}
if (length(files)<16L) for (i in seq_len(16L-length(files))) plot.new()
dev.off()
cat('Created S5D program maps; S5C not executed.\n')
