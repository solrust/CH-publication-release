# Reproduce Fig. 4B by combining the shared Fig. 3A UMAP background
# with Fig. 4B's *highlighted cells only* Source Data.
# Run from CH_publication_release root after export_Fig4B.R.

umap_file <- 'RNA_analysis/publication/source_data/Fig3A_umap_cells.tsv.gz'
clone_file <- 'RNA_clone_analysis/publication/source_data/Fig4B_example_clone_cells.tsv.gz'
for (f in c(umap_file, clone_file)) if (!file.exists(f)) stop('Missing ', f)

background_raw <- utils::read.delim(umap_file, check.names=FALSE, stringsAsFactors=FALSE)
# fig4b_find_umap_cols is shared with the exporter, but this plotter is standalone.
norm <- tolower(gsub('[^[:alnum:]]', '', names(background_raw)))
get_col <- function(keys, label) {
  m <- match(keys, norm, nomatch=0L)
  m <- m[m > 0L]
  if (!length(m)) stop('Cannot find ', label, ' in Fig3A; columns: ', paste(names(background_raw),collapse=', '))
  names(background_raw)[m[1L]]
}
x <- as.numeric(background_raw[[get_col(c('umap1','umapx','embedding1','emb1','dimension1','dim1','x'),'UMAP1')]])
y <- as.numeric(background_raw[[get_col(c('umap2','umapy','embedding2','emb2','dimension2','dim2','y'),'UMAP2')]])
selected <- utils::read.delim(clone_file, stringsAsFactors=FALSE)
expected_cols <- c('UMAP_1','UMAP_2','clone_number','reported_clone_size')
if (!all(expected_cols %in% names(selected))) stop('Fig4B Source Data missing required columns.')

outdir <- 'RNA_clone_analysis/publication/figures'
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
outfile <- file.path(outdir,'Fig4B_reproduced.png')
cols <- c('16'='#0000FF', '213'='#8000FF', '396'='#FFA500')
titles <- c('16'='myeloid biased','213'='lymphoid biased','396'='HSC biased')
grDevices::png(outfile, width=1500, height=500, res=150)
tryCatch({
  op <- par(mfrow=c(1,3),mar=c(1,1,4,1))
  for (number in names(cols)) {
    sub <- selected[as.character(selected$clone_number) == number, , drop=FALSE]
    plot(x,y,pch=19,cex=0.10,col='grey75',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',
         main=titles[[number]],font.main=1,cex.main=1.15)
    points(sub$UMAP_1,sub$UMAP_2,pch=21,cex=0.6,lwd=0.8,col=cols[[number]],bg='white')
    if (nrow(sub) > 0L) {
      legend('topright',legend=paste0('clone ',number,'\nsize = ',sub$reported_clone_size[1]),
             bty='n',cex=0.85)
    }
  }
  par(op)
},finally=grDevices::dev.off())
message('Created ',outfile)
