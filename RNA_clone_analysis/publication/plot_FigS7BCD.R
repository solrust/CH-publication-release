# Replot Supplementary Figure 7B-D from compact Source Data.
# Requires ggplot2; run from the CH_publication_release repository root.
local({
  if (!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2')
  library(ggplot2)
  pub <- 'RNA_clone_analysis/publication'
  src <- file.path(pub,'source_data')
  out <- file.path(pub,'figures'); dir.create(out,showWarnings=FALSE,recursive=TRUE)
  need <- c('FigS7B_dominant_clone_umis.tsv.gz','FigS7C_clone_sizes_by_sample.tsv',
            'FigS7D_clone_counts_by_sample.tsv')
  if (any(!file.exists(file.path(src,need)))) stop('Run export_FigS7BCD.R first')

  # Panel B: preserve the original 0.85x diagonal and five-UMI cutoff.
  b <- read.delim(gzfile(file.path(src,need[1])),check.names=FALSE)
  set.seed(239)
  b <- b[sample.int(nrow(b)),,drop=FALSE]
  plB <- ggplot(b,aes(total_clone_umis,dominant_clone_umis)) +
    geom_point(aes(color=library),alpha=0.01,size=2) +
    geom_abline(slope=0.85,linetype=2,color='black') +
    geom_vline(xintercept=5,linetype=2,color='black') +
    coord_cartesian(xlim=c(0,250),ylim=c(0,250),expand=FALSE) +
    labs(x='all clones, UMIs',y='dominant clone, UMIs') +
    theme_bw() + theme(text=element_text(size=14),legend.position='none')
  png(file.path(out,'FigS7B_reproduced.png'),width=700,height=700,res=175)
  print(plB);dev.off()

  # Panel C: logarithmic clone-size histogram, one facet per sample.
  cdat <- read.delim(file.path(src,need[2]),check.names=FALSE)
  cdat$genotype <- factor(cdat$genotype,levels=c('DNMT3A','TET2','WT'))
  pal <- c('DNMT3A'='#EF8F74','TET2'='#A8D69A','WT'='#777777')
  plC <- ggplot(cdat,aes(x=clone_size,fill=genotype)) +
    geom_histogram(bins=10,color='black',linewidth=0.12) +
    scale_x_log10() + scale_fill_manual(values=pal,drop=FALSE) +
    facet_wrap(~sample_id+expid,ncol=2,dir='v') +
    labs(x='size',y='count',fill=NULL) +
    theme_bw() + theme(strip.text=element_blank(),legend.position='top',
                       text=element_text(size=11))
  pdf(file.path(out,'FigS7C_reproduced.pdf'),width=7,height=10)
  print(plC);dev.off()

  # Panel D: sample-level clone counts with genotype and experiment indicators.
  d <- read.delim(file.path(src,need[3]),check.names=FALSE)
  xpos <- match(d$genotype,c('DNMT3A','WT','TET2'))
  if (anyNA(xpos)) stop('Unexpected genotype in S7D data')
  # Historical script horizontally separates WT observations from the two experiments.
  xpos[d$genotype=='WT' & d$expid=='123'] <- xpos[d$genotype=='WT' & d$expid=='123']-0.2
  xpos[d$genotype=='WT' & d$expid=='127'] <- xpos[d$genotype=='WT' & d$expid=='127']+0.2
  bg <- c('DNMT3A'='#F18D91','WT'='#686868','TET2'='#A2EA89')[d$genotype]
  pdf(file.path(out,'FigS7D_reproduced.pdf'),width=4,height=4)
  par(mar=c(5,5,1,1))
  plot(xpos,d$n_clones,xlim=c(0.8,3.2),ylim=c(0,max(d$n_clones)+1),
       pch=ifelse(d$expid=='123',21,24),bg=bg,col='black',cex=1.8,
       xlab='',ylab='scRNA-seq clones',xaxt='n',bty='n',cex.lab=1.3,cex.axis=1.3)
  axis(1,at=1:3,labels=c('DNMT3A','WT','TET2'),cex.axis=1.2)
  legend('bottomright',c('exp1','exp2'),pch=c(21,24),bty='n',cex=1)
  dev.off()
  message('Created S7B, S7C and S7D review plots in ',out)
})
