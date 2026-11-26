# Initial review plots from archived source-data tables, NOT model refits.
local({
  if(!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2')
  library(ggplot2)
  pub <- 'RNA_clone_analysis/publication'
  file <- file.path(pub,'source_data','Fig5GH_FigS8E_state_estimates.tsv')
  if(!file.exists(file)) stop('Run export_Fig5GH_S8E.R first.')
  out <- file.path(pub,'figures')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  tbl <- read.delim(file,check.names=FALSE,stringsAsFactors=FALSE)
  tbl$contrast <- factor(tbl$contrast,levels=c('TET2','DNMT3A'))
  tbl$state_col <- factor(tbl$state,levels=1:10)
  # Review palette only. Published Illustrator cell-state palette can be substituted later.
  state_cols <- setNames(grDevices::hcl.colors(10,'Dark 3'),as.character(1:10))
  make_effect <- function(prefix,outfile,label) {
    d <- data.frame(contrast=tbl$contrast,state=tbl$state,state_col=tbl$state_col,
                    wt=tbl[[paste0(prefix,'.wt')]],
                    mut=tbl[[paste0(prefix,'.mut')]],
                    wt_lo=tbl[[paste0(prefix,'.wt.down')]],
                    wt_hi=tbl[[paste0(prefix,'.wt.up')]],
                    mut_lo=tbl[[paste0(prefix,'.mut.down')]],
                    mut_hi=tbl[[paste0(prefix,'.mut.up')]])
    required <- c('wt','mut','wt_lo','wt_hi','mut_lo','mut_hi')
    if(any(!is.finite(as.matrix(d[required])) | as.matrix(d[required])<=0))
      stop('Nonpositive / nonfinite values in ',prefix,' estimates: cannot use log axes without review')
    pl <- ggplot(d,aes(wt,mut,color=state_col)) +
      geom_segment(aes(x=wt_lo,xend=wt_hi,y=mut,yend=mut),linewidth=0.6,alpha=0.65) +
      geom_segment(aes(x=wt,xend=wt,y=mut_lo,yend=mut_hi),linewidth=0.6,alpha=0.65) +
      geom_abline(slope=1,intercept=0,linetype='dashed',color='grey60') +
      geom_point(size=2.5) +
      geom_text(aes(label=state),size=3,vjust=-0.8,show.legend=FALSE,check_overlap=TRUE) +
      facet_wrap(~contrast,nrow=1,scales='free') +
      scale_x_log10() + scale_y_log10() +
      scale_color_manual(values=state_cols) +
      labs(x=paste(label,'WT'),y=paste(label,'mutant'),color='Cell state') +
      theme_bw(base_size=12)
    pdf(file.path(out,outfile),width=8.5,height=4.1)
    print(pl);dev.off()
  }
  make_effect('mean','Fig5G_mean_bias_review.pdf','Mean clonal bias,')
  make_effect('sd','Fig5H_variance_review.pdf','Coefficient of variation (archived),')
  # Empirical bootstrap P-values may be exactly zero. For display only,
  # place such points at 1e-5; the exported P-values remain unmodified.
  p <- tbl
  p$xplot <- ifelse(is.finite(p$mean.pv),pmax(p$mean.pv,1e-5),NA_real_)
  p$yplot <- ifelse(is.finite(p$sd.pv),pmax(p$sd.pv,1e-5),NA_real_)
  pl <- ggplot(p,aes(xplot,yplot,color=state_col)) +
    geom_vline(xintercept=0.01,linetype='dashed',color='grey50') +
    geom_hline(yintercept=0.01,linetype='dashed',color='grey50') +
    geom_point(size=2.5) +
    geom_text(aes(label=state),size=3,vjust=-0.8,show.legend=FALSE,check_overlap=TRUE) +
    facet_wrap(~contrast,nrow=1,scales='free') +
    scale_x_log10() + scale_y_log10() +
    scale_color_manual(values=state_cols) +
    labs(x='Mean-bias difference, empirical P value',
         y='Variance difference, empirical P value',color='Cell state') +
    theme_bw(base_size=12)
  pdf(file.path(out,'FigS8E_bootstrap_pvalues_review.pdf'),width=8.5,height=4.1)
  print(pl);dev.off()
  message('Review figures created. Fig5F requires a separate distribution-selection review.')
})
