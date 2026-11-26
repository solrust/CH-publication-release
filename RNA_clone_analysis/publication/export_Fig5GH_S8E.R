# Publication exports for Fig 5G/H and Supplementary Fig 8E.
# Run from the CH_publication_release root. Reads saved 100,000-bootstrap
# Poisson-lognormal results, not raw LARRY or original model fitting.
local({
  pub <- 'RNA_clone_analysis/publication'
  out <- file.path(pub,'source_data')
  val <- file.path(pub,'validation')
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  dir.create(val, recursive=TRUE, showWarnings=FALSE)
  fields <- c('state','mean.wt','mean.mut','mean.wt.down','mean.wt.up',
              'mean.mut.down','mean.mut.up', 'sd.wt','sd.mut',
              'sd.wt.down','sd.wt.up','sd.mut.down','sd.mut.up',
              'mean.pv','sd.pv')
  result <- list(); qa <- list()
  for (gen in c('TET2','DNMT3A')) {
    input <- file.path('RNA_clone_analysis','data','processed',
                       paste0('clone_patterns_',gen,'.RData'))
    if (!file.exists(input)) stop('Required saved checkpoint absent: ',input)
    e <- new.env(parent=baseenv())
    load(input, envir=e)
    if (!exists('dat',envir=e,inherits=FALSE)) stop('Missing dat in ',input)
    d <- e$dat
    missing <- setdiff(fields,colnames(d))
    if (length(missing)) stop(input,' lacks columns: ',paste(missing,collapse=', '))
    if(nrow(d)!=10L) stop('Expected 10 cell states in ', input, ', got ',nrow(d))
    if(!setequal(as.integer(as.character(d$state)),1:10))
      stop('Unexpected cell state IDs in ',input)
    if(exists('geno',envir=e,inherits=FALSE) && !identical(as.character(e$geno),gen))
      stop('Saved geno does not match filename: ',input)
    d <- d[,fields,drop=FALSE]
    d$contrast <- gen
    d$state <- as.integer(as.character(d$state))
    # Keep values as archived, including exact zeros in bootstrap P-values.
    result[[gen]] <- d[,c('contrast',fields),drop=FALSE]
    n.wt <- if(exists('bootsrap.wt0',envir=e,inherits=FALSE)) nrow(e$bootsrap.wt0) else NA_integer_
    n.mut <- if(exists('bootsrap.mut0',envir=e,inherits=FALSE)) nrow(e$bootsrap.mut0) else NA_integer_
    clone.cov <- if(exists('clone.cov',envir=e,inherits=FALSE)) e$clone.cov else NA_integer_
    qa[[gen]] <- data.frame(contrast=gen, states=nrow(d),
                            bootstrap_WT=n.wt, bootstrap_mutant=n.mut,
                            min_clone_coverage=clone.cov,
                            n_zero_mean_p=sum(d$mean.pv==0,na.rm=TRUE),
                            n_zero_variance_p=sum(d$sd.pv==0,na.rm=TRUE))
    message(gen,': ',nrow(d),' states; WT bootstraps=',n.wt,
            ', mutant bootstraps=',n.mut)
    rm(e)
    gc(verbose=FALSE)
  }
  dat <- do.call(rbind,result)
  rownames(dat) <- NULL
  if(anyDuplicated(paste(dat$contrast,dat$state))) stop('Repeated contrast-state pairs')
  write.table(dat,file.path(out,'Fig5GH_FigS8E_state_estimates.tsv'),
              sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  report <- do.call(rbind,qa)
  write.table(report,file.path(val,'Fig5GH_FigS8E_export_summary.tsv'),
              sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  print(report,row.names=FALSE)
  message('Source Data exported. Replot and compare against accepted panels before validation.')
})
