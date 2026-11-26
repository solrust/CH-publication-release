# Shared orchestration only. Never source a scientific helper during planning.
ch_early_run <- function(project, steps, default_stages, root, stages = NULL,
                         dry_run = TRUE, work_root = NULL, allow_refit = FALSE) {
  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, project, 'publication'))) stop('Wrong repository root.')
  if (is.null(stages)) stages <- default_stages
  if (!length(stages) || anyNA(stages) || any(!stages %in% names(steps)))
    stop('Unknown stage. Available: ', paste(names(steps), collapse=', '))
  stages <- unique(stages)
  plan <- do.call(rbind, lapply(stages, function(st) data.frame(
    project=project, stage=st, script=file.path(project,'publication',steps[[st]]),
    stringsAsFactors=FALSE)))
  plan$present <- file.exists(file.path(root, plan$script))
  plan$execution_type <- ifelse(plan$stage=='5J_plot','SAVED_SUMMARY_RENDER_ONLY',
    ifelse(plan$stage=='5J_refit','EXPLICIT_MODEL_REFIT','HISTORICAL_HELPERS_UNCHANGED'))
  print(plan, row.names=FALSE)
  if (!all(plan$present)) stop('Required helper missing. No analysis started.')
  if (isTRUE(dry_run)) {
    message('PLAN ONLY: no datasets loaded; no analyses run; no outputs changed.')
    if (project=='RNA_analysis') message('S5C and S6D: deferred by user; not scheduled.')
    return(invisible(plan))
  }
  if ('5J_refit' %in% stages && !isTRUE(allow_refit))
    stop('5J_refit requires allow_refit=TRUE and sufficient memory. Default uses saved summary.')
  if (is.null(work_root) || length(work_root)!=1L || !nzchar(work_root))
    stop('Execution requires a NEW work_root; publication baselines cannot be execution targets.')
  if (file.exists(work_root) || dir.exists(work_root)) stop('work_root must not already exist.')
  # work_root is a repository-shaped test workspace, not an output prefix.
  # Linux symlinks expose original read inputs. This is NOT an OS security sandbox.
  if (!dir.create(work_root, recursive=TRUE)) stop('Cannot create work_root.')
  work_root <- normalizePath(work_root, mustWork=TRUE)
  if (identical(root, work_root)) stop('Cannot execute into canonical root.')
  copy_one <- function(rel, required=TRUE) {
    from <- file.path(root,rel); to <- file.path(work_root,rel)
    if (!file.exists(from)) {
      if (required) stop('Missing input: ',rel)
      return(invisible(FALSE))
    }
    dir.create(dirname(to), recursive=TRUE, showWarnings=FALSE)
    if (!file.copy(from,to,overwrite=FALSE,copy.mode=TRUE)) stop('Cannot copy: ',rel)
    invisible(TRUE)
  }
  input_projects <- c('RNA_analysis','RNA_clone_analysis','Association_analysis','DNA_analysis')
  for (p in input_projects) {
    dir.create(file.path(work_root,p),recursive=TRUE,showWarnings=FALSE)
    for (d in c('data','scripts','R')) {
      from <- file.path(root,p,d)
      if (dir.exists(from) && !file.symlink(from,file.path(work_root,p,d)))
        stop('Cannot link input directory: ',from)
    }
    src <- file.path(root,p,'publication')
    dst <- file.path(work_root,p,'publication')
    dir.create(dst,recursive=TRUE,showWarnings=FALSE)
    # Copy immediate code/resources only, never historical output directories.
    if (dir.exists(src)) for (f in list.files(src,full.names=TRUE)) {
      if (!dir.exists(f) && grepl('\\.(R|r|tsv|csv|md|txt|json|yaml|yml)$', f))
        copy_one(file.path(p,'publication',basename(f)))
    }
  }
  seed <- character()
  # Small figure definitions travel with the code into the run workspace.
  if (project=='RNA_analysis' && '3DE' %in% stages)
    seed <- c(seed,'RNA_analysis/publication/resources/Fig3DE_published_pathway_ids.tsv')
  if (project=='RNA_analysis' && '3F' %in% stages)
    seed <- c(seed,'RNA_analysis/publication/resources/Fig3F_functional_categories.tsv')
  if (project=='RNA_clone_analysis' && '4B' %in% stages)
    seed <- c(seed,'RNA_analysis/publication/source_data/Fig3A_umap_cells.tsv.gz')
  if (project=='Association_analysis' && '4FH_S7F' %in% stages)
    seed <- c(seed,'Association_analysis/publication/resources/Fig4H_accepted_pathway_ids.tsv')
  if (project=='Association_analysis' && 'S8ABC' %in% stages && !'4M' %in% stages)
    seed <- c(seed,'Association_analysis/publication/source_data/Fig4M_HSC_module_projection_cells.tsv.gz')
  if (project=='Association_analysis' && '5J_plot' %in% stages)
    seed <- c(seed,'Association_analysis/publication/source_data/Fig5J_heritability_summary.tsv')
  for (f in seed) copy_one(f)
  logdir <- file.path(work_root,'execution_logs')
  dir.create(logdir)
  if (length(seed)) utils::write.table(data.frame(
    input=seed, provenance='COPIED_EXISTING_INPUT_NOT_RECOMPUTED',
    md5=unname(tools::md5sum(file.path(root,seed)))),
    file.path(logdir,'seeded_inputs.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
  rscript <- file.path(R.home('bin'),'Rscript')
  if (!file.exists(rscript)) stop('Rscript unavailable in current R installation.')
  quote_r <- function(x) paste(capture.output(dput(x)),collapse='\n')
  reports <- list()
  for (st in stages) {
    script <- file.path(logdir,paste0(st,'.R'))
    log <- file.path(logdir,paste0(st,'.log'))
    paths <- plan$script[plan$stage==st]
    # A new --vanilla process per stage avoids restored RStudio objects/namespaces.
    boot <- c(
      paste0('setwd(',quote_r(work_root),')'),
      'stopifnot(getRversion() == package_version("4.2.3"))',
      '.libPaths(unique(c(Sys.getenv("CH_R_LIB", unset = "/packages"), .libPaths())))',
      'library(Matrix, lib.loc=Sys.getenv("CH_R_LIB", unset = "/packages"))',
      'stopifnot(getNamespaceVersion("Matrix") == package_version("1.6-4"),',
      ' normalizePath(getNamespaceInfo("Matrix", "path")) == normalizePath(file.path(Sys.getenv("CH_R_LIB", unset = "/packages"), "Matrix")))',
      paste0('Sys.setenv(CH_DATA_ROOT=',quote_r(work_root),')'),
      'options(CH.s5d.reuse_loaded=FALSE, CH.fig5j.splits=300L)',
      'print(sessionInfo())',
      unlist(lapply(paths,function(f) c(
        paste0('message("RUN: ", ',quote_r(f),')'),
        paste0('source(',quote_r(f),', local=FALSE, chdir=FALSE)')))),
      'print(sessionInfo())')
    writeLines(boot,script)
    message('Isolated stage: ',st,'; log: ',log)
    status <- system2(rscript,c('--vanilla',shQuote(script)),stdout=log,stderr=log)
    reports[[st]] <- data.frame(stage=st,exit_status=status,
      status=if(status==0L) 'COMPLETED_NOT_NUMERICALLY_COMPARED' else 'FAILED',
      log=log,stringsAsFactors=FALSE)
    utils::write.table(do.call(rbind,reports),file.path(logdir,'execution_report.tsv'),
                      sep='\t',quote=TRUE,row.names=FALSE)
    if (status!=0L) stop('Stage failed: ',st,'. Baselines unchanged; inspect ',log)
  }
  message('Isolated execution finished. Compare outputs with baselines before acceptance.')
  invisible(list(work_root=work_root,plan=plan,execution=do.call(rbind,reports)))
}
