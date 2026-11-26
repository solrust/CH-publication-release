# Canonical orchestration: sourcing defines a function; it does NOT run analysis.
# The retained analytical helpers remain required runtime dependencies.
run_rna_publication <- function(root=getwd(), stages=NULL, dry_run=TRUE,
                              work_root=NULL, allow_refit=FALSE) {
  root <- normalizePath(root,mustWork=TRUE)
  h <- new.env(parent=globalenv())
  sys.source(file.path(root,'publication/tools/early_RNA_launcher.R'),envir=h)
  steps <- list(
    `3A_3B_S6AB` = c("export_Fig3A_FigS6AB.R", "export_Fig3B.R", "plot_Fig3A_Fig3B_FigS6AB.R"),
    `3DE` = c("export_Fig3DE.R", "plot_Fig3DE.R", "check_Fig3DE_published_pathways.R"),
    `3F` = c("export_Fig3F_only.R", "plot_Fig3F_only.R"),
    `S6C` = c("export_FigS6C.R", "plot_FigS6C.R"),
    `S5A` = c("export_FigS5A.R", "plot_FigS5A.R"),
    `S5B` = c("export_FigS5B.R", "plot_FigS5B.R"),
    `S5D` = c("export_FigS5D_only.R", "plot_FigS5D_only.R")
  )
  h$ch_early_run(project="RNA_analysis", steps=steps,
    default_stages=c("3A_3B_S6AB", "3DE", "3F", "S6C", "S5A", "S5B", "S5D"), root=root, stages=stages, dry_run=dry_run,
    work_root=work_root, allow_refit=allow_refit)
}
