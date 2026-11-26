# Canonical orchestration: sourcing defines a function; it does NOT run analysis.
# The retained analytical helpers remain required runtime dependencies.
run_rna_clone_publication <- function(root=getwd(), stages=NULL, dry_run=TRUE,
                              work_root=NULL, allow_refit=FALSE) {
  root <- normalizePath(root,mustWork=TRUE)
  h <- new.env(parent=globalenv())
  sys.source(file.path(root,'publication/tools/early_RNA_launcher.R'),envir=h)
  steps <- list(
    `4C` = c("export_plot_Fig4C.R"),
    `4B` = c("export_Fig4B.R", "plot_Fig4B.R"),
    `S7BCD` = c("export_FigS7BCD.R", "plot_FigS7BCD.R"),
    `5E_S8D` = c("export_Fig5E_FigS8D.R", "plot_Fig5E_FigS8D.R"),
    `5GH_S8E` = c("export_Fig5GH_S8E.R", "plot_Fig5GH_S8E.R")
  )
  h$ch_early_run(project="RNA_clone_analysis", steps=steps,
    default_stages=c("4C", "4B", "S7BCD", "5E_S8D", "5GH_S8E"), root=root, stages=stages, dry_run=dry_run,
    work_root=work_root, allow_refit=allow_refit)
}
