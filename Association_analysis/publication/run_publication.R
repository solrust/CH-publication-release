# Canonical orchestration: sourcing defines a function; it does NOT run analysis.
# The retained analytical helpers remain required runtime dependencies.
run_association_publication <- function(root=getwd(), stages=NULL, dry_run=TRUE,
                              work_root=NULL, allow_refit=FALSE) {
  root <- normalizePath(root,mustWork=TRUE)
  h <- new.env(parent=globalenv())
  sys.source(file.path(root,'publication/tools/early_RNA_launcher.R'),envir=h)
  steps <- list(
    `5C` = c("export_Fig5C.R", "plot_Fig5C.R"),
    `4FH_S7F` = c("export_Fig4FH_FigS7F.R", "audit_Fig4H_accepted_pathways.R", "curate_Fig4H_published_pathways.R", "plot_Fig4FH_FigS7F.R"),
    `4EG` = c("export_Fig4EG.R", "plot_Fig4EG.R"),
    `4I` = c("export_Fig4I.R", "reconstruct_Fig4I_GO_categories.R", "plot_Fig4I.R"),
    `4JK` = c("export_Fig4JK.R", "plot_Fig4JK.R", "plot_Fig4K_with_AUC.R"),
    `4L` = c("export_Fig4L.R", "plot_Fig4L.R"),
    `4M` = c("export_Fig4M.R", "plot_Fig4M.R"),
    `5B` = c("export_Fig5B.R", "plot_Fig5B.R"),
    `S8ABC` = c("export_FigS8AB.R", "plot_FigS8AB.R", "plot_FigS8C.R"),
    `5J_plot` = c("plot_Fig5J.R"),
    `5J_refit` = c("export_Fig5J.R", "plot_Fig5J.R")
  )
  h$ch_early_run(project="Association_analysis", steps=steps,
    default_stages=c("5C", "4FH_S7F", "4EG", "4I", "4JK", "4L", "4M", "5B", "S8ABC", "5J_plot"), root=root, stages=stages, dry_run=dry_run,
    work_root=work_root, allow_refit=allow_refit)
}
