# Source once, from the repository root, in a fresh R 4.2.3 session.
# Set CH_R_LIB to an absolute directory before starting R; /packages is the IRIS default.
ch_configure_environment <- function(root = getwd(), create_layout = TRUE) {
  root <- normalizePath(root, mustWork = TRUE)
  if (!file.exists(file.path(root, 'publication/tools/early_RNA_launcher.R')))
    stop('Use CH_publication_release as the working directory')
  lib <- Sys.getenv('CH_R_LIB', unset = '/packages')
  if (!dir.exists(lib)) stop('R library not found; set CH_R_LIB: ', lib)
  lib <- normalizePath(lib, mustWork = TRUE)
  Sys.setenv(CH_R_LIB = lib)
  stopifnot(getRversion() == package_version('4.2.3'))
  .libPaths(unique(c(lib, .libPaths())))
  matrix_path <- normalizePath(file.path(lib, 'Matrix'), mustWork = TRUE)
  if (!'Matrix' %in% loadedNamespaces()) library(Matrix, lib.loc = lib)
  if (getNamespaceVersion('Matrix') != package_version('1.6-4') ||
      normalizePath(getNamespaceInfo('Matrix', 'path'), mustWork = TRUE) != matrix_path)
    stop('Restart R and load Matrix 1.6-4 from CH_R_LIB before using archived data.')
  if (create_layout) {
    projects <- c('RNA_analysis','RNA_clone_analysis','Association_analysis','multiome',
                  'human_ch','DNA_analysis','clones','TF_overexpression','MYC_analysis',
                  'inflammatory_HSCs','indels_calling')
    # Git does not preserve empty historical figure directories. Some runners resolve
    # these to enforce their no-overwrite boundary. No historical figures are created.
    for (p in projects) if (dir.exists(file.path(root,p)))
      dir.create(file.path(root,p,'figures'), recursive=TRUE, showWarnings=FALSE)
  }
  message('R ', getRversion(), '; Matrix ', getNamespaceVersion('Matrix'), '; library ', lib)
  invisible(lib)
}
