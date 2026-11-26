# Supplementary Fig. 5D, stand-alone/resumable score export.
# Run from CH_publication_release root. Does NOT rerun S5C, Pagoda2 or DESeq2.
# By default reads the historical scrna.RData checkpoint.
# To reuse the loaded object `e` from a just-failed export_FigS5CD.R run:
#   options(CH.s5d.reuse_loaded = TRUE)
#   source('RNA_analysis/publication/export_FigS5D_only.R')
#   options(CH.s5d.reuse_loaded = NULL)

if (!requireNamespace('Matrix', quietly = TRUE)) stop('Install Matrix first.')
source('RNA_analysis/publication/project_gset_publication.R', local = TRUE)

root <- Sys.getenv('CH_DATA_ROOT', unset = '.')
data_dir <- file.path(root, 'RNA_analysis', 'data')
checkpoint <- file.path(data_dir, 'processed', 'scrna.RData')
program_file <- file.path(data_dir, 'genesets', 'processed', 'fate_markers.rds')
stopifnot(file.exists(program_file))

outdir <- 'RNA_analysis/publication/source_data'
valdir <- 'RNA_analysis/publication/validation'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(valdir, showWarnings = FALSE, recursive = TRUE)

if (isTRUE(getOption('CH.s5d.reuse_loaded', FALSE))) {
  if (!exists('e', envir = .GlobalEnv, inherits = FALSE)) {
    stop('No `e` environment in current session; disable reuse_loaded.')
  }
  src <- get('e', envir = .GlobalEnv, inherits = FALSE)
  message('Reusing loaded `e` environment from the previous export attempt.')
} else {
  stopifnot(file.exists(checkpoint))
  message('Loading ', checkpoint)
  src <- new.env(parent = globalenv())
  load(checkpoint, envir = src)
}
if (!is.environment(src) ||
    !all(vapply(c('p2counts', 'emb'), exists, logical(1),
                envir = src, inherits = FALSE))) {
  stop('Saved environment must contain p2counts and emb.')
}
mat <- src$p2counts
emb <- src$emb
stopifnot(length(dim(mat)) == 2L, ncol(emb) >= 2L,
          !is.null(rownames(mat)), !is.null(colnames(mat)),
          !is.null(rownames(emb)))
programs <- readRDS(program_file)

write_gz <- function(df, name) {
  con <- gzfile(file.path(outdir, name), open = 'wt')
  on.exit(close(con))
  utils::write.table(df, con, sep = '\t', quote = FALSE,
                     row.names = FALSE, na = 'NA')
}

# Do not include program 13: reproduces original 4_sc_analysis.R selection.
indices <- setdiff(seq_along(programs), 13L)
program_info <- vector('list', length(indices))
summary_rows <- vector('list', length(indices))
for (ii in seq_along(indices)) {
  k <- indices[ii]
  name <- names(programs)[k]
  if (is.null(name) || is.na(name) || !nzchar(name)) name <- paste0('program_', k)
  message('S5D program ', k, ': ', name)
  result <- tryCatch(project_gset_publication(programs[[k]], mat),
                     error = function(err) stop('Program ', k, ' (', name,
                                                ') failed: ', conditionMessage(err)))
  if (is.null(result)) stop('Program ', k, ' has no expressed genes.')
  vals <- result$av.fit
  if (length(vals) != nrow(mat)) stop('Incorrect score length for ', name)
  ids <- names(vals)
  if (is.null(ids) || !identical(ids, rownames(mat))) ids <- rownames(mat)
  ix <- match(ids, rownames(emb))
  keep <- !is.na(ix) & is.finite(vals)
  if (!any(keep)) stop('No cells match UMAP coordinates for ', name)
  y <- as.numeric(vals[keep])
  display_score <- sign(y) * (abs(y) * 1e8 + 1)^0.3
  df <- data.frame(program_index = k, program = name, cell_id = ids[keep],
                   umap_1 = emb[ix[keep], 1L], umap_2 = emb[ix[keep], 2L],
                   fitted_score = y, display_score = display_score,
                   stringsAsFactors = FALSE)
  fname <- sprintf('FigS5D_program_%02d_scores.tsv.gz', k)
  write_gz(df, fname)
  top <- names(sort(abs(result$gscore), decreasing = TRUE))
  program_info[[ii]] <- data.frame(program_index = k, program = name,
     rank = seq_along(top), gene = top,
     gene_loading = as.numeric(result$gscore[top]), stringsAsFactors = FALSE)
  summary_rows[[ii]] <- data.frame(program_index = k, program = name,
     genes_expressed = result$n_matched_expressed_genes,
     n_cells_exported = nrow(df), output_file = fname,
     stringsAsFactors = FALSE)
  rm(result, vals, df, display_score)
  gc(verbose = FALSE)
}
write.table(do.call(rbind, program_info),
   file.path(outdir, 'FigS5D_program_gene_loadings.tsv'),
   sep = '\t', row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, summary_rows),
   file.path(valdir, 'FigS5D_export_summary.tsv'),
   sep = '\t', row.names = FALSE, quote = FALSE)
cat('S5D: exported ', length(summary_rows),
    ' gene programs; S5C data were not changed.\n', sep = '')
