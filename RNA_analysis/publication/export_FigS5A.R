# Supplementary Fig. 5A: hashtag count pairs underlying six displayed plots.
# Run from CH_publication_release root. Does not run demultiplexing or modify inputs.

stopifnot(file.exists('RNA_analysis/scripts/1_hashtag_analysis.R'))
if (!requireNamespace('Matrix', quietly = TRUE)) stop('Install Matrix before running this script.')

root <- Sys.getenv('CH_DATA_ROOT', unset = '.')
data_dir <- file.path(root, 'RNA_analysis', 'data')
outdir <- 'RNA_analysis/publication/source_data'
valdir <- 'RNA_analysis/publication/validation'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(valdir, recursive = TRUE, showWarnings = FALSE)

plot_runs <- c(1L, 3L, 4L, 5L, 6L, 8L)  # original: omit runs 2 and 7
hash_ids <- list(c(1L, 2L), c(3L, 3L), c(1L, 3L), c(1L, 2L),
                 c(1L, 2L), c(3L, 1L), c(3L, 3L), c(1L, 2L))
input_paths <- vapply(plot_runs, function(i) {
  exp_id <- if (i <= 4L) '123' else '127'
  local_run <- if (i <= 4L) i else i - 4L
  file.path(data_dir, paste0('raw/exp', exp_id), 'hashtags',
            paste0('RNA_', local_run), 'Results', 'umi_count')
}, character(1))
missing_inputs <- unlist(lapply(input_paths, function(d) {
  f <- file.path(d, c('matrix.mtx.gz', 'barcodes.tsv.gz', 'features.tsv.gz'))
  f[!file.exists(f)]
}), use.names = FALSE)
if (length(missing_inputs)) {
  cat('Required files missing for S5A:\n', paste(missing_inputs, collapse='\n'), '\n')
  stop('S5A not exported. Locate archived hashtag count matrices or mark S5A as deferred.')
}

write_tsv_gz <- function(df, filename) {
  con <- gzfile(filename, open='wt')
  on.exit(close(con))
  write.table(df, file=con, sep='\t', quote=FALSE, row.names=FALSE, na='NA')
}

panels <- vector('list', length(plot_runs))
summaries <- vector('list', length(plot_runs))
for (j in seq_along(plot_runs)) {
  i <- plot_runs[j]
  experiment <- if (i <= 4L) '123' else '127'
  local_run <- if (i <= 4L) i else i - 4L
  cutoff <- if (i <= 4L) 150L else 90L
  d <- input_paths[j]
  mcon <- gzfile(file.path(d, 'matrix.mtx.gz'), open='rt')
  m <- tryCatch(Matrix::readMM(mcon), finally=close(mcon))
  barcode <- read.delim(gzfile(file.path(d, 'barcodes.tsv.gz')), header=FALSE,
                        stringsAsFactors=FALSE)[[1L]]
  features <- read.delim(gzfile(file.path(d, 'features.tsv.gz')), header=FALSE,
                         stringsAsFactors=FALSE)[[1L]]
  if (length(barcode) != ncol(m) || length(features) != nrow(m))
    stop('Mismatch between count matrix, features and barcodes in ', d)
  ids <- hash_ids[[i]]
  x <- as.numeric(m[ids[1], ])
  y <- as.numeric(m[ids[2], ])
  panels[[j]] <- data.frame(
    panel_order=j, experiment=experiment, RNA_run=local_run,
    original_run_index=i, barcode=barcode,
    hashtag_x=features[ids[1]], hashtag_y=features[ids[2]],
    x_count=x, y_count=y, threshold=cutoff,
    stringsAsFactors=FALSE
  )
  summaries[[j]] <- data.frame(
    panel_order=j, experiment=experiment, RNA_run=local_run,
    n_barcodes=length(barcode), threshold=cutoff,
    n_both_above=sum(x>cutoff & y>cutoff),
    n_x_only=sum(x>cutoff & y<cutoff),
    n_y_only=sum(x<cutoff & y>cutoff)
  )
  cat('S5A:', experiment, 'RNA', local_run, ':', length(barcode), 'barcodes\n')
  rm(m, barcode, features, x, y)
  gc(verbose=FALSE)
}
all_points <- do.call(rbind, panels)
write_tsv_gz(all_points, file.path(outdir, 'FigS5A_hashtag_count_pairs.tsv.gz'))
write.table(do.call(rbind, summaries),
            file.path(outdir, 'FigS5A_hashtag_panel_summary.tsv'),
            sep='\t', quote=FALSE, row.names=FALSE)
cat('S5A exported six panels. Source data: ', nrow(all_points), 'barcode pairs.\n')
