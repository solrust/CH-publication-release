# clones/publication/run_publication.R
# Figures 5L, S8H and S8I from historical Starcode count files.
# Execute from CH_publication_release; never edits historical scripts or data.

clones_pub_run <- function(root = getwd(), output_root = file.path(root, 'clones', 'publication'),
                           validate = TRUE) {
  stopifnot(dir.exists(root), dir.exists(file.path(root, 'clones', 'data', 'GS175')))
  if (!requireNamespace('ggplot2', quietly = TRUE)) stop('Install ggplot2 in the established R environment')
  if (!requireNamespace('confintr', quietly = TRUE)) stop('Install conflintr in the established R environment')
  root <- normalizePath(root)
  datadir <- file.path(root, 'clones', 'data', 'GS175')
  mdfile <- file.path(datadir, 'GS_175_larry in vitro 3.csv')
  md <- read.csv(mdfile, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(md) != 112L || !identical(as.integer(md$sample), seq_len(112L)))
    stop('Metadata must have sample numbers 1..112 in ascending order; inspect source before continuing')
  parts <- regexec('^(AAVS1|TET2)_([12])([AB])$', md$genotype)
  extracted <- regmatches(md$genotype, parts)
  if (any(lengths(extracted) != 4L)) stop('Unexpected genotype/replicate/technical split format')
  md$geno <- vapply(extracted, `[`, '', 2L)
  md$replicate <- vapply(extracted, `[`, '', 3L)
  md$technical <- vapply(extracted, `[`, '', 4L)
  stopifnot(identical(unique(md$geno), c('AAVS1', 'TET2')))
  files <- list.files(file.path(datadir, 'starcode'), full.names = TRUE, recursive = TRUE,
                      pattern = 'barcode_counts_s\\.txt$')
  ids <- suppressWarnings(as.integer(sub('.*175-([0-9]+)_.*', '\\1', basename(files))))
  if (length(files) != 112L || anyNA(ids) || anyDuplicated(ids) ||
      !setequal(ids, seq_len(112L))) stop('Expected exactly one Starcode count file for samples 1..112')
  files <- files[order(ids)]
  cat('clones: reading', length(files), 'historical Starcode tables\n')
  fls <- lapply(files, function(f) {
    x <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(x) < 2L) stop('Expected at least two Starcode columns in ', f)
    x <- data.frame(barcode = as.character(x[[1L]]), reads = as.numeric(x[[2L]]),
                    stringsAsFactors = FALSE)
    if (anyNA(x) || anyDuplicated(x[[1L]]) || any(x[[2L]] < 0))
      stop('Invalid or duplicate barcode counts in ', f)
    x
  })
  clone_ids <- unique(unlist(lapply(fls, `[[`, 1L), use.names = FALSE))
  counts <- matrix(0, nrow = length(clone_ids), ncol = 112L,
                   dimnames = list(clone_ids, paste0('sample_', seq_len(112L))))
  for (j in seq_along(fls)) counts[fls[[j]][[1L]], j] <- fls[[j]][[2L]]
  # Historical filter is >=20 total reads over all samples, before computing any fraction.
  counts <- counts[rowSums(counts) >= 20, , drop = FALSE]
  cat('clones: retained', nrow(counts), 'barcodes at >=20 total reads\n')

  # Reproduce exactly the historical sample selection and per-split cut_count=4.
  # Each condition must have ONE A and ONE B count column for target and bulk.
  calc_one <- function(day, repl, geno, pop) {
    getcol <- function(population, tech) {
      sel <- which(md$timepoint == day & md$replicate == repl &
                   md$geno == geno & md$population_short == population &
                   md$technical == tech)
      if (length(sel) != 1L) stop('Expected one sample, found ', length(sel), ': ',
                                  paste(day, repl, geno, population, tech, collapse = '/'))
      sel
    }
    ia <- getcol(pop, 'A'); ib <- getcol(pop, 'B')
    ba <- getcol('bulk', 'A'); bb <- getcol('bulk', 'B')
    target_a <- counts[, ia]; target_b <- counts[, ib]
    bulk_a <- counts[, ba]; bulk_b <- counts[, bb]
    keep <- (target_a + bulk_a >= 4) & (target_b + bulk_b >= 4)
    if (sum(keep) < 4L) stop('Fewer than four eligible clones for ', paste(day, repl, geno, pop))
    v1 <- target_a[keep] / (target_a[keep] + bulk_a[keep])
    v2 <- target_b[keep] / (target_b[keep] + bulk_b[keep])
    if (sd(v1) == 0 || sd(v2) == 0) stop('Undefined Pearson correlation for ', paste(day, repl, geno, pop))
    ci <- confintr::ci_cor(v1, v2, probs = c(0.05, 0.95), method = 'pearson', type = 'normal')
    r <- as.numeric(ci$estimate)[1L]; interval <- as.numeric(ci$interval)
    if (length(interval) != 2L || any(!is.finite(c(r, interval))))
      stop('Nonfinite/missing correlation/CI for ', paste(day, repl, geno, pop))
    if (validate && !isTRUE(all.equal(r, unname(stats::cor(v1, v2)), tolerance = 1e-10)))
      stop('Independent Pearson check failed for ', paste(day, repl, geno, pop))
    info <- data.frame(day = day, replicate = as.integer(repl), population = pop,
                       genotype = geno, n_clones = sum(keep), pearson_r = r,
                       ci90_lower = interval[1L], ci90_upper = interval[2L],
                       total_read_min = 20L, split_combined_read_min = 4L,
                       stringsAsFactors = FALSE)
    idx <- which(keep)
    pts <- data.frame(day = day, replicate = as.integer(repl), population = pop,
                      genotype = geno, clone_id = rownames(counts)[idx],
                      sample_A_target = md$sample[ia], sample_A_bulk = md$sample[ba],
                      sample_B_target = md$sample[ib], sample_B_bulk = md$sample[bb],
                      split_A_target_reads = target_a[idx], split_A_bulk_reads = bulk_a[idx],
                      split_B_target_reads = target_b[idx], split_B_bulk_reads = bulk_b[idx],
                      split_A_fraction = v1, split_B_fraction = v2,
                      stringsAsFactors = FALSE)
    list(summary = info, points = pts)
  }
  pops <- c('bulk CD34neg', 'early_My', 'early_Ery', 'early_undefined')
  combos <- expand.grid(day = c('day 5', 'day 10'), replicate = c('1', '2'),
                       population = pops, genotype = c('AAVS1', 'TET2'),
                       stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(combos)), function(i) {
    x <- combos[i, ]; calc_one(x$day, x$replicate, x$genotype, x$population)
  })
  summ <- do.call(rbind, lapply(res, `[[`, 'summary'))
  points <- do.call(rbind, lapply(res, `[[`, 'points'))
  stopifnot(nrow(summ) == 32L, all(summ$n_clones > 3),
            sum(summ$n_clones) == nrow(points))
  # 5L day5, repl1+2; S8I day10, repl1; S8H day5 repl1 erythroid scatter.
  plotted <- c('early_Ery', 'early_My', 'early_undefined')
  f5l <- subset(summ, day == 'day 5' & population %in% plotted)
  fs8i <- subset(summ, day == 'day 10' & replicate == 1L & population %in% plotted)
  fs8h <- subset(points, day == 'day 5' & replicate == 1L & population == 'early_Ery')
  stopifnot(nrow(f5l) == 12L, nrow(fs8i) == 6L,
            all(table(fs8h$genotype) > 3L))
  art_r <- sapply(c('AAVS1', 'TET2'), function(g) {
    z <- subset(summ, day == 'day 5' & replicate == 1L &
                population == 'early_Ery' & genotype == g)
    z$pearson_r
  })
  cat('S8H: Pearson r AAVS1=', sprintf('%.4f', art_r['AAVS1']),
      '; TET2=', sprintf('%.4f', art_r['TET2']),
      '; manuscript labels approximately 0.30 / 0.43\n', sep = '')
  if (abs(art_r['AAVS1'] - 0.30) > 0.025 || abs(art_r['TET2'] - 0.43) > 0.025)
    warning('S8H correlations differ from manuscript rounded labels; mark historical visual/numeric REVIEW')

  if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
  srcdir <- file.path(output_root, 'source_data')
  figdir <- file.path(output_root, 'figures')
  src_names <- c('Fig5L_lineage_heritability.tsv',
                 'FigS8H_erythroid_split_clones.tsv.gz',
                 'FigS8I_lineage_heritability.tsv',
                 'clones_all_heritability_correlations.tsv',
                 'clones_all_split_clone_fractions.tsv.gz')
  fig_names <- c('Fig5L_lineage_heritability_REVIEW.pdf',
                 'FigS8H_erythroid_split_scatter_REVIEW.pdf',
                 'FigS8I_lineage_heritability_REVIEW.pdf')
  target_files <- c(file.path(srcdir, src_names), file.path(figdir, fig_names))
  if (any(file.exists(target_files)))
    stop('Publication outputs already exist; refusing overwrite: ',
         paste(target_files[file.exists(target_files)], collapse = ', '))
  stage <- tempfile('.clones_stage_', tmpdir = output_root)
  dir.create(stage); on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  dir.create(file.path(stage, 'source_data'))
  dir.create(file.path(stage, 'figures'))
  put <- function(x, name) {
    file <- file.path(stage, 'source_data', name)
    if (grepl('\\.gz$', name)) {
      cc <- gzfile(file, open = 'wt'); on.exit(close(cc), add = TRUE)
      write.table(x, file = cc, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
    } else {
      write.table(x, file = file, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
    }
  }
  put(f5l, src_names[1L]); put(fs8h, src_names[2L]); put(fs8i, src_names[3L])
  put(summ, src_names[4L]); put(points, src_names[5L])

  # Original summary plots: three early lineages, AAVS1/TET2, 90% CIs.
  colors <- c(AAVS1 = 'grey50', TET2 = '#22A33A')
  labels <- c(early_Ery = 'early\nerythroid', early_My = 'early\nmyeloid',
              early_undefined = 'early\nuncommitted')
  plot_bars <- function(df, facet = FALSE) {
    df$population <- factor(df$population, levels = names(labels), labels = unname(labels))
    df$genotype <- factor(df$genotype, levels = c('AAVS1', 'TET2'))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = population, y = pearson_r, fill = genotype)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = .7), width = .64) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ci90_lower, ymax = ci90_upper),
                            position = ggplot2::position_dodge(width = .7), width = .17) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(x = NULL, y = 'Clonal correlation across replicates', fill = 'Genotype') +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(legend.position = 'top')
    if (facet) p <- p + ggplot2::facet_wrap(~replicate, nrow = 1,
                                          labeller = ggplot2::labeller(replicate = function(x) paste('replicate', x)))
    p
  }
  pdf(file.path(stage, 'figures', fig_names[1L]), width = 9.0, height = 4.0)
  print(plot_bars(f5l, facet = TRUE)); dev.off()
  # Publication S8H: early_Ery day 5 replicate 1, panels AAVS1/TET2.
  fs8h$genotype <- factor(fs8h$genotype, levels = c('AAVS1', 'TET2'))
  rlab <- setNames(sprintf('r = %.2f', art_r), names(art_r))
  fs8h$label <- rlab[as.character(fs8h$genotype)]
  p_scatter <- ggplot2::ggplot(fs8h, ggplot2::aes(split_A_fraction, split_B_fraction)) +
    ggplot2::geom_point(size = .9, shape = 21, fill = 'grey75', color = 'grey40', stroke = .2) +
    ggplot2::facet_wrap(~genotype, nrow = 1) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(x = 'Clonal bias, split 1', y = 'Clonal bias, split 2',
                  title = 'early_Ery bias; day 5, rep 1') +
    ggplot2::theme_classic(base_size = 12)
  pdf(file.path(stage, 'figures', fig_names[2L]), width = 9.0, height = 4.5)
  print(p_scatter + ggplot2::geom_text(data = unique(fs8h[c('genotype','label')]),
              ggplot2::aes(x = .96, y = .05, label = label), inherit.aes = FALSE,
              hjust = 1, size = 4, fontface = 'italic'))
  dev.off()
  pdf(file.path(stage, 'figures', fig_names[3L]), width = 5.0, height = 4.3)
  print(plot_bars(fs8i, facet = FALSE)); dev.off()

  staged <- c(file.path(stage, 'source_data', src_names), file.path(stage, 'figures', fig_names))
  if (any(!file.exists(staged)) || any(file.info(staged)$size < 100L))
    stop('Missing or empty staged output; no files promoted')
  dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  if (any(file.exists(target_files))) stop('Destination collision before promotion')
  for (i in seq_along(staged)) {
    if (!file.rename(staged[i], target_files[i])) stop('File promotion failed: ', target_files[i])
  }
  out <- data.frame(file = c(src_names, fig_names),
                    rows = c(nrow(f5l), nrow(fs8h), nrow(fs8i), nrow(summ), nrow(points), NA, NA, NA),
                    status = c(rep('INTERNAL_PEARSON_CHECK_PASS; HISTORICAL_REVIEW_PENDING', 5),
                               rep('PDF_GENERATED; HISTORICAL_VISUAL_REVIEW_PENDING', 3)))
  print(out, row.names = FALSE)
  invisible(list(status = 'INTERNAL_PEARSON_CHECK_PASS; HISTORICAL_VISUAL_REVIEW_PENDING',
                 tables = out, art_r = art_r))
}
