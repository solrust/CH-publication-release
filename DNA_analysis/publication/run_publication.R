# CH publication reproduction — DNA_analysis
# Initial checkpoint stage: Figure 2C and Supplementary Figure S3J.
# Run from CH_publication_release. No raw sequencing or model fitting.
# Input historical scripts are read-only.

dna_pub_preflight <- function(root = getwd()) {
  if (!dir.exists(file.path(root, "DNA_analysis"))) stop("Run from CH_publication_release root")
  needed <- c(
    "DNA_analysis/data/processed/DNA_clones.RData",
    "RNA_analysis/data/processed/hashtag_meta.rds")
  missing <- needed[!file.exists(file.path(root, needed))]
  if (length(missing)) stop("Missing required files: ", paste(missing, collapse = ", "))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  invisible(TRUE)
}

dna_pub_get_clone_counts <- function(ct, norm, h) {
  states <- c("6", "12", "16", "20", "24", "CD34", "Myeloid", "B", "T")
  stopifnot(identical(rownames(ct), rownames(norm)), all(states %in% names(norm)),
            all(c("genotype", "sample") %in% names(ct)))
  meta <- data.frame(sample = as.character(ct$sample),
                     genotype = as.character(ct$genotype),
                     stringsAsFactors = FALSE)
  if (anyNA(meta)) stop("Missing clone sample/genotype metadata")
  if (is.null(rownames(h)) || !"exp" %in% names(h)) stop("hashtag metadata lacks exp")
  if (any(!meta$sample %in% rownames(h))) stop("Unmatched sample IDs in hashtag metadata")
  meta$exp <- as.character(h[meta$sample, "exp"])
  if (anyNA(meta$exp)) stop("Missing experiment labels")
  detected <- as.matrix(norm[, states, drop = FALSE]) >= 1e-4
  detected[is.na(detected)] <- FALSE
  # One row per sample/genotype/experiment/timepoint or compartment.
  sample_key <- paste(meta$sample, meta$genotype, meta$exp, sep = "\r")
  sums <- rowsum(detected * 1L, group = sample_key, reorder = FALSE)
  keys <- strsplit(rownames(sums), "\r", fixed = TRUE)
  out <- data.frame(sample = vapply(keys, `[`, "", 1),
                    genotype = vapply(keys, `[`, "", 2),
                    experiment = vapply(keys, `[`, "", 3),
                    sums, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out <- out[order(out$experiment, out$sample), , drop = FALSE]
  list(counts = out, meta = meta, detected = detected, states = states)
}

dna_pub_source_tables <- function(ct, norm, h) {
  x <- dna_pub_get_clone_counts(ct, norm, h)
  states <- x$states
  counts <- x$counts
  # Fig 2C historical 'show.HSPC=TRUE': CD34-detected clones / all CD34 clones,
  # at detection threshold >= 1e-4, separately per sample and compartment.
  cd34 <- x$detected[, "CD34"]
  meta <- x$meta[cd34, , drop = FALSE]
  subset_mat <- x$detected[cd34, , drop = FALSE] * 1L
  if (!nrow(subset_mat)) stop("No CD34-positive clones found")
  keys <- paste(meta$sample, meta$genotype, meta$exp, sep = "\r")
  nums <- rowsum(subset_mat, keys, reorder = FALSE)
  ix <- match(paste(counts$sample, counts$genotype, counts$experiment, sep = "\r"), rownames(nums))
  if (anyNA(ix)) stop("Some samples lack CD34-selected clone counts")
  nums <- nums[ix, states, drop = FALSE]
  denom <- counts$CD34
  if (any(denom < 1L)) stop("CD34 denominator zero in sample")
  # Keyed exact equality is expected for the HSPC denominator.
  if (!identical(as.numeric(nums[, "CD34"]), as.numeric(denom)))
    stop("CD34 numerator/denominator mismatch")
  make_long <- function(mat, value_name) {
    ans <- do.call(rbind, lapply(states, function(st) {
      data.frame(sample = counts$sample, genotype = counts$genotype,
                 experiment = counts$experiment, compartment = st,
                 value = as.numeric(mat[, st]), stringsAsFactors = FALSE)
    }))
    names(ans)[names(ans) == "value"] <- value_name
    ans[order(ans$experiment, ans$sample,
              match(ans$compartment, states)), , drop = FALSE]
  }
  j <- make_long(counts[, states, drop = FALSE], "barcodes_detected")
  c <- make_long(nums, "shared_CD34_barcodes")
  # Avoid alignment ambiguity by joining with verified sample IDs.
  c$CD34_denominator <- denom[match(c$sample, counts$sample)]
  c$fraction_of_CD34_clones <- c$shared_CD34_barcodes / c$CD34_denominator
  stopifnot(nrow(j) == nrow(c), all(c$fraction_of_CD34_clones >= 0 & c$fraction_of_CD34_clones <= 1))
  list(S3J = j, `2C` = c)
}

dna_pub_run <- function(root = getwd(), panels = c("S3J", "2C"), validate_replay = TRUE) {
  dna_pub_preflight(root)
  panels <- match.arg(panels, c("S3J", "2C"), several.ok = TRUE)
  pub <- file.path(root, "DNA_analysis/publication")
  dest <- file.path(pub, "source_data")
  figs <- file.path(pub, "figures")
  paths <- c(if ("S3J" %in% panels) c(file.path(dest, "FigS3J_clone_detection_counts.tsv.gz"),
                                         file.path(figs, "FigS3J_clone_detection_counts_REVIEW.pdf")),
             if ("2C" %in% panels) c(file.path(dest, "Fig2C_HSPC_clone_sharing.tsv.gz"),
                                       file.path(figs, "Fig2C_HSPC_clone_sharing_REVIEW.pdf")))
  if (any(file.exists(paths))) stop("Refusing overwrite of existing output: ", paste(paths[file.exists(paths)], collapse = ", "))
  e <- new.env(parent = baseenv())
  load(file.path(root, "DNA_analysis/data/processed/DNA_clones.RData"), envir = e)
  stopifnot(all(c("clone.trends", "clone.trends.norm", "clone_state", "growth.info.nbinom", "rate") %in% ls(e)))
  h <- readRDS(file.path(root, "RNA_analysis/data/processed/hashtag_meta.rds"))
  tables <- dna_pub_source_tables(e$clone.trends, e$clone.trends.norm, h)
  if (validate_replay) {
    rerun <- dna_pub_source_tables(e$clone.trends, e$clone.trends.norm, h)
    stopifnot(isTRUE(all.equal(tables$S3J, rerun$S3J, tolerance = 1e-8)),
              isTRUE(all.equal(tables$`2C`, rerun$`2C`, tolerance = 1e-8)))
  }
  # Plot display is provisional, requiring approval against assembled artwork.
  gg <- asNamespace("ggplot2")
  palette <- c(WT = "#333333", DNMT3A = "#dc3838", TET2 = "#409b4c")
  stage <- tempfile(pattern = ".dna_pub_stage_", tmpdir = pub)
  dir.create(stage, recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  dir.create(file.path(stage, "source_data"))
  dir.create(file.path(stage, "figures"))
  for (panel in panels) {
    if (panel == "S3J") {
      dat <- tables$S3J
      tsv <- "FigS3J_clone_detection_counts.tsv.gz"
      pdfname <- "FigS3J_clone_detection_counts_REVIEW.pdf"
      plot_dat <- dat
      plot_dat$compartment <- factor(plot_dat$compartment,
                                     levels = c("6", "12", "16", "20", "24", "CD34", "Myeloid", "B", "T"))
      p <- gg$ggplot(plot_dat, gg$aes(x = compartment, y = barcodes_detected, fill = genotype)) +
        gg$geom_col(position = gg$position_dodge(width = 0.8), color = "black", linewidth = 0.1) +
        gg$facet_wrap(~experiment, ncol = 1, scales = "free_y") +
        gg$scale_fill_manual(values = palette) +
        gg$labs(x = NULL, y = "Number of clones detected", title = "S3J — checkpoint reproduction") +
        gg$theme_bw() + gg$theme(axis.text.x = gg$element_text(angle = 45, hjust = 1))
    } else {
      dat <- tables$`2C`
      tsv <- "Fig2C_HSPC_clone_sharing.tsv.gz"
      pdfname <- "Fig2C_HSPC_clone_sharing_REVIEW.pdf"
      plot_dat <- dat
      plot_dat$compartment <- factor(plot_dat$compartment,
                                     levels = c("6", "12", "16", "20", "24", "CD34", "Myeloid", "B", "T"))
      p <- gg$ggplot(plot_dat, gg$aes(x = compartment, y = fraction_of_CD34_clones, fill = genotype)) +
        gg$geom_col(position = gg$position_dodge2(), color = "black", linewidth = 0.1) +
        gg$facet_wrap(~experiment, ncol = 1) + gg$scale_fill_manual(values = palette) +
        gg$labs(x = NULL, y = "Fraction of CD34 barcodes shared", title = "2C — checkpoint reproduction") +
        gg$theme_bw() + gg$theme(axis.text.x = gg$element_text(angle = 45, hjust = 1))
    }
    con <- gzfile(file.path(stage, "source_data", tsv), "wt")
    tryCatch(write.table(dat, con, sep = "\t", quote = FALSE, row.names = FALSE),
             finally = close(con))
    grDevices::pdf(file.path(stage, "figures", pdfname), width = 9, height = 6)
    tryCatch(print(p), finally = grDevices::dev.off())
    cat("CHECKPOINT_REPLAY_PASS:", panel, nrow(dat), "Source Data rows; PDF generated\n")
  }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  dir.create(figs, recursive = TRUE, showWarnings = FALSE)
  for (p in paths) {
    sub <- if (grepl("/source_data/", p, fixed = TRUE)) "source_data" else "figures"
    staged <- file.path(stage, sub, basename(p))
    stopifnot(file.exists(staged), file.info(staged)$size > 0)
    if (!file.rename(staged, p)) stop("Could not promote staged output: ", p)
  }
  invisible(list(status = "CHECKPOINT_REPLAY_PASS; HISTORICAL_VISUAL_REVIEW_PENDING",
                 tables = vapply(tables[panels], nrow, integer(1))))
}

# DNA_analysis publication extension: Figure 2D/E/F/G and Supplementary S3K.
# Requires the already-installed dna_pub_run() for Figure 2C/S3J.
# Preserves original inputs and approved outputs. Figures are for visual review.

dna_pub_remaining <- function(root = getwd(),
                              panels = c("2D", "2E", "2F", "2G", "S3K"),
                              validate_replay = TRUE) {
  dna_pub_preflight(root)
  panels <- unique(match.arg(panels, c("2D", "2E", "2F", "2G", "S3K"), several.ok = TRUE))
  require_packages <- c("ggplot2", "RColorBrewer")
  if ("2F" %in% panels) require_packages <- c(require_packages, "pheatmap")
  absent <- require_packages[!vapply(require_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing packages: ", paste(absent, collapse = ", "))
  pub <- file.path(root, "DNA_analysis", "publication")
  src <- file.path(pub, "source_data")
  fig <- file.path(pub, "figures")
  filemap <- list(
    `2D` = list(tables = "Fig2D_clone_trajectories.tsv.gz", plots = c("Fig2D_growth_up_REVIEW.pdf", "Fig2D_growth_down_REVIEW.pdf")),
    `2E` = list(tables = "Fig2E_growth_CD34_groups.tsv.gz", plots = "Fig2E_growth_heterogeneity_REVIEW.pdf"),
    `2F` = list(tables = "Fig2F_lineage_heatmap_clones.tsv.gz", plots = "Fig2F_lineage_heatmap_REVIEW.pdf"),
    `2G` = list(tables = "Fig2G_growing_clones_by_sample.tsv", plots = "Fig2G_growing_clones_REVIEW.pdf"),
    `S3K` = list(tables = "FigS3K_lineage_fraction_by_growth.tsv.gz", plots = "FigS3K_lineage_fraction_by_growth_REVIEW.pdf")
  )
  destfiles <- unlist(lapply(filemap[panels], function(x) c(file.path(src, x$tables), file.path(fig, x$plots))), use.names = FALSE)
  if (any(file.exists(destfiles))) stop("Refusing to overwrite: ", paste(destfiles[file.exists(destfiles)], collapse = ", "))
  e <- new.env(parent = baseenv())
  load(file.path(root, "DNA_analysis", "data", "processed", "DNA_clones.RData"), envir = e)
  required <- c("clone.trends", "clone.trends.norm", "growth.info.nbinom", "clone_state", "rate", "sample_size")
  if (!all(required %in% ls(e))) stop("Incomplete DNA_clones.RData checkpoint")
  ct <- e$clone.trends; norm <- e$clone.trends.norm
  gf <- e$growth.info.nbinom; clone_state <- e$clone_state; rate <- e$rate
  sample_size <- e$sample_size
  stopifnot(identical(rownames(ct), rownames(norm)), all(c("genotype", "sample") %in% names(ct)),
            all(c("CD34", "Myeloid", "B", "T", "24") %in% names(norm)))
  h <- readRDS(file.path(root, "RNA_analysis", "data", "processed", "hashtag_meta.rds"))
  gg <- asNamespace("ggplot2")
  geno_colors <- c(WT = "#444444", DNMT3A = "#e33636", TET2 = "#49a55b")
  timepoints <- c(6, 12, 16, 20, 24)
  pvalues <- vapply(gf, function(z) as.numeric(z$coefs[2, 4]), numeric(1))
  candidates <- intersect(names(pvalues), rownames(ct))
  active <- candidates[is.finite(pvalues[candidates]) & pvalues[candidates] != 0 &
                         rowSums(as.matrix(ct[candidates, 1:5, drop = FALSE]) > 4, na.rm = TRUE) > 2]
  active <- intersect(active, names(rate))
  # Historically selected clone set for growth-rate boxplots and trends.
  get_growth <- function() {
    vals <- data.frame(clone_id = active, sample = as.character(ct[active, "sample"]),
                       genotype = as.character(ct[active, "genotype"]),
                       rate = as.numeric(rate[active]), pvalue = as.numeric(pvalues[active]),
                       CD34 = as.numeric(norm[active, "CD34"]),
                       Myeloid = as.numeric(norm[active, "Myeloid"]),
                       B = as.numeric(norm[active, "B"]), T = as.numeric(norm[active, "T"]),
                       stringsAsFactors = FALSE)
    vals <- vals[is.finite(vals$rate), , drop = FALSE]
    q <- quantile(vals$rate, probs = seq(0, 1, length.out = 6), na.rm = TRUE)
    if (anyDuplicated(q)) stop("Growth quintile breaks are duplicated; inspect historical rate selection")
    vals$rate_category <- as.integer(cut(vals$rate, breaks = q, include.lowest = TRUE))
    med <- tapply(vals$rate, vals$rate_category, median)
    vals$median_rate <- as.numeric(med[as.character(vals$rate_category)])
    vals$CD34_status <- ifelse(vals$CD34 > 1e-5, "CD34+", "CD34-")
    vals <- vals[order(vals$clone_id), , drop = FALSE]
    rownames(vals) <- NULL
    vals
  }
  need_growth <- any(c("2E", "S3K") %in% panels)
  growth <- if (need_growth) get_growth() else NULL
  stage <- tempfile(pattern = ".dna_remaining_", tmpdir = pub)
  dir.create(file.path(stage, "source_data"), recursive = TRUE)
  dir.create(file.path(stage, "figures"), recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  stage_tab <- function(panel, dat) {
    stopifnot(is.data.frame(dat), nrow(dat) > 0, !anyDuplicated(names(dat)))
    nm <- filemap[[panel]]$tables
    path <- file.path(stage, "source_data", nm)
    con <- if (grepl("\\.gz$", nm)) gzfile(path, "wt") else file(path, "wt")
    tryCatch(write.table(dat, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"),
             finally = close(con))
    cat("DERIVED:", panel, nrow(dat), "Source Data rows\n")
  }
  plot_file <- function(name, expr, width = 8, height = 5) {
    grDevices::pdf(file.path(stage, "figures", name), width = width, height = height)
    tryCatch(force(expr), finally = grDevices::dev.off())
  }
  output_counts <- integer(0)
  if ("2D" %in% panels) {
    # Selected example clone IDs are explicit in the historical plotting script.
    up <- c("12_TTTATGTCGACACAGGACCTGTGACTCTGTGTAGAGGAAA",
            "4_TGCATGCTGACACCCTACCCGGGAGGGCGTGACGAGGCGG",
            "10_GAGGTGGACTCAAGATACCTGAGATCAGGTACATAGTAAT",
            "1_TCCATGAAACCACTTAACGTATGAAATAGTTCAGAGTTTC")
    down <- c("2_GTTATGCGGTCAAGATACTAACGACTCCGTGATGAGATTC",
              "6_CCAATGACCACACGTGACCATAGAAAGTGTAAATAGTCTT",
              "11_CATTTGTGACCACGGTACTTAGGACGAAGTTACTAGGATA",
              "8_AAACTGCTCACACCAAACTACCGATAAAGTAAAGAGCCGG")
    missing_highlights <- setdiff(c(up, down), names(gf))
    if (length(missing_highlights)) stop("Historical Figure 2D highlighted clones missing from checkpoint: ", paste(missing_highlights, collapse = ", "))
    plot_ids <- unique(c(active, up, down))
    get_traj <- function(id) {
      sample_id <- as.character(ct[id, "sample"])
      if (!sample_id %in% rownames(sample_size)) stop("Sample not in sample_size: ", sample_id)
      denom <- as.numeric(sample_size[sample_id, 1:5])
      if (any(!is.finite(denom) | denom <= 0)) stop("Invalid sample-size denominator: ", sample_id)
      counts <- as.numeric(ct[id, 1:5])
      coef <- as.numeric(gf[[id]]$coefs[1:2, 1])
      if (any(!is.finite(coef))) stop("Invalid stored growth model coefficients: ", id)
      data.frame(clone_id = id, sample = sample_id, genotype = as.character(ct[id,"genotype"]),
                 category = if (id %in% up) "highlight_up" else if (id %in% down) "highlight_down" else "background",
                 time_weeks = timepoints, barcode_count = counts, sample_denominator = denom,
                 observed_frequency = counts / denom * 1e-2,
                 fitted_frequency = exp(coef[1] + coef[2] * timepoints) * 1e-2,
                 growth_rate = as.numeric(rate[id]), growth_pvalue = as.numeric(pvalues[id]),
                 stringsAsFactors = FALSE)
    }
    traj <- do.call(rbind, lapply(plot_ids, get_traj))
    rownames(traj) <- NULL
    if (validate_replay) stopifnot(nrow(traj) == 5L * length(plot_ids),
                                   all(is.finite(traj$fitted_frequency)), all(traj$sample_denominator > 0))
    stage_tab("2D", traj); output_counts["2D"] <- nrow(traj)
    for (kind in c("up", "down")) {
      examples <- if (kind == "up") up else down
      pool <- active[is.finite(rate[active]) & (if (kind == "up") rate[active] > 0 else rate[active] < 0)]
      set.seed(2207 + if (kind == "up") 1 else 2)
      background <- sample(pool, min(200, length(pool)))
      plot_file(paste0("Fig2D_growth_", kind, "_REVIEW.pdf"), {
        graphics::plot(timepoints, rep(0, 5), type = "n", xlim = c(5, 27),
                       ylim = c(0, if (kind == "up") 0.025 else 0.005),
                       xlab = "Weeks", ylab = "Clone frequency in PB", bty = "n", xaxt = "n")
        graphics::axis(1, at = timepoints)
        for (id in background) {
          z <- traj[traj$clone_id == id, ]
          graphics::lines(z$time_weeks, z$fitted_frequency, col = grDevices::adjustcolor("grey50", .3))
        }
        cols <- RColorBrewer::brewer.pal(4, "Set2")
        for (i in seq_along(examples)) {
          z <- traj[traj$clone_id == examples[i], ]
          graphics::lines(z$time_weeks, z$fitted_frequency, col = cols[i], lwd = 3)
          graphics::points(z$time_weeks, z$observed_frequency, pch = 21, bg = cols[i], cex = 1.1)
        }
        graphics::legend(if (kind == "up") "topleft" else "topright", legend = paste("clone", if (kind == "up") 1:4 else 5:8),
                         col = cols, pch = 19, bty = "n", cex = .85)
      }, width = 5.5, height = 4.5)
    }
  }
  if ("2E" %in% panels) {
    stage_tab("2E", growth); output_counts["2E"] <- nrow(growth)
    z <- growth
    z$genotype <- factor(z$genotype, levels = c("WT", "DNMT3A", "TET2"))
    z$CD34_status <- factor(z$CD34_status, levels = c("CD34-", "CD34+"))
    p <- gg$ggplot(z, gg$aes(x = genotype, y = rate)) +
      gg$geom_boxplot(color = "black", outlier.shape = NA, width = .65) +
      gg$geom_point(gg$aes(fill = genotype), shape = 21, size = 1.3, color = "black",
                    position = gg$position_jitter(width = .2, height = 0, seed = 123)) +
      gg$facet_grid(~ CD34_status) + gg$geom_hline(yintercept = 0, linetype = 2) +
      gg$scale_fill_manual(values = geno_colors) +
      gg$labs(x = "Clone genotype", y = "Clone growth rate") + gg$theme_bw() +
      gg$theme(legend.position = "none")
    plot_file(filemap[["2E"]]$plots, print(p), width = 9, height = 4.7)
  }
  if ("2F" %in% panels) {
    vars <- c("CD34", "Myeloid", "B", "T")
    ids <- intersect(rownames(ct)[rowSums(as.matrix(ct[, vars, drop = FALSE]), na.rm = TRUE) > 100], names(rate))
    if (length(ids) < 3) stop("Too few lineage-heatmap clones")
    x <- as.matrix(norm[ids, vars, drop = FALSE])
    status <- ifelse(rate[ids] > 0, "up", "down")
    qvals <- p.adjust(pvalues[ids], "BH")
    status[!is.finite(rate[ids]) | !is.finite(pvalues[ids]) | is.na(qvals) | qvals > .2] <- "uncertain"
    side <- x[, "CD34"] > 1e-5
    get_order <- function(ii) if (length(ii) < 2L) ii else ii[stats::hclust(stats::dist(log(x[ii, , drop = FALSE] * 1e6 + 1)))$order]
    ord <- c(get_order(which(side)), get_order(which(!side)))
    z <- data.frame(clone_id = ids, sample = as.character(ct[ids,"sample"]),
                    genotype = as.character(ct[ids,"genotype"]), growth_rate = as.numeric(rate[ids]),
                    growth_pvalue = as.numeric(pvalues[ids]), BH_adjusted_pvalue = as.numeric(qvals),
                    growth_category = status, CD34_positive = side,
                    CD34 = x[,"CD34"], Myeloid = x[,"Myeloid"], B = x[,"B"], T = x[,"T"],
                    heatmap_order = match(seq_along(ids), ord), stringsAsFactors = FALSE)
    z <- z[order(z$heatmap_order), , drop = FALSE]; rownames(z) <- NULL
    if (validate_replay) stopifnot(!anyDuplicated(z$clone_id), setequal(z$clone_id, ids),
                                   all(z$CD34_positive == (z$CD34 > 1e-5)))
    stage_tab("2F", z); output_counts["2F"] <- nrow(z)
    mat <- log(as.matrix(z[, vars, drop = FALSE]) * 1e6 + 1)
    rownames(mat) <- z$clone_id
    annotation <- data.frame(genotype = z$genotype, growth = z$growth_category, row.names = z$clone_id)
    annotation$genotype <- factor(annotation$genotype, levels = c("WT", "DNMT3A", "TET2"))
    annotation$growth <- factor(annotation$growth, levels = c("uncertain", "up", "down"))
    plot_file(filemap[["2F"]]$plots, {
      pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
                         show_rownames = FALSE, annotation_row = annotation,
                         annotation_colors = list(genotype = geno_colors,
                                                  growth = c(uncertain = "grey80", up = "#e33636", down = "#3765b0")),
                         color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100),
                         border_color = NA, fontsize_col = 11,
                         main = "Lineage contribution by growth phenotype")
    }, width = 9, height = 6)
  }
  if ("2G" %in% panels) {
    ids <- intersect(names(clone_state)[clone_state == "late"], rownames(ct))
    ids <- intersect(ids[rowSums(as.matrix(ct[ids, 1:5, drop = FALSE]) > 2, na.rm = TRUE) > 2], names(rate))
    ids <- ids[is.finite(rate[ids])]
    if (!length(ids)) stop("No persistent late clones for Figure 2G")
    keys <- as.character(ct[ids, "sample"])
    bysample <- split(ids, keys)
    z <- do.call(rbind, lapply(names(bysample), function(s) {
      sub <- bysample[[s]]
      if (!s %in% rownames(h)) stop("Unmatched sample in 2G: ", s)
      data.frame(sample = s, genotype = as.character(h[s, "genotype"]),
                 experiment = as.character(h[s, "exp"]),
                 persistent_late_clones = length(sub),
                 growing_clones = sum(rate[sub] >= 0),
                 frequent_24_week_clones = sum(norm[sub, "24"] >= 1e-5, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
    z <- z[order(z$experiment, z$sample), , drop = FALSE]
    rownames(z) <- NULL
    if (validate_replay) stopifnot(sum(z$persistent_late_clones) == length(ids),
                                   all(z$growing_clones <= z$persistent_late_clones))
    stage_tab("2G", z); output_counts["2G"] <- nrow(z)
    z$genotype <- factor(z$genotype, levels = c("DNMT3A", "WT", "TET2"))
    z$experiment <- factor(z$experiment)
    p <- gg$ggplot(z, gg$aes(x = genotype, y = growing_clones,
                            shape = experiment, fill = genotype)) +
      gg$geom_point(size = 3.4, position = gg$position_jitter(width = .08, height = 0, seed = 123),
                    color = "black") +
      gg$scale_fill_manual(values = geno_colors) +
      gg$labs(x = "Clone genotype", y = "Number of growing clones") + gg$theme_classic() +
      gg$theme(legend.position = "right")
    plot_file(filemap[["2G"]]$plots, print(p), width = 6, height = 4.5)
  }
  if ("S3K" %in% panels) {
    z <- do.call(rbind, lapply(c("CD34", "B", "Myeloid"), function(comp) {
      data.frame(clone_id = growth$clone_id, genotype = growth$genotype,
                 sample = growth$sample, growth_rate = growth$rate,
                 growth_quintile = growth$rate_category, median_growth_rate = growth$median_rate,
                 compartment = comp, clone_fraction = growth[[comp]], stringsAsFactors = FALSE)
    }))
    rownames(z) <- NULL
    stage_tab("S3K", z); output_counts["S3K"] <- nrow(z)
    z$genotype <- factor(z$genotype, levels = c("WT", "DNMT3A", "TET2"))
    z$compartment <- factor(z$compartment, levels = c("CD34", "B", "Myeloid"))
    z$growth_quintile <- factor(z$growth_quintile, levels = 1:5)
    plot_file(filemap[["S3K"]]$plots, {
      for (comp in levels(z$compartment)) {
        d <- droplevels(z[z$compartment == comp, , drop = FALSE])
        # Historical scale_y_log10 omits zeros from graphical geometry; keep them in Source Data.
        p <- gg$ggplot(d, gg$aes(x = as.numeric(growth_quintile), y = clone_fraction,
                                group = interaction(growth_quintile, genotype))) +
          gg$geom_boxplot(gg$aes(color = genotype), width = .6, outlier.shape = NA,
                          position = gg$position_dodge(width = .75)) +
          gg$geom_point(gg$aes(fill = genotype), shape = 21, size = .7, color = "black",
                        position = gg$position_jitterdodge(jitter.width = .15, jitter.height = 0,
                                                             dodge.width = .75, seed = 123)) +
          gg$geom_smooth(gg$aes(x = as.numeric(growth_quintile), y = clone_fraction,
                                 color = genotype, group = genotype), method = "lm", se = TRUE,
                         inherit.aes = FALSE) +
          gg$scale_y_log10() +
          gg$scale_x_continuous(breaks = 1:5,
                                labels = sprintf("%.2f", tapply(d$growth_rate, d$growth_quintile, median))) +
          gg$scale_color_manual(values = geno_colors) + gg$scale_fill_manual(values = geno_colors) +
          gg$labs(x = "Binned clone growth rate", y = paste0(comp, " clone fraction")) + gg$theme_bw()
        print(p)
      }
    }, width = 7, height = 5)
  }
  if (validate_replay) stopifnot(length(output_counts) == length(panels), all(output_counts > 0))
  dir.create(src, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig, recursive = TRUE, showWarnings = FALSE)
  for (path in destfiles) {
    sub <- if (grepl("/source_data/", path, fixed = TRUE)) "source_data" else "figures"
    temp <- file.path(stage, sub, basename(path))
    if (!file.exists(temp) || file.info(temp)$size <= 0) stop("Staged output missing/empty: ", temp)
  }
  for (path in destfiles) {
    sub <- if (grepl("/source_data/", path, fixed = TRUE)) "source_data" else "figures"
    temp <- file.path(stage, sub, basename(path))
    if (!file.rename(temp, path)) stop("Failed to promote output: ", path)
  }
  print(data.frame(panel = names(output_counts), Source_Data_rows = as.integer(output_counts),
                   status = "CHECKPOINT_DERIVATION_INTERNAL_PASS", row.names = NULL))
  invisible(list(status = "CHECKPOINT_DERIVATION_INTERNAL_PASS; HISTORICAL_VISUAL_REVIEW_PENDING",
                 tables = output_counts))
}


# DNA_analysis Fig2B publication extension (checkpoint only).
# Uses historical script 1_DNA_analysis.R: nine clone compartments, ln(frequency*1e6+1),
# Euclidean hierarchical clustering with complete linkage, saved clone_state annotations.
# Does not refit growth models or modify any historic inputs/figures.

dna_pub_fig2b <- function(root = getwd()) {
  dna_pub_preflight(root)
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("pheatmap required")
  pub <- file.path(root, "DNA_analysis", "publication")
  src <- file.path(pub, "source_data", "Fig2B_clone_dynamics.tsv.gz")
  pdf <- file.path(pub, "figures", "Fig2B_clone_dynamics_REVIEW.pdf")
  if (any(file.exists(c(src, pdf)))) stop("Figure 2B outputs already exist; refusing overwrite")
  stopifnot(dir.exists(dirname(src)), dir.exists(dirname(pdf)))
  e <- new.env(parent = baseenv())
  load(file.path(root, "DNA_analysis", "data", "processed", "DNA_clones.RData"), envir = e)
  stopifnot(all(c("clone.trends", "clone.trends.norm", "clone_state") %in% ls(e)))
  ct <- e$clone.trends
  norm <- e$clone.trends.norm
  states <- c("6", "12", "16", "20", "24", "CD34", "Myeloid", "B", "T")
  stopifnot(identical(rownames(ct), rownames(norm)), all(states %in% colnames(norm)))
  ids <- rownames(ct)[rownames(ct) %in% names(e$clone_state)]
  if (length(ids) < 100L) stop("Unexpectedly few selected clones; inspect clone_state IDs")
  status <- as.character(e$clone_state[ids])
  genotype <- as.character(ct[ids, "genotype"])
  if (anyNA(status) || anyNA(genotype)) stop("Missing clone annotation")
  freq <- as.matrix(norm[ids, states, drop = FALSE])
  storage.mode(freq) <- "double"
  if (anyNA(freq) || any(!is.finite(freq)) || any(freq < 0)) stop("Nonfinite or negative clone frequencies")
  transformed <- log(freq * 1e6 + 1)
  message("Fig2B: clustering ", nrow(transformed), " clones over ", ncol(transformed), " columns; no growth refit")
  # pheatmap historical default clustering: Euclidean distance, complete linkage.
  hc <- stats::hclust(stats::dist(transformed), method = "complete")
  hc_order <- hc$order
  stopifnot(length(hc_order) == length(ids), !anyDuplicated(hc_order))
  annotations <- data.frame(status = factor(status, levels = c("early", "late")),
                            genotype = factor(genotype, levels = c("WT", "DNMT3A", "TET2")),
                            row.names = ids)
  if (anyNA(annotations)) stop("Unknown historical status or genotype labels")
  record <- data.frame(clone_id = ids, sample = as.character(ct[ids, "sample"]),
                       genotype = genotype, status = status,
                       heatmap_order = match(seq_along(ids), hc_order),
                       freq, check.names = FALSE, stringsAsFactors = FALSE)
  record <- record[base::order(record$heatmap_order), , drop = FALSE]
  rownames(record) <- NULL
  stopifnot(nrow(record) == length(ids), !anyDuplicated(record$clone_id))
  # Stage outputs before publication. Copy each only after both are complete.
  stage <- tempfile(pattern = ".dna2b_", tmpdir = pub)
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  tabtmp <- file.path(stage, basename(src))
  pdftmp <- file.path(stage, basename(pdf))
  con <- gzfile(tabtmp, open = "wt")
  tryCatch(utils::write.table(record, con, sep = "\t", row.names = FALSE,
                              col.names = TRUE, quote = FALSE, na = "NA"),
           finally = close(con))
  if (!file.exists(tabtmp) || file.info(tabtmp)$size == 0) stop("Failed to stage Source Data")
  p <- pheatmap::pheatmap(
    transformed, cluster_rows = hc, cluster_cols = FALSE,
    show_rownames = FALSE, annotation_row = annotations,
    annotation_colors = list(
      genotype = c(WT = "#777777", DNMT3A = "#ee7777", TET2 = "#88bb77"),
      status = c(early = "grey80", late = "grey20")),
    color = grDevices::colorRampPalette(c("lightgrey", "darkred"))(100),
    legend_breaks = log(10^c(-6, -4, -2, -1) * 1e6 + 1),
    legend_labels = c("1e-6", "1e-4", "1e-2", "0.1"),
    border_color = NA, angle_col = 45,
    fontsize_col = 10, silent = TRUE)
  grDevices::pdf(pdftmp, width = 9, height = 6, useDingbats = FALSE)
  tryCatch({
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  }, finally = grDevices::dev.off())
  if (!file.exists(pdftmp) || file.info(pdftmp)$size < 1000) stop("Invalid staged Figure 2B PDF")
  if (file.exists(src) || file.exists(pdf)) stop("Outputs appeared during staging; refusing overwrite")
  if (!file.copy(tabtmp, src, overwrite = FALSE)) stop("Could not publish Source Data")
  if (!file.copy(pdftmp, pdf, overwrite = FALSE)) {
    unlink(src)
    stop("Could not publish figure PDF; removed new Source Data to avoid incomplete pair")
  }
  message("Fig2B: ", nrow(record), " clone rows; checkpoint-derived table and PDF generated")
  list(status = "CHECKPOINT_DERIVATION_INTERNAL_PASS; HISTORICAL_VISUAL_REVIEW_PENDING",
       table_rows = nrow(record), source_data = src, figure = pdf)
}
