# TF_overexpression/publication/run_publication.R
# CH publication, checkpoint-based reproduction: Fig. 7C-D, 7E;
# Supplementary Fig. 12B, 12E and Supplementary Fig. 13A-C
# Source: TF_overexpression/scripts/1_dna_processing.R, historical plots
#  figures/volcano_bias.pdf and figures/volcano_mixed_18.pdf.
# Run from CH_publication_release; do not source historical scripts end-to-end.
# Historical inputs are read-only. No new model fitting or reclustering.
# The heavy TF expression matrix is discarded as soon as cell.map.tf and TFn
# are recovered. scde is NOT required for checkpoint-based plotting.

.tf_assert_root <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  required <- file.path(root, "TF_overexpression", "data", "processed", "TF_trends.RData")
  if (!file.exists(required)) {
    stop("Cannot locate TF_trends.RData. Set 'root' to CH_publication_release.", call. = FALSE)
  }
  root
}

.tf_stage3_panels <- function() c("7E", "S12E", "S13A", "S13B", "S13C")

.tf_load_historical_matrix <- function() {
  # Read historic serialized sparse matrices only after pinning the Matrix
  # namespace. A wrong already-loaded namespace cannot be safely swapped.
  if (!"Matrix" %in% loadedNamespaces()) {
    if (!dir.exists(file.path(Sys.getenv("CH_R_LIB", unset = "/packages"), "Matrix"))) {
      stop("Historical Matrix library in CH_R_LIB not found. ",
           "Initialize the CH R 4.2.3 environment first.", call. = FALSE)
    }
    library(Matrix, lib.loc = Sys.getenv("CH_R_LIB", unset = "/packages"))
  }
  path <- normalizePath(getNamespaceInfo("Matrix", "path"), mustWork = TRUE)
  expected <- normalizePath(file.path(Sys.getenv("CH_R_LIB", unset = "/packages"), "Matrix"), mustWork = TRUE)
  if (!identical(path, expected) ||
      getNamespaceVersion("Matrix") < package_version("1.6-4")) {
    stop("Historical sparse matrix requires Matrix >=1.6-4 from ", expected,
         "; loaded ", path, " (", as.character(getNamespaceVersion("Matrix")),
         "). Restart R, then library(Matrix, lib.loc='/packages').", call. = FALSE)
  }
  invisible(TRUE)
}

tf_preflight <- function(panels = c("7C", "7D"), root = getwd()) {
  valid <- c("7C", "7D", "S12B", .tf_stage3_panels())
  if (!length(panels) || anyNA(panels) || any(!panels %in% valid)) {
    stop("Supported panels: ", paste(valid, collapse = ", "), ". Requested: ",
         paste(panels, collapse = ", "), call. = FALSE)
  }
  root <- .tf_assert_root(root)
  paths <- c(
    trends = file.path(root, "TF_overexpression/data/processed/TF_trends.RData"),
    invitro = file.path(root, "TF_overexpression/data/processed/TF_invitro.RData"),
    pagoda = file.path(root, "TF_overexpression/data/processed/p2_tfinvitro.rds"),
    dna = file.path(root, "TF_overexpression/data/raw/3_Count_161_HSC_filt.txt")
  )
  stage3 <- any(panels %in% .tf_stage3_panels())
  use_pagoda <- any(panels %in% c("7E", "S12E", "S13C"))
  needed <- c(if (any(panels %in% c("7C", "7D", "S12B"))) "trends",
              if (stage3) "invitro", if (use_pagoda) "pagoda",
              if ("S13B" %in% panels) "dna")
  input <- data.frame(input = needed, present = file.exists(paths[needed]),
                      MiB = round(file.info(paths[needed])$size / 1024^2, 3),
                      path = unname(paths[needed]), stringsAsFactors = FALSE)
  print(input, row.names = FALSE)
  if (!all(input$present)) {
    stop("Required checkpoint(s) missing: ",
         paste(input$path[!input$present], collapse = "; "), call. = FALSE)
  }
  if (stage3) .tf_load_historical_matrix()
  pkg <- unique(c("ggplot2", "ggrepel", if (stage3) "Matrix",
                  if (use_pagoda) "pagoda2", if ("7E" %in% panels) "RColorBrewer"))
  package_status <- data.frame(package = pkg,
                               available = vapply(pkg, requireNamespace, logical(1), quietly = TRUE))
  print(package_status, row.names = FALSE)
  if (!all(package_status$available)) {
    stop("Missing packages: ", paste(pkg[!package_status$available], collapse = ", "),
         call. = FALSE)
  }
  invisible(list(root = root, checkpoint = paths[["trends"]],
                 paths = paths, packages = package_status))
}

.tf_load_checkpoint <- function(path) {
  env <- new.env(parent = baseenv())
  loaded <- load(path, envir = env)
  need <- c("biastot", "tf.trends.no25.v1m", "bias_site", "barcode2tf")
  if (!all(need %in% loaded)) {
    stop("Checkpoint missing required objects: ",
         paste(setdiff(need, loaded), collapse = ", "), call. = FALSE)
  }
  if (!is.data.frame(env$biastot) || !is.data.frame(env$tf.trends.no25.v1m)) {
    stop("Expected biastot and tf.trends.no25.v1m to be data frames.", call. = FALSE)
  }
  list(bias = env$biastot, growth = env$tf.trends.no25.v1m,
       bias_site = env$bias_site, barcode2tf = env$barcode2tf)
}

.tf_source_7C <- function(ck) {
  # This is the 'WT' facet of historical figures/volcano_bias.pdf.
  # 'trend' was assigned in 1_dna_processing.R with BH FDR < 0.2
  # separately per comparison. Preserve the saved label, not a new cutoff.
  dat <- ck$bias
  required <- c("logFC", "pv", "TF", "trend", "comparison")
  if (!all(required %in% names(dat))) {
    stop("biastot lacks: ", paste(setdiff(required, names(dat)), collapse = ", "))
  }
  dat <- dat[!is.na(dat$comparison) & as.character(dat$comparison) == "WT", , drop = FALSE]
  if (!nrow(dat)) stop("No WT rows in biastot.", call. = FALSE)
  if (!all(as.character(dat$trend) %in% c("B", "myeloid", "no"))) {
    stop("Unexpected lineage-bias trend labels in checkpoint.", call. = FALSE)
  }
  data.frame(source_row = seq_len(nrow(dat)), TF = as.character(dat$TF),
             comparison = as.character(dat$comparison),
             logFC = as.numeric(dat$logFC), pvalue = as.numeric(dat$pv),
             trend = as.character(dat$trend),
             stringsAsFactors = FALSE)
}

.tf_source_7D <- function(ck) {
  # Historical 0–18 week growth panel uses tf.trends.no25.v1m,
  # and calls up/down at nominal p.t < 0.01 (no multiple-test correction).
  dat <- ck$growth
  required <- c("intercept", "t", "tDNMT3A", "tTET2", "p.intercept",
                "p.t", "p.tDNMT3A", "p.tTET2", "TF")
  if (!all(required %in% names(dat))) {
    stop("Growth checkpoint lacks: ",
         paste(setdiff(required, names(dat)), collapse = ", "), call. = FALSE)
  }
  trend <- rep("no", nrow(dat))
  trend[!is.na(dat$p.t) & dat$p.t < 1e-2 & dat$t > 0] <- "up"
  trend[!is.na(dat$p.t) & dat$p.t < 1e-2 & dat$t < 0] <- "down"
  data.frame(barcode = rownames(dat), TF = as.character(dat$TF),
             comparison = "WT", interval_weeks = "0-18",
             growth_rate = as.numeric(dat$t), pvalue = as.numeric(dat$p.t),
             trend = trend,
             intercept = as.numeric(dat$intercept),
             tDNMT3A = as.numeric(dat$tDNMT3A),
             tTET2 = as.numeric(dat$tTET2),
             p_intercept = as.numeric(dat$p.intercept),
             p_tDNMT3A = as.numeric(dat$p.tDNMT3A),
             p_tTET2 = as.numeric(dat$p.tTET2),
             stringsAsFactors = FALSE)
}

.tf_figure_7C <- function(dat) {
  plotdat <- dat[is.finite(dat$logFC) & is.finite(dat$pvalue) &
                   dat$pvalue > 0, , drop = FALSE]
  if (!nrow(plotdat)) stop("No plottable Fig7C rows.", call. = FALSE)
  plotdat$trend <- factor(plotdat$trend, levels = c("B", "myeloid", "no"))
  ggplot2::ggplot(plotdat,
                  ggplot2::aes(x = logFC, y = pvalue, label = TF, fill = trend)) +
    ggplot2::geom_point(color = "black", alpha = 0.8, size = 2, shape = 21) +
    ggrepel::geom_text_repel(max.overlaps = 20, size = 3, na.rm = TRUE, seed = 1) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-3, 3)) +
    ggplot2::scale_fill_manual(values = c(B = "darkgreen", myeloid = "darkorange", no = "grey"),
                               drop = FALSE) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
    ggplot2::labs(x = "B/My lineage bias, log", y = "p-value", fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}

.tf_figure_7D <- function(dat) {
  plotdat <- dat[is.finite(dat$growth_rate) & is.finite(dat$pvalue) &
                   dat$pvalue > 0, , drop = FALSE]
  if (!nrow(plotdat)) stop("No plottable Fig7D rows.", call. = FALSE)
  plotdat$trend <- factor(plotdat$trend, levels = c("down", "no", "up"))
  ggplot2::ggplot(plotdat,
                  ggplot2::aes(x = growth_rate, y = pvalue, label = TF, fill = trend)) +
    ggplot2::geom_point(color = "black", alpha = 0.8, size = 2, shape = 21) +
    ggrepel::geom_text_repel(max.overlaps = 15, size = 3, na.rm = TRUE, seed = 1) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-0.5, 0.5)) +
    ggplot2::scale_fill_manual(values = c(down = "blue", no = "grey", up = "red"),
                               drop = FALSE) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
    ggplot2::labs(x = "Growth rate (0-18 weeks)", y = "p-value", fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}


# ---- Supplementary Figure 12B ----
# Artwork labels are held below ONLY as an independent validation reference.
# The S12B selection and order are reconstructed automatically from the
# saved barcode-level results and reference-library barcode order.
# Each TF label can display multiple independent barcode-level estimates.
.tf_s12b_artwork_order <- function() {
  c("TCF23",
    "FXN",
    "HLF",
    "HOXC10",
    "ZSCAN21",
    "RPA2",
    "GTF2H3",
    "ZBTB43",
    "HEY1",
    "CDK1",
    "ID1",
    "ID3",
    "ID4",
    "IKZF4",
    "ETV7",
    "ZKSCAN3",
    "ZNF669",
    "ZNF710",
    "PAX3",
    "TFDP1",
    "PDCD2",
    "ZNF468",
    "ETV3",
    "TCF12",
    "YAP1",
    "CREB3L4",
    "ZNF528",
    "ZNF525",
    "SOHLH2",
    "FLI1",
    "PPARD",
    "MEIS2",
    "NR2C1",
    "ETS1",
    "CNBP",
    "TBP",
    "LDB2",
    "RNF2",
    "PMS1",
    "HMG20A",
    "MAEL",
    "PRRX2",
    "DMBX1",
    "HMBOX1",
    "HSFX2",
    "ZNF124",
    "HOXB3",
    "LHX6",
    "LHX2",
    "SHOX",
    "PIAS1",
    "ZNF155",
    "TEF",
    "GCFC2",
    "ZSCAN1",
    "PBX1",
    "PBX3",
    "CTCFL",
    "FOXB1",
    "CREM",
    "TCF7",
    "DBP",
    "LIN28B")
}


# Reconstructed rule independently reproduces the 63 TF names and their
# top-to-bottom order in approved Supplementary Fig. 12B:
# (1) retain TF names with >=1 barcode having nominal P<0.05 in the
#     lineage-bias OR growth arm (no FDR correction for selection);
# (2) order names by their earliest *tested* barcode position in the
#     original barcode2tf reference library, decreasing.
# This was inferred from the checkpoint and verified against the artwork;
# the original figure-assembly code has not been located.
.tf_s12b_reconstructed_order <- function(bias, growth, map) {
  required <- c("Sequence", "Gene")
  if (!all(required %in% names(map))) {
    stop("S12B: barcode reference requires Sequence and Gene.", call. = FALSE)
  }
  eligible <- function(df) {
    !is.na(df$TF) & nzchar(df$TF) &
      is.finite(df$pvalue) & df$pvalue < 0.05
  }
  selected <- union(unique(bias$TF[eligible(bias)]),
                    unique(growth$TF[eligible(growth)]))
  all_barcode <- c(as.character(bias$barcode), as.character(growth$barcode))
  all_tf <- c(as.character(bias$TF), as.character(growth$TF))
  reference_index <- match(all_barcode, as.character(map$Sequence))
  if (anyNA(reference_index) || anyNA(all_tf) ||
      !identical(all_tf, as.character(map$Gene[reference_index]))) {
    stop("S12B: original barcode reference does not uniquely map to all tested TFs.",
         call. = FALSE)
  }
  first_position <- tapply(reference_index, all_tf, min)
  if (!length(selected) || !all(selected %in% names(first_position))) {
    stop("S12B: missing selected TF(s) from reference barcode ordering.",
         call. = FALSE)
  }
  reconstructed <- names(sort(first_position[selected], decreasing = TRUE))
  artwork_reference <- .tf_s12b_artwork_order()
  if (!identical(reconstructed, artwork_reference)) {
    stop("S12B: reconstructed selection/order no longer matches approved artwork. ",
         "selected=", length(reconstructed),
         "; missing=", paste(setdiff(artwork_reference, reconstructed), collapse = ","),
         "; extra=", paste(setdiff(reconstructed, artwork_reference), collapse = ","),
         "; mismatched positions=",
         if (length(reconstructed) == length(artwork_reference))
           sum(reconstructed != artwork_reference) else NA_integer_,
         call. = FALSE)
  }
  cat("S12B selection: nominal P<0.05 in either barcode-level arm; ",
      length(reconstructed), " TF names; descending earliest tested ",
      "reference barcode position; independent artwork order check ",
      length(reconstructed), "/", length(artwork_reference), " PASS.\n", sep = "")
  reconstructed
}

.tf_source_S12B <- function(ck) {
  # Preserve separate barcode estimates for TFs with multiple barcodes.
  # Do not average by TF or select a representative barcode.
  # The approved artwork is used only for a post-reconstruction assertion.
  bias <- .tf_source_7C(ck)
  growth <- .tf_source_7D(ck)
  raw_bias <- ck$bias_site
  map <- ck$barcode2tf

  # Verify identity AND row order before recovering barcode sequences.
  if (is.null(rownames(raw_bias)) || anyDuplicated(rownames(raw_bias)) ||
      nrow(bias) != nrow(raw_bias) || ncol(raw_bias) < 2L ||
      !all(c("Sequence", "Gene") %in% names(map)) ||
      !isTRUE(all.equal(bias$logFC, as.numeric(raw_bias[, 1]),
                        check.attributes = FALSE)) ||
      !isTRUE(all.equal(bias$pvalue, as.numeric(raw_bias[, 2]),
                        check.attributes = FALSE))) {
    stop("S12B: bias_site and WT biastot rows do not align; " ,
         "refusing to infer barcode identity.", call. = FALSE)
  }
  bias$barcode <- rownames(raw_bias)
  mapped <- as.character(map$Gene[match(bias$barcode, map$Sequence)])
  if (anyNA(mapped) || !identical(bias$TF, mapped)) {
    stop("S12B: barcode2tf mapping does not match WT biastot TF names.",
         call. = FALSE)
  }
  if (is.null(growth$barcode) || anyNA(growth$barcode) ||
      anyDuplicated(growth$barcode)) {
    stop("S12B: growth barcode identifiers are missing or duplicated.",
         call. = FALSE)
  }
  # Independently infer the TF universe and exact display order.
  tf_order <- .tf_s12b_reconstructed_order(bias, growth, map)

  # BH-adjusted values are reported for audit but do NOT determine the
  # S12B colors (the manuscript's S12B legend states nominal P < 0.05).
  bias$FDR <- stats::p.adjust(bias$pvalue, "BH")
  growth$FDR <- stats::p.adjust(growth$pvalue, "BH")

  b <- bias[bias$TF %in% tf_order, , drop = FALSE]
  g <- growth[growth$TF %in% tf_order, , drop = FALSE]
  if (!nrow(b) || !nrow(g)) {
    stop("S12B: no barcode-level data for one artwork panel.", call. = FALSE)
  }
  label_class <- function(value, pval, pos, neg) {
    cls <- rep("no", length(value))
    hit <- is.finite(value) & is.finite(pval) & pval < 0.05
    cls[hit & value > 0] <- pos
    cls[hit & value < 0] <- neg
    cls
  }
  b_cls <- label_class(b$logFC, b$pvalue, "B", "myeloid")
  g_cls <- label_class(g$growth_rate, g$pvalue, "grow", "decrease")

  left <- data.frame(
    TF_position_top_to_bottom = match(b$TF, tf_order),
    TF = b$TF, arm = "lineage_bias", barcode = b$barcode,
    source_row = b$source_row,
    effect = b$logFC, nominal_pvalue = b$pvalue, BH_FDR = b$FDR,
    original_Fig7_class = b$trend,
    S12B_nominal_p005_class = b_cls,
    stringsAsFactors = FALSE)
  right <- data.frame(
    TF_position_top_to_bottom = match(g$TF, tf_order),
    TF = g$TF, arm = "growth_rate", barcode = g$barcode,
    source_row = match(g$barcode, rownames(ck$growth)),
    effect = g$growth_rate, nominal_pvalue = g$pvalue, BH_FDR = g$FDR,
    original_Fig7_class = g$trend,
    S12B_nominal_p005_class = g_cls,
    stringsAsFactors = FALSE)
  dat <- rbind(left, right)
  dat$tf_order_source <- "reconstructed_from_nominal_p_lt_0.05_union_and_descending_reference_barcode_position"
  dat$TF_selection_rule <- "at_least_one_barcode_nominal_P_lt_0.05_in_either_arm"
  dat$TF_order_rule <- "descending_earliest_tested_barcode_position_in_reference"
  dat$effect_definition <- ifelse(dat$arm == "lineage_bias",
                                  "B/My log-link model coefficient",
                                  "growth log-link time coefficient, weeks 0-18")

  # The position column now comes from the algorithm, not the artwork.
  # The published plotting windows exclude extreme/undefined effect values.
  # Do not discard them from Source Data; leave plot_include explicit.
  dat$plot_include <- is.finite(dat$effect) &
    ifelse(dat$arm == "lineage_bias",
           dat$effect >= -3 & dat$effect <= 2.5,
           dat$effect >= -0.7 & dat$effect <= 0.5)

  # Explicitly position each barcode as a separate short bar in its named
  # TF row. Order nonsignificant barcodes above more significant barcodes;
  # the precise within-TF artwork order has not been recovered.
  dat$bar_ymin <- NA_real_
  dat$bar_ymax <- NA_real_
  for (arm_now in c("lineage_bias", "growth_rate")) {
    for (tf_now in tf_order) {
      idx <- which(dat$arm == arm_now & dat$TF == tf_now)
      if (!length(idx)) next
      idx <- idx[order(dat$nominal_pvalue[idx], decreasing = TRUE,
                       na.last = TRUE)]
      n <- length(idx)
      row_center <- length(tf_order) + 1 - match(tf_now, tf_order)
      lane_height <- 0.78 / n
      for (j in seq_along(idx)) {
        center <- row_center + 0.39 - (j - 0.5) * lane_height
        dat$bar_ymin[idx[j]] <- center - 0.46 * lane_height
        dat$bar_ymax[idx[j]] <- center + 0.46 * lane_height
      }
    }
  }
  dat <- dat[order(dat$TF_position_top_to_bottom, dat$arm,
                   dat$nominal_pvalue, na.last = TRUE), , drop = FALSE]
  rownames(dat) <- NULL
  if (anyDuplicated(paste(dat$arm, dat$barcode, sep = ":"))) {
    stop("S12B: duplicate arm:barcode keys after export.", call. = FALSE)
  }
  for (arm_now in c("lineage_bias", "growth_rate")) {
    sub <- dat[dat$arm == arm_now, , drop = FALSE]
    duplicated_TFs <- unique(sub$TF[duplicated(sub$TF)])
    n_differ <- sum(sub$original_Fig7_class != sub$S12B_nominal_p005_class,
                    na.rm = TRUE)
    cat("S12B ", arm_now, ": ", nrow(sub), " barcode rows from ",
        length(unique(sub$TF)), " TF names; ", length(duplicated_TFs),
        " TFs have multiple barcodes; ",
        sum(sub$plot_include), " bars within the review plotting window; ",
        n_differ, " barcode classes differ from Fig7.\n", sep = "")
    if (length(duplicated_TFs)) {
      cat("  Multi-barcode TFs: ",
          paste(duplicated_TFs, collapse = ", "), "\n", sep = "")
    }
  }
  cat("S12B: 63 TFs selected algorithmically from the historical ",
      "checkpoint; one Source Data row = one barcode in one arm.\n", sep = "")
  dat
}

.tf_figure_S12B <- function(dat) {
  # Plot uses only the order reconstructed from the checkpoint's Source Data.
  tf_order <- unique(as.character(dat$TF[order(dat$TF_position_top_to_bottom)]))
  if (!length(tf_order) || anyNA(tf_order) ||
      length(tf_order) != length(unique(dat$TF))) {
    stop("S12B: cannot reconstruct plotting order from exported records.",
         call. = FALSE)
  }
  # A rectangle per barcode allows multiple barcodes in one TF-labelled row.
  # Using geom_col(aes(y=TF)) would overplot or combine independent values.
  common <- ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.title.y = ggplot2::element_blank())
  make_plot <- function(which_arm, title, axis_name, colors, x_breaks, x_range,
                        show_labels) {
    tbl <- dat[dat$arm == which_arm & dat$plot_include, , drop = FALSE]
    if (!nrow(tbl)) stop("S12B: no plot-ready barcodes in ", which_arm)
    tbl$class <- factor(tbl$S12B_nominal_p005_class,
                        levels = names(colors))
    p <- ggplot2::ggplot(tbl, ggplot2::aes(
      xmin = pmin(0, effect), xmax = pmax(0, effect),
      ymin = bar_ymin, ymax = bar_ymax, fill = class)) +
      ggplot2::geom_rect(color = "black", linewidth = 0.15) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.35) +
      ggplot2::scale_y_continuous(
        limits = c(0.5, length(tf_order) + 0.5),
        breaks = seq_along(tf_order), labels = rev(tf_order),
        expand = c(0, 0)) +
      ggplot2::scale_x_continuous(limits = x_range, breaks = x_breaks) +
      ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
      ggplot2::labs(x = axis_name, title = title) + common
    if (!show_labels) {
      p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                              axis.ticks.y = ggplot2::element_blank())
    }
    p
  }
  left <- make_plot(
    "lineage_bias", "Lineage bias", "B/My bias, log",
    c(B = "darkgreen", myeloid = "darkorange", no = "grey75"),
    c(-2, -1, 0, 1, 2), c(-3, 2.5), TRUE)
  right <- make_plot(
    "growth_rate", "Growth", "Growth rate, 0-18 weeks",
    c(grow = "red", decrease = "royalblue4", no = "grey75"),
    c(-0.6, -0.4, -0.2, 0, 0.2, 0.4), c(-0.7, 0.5), FALSE)
  list(left = left, right = right)
}

.tf_write_figure_S12B <- function(plots, path) {
  grDevices::pdf(path, width = 8.5, height = 13.5, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 1, ncol = 2, widths = grid::unit(c(0.57, 0.43), "null"))))
  print(plots$left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plots$right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
}


# ---- In-vitro TF cell assignment and saved UMAP (7E/S12E/S13A-C) ----
# Original analysis: scripts/4_tf_invitro_processing.R, lines 140-217;
# scripts/5_tf_invitro_analysis.R, lines 11-33.
.tf_load_stage3 <- function(paths, panels) {
  .tf_load_historical_matrix()
  e <- new.env(parent = baseenv())
  loaded <- load(paths[["invitro"]], envir = e)
  if (!all(c("cell.map.tf", "TFn") %in% loaded)) {
    stop("TF_invitro.RData missing cell.map.tf or TFn.", call. = FALSE)
  }
  m <- e$cell.map.tf
  TFn <- e$TFn
  need <- c("cell", "TFbarcode", "TF", "dominant", "total", "id")
  if (!all(need %in% names(m)) || anyDuplicated(as.character(m$cell))) {
    stop("Unexpected cell.map.tf schema or duplicated cell IDs.", call. = FALSE)
  }
  observed <- as.logical(m$id)
  threshold <- is.finite(m$dominant) & is.finite(m$total) &
    m$dominant > 10 & m$dominant > 0.9 * m$total
  if (anyNA(observed) || !identical(unname(observed), unname(threshold))) {
    stop("Saved TF assignments do not match historical dominant>10 and >90% rule.",
         call. = FALSE)
  }
  # Release the approximately 0.9-GiB raw count matrix before loading p2.
  rm(list = setdiff(ls(e), c("cell.map.tf", "TFn")), envir = e)
  rm(e)
  invisible(gc(verbose = FALSE))
  cat("Saved TF assignments: ", nrow(m), " amplicon cells; ",
      sum(observed), " pass the historical dominance threshold.\n", sep = "")

  use_pagoda <- any(panels %in% c("7E", "S12E", "S13C"))
  emb <- NULL
  clusters <- NULL
  if (use_pagoda) {
    p2 <- readRDS(paths[["pagoda"]])
    emb <- p2$embeddings$PCA$UMAP
    clusters <- p2$clusters$PCA$leiden
    if (is.null(dim(emb)) || ncol(emb) != 2L ||
        anyDuplicated(rownames(emb)) || is.null(names(clusters)) ||
        anyNA(match(rownames(emb), names(clusters)))) {
      stop("Invalid or unaligned Pagoda2 UMAP and cluster IDs.", call. = FALSE)
    }
    clusters <- as.character(clusters[rownames(emb)])
    if (anyNA(clusters)) stop("NA in aligned Pagoda2 clusters.", call. = FALSE)
    cat("Pagoda2: ", nrow(emb), " UMAP cells; ",
        sum(rownames(emb) %in% m$cell), " with amplicon records; ",
        sum(!rownames(emb) %in% m$cell), " without amplicon records.\n", sep = "")
    rm(p2)
    invisible(gc(verbose = FALSE))
  }
  dna <- NULL
  if ("S13B" %in% panels) {
    dna <- utils::read.table(paths[["dna"]], header = TRUE,
                             stringsAsFactors = FALSE)
    if (!all(c("ORFID", "TF", "plasmid", "HSC") %in% names(dna)) ||
        any(!is.finite(dna$HSC)) || any(dna$HSC < 0)) {
      stop("Invalid S13B DNA frequency table.", call. = FALSE)
    }
    if (anyDuplicated(dna$TF)) {
      cat("S13B: DNA table contains ", sum(duplicated(dna$TF)),
          " repeated TF name(s); preserving the original ORF-level rows ",
          "and the historical RNA denominator.\n", sep = "")
    }
  }
  list(m = m, TFn = TFn, emb = emb, clusters = clusters, dna = dna)
}

.tf_UMAP_cells <- function(v) {
  if (is.null(v$emb)) stop("Panel requires Pagoda2 UMAP.", call. = FALSE)
  ids <- rownames(v$emb)
  idx <- match(ids, as.character(v$m$cell))
  ok <- !is.na(idx)
  assigned <- rep(FALSE, length(ids))
  assigned[ok] <- as.logical(v$m$id[idx[ok]])
  TF <- rep(NA_character_, length(ids))
  barcode <- rep(NA_character_, length(ids))
  TF[ok] <- as.character(v$m$TF[idx[ok]])
  barcode[ok] <- as.character(v$m$TFbarcode[idx[ok]])
  data.frame(cell = ids, UMAP1 = as.numeric(v$emb[, 1]),
             UMAP2 = as.numeric(v$emb[, 2]),
             cluster = as.character(v$clusters),
             has_amplicon_record = ok, assigned = assigned,
             TF = TF, TFbarcode = barcode, stringsAsFactors = FALSE)
}

.tf_source_7E <- function(v) {
  d <- .tf_UMAP_cells(v)
  assigned_tf <- ifelse(d$assigned, d$TF, NA_character_)
  counts <- sort(table(assigned_tf), decreasing = TRUE)
  if (length(counts) < 10L || names(counts)[1] != "GFP") {
    stop("Fig7E: expected GFP and at least nine assigned TFs on UMAP.",
         call. = FALSE)
  }
  major <- names(counts)[2:10]
  d$TF_display <- ifelse(!is.na(assigned_tf) & assigned_tf %in% major,
                         assigned_tf, NA_character_)
  d$display_rank <- match(d$TF_display, major)
  d$display_status <- ifelse(!d$has_amplicon_record, "no_amplicon_record",
                      ifelse(!d$assigned, "unassigned_by_threshold",
                      ifelse(is.na(d$TF_display), "assigned_other_including_GFP",
                             "major_perturbation")))
  if (anyNA(d$UMAP1) || anyNA(d$UMAP2)) {
    stop("Fig7E: missing UMAP coordinates.", call. = FALSE)
  }
  cat("Fig7E displayed TF rank 2-10 (after GFP): ",
      paste(major, collapse = ", "), "; ",
      sum(!is.na(d$TF_display)), " highlighted UMAP cells.\n", sep = "")
  d
}

.tf_figure_7E <- function(dat) {
  major <- unique(dat$TF_display[order(dat$display_rank, na.last = TRUE)])
  major <- major[!is.na(major)]
  if (length(major) != 9L) stop("Fig7E: expected nine major TFs.", call. = FALSE)
  colors <- stats::setNames(RColorBrewer::brewer.pal(9, "Set1"), major)
  sel <- dat[!is.na(dat$TF_display), , drop = FALSE]
  sel$TF_display <- factor(sel$TF_display, levels = major)
  ggplot2::ggplot() +
    ggplot2::geom_point(data = dat,
      ggplot2::aes(x = UMAP1, y = UMAP2),
      color = "grey85", size = 0.32) +
    ggplot2::geom_point(data = sel,
      ggplot2::aes(x = UMAP1, y = UMAP2, color = TF_display), size = 0.75) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, color = "Perturbation") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.position = "right") +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 3)))
}

.tf_source_S12E <- function(v) {
  d <- .tf_UMAP_cells(v)
  if (!setequal(unique(d$cluster), as.character(1:10))) {
    stop("FigS12E: expected ten stored clusters 1-10.", call. = FALSE)
  }
  d
}

.tf_figure_S12E <- function(dat) {
  dat$cluster <- factor(dat$cluster, levels = as.character(1:10))
  centers <- stats::aggregate(cbind(UMAP1, UMAP2) ~ cluster, dat,
                              FUN = stats::median)
  pal <- stats::setNames(grDevices::rainbow(10), as.character(1:10))
  ggplot2::ggplot(dat, ggplot2::aes(UMAP1, UMAP2, color = cluster)) +
    ggplot2::geom_point(size = 0.32, alpha = 0.65, shape = 1) +
    ggplot2::geom_label(data = centers,
      ggplot2::aes(x = UMAP1, y = UMAP2, label = cluster),
      color = "black", fill = "white", label.size = 0, size = 3.5,
      inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, color = "Cluster") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::guides(color = "none")
}

.tf_source_S13A <- function(v) {
  m <- v$m
  a <- as.logical(m$id)
  data.frame(cell = as.character(m$cell),
             TFbarcode = as.character(m$TFbarcode), TF = as.character(m$TF),
             dominant_TF_barcode_UMIs = as.numeric(m$dominant),
             total_TF_barcode_UMIs = as.numeric(m$total),
             dominant_fraction = ifelse(m$total > 0, m$dominant / m$total, NA_real_),
             TF_assigned = a, dominant_min_exclusive = 10L,
             dominant_fraction_min_exclusive = 0.90,
             stringsAsFactors = FALSE)
}

.tf_figure_S13A <- function(dat) {
  # The *historical* plot and approved artwork have no dashed threshold lines.
  # Caption must not claim that such lines are displayed. Thresholds in Source Data.
  dat$TF_assigned <- factor(dat$TF_assigned, levels = c(FALSE, TRUE))
  ggplot2::ggplot(dat, ggplot2::aes(total_TF_barcode_UMIs,
                                   dominant_TF_barcode_UMIs,
                                   color = TF_assigned)) +
    ggplot2::geom_point(size = 0.45, alpha = 0.75) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "black")) +
    ggplot2::labs(title = "TF barcode read counts", x = "total",
                  y = "dominant", color = "id") +
    ggplot2::theme_bw(base_size = 10)
}

.tf_source_S13B <- function(v) {
  m <- v$m
  dna <- v$dna
  if (is.null(dna)) stop("S13B: DNA checkpoint not loaded.", call. = FALSE)
  # Preserve 4_tf_invitro_processing.R line 167 exactly:
  # rna.tf[dna.tf$TF] / sum(rna.tf[dna.tf$TF], na.rm=TRUE).
  # TF names can occur on multiple ORF rows; denominator includes repeats.
  rna.tf <- table(as.character(m$TF[as.logical(m$id)]))
  rna_match <- as.numeric(rna.tf[as.character(dna$TF)])
  historical_denominator <- sum(rna_match, na.rm = TRUE)
  all_assigned_denominator <- sum(as.numeric(rna.tf))
  if (historical_denominator <= 0 || sum(dna$HSC) <= 0) {
    stop("S13B: zero frequency denominator.", call. = FALSE)
  }
  d <- data.frame(ORFID = as.character(dna$ORFID),
    TF = as.character(dna$TF), DNA_HSC_reads = as.numeric(dna$HSC),
    DNA_plasmid_reads = as.numeric(dna$plasmid),
    RNA_assigned_cell_count = rna_match,
    DNA_frequency_HSC = dna$HSC / sum(dna$HSC),
    RNA_frequency_historical = rna_match / historical_denominator,
    RNA_frequency_all_assigned_audit = rna_match / all_assigned_denominator,
    DNA_HSC_denominator = sum(dna$HSC),
    RNA_historical_denominator = historical_denominator,
    RNA_all_assigned_denominator = all_assigned_denominator,
    TF_has_multiple_DNA_rows = duplicated(dna$TF) | duplicated(dna$TF, fromLast = TRUE),
    stringsAsFactors = FALSE)
  d$plot_include <- is.finite(d$DNA_frequency_HSC) & d$DNA_frequency_HSC > 0 &
    is.finite(d$RNA_frequency_historical) & d$RNA_frequency_historical > 0
  cat("FigS13B: ", nrow(d), " DNA ORF rows; ", sum(d$plot_include),
      " plottable RNA/DNA TF frequencies; historical RNA denominator=",
      historical_denominator, ", all-assigned denominator=",
      all_assigned_denominator, ".\n", sep = "")
  d
}

.tf_figure_S13B <- function(dat) {
  d <- dat[dat$plot_include, , drop = FALSE]
  if (!nrow(d)) stop("S13B: no log-scale-plottable TFs.", call. = FALSE)
  # Train the DNA x-axis using all positive ORF frequencies.
  # Continue plotting only complete RNA/DNA observations.
  x_all <- dat$DNA_frequency_HSC
  x_all <- x_all[is.finite(x_all) & x_all > 0]

  if (!length(x_all))
    stop("S13B: no positive DNA frequencies.")

  x_min <- min(x_all)
  ggplot2::ggplot(d, ggplot2::aes(DNA_frequency_HSC,
                                 RNA_frequency_historical, label = TF)) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggrepel::geom_text_repel(size = 2.5, seed = 1, max.overlaps = Inf) +
    ggplot2::expand_limits(x = x_min) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(title = "in vitro TF frequencies", x = "frequency, DNA",
                  y = "frequency, RNA") +
    ggplot2::theme_bw(base_size = 10)
}

.tf_source_S13C <- function(v) {
  if (is.null(v$emb)) stop("S13C: requires UMAP checkpoint for auditing.", call. = FALSE)
  tf <- names(v$TFn)
  cells <- as.integer(v$TFn)
  if (!length(tf) || length(tf) != length(cells) || anyNA(cells)) {
    stop("S13C: invalid saved TFn count table.", call. = FALSE)
  }
  m <- v$m
  selected <- as.logical(m$id) & as.character(m$cell) %in% rownames(v$emb)
  embedded_counts <- table(as.character(m$TF[selected]))
  all_counts <- table(as.character(m$TF[as.logical(m$id)]))
  d <- data.frame(TF = tf, assigned_cells_in_p2_saved_TFn = cells,
                  assigned_cells_in_UMAP_audit = as.integer(embedded_counts[tf]),
                  assigned_cells_in_all_amplicons = as.integer(all_counts[tf]),
                  stringsAsFactors = FALSE)
  d$assigned_cells_in_UMAP_audit[is.na(d$assigned_cells_in_UMAP_audit)] <- 0L
  d$assigned_cells_in_all_amplicons[is.na(d$assigned_cells_in_all_amplicons)] <- 0L
  d$saved_vs_UMAP_count_equal <- d$assigned_cells_in_p2_saved_TFn ==
    d$assigned_cells_in_UMAP_audit
  d$plot_include <- d$assigned_cells_in_p2_saved_TFn > 2L
  d$historical_rank <- seq_len(nrow(d))
  cat("FigS13C: ", nrow(d), " TF names in saved TFn; ",
      sum(d$plot_include), " with >2 cells; ",
      sum(!d$saved_vs_UMAP_count_equal),
      " TFs have saved-vs-UMAP count differences (audit column).\n", sep = "")
  d
}

.tf_figure_S13C <- function(dat) {
  d <- dat[dat$plot_include, , drop = FALSE]
  if (!nrow(d)) stop("S13C: no TFs with >2 cells.", call. = FALSE)
  d <- d[order(-d$assigned_cells_in_p2_saved_TFn, d$historical_rank), , drop = FALSE]
  d$position <- seq_len(nrow(d))
  ggplot2::ggplot(d) +
    ggplot2::geom_rect(ggplot2::aes(
      xmin = position - .38, xmax = position + .38,
      ymin = 1, ymax = assigned_cells_in_p2_saved_TFn),
      fill = "grey75", color = "black", linewidth = .22) +
    ggplot2::scale_x_continuous(breaks = d$position, labels = d$TF,
                                expand = ggplot2::expansion(add = .6)) +
    ggplot2::scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500)) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "cells") +
    ggplot2::theme_classic(base_size = 10)
}

.tf_write_table <- function(dat, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(dat, con, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
}

.tf_write_figure <- function(plot, path, width = 5.5, height = 4.4) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}


# ===== Figure S12F: restored from saved Pagoda2 checkpoint =====
# Original: TF_overexpression/scripts/4_tf_invitro_processing.R, lines 30-74.
# This stage does NOT recompute UMAP, clustering, or differential genes.
# Historic cell sampling lacked a seed; here seed=712006 is a *new* recorded
# deterministic choice. The original sampled 1,252-cell identities are unknown.
.tf_s12f_checkpoint_genes <- function() {
  # Verified selection/order from original checkpoint, 67 marker genes.
  # This is a validation reference only; selection comes from diffgenes.
  c("SPAG17", "MAF", "HLF", "AVP", "MEIS1", "MCM4", "GINS2", "MSH6",
    "UNG", "MCM3", "HELLS", "MCM6", "PCNA", "HIST1H1B", "HIST1H4C",
    "HIST1H1D", "HIST1H3D", "HIST1H2AH", "HIST1H1E", "HIST1H1C",
    "HIST1H3G", "CENPF", "HMMR", "ASPM", "DLGAP5", "TUBB4B", "TOP2A",
    "CKS2", "CCNB1", "MKI67", "PTTG1", "IGLL1", "HMGN2", "PTMS",
    "SNHG3", "ILF3-DT", "DNTT", "HLA-DRB1", "XIST", "RNASET2",
    "FAM118A", "AP002387.2", "ZNF667-AS1", "SCN3A", "PRSS2", "LST1",
    "LY86", "MS4A6A", "SAMHD1", "LYZ", "NKG7", "TYROBP", "CLEC10A",
    "IRF7", "LGALS1", "SRGN", "LMO4", "HDC", "PRG2", "S100A6",
    "CPA3", "ANXA1", "ALOX5AP", "FCER1A", "PPM1F", "ATP8B4",
    "AC004687.1")
}

.tf_s12f_markers <- function(obj) {
  dg <- obj$diffgenes$PCA[[1]]
  grp <- obj$clusters$PCA[[1]]
  counts <- obj$counts
  if (length(dg) != 10L || length(grp) != 20459L ||
      is.null(names(grp)) || anyDuplicated(names(grp)) ||
      !all(names(grp) %in% rownames(counts))) {
    stop("S12F: unexpected Pagoda2 checkpoint gene/cluster structure.",
         call. = FALSE)
  }
  marker_genes <- lapply(dg, function(de) {
    if (!all(c("Z", "M") %in% colnames(de))) {
      stop("S12F: missing M or Z from saved differential-gene table.")
    }
    de <- de[!is.na(de$Z) & !is.na(de$M) &
             abs(de$Z) > 8 & de$M > 1, , drop = FALSE]
    if (!nrow(de)) return(character())
    rownames(de)[head(order(-de$Z), min(8L, nrow(de)))]
  })
  hsc <- c("HLF", "PROCR", "HOPX", "AVP", "MEIS1")
  bprog <- c("CD79A", "MS4A1", "VPREB1/3", "DNTT", "RAG1/2")
  dc <- c("HLA-DRA", "CD74", "LST1", "FCER1A", "CLEC10A", "IRF7", "GZMB")
  genes <- unique(c(unlist(marker_genes, use.names = FALSE), hsc, bprog, dc))
  genes <- genes[!is.na(genes) & genes %in% colnames(counts)]
  y <- do.call(cbind, lapply(dg, function(de) {
    as.numeric(de[genes, "M"])
  }))
  rownames(y) <- genes
  max_cl <- apply(y, 1L, function(x) which.max(x)[1L])
  remove <- c(
    hsc[!max_cl[hsc] %in% c(1, 6)],
    bprog[!max_cl[bprog] %in% 5],
    dc[!max_cl[dc] %in% 7:9]
  )
  selected <- names(sort(max_cl[!is.na(max_cl) &
                                  !names(max_cl) %in% remove]))
  expected <- .tf_s12f_checkpoint_genes()
  if (length(selected) != 67L || !identical(selected, expected)) {
    stop("S12F: saved 67-gene marker selection/order differs from historical ",
         "validated checkpoint; inspect selected vs reference first. Missing: ",
         paste(setdiff(expected, selected), collapse = ", "),
         "; extra: ", paste(setdiff(selected, expected), collapse = ", "),
         call. = FALSE)
  }
  # Export the per-cluster DE statistics underlying marker selection.
  stat_list <- lapply(seq_along(dg), function(k) {
    de <- dg[[k]]
    data.frame(
      gene = selected, marker_order = seq_along(selected),
      cluster_test_index = as.integer(k),
      M = as.numeric(de[selected, "M"]),
      Z = as.numeric(de[selected, "Z"]),
      gene_max_M_cluster_index = as.integer(max_cl[selected]),
      stringsAsFactors = FALSE
    )
  })
  stats <- do.call(rbind, stat_list)
  stopifnot(nrow(stats) == 670L, !anyDuplicated(paste(stats$gene,
                     stats$cluster_test_index, sep = ":")))
  list(selected = selected, stats = stats, groups = grp)
}

.tf_s12f_sampling <- function(groups, seed = 712006L) {
  # Preserve the *historical formula and group iteration order*.
  # sample() coerces the real-valued requested size to integer.
  counts <- table(groups)
  minimum <- min(counts)
  set.seed(seed)
  sampled <- unlist(lapply(unique(groups), function(x) {
    group <- groups[groups %in% x]
    sample(group, minimum + (length(group) - minimum)^0.6)
  }))
  if (!length(sampled) || is.null(names(sampled)) ||
      anyNA(names(sampled)) || anyDuplicated(names(sampled)) ||
      length(sampled) != 1252L) {
    stop("S12F: unexpected sampled cell IDs/count. ",
         "Do not proceed without verifying Pagoda2 cluster names.", call. = FALSE)
  }
  sizes <- c(`1`=222L, `2`=176L, `3`=137L, `4`=140L, `5`=97L,
             `6`=177L, `7`=116L, `8`=76L, `9`=41L, `10`=70L)
  actual <- table(factor(as.character(sampled), levels = names(sizes)))
  if (!identical(as.integer(actual), as.integer(sizes))) {
    stop("S12F: sampling formula no longer yields the inspected per-cluster sizes.",
         call. = FALSE)
  }
  sampled
}

.tf_s12f_write_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  tryCatch(utils::write.table(x, con, sep = "\t", quote = FALSE,
                              row.names = FALSE, na = "NA"),
           finally = close(con))
  back <- utils::read.delim(path, check.names = FALSE,
                            stringsAsFactors = FALSE)
  if (!identical(dim(back), dim(x)) || !identical(names(back), names(x))) {
    stop("S12F Source Data roundtrip failed: ", path, call. = FALSE)
  }
  invisible(path)
}

.tf_s12f_pdf <- function(obj, genes, sampled, path) {
  grDevices::pdf(path, height = 6, width = 3, useDingbats = FALSE)
  tryCatch({
    obj$plotGeneHeatmap(genes = genes, groups = sampled,
                        useRaster = FALSE, box = TRUE)
  }, finally = grDevices::dev.off())
  if (!file.exists(path) || is.na(file.info(path)$size) ||
      file.info(path)$size < 1000L) {
    stop("S12F: PDF missing or empty after native Pagoda2 rendering.",
         call. = FALSE)
  }
  invisible(path)
}

reproduce_FigS12F <- function(
  root = getwd(),
  output_root = file.path(root, "TF_overexpression/publication"),
  obj = NULL,
  seed = 712006L
) {
  if (!requireNamespace("pagoda2", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE)) {
    stop("S12F requires available pagoda2 and Matrix; no installation attempted.",
         call. = FALSE)
  }
  if (getNamespaceVersion("Matrix") < package_version("1.6-4") ||
      normalizePath(getNamespaceInfo("Matrix", "path")) !=
      normalizePath(file.path(Sys.getenv("CH_R_LIB", unset = "/packages"), "Matrix"))) {
    stop("S12F: load historical Matrix 1.6-4 from /packages in a fresh R session.",
         call. = FALSE)
  }
  root <- normalizePath(root, mustWork = TRUE)
  input <- file.path(root, "TF_overexpression/data/processed/p2_tfinvitro.rds")
  if (is.null(obj)) {
    if (!file.exists(input)) stop("S12F: missing saved Pagoda2 checkpoint: ", input)
    cat("Loading saved Pagoda2 checkpoint: ", input, "\n", sep = "")
    obj <- readRDS(input)
  } else {
    cat("S12F: using supplied Pagoda2 object; no checkpoint reload.\n")
  }
  if (!is.function(obj$plotGeneHeatmap)) {
    stop("S12F: checkpoint has no plotGeneHeatmap method.", call. = FALSE)
  }
  marker <- .tf_s12f_markers(obj)
  # Sampled cells are deterministic for R 4.2.3 under the documented seed.
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  sampled <- .tf_s12f_sampling(marker$groups, seed = as.integer(seed))
  cells <- names(sampled)
  if (!all(cells %in% rownames(obj$counts))) {
    stop("S12F: sampled cell IDs absent from saved expression matrix.")
  }
  cells_df <- data.frame(
    cell_id = cells, cluster = as.character(sampled),
    input_sampling_order = seq_along(cells),
    sampling_seed = as.integer(seed),
    stringsAsFactors = FALSE)
  values <- as.matrix(obj$counts[cells, marker$selected, drop = FALSE])
  if (!identical(dim(values), c(1252L, 67L))) {
    stop("S12F: unexpected sampled expression matrix dimensions.")
  }
  expr_df <- data.frame(
    cell_id = rep(cells, times = length(marker$selected)),
    cluster = rep(as.character(sampled), times = length(marker$selected)),
    input_sampling_order = rep(seq_along(cells), times = length(marker$selected)),
    gene = rep(marker$selected, each = length(cells)),
    marker_order = rep(seq_along(marker$selected), each = length(cells)),
    saved_pagoda2_value = as.numeric(values),
    stringsAsFactors = FALSE)
  stopifnot(nrow(expr_df) == 83884L, nrow(marker$stats) == 670L,
            nrow(cells_df) == 1252L)

  output_root <- normalizePath(output_root, mustWork = FALSE)
  historical <- normalizePath(file.path(root, "TF_overexpression/figures"),
                              mustWork = FALSE)
  if (identical(output_root, historical) ||
      startsWith(output_root, paste0(historical, "/"))) {
    stop("S12F: refusing to write into original figures directory.")
  }
  specs <- c(
    "source_data/FigS12F_marker_cluster_statistics.tsv.gz",
    "source_data/FigS12F_sampled_cells.tsv.gz",
    "source_data/FigS12F_sampled_cell_expression.tsv.gz",
    "figures/FigS12F_marker_heatmap_REVIEW.pdf"
  )
  targets <- file.path(output_root, specs)
  if (any(file.exists(targets))) {
    stop("S12F: existing approved output(s) protected:\n",
         paste(targets[file.exists(targets)], collapse = "\n"),
         "\nUse a fresh output_root for an independent rerun.", call. = FALSE)
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  staged <- tempfile(".S12F_stage_", tmpdir = output_root)
  dir.create(staged)
  on.exit(unlink(staged, recursive = TRUE), add = TRUE)
  for (d in c("source_data", "figures"))
    dir.create(file.path(staged, d), recursive = TRUE, showWarnings = FALSE)
  staged_paths <- file.path(staged, specs)
  .tf_s12f_write_gz(marker$stats, staged_paths[1L])
  .tf_s12f_write_gz(cells_df, staged_paths[2L])
  .tf_s12f_write_gz(expr_df, staged_paths[3L])
  .tf_s12f_pdf(obj, marker$selected, sampled, staged_paths[4L])
  for (i in seq_along(targets)) {
    dir.create(dirname(targets[i]), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(targets[i]) || !file.rename(staged_paths[i], targets[i])) {
      stop("S12F: could not finalize output: ", targets[i], call. = FALSE)
    }
  }
  cat("S12F: selected=", length(marker$selected),
      " historical marker genes; sampled=", length(sampled),
      " cells across 10 clusters; numeric Source Data=670 marker-cluster, ",
      "1252 cells, 83884 cell-gene rows.\n", sep = "")
  cat("S12F: deterministic seed=", seed,
      "; historical sampled cell identities unknown.\n", sep = "")
  cat("S12F: source-data roundtrip PASS; historical visual comparison REVIEW.\n")
  invisible(list(panel = "S12F", marker_genes = marker$selected,
                 n_cells = length(sampled), output_root = output_root,
                 source_data = targets[1:3], pdf = targets[4], seed = seed))
}
# ===== End S12F checkpoint reproduction extension =====


# ===== TF pathway enrichment from archived in-vitro DE checkpoint =====
# Original: scripts/6_tf_invitro_de.R lines ~124 and 190-230.
# Source DE genes and fgsea inputs are checkpointed; no DESeq2 refitting.
# IMPORTANT: The historical Fig7F matrix combines fgsea results by ROW INDEX,
# then assigns first-TF pathway names. This may mislabel columns if fgsea
# results differ in order. The audit below tests that assumption rather than
# silently using positional alignment as if it were pathway-keyed.
# Per-TF RDS checkpoints allow restart without recomputing completed TFs.

.tf_pathway_inputs <- function(root = getwd()) {
  root <- normalizePath(root, mustWork = TRUE)
  paths <- c(
    DE = file.path(root, "TF_overexpression/data/processed/TF_invitro_DE.RData"),
    genesets = file.path(root, "RNA_analysis/data/genesets/processed/gene_sets_use.rds"),
    TF_counts = file.path(root, "TF_overexpression/publication/source_data/FigS13C_TF_cells_per_perturbation.tsv.gz")
  )
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("TF pathway checkpoints missing: ", paste(missing, collapse = "; "),
         call. = FALSE)
  }
  pkgs <- c("fgsea", "Matrix")
  pk_ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!pk_ok)) {
    stop("Missing pathway-enrichment R packages: ",
         paste(pkgs[!pk_ok], collapse = ", "),
         "; no installation attempted.", call. = FALSE)
  }
  if (getNamespaceVersion("Matrix") < package_version("1.6-4") ||
      normalizePath(getNamespaceInfo("Matrix", "path")) !=
        normalizePath(file.path(Sys.getenv("CH_R_LIB", unset = "/packages"), "Matrix"))) {
    stop("Load historical Matrix 1.6-4 from /packages in a fresh R session before fgsea.",
         call. = FALSE)
  }
  list(paths = paths, md5 = unname(tools::md5sum(paths)),
       fgsea_version = as.character(utils::packageVersion("fgsea")))
}

.tf_pathway_prepare_de <- function(paths) {
  e <- new.env(parent = baseenv())
  load(paths[["DE"]], envir = e)
  if (!all(c("logFC.tf", "pv.tf") %in% ls(e))) {
    stop("Missing logFC.tf or pv.tf from saved TF_invitro_DE.RData.", call. = FALSE)
  }
  lfc <- e$logFC.tf
  pv <- e$pv.tf
  if (!identical(dim(lfc), dim(pv)) ||
      !identical(rownames(lfc), rownames(pv)) ||
      !identical(colnames(lfc), colnames(pv)) ||
      is.null(rownames(lfc)) || is.null(colnames(lfc)) ||
      anyDuplicated(rownames(lfc)) || anyDuplicated(colnames(lfc))) {
    stop("TF fgsea: the saved logFC and P matrices are not uniquely aligned.",
         call. = FALSE)
  }
  if (any(!is.finite(lfc)) || any(!is.finite(pv)) ||
      any(pv <= 0 | pv > 1)) {
    stop("TF fgsea: nonfinite/zero/out-of-range DE statistics; preserve historical formula and inspect first.",
         call. = FALSE)
  }
  if (!identical(dim(lfc), c(11176L, 38L))) {
    stop("TF fgsea: unexpected saved DE dimensions; inspect input revision first.",
         call. = FALSE)
  }
  tf_counts <- utils::read.delim(paths[["TF_counts"]], check.names = FALSE,
                                stringsAsFactors = FALSE)
  required <- c("TF", "assigned_cells_in_p2_saved_TFn")
  if (!all(required %in% names(tf_counts)) || anyDuplicated(tf_counts$TF)) {
    stop("TF fgsea: missing or duplicated source TF-cell counts.", call. = FALSE)
  }
  selected <- intersect(colnames(lfc),
    tf_counts$TF[!is.na(tf_counts$assigned_cells_in_p2_saved_TFn) &
                 tf_counts$assigned_cells_in_p2_saved_TFn >= 4])
  if (length(selected) != 25L) {
    stop("TF fgsea: expected 25 eligible TF comparisons, got ", length(selected),
         "; inspect eligibility threshold and saved cell counts first.", call. = FALSE)
  }
  gene_sets <- readRDS(paths[["genesets"]])
  if (!is.list(gene_sets) || is.null(names(gene_sets)) ||
      anyNA(names(gene_sets)) || any(!nzchar(names(gene_sets))) ||
      anyDuplicated(names(gene_sets))) {
    stop("TF fgsea: expected uniquely named list of gene sets.", call. = FALSE)
  }
  cat("TF fgsea inputs: ", nrow(lfc), " genes, ", ncol(lfc),
      " DE comparisons, ", length(selected), " eligible TFs, ",
      length(gene_sets), " gene sets.\n", sep = "")
  cat("Selected TF order: ", paste(selected, collapse = ", "), "\n", sep = "")
  list(lfc = lfc, pvalue = pv, tfs = selected, pathways = gene_sets)
}

.tf_pathway_save_rds_atomic <- function(object, target) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".fgsea_tmp_", tmpdir = dirname(target), fileext = ".rds")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (file.exists(target) || !file.rename(tmp, target)) {
    stop("Cannot commit TF fgsea checkpoint: ", target,
         "; existing checkpoints remain untouched.", call. = FALSE)
  }
  invisible(target)
}

# Use max_tf = 2L for a small alignment pilot (saved for reuse), then omit
# max_tf to compute the remaining TFs. source_data/ and figures/ are not touched.
tf_pathway_prepare <- function(
  root = getwd(),
  output_root = file.path(root, "TF_overexpression/publication"),
  max_tf = 2L,
  seed = 713100L
) {
  if (!is.null(max_tf) &&
      (length(max_tf) != 1L || is.na(max_tf) || max_tf < 1L ||
       max_tf != as.integer(max_tf))) {
    stop("max_tf must be NULL (all eligible TFs) or a positive integer.",
         call. = FALSE)
  }
  input <- .tf_pathway_inputs(root)
  dat <- .tf_pathway_prepare_de(input$paths)
  if (is.null(max_tf)) max_tf <- length(dat$tfs)
  selected <- head(dat$tfs, as.integer(max_tf))
  fgsea_dir <- file.path(output_root, "intermediate/TF_fgsea_pathways")
  dir.create(fgsea_dir, recursive = TRUE, showWarnings = FALSE)
  # Preserve the caller's random state when a pilot or full run is executed.
  old_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_rng) saved_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_rng) assign(".Random.seed", saved_rng, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  res_list <- vector("list", length(selected))
  names(res_list) <- selected
  for (tf in selected) {
    i <- match(tf, dat$tfs)
    tf_seed <- as.integer(seed + i)
    target <- file.path(fgsea_dir, paste0("TF_", tf, ".rds"))
    if (file.exists(target)) {
      item <- readRDS(target)
      if (!is.list(item) || !identical(item$TF, tf) ||
          !identical(item$seed, tf_seed) ||
          !identical(item$input_md5, input$md5) ||
          !identical(item$fgsea_version, input$fgsea_version) ||
          !is.data.frame(item$results)) {
        stop("Saved TF fgsea checkpoint metadata differs: ", target,
             ". Review rather than overwriting it.", call. = FALSE)
      }
      cat("REUSE fgsea TF=", tf, " (", nrow(item$results), " pathways)\n", sep = "")
    } else {
      cat("COMPUTE fgsea TF=", tf, " (", i, "/", length(dat$tfs), ")\n", sep = "")
      # Historical statistic and fgsea call: original 6_tf_invitro_de.R
      sc <- sign(dat$lfc[, tf]) * -log(dat$pvalue[, tf])
      names(sc) <- rownames(dat$lfc)
      if (any(!is.finite(sc))) stop("TF fgsea input ranks nonfinite: ", tf)
      set.seed(tf_seed)
      result <- fgsea::fgsea(pathways = dat$pathways, stats = sc,
                             minSize = 15, maxSize = 1000)
      if (!all(c("pathway", "NES", "pval", "padj", "size") %in% names(result)) ||
          anyDuplicated(as.character(result$pathway))) {
        stop("TF fgsea returned invalid or duplicated pathway records: ", tf)
      }
      item <- list(TF = tf, seed = tf_seed,
                   fgsea_version = input$fgsea_version,
                   input_md5 = input$md5,
                   rank_definition = "sign(DE log2FC) * -ln(nominal DE P)",
                   results = result)
      .tf_pathway_save_rds_atomic(item, target)
      cat("SAVED fgsea TF=", tf, " (", nrow(result), " pathways)\n", sep = "")
    }
    res_list[[tf]] <- item$results
  }
  # Important audit: historical 7F uses positional cbind of independent
  # pathway-enrichment tables. Distinct fgsea tables may have different orders.
  reference <- as.character(res_list[[1L]]$pathway)
  audit <- do.call(rbind, lapply(names(res_list), function(tf) {
    pathways <- as.character(res_list[[tf]]$pathway)
    npos <- if (length(pathways) == length(reference)) {
      sum(pathways == reference)
    } else {
      NA_integer_
    }
    data.frame(TF = tf, n_pathways = length(pathways),
               n_shared_with_first = sum(pathways %in% reference),
               matched_at_same_position = npos,
               exact_pathway_order = identical(reference, pathways),
               stringsAsFactors = FALSE)
  }))
  audit_path <- file.path(fgsea_dir, "pathway_order_audit.tsv")
  utils::write.table(audit, audit_path, sep = "\t", row.names = FALSE,
                     quote = FALSE, na = "NA")
  cat("\n=== Historical Fig7F positional pathway alignment audit ===\n")
  print(audit, row.names = FALSE)
  cat("Audit: ", audit_path, "\n", sep = "")
  if (any(!audit$exact_pathway_order)) {
    cat("REVIEW: fgsea pathway row orders differ across TFs. The historical ",
        "7F/S13E positional cbind may assign incorrect pathway labels. ",
        "Do not approve the published heatmap mapping without investigating.\n", sep = "")
  } else {
    cat("PASS: all checked fgsea pathway rows share identical order.\n")
  }
  invisible(list(n_completed = length(selected), n_total = length(dat$tfs),
                 tf_order = selected, audit = audit,
                 checkpoint_dir = fgsea_dir,
                 input_md5 = input$md5,
                 fgsea_version = input$fgsea_version))
}
# ===== End TF pathway preparation extension =====


# ==== Reproduce Fig7F, FigS13E and FigS13F from the 25 saved TF fgsea checkpoints ====
# Historical: scripts/6_tf_invitro_de.R (Fig7F/S13E) and
# scripts/8_revision.R (S13F). No fgsea or DESeq2 is rerun by this exporter.
# The original Figure 7F pathways were aggregated by row position; the
# checkpoint audit must certify exact pathway ordering in all 25 TFs.

.tf_pathway_figures_load <- function(root, output_root) {
  root <- normalizePath(root, mustWork = TRUE)
  check <- .tf_pathway_inputs(root) # current input hashes and fgsea version
  directory <- file.path(output_root, "intermediate/TF_fgsea_pathways")
  audit_file <- file.path(directory, "pathway_order_audit.tsv")
  if (!file.exists(audit_file)) {
    stop("Missing 25-TF pathway-order audit; run tf_pathway_prepare(max_tf=NULL) first.", call. = FALSE)
  }
  audit <- utils::read.delim(audit_file, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("TF", "n_pathways", "matched_at_same_position", "exact_pathway_order")
  if (!all(required %in% names(audit)) || nrow(audit) != 25L ||
      anyDuplicated(audit$TF) || anyNA(audit$exact_pathway_order) ||
      !all(audit$exact_pathway_order) || length(unique(audit$n_pathways)) != 1L ||
      any(audit$matched_at_same_position != audit$n_pathways)) {
    stop("Full 25-TF audit is not a valid exact alignment PASS.", call. = FALSE)
  }
  tfs <- as.character(audit$TF)
  checkpoint <- vector("list", length(tfs)); names(checkpoint) <- tfs
  for (tf in tfs) {
    path <- file.path(directory, paste0("TF_", tf, ".rds"))
    if (!file.exists(path)) stop("Missing checkpoint: ", path, call. = FALSE)
    item <- readRDS(path)
    if (!identical(item$TF, tf) ||
        !identical(item$input_md5, check$md5) ||
        !identical(item$fgsea_version, check$fgsea_version) ||
        !is.data.frame(item$results)) {
      stop("TF fgsea checkpoint metadata mismatch: ", path, call. = FALSE)
    }
    needed <- c("pathway", "pval", "padj", "NES", "ES", "size")
    if (!all(needed %in% names(item$results))) {
      stop("TF fgsea checkpoint is missing required statistical columns: ", tf,
           call. = FALSE)
    }
    checkpoint[[tf]] <- as.data.frame(item$results)
  }
  ref <- as.character(checkpoint[[1L]]$pathway)
  if (!length(ref) || anyNA(ref) || anyDuplicated(ref) ||
      !all(vapply(checkpoint, function(x) identical(as.character(x$pathway), ref),
                  logical(1)))) {
    stop("Checkpoint pathway order inconsistent despite saved audit; stop.", call. = FALSE)
  }
  cat("Pathway export: ", length(checkpoint), " TFs x ", length(ref),
      " identical pathway IDs; no fgsea recomputation.\n", sep = "")
  list(tf_order = tfs, results = checkpoint, pathways = ref,
       input_md5 = check$md5, fgsea_version = check$fgsea_version)
}

.tf_pathway_figures_scores <- function(obj) {
  v <- obj$results
  mat <- do.call(cbind, lapply(v, function(d) {
    pv <- as.numeric(d$pval)
    nes <- as.numeric(d$NES)
    if (any(!is.na(pv) & (pv <= 0 | pv > 1))) {
      stop("fgsea contains zero/out-of-range pathway P values; " ,
           "cannot take historical -ln(P) without review.", call. = FALSE)
    }
    -log(pv) * sign(nes)
  }))
  dimnames(mat) <- list(obj$pathways, obj$tf_order)
  # Original 6_tf_invitro_de.R uses: tf.enrich.go[is.na(tf.enrich.go)] <- 0
  missing <- sum(!is.finite(mat))
  if (missing) {
    cat("Historical 7F display rule: ", missing,
        " undefined pathway scores replaced by 0 in heatmap only; ",
        "raw fgsea values remain NA in Source Data.\n", sep = "")
  }
  mat[is.na(mat)] <- 0
  if (any(!is.finite(mat))) stop("Nonfinite heatmap scores remain.", call. = FALSE)
  selected <- apply(abs(mat), 1L, max) > 9 # historical Fig7F rule
  if (sum(selected) < 5L) {
    stop("Historical >9 pathway threshold retained fewer than 5 pathways; " ,
         "cannot form the five pathway clusters.", call. = FALSE)
  }
  m <- mat[selected, , drop = FALSE]
  # Historical: as.dist(1 - cor(t(tf.enrich.go[ind1,])))
  cor_m <- stats::cor(t(m))
  if (any(!is.finite(cor_m))) {
    stop("Undefined correlation in historical Fig7F clustering.", call. = FALSE)
  }
  dist_m <- stats::as.dist(1 - cor_m)
  tree <- stats::hclust(dist_m)
  cl <- factor(stats::cutree(tree, k = 5L))
  names(cl) <- rownames(m)
  col.tf.go <- RColorBrewer::brewer.pal(length(unique(cl)), "Set2")
  names(col.tf.go) <- unique(as.character(cl))
  cat("Fig7F: ", sum(selected), " pathways exceed historical |signed -ln(P)| > 9; ",
      "five historical correlation clusters.\n", sep = "")
  list(full = mat, selected = selected, subset = m, cluster = cl,
       distance = dist_m, tree = tree, cluster_colors = col.tf.go)
}

.tf_pathway_export_all <- function(obj, scores) {
  long <- lapply(obj$tf_order, function(tf) {
    x <- obj$results[[tf]]
    data.frame(TF = tf, pathway = as.character(x$pathway),
               ES = as.numeric(x$ES), NES = as.numeric(x$NES),
               pvalue = as.numeric(x$pval), padj = as.numeric(x$padj),
               pathway_size = as.integer(x$size),
               signed_minus_ln_pathway_P = as.numeric(scores$full[, tf]),
               selected_Fig7F = as.logical(scores$selected),
               Fig7F_cluster = as.character(scores$cluster[as.character(x$pathway)]),
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, long); rownames(d) <- NULL
  if (nrow(d) != length(obj$tf_order) * length(obj$pathways) ||
      anyDuplicated(paste(d$TF, d$pathway, sep = "\t"))) {
    stop("Invalid TF/pathway Source Data keys or row count.", call. = FALSE)
  }
  d
}

.tf_pathway_export_7F <- function(all, scores) {
  # One row per (TF, selected pathway): export exactly the values used by the PDF.
  d <- all[all$selected_Fig7F, , drop = FALSE]
  d$matrix_clustering_distance <- "1_minus_Pearson_correlation_across_TFs"
  d$pathway_selection <- "max_abs_signed_minus_ln_nominal_P_gt_9"
  d
}

.tf_pathway_figure_7F <- function(scores, file) {
  ann <- data.frame(clust = scores$cluster)
  rownames(ann) <- names(scores$cluster)
  grDevices::pdf(file, width = 5, height = 5, useDingbats = FALSE)
  tryCatch({
    # Reproduce historical figure settings, including signed -ln(P),
    # fixed [-30,30] palette, and correlation distance-based pathway tree.
    p <- pheatmap::pheatmap(
      t(scores$subset), show_colnames = FALSE, main = "Pathway enrichment",
      treeheight_col = 10, treeheight_row = 20,
      clustering_distance_cols = scores$distance,
      annotation_col = ann,
      annotation_colors = list(clust = scores$cluster_colors),
      color = grDevices::colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-30, 30, length.out = 100),
      silent = TRUE)
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  }, finally = grDevices::dev.off())
  if (!file.exists(file) || file.info(file)$size < 1000) {
    stop("Fig7F PDF missing or empty after rendering.", call. = FALSE)
  }
  invisible(file)
}

.tf_pathway_pheno_groups <- function(root) {
  env <- new.env(parent = baseenv())
  load(file.path(root, "TF_overexpression/data/processed/TF_trends.RData"), envir = env)
  growth <- env$tf.trends.no25.v1m
  bias <- env$biastot
  if (!all(c("TF", "p.t", "t") %in% names(growth)) ||
      !all(c("TF", "trend") %in% names(bias))) {
    stop("TF phenotype checkpoint format changed.", call. = FALSE)
  }
  q <- stats::p.adjust(growth$p.t, method = "BH")
  list(up = unique(as.character(growth$TF[!is.na(q) & q < 0.1 & growth$t > 0])),
       down = unique(as.character(growth$TF[!is.na(q) & q < 0.1 & growth$t < 0])),
       myeloid = unique(as.character(bias$TF[as.character(bias$trend) == "myeloid"])))
}

.tf_pathway_S13F_patterns <- function() {
  # Exact four textual selections in scripts/8_revision.R; resolve via the
  # archived pathway names; stop if a pattern maps to zero or multiple paths.
  c("MYC_TARGETS_V1", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "Signaling by Receptor Tyrosine Kinases")
}

.tf_pathway_export_S13F <- function(obj, all, root) {
  ph <- .tf_pathway_pheno_groups(root)
  patterns <- .tf_pathway_S13F_patterns()
  resolved <- vapply(patterns, function(pattern) {
    hit <- grep(pattern, obj$pathways, fixed = TRUE, value = TRUE)
    if (length(hit) != 1L) {
      stop("S13F original pathway label '", pattern,
           "' maps to ", length(hit), " saved pathways: ",
           paste(head(hit, 6), collapse = "; "),
           ". Inspect rather than guess.", call. = FALSE)
    }
    hit
  }, character(1))
  out <- lapply(seq_along(patterns), function(i) {
    d <- all[all$pathway == resolved[i], , drop = FALSE]
    d$panel_position <- i
    d$historical_pathway_pattern <- patterns[i]
    d$growth_class <- ifelse(d$TF %in% ph$up, "up",
                             ifelse(d$TF %in% ph$down, "down", "stable"))
    d$lineage_class <- ifelse(d$TF %in% ph$myeloid, "myeloid", "none")
    d$significant_FDR_005 <- !is.na(d$padj) & d$padj < 0.05
    d
  })
  d <- do.call(rbind, out)
  rownames(d) <- NULL
  if (nrow(d) != 4L * length(obj$tf_order)) {
    stop("S13F expected one result per TF per four pathways.", call. = FALSE)
  }
  cat("S13F: resolved all 4 historical pathway names uniquely; ",
      nrow(d), " TF/pathway observations.\n", sep = "")
  d
}

.tf_pathway_figure_S13F <- function(d, file) {
  base <- lapply(seq_len(4L), function(i) {
    x <- d[d$panel_position == i & is.finite(d$NES) &
             is.finite(d$pvalue) & d$pvalue > 0, , drop = FALSE]
    if (!nrow(x)) stop("No plottable S13F observations for panel ", i, call. = FALSE)
    x$growth_class <- factor(x$growth_class,
                             levels = c("stable", "up", "down"))
    x$lineage_class <- factor(x$lineage_class,
                              levels = c("none", "myeloid"))
    ggplot2::ggplot(x, ggplot2::aes(x = NES, y = pvalue, label = TF)) +
      ggplot2::geom_point(ggplot2::aes(shape = growth_class,
                                    fill = lineage_class,
                                    alpha = significant_FDR_005),
                          color = "black", size = 5) +
      ggplot2::scale_shape_manual(values = c(stable = 21, up = 24, down = 23)) +
      ggplot2::scale_fill_manual(values = c(none = "grey60", myeloid = "orange")) +
      ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.2, `TRUE` = 1),
                                  guide = "none") +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::scale_y_log10() +
      ggrepel::geom_text_repel(max.overlaps = 10, box.padding = 0.5,
                               point.padding = 0.3, seed = 1, size = 3) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Normalized enrichment score", y = "P value",
                    title = unique(x$historical_pathway_pattern)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5))
  })
  grDevices::pdf(file, width = 10, height = 8, useDingbats = FALSE)
  tryCatch(print(cowplot::plot_grid(plotlist = base,
                 labels = c("A", "B", "C", "D"), ncol = 2,
                 nrow = 2, scale = 0.9)),
           finally = grDevices::dev.off())
  if (!file.exists(file) || file.info(file)$size < 1000) {
    stop("S13F PDF missing or empty.", call. = FALSE)
  }
  invisible(file)
}

.tf_pathway_figure_S13E <- function(scores, root, file) {
  helper <- file.path(root, "RNA_analysis/R/draw_go.R")
  if (!file.exists(helper)) {
    stop("S13E requires original RNA_analysis/R/draw_go.R helper; " ,
         "historical helper absent. Do not substitute another function.", call. = FALSE)
  }
  # Reuse the original draw_go function and its documented 6_tf_invitro_de.R
  # arguments. Source into a private environment with access to R attachments.
  env <- new.env(parent = .GlobalEnv)
  sys.source(helper, envir = env)
  if (!is.function(env$draw_go)) {
    stop("Historical helper does not define draw_go().", call. = FALSE)
  }
  suppressPackageStartupMessages(library(ggplot2))
  ann <- data.frame(clust = scores$cluster)
  rownames(ann) <- names(scores$cluster)
  pl <- env$draw_go(scores$subset, rownames(scores$subset),
                    ann, scores$cluster_colors,
                    n.terms = 10, only.go = FALSE,
                    term.length = 70, cut.pval = 40)
  grDevices::pdf(file, width = 12, height = 14, useDingbats = FALSE)
  tryCatch(print(pl), finally = grDevices::dev.off())
  if (!file.exists(file) || file.info(file)$size < 1000) {
    stop("S13E PDF missing or empty.", call. = FALSE)
  }
  invisible(file)
}

.tf_pathway_write_gz <- function(d, file) {
  con <- gzfile(file, open = "wt")
  tryCatch(utils::write.table(d, con, sep = "\t", row.names = FALSE,
                              quote = FALSE, na = "NA"), finally = close(con))
  x <- utils::read.delim(file, check.names = FALSE)
  if (!identical(dim(x), dim(d)) || !identical(names(x), names(d))) {
    stop("Pathway Source Data dimensions/columns roundtrip failed: ", file,
         call. = FALSE)
  }
  invisible(file)
}

# Export plot inputs and review PDFs from SAVED 25 TF checkpoints only.
# Default all 3, but a preliminary run may request c("7F", "S13F")
# before verifying RNA_analysis/R/draw_go.R for S13E.
tf_pathway_export <- function(
  panels = c("7F", "S13E", "S13F"),
  root = getwd(),
  output_root = file.path(root, "TF_overexpression/publication")
) {
  valid <- c("7F", "S13E", "S13F")
  if (!length(panels) || anyNA(panels) || any(!panels %in% valid) ||
      anyDuplicated(panels)) {
    stop("Choose unique pathway figure panels from 7F, S13E, S13F.", call. = FALSE)
  }
  needed <- c("pheatmap", "RColorBrewer", "ggplot2", "ggrepel", "cowplot")
  needed <- if (identical(panels, "7F")) c("pheatmap", "RColorBrewer") else needed
  ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) {
    stop("Missing plotting packages: ", paste(needed[!ok], collapse = ", "),
         "; no package installation performed.", call. = FALSE)
  }
  root <- normalizePath(root, mustWork = TRUE)
  if ("S13E" %in% panels &&
      !file.exists(file.path(root, "RNA_analysis/R/draw_go.R"))) {
    stop("Missing original S13E draw_go.R helper; inspect historic paths first.",
         call. = FALSE)
  }
  out <- normalizePath(output_root, mustWork = FALSE)
  approved_hist <- normalizePath(file.path(root, "TF_overexpression/figures"),
                                 mustWork = TRUE)
  if (identical(out, approved_hist) ||
      startsWith(out, paste0(approved_hist, "/"))) {
    stop("Refusing to overwrite historical TF figures.", call. = FALSE)
  }
  targets <- list(
    `7F` = c("Fig7F_all_TF_pathway_enrichment.tsv.gz",
             "Fig7F_selected_pathway_heatmap_values.tsv.gz",
             "Fig7F_TF_pathway_heatmap_REVIEW.pdf"),
    `S13E` = c("FigS13E_pathway_category_candidates.tsv.gz",
               "FigS13E_TF_pathway_examples_REVIEW.pdf"),
    `S13F` = c("FigS13F_selected_pathway_associations.tsv.gz",
               "FigS13F_TF_pathway_associations_REVIEW.pdf")
  )
  outputfiles <- unlist(lapply(panels, function(x) {
    c(file.path(out, "source_data", grep("\\.tsv\\.gz$", targets[[x]], value = TRUE)),
      file.path(out, "figures", grep("\\.pdf$", targets[[x]], value = TRUE)))
  }), use.names = FALSE)
  if (any(file.exists(outputfiles))) {
    stop("Existing pathway figure outputs are protected: ",
         paste(outputfiles[file.exists(outputfiles)], collapse = "; "),
         call. = FALSE)
  }
  obj <- .tf_pathway_figures_load(root, file.path(root, "TF_overexpression/publication"))
  scores <- .tf_pathway_figures_scores(obj)
  all <- .tf_pathway_export_all(obj, scores)
  # Build S13F selected TF sets from the saved in-vivo statistics, no refitting.
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".tf_pathway_export_", tmpdir = out)
  dir.create(file.path(stage, "source_data"), recursive = TRUE)
  dir.create(file.path(stage, "figures"), recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  report <- list()
  for (panel in panels) {
    if (panel == "7F") {
      file1 <- file.path(stage, "source_data", targets[[panel]][1])
      file2 <- file.path(stage, "source_data", targets[[panel]][2])
      fig <- file.path(stage, "figures", targets[[panel]][3])
      .tf_pathway_write_gz(all, file1)
      .tf_pathway_write_gz(.tf_pathway_export_7F(all, scores), file2)
      .tf_pathway_figure_7F(scores, fig)
      rows <- nrow(all)
    } else if (panel == "S13E") {
      file1 <- file.path(stage, "source_data", targets[[panel]][1])
      fig <- file.path(stage, "figures", targets[[panel]][2])
      d <- .tf_pathway_export_7F(all, scores)
      d$caption_scope <- "candidate_pathways_by_fig7f_cluster; helper_selects_examples"
      .tf_pathway_write_gz(d, file1)
      .tf_pathway_figure_S13E(scores, root, fig)
      rows <- nrow(d)
    } else if (panel == "S13F") {
      file1 <- file.path(stage, "source_data", targets[[panel]][1])
      fig <- file.path(stage, "figures", targets[[panel]][2])
      d <- .tf_pathway_export_S13F(obj, all, root)
      .tf_pathway_write_gz(d, file1)
      .tf_pathway_figure_S13F(d, fig)
      rows <- nrow(d)
    }
    report[[panel]] <- data.frame(panel = panel, n_rows = rows,
                                  status = "GENERATED; VISUAL_REVIEW_PENDING")
    cat("GENERATED ", panel, ": ", rows,
        " numerical rows; PDF created; review against historical figure pending.\n",
        sep = "")
  }
  # Commit all panels only after every table and PDF rendered successfully.
  staged <- list.files(stage, recursive = TRUE, full.names = TRUE)
  staged <- staged[!file.info(staged)$isdir]
  for (f in staged) {
    relative <- substring(f, nchar(stage) + 2L)
    dest <- file.path(out, relative)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(dest) || !file.rename(f, dest)) {
      stop("Unable to commit fully rendered pathway export: ", dest,
           call. = FALSE)
    }
  }
  invisible(do.call(rbind, report))
}
# ==== End checkpoint-based Fig7F/S13E/S13F export ====

run_tf_publication <- function(panels = c("7C", "7D"),
                               root = getwd(),
                               output_root = file.path(root, "TF_overexpression/publication")) {
  # S12F uses the archived Pagoda2 checkpoint. Keep it separate from the
  # other panels to avoid reloading the large TF cell-assignment checkpoint.
  if (length(panels) == 1L && identical(as.character(panels), "S12F")) {
    return(reproduce_FigS12F(root = root, output_root = output_root))
  }
  if ("S12F" %in% as.character(panels)) {
    stop("Run S12F separately: run_tf_publication(panels = 'S12F').",
         call. = FALSE)
  }

  pf <- tf_preflight(panels = panels, root = root)
  panels <- unique(panels)
  old <- c("7C", "7D", "S12B")
  ck <- if (any(panels %in% old)) .tf_load_checkpoint(pf$checkpoint) else NULL
  inv <- if (any(panels %in% .tf_stage3_panels()))
    .tf_load_stage3(pf$paths, panels) else NULL
  dirs <- c(file.path(output_root, "source_data"),
            file.path(output_root, "figures"))
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  spec <- list(
    `7C` = list(tab = "Fig7C_TF_lineage_bias.tsv.gz",
                pdf = "Fig7C_TF_lineage_bias_REVIEW.pdf",
                data = .tf_source_7C, plot = .tf_figure_7C),
    `7D` = list(tab = "Fig7D_TF_growth_rate.tsv.gz",
                pdf = "Fig7D_TF_growth_rate_REVIEW.pdf",
                data = .tf_source_7D, plot = .tf_figure_7D),
    `S12B` = list(tab = "FigS12B_TF_functional_effects.tsv.gz",
                  pdf = "FigS12B_TF_functional_effects_REVIEW.pdf",
                  data = .tf_source_S12B, plot = .tf_figure_S12B,
                  writer = .tf_write_figure_S12B),
    `7E` = list(tab = "Fig7E_TF_major_UMAP_cells.tsv.gz",
                pdf = "Fig7E_TF_major_UMAP_REVIEW.pdf",
                data = .tf_source_7E, plot = .tf_figure_7E, dims = c(6, 4)),
    `S12E` = list(tab = "FigS12E_TF_clusters_UMAP_cells.tsv.gz",
                  pdf = "FigS12E_TF_clusters_UMAP_REVIEW.pdf",
                  data = .tf_source_S12E, plot = .tf_figure_S12E, dims = c(5.4, 5.4)),
    `S13A` = list(tab = "FigS13A_TF_barcode_cell_counts.tsv.gz",
                  pdf = "FigS13A_TF_barcode_assignment_REVIEW.pdf",
                  data = .tf_source_S13A, plot = .tf_figure_S13A, dims = c(5.5, 5)),
    `S13B` = list(tab = "FigS13B_TF_DNA_RNA_frequencies.tsv.gz",
                  pdf = "FigS13B_TF_DNA_RNA_frequencies_REVIEW.pdf",
                  data = .tf_source_S13B, plot = .tf_figure_S13B, dims = c(5.5, 5)),
    `S13C` = list(tab = "FigS13C_TF_cells_per_perturbation.tsv.gz",
                  pdf = "FigS13C_TF_cells_per_perturbation_REVIEW.pdf",
                  data = .tf_source_S13C, plot = .tf_figure_S13C, dims = c(5.5, 6.0))
  )
  targets <- unlist(lapply(panels, function(panel) {
    c(file.path(dirs[1], spec[[panel]]$tab),
      file.path(dirs[2], spec[[panel]]$pdf))
  }), use.names = FALSE)
  if (any(file.exists(targets))) {
    stop("Refusing to overwrite outputs: ",
         paste(targets[file.exists(targets)], collapse = "; "),
         ". Use a fresh output_root for independent validation.", call. = FALSE)
  }
  result <- lapply(panels, function(panel) {
    s <- spec[[panel]]
    dat <- s$data(if (panel %in% old) ck else inv)
    plot <- s$plot(dat)
    table_path <- file.path(dirs[1], s$tab)
    figure_path <- file.path(dirs[2], s$pdf)
    .tf_write_table(dat, table_path)
    if (is.null(s$writer)) {
      if (is.null(s$dims)) .tf_write_figure(plot, figure_path)
      else .tf_write_figure(plot, figure_path, width = s$dims[1], height = s$dims[2])
    } else {
      s$writer(plot, figure_path)
    }
    if (panel %in% c("S12B", .tf_stage3_panels())) {
      cat("GENERATED, awaiting visual and legend review: Fig ", panel,
          "; ", nrow(dat), " Source Data rows; 1 PDF.\n", sep = "")
    } else {
      cat("PASS: Fig ", panel, "; ", nrow(dat),
          " Source Data rows; 1 PDF.\n", sep = "")
    }
    data.frame(panel = panel, rows = nrow(dat),
               source_data = table_path, figure = figure_path,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}

# === Figure 7H/7I checkpoint-based gene-to-clone associations ===
# Based on the B scatterplots in TF_overexpression/scripts/8_revision.R.
# Loads the saved DESeq2 consensus signatures and archived clone association
# checkpoint. NO DESeq2 refit or other upstream analysis is performed.
# Figure 7H: myeloid bias, consensus res_myeloid, BH padj < 0.20.
# Figure 7I: growth rate, consensus res_growth, nominal pvalue < 0.01.
# Values are GENE-level TF consensus log2FC versus a gene-level clone-
# phenotype association score; one row is NOT one TF perturbation or clone.

.tf_fig7hi_specs <- function() {
  list(
    `7H` = list(result = "res_myeloid", phenotype = "myeloid bias",
                field = "padj", cutoff = 0.20,
                genes = "Fig7H_TF_myeloid_bias_gene_associations.tsv.gz",
                summary = "Fig7H_TF_myeloid_bias_correlation_summary.tsv",
                figure = "Fig7H_TF_myeloid_bias_REVIEW.pdf",
                expected_n = 96L, expected_r = 0.36201896,
                expected_p = 0.0002897131),
    `7I` = list(result = "res_growth", phenotype = "growth rate",
                field = "pvalue", cutoff = 0.01,
                genes = "Fig7I_TF_growth_rate_gene_associations.tsv.gz",
                summary = "Fig7I_TF_growth_rate_correlation_summary.tsv",
                figure = "Fig7I_TF_growth_rate_REVIEW.pdf",
                expected_n = 43L, expected_r = 0.04104002,
                expected_p = 0.7938632353)
  )
}

.tf_fig7hi_load <- function(root) {
  f1 <- file.path(root, "TF_overexpression/data/processed/TF_consensus_signatures")
  f2 <- file.path(root,
    "Association_analysis/data/processed/use/association_summary.RData")
  if (!file.exists(f1) || !file.exists(f2)) {
    stop("Missing Figure 7H/I checkpoint(s): ",
         paste(c(f1, f2)[!file.exists(c(f1, f2))], collapse = "; "),
         call. = FALSE)
  }
  consensus <- readRDS(f1)
  e <- new.env(parent = baseenv())
  load(f2, envir = e)
  if (!is.list(consensus) ||
      !all(c("res_myeloid", "res_growth") %in% names(consensus)) ||
      !exists("test.basic.list", envir = e, inherits = FALSE) ||
      !all(c("myeloid bias", "growth rate") %in% names(e$test.basic.list))) {
    stop("Unexpected Figure 7H/I consensus or association checkpoint structure.",
         call. = FALSE)
  }
  list(consensus = consensus, phenotypes = e$test.basic.list,
       source_consensus = f1, source_association = f2)
}

.tf_fig7hi_build <- function(panel, inputs) {
  s <- .tf_fig7hi_specs()[[panel]]
  de <- as.data.frame(inputs$consensus[[s$result]])
  ph <- as.data.frame(inputs$phenotypes[[s$phenotype]])
  if (!all(c("log2FoldChange", "pvalue", "padj") %in% names(de)) ||
      !"pheno" %in% names(ph) ||
      is.null(rownames(de)) || is.null(rownames(ph)) ||
      anyDuplicated(rownames(de)) || anyDuplicated(rownames(ph))) {
    stop("Missing statistics or nonunique gene identifiers for ", panel,
         call. = FALSE)
  }
  genes <- intersect(rownames(de), rownames(ph))
  if (!length(genes)) stop("No shared genes in ", panel, call. = FALSE)
  selection_p <- as.numeric(de[genes, s$field])
  log2fc <- as.numeric(de[genes, "log2FoldChange"])
  pheno_score <- as.numeric(ph[genes, "pheno"])
  selected <- !is.na(selection_p) & selection_p < s$cutoff
  finite <- is.finite(log2fc) & is.finite(pheno_score)
  plotted <- selected & finite
  if (sum(plotted) < 3L) stop("Too few finite gene pairs for ", panel,
                             call. = FALSE)
  d <- data.frame(
    panel = panel, gene = genes,
    TF_consensus_log2FoldChange = log2fc,
    TF_consensus_pvalue = as.numeric(de[genes, "pvalue"]),
    TF_consensus_padj = as.numeric(de[genes, "padj"]),
    clone_phenotype_GE_shift = pheno_score,
    selection_P_field = s$field, selection_P_cutoff = s$cutoff,
    selected_by_signature_P = selected,
    finite_gene_pair = finite,
    plotted_in_Fig7HI = plotted,
    stringsAsFactors = FALSE)
  for (field in c("pvalue", "padj")) {
    if (field %in% names(ph)) {
      d[[paste0("clone_phenotype_", field)]] <- as.numeric(ph[genes, field])
    }
  }
  fit <- stats::cor.test(log2fc[plotted], pheno_score[plotted],
                         method = "pearson")
  r <- as.numeric(fit$estimate)
  p <- fit$p.value
  tab <- data.frame(
    panel = panel, phenotype = s$phenotype,
    model_comparison = s$result,
    selection_P_field = s$field, selection_P_cutoff = s$cutoff,
    n_shared_genes = length(genes),
    n_significant_signature_genes = sum(selected),
    n_finite_selected_pairs = sum(plotted),
    correlation_method = "Pearson_two_sided",
    pearson_r = r, pvalue = p,
    provenance = "saved_TF_consensus_and_clone_association_checkpoints",
    stringsAsFactors = FALSE)
  # Historical controls were verified in the original IRIS session.
  # Halt on statistical drift rather than publishing a different result.
  if (length(genes) != 11093L || sum(selected) != s$expected_n ||
      sum(plotted) != s$expected_n ||
      abs(r - s$expected_r) > 1e-5 ||
      !isTRUE(all.equal(p, s$expected_p, tolerance = 1e-4))) {
    stop("Figure ", panel, " drifted from validated checkpoint statistics: ",
         "shared=", length(genes), ", selected=", sum(selected),
         ", plotted=", sum(plotted),
         ", Pearson r=", signif(r, 8), ", P=", signif(p, 8),
         ". Do not replace approved Source Data.", call. = FALSE)
  }
  cat(sprintf("Figure %s: %d shared genes; %d selected; Pearson r=%.8f, P=%.8g; historical numerical check PASS.\n",
              panel, length(genes), sum(plotted), r, p))
  list(all_genes = d, summary = tab)
}

.tf_fig7hi_plot <- function(x, panel) {
  d <- x$all_genes
  d <- d[d$plotted_in_Fig7HI, , drop = FALSE]
  s <- .tf_fig7hi_specs()[[panel]]
  stat <- x$summary
  label <- sprintf("r = %.2f\np = %s", stat$pearson_r,
                   format(stat$pvalue, scientific = TRUE, digits = 2))
  ggplot2::ggplot(d, ggplot2::aes(
      x = TF_consensus_log2FoldChange,
      y = clone_phenotype_GE_shift, label = gene)) +
    ggplot2::geom_point(color = "black", size = 1.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.45) +
    ggrepel::geom_text_repel(max.overlaps = 5, box.padding = 0.5,
          point.padding = 0.3, segment.color = "grey50", size = 3,
          seed = 1) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = label,
          hjust = -0.08, vjust = 1.25, size = 4.1, color = "red3") +
    ggplot2::labs(
      x = "In vitro TF perturbation score (log2FC)",
      y = paste0("Clones with ", s$phenotype, ", GE shift")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 12))
}

.tf_fig7hi_write_gz <- function(d, path) {
  con <- gzfile(path, open = "wt")
  tryCatch(utils::write.table(d, con, sep = "\t", row.names = FALSE,
                              quote = FALSE, na = "NA"), finally = close(con))
  read_con <- gzfile(path, open = "rt")
  back <- tryCatch(
    utils::read.delim(read_con, check.names = FALSE,
                      stringsAsFactors = FALSE),
    finally = close(read_con))
  if (!identical(dim(back), dim(d)) ||
      !identical(names(back), names(d)) ||
      !isTRUE(all.equal(d, back, tolerance = 1e-8,
                        check.attributes = FALSE))) {
    stop("Figure 7H/I numerical Source Data round-trip mismatch: ",
         path, call. = FALSE)
  }
}

.tf_fig7hi_write_summary <- function(d, path) {
  utils::write.table(d, path, sep = "\t", row.names = FALSE,
                     quote = FALSE, na = "NA")
  back <- utils::read.delim(path, check.names = FALSE,
                            stringsAsFactors = FALSE)
  if (!isTRUE(all.equal(d, back, tolerance = 1e-8,
                        check.attributes = FALSE))) {
    stop("Figure 7H/I correlation summary round-trip mismatch.", call. = FALSE)
  }
}

# Keeps all existing panel functions in run_publication.R unchanged.
# Produces two gene-level tables, two correlation-summary tables and two PDFs.
tf_association_export <- function(
    panels = c("7H", "7I"), root = getwd(),
    output_root = file.path(root, "TF_overexpression/publication")) {
  if (!length(panels) || anyNA(panels) ||
      !all(panels %in% c("7H", "7I")) || anyDuplicated(panels)) {
    stop("Specify unique panels from '7H' and '7I'.", call. = FALSE)
  }
  pkgs <- c("ggplot2", "ggrepel")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing package(s): ", paste(missing, collapse = ", "),
         "; no installation attempted.", call. = FALSE)
  }
  root <- normalizePath(root, mustWork = TRUE)
  if (!file.exists(file.path(root, "TF_overexpression/publication/run_publication.R"))) {
    stop("Run from CH_publication_release; TF publication runner absent.",
         call. = FALSE)
  }
  output_root <- normalizePath(output_root, mustWork = FALSE)
  hist <- normalizePath(file.path(root, "TF_overexpression/figures"),
                        mustWork = TRUE)
  if (identical(output_root, hist) ||
      startsWith(output_root, paste0(hist, "/"))) {
    stop("Refusing to write into original historical figures.", call. = FALSE)
  }
  spec <- .tf_fig7hi_specs()
  targets <- unlist(lapply(panels, function(panel) {
    ss <- spec[[panel]]
    c(file.path(output_root, "source_data", ss$genes),
      file.path(output_root, "source_data", ss$summary),
      file.path(output_root, "figures", ss$figure))
  }), use.names = FALSE)
  if (any(file.exists(targets))) {
    stop("Output already exists; refusing overwrite: ",
         paste(targets[file.exists(targets)], collapse = "; "),
         call. = FALSE)
  }
  inputs <- .tf_fig7hi_load(root)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  result <- vector("list", length(panels))
  names(result) <- panels
  for (panel in panels) {
    ss <- spec[[panel]]
    stage <- tempfile(paste0(".stage_Fig", panel, "_"), tmpdir = output_root)
    dir.create(stage)
    committed <- character()
    tryCatch({
      z <- .tf_fig7hi_build(panel, inputs)
      st_genes <- file.path(stage, ss$genes)
      st_summary <- file.path(stage, ss$summary)
      st_figure <- file.path(stage, ss$figure)
      .tf_fig7hi_write_gz(z$all_genes, st_genes)
      .tf_fig7hi_write_summary(z$summary, st_summary)
      grDevices::pdf(st_figure, width = 6, height = 5.4,
                     useDingbats = FALSE)
      tryCatch(print(.tf_fig7hi_plot(z, panel)),
               finally = grDevices::dev.off())
      if (!is.finite(file.info(st_figure)$size) ||
          file.info(st_figure)$size < 1000) {
        stop("Figure ", panel, " PDF empty or missing.", call. = FALSE)
      }
      dst <- c(file.path(output_root, "source_data", ss$genes),
               file.path(output_root, "source_data", ss$summary),
               file.path(output_root, "figures", ss$figure))
      dir.create(dirname(dst[1L]), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(dst[3L]), recursive = TRUE, showWarnings = FALSE)
      src <- c(st_genes, st_summary, st_figure)
      for (i in seq_along(src)) {
        if (file.exists(dst[i]) || !file.rename(src[i], dst[i])) {
          stop("Cannot commit output: ", dst[i], call. = FALSE)
        }
        committed <- c(committed, dst[i])
      }
      result[[panel]] <- data.frame(
        panel = panel, shared_genes = nrow(z$all_genes),
        selected_genes = sum(z$all_genes$plotted_in_Fig7HI),
        pearson_r = z$summary$pearson_r, pvalue = z$summary$pvalue,
        source_data = dst[1L], correlation_summary = dst[2L],
        figure = dst[3L], status = "NUMERICAL_PASS; VISUAL_REVIEW_PENDING",
        stringsAsFactors = FALSE)
      cat("GENERATED Fig", panel, ": ", nrow(z$all_genes),
          " gene-level rows, 1 correlation summary, 1 PDF.\n", sep = "")
    }, error = function(err) {
      if (length(committed)) unlink(committed)
      stop(err)
    }, finally = unlink(stage, recursive = TRUE))
  }
  invisible(do.call(rbind, result))
}

# === Supplementary Figure 13G: TF response vs HSC module correlations ===
# Reproduces TF_overexpression/scripts/6_tf_invitro_de.R, section
# "compare TF-induced changes with clone signatures".
# No DESeq2 or fgsea refitting. The 38-column saved logFC.tf/pv.tf matrices
# and the seven columns of gsea.clone$gscore.clone are read as-is.
# Original: compos.tf = sign(logFC.tf) * -log(pv.tf)
# Selected genes: intersect row identifiers, then any TF nominal p < 1e-3.
# Correlation: stats::cor(compos.tf[selected,], gscore[selected,])
# (default Pearson, NA policy 'everything'), 38 TFs x 7 modules.

.tf_s13g_load <- function(root) {
  a <- file.path(root, 'TF_overexpression/data/processed/TF_invitro_DE.RData')
  b <- file.path(root, 'Association_analysis/data/processed/use/association_summary.RData')
  if (!file.exists(a) || !file.exists(b)) {
    stop('S13G missing saved input: ', paste(c(a,b)[!file.exists(c(a,b))], collapse='; '), call.=FALSE)
  }
  e1 <- new.env(parent=baseenv()); e2 <- new.env(parent=baseenv())
  load(a, envir=e1); load(b, envir=e2)
  if (!all(c('logFC.tf','pv.tf') %in% ls(e1, all.names=TRUE)) ||
      !exists('gsea.clone', envir=e2, inherits=FALSE) ||
      !'gscore.clone' %in% names(e2$gsea.clone)) {
    stop('S13G missing expected saved DE or clone-score objects.', call.=FALSE)
  }
  lfc <- as.matrix(e1$logFC.tf)
  pv <- as.matrix(e1$pv.tf)
  scores <- as.matrix(e2$gsea.clone$gscore.clone)
  if (!identical(dim(lfc), dim(pv)) ||
      !identical(rownames(lfc), rownames(pv)) ||
      !identical(colnames(lfc), colnames(pv)) ||
      !is.numeric(lfc) || !is.numeric(pv) || !is.numeric(scores) ||
      anyDuplicated(rownames(lfc)) || anyDuplicated(rownames(scores)) ||
      anyDuplicated(colnames(lfc)) || anyDuplicated(colnames(scores))) {
    stop('S13G dimension, name or numeric-type inconsistency.', call.=FALSE)
  }
  if (any(pv < 0 | pv > 1, na.rm=TRUE) || any(pv == 0, na.rm=TRUE)) {
    stop('S13G missing/zero/invalid P values require a historical semantics audit.', call.=FALSE)
  }
  list(lfc=lfc, pv=pv, scores=scores, input_de=a, input_modules=b)
}

.tf_s13g_calculate <- function(x) {
  common <- intersect(rownames(x$lfc), rownames(x$scores))
  p <- x$pv[common, , drop=FALSE]
  selected <- rowSums(is.finite(p) & p < 1e-3) > 0
  genes <- common[selected]
  if (!length(genes)) stop('S13G has no selected genes.', call.=FALSE)
  lfc <- x$lfc[genes, , drop=FALSE]
  p <- x$pv[genes, , drop=FALSE]
  module <- x$scores[genes, , drop=FALSE]
  # Preserve historical rule without fallback imputation or modified ordering.
  ranking <- sign(lfc) * -log(p)
  corr <- stats::cor(ranking, module, method='pearson', use='everything')
  if (!identical(dim(corr), c(ncol(x$lfc), ncol(x$scores))) ||
      any(!is.finite(corr))) {
    stop('S13G produced undefined correlation(s). Check missing gene scores.', call.=FALSE)
  }
  n <- length(genes)
  selected_tab <- data.frame(
    gene=genes, n_TFs_with_nominal_P_lt_0.001=rowSums(is.finite(p) & p < 1e-3),
    minimum_TF_nominal_P=apply(p, 1, min, na.rm=TRUE),
    selection_rule='any_TF_nominal_P_lt_0.001', stringsAsFactors=FALSE)
  tf_tab <- data.frame(
    gene=rep(genes, ncol(lfc)), TF=rep(colnames(lfc), each=n),
    TF_log2FoldChange=as.numeric(lfc), TF_nominal_pvalue=as.numeric(p),
    TF_signed_neg_ln_P=as.numeric(ranking), stringsAsFactors=FALSE)
  module_tab <- data.frame(
    gene=rep(genes, ncol(module)),
    HSC_module=rep(colnames(module), each=n),
    module_gene_score=as.numeric(module), stringsAsFactors=FALSE)
  corr_tab <- data.frame(
    TF=rep(rownames(corr), times=ncol(corr)),
    HSC_module=rep(colnames(corr), each=nrow(corr)),
    pearson_r=as.numeric(corr),
    n_selected_genes=n,
    correlation_statistic='Pearson_cor_signed_negative_ln_TF_P_and_HSC_module_gene_score',
    stringsAsFactors=FALSE)
  if (n != 921L || nrow(corr_tab) != 266L ||
      ncol(x$lfc) != 38L || ncol(x$scores) != 7L) {
    stop('S13G drifted from IRIS preflight counts: genes=',n,
         ', TFs=', ncol(x$lfc), ', modules=', ncol(x$scores), call.=FALSE)
  }
  cat(sprintf('S13G: %d common genes; %d selected; %d TFs x %d modules = %d correlations.\n',
              length(common), n, nrow(corr), ncol(corr), length(corr)))
  list(selected=selected_tab, tf=tf_tab, module=module_tab,
       correlations=corr_tab, matrix=corr)
}

.tf_s13g_paths <- function(output_root) {
  list(
    selected=file.path(output_root, 'source_data/FigS13G_selected_genes.tsv.gz'),
    tf=file.path(output_root, 'source_data/FigS13G_TF_gene_scores.tsv.gz'),
    module=file.path(output_root, 'source_data/FigS13G_HSC_module_gene_scores.tsv.gz'),
    correlations=file.path(output_root, 'source_data/FigS13G_TF_module_correlations.tsv.gz'),
    figure=file.path(output_root, 'figures/FigS13G_TF_module_correlations_REVIEW.pdf')
  )
}

.tf_s13g_write_table <- function(d, file) {
  z <- gzfile(file, open='wt')
  tryCatch(utils::write.table(d, z, sep='\t', row.names=FALSE,
                              quote=FALSE, na='NA'), finally=close(z))
  z <- gzfile(file, open='rt')
  back <- tryCatch(utils::read.delim(z, check.names=FALSE,
                         stringsAsFactors=FALSE), finally=close(z))
  if (!identical(dim(d), dim(back)) || !identical(names(d), names(back)) ||
      !isTRUE(all.equal(d, back, tolerance=1e-8,
                        check.attributes=FALSE))) {
    stop('S13G source-data round-trip mismatch for ', file, call.=FALSE)
  }
}

.tf_s13g_pdf <- function(corr, file) {
  if (!requireNamespace('pheatmap', quietly=TRUE)) {
    stop('S13G requires existing package pheatmap.', call.=FALSE)
  }
  limit <- max(abs(corr))
  if (!is.finite(limit) || limit == 0) stop('S13G correlation scale invalid.', call.=FALSE)
  p <- pheatmap::pheatmap(corr, angle_col=45,
                         main='TF perturbation response and clone phenotypes',
                         fontsize=10,
                         breaks=seq(-limit,limit,length.out=100),
                         silent=TRUE)
  grDevices::pdf(file, width=5, height=8, useDingbats=FALSE)
  tryCatch(grid::grid.draw(p$gtable), finally=grDevices::dev.off())
  if (!file.exists(file) || file.info(file)$size < 1000) {
    stop('S13G PDF absent or empty: ', file, call.=FALSE)
  }
}

# One command exports FOUR fully reconstructive numerical tables + heatmap PDF.
# validate_replay=TRUE runs a second clean checkpoint export into an isolated
# subdirectory, then checks all four tables. It does not refit models.
tf_s13g_export <- function(root=getwd(),
    output_root=file.path(root,'TF_overexpression/publication'),
    validate_replay=TRUE) {
  root <- normalizePath(root, mustWork=TRUE)
  output_root <- normalizePath(output_root, mustWork=FALSE)
  if (!file.exists(file.path(root,'TF_overexpression/publication/run_publication.R'))) {
    stop('Start in CH_publication_release.', call.=FALSE)
  }
  historical <- normalizePath(file.path(root,'TF_overexpression/figures'),
                              mustWork=TRUE)
  if (identical(output_root,historical) ||
      startsWith(output_root,paste0(historical,'/'))) {
    stop('Refusing to write into historical figures.', call.=FALSE)
  }
  paths <- .tf_s13g_paths(output_root)
  if (any(file.exists(unlist(paths)))) {
    stop('Existing S13G output(s) preserved. Choose a fresh output_root.',
         call.=FALSE)
  }
  # Precheck replay destination before beginning the approved main export.
  replay_root <- file.path(root,
        'TF_overexpression/publication/validation/S13G_clean')
  if (validate_replay && dir.exists(replay_root)) {
    stop('Existing S13G clean replay directory; review it before retrying.',call.=FALSE)
  }
  pkgs <- c('pheatmap')
  missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop('Missing S13G package(s): ',paste(missing,collapse=', '))
  x <- .tf_s13g_load(root)
  data <- .tf_s13g_calculate(x)
  rm(x); invisible(gc())
  dir.create(file.path(output_root,'source_data'),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(output_root,'figures'),recursive=TRUE,showWarnings=FALSE)
  for (field in c('selected','tf','module','correlations')) {
    .tf_s13g_write_table(data[[field]], paths[[field]])
  }
  .tf_s13g_pdf(data$matrix, paths$figure)
  cat('S13G Source Data: ',nrow(data$selected),' selected genes, ',
      nrow(data$tf),' TF-gene rows, ',nrow(data$module),' module-gene rows, ',
      nrow(data$correlations),' TF-module Pearson correlations.\n',sep='')
  cat('S13G PDF generated; historical visual agreement: REVIEW.\n')
  if (validate_replay) {
    # Recalculate from saved checkpoints in an isolated root; no stale in-memory
    # objects and no copied numerical Source Data.
    replay <- tf_s13g_export(root=root, output_root=replay_root,
                             validate_replay=FALSE)
    report <- do.call(rbind, lapply(c('selected','tf','module','correlations'),
        function(field) {
      f <- function(path) {
        z <- gzfile(path,open='rt')
        tryCatch(utils::read.delim(z,check.names=FALSE,
                   stringsAsFactors=FALSE),finally=close(z))
      }
      a <- f(paths[[field]]); b <- f(replay$paths[[field]])
      ok <- identical(dim(a),dim(b)) && identical(names(a),names(b)) &&
        isTRUE(all.equal(a,b,tolerance=1e-8,check.attributes=FALSE))
      data.frame(table=basename(paths[[field]]),rows=nrow(a),
                 status=if(ok) 'PASS' else 'REVIEW',stringsAsFactors=FALSE)
    }))
    report_file <- file.path(root,
      'TF_overexpression/publication/validation/S13G_clean_comparison.tsv')
    utils::write.table(report,report_file,sep='\t',quote=FALSE,row.names=FALSE)
    print(report,row.names=FALSE)
    if (!all(report$status=='PASS')) stop('S13G clean numerical replay mismatch.')
    cat('PASS: S13G 4/4 numerical checkpoint comparisons.\n')
  }
  invisible(list(paths=paths,statistics=data$correlations,
                 selected_genes=nrow(data$selected),
                 status='NUMERICAL_PASS; HISTORICAL_VISUAL_REVIEW_PENDING'))
}


# === Supplementary Figure 13D: DE volcanoes from saved matrices ===
# Historical source: TF_overexpression/scripts/6_tf_invitro_de.R,
# section `fdr.tf <- apply(pv.tf, 2, ...)` through
# `figures/invitro/TF_response_scatterplots.pdf`.
# 25 TFs with at least 4 saved assigned cells; no DESeq2 refitting.
# x = shrunken log2FC, y = nominal P (log10), red if BH FDR < 0.30,
# label genes if BH FDR <= 0.20, max.overlaps = 15.

.tf_s13d_inputs <- function(root) {
  de_file <- file.path(root, 'TF_overexpression/data/processed/TF_invitro_DE.RData')
  count_file <- file.path(root, 'TF_overexpression/publication/source_data/FigS13C_TF_cells_per_perturbation.tsv.gz')
  for (p in c(de_file, count_file)) {
    if (!file.exists(p)) stop('S13D missing prerequisite: ', p, call.=FALSE)
  }
  e <- new.env(parent=baseenv())
  load(de_file, envir=e)
  if (!all(c('logFC.tf','pv.tf') %in% ls(e, all.names=TRUE))) {
    stop('S13D: missing saved logFC.tf or pv.tf.', call.=FALSE)
  }
  lfc <- as.matrix(e$logFC.tf); pv <- as.matrix(e$pv.tf)
  if (!identical(dim(lfc), dim(pv)) ||
      !identical(rownames(lfc), rownames(pv)) ||
      !identical(colnames(lfc), colnames(pv)) ||
      is.null(rownames(lfc)) || is.null(colnames(lfc)) ||
      anyDuplicated(rownames(lfc)) || anyDuplicated(colnames(lfc))) {
    stop('S13D: incompatible differential-expression matrices.', call.=FALSE)
  }
  if (any(!is.finite(pv)) || any(pv <= 0 | pv > 1)) {
    stop('S13D: zero, missing or invalid nominal P requires separate audit.', call.=FALSE)
  }
  tfc <- utils::read.delim(count_file, check.names=FALSE, stringsAsFactors=FALSE)
  if (!all(c('TF','assigned_cells_in_p2_saved_TFn') %in% names(tfc)) ||
      anyDuplicated(tfc$TF)) stop('S13D: invalid S13C TF count table.', call.=FALSE)
  tf_n <- tfc$assigned_cells_in_p2_saved_TFn[match(colnames(lfc), tfc$TF)]
  if (anyNA(tf_n)) stop('S13D: missing DE TF(s) in S13C table.', call.=FALSE)
  keep <- which(tf_n >= 4)
  expected <- c('GFP','ID3','ID4','ID1','DBP','TEF','HLF','CNBP','ATF2',
                'HOXB1','HMX1','TAL1','MEIS3','ZSCAN1','LDB2','CDK1',
                'HMBOX1','JUN','PBX3','SIX3','TCF7','BARHL2','TBX15','IRX3','NFYC')
  if (!identical(colnames(lfc)[keep], expected) || nrow(lfc) != 11176L) {
    stop('S13D eligible TFs or DE feature universe differs from historical preflight.',
         call.=FALSE)
  }
  list(lfc=lfc[,keep,drop=FALSE], pv=pv[,keep,drop=FALSE],
       counts=as.integer(tf_n[keep]), source_de=de_file, source_counts=count_file)
}

.tf_s13d_tables <- function(x) {
  ng <- nrow(x$lfc); nt <- ncol(x$lfc)
  padj <- apply(x$pv, 2, stats::p.adjust, method='BH')
  stopifnot(identical(dim(padj), dim(x$lfc)))
  df <- data.frame(
    gene=rep(rownames(x$lfc), times=nt),
    TF=rep(colnames(x$lfc), each=ng),
    log2FoldChange=as.numeric(x$lfc),
    pvalue=as.numeric(x$pv),
    BH_FDR=as.numeric(padj),
    assigned_cells=rep(x$counts, each=ng),
    stringsAsFactors=FALSE)
  df$highlight_FDR_lt_0.30 <- is.finite(df$BH_FDR) & df$BH_FDR < 0.30
  df$label_FDR_le_0.20 <- is.finite(df$BH_FDR) & df$BH_FDR <= 0.20
  df$plot_include <- is.finite(df$log2FoldChange) &
    is.finite(df$pvalue) & df$pvalue > 0
  tf_order <- colnames(x$lfc)
  summary <- do.call(rbind, lapply(tf_order, function(tf) {
    z <- df[df$TF == tf, , drop=FALSE]
    data.frame(TF=tf, n_assigned_cells=z$assigned_cells[1],
       n_tested_genes=nrow(z),
       n_FDR_lt_0.30=sum(z$highlight_FDR_lt_0.30),
       n_FDR_le_0.20=sum(z$label_FDR_le_0.20),
       n_plottable=sum(z$plot_include), stringsAsFactors=FALSE)
  }))
  rownames(summary) <- NULL
  stopifnot(nrow(df)==279400L, nrow(summary)==25L)
  list(genes=df, summary=summary, TF_order=tf_order)
}

.tf_s13d_paths <- function(out) {
  list(genes=file.path(out,'source_data/FigS13D_TF_DE_volcano_genes.tsv.gz'),
       summary=file.path(out,'source_data/FigS13D_TF_DE_volcano_summary.tsv'),
       pdf=file.path(out,'figures/FigS13D_TF_DE_volcano_REVIEW.pdf'))
}

.tf_s13d_write <- function(dat, path) {
  if (grepl('\\.gz$', path)) {
    con <- gzfile(path, 'wt')
    tryCatch(utils::write.table(dat, con, sep='\t', quote=FALSE,
                                row.names=FALSE, na='NA'), finally=close(con))
  } else {
    utils::write.table(dat, path, sep='\t', quote=FALSE,
                       row.names=FALSE, na='NA')
  }
  back <- utils::read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
  if (!identical(dim(dat), dim(back)) || !identical(names(dat),names(back)) ||
      !isTRUE(all.equal(dat, back, tolerance=1e-8,
                        check.attributes=FALSE))) {
    stop('S13D Source Data round-trip mismatch: ', path, call.=FALSE)
  }
}

.tf_s13d_plot <- function(tables, path) {
  if (!requireNamespace('ggplot2', quietly=TRUE) ||
      !requireNamespace('ggrepel', quietly=TRUE)) {
    stop('S13D needs existing ggplot2 and ggrepel.', call.=FALSE)
  }
  pdat <- tables$genes[tables$genes$plot_include, , drop=FALSE]
  pdat$TF <- factor(pdat$TF, levels=tables$TF_order)
  labels <- pdat[pdat$label_FDR_le_0.20, , drop=FALSE]
  g <- ggplot2::ggplot(pdat,
       ggplot2::aes(x=log2FoldChange,y=pvalue,
                    color=highlight_FDR_lt_0.30)) +
    ggplot2::geom_point(size=0.2, alpha=1) +
    ggrepel::geom_text_repel(data=labels,
         ggplot2::aes(label=gene), max.overlaps=15, seed=1,
         size=2.3, show.legend=FALSE) +
    ggplot2::theme_bw() + ggplot2::scale_y_log10() +
    ggplot2::labs(x='logFC',y='p-value, log') +
    ggplot2::scale_color_manual(
      values=c('TRUE'='darkred','FALSE'='grey'),drop=FALSE) +
    ggplot2::theme(text=ggplot2::element_text(size=12),
                   legend.position='none') +
    ggplot2::facet_wrap(~TF,scales='free',ncol=5)
  grDevices::pdf(path,width=14,height=14,useDingbats=FALSE)
  tryCatch(print(g),finally=grDevices::dev.off())
  if (!file.exists(path) || file.info(path)$size < 1000) {
    stop('S13D plot missing/empty: ',path,call.=FALSE)
  }
}

# Checkpoint-level numeric reproduction. Safe default: rejects pre-existing targets.
# Source Data includes ALL 11,176 genes for each of 25 TFs, not just labels.
tf_s13d_export <- function(root=getwd(),
    output_root=file.path(root,'TF_overexpression/publication'),
    validate_replay=TRUE) {
  root <- normalizePath(root,mustWork=TRUE)
  output_root <- normalizePath(output_root,mustWork=FALSE)
  if (!file.exists(file.path(root,'TF_overexpression/publication/run_publication.R'))) {
    stop('S13D: run from CH_publication_release.',call.=FALSE)
  }
  paths <- .tf_s13d_paths(output_root)
  if (any(file.exists(unlist(paths)))) {
    stop('S13D outputs already exist; preserve them or choose a fresh output_root.',
         call.=FALSE)
  }
  replay_root <- file.path(root,'TF_overexpression/publication/validation/S13D_clean')
  if (validate_replay && dir.exists(replay_root)) {
    stop('S13D clean directory already exists; choose a clean run before retry.',call.=FALSE)
  }
  source <- .tf_s13d_inputs(root)
  data <- .tf_s13d_tables(source)
  rm(source); invisible(gc())
  dir.create(file.path(output_root,'source_data'),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(output_root,'figures'),recursive=TRUE,showWarnings=FALSE)
  .tf_s13d_write(data$genes, paths$genes)
  .tf_s13d_write(data$summary, paths$summary)
  .tf_s13d_plot(data, paths$pdf)
  cat('GENERATED S13D: ',nrow(data$genes),' TF-gene Source Data rows, ',
      nrow(data$summary),' TF summaries, one 25-panel PDF.\n',sep='')
  if (validate_replay) {
    # Independent checkpoint reload; the replay is not a DESeq2 refit.
    secondary <- .tf_s13d_tables(.tf_s13d_inputs(root))
    comparison <- do.call(rbind, lapply(c('genes','summary'),function(key) {
      z <- utils::read.delim(paths[[key]],check.names=FALSE,
                             stringsAsFactors=FALSE)
      check <- secondary[[key]]
      ok <- identical(dim(z),dim(check)) && identical(names(z),names(check)) &&
        isTRUE(all.equal(z,check,tolerance=1e-8,check.attributes=FALSE))
      data.frame(table=basename(paths[[key]]),rows=nrow(z),
                 status=if (ok) 'PASS' else 'REVIEW',stringsAsFactors=FALSE)
    }))
    report <- file.path(root,'TF_overexpression/publication/validation/S13D_clean_comparison.tsv')
    dir.create(dirname(report),recursive=TRUE,showWarnings=FALSE)
    if (file.exists(report)) stop('Existing S13D validation report preserved.',call.=FALSE)
    utils::write.table(comparison,report,sep='\t',quote=FALSE,row.names=FALSE)
    print(comparison,row.names=FALSE)
    if (!all(comparison$status=='PASS')) stop('S13D checkpoint replay mismatch.',call.=FALSE)
    cat('PASS: S13D two numerical tables independently re-exported from checkpoint.\n')
  }
  invisible(list(status='NUMERICAL_PASS; VISUAL_REVIEW_PENDING',
                 source_data=c(paths$genes,paths$summary),figure=paths$pdf,
                 TFs=data$TF_order,n_genes=nrow(data$genes)))
}

