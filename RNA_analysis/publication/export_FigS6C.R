#!/usr/bin/env Rscript
# CH manuscript: Supplementary Figure 6C Source Data.
# Inputs: archived DESeq2 results in mutant_de.RData only.
# Does not rerun DESeq2, fgsea, Pagoda2, or load scrna.RData.
# Run from the CH_publication_release repository root.

root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "RNA_analysis", "scripts", "5_mutant_DE.R"))) {
  stop("Please run from the CH_publication_release repository root.")
}

data_root <- Sys.getenv("CH_DATA_ROOT", unset = root)
input_file <- file.path(data_root, "RNA_analysis", "data", "processed", "mutant_de.RData")
if (!file.exists(input_file)) {
  stop("Input missing: ", input_file,
       "\nSet CH_DATA_ROOT to the directory containing RNA_analysis if needed.")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("Install/load DESeq2 to read the archived DESeqResults objects.")
}

outdir <- file.path(root, "RNA_analysis", "publication", "source_data")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ref <- new.env(parent = globalenv())
loaded <- load(input_file, envir = ref)
if (!("de.mut2wt.shrink" %in% loaded)) stop("Archived file lacks de.mut2wt.shrink")
n_states <- length(ref$de.mut2wt.shrink)
if (n_states != 10L) stop("Expected 10 states; found ", n_states)

states <- vector("list", n_states)
for (cl in seq_len(n_states)) {
  contrast <- ref$de.mut2wt.shrink[[cl]]
  if (!all(c("DNMT3A", "TET2") %in% names(contrast))) {
    stop("Missing genotype contrast in state ", cl)
  }
  dn <- as.data.frame(contrast$DNMT3A)
  te <- as.data.frame(contrast$TET2)
  cols <- c("log2FoldChange", "pvalue", "padj")
  if (!all(cols %in% names(dn)) || !all(cols %in% names(te))) {
    stop("Missing DESeq2 columns in state ", cl)
  }
  gene_ids <- rownames(dn)
  if (is.null(gene_ids) || anyDuplicated(gene_ids) ||
      !identical(gene_ids, rownames(te))) {
    stop("DESeq2 gene identities/order do not match in state ", cl)
  }

  # Assign gene row names, as in the original code. When ten data.frames are
  # rbind-ed, base R disambiguates repeated row names. This preserves the
  # historical ggplot text-label candidates without changing the gene_id.
  states[[cl]] <- data.frame(
    cluster = rep.int(cl, nrow(dn)),
    gene_id = gene_ids,
    log2FoldChange_DNMT3A = dn$log2FoldChange,
    log2FoldChange_TET2 = te$log2FoldChange,
    pvalue_DNMT3A = dn$pvalue,
    pvalue_TET2 = te$pvalue,
    padj_DNMT3A = dn$padj,
    padj_TET2 = te$padj,
    row.names = gene_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  message("Prepared cluster ", cl, ": ", nrow(states[[cl]]), " genes")
}

full <- do.call(rbind, states)
full$historical_plot_rowname <- rownames(full)

# Replicate 5_mutant_DE.R lines ~287-294, with explicit NA handling.
# A call is assigned a color ONLY if both relevant FDRs are finite.
p_dn <- full$padj_DNMT3A
p_te <- full$padj_TET2
valid <- is.finite(p_dn) & is.finite(p_te)
full$plot_group <- NA_character_
full$plot_group[valid & p_dn < 0.1 & p_te < 0.1] <- "both"
full$plot_group[valid & p_dn < 0.1 & p_te > 0.1] <- "DNMT3A"
full$plot_group[valid & p_dn > 0.1 & p_te < 0.1] <- "TET2"

# IMPORTANT HISTORICAL DETAIL: in 5_mutant_DE.R lines ~296-305, BOTH
# label-selection clauses use DNMT3A padj, including the clause testing
# abs(TET2 log2FoldChange) > 0.2. We preserve that rule here intentionally.
# It should NOT be silently changed to TET2 padj for publication reproduction.
label_rule <- (
  is.finite(p_dn) & p_dn < 0.01 &
    is.finite(full$log2FoldChange_DNMT3A) &
    abs(full$log2FoldChange_DNMT3A) > 0.4
) | (
  is.finite(p_dn) & p_dn < 0.01 &
    is.finite(full$log2FoldChange_TET2) &
    abs(full$log2FoldChange_TET2) > 0.2
)
full$plot_label_candidate <- ifelse(
  label_rule & !is.na(full$plot_group),
  full$historical_plot_rowname,
  NA_character_
)

# Original code used Pearson correlation (default cor()) for each cluster,
# with use="na.or.complete", then rounded the displayed text to two digits.
correlations <- do.call(rbind, lapply(seq_len(n_states), function(cl) {
  x <- full[full$cluster == cl, , drop = FALSE]
  good <- is.finite(x$log2FoldChange_DNMT3A) &
    is.finite(x$log2FoldChange_TET2)
  if (sum(good) < 3L) stop("Too few paired finite logFC values in cluster ", cl)
  r <- stats::cor(x$log2FoldChange_DNMT3A,
                  x$log2FoldChange_TET2,
                  method = "pearson", use = "na.or.complete")
  data.frame(
    cluster = cl,
    n_genes_total = nrow(x),
    n_genes_paired_finite = sum(good),
    pearson_r = as.numeric(r),
    displayed_r = round(as.numeric(r), 2),
    label_x = min(x$log2FoldChange_DNMT3A, na.rm = TRUE) + 0.2,
    label_y = min(x$log2FoldChange_TET2, na.rm = TRUE) + 0.05,
    n_colored_points = sum(!is.na(x$plot_group)),
    n_label_candidates = sum(!is.na(x$plot_label_candidate)),
    stringsAsFactors = FALSE
  )
}))

write_tsv <- function(x, dest) {
  utils::write.table(x, file = dest, sep = "\t", row.names = FALSE,
                     col.names = TRUE, quote = FALSE, na = "NA")
}

outfile <- file.path(outdir, "FigS6C_DE_concordance.tsv.gz")
con <- gzfile(outfile, open = "wt", compression = 6)
tryCatch(
  utils::write.table(full, file = con, sep = "\t", row.names = FALSE,
                     col.names = TRUE, quote = FALSE, na = "NA"),
  finally = close(con)
)
write_tsv(correlations, file.path(outdir, "FigS6C_correlations.tsv"))
write_tsv(data.frame(
  metric = c("cell_states", "gene_state_rows", "plot_points_with_paired_finite_logFC",
             "colored_points", "gene_label_candidates"),
  value = c(n_states, nrow(full),
            sum(is.finite(full$log2FoldChange_DNMT3A) &
                is.finite(full$log2FoldChange_TET2)),
            sum(!is.na(full$plot_group)),
            sum(!is.na(full$plot_label_candidate)))
), file.path(outdir, "FigS6C_export_summary.tsv"))

writeLines(c(
  "Supplementary Figure 6C: publication reproduction provenance",
  "Input: RNA_analysis/data/processed/mutant_de.RData",
  "Script basis: RNA_analysis/scripts/5_mutant_DE.R, DE-scatterplot section (~lines 285-329)",
  "Each row: one gene in one cell state, compared across DNMT3A and TET2 contrasts.",
  "X: DESeq2 normal-shrunken log2FoldChange, DNMT3A vs WT.",
  "Y: DESeq2 normal-shrunken log2FoldChange, TET2 vs WT.",
  "Color groups: both FDR<0.1; DNMT3A-only if DNMT3A FDR<0.1 and TET2 FDR>0.1; vice versa for TET2-only.",
  "Unclassified points: both not significant or either adjusted p-value missing or on cutoff.",
  "Historical label rule: DNMT3A padj<0.01 AND (abs(DNMT3A log2FC)>0.4 OR abs(TET2 log2FC)>0.2); colored points only.",
  "IMPORTANT: historical code uses DNMT3A FDR for BOTH label-selection clauses; rule preserved deliberately.",
  "plot_label_candidate indicates eligibility, NOT necessarily text displayed after ggrepel collision filtering.",
  "historical_plot_rowname retains original rbind-disambiguated row names; gene_id is the unmodified gene symbol/identifier.",
  "Pearson correlations use cor(method='pearson',use='na.or.complete'); displayed values are rounded to 2 decimals.",
  "Tables provide raw and adjusted DE p-values in addition to the plotted log2FC values.",
  "The export loads archived results; no DESeq2, fgsea, or Pagoda2 calculations are rerun.",
  "A separate plot script reads ONLY these TSV outputs; compare its PDF with the accepted Fig. S6C artwork."
), file.path(outdir, "FigS6C_provenance.txt"))

print(correlations[, c("cluster", "n_genes_paired_finite", "displayed_r",
                       "n_colored_points", "n_label_candidates")], row.names = FALSE)
message("Exported Source Data to ", outdir)
