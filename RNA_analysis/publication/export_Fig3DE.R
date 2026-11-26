#!/usr/bin/env Rscript
# CH publication: canonical Figure 3D/E numerical Source Data.
# Uses accepted-figure pathway selection, archived DESeq2 and fgsea results.
# This script does NOT refit DESeq2 or recompute fgsea.
# Run from CH_publication_release root.
root <- normalizePath(getwd(), mustWork = TRUE)
pub <- file.path(root, "RNA_analysis", "publication")
stopifnot(file.exists(file.path(root, "RNA_analysis", "scripts", "6_mutant_analysis.R")))
manifest_file <- file.path(pub, "resources", "Fig3DE_published_pathway_ids.tsv")
if (!file.exists(manifest_file)) stop("Missing accepted pathway manifest: ", manifest_file)
selection <- utils::read.delim(manifest_file, stringsAsFactors = FALSE)
stopifnot(all(c("panel", "genotype", "published_order_top_to_bottom",
                "pathway_id", "published_label", "published_direction") %in% names(selection)))
input_root <- Sys.getenv("CH_DATA_ROOT", unset = root)
input_path <- file.path(input_root, "RNA_analysis", "data", "processed", "mutant_de.RData")
if (!file.exists(input_path)) stop("Archived data missing: ", input_path)
if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 is required to read saved results")
outdir <- file.path(pub, "source_data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ref <- new.env(parent = globalenv())
loaded <- load(input_path, envir = ref)
needed <- c("de.mut2wt.shrink", "fgsea.dnmt3a.clust", "fgsea.tet2.clust")
if (!all(needed %in% loaded)) stop("Missing archived objects: ", paste(setdiff(needed,loaded),collapse=", "))
write_tsv <- function(x, filename) {
  utils::write.table(x, filename, sep = "\t", row.names = FALSE,
                     col.names = TRUE, quote = FALSE, na = "NA")
}
summary <- list()
for (geno in c("DNMT3A", "TET2")) {
  panel <- if (geno == "DNMT3A") "Fig3D" else "Fig3E"
  de <- as.data.frame(ref$de.mut2wt.shrink[[1]][[geno]])
  stopifnot(!is.null(rownames(de)), all(c("log2FoldChange", "pvalue", "padj") %in% names(de)))
  genes <- data.frame(gene_id = rownames(de), de, check.names=FALSE, stringsAsFactors=FALSE)
  if (anyDuplicated(genes$gene_id)) stop("Duplicate genes: ", panel)
  genes$volcano_label_candidate <- ifelse(is.finite(genes$pvalue) & genes$pvalue <= 1e-3,
                                          genes$gene_id, NA_character_)
  genes$volcano_raw_pvalue <- genes$pvalue
  write_tsv(genes, file.path(outdir, paste0(panel, "_genes.tsv")))

  sel <- selection[selection$panel == panel & selection$genotype == geno, , drop=FALSE]
  sel <- sel[order(sel$published_order_top_to_bottom), , drop=FALSE]
  n <- nrow(sel)
  if (n == 0L || !identical(as.integer(sel$published_order_top_to_bottom),seq_len(n)) ||
      anyDuplicated(sel$pathway_id)) stop("Invalid pathway selection: ", panel)
  if (!all(sel$published_direction %in% c("positive", "negative"))) stop("Invalid directions: ", panel)
  old <- as.data.frame(if (geno == "DNMT3A") ref$fgsea.dnmt3a.clust[[1]] else ref$fgsea.tet2.clust[[1]])
  old$pathway <- as.character(old$pathway)
  old$pathway_id <- sub(":\\s+.*$", "", old$pathway, perl=TRUE)
  if (anyDuplicated(old$pathway_id)) stop("Duplicate archived IDs: ", panel)
  ix <- match(sel$pathway_id, old$pathway_id)
  if (anyNA(ix)) stop("Missing published pathways: ",panel, ": ", paste(sel$pathway_id[is.na(ix)],collapse=", "))
  fg <- old[ix, , drop=FALSE]
  if (any(!is.finite(fg$pval) | fg$pval <= 0) || any(!is.finite(fg$NES)))
    stop("Non-finite pathway values: ", panel)
  signname <- ifelse(fg$NES > 0, "positive", ifelse(fg$NES < 0,"negative","zero"))
  if (!identical(signname, sel$published_direction)) stop("Enrichment direction mismatch: ",panel)
  pathways <- data.frame(
    pathway = paste(sel$pathway_id, sel$published_label, sep=": "),
    pathway_id = sel$pathway_id,
    archived_pathway = fg$pathway,
    published_order_top_to_bottom = seq_len(n),
    plot_order = n + 1L - seq_len(n),
    size = if ("size" %in% names(fg)) fg$size else rep(NA_integer_,n),
    ES=fg$ES, NES=fg$NES, pval=fg$pval, padj=fg$padj,
    signed_neglog10_raw_pvalue = -log10(fg$pval)*sign(fg$NES),
    enrichment_direction=signname,
    selection_source="Accepted Fig. 3 Illustrator artwork (manually transcribed and checked)",
    stringsAsFactors=FALSE
  )
  if ("log2err" %in% names(fg)) pathways$log2err <- fg$log2err
  write_tsv(pathways, file.path(outdir,paste0(panel,"_pathways.tsv")))
  summary[[panel]] <- data.frame(
    panel=panel, genotype=geno,cell_state="HSC / cluster 1",
    n_genes=nrow(genes), n_label_candidates=sum(!is.na(genes$volcano_label_candidate)),
    n_displayed_pathways=nrow(pathways),
    n_genes_with_finite_pvalue=sum(is.finite(genes$pvalue)),
    n_genes_with_finite_padj=sum(is.finite(genes$padj)),
    input_file="RNA_analysis/data/processed/mutant_de.RData",
    stringsAsFactors=FALSE
  )
  message(panel, ": ",nrow(genes), " genes and ",nrow(pathways)," published pathways")
}
write_tsv(do.call(rbind,summary), file.path(outdir,"Fig3DE_export_summary.tsv"))
writeLines(c(
  "Figure 3D/E canonical Source Data (accepted artwork selection)",
  "Input: RNA_analysis/data/processed/mutant_de.RData",
  "Pathway manifest: RNA_analysis/publication/resources/Fig3DE_published_pathway_ids.tsv",
  "HSC / cluster 1; DNMT3A and TET2 vs WT",
  "Volcano x: DESeq2 normal-shrunken log2FoldChange; y: raw p-value displayed on log10 scale",
  "Volcano label candidate: raw DE p-value <= 0.001",
  "Pathway bar length: sign(NES) * -log10(raw fgsea pval); padj is included for interpretation",
  "Published representative pathway membership/order explicitly preserved",
  "Regenerates Source Data from saved results; no DE or fgsea computations",
  "Export and final figure still require release/legend review"
),file.path(outdir,"Fig3DE_provenance.txt"))
message("Wrote canonical Fig3D/E tables into: ",outdir)
