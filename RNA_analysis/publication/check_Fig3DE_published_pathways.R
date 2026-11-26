#!/usr/bin/env Rscript
# Compare pathway labels in accepted Illustrator Fig. 3D/E with archived
# fgsea results and the initial published Source Data candidates.
# No differential expression or fgsea calculations are performed.
# Run from CH_publication_release root.

root <- normalizePath(getwd(), mustWork = TRUE)
publication_dir <- file.path(root, "RNA_analysis", "publication")
selection_file <- file.path(publication_dir, "resources", "Fig3DE_published_pathway_ids.tsv")
stopifnot(file.exists(selection_file))
selected <- read.delim(selection_file, stringsAsFactors = FALSE,
                       check.names = FALSE)
stopifnot(!anyDuplicated(paste(selected$panel, selected$pathway_id)),
          all(selected$published_direction %in% c("positive", "negative")))

input_root <- Sys.getenv("CH_DATA_ROOT", unset = root)
input_file <- file.path(input_root, "RNA_analysis", "data", "processed", "mutant_de.RData")
if (!file.exists(input_file)) stop("Missing file: ", input_file)
if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 package not available")

archived <- new.env(parent = globalenv())
load(input_file, envir = archived)
required <- c("fgsea.dnmt3a.clust", "fgsea.tet2.clust")
if (!all(vapply(required, exists, logical(1), envir=archived, inherits=FALSE))) {
  stop("Archived object is missing fgsea results")
}

validation_dir <- file.path(publication_dir, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

check_rows <- list()
for (panel in c("Fig3D", "Fig3E")) {
  genotype <- if (panel == "Fig3D") "DNMT3A" else "TET2"
  rows <- selected[selected$panel == panel, , drop = FALSE]
  rows <- rows[order(rows$published_order_top_to_bottom), , drop = FALSE]
  fg <- if (genotype == "DNMT3A") archived$fgsea.dnmt3a.clust[[1]] else archived$fgsea.tet2.clust[[1]]
  fg <- as.data.frame(fg)
  fg$pathway <- as.character(fg$pathway)
  # Match by stable GO/MSigDB/Reactome IDs, not potentially revised names.
  fg$pathway_id <- sub(":[[:space:]].*$", "", fg$pathway)
  if (anyDuplicated(fg$pathway_id)) stop("Duplicate pathway IDs for ", genotype)
  i <- match(rows$pathway_id, fg$pathway_id)

  old_file <- file.path(publication_dir, "source_data", paste0(panel, "_pathways.tsv"))
  if (file.exists(old_file)) {
    old <- read.delim(old_file, stringsAsFactors=FALSE)
    old$pathway_id <- sub(":[[:space:]].*$", "", as.character(old$pathway))
    rows$included_in_first_export <- rows$pathway_id %in% old$pathway_id
    count_first <- nrow(old)
  } else {
    rows$included_in_first_export <- NA
    count_first <- NA_integer_
  }

  rows$archived_pathway <- fg$pathway[i]
  rows$ES <- fg$ES[i]
  rows$NES <- fg$NES[i]
  rows$pval <- fg$pval[i]
  rows$padj <- fg$padj[i]
  rows$signed_neglog10_raw_pvalue <- -log10(rows$pval) * sign(rows$NES)
  rows$found_in_archived_fgsea <- !is.na(i)
  rows$direction_matches_accepted <- !is.na(rows$NES) &
    ifelse(rows$NES > 0, "positive", "negative") == rows$published_direction
  check_rows[[panel]] <- rows
  message(panel, ": ", nrow(rows), " pathways in accepted figure; ",
          count_first, " pathways in first export; ",
          sum(!rows$found_in_archived_fgsea), " IDs absent from archived fgsea; ",
          sum(!rows$direction_matches_accepted), " direction mismatches")
}

out <- do.call(rbind, check_rows)
outfile <- file.path(validation_dir, "Fig3DE_published_pathway_check.tsv")
write.table(out, outfile, sep="\t", quote=FALSE,
            row.names=FALSE, na="NA")
message("Wrote: ", outfile)
if (any(!out$found_in_archived_fgsea)) {
  warning("At least one final-figure pathway was not found in the archived results")
}
if (any(!out$direction_matches_accepted)) {
  warning("At least one final-figure pathway has an unexpected enrichment direction")
}
