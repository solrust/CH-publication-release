### select DE gene set
cl <- 1
tbl <- de.env$de.mut2wt.shrink[[cl]]$TET2
gset <- rownames(tbl)[!is.na(tbl$padj) & tbl$padj < 0.2 & abs(tbl$log2FoldChange) > 0.1]
str(gset)

### identify neighboring peaks
peaks_gr <- granges(obj.combined)
annotations <- Annotation(obj.combined)
promoters <- promoters(annotations, upstream = 1000, downstream = 1000)

# Filter the tss_window object first
specific_tss_window <- promoters[promoters$gene_name %in% gset]
specific_overlaps <- findOverlaps(peaks_gr, specific_tss_window)
specific_peaks <- peaks_gr[unique(queryHits(specific_overlaps))]
length(specific_peaks)


