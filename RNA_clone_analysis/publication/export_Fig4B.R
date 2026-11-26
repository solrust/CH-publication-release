# Figure 4B: export only the three highlighted clone cells.
# Reuse the canonical UMAP background from Figure 3A Source Data.
# Run from the CH_publication_release repository root:
#   source('RNA_clone_analysis/publication/export_Fig4B.R')
# No scrna.RData load or UMAP recomputation is required.

fig3a_file <- 'RNA_analysis/publication/source_data/Fig3A_umap_cells.tsv.gz'
clone_file <- 'RNA_clone_analysis/data/processed/clone_info.RData'
outdir <- 'RNA_clone_analysis/publication/source_data'
valdir <- 'RNA_clone_analysis/publication/validation'
for (f in c(fig3a_file, clone_file)) {
  if (!file.exists(f)) stop('Missing required input: ', f, '. Run from repository root.')
}
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(valdir, recursive=TRUE, showWarnings=FALSE)

# Support a small number of conventional coordinate column names. Fail loudly otherwise.
fig4b_find_umap_cols <- function(dat) {
  norm <- tolower(gsub('[^[:alnum:]]', '', names(dat)))
  get_one <- function(choices, purpose) {
    i <- match(choices, norm, nomatch=0L)
    i <- i[i > 0L]
    if (!length(i)) {
      stop('Cannot identify ', purpose, ' column in Fig3A_umap_cells.tsv.gz. ',
           'Available columns: ', paste(names(dat), collapse=', '))
    }
    names(dat)[i[1L]]
  }
  list(
    cell_id=get_one(c('cellid', 'cellbarcode', 'barcode', 'cell', 'cellname'), 'cell ID'),
    x=get_one(c('umap1', 'umapx', 'embedding1', 'emb1', 'dimension1', 'dim1', 'x'), 'UMAP-1'),
    y=get_one(c('umap2', 'umapy', 'embedding2', 'emb2', 'dimension2', 'dim2', 'y'), 'UMAP-2')
  )
}

message('Reading Fig. 3A Source Data as the shared UMAP coordinate reference.')
umap_raw <- utils::read.delim(fig3a_file, check.names=FALSE, stringsAsFactors=FALSE)
cols <- fig4b_find_umap_cols(umap_raw)
umap <- data.frame(
  cell_id=as.character(umap_raw[[cols$cell_id]]),
  UMAP_1=suppressWarnings(as.numeric(umap_raw[[cols$x]])),
  UMAP_2=suppressWarnings(as.numeric(umap_raw[[cols$y]])),
  stringsAsFactors=FALSE
)
if (anyNA(umap$cell_id) || any(!nzchar(umap$cell_id)) || anyDuplicated(umap$cell_id)) {
  stop('Fig. 3A cell IDs are missing, empty or duplicated.')
}
if (any(!is.finite(umap$UMAP_1) | !is.finite(umap$UMAP_2))) {
  stop('Fig. 3A UMAP coordinates contain missing or nonfinite values.')
}

cl <- new.env(parent=baseenv())
load(clone_file, envir=cl)
if (!all(c('cell.clone','clone.info') %in% ls(cl))) {
  stop('clone_info.RData must contain cell.clone and clone.info.')
}
if (!'cloneid' %in% names(cl$cell.clone) || is.null(rownames(cl$cell.clone))) {
  stop('cell.clone must have cell IDs as row names and a cloneid field.')
}

clone_numbers <- c(16L, 213L, 396L)
clone_descriptions <- c('myeloid biased', 'lymphoid biased', 'HSC biased')
if (nrow(cl$clone.info) < max(clone_numbers)) stop('clone.info lacks expected clone rows.')
clone_ids <- rownames(cl$clone.info)[clone_numbers]
if (anyNA(clone_ids) || anyDuplicated(clone_ids)) stop('Invalid example clone IDs.')

assigned <- data.frame(
  cell_id=rownames(cl$cell.clone),
  clone_id=as.character(cl$cell.clone$cloneid),
  stringsAsFactors=FALSE
)
assigned <- assigned[assigned$clone_id %in% clone_ids, , drop=FALSE]
if (anyDuplicated(assigned$cell_id)) stop('Selected cell IDs appear more than once in cell.clone.')
assigned$clone_number <- clone_numbers[match(assigned$clone_id, clone_ids)]
# Keep only plotted cells: identifiers must be in the canonical Fig. 3A UMAP.
where <- match(assigned$cell_id, umap$cell_id)
visible <- !is.na(where)
selected <- assigned[visible, , drop=FALSE]
if (!nrow(selected)) stop('No selected clone cells match Fig. 3A cell IDs; check ID conventions.')
selected$UMAP_1 <- umap$UMAP_1[where[visible]]
selected$UMAP_2 <- umap$UMAP_2[where[visible]]
selected$genotype <- as.character(cl$clone.info$genotype[clone_numbers[match(selected$clone_id, clone_ids)]])
selected$phenotype_label <- clone_descriptions[match(selected$clone_id, clone_ids)]
selected$reported_clone_size <- as.integer(cl$clone.info$size[clone_numbers[match(selected$clone_id, clone_ids)]])
selected <- selected[order(selected$clone_number, match(selected$cell_id, umap$cell_id)),
                     c('cell_id','UMAP_1','UMAP_2','clone_number','clone_id','genotype','phenotype_label','reported_clone_size'),
                     drop=FALSE]

summary <- data.frame(
  clone_number=clone_numbers,
  phenotype_label=clone_descriptions,
  clone_id=clone_ids,
  genotype=as.character(cl$clone.info$genotype[clone_numbers]),
  n_cells_reported=as.integer(cl$clone.info$size[clone_numbers]),
  n_cells_assigned_in_checkpoint=as.integer(vapply(clone_ids, function(id) sum(assigned$clone_id == id), integer(1))),
  n_cells_on_Fig3A_umap=as.integer(vapply(clone_ids, function(id) sum(selected$clone_id == id), integer(1))),
  stringsAsFactors=FALSE
)
summary$n_assigned_not_on_umap <- summary$n_cells_assigned_in_checkpoint - summary$n_cells_on_Fig3A_umap

outfile <- file.path(outdir, 'Fig4B_example_clone_cells.tsv.gz')
con <- gzfile(outfile, open='wt')
tryCatch(utils::write.table(selected, con, sep='\t', quote=FALSE, row.names=FALSE),
         finally=close(con))
utils::write.table(summary, file.path(valdir,'Fig4B_clone_summary.tsv'),
                   sep='\t', quote=FALSE, row.names=FALSE)
print(summary[, c('clone_number','genotype','n_cells_reported',
                  'n_cells_assigned_in_checkpoint','n_cells_on_Fig3A_umap')],
      row.names=FALSE)
message('Wrote ', nrow(selected), ' highlighted cell rows to ', outfile)
if (any(summary$n_cells_reported != summary$n_cells_on_Fig3A_umap)) {
  message('NOTE: Some reported clone cells are not shown on the Fig. 3A embedding. ',
          'Do not substitute reported clone sizes for plotted cell counts.')
}
