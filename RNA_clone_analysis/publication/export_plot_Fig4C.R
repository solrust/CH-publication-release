# Figure 4C: clone phenotype distributions (publication correction)
# Run from the CH_publication_release repository root.
# Input: RNA_clone_analysis/data/processed/clone_info.RData
# This script deliberately does not source the historical full analysis.
# Correction: HSC bias = state 1, myeloid bias = state 8, lymphoid bias = state 4.

input <- file.path('RNA_clone_analysis', 'data', 'processed', 'clone_info.RData')
if (!file.exists(input)) stop('Missing input: ', input, '; run from the repository root.')

outroot <- file.path('RNA_clone_analysis', 'publication')
src_dir <- file.path(outroot, 'source_data')
fig_dir <- file.path(outroot, 'figures')
val_dir <- file.path(outroot, 'validation')
for (d in c(src_dir, fig_dir, val_dir)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

e <- new.env(parent = emptyenv())
load(input, envir=e)
required <- c('clone.info', 'clone.states', 'clone.states2')
if (!all(required %in% ls(e))) stop('Missing objects: ', paste(setdiff(required, ls(e)), collapse=', '))
if (!all(c('1', '4', '8') %in% colnames(e$clone.states2)))
  stop('Missing state columns 1, 4, or 8 in clone.states2')

# Match the historical panel's >=20-cells-per-clone selection.
keep <- rownames(e$clone.states)[rowSums(e$clone.states) >= 20]
if (!all(keep %in% rownames(e$clone.states2)) || !all(keep %in% rownames(e$clone.info)))
  stop('Clone IDs are inconsistent across the saved objects.')

src <- data.frame(
  clone_id = keep,
  genotype = as.character(e$clone.info[keep, 'genotype']),
  clone_size = as.numeric(e$clone.info[keep, 'size']),
  myeloid_bias_state8 = as.numeric(e$clone.states2[keep, '8']),
  lymphoid_bias_state4 = as.numeric(e$clone.states2[keep, '4']),
  HSC_bias_state1 = as.numeric(e$clone.states2[keep, '1']),
  stringsAsFactors = FALSE
)

stopifnot(nrow(src) > 0, !anyDuplicated(src$clone_id))
counts <- as.data.frame(table(factor(src$genotype, levels=c('WT','DNMT3A','TET2'))), stringsAsFactors=FALSE)
colnames(counts) <- c('genotype', 'n_clones')

write.table(src, file.path(src_dir, 'Fig4C_clone_phenotypes.tsv'),
            sep='\t', quote=FALSE, row.names=FALSE)
write.table(counts, file.path(val_dir, 'Fig4C_clone_counts.tsv'),
            sep='\t', quote=FALSE, row.names=FALSE)

# Numerical size audit for the three examples printed in accepted Fig. 4B.
example_indices <- c(16L, 213L, 396L)
examples <- data.frame(
  figure_label = paste('clone', example_indices),
  clone_index_in_saved_object = example_indices,
  clone_id = rownames(e$clone.info)[example_indices],
  genotype = as.character(e$clone.info[example_indices, 'genotype']),
  n_cells = as.integer(e$clone.info[example_indices, 'size']),
  stringsAsFactors=FALSE
)
write.table(examples, file.path(val_dir, 'Fig4B_example_clone_sizes.tsv'),
            sep='\t', quote=FALSE, row.names=FALSE)

cat('Fig4C clones by genotype:\n')
print(counts, row.names=FALSE)
cat('Fig4B example clone sizes:\n')
print(examples[, c('figure_label','genotype','n_cells')], row.names=FALSE)
cat('Fig4C HSC and myeloid biases differ for ',
    sum(src$HSC_bias_state1 != src$myeloid_bias_state8, na.rm=TRUE),
    ' of ', nrow(src), ' clones.\n', sep='')

if (!requireNamespace('ggplot2', quietly=TRUE)) {
  message('Package ggplot2 unavailable: Source Data exported; skipping plot.')
} else {
  n <- nrow(src)
  plot_data <- data.frame(
    variable = factor(rep(c('Size','My bias','Ly bias','HSC bias'), each=n),
                      levels=c('Size','My bias','Ly bias','HSC bias')),
    value = c(src$clone_size, src$myeloid_bias_state8,
              src$lymphoid_bias_state4, src$HSC_bias_state1),
    genotype = factor(rep(src$genotype, times=4), levels=c('WT','DNMT3A','TET2')),
    stringsAsFactors=FALSE
  )
  pl <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = "", y = value)
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      color = "black"
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(fill = genotype),
      width = 0.4,
      height = 0,
      seed = 239,
      shape = 21,
      size = 2,
      stroke = 0.5,
      color = "black"
    ) +
    ggplot2::facet_wrap(
      ~variable,
      scales = "free_y",
      nrow = 1
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        WT = adjustcolor("black", alpha.f = 0.5),
        DNMT3A = adjustcolor("red", alpha.f = 0.5),
        TET2 = adjustcolor("green", alpha.f = 0.5)
      )
    ) +
    ggplot2::labs(
      x = "",
      y = "clonal phenotype",
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_blank()
    )
  outfile <- file.path(fig_dir, 'Fig4C_reproduced_corrected.pdf')
  ggplot2::ggsave(outfile, plot=pl, width=9, height=4, units='in', device='pdf')
  cat('Wrote ', outfile, '\n', sep='')
}
cat('Fig4C export complete. Source data: ', nrow(src), ' clones.\n', sep='')
