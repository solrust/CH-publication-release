# Fig. 4K review plot from archived ROC coordinates and AUC summary.
# This script modifies PRESENTATION ONLY. It does not recompute ROC curves,
# AUCs, gene-signature membership, association scores, or Source Data.
# Run from CH_publication_release repository root.

pub <- 'Association_analysis/publication'
source_dir <- file.path(pub, 'source_data')
val_dir <- file.path(pub, 'validation')
fig_dir <- file.path(pub, 'figures')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace('ggplot2', quietly = TRUE)) stop('Package ggplot2 is required')

roc_file <- file.path(source_dir, 'Fig4K_ROC_coordinates.tsv.gz')
auc_file <- file.path(val_dir, 'Fig4K_ROC_AUC_summary.tsv')
for (f in c(roc_file, auc_file)) if (!file.exists(f)) stop('Missing Fig4K input: ', f)
roc <- utils::read.delim(roc_file, check.names = FALSE)
auc <- utils::read.delim(auc_file, check.names = FALSE)

needed_roc <- c('signature', 'phenotype', 'false_positive_rate', 'sensitivity')
needed_auc <- c('signature', 'phenotype', 'auc', 'n_genes')
if (!all(needed_roc %in% names(roc))) stop('ROC input columns missing')
if (!all(needed_auc %in% names(auc))) stop('AUC input columns missing')
if (anyNA(roc[, needed_roc]) || anyNA(auc[, needed_auc])) stop('Missing ROC/AUC values')

signatures <- c('HSC dormancy', 'HSC label-retention', 'HSC output')
phenotypes <- c('growth rate', 'HSC bias', 'clone size')
keys <- expand.grid(signature = signatures, phenotype = phenotypes, stringsAsFactors = FALSE)
key <- function(s,p) paste(s,p,sep='\r')
if (anyDuplicated(key(auc$signature,auc$phenotype))) stop('Duplicate AUC row')
if (!setequal(key(auc$signature,auc$phenotype), key(keys$signature,keys$phenotype)))
  stop('Expected exactly one AUC for each of the 3 signatures x 3 phenotypes')
if (any(!is.finite(auc$auc) | auc$auc < 0 | auc$auc > 1)) stop('Invalid AUC')
if (!setequal(unique(roc$signature), signatures) ||
    !setequal(unique(roc$phenotype), phenotypes)) stop('Unexpected ROC categories')
if (any(!is.finite(roc$false_positive_rate) | !is.finite(roc$sensitivity))) stop('Non-finite ROC points')
if (any(roc$false_positive_rate < 0 | roc$false_positive_rate > 1 |
        roc$sensitivity < 0 | roc$sensitivity > 1)) stop('ROC points outside [0,1]')

roc$signature <- factor(roc$signature, levels = signatures)
roc$phenotype <- factor(roc$phenotype, levels = phenotypes)
auc$signature <- factor(auc$signature, levels = signatures)
auc$phenotype <- factor(auc$phenotype, levels = phenotypes)
# Fixed positions ensure each panel contains its three matching AUC values.
auc$label <- sprintf('%s (AUC = %.2f)', as.character(auc$phenotype), auc$auc)
auc$label_y <- c('growth rate'=0.205, 'HSC bias'=0.135, 'clone size'=0.065)[as.character(auc$phenotype)]
auc$label_x <- 0.978

backdrop <- data.frame(signature = factor(signatures, levels=signatures),
                       xmin = 0.40, xmax = 0.996, ymin = 0.025, ymax = 0.255)
colors <- c('growth rate' = '#EF756C', 'HSC bias' = '#20AF59', 'clone size' = '#65A2E6')
p <- ggplot2::ggplot(roc, ggplot2::aes(x=false_positive_rate, y=sensitivity,
                                       color=phenotype, group=phenotype)) +
  ggplot2::geom_abline(intercept=0, slope=1, linetype='dashed', color='grey60', linewidth=0.45) +
  ggplot2::geom_path(linewidth=0.85) +
  ggplot2::geom_rect(data=backdrop,
                     mapping=ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                     inherit.aes=FALSE, fill='white', color='grey65', linewidth=0.25,
                     alpha=0.94) +
  ggplot2::geom_text(data=auc, ggplot2::aes(x=label_x, y=label_y,
                                         label=label, color=phenotype),
                     inherit.aes=FALSE, hjust=1, size=2.95, show.legend=FALSE) +
  ggplot2::facet_wrap(~signature, nrow=1) +
  ggplot2::scale_color_manual(values=colors, breaks=phenotypes) +
  ggplot2::coord_equal(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
  ggplot2::labs(x='1 - specificity', y='Sensitivity') +
  ggplot2::theme_bw(base_size=12) +
  ggplot2::theme(legend.position='none', strip.text=ggplot2::element_text(size=11))

out <- file.path(fig_dir, 'Fig4K_reproduced_with_AUC.pdf')
ggplot2::ggsave(filename=out, plot=p, width=12.8, height=4.5, useDingbats=FALSE)
cat('Created ', out, '\n', sep='')
print(auc[,c('signature','phenotype','auc','n_genes')], row.names=FALSE)
cat('The accepted Fig. 4K already includes these rounded AUCs;\n',
    'this plot is a clearer REPRODUCTION, not a numerical correction.\n',sep='')
