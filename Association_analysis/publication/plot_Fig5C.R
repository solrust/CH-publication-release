# Plot Fig. 5C from the exported table only.
# Run from the CH_publication_release repository root.
if (!requireNamespace('ggplot2', quietly=TRUE)) stop('Install ggplot2')
base <- 'Association_analysis/publication'
f <- file.path(base,'source_data','Fig5C_gene_association_counts.tsv')
if (!file.exists(f)) stop('Run export_Fig5C.R first')
df <- read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
if (!all(c('phenotype','term','n_genes_BH_FDR_le_0p5') %in% names(df))) stop('Unexpected count table columns')
phen.order <- names(sort(tapply(df$n_genes_BH_FDR_le_0p5, df$phenotype, sum),decreasing=TRUE))
df$phenotype <- factor(df$phenotype,levels=rev(phen.order))
df$term <- factor(df$term, levels=c('phenotype','phenotype:DNMT3A','phenotype:TET2'))
p <- ggplot2::ggplot(df,ggplot2::aes(x=n_genes_BH_FDR_le_0p5,y=phenotype,fill=term)) +
  ggplot2::geom_col(position=ggplot2::position_dodge(width=0.8), width=0.75,color='black',linewidth=0.2) +
  ggplot2::scale_fill_manual(values=c('phenotype'='#999999','phenotype:DNMT3A'='#E69F00','phenotype:TET2'='#56B4E9')) +
  ggplot2::theme_bw(base_size=13) +
  ggplot2::labs(x='Significant gene associations (BH FDR <= 0.5)',y='Clone phenotype',fill=NULL)
dir.create(file.path(base,'figures'),recursive=TRUE,showWarnings=FALSE)
ggplot2::ggsave(file.path(base,'figures','Fig5C_reproduced_review.pdf'),plot=p,width=9,height=5.5)
message('Created Fig5C review plot from Source Data.')
