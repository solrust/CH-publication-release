library(DESeq2)
library(ggrepel)

tbl <- read.table('data/raw/Exp3_Count_all_HSC.txt', header = T)

# comments:
# examine HSC/bulk ratio; use edge/deseq to assess mutant to control stat significance.
# for each TF: logFC of HSC increase + p-value
# comparison with other data?

counts <- tbl[,-c(1,2)]
rownames(counts) <- tbl$ORFID
metadat <- data.frame(colnames(counts)); colnames(metadat) <- 'id'
metadat$genotype <- 'WT'
ind <- grepl('D', metadat$id)
metadat$genotype[ind] <- 'DNMT3A' #  unlist(lapply(strsplit(metadat$id[ind], '_'), function(x) x[2]))
ind <- grepl('T', metadat$id)
metadat$genotype[ind] <- 'TET2' #  unlist(lapply(strsplit(metadat$id[ind], '_'), function(x) x[2]))
metadat$genotype <- factor(metadat$genotype, levels = c("WT", "DNMT3A", "TET2"))
metadat$state <- 'bulk'
metadat$state[grepl('HSC', metadat$id)] <- 'HSC'
dds_mut <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadat,
                                  design = ~ genotype * state)#~ genotype + state + geontype:state)
dds_mut <- DESeq(dds_mut)
resultsNames(dds_mut)
resd <- lfcShrink(dds_mut,
          #coef = "genotype_TET2_vs_WT",
          coef = "state_HSC_vs_bulk",
          #coef = "genotypeTET2.stateHSC",
          #coef = "genotypeDNMT3A.stateHSC",
          type="apeglm"
          #type="normal"
)
resd$TF <- tbl$TF
#head(resd)
ggplot(resd, aes(log2FoldChange, pvalue, label = TF)) + geom_point(size = 4, shape = 21, fill = 'grey') +
  scale_y_log10() + geom_hline(yintercept = 0.1, linetype = 2) +
  geom_text_repel(
    #data = genes_to_label,
    #aes(label = TF),
    size = 5,                  # Text size
    box.padding = 0.5,         # How far to stay away from neighboring text boxes
    point.padding = 0.3,       # How far to stay away from the actual data point
    segment.color = "black",   # Color of the line connecting text to point
    segment.size = 0.5,        # Thickness of the connector line
    max.overlaps = 20          # Excludes text if it's too crowded (set to Inf to force all)
  ) + theme_bw()
head(resd[order(resd$pvalue),])
resd[resd$TF %in% c('IL18', 'RNF24'),]


env.mut <- new.env()
load(file = '../RNA_analysis/data/processed/mutant_de.RData', envir = env.mut)
head(env.mut$fgsea.tet2.score[,1])
ind <- grepl('IL', rownames(env.mut$fgsea.tet2.score))
table(ind)
env.mut$fgsea.tet2.score[ind,]

ind <- grepl('IL', rownames(env.mut$fgsea.dnmt3a.score))
table(ind)
env.mut$fgsea.dnmt3a.score[ind,]
