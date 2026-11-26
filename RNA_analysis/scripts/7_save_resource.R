library(openxlsx)
#
ll <- lapply(1:length(de.mut2wt.shrink), function(cl)
  data.frame(de.mut2wt.shrink[[cl]]$DNMT3A, DE = ifelse(!is.na(de.mut2wt.shrink[[cl]]$DNMT3A$padj) & de.mut2wt.shrink[[cl]]$DNMT3A$padj < 0.1, 'yes', 'no') ))
names(ll) <- paste('cluster', 1:length(ll))
write.xlsx(ll, '../resource/DNMT3A_DE.xlsx', row.names = TRUE, col.names = TRUE)
#
ll <- lapply(1:length(de.mut2wt.shrink), function(cl)
  data.frame(de.mut2wt.shrink[[cl]]$TET2, DE = ifelse(!is.na(de.mut2wt.shrink[[cl]]$TET2$padj) & de.mut2wt.shrink[[cl]]$TET2$padj < 0.1, 'yes', 'no') ))
names(ll) <- paste('cluster', 1:length(ll))
write.xlsx(ll, '../resource/TET2_DE.xlsx', row.names = TRUE, col.names = TRUE)


### save metadata
load(file = 'data/processed/scrna.RData')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
hash.match <- readRDS(file = "data/processed/hashtag_meta.rds")

meta <- cbind(hash.match[sampleid, c('exp', 'genotype')], sample = sampleid)
rownames(meta) <- names(sampleid)
meta$clusters <- NA
meta[names(clust.final),]$clusters <- clust.final
meta[,c('emb1', 'emb2')] <- NA
meta[rownames(emb), c('emb1', 'emb2')] <- emb
meta$clone <- NA
meta[rownames(cell.clone),]$clone <- cell.clone$cloneid
write.csv(meta, file = 'data/processed/metainfo.csv', quote = FALSE)

### save cluster-level summary table
library(dplyr)
library(tidyr)
meta <- read.csv(file = 'data/processed/metainfo.csv')
scrna.env <- new.env()
load(file = 'data/processed/scrna.RData', envir = scrna.env)
clust.labels <- c(`1`="HSC", `2`="LMPP-1", `3`="MLP", `4`="Pro-B", `5`="pDC",
                  `6`="LMPP-2", `7`="early GMP", `8`="late GMP", `9`="EMPP/MkP", `10`="EryP")

meta$cluster_final <- scrna.env$clust.final[meta$X]
meta$population <- clust.labels[as.character(meta$cluster_final)]

cell.summary <- meta %>%
  filter(!is.na(cluster_final)) %>%
  count(exp, sample, genotype, cluster_final, population, name="n_cells") %>%
  complete(nesting(exp, sample, genotype),
           cluster_final=as.numeric(names(clust.labels)),
           fill=list(n_cells=0)) %>%
  mutate(population=clust.labels[as.character(cluster_final)]) %>%
  group_by(exp, sample, genotype) %>%
  mutate(total_cells_sample=sum(n_cells),
         fraction_of_sample=n_cells/total_cells_sample) %>%
  ungroup() %>%
  arrange(exp, sample, cluster_final)

head(cell.summary)

cell.summary <- cell.summary %>%
  select(exp, sample, genotype, cluster_final, population,
         n_cells, total_cells_sample, fraction_of_sample)

write.csv(cell.summary, "../resource/scRNA_cell_counts.csv", row.names=FALSE)
