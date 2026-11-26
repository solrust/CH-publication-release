library(fgsea)
library(parallel)
source('R/gsea_clone_summary.R')
gene_sets <- readRDS('../RNA_analysis/data/genesets/processed/gene_sets_use.rds')

### choose folder to save the results
save_folder <- 'default'
data_folder <- paste('data/processed/', save_folder, '/', sep = '')
if (!file.exists(data_folder)) dir.create(path_folder)

###### test pathway-mediated association of mutation with clone phenotypes
pheno_names_new <- c('growth rate', 'clone size', 'HSC bias', 'HSPC bias', 'lymphoid bias', 'dendritic bias', 'pre-B bias', 'myeloid bias', 'erythroid bias')
pheno_names_new <- c('growth rate', 'clone size', 'HSC bias', 'lymphoid bias', 'dendritic bias', 'myeloid bias', 'erythroid bias')
test.basic.list <- lapply(#colnames(clone.pheno),
                          pheno_names_new,
                          function(nm){
    fl <- paste(data_folder, 'associations_', nm, '.rds', sep = '')
    test.basic.nm <- readRDS(fl)
    test.basic.nm[[1]]
  })
names(test.basic.list) <- pheno_names_new#[c(1,2,3,5,6,8,9)]
#
test.inter.list <- lapply(#colnames(clone.pheno),
                          pheno_names_new,
                          function(nm){
    fl <- paste(data_folder, 'associations_', nm, '.rds', sep = '')
    test.basic.nm <- readRDS(fl)
    test.basic.nm[[2]]
})
names(test.inter.list) <- pheno_names_new[c(1,2,3,5,6,8,9)]
### estimate gene set enrichment of clone-specific DE
gsea.clone <- gsea_clone_summary(test.basic.list, gene_sets)

# prioritize non-redundant pathways
mainPathways.clones <- mclapply(1:length(test.basic.list), function(cl){
  print(cl)
  gsc <- gsea.clone$gscore.clone[,cl]
  fgseaRes <- gsea.clone$gsea.clone[[cl]]
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
                                        gene_sets[ nchar(names(gene_sets)) < 60],
                                        gsc[!is.na(gsc)],
                                        pval.threshold = 0.001)
  mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
  return(mainPathways)
}, mc.cores = length(test.basic.list))
str(mainPathways.clones)

fl <- paste(data_folder, 'association_summary.RData', sep = '')
save(test.basic.list, test.inter.list, gsea.clone, mainPathways.clones,
     #save_folder,
     file = fl)

