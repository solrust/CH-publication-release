source('R/association_routines.R')
source('R/visualization_routines.R')
library(parallel)
library(Matrix)
library(MASS)

load('../RNA_analysis/data/processed/scrna.RData')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
load('../DNA_analysis/data/processed/DNA_clones.RData')

### estimate clone-specific count matrices for each leiden cluster
clusts <- p2clusts$PCAncc$leiden
clone.counts <- mclapply(1:length(unique(clusts)), function(cl){
  cls <- intersect(rownames(cell.clone),
                   names(clusts)[clusts == cl])
  str(cls)
  x <- aggregate(t(as.matrix(cd[, cls])),
                 by = list(cell.clone[cls,]$cloneid), sum)

  rownames(x) <- x[,1]
  x <- x[,-1]
  x <- Matrix::Matrix(t(as.matrix(x)))
  str(x)
  return(x)
}, mc.cores = 10)

### specify matrix of clone phenotypes:
clons <- rownames(clone.states2)
summary.size <- aggregate(clone.info[clons,]$size, by = list(clone.info[clons,]$id), sum)
rownames(summary.size) <- summary.size[,1]; summary.size <- summary.size[,-1,drop = F]
clone.pheno <- cbind('growth rate' = rate[clons],
                     'size' = log(clone.info[clons,]$size / summary.size[clone.info[clons,]$id,]),
                     clone.states2[clons, c(1, 4, 5, 8, 10)]
)
rownames(clone.pheno) <- clons
colnames(clone.pheno)[-c(1:2)] <-  c('HSPC', 'Lymphoid', 'DC', 'Myeloid', 'Erythroid') #setdiff(coarse_ann, 'pre-B')
#
pheno_names_new <- c('growth rate', 'clone size', 'HSC bias', 'HSPC bias', 'lymphoid bias', 'dendritic bias', 'pre-B bias', 'myeloid bias', 'erythroid bias')
colnames(clone.pheno) <-  pheno_names_new[c(1,2,3,5,6,8,9)]
head(clone.pheno)


# fix count matrix
#counts <- old.env$clone.counts[[11]]
counts <- clone.counts[[3]]
#counts <- old.env$clone.counts.hspc[[3]]
# remove low covered clones
clons <- intersect(colnames(counts), rownames(clone.info))
str(clons)
clons <- clons[clone.info[clons, ]$size > 50 &  colSums(counts[,clons]) > 0 ]
str(clons)
clons <- intersect(clons, names(rate))
str(clons)
zrate <- (rate[clons] - mean(rate[clons])) / sd(rate[clons])
clons <- clons[abs(zrate) < 1] # remove rate-outlier clones
#clons <- clons[rate[clons] >= 0]
str(clons)
# remove lowly-expressed genes
counts <- counts[, clons]
exn <- rowSums(counts > 0)
counts <- counts[exn >= 20, ]
fpm <- t(t(counts) / colSums(counts))
dim(fpm)


save_folder <- 'default'
data_folder <- paste('data/processed/', save_folder, '/', sep = '')
figure_folder <- paste('figures/', save_folder, '/', sep = '')
if (!file.exists(data_folder)) dir.create(path_folder)
if (!file.exists(figure_folder)) dir.create(figure_folder)

for (i in 1:ncol(clone.pheno)){
  #i <- 1
  pheno.state <- clone.pheno[clons,i]
  pheno.name <- colnames(clone.pheno)[i]
  str(pheno.state)
  str(pheno.name)

  ### run basic association model
  ind <- TRUE
  test.basic <- fit.nb(t(counts[,clons][,ind]), clone.info[clons,][ind,], pheno.state[ind], #clusts = cls,
                       do.batch = T, do.interactions = FALSE, nc.cores = 10)
  colnames(test.basic) <- c('pheno', 'exp', 'DNMT3A', 'TET2', 'p.pheno', 'p.exp', 'p.DNMT3A', 'p.TET2')
  sum(p.adjust(test.basic[,5], 'BH') < 0.2); sum(p.adjust(test.basic[,6], 'BH') < 0.2); sum(p.adjust(test.basic[,7], 'BH') < 0.2)
  sum(p.adjust(test.basic[,5], 'BH') < 0.5); sum(p.adjust(test.basic[,6], 'BH') < 0.5); sum(p.adjust(test.basic[,7], 'BH') < 0.5)
  de <- test.basic[,1]; pv <- test.basic[,5]
  plot_associations(de, pv, pheno.name)

  ### run association model with mutation-interactions
  test.interactions <- fit.nb(t(counts[ind,clons]), clone.info[clons,][ind,], pheno.state[ind], #clusts = cls,
                             do.batch = TRUE, do.interactions = TRUE, nc.cores = 10)
  colnames(test.interactions) <- c('pheno','DNMT3A', 'TET2', 'exp', 'pheno:DNMT3A', 'pheno:TET2', 'p.pheno','p.DNMT3A', 'p.TET2', 'p.exp', 'p.pheno:DNMT3A', 'p.pheno:TET2')
  sum(p.adjust(test.interactions[,7], 'BH') < 0.2, na.rm = T); sum(p.adjust(test.interactions[,8], 'BH') < 0.2, na.rm = T); sum(p.adjust(test.interactions[,9], 'BH') < 0.2, na.rm = T); sum(p.adjust(test.interactions[,10], 'BH') < 0.2, na.rm = T);
  sum(p.adjust(test.interactions[,11], 'BH') < 0.2, na.rm = T); sum(p.adjust(test.interactions[,12], 'BH') < 0.2, na.rm = T)

  fl <- paste(data_folder, 'associations_', pheno.name,'.rds', sep = '')
  saveRDS(list(test.basic, test.interactions), file = fl)

  ### select significant genes (associated with pheno or pheno::geno)
  bh.cutoff <- 0.5
  pdj.dnmt3 <- p.adjust(test.interactions[,11], 'BH')
  pdj.tet2 <- p.adjust(test.interactions[,12], 'BH')
  ind <- p.adjust(test.basic[,5], 'BH') < bh.cutoff |
         p.adjust(test.interactions[,11], 'BH') < bh.cutoff |
         p.adjust(test.interactions[,12], 'BH') < bh.cutoff
  sum(ind)
  gg0 <- names(sort(test.basic[ind,5]))
  str(gg0)

  ### save scatterplot of gene associations for significant genes
  fl <- paste(figure_folder, 'association_', pheno.name, '_genes' ,'.pdf', sep = '')
  pdf(fl, width = 6, height = 4)
  for(gn in gg0[1:min(100, length(gg0))] ){
    print(gn)
    #gn <- names(sort(pv))[1]
    print(show_gene(as.matrix(t(fpm[,clons])), gn,  pheno.state, clone.info[clons,], sign.clones = TRUE,
                    xname = pheno.name, title = genes.ital[[gn]], text.size = 20))
  }
  dev.off()
}

