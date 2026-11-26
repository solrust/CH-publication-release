library(ggplot2)
library(scales)
library(RColorBrewer)

file_path <- '../Association_analysis/data/processed/use/association_summary.RData'
ass.env <- new.env()
load(file_path, envir = ass.env)
env.dna_clones <- new.env()
load('../DNA_analysis/data/processed/DNA_clones.RData', envir = env.dna_clones)

### examine expression of MYC and MYC targets across HSC subpopulations.
hsc.env <- new.env()
load('../Association_analysis/data/processed/hsc/hsc_axes.RData', env = hsc.env)
rna.env <- new.env()
load('../RNA_analysis/data/processed/scrna.RData', env = rna.env)


### revisit correlation of HSC phenotypes and hscPCs
env.clones <- new.env()
load('../RNA_clone_analysis/data/processed/clone_info.RData', envir = env.clones)
head(env.clones$cell.clone)
str(hsc.env$prr$x)
head(env.clones$clone.info)
head(env.clones$clone.states2)


### specify matrix of clone phenotypes:
clons <- rownames(env.clones$clone.states2)
summary.size <- aggregate(env.clones$clone.info[clons,]$size, by = list(env.clones$clone.info[clons,]$id), sum)
rownames(summary.size) <- summary.size[,1]; summary.size <- summary.size[,-1,drop = F]
clone.pheno <- cbind('growth rate' = env.dna_clones$rate[clons],
                     'size' = log(env.clones$clone.info[clons,]$size / summary.size[env.clones$clone.info[clons,]$id,]),
                     env.clones$clone.states2[clons, c(1, 4, 5, 8, 10)]
)
rownames(clone.pheno) <- clons
colnames(clone.pheno)[-c(1:2)] <-  c('HSPC', 'Lymphoid', 'DC', 'Myeloid', 'Erythroid') #setdiff(coarse_ann, 'pre-B')
#
pheno_names_new <- c('growth rate', 'clone size', 'HSC bias', 'HSPC bias', 'lymphoid bias', 'dendritic bias', 'pre-B bias', 'myeloid bias', 'erythroid bias')
colnames(clone.pheno) <-  pheno_names_new[c(1,2,3,5,6,8,9)]
head(clone.pheno)
clone.pheno <- clone.pheno[env.clones$cell.clone$cloneid,]
rownames(clone.pheno) <- rownames(env.clones$cell.clone)
clone.pheno <- apply(clone.pheno, 2, function(x) rank(x)/length(x))


# choose HSCs
cls <- intersect(rownames(hsc.env$prr$x), rownames(clone.pheno))
str(cls)

# smooth clone phenotypes
ind.pc <- 1:10
ft <- lm((clone.pheno[cls, ]) ~ hsc.env$prr$x[cls, ind.pc])
clone.pheno.smoothed <- ft$fitted.values


png('figures/hsc/hsc_clonal_density.png', width = 3 * 150, height = 3 * 150)
ii <- 1:ncol(clone.pheno.smoothed)
nms <- colnames(clone.pheno)
par(mfrow=c(3, 3), mar = rep(2, 4), oma = c(0, 0, 3, 0))
for (i in 1:length(ii)){
  scor <- clone.pheno.smoothed[,ii[i]]
  cols <- colorRampPalette(c("grey95", "darkred" ))(n = 10)[as.numeric(cut(rank(scor), 10))]
  plot(rna.env$emb[cls,], pch = 21, cex = 1, bg = adjustcolor(cols, 0.5), col = adjustcolor('black', 0.3),xaxt = 'n', yaxt = 'n',
       xlim = c(-2, 5), ylim = c(2, 10), cex.main = 2, font.main = 1, main = nms[i])
}
mtext("clonal phenotype densities", side = 3, line = 1, outer = TRUE, cex = 2, font = 2)
dev.off()


# estimate enrichment of clone phenotype densitis along PCs
cor.permut <- lapply(1:10000, function(i){
  #print(i)
  cls1 <- sample(cls)
  cor(clone.pheno[cls1, ii], hsc.env$prr$x[cls, 1:3])
})
cor_array <- simplify2array(cor.permut)
cor_means <- apply(cor_array, c(1, 2), mean)
cor_sd <- apply(cor_array, c(1, 2), sd)
x <- cor(clone.pheno[cls, ii], hsc.env$prr$x[cls, 1:3])
y <- (x-cor_means)/cor_sd
y <- t(y)

png('figures/hsc/hsc_phenotype_density_pca.png', width = 3 * 150, height = 2 * 150)
par(mfrow=c(2,3), mar = c(2, 2, 0, 2))
col.types <- c('RdBu', 'PiYG', 'PRGn')
nms.types <- c('growth', 'output', 'unknown')
reord.pnames <- c("lymphoid bias", "clone size", "growth rate", "erythroid bias", "myeloid bias", "HSC bias", "dendritic bias")
reord.pnames <- rev(reord.pnames)
for (i in 1:3){
  #cols <-  colorRampPalette(c("blue","grey90", "red" ))(n = 10)[as.numeric(cut(rank(prr$x[,i]), 10))]
  cols <-  brewer.pal(10, col.types[i])[as.numeric(cut(rank(-hsc.env$prr$x[cls,i]), 10))]
  plot(rna.env$emb[cls,], pch = 19, cex = 0.5, col = adjustcolor(cols, 0.5), bty = 'n', xaxt = 'n', yaxt = 'n',
       #main = nms.types[i],
       font.main=1, cex.main = 2, xlim = c(-2, 5), ylim = c(3, 9))
}
xmax <- 8#max(abs(y))
for (i in 1:3){
  #cols <-  brewer.pal(10, col.types[i])[as.numeric(cut(-y[i,], 10))]
  cols <-  brewer.pal(10, col.types[i])[as.numeric(cut(-y[i,],
                                                       breaks = seq(-max(abs(y[i,])+1e-1), max(abs(y[i,])+1e-1), length.out = 10)))]
  barplot(y[i,reord.pnames], horiz = T, col =  cols[match(reord.pnames, colnames(y))], cex.axis = 2,
          xlim = c(-xmax,xmax),
          names.arg = '')
}
dev.off()


