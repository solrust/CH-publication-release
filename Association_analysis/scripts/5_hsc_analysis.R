load('../RNA_analysis/data/processed/scrna.RData')
load('../RNA_analysis/data/processed/mutant_de.RData')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
load('../TF_overexpression/data/processed/TF_trends.RData')
source('R/association_routines.R')
library(Matrix)
library(RColorBrewer)
library(parallel)
library(Matrix)
library(MASS)
library(glmmTMB)
library(ggplot2)
library(data.table)
library(gridExtra)
hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
tfs.list <- readRDS('../RNA_analysis/data/genesets/processed/tfs.rds')
tfs.all <- unlist(tfs.list)
save_folder <- 'default'
data_folder <- paste('data/processed/', save_folder, '/', sep = '')
load(paste(data_folder, 'association_summary.RData', sep = ''))

load(file = 'data/processed/hsc/hsc_axes.RData')

######### run analysis
do.small.set <- F
tfs.lfc <- do.call(cbind, lapply(test.pc.mutant, function(t) t$test.pc[,1] ))
tfs.p <- do.call(cbind, lapply(test.pc.mutant, function(t) t$test.pc[,5] ))
tfs.lfc <- tfs.lfc[intersect(tfs.all, rownames(tfs.lfc)),]
tfs.p <- tfs.p[intersect(tfs.all, rownames(tfs.p)),]
tshow <- apply(tfs.p, 1, min) < 1e-5 &
  rank(-apply(abs(tfs.lfc), 1, max))  <= ifelse(do.small.set, 40, 150)
sum(tshow)
x <- tfs.lfc[tshow,]
rownames(x) <- unlist(lapply(rownames(x), function(y) paste(substr(y,1,1),substr(tolower(y), 2, nchar(y)), sep = '') ))
pheat <- pheatmap::pheatmap(x, cluster_cols = F, breaks = seq(-0.2, 0.2, length.out = 100), cellheight = 22, angle_col = 0, fontsize = 20,
                            labels_col = c('PC1', 'PC2', 'PC3'))
#
pdf(paste('figures/hsc/hsc_TFs', ifelse(!do.small.set, '_large', '') ,'.pdf', sep = ''),
    height = ifelse(do.small.set, 10, 10 * sum(tshow) / 31), width = 5)
pheat
dev.off()

### plot genotype cell enrichment across PC axes
x <- data.frame(sample = sampleid[rownames(prr$x)], hash.match[sampleid[rownames(prr$x)],])
y <- do.call(rbind, lapply(1:3, function(i){
  ff <- 0.5
  ind1 <- rank(prr$x[,i]) <= ff * nrow(prr$x)
  n1 <- table(x$sample[ind1])
  ind2 <- rank(-prr$x[,i]) <= ff * nrow(prr$x)
  n2 <- table(x$sample[ind2])
  #data.frame(n1 = n1, n2 = n2, sample = names(n1), pc = rep(paste('PC',i,sep=''), length(n1)))
  as.data.frame(cbind(n1 = n1, n2 = n2, sample = names(n1), pc = rep(paste('PC',i,sep=''), length(n1)) ))
}))
y[,1] <- as.numeric(y[,1]); y[,2] <- as.numeric(y[,2])
y <- cbind(y, hash.match[y$sample, c('genotype', 'exp')])
y$facet <- paste(y$exp, y$pc, sep = '_')
#ord <- order(y$genotype, y$n2 / y$n1)
#head(y[ord,])
z <- unique(y$facet); names(z) <- z; z <- substr(z, 5, nchar(z))
pl <- ggplot(y, aes(rownames(y), (n2 / n1), fill = genotype)) +
  theme_bw() + geom_col(col = 'black', lwd = 0.1) + labs(x = 'sample', y = 'cell ratio, 50% top to bottom PC ')  +
  scale_fill_manual(values = col.genotype) +
  theme(text = element_text(size = 20), legend.title = element_blank(), axis.text.x = element_blank()) +
  #facet_wrap(~ exp + pc, scales = 'free_x', ) + theme(strip.text.x  = element_blank())
  facet_wrap(~ facet, scales = 'free_x', labeller = labeller(facet = z)) #+ theme(strip.text.x  = element_blank())
pl
#
pdf('figures/hsc/hsc_density_bias.pdf', width = 5, height = 5)
pl
dev.off()


## check stat sign:
y$ratio <- y$n2 / y$n1
y$total <- y$n1 + y$n2
do.call(rbind, lapply(unique(y$pc), function(pc){
  unlist(lapply(unique(y$exp), function(ex){
    y1 <- y[y$pc == pc & y$exp == ex,]

    ft <- MASS::glm.nb(y1$n1 ~ y1$genotype + offset(log(y1$total)))
    sm <- summary(ft)
    sm$coefficients[2,4]

    ft <-  glmmTMB(y1$n1 ~ y1$genotype + offset(log(y1$total)),
                   family = nbinom1(link = 'log'))
    sm <- summary(ft)
    sm$coefficients$cond[2,4]
  }))
}))


### visulalize TF expression and TF perturbation
head(tf.trends.no25.v1m[,c('t', 'p.t')])
bias_site1 <- data.frame(bias_site,  TF = barcode2tf$Gene[match(rownames(bias_site), barcode2tf$Sequence)])
#
gg <- unique(c(rownames(tf.trends.no25.v1m), rownames(bias_site1)))
str(gg)
#tf.features <- data.frame(cbind(tf.trends.no25.v1m[gg, c('t', 'p.t')], bias_save[gg, c('My-B', 'p.My-B')]))
tf.features <- cbind(tf.trends.no25.v1m[gg, c('t', 'p.t')], bias_site1[gg, 1:2]); colnames(tf.features)[3:4] <- c('My-B', 'p.My-B')
tf.features <- data.frame(tf.features)
tf.features$My.B[abs(tf.features$My.B) > 10 ] <- NA
rownames(tf.features) <- gg
tf.features$TF <- barcode2tf$Gene[match(rownames(tf.features), barcode2tf$Sequence)]
head(tf.features)
tf.features1 <- data.frame(rbind(as.matrix(cbind(tf.features[,c('t', 'p.t', 'TF')], 't')), as.matrix(cbind(tf.features[,c('My.B', 'p.My.B', 'TF')], 'bias') )))
head(tf.features1)
colnames(tf.features1) <- c('pheno', 'p', 'TF', 'name')
tf.features1[,1] <- as.numeric(tf.features1[,1]); tf.features1[,2] <- as.numeric(tf.features1[,2])
tf.features1$TF <- factor(tf.features1$TF, level = unique(tf.features1$TF))
str(tf.features1$TF )
pcut <- 5e-2
tf.features1$trend <- 'no'
tf.features1$trend[tf.features1$p < pcut & tf.features1$pheno > 0 & tf.features1$name == 't'] <- 'grow'
tf.features1$trend[tf.features1$p < pcut & tf.features1$pheno < 0 & tf.features1$name == 't'] <- 'decrease'
tf.features1$trend[tf.features1$p < pcut & tf.features1$pheno > 0 & tf.features1$name == 'bias'] <- 'B'
tf.features1$trend[tf.features1$p < pcut & tf.features1$pheno < 0 & tf.features1$name == 'bias'] <- 'My'
table(tf.features1$trend)
ind <- tf.features1$TF %in% unique(tf.features1$TF[which(tf.features1$p < 5e-2)])
tf.features1$trend <- factor(tf.features1$trend, levels = c('B', 'My', 'grow', 'no', 'decrease'))
pl1 <- ggplot(tf.features1[ind,])  + theme_bw() +
  geom_bar(aes(pheno, TF, fill = trend), col = 'black', lwd = 0.2, position = position_dodge(),  stat = 'identity',
                                                           orientation = 'y') +
  scale_fill_manual(values = c('decrease' = 'blue', 'grow' = 'red', 'My' = 'orange', 'B' = 'darkgreen', 'no' = 'grey')) +
  #scale_fill_manual(values = c('down' = 'Blue Label', 'up' = 'Red Label', 'My' = 'Orange Label', 'B' = 'Green Label', 'no' = 'Grey Label'), name = "Your Legend Title") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x='', y = '') +
  facet_wrap(~name, scales = 'free_x') + theme(strip.text = element_blank())
pl1
head(tf.features1[ind,])
###
tshow <- as.character(unique(tf.features1$TF[ind]))
str(tshow)
#
tfs.lfc <- do.call(cbind, lapply(test.pc.mutant, function(t) t$test.pc[,1] ))
tfs.p <- do.call(cbind, lapply(test.pc.mutant, function(t) t$test.pc[,5] ))
tadd <- setdiff(tshow, rownames(tfs.lfc))
str(tadd)
na_row <- as.data.frame(matrix(NA, nrow = length(tadd), ncol = 3))
rownames(na_row) <- tadd; colnames(na_row) <- colnames(tfs.lfc)
tfs.lfc <- rbind(tfs.lfc, na_row)
tfs.p <- rbind(tfs.p, na_row)
head(tfs.lfc)
tail(tfs.lfc)
tfs.lfc$TF <- tfs.p$TF <- factor(rownames(tfs.lfc), level = tshow)
#
x1 <- melt(tfs.lfc[tshow,], na.rm = F, id.vars = "TF"); x2 <- melt(tfs.p[tshow,], na.rm = F, id.vars = "TF")
x1 <- cbind(x1, p = x2$value)
#x1$Var1 <- factor(x1$Var1, level = tshow)
x1$state <-  ifelse(x1$value > 0, 'up', 'down'); x1$state[x1$p > 1e-2] <- ''
pl2 <- ggplot(x1) + theme_bw() + geom_point(aes(variable, TF, size = abs(value), col = state ))  +
  theme(axis.text.y = element_blank(), #axis.text.x = element_text(size = 15),
        legend.position = "bottom", legend.title = element_blank()) +
  labs(y = '', x = '') + guides(size = FALSE) +
  scale_color_manual(values = c('up' = adjustcolor('red', 0.5), 'down' = adjustcolor('blue', 0.5), '' = 'grey'))
pl2
#
grid.arrange(grobs = list(pl1, pl2), ncol = 2, widths = c(3, 1))

#
pdf('figures/hsc/TF_functional_expression.pdf', width = 6, height = 8)
grid.arrange(grobs = list(pl1, pl2), ncol = 2, widths = c(3, 1))
dev.off()


