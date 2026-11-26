source('R/regression_time.R')
library(ggplot2)
library(parallel)
library(ggrepel)
library(data.table)
col.genotype <- adjustcolor(c('black', 'red', 'green'), 0.5)
names(col.genotype) <- c('WT', 'DNMT3A', 'TET2')


barcode2tf <- read.csv('data/raw/barcode list Zhang copy.txt', sep = '')
# check that Zhang's barodes are well separated between each other
seqdist <- adist(barcode2tf$Sequence)
diag(seqdist) <- NA
y <- apply(seqdist, 1, function(a) min(a, na.rm = T))
qplot(y, xlim = c(0, max(y)), xlab = 'min dist to other Zhang\'s barcode', main = paste('Min dist is', min(y)))
gfp.barcode <- barcode2tf$Sequence[barcode2tf$Gene == 'GFP']

sample.list1 <- read.csv('data/raw/sample list.csv', sep = ',')
sample.list1$Treatment <- sample.list1$X
sample.list1 <- sample.list1[,1:4]
##
conv.s2g <- sample.list1[sample.list1$Timepoint == '6 weeks', c(2,4)]
rownames(conv.s2g) <- conv.s2g$Mouse
conv.s2g <- conv.s2g[, 2, drop = F]
##
sample.list1 <- sample.list1[sample.list1$Treatment != '' | sample.list1$Sample.name != '',] # remove empty rows
sample.list1$Treatment <- conv.s2g[sample.list1$Mouse,]
sample.list1$Timepoint[grepl('18w', sample.list1$Sample.name)] <- '18 weeks'
sample.list1$Timepoint[grepl('PB 25 weeks', sample.list1$Timepoint)] <- '25 weeks'
sample.list1$site <- sample.list1$Timepoint
sample.list1$site[sample.list1$site %in% c('6 weeks', '12 weeks', '18 weeks', '25 weeks')] <- 'PB'
rownames(sample.list1) <- sample.list1$Sample.name
sample.list1$tm <- as.numeric(substr(sample.list1$Timepoint, 1, 2))
sample.list <- sample.list1


counts.tf <- do.call(cbind, mclapply(sample.list1$Sample.name, function(nm){
  print(nm)
  nm <- paste(unlist(strsplit(nm, split = '\\.'))[1], '_s.txt', sep = '')
  fl <- paste("data/raw/in vivo/", nm, sep = '')
  if (!file.exists(fl)){
    fl <- paste("data/raw/", nm, sep = '')
  }
  x <- read.csv(fl, sep = '\t', header = F)
  x <- x[,1:2]
  colnames(x) <- c('barcode', 'counts')
  ### match starcode barcodes to Zhang list
  seqdist <- adist(barcode2tf$Sequence, x$barcode)
  abund <- (seqdist <= 1) %*% x$counts
  #abund[abund == 0] <- NA
  return(abund)
}, mc.cores = 1))
rownames(counts.tf) <- barcode2tf$Sequence
colnames(counts.tf) <- sample.list1$Sample.name



counts.tf.vitro <- do.call(cbind, mclapply(unique(sample.list1$Treatment), function(nm){
  print(nm)
  nm <- paste(unlist(strsplit(nm, split = '\\.'))[1], '_s.txt', sep = '')
  fl <- paste("data/raw/in vitro/", nm, sep = '')
  x <- read.csv(fl, sep = '\t', header = F)
  x <- x[,1:2]
  colnames(x) <- c('barcode', 'counts')
  ### match starcode barcodes to Zhang list
  seqdist <- adist(barcode2tf$Sequence, x$barcode)
  abund <- (seqdist <= 1) %*% x$counts
  #abund[abund == 0] <- NA
  return(abund)
}, mc.cores = 10))
rownames(counts.tf.vitro) <- barcode2tf$Sequence
colnames(counts.tf.vitro) <- unique(sample.list1$Treatment)

# merge in vitro and in vivo data
x <- as.data.frame(do.call(rbind, lapply(1:4, function(i){
  nm <- paste(unlist(strsplit(unique(sample.list$Treatment)[i], split = '\\.'))[1], '_s.txt', sep = '')
  c(nm, paste(unique(sample.list$Treatment)[i], 'vitro', sep = '_'), '0 weeks', unique(sample.list$Treatment)[i], 'in vitro', 0)
})))
colnames(x) <- colnames(sample.list)
sample.list <- rbind(sample.list, x)
sample.list$t <- as.numeric(sample.list$t)
counts.tf <- cbind(counts.tf, counts.tf.vitro)
colnames(counts.tf) <- rownames(sample.list)
expr <- t(t(counts.tf) / colSums(counts.tf))

### check that HSC and other samples express HSC_list and MUT_list:
hsc.list <- read.csv('data/raw/HSC_list.txt', sep = '\t', header = F)
hsc.list1 <- read.csv('data/raw/HSC plasmid.txt', sep = '\t', header = F)
x <- rowSums(counts.tf[, sample.list$Treatment != 'HSC'])
flag <- barcode2tf$Gene %in% hsc.list$V2
qplot(x, log = 'x', fill = flag, alpha = I(0.2), geom = 'density')

#
mut.list <- read.csv('data/raw/MUT_list.txt', sep = '\t', header = F)
mut.list1 <- read.csv('data/raw/mut plasmid.txt', sep = '\t', header = F)
x <- rowSums(counts.tf[, sample.list$Treatment != 'HSC'])
flag <- barcode2tf$Gene %in% mut.list$V2
qplot(x, log = 'x', fill = flag, alpha = I(0.2), geom = 'density')

x1 <- rowMeans(counts.tf[, sample.list$Treatment == 'HSC'])
x2 <- rowMeans(counts.tf[, sample.list$Treatment != 'HSC'])
flag <- rep('NA', nrow(barcode2tf))
flag[barcode2tf$Gene %in% mut.list1$V2] <- 'mut'
flag[barcode2tf$Gene %in% hsc.list1$V2] <- 'hsc'
flag[barcode2tf$Gene %in% mut.list1$V2 & barcode2tf$Gene %in% hsc.list1$V2] <- 'both'
lbl <- barcode2tf$Gene; lbl[x1 < 1e+2 & x2 < 1e+2] <- NA
pl <- ggplot(, aes(x1, x2, col = flag, label = lbl)) + theme_bw() + geom_point() + scale_x_log10() + scale_y_log10() +
  geom_text_repel(size = I(2), max.overlaps = 30) + labs(x = 'mean counts, HSC samples', y = 'mean counts, non-HSC samples')  + geom_abline()
pl
png('figures/counts_HSC_non-HSC.png', width = 700, height = 700)
pl
dev.off()
### explore count trends
ind <- sample.list$Timepoint %in% c('0 weeks', '6 weeks', '12 weeks', '18 weeks', '25 weeks')
barabund <- rowSums(counts.tf[,ind])
barabund <- apply(counts.tf[,ind], 1, function(x) sum(x > 1e+1) )
qplot(barabund, log = 'x')
ind1 <- barabund > 6#1e+2
sum(ind1)


###### run TF time-series analysis
include.samples <- sample.list$Timepoint %in% c('0 weeks', '6 weeks', '12 weeks', '18 weeks', '25 weeks')
table(include.samples)
cov.high <- rowSums(counts.tf[, include.samples] > 1e+1)
include.barcodes <- cov.high > 12
table(include.barcodes)

### run mixed time-series model without HSC population
include.samples.nohcs <- include.samples & sample.list$Treatment != 'HSC'
table(include.samples.nohcs)
include.samples.hcs <- include.samples & sample.list$Treatment == 'HSC'
table(include.samples.hcs)
####
tf.trends.v1m <- run_reg(counts.tf, treatment = NULL, perm = FALSE, meth = 'nbinom1', mixed.treatment = TRUE, run_mixed_model = TRUE,
                         include.barcodes = include.barcodes, include.samples = include.samples.nohcs, time.points = NA,
                         nc.cores = 10)
tf.trends.hsc <- run_reg(counts.tf, treatment = NULL, perm = FALSE, meth = 'nbinom1',
                         mixed.treatment = FALSE, run_mixed_model = TRUE,
                         include.barcodes = include.barcodes, include.samples = !include.samples.hcs, time.points = NA,
                         nc.cores = 10)
rownames(tf.trends.v1m) <- rownames(tf.trends.hsc) <- rownames(counts.tf)[include.barcodes]

### run mixed time-series model without HSC population and last time point
include.samples.nohcs.no25 <- include.samples & sample.list$Treatment != 'HSC' & sample.list$t != 25
table(include.samples.nohcs.no25)
cov.high.no25 <- rowSums(counts.tf[, include.samples.nohcs.no25] > 1e+1)
include.barcodes.no25 <- cov.high.no25 > 10
table(include.barcodes.no25)
tf.trends.no25.v1m <- run_reg(counts.tf, treatment = NULL, perm = FALSE, meth = 'nbinom1', mixed.treatment = TRUE, run_mixed_model = TRUE,
                              include.barcodes = include.barcodes, include.samples = include.samples.nohcs.no25, time.points = NA,
                              nc.cores = 10)
rownames(tf.trends.no25.v1m) <- rownames(counts.tf)[include.barcodes]

colnames(tf.trends.hsc) <- c('intercept', 't', 'tDNMT3A', 'tTET2', 'p.intercept', 'p.t', 'p.tDNMT3A', 'p.tTET2')
colnames(tf.trends.v1m) <- c('intercept', 't', 'tDNMT3A', 'tTET2', 'p.intercept', 'p.t', 'p.tDNMT3A', 'p.tTET2')
colnames(tf.trends.no25.v1m) <- c('intercept', 't', 'tDNMT3A', 'tTET2', 'p.intercept', 'p.t', 'p.tDNMT3A', 'p.tTET2')
rownames(tf.trends.v1m) <- rownames(tf.trends.no25.v1m) <- rownames(counts.tf)[include.barcodes]
tf.trends.v1m <- data.frame(tf.trends.v1m)
tf.trends.v1m$TF <- barcode2tf$Gene[match(rownames(counts.tf)[include.barcodes], barcode2tf$Sequence)]
tf.trends.hsc <- data.frame(tf.trends.hsc)
tf.trends.hsc$TF <- barcode2tf$Gene[match(rownames(counts.tf)[include.barcodes], barcode2tf$Sequence)]
tf.trends.no25.v1m <- data.frame(tf.trends.no25.v1m)
tf.trends.no25.v1m$TF <- barcode2tf$Gene[match(rownames(counts.tf)[include.barcodes], barcode2tf$Sequence)]

#saveRDS(tf.trends.v1m, '/home/solrust/CH/draft/data/pooledTF_trends.rds')
#saveRDS(tf.trends.no25.v1m, '/home/solrust/CH/draft/data/pooledTF_trends_wo25.rds')

tfnames <-  barcode2tf$Gene[match(rownames(counts.tf)[include.barcodes], barcode2tf$Sequence)]
tfsign <- do.call(rbind, lapply(0:2, function(i){
  trnd <- sign(tf.trends.v1m[,2+i]) + sign(tf.trends.no25.v1m[,2+i])
  trnd[trnd > 0] <- 'up'; trnd[trnd < 0] <- 'down'; trnd[trnd == 0] <- 'mixed'
  tfsign <- data.frame(pv1 = tf.trends.v1m[,6+i], pv2 = tf.trends.no25.v1m[,6+i], TF = tfnames,
                       trend = trnd, comparison = as.factor(ifelse(i==0, 'WT', ifelse(i==1, 'DNMT3A-WT', 'TET2-WT'))))
}))
cutx <- 1e-2; cuty <- 1e-2
pl <- ggplot(tfsign, aes(pv1, pv2, label = TF, fill = trend)) + geom_point(col = 'black', alpha = I(0.8), size = 2, shape = 21) + theme_bw() +
  #theme(legend.position = c(0.1, 0.9), legend.text = element_text(size = 20),  legend.title = element_text(size = 0)) +
  theme(text = element_text(size = 20))  +
  theme(legend.position = 'top',
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 0)) +
  geom_text_repel(max.overlaps = 10) +
  #scale_x_log10(limits = c(min(c(tfsign$pv1, tfsign$pv2), na.rm = T) / 1, 1)) + scale_y_log10(limits = c(min(c(tfsign$pv1, tfsign$pv2), na.rm = T) / 1, 1)) +
  scale_x_log10() + scale_y_log10() +
  scale_fill_manual(values = c(up = 'red', down = 'blue', mixed = 'grey')) +
  labs(x = '0-25 weeks, p-value', y = '0-18 weeks, p-value') +
  geom_vline(xintercept = cutx, linetype = 2, size = 0.5) + geom_hline(yintercept = cuty, linetype = 2, size = 0.5) +
  geom_abline(linetype = 2, size = 0.1) + facet_wrap(~comparison,scales='free')
pl


pdf('figures/pvalues_mixed_18vs25.pdf', width = 11, height = 4)
pl
dev.off()


###
for (incl.25 in c(F, T)){
#incl.25 <- F
tf.tr <- tf.trends.no25.v1m
if (incl.25) tf.tr <- tf.trends.v1m
tf.all <- do.call(rbind, lapply(0:2, function(i){
  trnd <- sign(tf.tr[,2+i]) * as.numeric(tf.tr[,6+i] < 1e-2)
  trnd[trnd > 0] <- 'up'; trnd[trnd < 0] <- 'down'; trnd[trnd == 0] <- 'no'
  tfsign <- data.frame(logFC = tf.tr[,2+i], pv = tf.tr[,6+i], TF = tfnames,
                       trend = trnd, comparison = as.factor(ifelse(i==0, 'WT', ifelse(i==1, 'DNMT3A-WT', 'TET2-WT'))))
}))
pl <- ggplot(tf.all, aes(logFC, pv, label = TF, fill = trend)) + geom_point(col = 'black', alpha = I(0.8), size = 2, shape = 21) + theme_bw() +
  theme(text = element_text(size = 20))  +
  theme(legend.position = 'top', legend.text = element_text(size = 20),  legend.title = element_text(size = 0)) +
  geom_text_repel(max.overlaps = 15) +
  scale_y_log10() + xlim(-0.5, 0.5) +
  scale_fill_manual(values = c(up = 'red', down = 'blue', no = 'grey')) +
  labs(x = paste('growth rate 0-',ifelse(incl.25, '25', '18'), ', weeks', sep = ''), y = 'p-value') +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +#+ geom_hline(yintercept = cuty, linetype = 2, size = 0.5) +
  facet_wrap(~comparison, scales='free')
pl

pdf(paste('figures/volcano_mixed_', ifelse(incl.25, '25', '18'),'.pdf', sep = ''),
    width = 11, height = 4)
print(pl)
dev.off()
}

##### plot TF-specific trends
size.tf <- colSums(counts.tf[include.barcodes,])
#size.tf <- counts.tf[barcode2tf$Sequence[barcode2tf$Gene == 'GFP'], ]
pdf('figures/time_trends_mixed.pdf')
tf.trnd <- tf.trends.v1m
tm <- seq(0, 25, length.out = 100)
#tf.trnd <- tf.trends.no25.v1
for (brcd in rownames(tf.trends.v1m)[order(tf.trends.v1m[,6])] ){
  #brcd <- rownames(tf.trnd)[order(tf.trnd[,6])][1:10][1]
  #brcd <- barcode2tf$Sequence[barcode2tf$Gene == 'GFP'][1]
  print(brcd)
  tfsample <- cbind(sample.list, counts = counts.tf[brcd, ] / size.tf)
  tfname <- barcode2tf$Gene[barcode2tf$Sequence == brcd]
  tfname

  datfit <- rbind(data.frame(tm, trend = exp(tf.trnd[brcd,1] - 0 + tf.trnd[brcd,2] * tm), Treatment = 'AAVS1'),
                  data.frame(tm, trend = exp(tf.trnd[brcd,1] - 0 + (tf.trnd[brcd,2] + tf.trnd[brcd,3]) * tm), Treatment = 'DNMT3A'),
                  data.frame(tm, trend = exp(tf.trnd[brcd,1] - 0 + (tf.trnd[brcd,2] + tf.trnd[brcd,4]) * tm), Treatment = 'TET2'))

  pp <- ggplot(data = subset(tfsample, include.samples),
               aes(t, counts, col = Treatment, linetype = Mouse)) +
    scale_x_continuous(breaks = sort(unique(subset(tfsample, include.samples)$t))) +
    scale_y_log10() +
    xlab('time') + ylab('frequency') +
    geom_line(data = datfit, aes(tm, trend, color = Treatment),
              alpha = I(0.5), size = 2, linetype = 1) +
    geom_line(alpha = I(0.9)) + geom_point(alpha = I(1), size = 1) +
    scale_color_manual(values = c('AAVS1' = as.character(col.genotype['WT']), 'DNMT3A' = as.character(col.genotype['DNMT3A']), 'TET2' = as.character(col.genotype['TET2']), 'HSC' = 'grey' )) +
    scale_linetype_manual(values = rep(1, length(unique(subset(tfsample, include.samples)$Mouse))), guide = 'none') +
    annotate(geom = "text", x = 3, y = max(tfsample$counts), label= tfname, size = I(5)) +
    theme_bw() + theme(text = element_text(size = 20))
  print(pp)
}
dev.off()

####### examine lineage biases
smpls <- sample.list$site %in% c('Spleen B cell', 'SP myeloid') & sample.list$Treatment %in% c('AAVS1', 'DNMT3A', 'TET2')
sum(smpls)
gg <- rowSums(counts.tf[,smpls] > 0) > 4
sum(gg)
# build permut matrix
indx <- rep(NA, sum(smpls))
for (msid in unique(sample.list[smpls,]$Mouse)){
  xx <- which(sample.list[smpls,]$Mouse == msid)
  indx[xx] <- sample(xx)
}
counts.tf.mut <- counts.tf[,smpls][,indx]
colnames(counts.tf.mut) <- colnames(counts.tf)[smpls]
#
do.perm <- FALSE
bias_site <- do.call(rbind, mclapply(rownames(counts.tf[gg,]), function(gn){
  y <- counts.tf[gn,smpls]
  logsize1 <- log(colSums(counts.tf[,smpls]))
  if (do.perm){
    y <- counts.tf.mut[gn,]
    logsize1 <- log(colSums(counts.tf.mut[,]))
  }
  sitt <- sample.list[smpls, c('site')]
  trtmt <- sample.list[smpls, c('Treatment')]
  smpl <- sample.list[smpls, c('Mouse')]
  ft <- glmmTMB(y ~ offset(logsize1) + sitt + trtmt,
                family=nbinom1(link = 'log'))
  c(summary(ft)$coefficients$cond[2,c(1,4)],
    summary(ft)$coefficients$cond[3,c(1,4)],
    summary(ft)$coefficients$cond[4,c(1,4)])
}, mc.cores = 10))
rownames(bias_site) <- rownames(counts.tf[gg,])

bias_site_geno <- do.call(rbind, mclapply(rownames(counts.tf[gg,]), function(gn){
  y <- counts.tf[gn,smpls]
  logsize1 <- log(colSums(counts.tf[,smpls]))
  if (do.perm){
    y <- counts.tf.mut[gn,]
    logsize1 <- log(colSums(counts.tf.mut[,]))
  }
  sitt <- sample.list[smpls, c('site')]
  trtmt <- sample.list[smpls, c('Treatment')]
  smpl <- sample.list[smpls, c('Mouse')]
  ft <- glmmTMB(y ~ offset(logsize1) + sitt + trtmt + sitt:trtmt,
                family=nbinom1(link = 'log'))
  c(summary(ft)$coefficients$cond[2,c(1,4)],
    summary(ft)$coefficients$cond[3,c(1,4)],
    summary(ft)$coefficients$cond[4,c(1,4)],
    summary(ft)$coefficients$cond[5,c(1,4)],
    summary(ft)$coefficients$cond[6,c(1,4)])
}, mc.cores = 10))
rownames(bias_site_geno) <- rownames(counts.tf[gg,])

tfnms <- barcode2tf$Gene[match(rownames(bias_site), barcode2tf$Sequence)]
trnd <- sign(bias_site[,1]) * as.numeric(p.adjust(bias_site[,2], 'BH') < 2e-1)
trnd[trnd > 0] <- 'B'; trnd[trnd < 0] <- 'myeloid'; trnd[trnd == 0] <- 'no'
bias1 <- data.frame('logFC' = bias_site[,1], 'pv' = bias_site[,2],
                    TF = tfnms, trend  = trnd, comparison = as.factor('WT'))
#
trnd <- sign(bias_site_geno[,7]) * as.numeric(p.adjust(bias_site_geno[,8], 'BH') < 2e-1)
trnd[trnd > 0] <- 'B'; trnd[trnd < 0] <- 'myeloid'; trnd[trnd == 0] <- 'no'
bias2 <- data.frame('logFC' = bias_site_geno[,7], 'pv' = bias_site_geno[,8],
                    TF = tfnms, trend  = trnd, comparison = as.factor('DNMT3A-WT'))
#
trnd <- sign(bias_site_geno[,9]) * as.numeric(p.adjust(bias_site_geno[,10], 'BH') < 2e-1)
trnd[trnd > 0] <- 'B'; trnd[trnd < 0] <- 'myeloid'; trnd[trnd == 0] <- 'no'
bias3 <- data.frame('logFC' = bias_site_geno[,9], 'pv' = bias_site_geno[,10],
                    TF = tfnms, trend  = trnd, comparison = as.factor('TET2-WT'))
biastot <- rbind(bias1, bias2, bias3)
biastot$trend[is.na(biastot$trend)] <- 'no'


pl <- ggplot(biastot, aes(logFC, pv, label = TF, fill = trend)) + geom_point(col = 'black', alpha = I(0.8), size = 2, shape = 21) + theme_bw() +
  theme(text = element_text(size = 20))  +
  theme(legend.position = 'top', legend.text = element_text(size = 20),  legend.title = element_text(size = 0)) +
  geom_text_repel(max.overlaps = 20) +
  scale_y_log10() + xlim(-3, 3) +
  scale_fill_manual(values = c(B = 'darkgreen', myeloid = 'darkorange', no = 'grey')) +
  labs(x = 'B/My bias, log', y = 'p-value') +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +#+ geom_hline(yintercept = cuty, linetype = 2, size = 0.5) +
  facet_wrap(~comparison, scales='free')
pl

pdf('figures/volcano_bias.pdf', width = 11, height = 4)
pl
dev.off()

### print per-gene bias
ind1 <- sample.list[smpls,]$site == 'Spleen B cell'
ind2 <- sample.list[smpls,]$site == 'SP myeloid'
fctr <- sample.list$Treatment[smpls][sample.list[smpls,]$site == 'Spleen B cell']
sz1 <- colSums(counts.tf[gg, smpls][,ind1])
x1 <- as.numeric(t(t(counts.tf[gg, smpls][,ind1]) / sz1))
sz2 <- colSums(counts.tf[gg, smpls][,ind2])
x2 <- as.numeric(t(t(counts.tf[gg, smpls][,ind2])/sz2))

cnts <- t(t(counts.tf[, smpls]) / colSums(counts.tf)[smpls])
tbx <- melt(cnts)
tbx <- cbind(tbx, site = sample.list[as.character(tbx$Var2),c('site')],
             mouse = sample.list[as.character(tbx$Var2),c('Mouse')],
             genotype = sample.list[as.character(tbx$Var2),c('Treatment')])
freq.bias <- dcast(tbx[,-2], Var1+mouse+genotype ~ ..., sum)
freq.bias$TF <- barcode2tf$Gene[match(barcode2tf$Sequence, freq.bias$barcode)]
colnames(freq.bias) <- c('barcode', 'mouse', 'genotype', 'myeloid', 'B')

pdf('figures/bias_genes.pdf', width = 6, height = 4)
for (barid in rownames(bias_site)[order(bias_site[,2])] ){
  alpp <- ifelse(as.character(freq.bias$barcode) == barid, 1, 0.1)
  cl <- ifelse(as.character(freq.bias$barcode) == barid, 'black', 'white')
  tf.name <- 'GFP'
  tf.name <- barcode2tf$Gene[barcode2tf$Sequence == barid]
  pl <- ggplot(freq.bias, aes(myeloid, B, fill = genotype)) +
    geom_point(col = cl, shape = 21, size = 3, alpha = I(alpp)) +
    geom_point(data = subset(freq.bias, barcode == barid), aes(myeloid, B, fill = genotype),
               shape = 21, size = 3) +
    theme_bw() + geom_abline(linetype = 2) + scale_x_log10() + scale_y_log10() +
    annotate("text", x = 1e-4, y = 0.5, label= tf.name, size = 6) +
    labs(x = 'frequency, Myeloid', y = 'frequency, B cells') +  theme(text = element_text(size = 20))
  print(pl)
}
dev.off()

save(barcode2tf, sample.list1, tf.trends.no25.v1m, tf.trends.v1m, bias_site, bias_site_geno, biastot,
           file = 'data/processed/TF_trends.RData')


