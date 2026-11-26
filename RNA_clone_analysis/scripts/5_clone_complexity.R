library(ggplot2)
library(ggrepel)
library(gridExtra)
library(parallel)
library(data.table)
load('../RNA_analysis/data/processed/scrna.RData')
load('data/processed/clone_info.RData')
hash.annot <- readRDS(file = "../RNA_analysis/data/processed/hashtag_info.rds")
hash.match <- readRDS(file = "../RNA_analysis/data/processed/hashtag_meta.rds")
expids <- c('exp123', 'exp127')
smpls <- c('RNA1', 'RNA2', 'RNA3', 'RNA4')

library(PLNmodels)
logNorm <- function(mlog, sdlog, len = 1000){
  sq <- exp(seq(log(1e-4), log(1), length.out = len))
  f <- exp(- (log(sq) - mlog)^2 / sdlog^2 / 2) / (sq * sdlog * sqrt(2 * pi))
  return(list(x = sq, f = f))
}

sz <- rowSums(clone.states)
#offset.mat <-as.matrix(summary.states.norm[clone.info[,]$id,]) * sz

#
geno <- 'DNMT3A'
clone.cov <- 20
ind0 <- sz >= clone.cov & clone.info$expid == ifelse(geno == 'TET2', '127', '123')

# prepare data for PLM modeling
cl <- 1:10
dt <- prepare_data(
  counts     = clone.states[ind0, cl, drop = F],
  covariates = clone.info[ind0, c('id', 'genotype'), drop = F],
  offset     = sz[ind0],
  #offset = offset.mat
)
dt$genotype <- factor(dt$genotype, levels = c('WT', geno))
dt$id <- as.character(dt$id)

#
fate.pln0 <- PLN(Abundance ~ 1 + id + offset(log(Offset)), data = dt, control = PLN_param(covariance = 'diagonal', trace = 0))
fate.pln.wt <- PLN(Abundance ~ 1 + id + offset(log(Offset)), subset = !dt$genotype == geno, data = dt, control = PLN_param(covariance = 'diagonal', trace = 0))
fate.pln.mut <- PLN(Abundance ~ 1 + id + offset(log(Offset)), subset = dt$genotype == geno, data = dt, control = PLN_param(covariance = 'diagonal', trace = 0))
#
dat0 <- data.frame(cbind(mean.wt = fate.pln.wt$model_par$B[1,] + colMeans(rbind(0, fate.pln.wt$model_par$B[-1,])), mean.mut = fate.pln.mut$model_par$B[1,] + colMeans(rbind(0, fate.pln.mut$model_par$B[-1,])), sd.wt = diag(fate.pln.wt$model_par$Sigma), sd.mut = diag(fate.pln.mut$model_par$Sigma), state = 1:nrow(fate.pln.wt$model_par$Sigma)))
dat <- dat0
dat$mean.wt <- exp(dat$mean.wt + 0.5 * dat$sd.wt)
dat$mean.mut <- exp(dat$mean.mut + 0.5 * dat$sd.mut)
dat$sd.wt <- exp(dat$sd.wt) - 1
dat$sd.mut <- exp(dat$sd.mut) - 1
# run bootsrap
bootsrap.wt0 <- do.call(rbind, mclapply(1:100000, function(i){
  ind <- which(!dt$genotype == geno)
  ind <- sample(ind, replace = T)
  fate.pln <- PLN(Abundance ~ 1 + id + offset(log(Offset)), subset = ind, data = dt, control = PLN_param(covariance = 'diagonal', trace = 0))
  c(fate.pln$model_par$B[1,] + colMeans(rbind(0, fate.pln$model_par$B[-1,])), diag(fate.pln$model_par$Sigma))
}, mc.cores = 20))
#
bootsrap.mut0 <- do.call(rbind, mclapply(1:100000, function(i){
  ind <- which(dt$genotype == geno)
  ind <- sample(ind, replace = T)
  fate.pln <- PLN(Abundance ~ 1 + id + offset(log(Offset)), subset = ind, data = dt, control = PLN_param(covariance = 'diagonal', trace = 0))
  c(fate.pln$model_par$B[1,] + colMeans(rbind(0, fate.pln$model_par$B[-1,])), diag(fate.pln$model_par$Sigma))
}, mc.cores = 20))
# estimate mean and CV of fate bias
bootsrap.wt <- cbind(exp(bootsrap.wt0[,1:10] + 0.5 * bootsrap.wt0[,11:20]), exp(bootsrap.wt0[,11:20]) - 1)
bootsrap.mut <- cbind(exp(bootsrap.mut0[,1:10] + 0.5 * bootsrap.mut0[,11:20]), exp(bootsrap.mut0[,11:20]) - 1)
#
qplot(dat$mean.wt, colMeans(bootsrap.wt[,1:10])) + geom_abline(linetype = 2)
qplot(dat$mean.mut, colMeans(bootsrap.mut[,1:10])) + geom_abline(linetype = 2)
qplot(dat$sd.wt, colMeans(bootsrap.wt[,11:20])) + geom_abline(linetype = 2)
qplot(dat$sd.mut, colMeans(bootsrap.mut[,11:20])) + geom_abline(linetype = 2)
#
dat$mean.wt.up <- apply(bootsrap.wt[,1:10], 2, function(x) quantile(x, probs = 0.95))
dat$mean.wt.down <- apply(bootsrap.wt[,1:10], 2, function(x) quantile(x, probs = 0.05))
dat$mean.mut.up <- apply(bootsrap.mut[,1:10], 2, function(x) quantile(x, probs = 0.95))
dat$mean.mut.down <- apply(bootsrap.mut[,1:10], 2, function(x) quantile(x, probs = 0.05))
dat$sd.wt.up <- apply(bootsrap.wt[,11:20], 2, function(x) quantile(x, probs = 0.95))
dat$sd.wt.down <- apply(bootsrap.wt[,11:20], 2, function(x) quantile(x, probs = 0.05))
dat$sd.mut.up <- apply(bootsrap.mut[,11:20], 2, function(x) quantile(x, probs = 0.95))
dat$sd.mut.down <- apply(bootsrap.mut[,11:20], 2, function(x) quantile(x, probs = 0.05))
pvs <- pmin(colMeans(bootsrap.wt - bootsrap.mut < 0), 1 - colMeans(bootsrap.wt - bootsrap.mut < 0)) * 2
dat$mean.pv <- pvs[1:10]
dat$sd.pv <- pvs[11:20]

save(bootsrap.wt0, bootsrap.wt, bootsrap.mut0, bootsrap.mut,
           dat, fate.pln.wt, fate.pln.mut,
           geno , clone.cov, ind0, file = paste('data/processed/clone_patterns_', geno, '.RData', sep = '') )


geno <- 'DNMT3A'
load(file = paste('data/processed/clone_patterns_', geno, '.RData', sep = ''))

# mean scatterplot
pl <- ggplot(data = dat) + geom_point(aes(mean.wt, mean.mut), size = I(2 - 2 * log10(dat$mean.pv + 1e-2)),
                                      fill = adjustcolor(clust.final.cols, 1),
                                      shape = 21, stroke = 1,
                                      col =  ifelse(dat$mean.pv < 0.05, 'black', 'white')
) + #, size = I(10)
  geom_abline(linetype = 2) +
  geom_errorbar(aes(x = mean.wt, y = mean.mut, xmin = mean.wt.down, xmax = mean.wt.up), width = 0.05, size = 0.6 , col = clust.final.cols) +
  geom_errorbar(aes(x = mean.wt, y = mean.mut, ymin = mean.mut.down, ymax = mean.mut.up), width = 0.05, size = 0.6, col = clust.final.cols) +
  #lims(x = c(0, max(dat$sd.wt)), y = c(0, max(dat$sd.mut))) +
  labs(x = paste('mean bias', ', WT', sep = ''), y = paste('mean bias', ',', geno, sep = '')) +
  scale_x_log10() + scale_y_log10() +
  #theme(text = element_text(size = 20)) +
  geom_text(aes(mean.wt, mean.mut, label = state), size = I(5), check_overlap = TRUE) + theme_bw() + theme(text = element_text(size = 20))
pl
#
pdf(paste('figures/fate_mean_', geno, '.pdf', sep = ''), width = 5, height = 5)
pl
dev.off()
# variance scatterplot
pl <- ggplot(data = dat) + geom_point(aes(sd.wt, sd.mut), size = I(2 - 2 * log10(dat$sd.pv + 1e-2)),
                                      fill = adjustcolor(clust.final.cols, 1),
                                      shape = 21, stroke = 1,
                                      col =  ifelse(dat$sd.pv < 0.05, 'black', 'white')
                                      ) + #, size = I(10)
   geom_abline(linetype = 2) +
  geom_errorbar(aes(x = sd.wt, y = sd.mut, xmin = sd.wt.down, xmax = sd.wt.up), width = 0.05, size = 0.6 , col = clust.final.cols) +
  geom_errorbar(aes(x = sd.wt, y = sd.mut, ymin = sd.mut.down, ymax = sd.mut.up), width = 0.05, size = 0.6, col = clust.final.cols) +
   #lims(x = c(0, max(dat$sd.wt)), y = c(0, max(dat$sd.mut))) +
  labs(x = paste('coefficient of variation', ', WT', sep = ''), y = paste('coefficient of variation', ',', geno, sep = '')) +
  scale_x_log10() + scale_y_log10() +
  #theme(text = element_text(size = 20)) +
  geom_text(aes(sd.wt, sd.mut, label = state), size = I(5), check_overlap = TRUE) + theme_bw() + theme(text = element_text(size = 20))
pl
#
pdf(paste('figures/fate_variance_', geno, '.pdf', sep = ''), width = 5, height = 5)
pl
dev.off()
#
#show individual logNormal distributions
pdf(paste('figures/clone_bias_distribution_', geno, '.pdf', sep = ''), height = 4, width = 4)
for (i in 1:length(unique(clust.final))){
  print(i)
  #i <- 4
  m <- dat0$mean.wt[i]; ssd <- sqrt(dat0$sd.wt[i])
  x <- rlnorm(100000, meanlog = m, sdlog = ssd)
  y <- logNorm(m, ssd)
  #
  m1 <- dat0$mean.mut[i]; ssd1 <- sqrt(dat0$sd.mut[i])
  x1 <- rlnorm(100000, meanlog = m1, sdlog = ssd1)
  y1 <- logNorm(m1, ssd1)
#
  xlm <- min(max(quantile(x, probs = 0.98), quantile(x1, probs = 0.98)), 1)
pl <- ggplot() + theme_bw() + #geom_density(aes(x, col = 'WT'), size = I(0.02)) +
  #geom_density(aes(x1, col = geno), size = I(0.02)) +
  labs(x = paste('clone bias, cluster', i), y = '') + #, title = paste('cluster', i)) +
  geom_line(aes(y$x, y$f, col = 'WT'), size = I(3), alpha = I(0.5), linetype = 1) +
  geom_line(aes(y1$x, y1$f, col = geno), size = I(3), alpha = I(0.5), linetype = 1) +
  coord_cartesian(xlim = c(0.00, xlm)) +
  #scale_x_log10() +
  geom_errorbar(aes(x = dat$mean.wt[i], y = max(y$f), xmin = dat$mean.wt[i] - sd(x), xmax = dat$mean.wt[i] + sd(x), col = 'WT'), width = 2, size = I(1), linetype = 1) +
  geom_errorbar(aes(x = dat$mean.mut[i], y = max(y1$f), xmin = mean(x1) - sd(x1), xmax = mean(x1) + sd(x1), col = geno), width = 2, size = I(1), linetype = 1) +
  geom_vline(aes(xintercept = dat$mean.wt[i], col = 'WT'), linetype = 2, size = I(1)) +
  geom_vline(aes(xintercept = dat$mean.mut[i], col = geno), linetype = 2, size = I(1)) +
  scale_color_manual(values = col.genotype[c('WT', geno)]) +
  theme(text = element_text(size = 20), legend.position = c(0.8, 0.5), legend.title = element_blank())
  print(pl)
}
dev.off()

# mean-varaince scatterplot
dat1 <- dat
dat1$state <- as.character(dat1$state)
shape.values <- c('WT' = 21, 'DNMT3A' = 22)
if (geno == 'TET2') shape.values <- c('WT' = 21, 'TET2' = 22)
pl <- ggplot(data = dat1) +
  #geom_point(aes(mean.wt, sd.wt, shape = 'WT', fill = state), #), fill = clust.final.cols,
  #           size = I(4), stroke = 0.1) +
  #geom_point(aes(mean.mut, sd.mut, shape = geno), fill = clust.final.cols, size = I(4), stroke = 0.1) +
  geom_point(aes(mean.wt, sd.wt, fill = as.character(state)), shape = 21, size = I(4), stroke = 0.1) +
  geom_point(aes(mean.mut, sd.mut, fill = as.character(state)), shape = 22, size = I(4), stroke = 0.1) +
  geom_point(aes(mean.wt, sd.wt, shape = 'WT'), size = I(4), stroke = 0.3) +
  geom_point(aes(mean.mut, sd.mut, shape = geno), size = I(4), stroke = 0.3) +
  geom_errorbar(aes(x = mean.wt, y = sd.wt, xmin = mean.wt.down, xmax = mean.wt.up, col = as.character(1:10)), width = 0.0, size = 0.5) +
  geom_errorbar(aes(x = mean.mut, y = sd.mut, xmin = mean.mut.down, xmax = mean.mut.up, col = as.character(1:10)), width = 0.0, size = 0.5) +
  geom_errorbar(aes(x = mean.wt, y = sd.wt, ymin = sd.wt.down, ymax = sd.wt.up, col = as.character(1:10)), width = 0.0, size = 0.5) +
  geom_errorbar(aes(x = mean.mut, y = sd.mut, ymin = sd.mut.down, ymax = sd.mut.up, col = as.character(1:10)), width = 0.0, size = 0.5) +
  #geom_text_repel(aes(mean.mut, sd.mut, label = state, col = state), size = I(5)) +
  geom_segment(aes(x = mean.wt, y = sd.wt, xend = mean.mut, yend = sd.mut),
                   arrow = arrow(length = unit(0.3, "cm")),  size = 0.5, color = "black", linetype = 1
                   #color = clust.final.cols
  ) +
  scale_fill_manual(values = clust.final.cols) +
  scale_shape_manual(values = shape.values)  +
  scale_color_manual(values = clust.final.cols) +
  labs(x = 'mean', y = 'coefficient of variation', title = 'clone fate patterns') +
  scale_x_log10() + scale_y_log10() +
  theme_bw() + theme(text = element_text(size = 20), #legend.position = c(0.5, 0.2),
                     legend.title = element_blank()) + guides(color = FALSE)
pl
#
#
pdf(paste('figures/clone_bias_patterns_', geno, '.pdf', sep = ''), height = 4, width = 5)
pl
dev.off()

pl <- ggplot(dat) + geom_point(aes(mean.pv, sd.pv, col = as.character(state)), size = I(6)) + theme_bw() +
  scale_x_log10() + scale_y_log10() +
  labs(x = 'mean difference, p-value', y = 'CV difference, p-value', title = paste(geno, 'vs WT')) +
  geom_text(aes(x = mean.pv, y = sd.pv, label = as.character(state)), size = I(5)) +
  scale_color_manual(values = clust.final.cols) +
  geom_vline(aes(xintercept = 1e-2), linetype = 2) + geom_hline(aes(yintercept = 1e-2), linetype = 2) +
  theme(text = element_text(size = 20), legend.position = 'none')
pl

pdf(paste('figures/clone_bias_pvalue_', geno, '.pdf', sep = ''), height = 4, width = 4)
pl
dev.off()

