library(ggplot2)
library(data.table)
#library(glmmTMB)
library(parallel)
library(RColorBrewer)

load('data/processed/DNA_clones.RData')
hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
col.genotype <- adjustcolor(c('black', 'red', 'green'), 0.5)
names(col.genotype) <- c('WT', 'DNMT3A', 'TET2')

### compared DNA-seq-derived and scRNA-seq-derived clone abundances
eids <- as.numeric(as.factor(paste(clone.rna.trends[,2], clone.rna.trends[,3], sep = '_')))
cols <- rainbow(length(unique(eids)))[eids]
do.cd34 <- TRUE
par2 <- c(2, 5)
if (do.cd34) par2 <-  c(2,4)
png(ifelse(do.cd34, "figures/comparison_clone_dna.png", "figures/comparison_clone_dna_rna.png"),
    width = 150 * ifelse(do.cd34, 4, 5), height = 150 * 2)
par(mfrow = par2, mar = c(5,5,0.2,0.2))
cex <- 2
cls <- unique(c(rownames(clone.rna.trends), rownames(clone.trends)))
x <- rep(0, length(cls))
names(x) <- cls
if (!do.cd34){
  x[rownames(clone.rna.trends)] <- clone.rna.trends[,1]
  inds <- 1:9
}else{
  x[rownames(clone.trends)] <- clone.trends.norm[, 'CD34']
  inds <- c(1:5, 7:9)
}
for (i in inds){
  y <- rep(0, length(cls))
  names(y) <- cls
  y[rownames(clone.trends)] <- clone.trends.norm[, i]
  smpl <- sample(1:length(x))
  print(paste(i, cor(log(x[smpl] + min(x[x>0])*8e-1), log(y[smpl] + min(y[y>0])*8e-1), method = 'pearson')))
  plot(x[smpl] + min(x[x>0])*8e-1, y[smpl] + min(y[y>0])*8e-1,   cex.lab = cex, cex.axis = cex, pch = 21, cex = 1,
       log = 'xy',
       xlab = ifelse(!do.cd34, 'scRNA-seq', 'CD34, DNA'), ylab = paste(colnames(clone.trends)[i], ', DNA', sep = ''),
       xaxt = 'n', yaxt = 'n', xlim = c(min(x[x>0])*8e-1, 1), ylim = c(min(y[y>0])*8e-1, 1),
       bg = adjustcolor(cols[smpl], 0.4), col = adjustcolor('black', 0.2))
  axis(1, at = c(1e-3, 1e-1), labels = c(1e-3, 1e-1), cex.axis = 1.5)
  axis(2, at = c(1e-3, 1e-1), labels = c(1e-3, 1e-1), cex.axis = 1.5)
}
if (!do.cd34){
y <- rep(0, length(cls)); names(y) <- cls
#y[rownames(clone.trends)] <- rowMeans(clone.trends.norm[,1:9])
y[rownames(clone.trends)] <- rowSums(clone.trends[,c( 'Myeloid', 'B', 'T')])
y <- y / sum(y)
smpl <- sample(1:length(x))
plot(x[smpl] + min(x[x>0])*5e-1, y[smpl] + min(y[y>0])*5e-1,   cex.lab = cex, cex.axis = cex, pch = 21, cex = 1, log = 'xy',
     xlab = ifelse(!do.cd34, 'scRNA-seq', 'CD34, DNA'), ylab = paste('total, DNA', sep = ''),
     xaxt = 'n', yaxt = 'n',
     bg = adjustcolor(cols[smpl], 0.4), col = adjustcolor('black', 0.2))
axis(1, at = c(1e-3, 1e-1), labels = c(1e-3, 1e-1), cex.axis = 1.5)
axis(2, at = c(1e-3, 1e-1), labels = c(1e-3, 1e-1), cex.axis = 1.5)
}
dev.off()



### show time-series heatmap
inds <- rownames(clone.trends) %in% names(clone_state)
xx <- clone.trends.norm[rownames(clone.trends.norm) %in% names(clone_state), 1:9]
xx <- xx * 1e+6 + 1
lbl <- 10^c(-6,-4,-2,-1)
annot <- cbind(clone.trends[inds, c('genotype'),drop=F], clone_state)
annot <- cbind(clone_state, clone.trends[inds, c('genotype'),drop=F])
colnames(annot)[1] <- 'status'
pheat <- pheatmap::pheatmap( log(xx[,])^(1/1),
                             color = colorRampPalette(c('lightgrey','darkred'))(100),
                    fontsize = 25, #fontsize_row = 20, fontsize_col = 20,
                    angle_col = 45, cellwidth = 50,
                    cluster_cols = F, cluster_rows = T,
                    show_rownames = F,
                    annotation_row =  annot,
                    annotation_names_row = F,
                    legend_breaks = log(lbl * 1e+6 + 1),
                    legend_labels = lbl,
                    annotation_colors = list(genotype = c(WT = adjustcolor("black", 0.5), DNMT3A = adjustcolor("red", 0.5), TET2 = adjustcolor('green', 0.5) ),
                                             status = c(early = 'grey80', late = 'grey20')),
                    treeheight_row = 5, annotation_legend=T)

png('figures/dna_clone_heatmap_recolor.png', width = 800, height = 400)
pheat
dev.off()

### estimate number of dominant clones at every time point:
inds <- rownames(clone.trends) %in% names(clone_state)
xx <- as.data.frame(cbind(peak = as.numeric(colnames(clone.trends.norm)[1:5][apply(clone.trends.norm[inds,1:5], 1, which.max)]),
                          clone.trends[inds, c('genotype', 'sample')], exp = clone.trends[inds,]$`hash.match$exp[[i]]`))
xx$count <- 1
#clone_dominant <- clone_dominant[,names(col.genotype)]
#rownames(clone_dominant) <- tm
yy <- aggregate(count ~ peak + genotype + sample + exp, data = xx, sum)
pl <- ggplot(yy, aes(x = as.factor(peak), y = count, col = genotype, shape = exp)) +
  scale_y_log10() +
  scale_color_manual(values = col.genotype) +
  labs(x = 'time', y = 'number of peaked clones') +
  geom_boxplot(position=position_dodge(0.5), width = 0.3, alpha = 1, lwd = 0.3) +
  geom_point(position=position_dodge(0.5), size = 3, alpha = 0.5) +
  theme_bw() + theme(text = element_text(size = 20), legend.position = c(0.7, 0.8), legend.box = "horizontal")
pl
#
pdf('figures/time_clone_peaks.pdf', width = 6, height = 4)
pl
dev.off()

### estimate number of clones
y <- aggregate(as.matrix(clone.trends.norm[,1:9] >= 1e-4) ~ genotype + sample + `hash.match$exp[[i]]`,
               data = clone.trends[,], sum)
colnames(y)[3] <- 'exp'
y$exp[y$exp == '123'] <- 'exp 1'; y$exp[y$exp == '127'] <- 'exp 2'

z <- reshape2::melt(y, id = 1:3)
#
pl <- ggplot(z, aes(x = variable, y = value)) +
  labs(x = '', y = 'barcodes detected') + #ylim(0,1) +
  geom_col(aes(fill = genotype), col = 'black', lwd = 0.1, position = position_dodge2(), width = 0.9) +
  scale_fill_manual(values = col.genotype) +
  theme_bw() + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 30, vjust = 0.6),
                     legend.title = element_text(size = 0), legend.position = c(0.8, 0.8), legend.box = "vertical") +
  facet_wrap(~ exp, nrow = 2)  +
  theme(strip.text = element_blank())
pl

pdf("/home/solrust/CH/draft/figures/dna_clones/clone_number.pdf", width = 5, height = 5)
pl
dev.off()

t.test(z[z$variable == '6' & z$exp == 'exp 1' & z$genotype == 'WT',]$value,
       z[z$variable == '6' & z$exp == 'exp 1' & z$genotype != 'WT',]$value
)


###
show.HSPC <- F
#for (show.HSPC in c(T, F)){
show.HSPC <- T
print(show.HSPC)
cutfrac <- 1e-4#1e-4
ind1 <- clone.trends.norm$CD34 >= cutfrac
y1 <- aggregate(as.matrix(clone.trends.norm[ind1,1:9] >= cutfrac) ~ genotype + sample + `hash.match$exp[[i]]`,
                data = clone.trends[ind1,], sum)
colnames(y1)[3] <- 'exp'
y1$exp[y1$exp == '123'] <- 'exp 1'; y1$exp[y1$exp == '127'] <- 'exp 2'
#
if (show.HSPC){
  y1[,-c(1:3)] <- y1[,-c(1:3)] / y$CD34
  yname <- '% barcodes shared of HSPC'
}else{
  y1[,-c(1:3)] <- y1[,-c(1:3)] / y[,-c(1:3)]
  yname <- '% barcodes shared with HSPC'
}
z <- reshape2::melt(y1, id = 1:3)
#
pl <- ggplot(z, aes(x = variable, y = value)) +
  labs(x = '', y = '% barcodes shared with HSPC') + #ylim(0,1) +
  geom_col(aes(fill = genotype), col = 'black', lwd = 0.1, position = position_dodge2(), width = 0.9) +
  scale_fill_manual(values = col.genotype) +
  theme_bw() + theme(text = element_text(size = 15), axis.text.x = element_text(angle = 30, vjust = 0.6),
                     legend.title = element_text(size = 0), legend.position = 'none', #c(0.15, 0.85),
                     legend.box = "vertical") +
  facet_wrap(~ exp, nrow = 2)  +
  theme(strip.text = element_blank())
print(pl)

pdf(paste("figures/clone_intersection_", ifelse(!show.HSPC, 'with', 'of'),"_HSPC.pdf", sep = ''), width = 5, height = 5)
pl
dev.off()
#}


# estimate trends
trend_list <- lapply(names(growth.info.nbinom), function(nm){
  trend <- clone.trends[nm, 1:5]/sample_size[clone.trends[nm,]$sample, 1:5]
  ftcoef <- growth.info.nbinom[[nm]]$coefs[,1]
  pred <- exp(ftcoef[1] + ftcoef[2] * tm)
  list(trend = trend, pred = pred, pv = growth.info.nbinom[[nm]]$coefs[2,4])
})
names(trend_list) <- names(growth.info.nbinom)

### visualize a few clones
clones.late <- names(which(clone_state == 'late'))
tm <- c(6,12, 16, 20, 24)
pvs <- unlist(lapply(growth.info.nbinom, function(x) x$coefs[2,4]))
cls <- clones.late[pvs[clones.late] & rowSums(clone.trends[clones.late,1:5] > 4) > 2]
cls <- names(pvs)[pvs & rowSums(clone.trends[names(pvs),1:5] > 4) > 2]
cls <- intersect(cls, names(growth.info.nbinom))
str(cls)
cls.up <- cls[rate[cls] > 0]
str(cls.up)
cls.down <- cls[rate[cls] < 0]
str(cls.down)

n.clones <- 4
fctr <- 1e-2
cols <- adjustcolor(brewer.pal(n.clones, ifelse(T, "Set2", "Set3") ), 1)

cls.up.show <- c("12_TTTATGTCGACACAGGACCTGTGACTCTGTGTAGAGGAAA", "4_TGCATGCTGACACCCTACCCGGGAGGGCGTGACGAGGCGG", "10_GAGGTGGACTCAAGATACCTGAGATCAGGTACATAGTAAT", "1_TCCATGAAACCACTTAACGTATGAAATAGTTCAGAGTTTC")
cls.down.show <- c("2_GTTATGCGGTCAAGATACTAACGACTCCGTGATGAGATTC", "6_CCAATGACCACACGTGACCATAGAAAGTGTAAATAGTCTT", "11_CATTTGTGACCACGGTACTTAGGACGAAGTTACTAGGATA", "8_AAACTGCTCACACCAAACTACCGATAAAGTAAAGAGCCGG")

do.up <- F
if (do.up){
  cls <- cls.up
  cls.show <- cls.up.show
  #nms.show <- paste('clone ', substr(cls.up.show, 1, 6))
  nms.show <- substr(cls.up.show, 1, 6)
}else{
  cls <- cls.down
  cls.show <- cls.down.show
  #nms.show <- paste('clone ',4:6)
  nms.show <- substr(cls.down.show, 1, 6)
}

pdf(paste('figures/clone_', ifelse(do.up, 'up', 'down') ,'_examples_v1.pdf', sep = ''),
    width = 4, height = 4)
par(mfrow=c(1,1), mar = c(4.5,4.5,0,0))
plot(tm, rep(0, length(tm)), ylim = c(0, ifelse(do.up, 0.025, 0.005)),
     xlim = c(tm[1]-1, 27), bty = 'n', col = 'white', cex.lab = 2, cex.axis = 1.5, xaxt = 'n',
     cex = cex, pch = 21, lwd = 1.5, xlab = 'time', ylab = 'frequency',
     #main = ifelse(do.up, 'relative growth', 'relative reduction'),
     font.main = 1, cex.main = cex,
)
axis(1, at = c(0, tm), labels = c(0, tm), cex.axis = 1.5)
for (i in sample(1:length(cls), 200)){
  nm <- cls[i]
  points(tm, trend_list[[nm]]$pred * fctr, col = adjustcolor('grey', 0.5), type = c('l'), lwd = 1, cex = 1)
}
for (i in 1:length(cls.show)){
  nm <- cls.show[i]
  points(tm, trend_list[[nm]]$pred * fctr, col = cols[i], type = c('l'), cex = cex, pch = 21, lwd = 6)
  points(tm, trend_list[[nm]]$trend * fctr, bg = cols[i], col = 'black', type = c('p'), cex = cex*0.5, pch = 21, lwd = 1)
}
for(i in 1:length(cls.show)){
  nm <- cls.show[i]
  j <- ifelse(do.up, 5, 1)
  shift <- ifelse(do.up, 2, -0.0)
  text(tm[j] + shift, trend_list[[nm]]$pred[j] * fctr + 0.02 * fctr, round(rate[nm], 2),
       col = cols[i], cex = 1.0)
}
pos.legend <- ifelse(do.up, 'topleft', 'topright')
legend(pos.legend, title = 'clone', title.cex = 1.5, cex = 1.0, title.adj = 0.2,
       nms.show, col = cols, pch = 19, pt.cex = 1.5, bty = 'n')
dev.off()

### visualize distribution of growth rates across CD34+ and CD34- clones
do.revised <- TRUE
if (!do.revised){
  #tt <- data.frame(rate = rate[cls], CD34 = clone.trends.norm[cls,]$CD34, genotype = clone.trends[cls, ]$genotype)
  tt <- data.frame(rate = rate[cls], clone.trends.norm[cls, c('CD34', 'Myeloid', 'B', 'T')], genotype = clone.trends[cls, ]$genotype)
  rate_categories <- cut(tt$rate, breaks = quantile(tt$rate, probs = seq(0, 1, length.out = 6), na.rm = TRUE),  include.lowest = T)
  tt$rate_scores <- (as.numeric(rate_categories))
  tt$rate_categories <- as.factor(as.numeric(rate_categories))
  text.scale <- 25
  pl <- ggplot(tt, aes(x = CD34 > 1e-5, y = rate)) + geom_boxplot(col = 'black') + geom_jitter(aes(bg = genotype), size = I(3), pch = 21, col = 'black') + theme_bw() +
    labs(y = 'clone growth rate', x = '') + scale_fill_manual(values = col.genotype) + geom_hline(yintercept = 0, linetype = 2, lwd = 1) +
    scale_x_discrete(labels = c("FALSE" = "CD34-", "TRUE" = "CD34+")) +
    theme(axis.title = element_text(size = text.scale), axis.text = element_text(size = text.scale),
          legend.title = element_blank(),  legend.text = element_text(size = text.scale), legend.position = c(0.8, 0.15))
  pl
  pdf('figures/growth_heterogeneity_boxplot.pdf')
  pl
  dev.off()
}else{
  tt <- data.frame(
    rate = rate[cls],
    clone.trends.norm[cls, c("CD34", "Myeloid", "B", "T")],
    genotype = clone.trends[cls, ]$genotype
  )

  # Preserve the original growth-rate categories used downstream
  rate_categories <- cut(
    tt$rate,
    breaks = quantile(
      tt$rate,
      probs = seq(0, 1, length.out = 6),
      na.rm = TRUE
    ),
    include.lowest = TRUE
  )

  tt$rate_scores <- as.numeric(rate_categories)
  tt$rate_categories <- as.factor(as.numeric(rate_categories))
  median_rate_lookup <- tapply(
    tt$rate,
    tt$rate_categories,
    median,
    na.rm = TRUE
  )

  # Add the corresponding median rate to every row
  tt$median_rate <- unname(
    median_rate_lookup[as.character(tt$rate_categories)]
  )

  # Add CD34-status variable for plotting
  tt$CD34_status <- ifelse(
    tt$CD34 > 1e-5,
    "CD34+",
    "CD34-"
  )

  tt$CD34_status <- factor(
    tt$CD34_status,
    levels = c("CD34-", "CD34+")
  )

  # Set genotype order
  tt$genotype <- factor(
    tt$genotype,
    levels = c("WT", "DNMT3A", "TET2")
  )

  # Validate labels before plotting
  if (anyNA(tt$genotype)) {
    stop(
      "Unexpected genotype labels. Labels present before factor conversion: ",
      paste(unique(clone.trends[cls, ]$genotype), collapse = ", ")
    )
  }

  print(table(tt$CD34_status, tt$genotype, useNA = "ifany"))

  text.scale <- 25

  pl <- ggplot(
    tt,
    aes(x = genotype, y = rate)
  ) +
    geom_boxplot(
      color = "black",
      outlier.shape = NA,
      width = 0.65
    ) +
    geom_point(
      aes(fill = genotype),
      shape = 21,
      size = 2,
      color = "black",
      position = position_jitter(
        width = 0.2,
        height = 0,
        seed = 123
      )
    ) +
    facet_grid(
      cols = vars(CD34_status),
      scales = "fixed"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = 2,
      linewidth = 1
    ) +
    scale_fill_manual(
      values = col.genotype,
      drop = FALSE
    ) +
    labs(
      y = "clone growth rate",
      x = NULL
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = text.scale),
      axis.text = element_text(size = text.scale),
      strip.text = element_text(size = text.scale),
      strip.background = element_rect(fill = "white"),
      panel.spacing.x = grid::unit(1, "lines"),
      legend.position = "none"
    )

  print(pl)

  ggsave(
    filename = "figures/growth_heterogeneity_boxplot.pdf",
    plot = pl,
    width = 10,
    height = 5
  )
}

## visualize dependence of CD34+/B/T/Myeloid on growth rate
do.revised
if (!do.revised){
  vars <- c('CD34', 'Myeloid', 'B', 'T')
  file.nm <- paste('figures/growth_lineages_boxplot.pdf', sep = '')
  pdf(file.nm)
  for (var in vars){
    pl <- ggplot(tt[tt$rate >= -1000000 & !is.na(tt$rate_categories),], aes(x = rate_categories, y = .data[[var]])) + geom_boxplot(col = 'black') +
      geom_jitter(aes(bg = genotype), size = I(3), pch = 21, col = 'black') +
      geom_smooth(aes(rate_scores, .data[[var]], col = genotype), method = "lm") +
      scale_y_log10() +
      theme_bw() +
      labs(y = paste(var, '+ fraction', sep = ''), x = 'clone growth rate') +
      scale_colour_manual(values = col.genotype) + scale_fill_manual(values = col.genotype) +
      #scale_x_discrete(labels = round(sort(rate_unique), 2)) +
      theme(axis.title = element_text(size = 30), axis.text = element_text(size = text.scale),
            legend.title = element_blank(),  legend.text = element_text(size = text.scale),
            legend.position = c(0.2, 0.9))
    print(pl)
  }
dev.off()
}else{
  vars <- c("CD34", "Myeloid", "B", "T")
  file.nm <- "figures/growth_lineages_boxplot_by_genotype.pdf"

  tt$genotype <- factor(
    tt$genotype,
    levels = c("WT", "DNMT3A", "TET2")
  )

  # Median-rate labels in rate-category order
  x_breaks <- seq_along(median_rate_lookup)
  x_labels <- sprintf("%.2f", as.numeric(median_rate_lookup))

  pdf(file.nm, width = 8, height = 5)

  for (var in vars) {

    plot.df <- tt[
      tt$rate >= -1000000 &
        !is.na(tt$rate_categories),
    ]

    pl <- ggplot(
      plot.df,
      aes(
        x = rate_scores,
        y = .data[[var]],
        group = interaction(rate_scores, genotype)
      )
    ) +
      geom_boxplot(
        aes(color = genotype),
        position = position_dodge(width = 0.75),
        outlier.shape = NA,
        width = 0.6
      ) +
      geom_point(
        aes(fill = genotype),
        shape = 21,
        size = 2,
        color = "black",
        position = position_jitterdodge(
          jitter.width = 0.15,
          jitter.height = 0,
          dodge.width = 0.75,
          seed = 123
        )
      ) +
      geom_smooth(
        aes(
          x = rate_scores,
          y = .data[[var]],
          color = genotype,
          group = genotype
        ),
        method = "lm",
        se = TRUE,
        inherit.aes = FALSE
      ) +
      scale_y_log10() +
      scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels
      ) +
      theme_bw() +
      labs(
        y = paste0(var, "+ fraction"),
        x = "clone growth rate"
      ) +
      scale_colour_manual(values = col.genotype) +
      scale_fill_manual(values = col.genotype) +
      theme(
        axis.title = element_text(size = 30),
        axis.text = element_text(size = text.scale),
        legend.title = element_blank(),
        legend.text = element_text(size = text.scale),
        legend.position = c(0.2, 0.85)
      )

    print(pl)
  }

  dev.off()
}

### estimate number of dominant clones at every time point:
# choose clones seen in at least 2 time points
#clones.persist <- clones.late[rowSums(clone.trends[clones.late,1:5] > 2) > 2]
clones.persist <- clones.late[rowSums(clone.trends[clones.late,1:5] > 2) > 2]
str(clones.persist)

show_increasing <- 1 #-1
for (show_increasing in c(-1, 1)){
clone_grow <- table(as.data.frame(cbind(rate[clones.persist]*(show_increasing) >= 0.0, clone.trends[clones.persist,c('sample')] )))
clone_freq <- table(as.data.frame(cbind((clone.trends.norm[clones.persist, '24'] - 1e-5)*sign(show_increasing) >= 0, clone.trends[clones.persist,c('sample')] )))
clone_stat <- as.data.frame(cbind(clone_grow[2,], clone_freq[2,]))
#clone_growing <- table(as.data.frame(cbind(clone.trends.norm[, '24'] > 1e-4, clone.trends[,c('sample')] )))
clone_stat <- as.data.frame(cbind(clone_stat, hash.match[rownames(clone_stat),]$genotype, hash.match[rownames(clone_stat),]$exp))
colnames(clone_stat) <- c('rate', 'freq', 'geno', 'exp')
clone_stat
clone_stat_up <- clone_stat
#
nm <- ifelse(show_increasing == 1, 'increase', 'decrease')
pdf(paste('figures/clone_numbers_', nm,'.pdf', sep = ''), width = 7, height = 3.5)
cex <- 1.4
par(mfrow=c(1,2), mar = c(5,5,1,1))
nums0 <- nums <- as.numeric(factor(clone_stat$geno, levels = c('DNMT3A', 'WT', 'TET2')))
ind1 <- clone_stat$geno == 'WT' & hash.match[rownames(clone_stat),]$exp == '123'
nums[ind1] <- nums[ind1] - 0.2
ind2 <- clone_stat$geno == 'WT' & hash.match[rownames(clone_stat),]$exp == '127'
nums[ind2] <- nums[ind2] + 0.2
plot(nums, clone_stat$rate, cex = 2, cex.lab = cex, ylab = ifelse(!show_increasing, 'clones, decreasing', 'clones, growing'),
     cex.axis = cex, xlab = '', bty = 'n',
     xlim = c(0.8, 3.2), ylim = c(0, max(clone_stat$rate)+1),
     bg = col.genotype[clone_stat$geno], col = 'black', pch = ifelse(clone_stat$exp == '123', 21, 24), xaxt = 'n')
legend('bottomright', c('exp1', 'exp2'), pch = c(21, 24), cex = cex, bty = 'n')
axis(1, labels = c(as.factor(clone_stat$geno)), at = nums0,
     cex.axis = cex)
#
plot(nums, clone_stat$freq, cex = 2, cex.lab = cex, cex.axis = cex, ylab = ifelse(!show_increasing, 'clones, rare', 'clones, frequent'),
     xlab = '', bty = 'n',
     xlim = c(0.8, 3.2), ylim = c(0, max(clone_stat$freq)+1), bg = col.genotype[clone_stat$geno], col = 'black', pch = ifelse(clone_stat$exp == '123', 21, 24), xaxt = 'n')
#legend('topright', c('exp 1', 'exp 2'), pch = c(21, 24), cex = 1.5)
axis(1, labels = c(as.factor(clone_stat$geno)), at = nums0, srt = 45,
     cex.axis = cex)
dev.off()
}

# estiamte p-values
ix <- hash.match[rownames(clone_stat),]$exp %in% '123'
swt <- clone_stat[ix & clone_stat$geno == 'WT',]$rate
smut <- clone_stat[ix & clone_stat$geno != 'WT',]$rate
#wilcox.test(swt, smut)
t.test(swt, smut)
swt
smut

### show compartment heatmap
inds <- clones.late
str(inds)
inds <- rownames(clone.trends)[rowSums(clone.trends[,c("CD34", "Myeloid", "B", "T")]) > 100]
str(inds)
inds <- intersect(inds, names(rate))
str(inds)
xx <- clone.trends.norm[inds, c("CD34", "Myeloid", "B", "T")]
xx <- xx * 1e+6 + 1
lbl <- 10^c(-6,-4,-2,-1)
annot <- cbind(clone.trends[inds, c('genotype'),drop=F], ifelse(rate[inds] > 0 , 'up', 'down'), clone.trends[inds, c('sample'),drop=F])
colnames(annot) <- c('genotype', 'growth', 'sample')
padj <- p.adjust(pvs[rownames(annot)], 'BH')
annot$growth[ padj > 0.2 ] <- 'uncertain'
inds1 <- clone.trends.norm[inds,]$CD34 > 1e-5
hcl1 <- hclust(dist(log(xx[inds1,])))
hcl2 <- hclust(dist(log(xx[!inds1,])))
ord <- c(which(inds1)[hcl1$order],
         which(!inds1)[hcl2$order])
pheat <- pheatmap::pheatmap( log(xx[ord,]),#[clone.trends.norm[inds,]$CD34 > 0e-5][ord,])^(1/1),
                    fontsize = 25, #fontsize_row = 20, fontsize_col = 20,
                    angle_col = 45, cellwidth = 50,
                    cluster_cols = F, cluster_rows = F,
                    show_rownames = F,
                    annotation_row =  annot[ord, -3],
                    annotation_names_row = F,
                    legend_breaks = log(lbl * 1e+6 + 1),
                    legend_labels = lbl,
                    annotation_colors = list(genotype = c(WT = adjustcolor("black", 0.5), DNMT3A = adjustcolor("red", 0.5), TET2 = adjustcolor('green', 0.5) ),
                                             growth = c(uncertain = 'grey80', up = 'red', down = ' blue')
                    ),
                    treeheight_row = 5, annotation_legend=T)
#
pdf('figures/dna_clone_heatmap_lineages.pdf', width = 5*2, height = 5*1)
pheat
dev.off()

mean(xx[annot$growth == 'up',]$CD34 > 1e+2)
mean(xx[annot$growth == 'up',]$Myeloid > 1e+2 & xx[annot$growth == 'up',]$B > 1e+2)
#
mean(xx[annot$growth == 'down',]$CD34 > 1e+2)
mean(xx[annot$growth == 'down',]$Myeloid > 1e+2 & xx[annot$growth == 'down',]$B > 1e+2)

### plot distinct behaviour of CD34+ and CD34- mutant clones
yy <- data.frame(xx, annot, exp = hash.match[annot$sample,]$exp)
z1 <- aggregate( cbind(CD34 < 5 & B > 10, CD34 < 5 & Myeloid > 10) ~ genotype + sample + exp, data = yy, sum)
z1.tot <- aggregate( cbind(CD34 < 5 & B > -10, CD34 < 5 & Myeloid > -10) ~ genotype + sample + exp, data = yy, sum)
z2 <- aggregate( cbind(CD34 >= 5 & B > 10, CD34 >= 5 & Myeloid > 10) ~ genotype + sample + exp, data = yy, sum)
z2.tot <- aggregate( cbind(CD34 >= 5 & B > -10, CD34 >= 5 & Myeloid > -10) ~ genotype + sample + exp, data = yy, sum)
#
z1[,c('V1', 'V2')] <- z1[,c('V1', 'V2')]  / z1.tot[,c('V1', 'V2')]
z2[,c('V1', 'V2')] <- z2[,c('V1', 'V2')]  / z2.tot[,c('V1', 'V2')]
#
z <- rbind(cbind(z1, compartment = 'CD34-') , cbind(z2, compartment = 'CD34+'))
pl <- ggplot() + geom_point(aes(V1, V2, col = genotype, shape = exp), data = z, size = 4) +
  scale_color_manual(values = col.genotype) +
  labs(x = 'B+ clones', y = 'Myeloid+ clones') + theme_bw() + theme(text = element_text(size = 20)) +
  facet_grid(~ compartment)
pl
#
pdf('figures/comparison_B_My_fraction.pdf', width = 6, height = 3)
pl
dev.off()

### show lineage correlations
inds1 <- inds[clone.trends.norm[inds,]$CD34 > 0e-5]
cr1 <- cor(clone.trends.norm[inds1,c('T', 'B', 'Myeloid')], as.numeric(rate[inds1] > 0),
           method = 'spearman')
inds1 <- inds[clone.trends.norm[inds,]$CD34 == 0e-5]
cr2 <- cor(clone.trends.norm[inds1,c('T', 'B', 'Myeloid')], as.numeric(rate[inds1] > 0),
           method = 'spearman')
cr <- cbind(cr1, cr2)
colnames(cr) <- c('CD34+', 'CD34-')
pheat <- pheatmap::pheatmap( t(cr),
                    main = 'growth vs lineage',
                    fontsize = 20, fontsize_row = 20, fontsize_col = 20,
                    angle_col = 45, cellwidth = 50,
                    cluster_cols = F, cluster_rows = F,
                    show_rownames = T,
                    annotation_names_row = F,
                    #legend_breaks = log(lbl * 1e+6 + 1),
                    #legend_labels = lbl,
                    #annotation_colors = list(genotype = c(WT = adjustcolor("black", 0.5), DNMT3A = adjustcolor("red", 0.5), TET2 = adjustcolor('green', 0.5) ),
                    #                         growth = c(uncertain = 'grey80', up = 'red', down = ' blue')
                    #),
                    treeheight_row = 5, annotation_legend=T)
#
png('figures/dna_clone_heatmap_lineage_cor.png', width = 600, height = 300)
pheat
dev.off()

### examine lineage bias
bias.mat  <- data.frame(clone.trends.norm[clones.late, c('Myeloid', 'B', 'T')], exp = clone.trends[clones.late,]$`hash.match$exp[[i]]`,
                        geno = clone.trends[clones.late,]$genotype, CD34 = ifelse(clone.trends.norm[clones.late, ]$CD34 > 0, 'CD34+', 'CD34-'))
pl <- ggplot(bias.mat, aes(Myeloid, B, col = geno, shape = CD34)) + theme_bw() + geom_point(size = 3) +
  scale_x_log10() + scale_y_log10() + geom_abline(linetype = 2) +
  scale_shape_manual(values = c('CD34+' = 19, 'CD34-' = 4)) +
  scale_color_manual(values = col.genotype) +
  labs(x = 'Myeloid frequency', y = 'B frequency') + facet_wrap(~ exp , #scales = 'free',
                                                                nrow = 2, labeller = labeller(exp = c("123" = "DNMT3A", "127" = "TET2"))) +
  theme(text = element_text(size = 30), legend.title = element_blank())
pl

pdf('figures/lineage_bias_My_B.pdf')#, width = 6, height = 8)
pl
dev.off()

