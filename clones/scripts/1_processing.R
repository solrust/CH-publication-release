library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(tidyr)
library(confintr)
#library(reshape2)

# TODO: change pathway after cleaning up the raw data folder
path <- "data/GS175/"

### upload metadata
metadata <- data.frame(genotype = c(rep('AAVS1', 6), rep('TET2', 6), rep('AAVS1', 8), rep('TET2', 8)),
                       replicate = c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4) ),
                       population = c('ery progenitors', 'neu progenitors', 'bulk', 'ery progenitors', 'neu progenitors', 
                                      'bulk', 'ery progenitors', 'neu progenitors', 'bulk', 'ery progenitors', 'neu progenitors',
                                      'bulk', 'early my progentors', 'early ery progentors', 'committed my progentors', 'committed ery progentors',
                                      'early my progentors', 'early ery progentors', 'committed my progentors', 'committed ery progentors', 'early my progentors',
                                      'early ery progentors', 'committed my progentors', 'committed ery progentors', 'early my progentors',
                                      'early ery progentors', 'committed my progentors', 'committed ery progentors'))
rownames(metadata) <-  paste(metadata$genotype, metadata$replicate, 1:nrow(metadata), sep = '_')
metadata <- read.csv(paste(path, 'GS_175_larry in vitro 3.csv', sep = ''))
metadata <- metadata[,-1]
x <- strsplit(metadata$genotype, '_')
metadata$genotype <- unlist(lapply(x, function(y) y[[1]]))
metadata$replicate <- unlist(lapply(x, function(y) substr(y[[2]],1,1)))
metadata$technical <- unlist(lapply(x, function(y) substr(y[[2]],2,2)))
rownames(metadata) <- paste(metadata$genotype, 1:nrow(metadata), sep = '_')
head(metadata)

### upload barcode counts
files_to_upload <- list.files(path = paste(path, 'starcode/', sep = ''), recursive = TRUE, full.names = TRUE)
#files_to_upload <- list.files(path = paste(path, 'txt/', sep = ''), recursive = TRUE, full.names = TRUE)
files_num <- as.numeric(unlist(lapply(files_to_upload, function(x){
  y <- unlist(strsplit(x, '175-'))[2]
  unlist(strsplit(y, '_'))[1]
})))
# clone names
fls <- lapply(files_to_upload[order(files_num)], function(fl) read.table(fl))
clone_names <- unlist(lapply(fls, function(dat) dat$V1))
clone_names <- unique(unique(clone_names))
str(clone_names)
# raw clone counts
counts_raw <- data.frame(matrix(0, nrow = length(clone_names), ncol = length(fls)))
rownames(counts_raw) <- clone_names
colnames(counts_raw) <- rownames(metadata) 
head(counts_raw)
for (i in 1:ncol(counts_raw)) counts_raw[fls[[i]]$V1 , i] <- fls[[i]]$V2
qplot(rowSums(counts_raw), log = 'x')
# filter out low count clones
clone.coverage <- 20
ind <- rowSums(counts_raw) >= clone.coverage 
table(ind)
counts <- counts_raw[ind, ]
pheat <- pheatmap::pheatmap(log(counts[rowSums(counts) >= 100, ]+1), cluster_cols = TRUE, show_rownames = FALSE,
                            annotation_col = metadata[,c(1,2,3,5,6)])
#
pdf('figures/heatmap_clones.pdf')
pheat
dev.off()


### estimate clonal correlations
estimate_clonal_cor <- function(counts, metadata, popl, pop_background, geno, day, repl, cut_count = 4){
  ind <- metadata$population_short == popl & metadata$genotype == geno & 
    metadata$replicate %in% repl & metadata$timepoint %in% day
  ind1 <- metadata$population_short == pop_background & metadata$genotype == geno & 
    metadata$replicate %in% repl & metadata$timepoint %in% day
  
  x1 <- counts[, ind]
  x2 <- counts[, ind1]
  ii <- x1[,1] + x2[,1] >= cut_count & x1[,2] + x2[,2] >= cut_count 
  
  r1 <- x1[ii,1] / (x1[ii,1] + x2[ii,1])
  r2 <- x1[ii,2] / (x1[ii,2] + x2[ii,2])
  
  corest <- ci_cor(r1, r2, probs = c(0.05, 0.95), method = c("pearson"), type = c("normal"))
  list(estimate = corest$estimate, interval = corest$interval, r1 = r1, r2 = r2)
}
# a routine to restructure 'crop' table for plotting
restructure_table <- function(crpop){
  df <- data.frame(group = rownames(crpop),
                   mean1 = crpop[,1], lower1 = crpop[,2], upper1 = crpop[,3], mean2 = crpop[,4], lower2 = crpop[,5],  upper2 = crpop[,6])
  # Reshape for ggplot
  df_long <- df %>% pivot_longer(cols = -group, names_to = c("stat", "group_id"), names_pattern = "(mean|lower|upper)(\\d)") %>%
    pivot_wider(names_from = stat, values_from = value)
  df_long$group_id <- factor(df_long$group_id, labels = c("WT", "TET2"))
  return(df_long)
}
# a routine to plot heritability correlation barplots
plot_herit_barplots <- function(df_long){
  ggplot(df_long, aes(x = group, y = mean, fill = group_id)) +
    geom_col(position = position_dodge(0.6), width = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.6), width = 0.2) +
    scale_fill_manual(values = c("grey50", "firebrick")) +
    labs(x = "", y = "Clonal correlation", fill = "Genotype") +
    ggtitle(paste("Heritability of lineage bias, replicate", repl)) +
    theme_minimal(base_size = 20) +  # base font size
    theme(legend.position = c(0.4, 0.95),  # top-right corner inside the plot
          legend.justification = c("right", "top"),
          legend.background = element_rect(fill = "white", color = "black"),
          plot.title = element_text(size = 16, face = "bold"),
          #axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.title = element_text(size = 13), legend.text = element_text(size = 12)
    ) 
}
# a routine to plot scatterplot of clonal bias correlations between technical splits
plot_clonal_cor <- function(corest, popl, geno, day, repl){
  r1 <- corest$r1; r2 <- corest$r2
  r_value <- round(corest$estimate, 2)
  r_label <- paste0("r = ", r_value)
  pl <- ggplot(, aes(x = r1, y = r2)) + geom_point(size = 5, shape = 21, fill = 'grey') +
    labs(x = paste(geno, 'clonal bias, split 1'), y = paste(geno, 'clonal bias, split 2'),
         title = paste(popl, ' bias; ', day, ', ', 'rep ', repl, sep = '')) +
    annotate("text", x = Inf, y = -Inf, 
             label = r_label, 
             hjust = 1.1, vjust = -2, # Shifts text slightly inward from the borders
             size = 10, fontface = "italic", color = "black") + 
    theme_bw() + theme(axis.title = element_text(size = 20), title = element_text(size = 15))
  return(pl)
}

pops <- c("bulk CD34neg", "early_My", "early_Ery", "early_undefined")
pop_background <- 'bulk'

for (repl in c(1, 2)){
  for (day in c('day 5', 'day 10')){
    path_plot <- paste('figures/clonal_bias_scatterplots_', gsub(" ", "_", day), '_repl', repl,'.pdf', sep = '')
    pdf(path_plot)
    for (popl in pops){
      for (geno in unique(metadata$genotype)){
        corest <- estimate_clonal_cor(counts, metadata, popl, pop_background, geno, day, repl)
        pl <- plot_clonal_cor(corest, popl, geno, day, repl)
        print(pl)
      }
    }
    dev.off()
  }
}

repl <- c(1); day <- c('day 5')
crpop <- do.call(rbind, lapply(pops, function(popl){
  unlist(lapply(unique(metadata$genotype), function(geno){
    corest <- estimate_clonal_cor(counts, metadata, popl, pop_background, geno, day, repl)
    c(corest$estimate, corest$interval)
  }))
}))
rownames(crpop) <- pops
# plot clonal correlations
df_long <- restructure_table(crpop)
pl <- plot_herit_barplots(df_long)
pl
crpop

# plot all pairs of day and replicate
for (repl in c(1, 2)){
  for (day in c('day 5', 'day 10')){
    crpop <- do.call(rbind, lapply(pops, function(popl){
      unlist(lapply(unique(metadata$genotype), function(geno){
        corest <- estimate_clonal_cor(counts, metadata, popl, pop_background, geno, day, repl)
        c(corest$estimate, corest$interval)
      }))
    }))
    rownames(crpop) <- pops
    #
    df_long <- restructure_table(crop)
    pl <- plot_herit_barplots(df_long)
    #
    path_plot <- paste('figures/heritable_cor_repl_', gsub(" ", "_", day), '_repl', repl,'.pdf', sep = '')
    pdf(path_plot)
    print(pl)
    dev.off()
  }
}
