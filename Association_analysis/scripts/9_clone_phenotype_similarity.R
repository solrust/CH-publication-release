library(ggplot2)
library(ggrastr)
library(tidyr)
library(dplyr)

sz <- clone.info[rownames(clone.pheno),]$size
clone.outlier <-  abs(clone.pheno$`growth rate`) > 1 | sz <= 50
clone.info.filter <- clone.pheno[!clone.outlier,1:7]
clone.info.filter[,3:7] <- log(clone.info.filter[,3:7]+1)
xx <- cor(clone.info.filter, method = 'pearson', use = 'na.or.complete')
#
dist_matrix <- as.dist(1 - xx)
custom_hclust <- hclust(dist_matrix, method = "complete")
diag(xx) <- NA
#
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
max_val <- 0.5#max(abs(xx), na.rm = TRUE)
my_breaks <- seq(-max_val, max_val, length.out = length(my_colors) + 1)
#
pheat <- pheatmap::pheatmap(xx,
                   cluster_rows = custom_hclust, cluster_cols = custom_hclust,
                   color = my_colors, breaks = my_breaks, main = 'Similarity of clone phenotypes')
pdf('figures/default/clone_phenotypes_similarity.pdf', width = 3.5, height = 3.5)
pheat
dev.off()

# 1. Select your target columns
cols_to_plot <- c("growth rate", "clone size", "HSC bias", "erythroid bias", "myeloid bias", "lymphoid bias")
# 2. Create a clean base data frame with an explicit ID column first
tbl_with_id <- clone.pheno[!clone.outlier, cols_to_plot]
pdf('figures/default/clone_phenotypes_pairwise.pdf', width = 4, height = 4)
for (phen1 in cols_to_plot){
  for (phen2 in setdiff(cols_to_plot, phen1)){
    r_val <- cor(tbl_with_id[[phen1]], tbl_with_id[[phen2]], use = "pairwise.complete.obs")
    stat_label <- sprintf("r = %.2f", r_val)
    pl <- ggplot(tbl_with_id, aes(x = .data[[phen2]], y = .data[[phen1]])) +
      geom_point(size = 4, alpha = 0.5, col = "black") +
      geom_smooth(method = "lm", formula = y ~ x, color = "blue", fill = "grey80", alpha = 0.5, se = TRUE) +
      geom_text(
        aes(x = Inf, y = Inf, label = stat_label),
        color = "darkred",
        fontface = "bold",
        size = 6,
        hjust = 1.1,  # Push slightly left from the absolute edge
        vjust = 1.5   # Push slightly down from the absolute top
      ) +
      theme_bw() + labs(x = phen2, y = phen1) +
      theme(text = element_text(size = 20), panel.grid.minor = element_blank())
    print(pl)
  }
}
dev.off()
