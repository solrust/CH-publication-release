source('R/gsea_clone_summary.R')
source('R/association_routines.R')
source('../RNA_analysis/R/plotEnrichment_new.R')
source('../RNA_analysis/R/draw_go.R')
library(parallel)
library(Matrix)
library(MASS)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(ggtext)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(data.table)
library(plyr)


load('../RNA_analysis/data/processed/scrna.RData')
load('../RNA_analysis/data/processed/mutant_de.RData')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
gene_sets <- readRDS('../RNA_analysis/data/genesets/processed/gene_sets_use.rds')
de.env <- new.env()
load('../RNA_analysis/data/processed/mutant_de.RData', envir = de.env)
fgsea.mut.score <- de.env$fgsea.all.score[,c('DNMT3A 1', 'TET2 1')]
colnames(fgsea.mut.score) <- c('DNMT3A', 'TET2')
rm(de.env)

### choose folder to access and save the results
save_folder <- 'default'
data_folder <- paste('data/processed/', save_folder, '/', sep = '')
figure_folder <- paste('figures/', save_folder, '/', sep = '')
if (!file.exists(data_folder)) dir.create(path_folder)
if (!file.exists(figure_folder)) dir.create(figure_folder)
### load the results of association analysis
load(paste(data_folder, 'association_summary.RData', sep = ''))



### save basic statistics on the number of significant genes
do.adj <- TRUE
padj.pheno <- do.call(cbind, lapply(test.basic.list, function(test.basic){
  if (do.adj){
    p.adjust(test.basic[,'p.pheno'], 'BH')
  }else{
    test.basic[,'p.pheno']
  }
}))
padj.pheno.DNMT3A <- do.call(cbind, lapply(test.inter.list, function(test.inter){
  if (do.adj){
    p.adjust(test.inter[,'p.pheno:DNMT3A'], 'BH')
  }else{
    test.inter[,'p.pheno:DNMT3A']
  }
}))
padj.pheno.TET2 <- do.call(cbind, lapply(test.inter.list, function(test.inter){
  if (do.adj){
    p.adjust(test.inter[,'p.pheno:TET2'], 'BH')
  }else{
    test.inter[,'p.pheno:TET2']
  }
}))
#
padj.cutoff <- 0.5
#padj.cutoff <- 0.001
padj.summary <- cbind(pheno = colSums(padj.pheno <= padj.cutoff, na.rm = T),
                      DNMT3A = colSums(padj.pheno.DNMT3A <= padj.cutoff, na.rm = T),
                      TET2 = colSums(padj.pheno.TET2 <= padj.cutoff, na.rm = T)
)
padj.summary
brks <- exp(seq(0, log(max(padj.summary)+1), length.out = 100))-1
pheat <- pheatmap::pheatmap(t(padj.summary),
                   #main = 'gene associations',
                   legend = FALSE,
                   angle_col =  315, fontsize = 24, treeheight_row  = F, treeheight_col = F,
                   #labels_row = c('DE genes', 'interactions\nclone-DNMT3A', 'interactions\nclone-TET2'),
                   labels_row = c('DE genes', 'DNMT3A:pheno', 'TET2:pheno'),
                   color = adjustcolor(colorRampPalette(c("grey95", "darkorange", "darkred" ))(n = 100), 0.8),
                   breaks = brks, display_numbers = T, fontsize_number = 24, number_format = "%.0f", number_color = 'black',
                   border_color = "grey50", #border_width = 35, #border_lty = "dashed"
)
pdf(paste(figure_folder, 'associations_summary_heatmap.pdf', sep = ''), width = 7, height = 3)
#width = 4, height = 5)
pheat
dev.off()

# 1. Convert your summary table to a data.table and capture the row names
# (If it's already a data.table with a column name, adjust 'keep.rownames' accordingly)
dt_plot <- as.data.table(padj.summary, keep.rownames = "Phenotype")

# 2. Rename 'pheno' to 'WT' to match your preferred plot labels
if ("pheno" %in% colnames(dt_plot)) {
  setnames(dt_plot, "pheno", "WT")
}

# 3. Sort the phenotypes by total count (highest to lowest) to keep the plot clean
dt_plot[, Total := WT + DNMT3A + TET2]
setorder(dt_plot, -Total)
dt_plot[, Phenotype := factor(Phenotype, levels = Phenotype)]

# 4. Melt the wide table into a long format required by ggplot2
df_long <- melt(
  dt_plot[, .(Phenotype, WT, DNMT3A, TET2)],
  id.vars       = "Phenotype",
  variable.name = "Genotype",
  value.name    = "Count"
)

# 5. Lock in the genotype factor order so it reads left-to-right consistently
df_long[, Genotype := factor(Genotype, levels = c("WT", "DNMT3A", "TET2"))]

# 6. Build the Grouped Bar Plot
pl <- ggplot(df_long, aes(x = Count, y = Phenotype, fill = Genotype)) +
  # position_dodge() groups the bars together; adjusting 'width' changes the spacing
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    color    = "black",
    width    = 0.75
  ) +
  theme_bw() +
  scale_fill_manual(values = c("WT" = "#999999", "DNMT3A" = "#E69F00", "TET2" = "#56B4E9"),
                    labels = c("WT" = "phenotype", "DNMT3A" = "pheno:DNMT3A", "TET2" = "pheno:TET2")) +
  labs(
    y     = "Phenotype",
    x     = "Significant genes (FDR < 0.5)",
    fill  = "",
  ) +
  theme(
    text               = element_text(size = 20),
    axis.text.x        = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title         = element_text(face = "bold", hjust = 0.5),
    legend.position    = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

pdf(paste(figure_folder, 'association_barplot_summary.pdf', sep = ''), width = 7, height = 4)
pl
dev.off()


pcut <- 1e-3
for(cl in 1:length(test.basic.list)){
  print(cl)
  #cl <- 3
  res <- as.data.frame(test.basic.list[[cl]])
  res$TF <- NA#as.factor(rownames(res) %in% tfs & res$p.pheno < pcut)
  lbl <- rownames(res)
  #lbl <- unlist(lapply(lbl, function(x) paste(substr(x,1,1), substr(tolower(x),2,nchar(x)), sep = '') ))
  lbl[res$p.pheno > pcut] <- ''
  names(lbl) <- rownames(res)
  pl1 <- ggplot(res, aes(pheno / log(2), p.pheno, label = lbl,
                         #col = TF,
                         alpha = I(ifelse(p.pheno > pcut, 0.1, 1)))) +
    geom_point(size = I(2), aes(col = res$TF)) + geom_vline(xintercept = 0, linetype = 2) +
    geom_text_repel(box.padding = 0.5, size = 6,  max.overlaps = 15, fontface = 'italic') +
    labs(x = paste('fold change, log2'), y = 'p-value', title = names(test.basic.list)[cl]) +
    scale_y_log10(limits = c(min(res$p.pheno), 1)) + theme_bw() +
    scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red')) +
    theme(text = element_text(size = 20), axis.title = element_text(size = 30), axis.text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5), legend.position =  'none',
          plot.margin = margin(rep(2, 4), "pt"))
  pl1


  mainPathways <- mainPathways.clones[[cl]]
  fgseaRes <- gsea.clone$gsea.clone[[cl]]

  ind <- fgseaRes$pathway %in% mainPathways
  # a different pathway selection for HSC phenotype
  if (cl == 3){
    ind <- rank(fgseaRes$padj) < 10 & nchar(fgseaRes$pathway) < 60
  }else{
    #ind <- fgseaRes$pathway %in% mainPathways & rank(fgseaRes$padj) <= 15
    ind[ rank(fgseaRes[ind,]$padj) >= 15 ] <- FALSE
  }
  if (sum(ind) > 0){
    mm <- max(abs(-log10(fgseaRes$pval) * sign(fgseaRes$NES)))
    pl2 <- ggplot(fgseaRes[ind,],
                  aes(x = factor(pathway, level = fgseaRes[ind,][order(-sign(NES) * log(pval))]$pathway),
                      y = -log10(pval) * sign(NES), fill = as.factor(sign(NES)) )
    ) +
      ylab('p-value, log10') + geom_bar(position = position_dodge(), stat = "identity", col = adjustcolor('black', 1)) +
      scale_fill_manual(values = c('-1' = adjustcolor('blue', 0.2), '1' = adjustcolor('red', 0.2))) +
      ylim(-mm, mm) +
      coord_flip() + theme_bw() +
      theme(text = element_text(size = 20),
            axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 20), legend.position = 'none',#)+
            axis.title = element_text(size = 25),
            #axis.title.x = element_text(margin = margin(t = 20, b = 0, unit = "pt")),
            #title = element_text(margin = margin(t = 50, b = 50, unit = "pt")),
            #axis.title.y = element_text(margin = margin(r = 50, l = 50, unit = "pt")),
            #plot.margin = margin(0, 0, 0, 0, "pt")
            plot.margin = margin(30, rep(2,3), "pt")
      ) + scale_x_discrete(name = "", position = "top")
    pl2
  }else{
    pl2 <- NULL
  }

  pdf(paste(figure_folder, 'scatterplot_new_', names(test.basic.list)[cl],'.pdf', sep = ''),
      width = 12*(2.5/2), height = 8)
  grid.arrange(arrangeGrob(pl1), arrangeGrob(pl2), widths = c(0.9,  1.6))#, top = paste('cluster', cl))
  dev.off()
}

### show some GO enrichment plots
# plot MYC depletion
cl <- 1
gsc <- gsea.clone$gscore.clone[,cl]
fgseaRes <- gsea.clone$gsea.clone[[cl]]
path <- 'M5926: HALLMARK_MYC_TARGETS_V1'; nm <- 'MYC targets'
pl <- plotEnrichment_new(gene_sets[[path]], gsc[!is.na(gsc)], ticksSize = 0.2, spreadES = 0.4, lsize = 2) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) + labs(title = nm) #+ geom_line(aes(size = I(3)))
pl
#
fl <- paste(figure_folder, 'enrichmnt_growth_myc.pdf', sep = '')
pdf(fl, width = 4, height = 4)
pl
dev.off()


### choose overdispersed/significant GO pathways for analysis of clone/mutant features
do.go <- T
if (do.go){
  mat <- matcl <- gsea.clone$gsea.score
  matmut <- fgsea.mut.score
  feature.var <-  rownames(mat)[apply(mat, 1, sd) > 7]
  feature.var <- intersect(feature.var, rownames(matmut))
  #feature.var <-  rownames(mat)[apply(abs(mat), 1, max) > 10]
  #feature.var <- feature.var0 <- intersect(feature.var, mainPathways)
  feature.var0 <- feature.var #<- old.env$feature.var
  str(feature.var)
  n.cl <- 4
}else{
  mat <- matcl <- gsea.clone$gscore.clone
  matmut <- gscore.mut
  feature.var <-  rownames(mat)[apply(mat, 1, sd) > 4]
  n.cl <- 2
}
head(sort(apply(mat, 1, sd), decreasing = T))
str(feature.var)
#feature.var0 <- feature.var <- feature.var[setype[feature.var] == 'GO']
table(grepl('GO', feature.var))

### cluster GO terms with similar gene sets
str(feature.var0)
gsl <- gene_sets[feature.var0]
gs.dst <- do.call(rbind, lapply(gsl, function(gs1){
  unlist(lapply(gsl, function(gs2) length(intersect(gs1, gs2)) / min(length(gs1), length(gs2)) ))
}))
gs.dst1 <- 1 - gs.dst
diag(gs.dst1) <- 0
#
ig <- igraph::graph_from_adjacency_matrix( gs.dst > 0.15, weighted=TRUE, mode="undirected")
lecl <- igraph::cluster_leiden(ig, objective_function = 'modularity',
                               resolution_parameter = 1.5)#1.5
clsets <- as.factor(lecl$membership)
table(clsets)
clsets[clsets == 4] <- 3
clsets <- factor(clsets, levels = c(1,2,3,5), labels = 1:4)
names(clsets) <- rownames(gs.dst)
#
clsets.dat <- as.data.frame(clsets); colnames(clsets.dat) <- 'GO main'
cols.goseq <- brewer.pal(length(unique(clsets)), "Set3")
names(cols.goseq) <- unique(clsets)#1:length(unique(clsets))
pheat <- pheatmap::pheatmap( gs.dst[order(clsets), order(clsets)],
                    #gs.dst[hcl$order, hcl$order],
                    main = 'Intersection of GO terms', show_colnames  = F, show_rownames  = F,
                    cluster_rows = F, cluster_cols  = F, cellwidth = 1.5, cellheight  = 1.5,
                    treeheight_row = 0, treeheight_col = 0,
                    border_color = NA,#adjustcolor('black', 0.01),
                    annotation_row = clsets.dat, annotation_col = clsets.dat,
                    annotation_colors = list(#'GO terms' = cols.go,
                                             'GO main' = cols.goseq),
                    annotation_legend = F,
)
#
fl <- paste(figure_folder, 'pathways_merge.pdf', sep = '')
pdf(fl, width = 5, height = 5)
pheat
dev.off()


### visualize representative of clustered GO terms with similar gene sets
fl <- paste(figure_folder, 'pathways_merge_representative.pdf', sep = '')
pdf(fl, width = 12, height = 14)
pl <- draw_go(matcl, feature.var0, as.data.frame(clsets), cols.goseq,
              n.terms = 10, only.go = FALSE, term.length = 70, cut.pval = 40)
pl
dev.off()

### annotate names of clustered GO terms
#go.main.names <- c('lymphocyte activation', 'translation', 'cell cycle', 'metabolism', '')
go.main.names <- c('translation', 'immune activation', 'metabolism', 'cell cycle')
clsets.dat1 <- clsets.dat
clsets.dat1[,1] <- go.main.names[clsets.dat1[,1]]
#clsets.dat1 <- clsets.dat1[!is.na(clsets.dat1[,1]),,drop = F]
cols.goseq1 <- cols.goseq
names(cols.goseq1) <- go.main.names[as.numeric(names(cols.goseq))]

### show correlation of GO terms across clones
dst <- as.dist(1-cor(t(mat[feature.var,])))
hcl <- hclust(dst)
dst1 <- as.dist(1-cor((mat[feature.var,])))
#dst1 <- dist(t(mat[feature.var,]))
hcl1 <- hclust(dst1)
setclust <- as.factor(cutree(hcl, n.cl))
setclust <- as.data.frame(setclust)
colnames(setclust) <- ifelse(do.go, 'GO terms', 'genes')
if (do.go){
  setclust.go <- setclust
}else{
  setclust.gene <- setclust
}
cols.go <- rainbow(length(unique(setclust[,1])))
cols.go <- brewer.pal(length(unique(setclust[,1])), "Set2")
names(cols.go) <- unique(setclust[,1])

pheat <- pheatmap::pheatmap(t(cbind(matcl[feature.var, rev(hcl1$order)] ,
                                  matmut[feature.var,])
                              ),
                            show_colnames  = F,
                          #annotation_col = setclust,
                          annotation_col = cbind(clsets.dat1, setclust),
                          main = ifelse(do.go, 'associated GO terms','associated genes'),
                          cluster_rows = F, treeheight_col = 10,
                          #annotation_legend = F,
                          gaps_row = 7,
                          fontsize = 18,
                          breaks = seq(-25, 25, length.out = 100),
                          clustering_distance_cols = dst,
                          annotation_colors = list('GO terms' = cols.go, 'GO main' = cols.goseq1, 'genes' = cols.go),
                          box = TRUE, border_color = "red", border_width = 2, border_lty = "dashed"
)
#
fl <- paste(figure_folder, 'sets_heatmap.pdf', sep = '')
pdf(fl, width = 8, height = 5)
pheat
dev.off()

### plot scatterplot of growth vs HSC-bias
colnames(gsea.clone$gsea.score) <- colnames(gsea.clone$gscore.clone)
tbl <- data.frame(gsea.clone$gsea.score)
tbl$GOterms <- NA
tbl[rownames(clsets.dat1),]$GOterms <- clsets.dat1$`GO main`
tbl <- tbl[order(is.na(tbl$GOterms), decreasing = TRUE), ]
pl <- ggplot(tbl, aes(x = growth.rate, y = HSC.bias, col = GOterms)) +
  geom_point(size = 3, alpha = I(0.75)) + theme_bw() +
  #geom_point_raster(size = 2, alpha = I(0.3), dpi = 300) +
  #rasterise(geom_point(size = 2, alpha = I(0.3)), dpi = 300) +
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) +
  xlim(-30, 30) + ylim(-13, 13) + labs(x = 'growing clones, log pathway significance',
                                       y = 'HSC-biased clones, log pathway significance') +
  theme(text = element_text(size = 15),
        legend.position = c(0.99, 0.99),
        legend.title = element_blank(),
        # Anchors the legend by its top-right corner to prevent clipping
        legend.justification = c("right", "top"),
        # Optional: Adds a clean white background with a subtle border to make it readable
        legend.background = element_rect(fill = "white", color = "gray80", size = 0.2)
  )
pdf('figures/default/scatterplot_growth_vs_output.pdf', width = 5, height = 5)
plot(pl)
dev.off()

### plot representative cluster-specific GO terms
if (do.go){
  pl <- draw_go(matcl, feature.var0, setclust.go, cols.go, n.terms = 15, only.go = FALSE, term.length = 70, cut.pval = 40)
  fl <- paste(figure_folder, 'clone_pathways_representatives.pdf', sep = '')
  pdf(fl,  width = 12, height = 14)
  pl
  dev.off()
}


### upload and process AML-related gene sets
aml_sets <- readRDS("../RNA_analysis/data/genesets/processed/aml_sets.rds")
aml_sign <- aml_sets$aml_sign
kentsis_sign <- aml_sets$kentsis_sign
#
hsc_output <- readxl::read_excel("../RNA_analysis/data/genesets/processed/HSC_output_RF.xlsx")
hsc_output <- as.data.frame(hsc_output)
colnames(hsc_output) <- hsc_output[1,]
hsc_output <- hsc_output[-1,]
hsc_output[,2] <- as.numeric(hsc_output[,2])
hsc_output[,3] <- as.numeric(hsc_output[,3])
hsc_output[,4] <- as.numeric(hsc_output[,4])
head(hsc_output)
#
hsc_output_sign <- hsc_output$gene_names[abs(hsc_output$logfoldchange) > 2 &
                                           hsc_output$`adj.p-value` < 1e-3]
str(hsc_output_sign)
qplot(data = hsc_output, x = logfoldchange, y = `adj.p-value`, log = 'y', col = hsc_output$gene_names %in% hsc_output_sign)

#
gene_list_hsc <- list('AML' = list('up' = rownames(aml_sign)[aml_sign$log2FC > 0 & aml_sign$logP.sum > -log(1e-5) & aml_sign$n.patients.sign >= 2],
                                   'down' = rownames(aml_sign)[aml_sign$log2FC < 0 & aml_sign$logP.sum > -log(1e-5) & aml_sign$n.patients.sign >= 2]),
                      'dormancy' = list('dormant' = rownames(aml_sign)[aml_sign$Zhang == 'dHSC'],
                                        'active' = rownames(aml_sign)[aml_sign$Zhang == 'aHSC']),
                      'label-retention' = list('positive' = names(kentsis_sign)[kentsis_sign == 2],
                                               'negative' = names(kentsis_sign)[kentsis_sign == -2]),
                      'HSC output' = list(#high = rownames(aml_sets$aml_sign)[aml_sets$aml_sign$RF == 'High'],
                                          #low = rownames(aml_sets$aml_sign)[aml_sets$aml_sign$RF == 'dHSC'],
                                            'high' = toupper(hsc_output_sign[hsc_output[hsc_output$gene_names%in%hsc_output_sign,]$logfoldchange < 0]),
                                            'low' = toupper(hsc_output_sign[hsc_output[hsc_output$gene_names%in%hsc_output_sign,]$logfoldchange > 0])
                                          )

)
str(gene_list_hsc)
#
de.mut2wt <- de.mut2wt.shrink
mut.cut.padj <- 3e-1; mut.cut.lfc <- 0.0
#mut.cut.padj <- 1e-2; mut.cut.lfc <- 0.3
gene_list_hsc <- c(list('DNMT3A' = list('up' = rownames(de.mut2wt[[1]]$DNMT3A)[!is.na(de.mut2wt[[1]]$DNMT3A$padj) & de.mut2wt[[1]]$DNMT3A$log2FoldChange >= mut.cut.lfc & de.mut2wt[[1]]$DNMT3A$padj < mut.cut.padj],
                                        'down' = rownames(de.mut2wt[[1]]$DNMT3A)[!is.na(de.mut2wt[[1]]$DNMT3A$padj) & de.mut2wt[[1]]$DNMT3A$log2FoldChange <= -mut.cut.lfc & de.mut2wt[[1]]$DNMT3A$padj < mut.cut.padj]),
                        'TET2' = list('up' = rownames(de.mut2wt[[1]]$TET2)[!is.na(de.mut2wt[[1]]$TET2$padj) & de.mut2wt[[1]]$TET2$log2FoldChange >= mut.cut.lfc & de.mut2wt[[1]]$TET2$padj < mut.cut.padj],
                                      'down' = rownames(de.mut2wt[[1]]$TET2)[!is.na(de.mut2wt[[1]]$TET2$padj) & de.mut2wt[[1]]$TET2$log2FoldChange <= -mut.cut.lfc & de.mut2wt[[1]]$TET2$padj < mut.cut.padj])),
                   gene_list_hsc
)
str(gene_list_hsc)


### phenotype-specific signatures
gscore.mat <- data.frame(gsea.clone$gscore.clone)
colnames(gscore.mat) <- colnames(gsea.clone$gscore.clone)
gscore.mat$gene <- rownames(gscore.mat)
gscore.mat <- melt(as.data.table(gscore.mat), id.vars = c('gene'))#,id.vars = c('growth.rate', 'clone.size'))
colnames(gscore.mat) <- c('gene', 'phenotype', 'score')
head(gscore.mat)
#
for (i in 1:length(gene_list_hsc)){
  vec <- rep(NA, nrow(gscore.mat))
  vec[gscore.mat$gene %in% gene_list_hsc[[i]][[1]]] <- names(gene_list_hsc[[i]])[1]
  vec[gscore.mat$gene %in% gene_list_hsc[[i]][[2]]] <- names(gene_list_hsc[[i]])[2]
  gscore.mat <- cbind(gscore.mat, vec)
  colnames(gscore.mat)[ncol(gscore.mat)] <- names(gene_list_hsc)[i]
}
str(gscore.mat)
#
# TODO:
# .
# .
for (show.mutants in c(T, F)){
  if (show.mutants){
    tmp=melt(gscore.mat[,c(1:3, 4:6)], id.vars = c('gene','score', 'phenotype'))
  }else{
    tmp=melt(gscore.mat[,c(1:3, 7:9)], id.vars = c('gene','score', 'phenotype'))
  }
  tmp1=unique(subset(tmp,!is.na(value),select=c(value,variable)))
  tmp1=ddply(tmp1,.(variable), mutate, k=0:1)
  #
  #clone_pheno <- 'HSC bias' # 'HSC bias' growth rate'
  #clone_pheno <- 'growth rate' # 'HSC bias' growth rate'
  clone_pheno <- c('growth rate', 'erythroid bias')
  pl <- ggplot(subset(tmp[tmp$phenotype %in% clone_pheno, drop = F],!is.na(value)),
           #subset(tmp, phenotype %in% clone_pheno & !is.na(value) & variable == var),
               aes(x = value,
                   #x = phenotype,
                   #x = genotype,
                   y = score, color = value)) +
    ylab('growth gene score') + theme_bw() +
    #ylim(-ymm, ymm) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme(
      #axis.title.x = element_blank(),
      #axis.text.x = element_blank(),
      #text = element_text(size = 20),
      legend.title  = element_blank(), legend.position = "none")+#, legend.position = "bottom") + #, legend.position = "none") +
    geom_boxplot(width = 0.8, outlier.shape=NA, alpha = 1) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1, alpha = 0.5) +
    geom_text(data=tmp1,aes(1.1,-6+k*1.5,label=value),hjust=0, size = 6)+
    geom_point(data=tmp1,aes(1,-6+k*1.5))+
    facet_grid(~variable + phenotype) #+ coord_flip()
  print(pl)
  #
  if (show.mutants){
    fl <- paste(figure_folder, 'HSC_growth_mutants.pdf', sep = '')
  }else{ fl <- paste(figure_folder, 'HSC_growth_aml.pdf', sep = '')}

  pdf(fl, width = 6, height = 3)
  print(pl)
  dev.off()
}

for (show.mutants in c(T, F)){
  show.mutants <- FALSE
  if (show.mutants){
    tmp=melt(gscore.mat[,c(1:3, 4:6)], id.vars = c('gene','score', 'phenotype'))
  }else{
    tmp=melt(gscore.mat[,c(1:3, 7:9)], id.vars = c('gene','score', 'phenotype'))
  }
  tmp1=unique(subset(tmp,!is.na(value),select=c(value,variable)))
  tmp1=ddply(tmp1,.(variable), mutate, k=0:1)
  #
  clone_pheno <- 'HSC bias' # 'HSC bias' growth rate'
  clone_pheno <- c('growth rate')
  alims <- range(tmp[tmp$phenotype == clone_pheno,]$score)
  qlims <- quantile(tmp[tmp$phenotype == clone_pheno,]$score, probs = seq(0, 1, length.out = 100))
  alims <- c(qlims[2], qlims[length(qlims) - 1])
  pls <- lapply(unique(tmp$variable), function(var){
    tmp_subset <- subset(tmp, phenotype %in% clone_pheno & !is.na(value) & variable == var)
    pval <- wilcox.test(score ~ value, data = tmp_subset)$p.value
    p_label <- paste0("p < ", signif(pval, 1))
    pl <- ggplot(tmp_subset,
                 aes(x = value, y = score, color = value)) +
      theme_bw() + labs(title = var) + ylim(alims) +
      scale_x_discrete(guide = guide_axis(angle = 0)) +
      theme(
        plot.title.position = "panel",
        plot.title = element_text(
          hjust = 0,
          #margin = margin(
          #  t = 8,
          #  r = 0,
          #  b = -35,   # negative pushes title into plot area
          #  l = 5
          #),
          size = 25
        ),
        axis.title = element_blank(),
        text = element_text(size = 25),
        legend.position = "none")+#, legend.position = "bottom") + #, legend.position = "none") +
      geom_boxplot(width = 0.8, outlier.shape=NA, alpha = 1) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_point(aes(fill = value), position=position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8),
                 #size = 1, alpha = 0.5
                 shape = 21, color = 'black',
                 size = 4, alpha = 0.5,
                 )
    pl <- pl +   annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = p_label,
      hjust = 1.2,
      vjust = -0.5,
      size = 5
    )
    pl
      #+
  })
  #
  require(cowplot)
  # keep y-axis only on first plot
  pls[[2]] <- pls[[2]] + theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
  pls[[3]] <- pls[[3]] + theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
  aligned <- align_plots(pls[[1]], pls[[2]], pls[[3]], align = "hv", axis = "l")
  pdf(paste(figure_folder, clone_pheno, '_gene_classes.pdf', sep = ''), width = 2.5*4, height = 4)
  grid.arrange(grobs = aligned, nrow = 1,
               left = textGrob(paste("Clone ", clone_pheno,"\n HSC gene Z-score", sep = ''),
                               rot = 90, y = 0.5,  gp = gpar(fontsize = 25)),
               bottom = textGrob("genes", gp = gpar(fontsize = 25)))
  dev.off()
}

### show AUCs
vars <- c( 'dormancy', 'HSC output', 'label-retention')
clone_features_show <- c('growth rate', 'HSC bias', 'clone size')
tmp1 <- tmp[tmp$variable == var & tmp$phenotype %in% clone_features_show,]
# function to prepare data structure for roc curve visualization
prepare_auc_params <- function(tmp1){
  require(dplyr)
  df <- tmp1 %>%
    filter(!is.na(value)) %>%
    mutate(value = factor(value))
  table(df$value, df$phenotype)
  require(pROC)
  roc_list <- df %>%
    group_by(phenotype) %>%
    group_map(~{
      roc_obj <- roc(.x$value, .x$score, quiet = TRUE)
      list(
        phenotype = .y$phenotype,   # ✅ correct way
        roc = roc_obj,
        auc = auc(roc_obj)
      )
    })
  auc_df <- do.call(rbind, lapply(roc_list, function(x) {
    data.frame(phenotype = x$phenotype, AUC = as.numeric(x$auc))
  }))
  #plot
  roc_plot_df <- do.call(rbind, lapply(roc_list, function(x) {
    data.frame(
      specificity = x$roc$specificities,
      sensitivity = x$roc$sensitivities,
      phenotype = x$phenotype
    )
  }))
  #
  auc_labels <- paste0(auc_df$phenotype, " (AUC=", round(auc_df$AUC, 2), ")")
  return(list(roc_plot_df = roc_plot_df, auc_labels = auc_labels))
}


pls <- lapply(vars, function(var){
  tmp1 <- tmp[tmp$variable == var & tmp$phenotype %in% clone_features_show,]
  roc_stat <- prepare_auc_params(tmp1)
  pl <- ggplot(roc_stat$roc_plot_df, aes(x = 1 - specificity, y = sensitivity, color = phenotype)) +
    #geom_point(size = 0.5) +
    geom_path(size = 3) +
    geom_abline(linetype = "dashed") +
    theme_bw() +
    scale_color_discrete(labels = roc_stat$auc_labels) +
    labs(title = paste(var,'signature'), color = "HSC gene sets") +
    theme(text = element_text(size = 30),
          legend.position = c(0.95, 0.05), legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  return(pl)
})
pdf('figures/default/roc_curve_features.pdf', width = 21, height = 7)
grid.arrange(pls[[1]], pls[[2]], pls[[3]], nrow = 1)#, top = paste('cluster', cl))
dev.off()


roc_stat <- prepare_auc_params(tmp1)
pl <- ggplot(roc_stat$roc_plot_df, aes(x = 1 - specificity, y = sensitivity, color = phenotype)) +
  #geom_point(size = 0.5) +
  geom_path(size = 3) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  scale_color_discrete(labels = roc_stat$auc_labels) +
  labs(title = var, color = "HSC gene sets") +
  theme(text = element_text(size = 30),
        legend.position = c(0.95, 0.05), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
pl

### draw p-values of associations
gscore.mat$DNMT3A <- factor(gscore.mat$DNMT3A, levels = c('up', 'down'))
pstat <- do.call(cbind, lapply(c(4:6),
                               #c(7:8),
                               function(i){
                                 print(i)
                                 gscore.mat1 <- subset(gscore.mat, !is.na(gscore.mat[,i]))

                                 #if (do.mutant){
                                 #  vals <- gscore.mat1$genotype
                                 #}else{
                                   vals <- gscore.mat1$phenotype
                                 #}
                                 pvs <- unlist(lapply(unique(vals), function(val){
                                   print(val)
                                   gscore.mat1d <- gscore.mat1[vals == val,]
                                   x1 <- gscore.mat1d$score[gscore.mat1d[,i] == unique(gscore.mat1d[,i])[1]]
                                   x2 <- gscore.mat1d$score[gscore.mat1d[,i] == unique(gscore.mat1d[,i])[2]]
                                   #
                                   #x1 <- gscore.mat1d$score[gscore.mat1d[,i] == 'up']
                                   #x2 <- gscore.mat1d$score[gscore.mat1d[,i] == 'down']
                                   #
                                   wtest <- wilcox.test(x1, x2)
                                   wtest <- t.test(x1, x2)
                                   wtest$p.value * sign(mean(x1) - mean(x2))
                                 }))
                                 names(pvs) <- unique(vals)
                                 return(pvs)
                               }))
#colnames(pstat) <- colnames(gscore.mat)[4:ncol(gscore.mat)]
colnames(pstat) <- colnames(gscore.mat)[4:6]
pstat.sign <- matrix('', ncol = 3, nrow = nrow(pstat))
pstat.sign[abs(pstat) < 5e-2] <- '*'
pstat.sign[abs(pstat) < 1e-3] <- '**'
pstat.sign[abs(pstat) < 5e-4] <- '***'
xx <- -log10(abs(pstat[,1:3])) * sign(pstat[,1:3])
pheat <- pheatmap::pheatmap(xx,
                            border_color = 'black', gaps_col = rep(2, 3),
                            color = colorRampPalette(c('darkblue','grey90', 'darkred'))(100),
                            breaks = seq(-max(abs(xx)), max(abs(xx)), length.out = 100),
                            cluster_rows = F, cluster_cols = F, fontsize = 20,
                            #legend = F,
                            legend_breaks = c(-4, -2, 0, 2, 4),
                            legend_labels = c('0.001', '0.01', '1', '0.01', '0.001'),
                            display_numbers = pstat.sign,
                            angle_col = 0, fontsize_number = 50,
)
#
fl <- paste(figure_folder, 'mutant_effects_clones.pdf', sep = '')
pdf(fl, width = 6, height = 5)
pheat
dev.off()


