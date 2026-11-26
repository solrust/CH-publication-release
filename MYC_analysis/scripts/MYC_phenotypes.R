library(ggplot2)
library(scales)
library(RColorBrewer)

file_path <- '../Association_analysis/data/processed/use/association_summary.RData'
ass.env <- new.env()
load(file_path, envir = ass.env)
env.dna_clones <- new.env()
load('../DNA_analysis/data/processed/DNA_clones.RData', envir = env.dna_clones)

### estimate correlation of MYC expression and MYC targets with clonal behaviors

myc_path_names <- c('HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MYC_TARGETS_V2')#, 'MYC')
myc_in_pheno <- do.call(rbind, lapply(myc_path_names, function(pathn){
  myc_in_pheno <- do.call(rbind, lapply(ass.env$gsea.clone$gsea.clone, function(tbl) tbl[grepl(pathn, tbl$pathway), 1:6]))
  myc_in_pheno <- as.data.frame(myc_in_pheno)
  myc_in_pheno$phenotype <- names(ass.env$test.basic.list)
  myc_in_pheno$pathway <- pathn
  return(myc_in_pheno)
}))
myc_in_pheno$logP <- pmin(-log10(myc_in_pheno$pval), 10)

pl <- ggplot(myc_in_pheno, aes(y = phenotype, x = pathway, size = logP, fill = NES)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  #scale_size_continuous(range = c(0, 20)) +
  scale_size_continuous(range = c(0, 20),
                        limits = c(0, 10),
                        breaks = c(1, 2, 5, 10),
  )+
  theme_minimal() + labs(y = "clonal phenotype", size = "logP", color = "NES") +
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle =30, vjust = 1, hjust = 1))
pl

pdf('figures/HSC_phenotypes_MYC.pdf')#, width = 6, height = 6)
pl
dev.off()


### upload and estimate significance of MYC motif in HSC phenotypes
file_path_motif <- '../resource/motif_associations.xlsx'
sheet_names <- readxl::excel_sheets('../resource/motif_associations.xlsx')
tbl_motifs <- lapply(sheet_names, function(x) {
  y <- readxl::read_excel(file_path_motif, sheet = x)
  y <- as.data.frame(y)
  rownames(y) <- y[,1]
  y <- y[,-1]
})
names(tbl_motifs) <- sheet_names
head(tbl_motifs)

myc_motif_pheno <- do.call(rbind, lapply(sheet_names, function(nm){
  cbind(tbl_motifs[[nm]]['MYC', , drop = F], 'phenotype' = nm)
}))
myc_motif_pheno$motif <- 'MYC motif'
myc_motif_pheno <- myc_motif_pheno[myc_motif_pheno$phenotype != 'mk bias',]

pl1 <- ggplot(myc_motif_pheno, aes(y = phenotype, x = motif, size = -log(pvalue), fill = log2FoldChange)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(0, 20),
                        limits = c(0, 10),
                        breaks = c(1, 2, 5, 10),
  )+
  theme_minimal() + labs(y = "clonal phenotype", size = "-logP", color = "log2FoldChange") +
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle =30, vjust = 1, hjust = 1))
pl1

pdf('figures/HSC_phenotypes_MYC_motif.pdf')#, width = 6, height = 6)
pl1
dev.off()

### examine expression of MYC and MYC targets across HSC subpopulations.
hsc.env <- new.env()
load('../Association_analysis/data/processed/hsc/hsc_axes.RData', env = hsc.env)
rna.env <- new.env()
load('../RNA_analysis/data/processed/scrna.RData', env = rna.env)


cls <- rownames(hsc.env$prr$x)
# MYC gene expression
ft <- lm(rna.env$p2counts[cls, 'MYC'] ~ hsc.env$prr$x[, 1:40])
myc_expr <- ft$fitted.values

# MYC target gene expression
myc_targets <- unlist(ass.env$gsea.clone$gsea.clone[[1]][grepl(myc_path_names[1], ass.env$gsea.clone$gsea.clone[[1]]$pathway),]$leadingEdge)
str(myc_targets)
myc_targets_over <- intersect(rownames(hsc.env$prr$rotation), myc_targets)
str(myc_targets_over)
myc_targets_expr <- hsc.env$prr$x %*%  colMeans(hsc.env$prr$rotation[myc_targets_over, ])
myc_targets_expr <- myc_targets_expr[cls,1]
str(myc_targets_expr)

png('figures/MYC_hsc.png', width = 6 * 150, height = 3 * 150)
scor <- myc_expr
par(mfrow=c(1,2), mar = rep(2, 4))
cols <- colorRampPalette(c("blue","grey90", "red" ))(n = 10)[as.numeric(cut(rank(scor), 10))]
plot(rna.env$emb[cls,], pch = 19, cex = 0.5, col = adjustcolor(cols, 0.5), xaxt = 'n', yaxt = 'n',
     xlim = c(-2, 5), ylim = c(2, 10), cex.main = 2, font.main = 1, main = 'MYC')
#
scor <- myc_targets_expr
cols <- colorRampPalette(c("blue","grey90", "red" ))(n = 10)[as.numeric(cut(rank(scor), 10))]
plot(rna.env$emb[cls,], pch = 19, cex = 0.5, col = adjustcolor(cols, 0.5), xaxt = 'n', yaxt = 'n',
     xlim = c(-2, 5), ylim = c(2, 10), cex.main = 2, font.main = 1, main = 'MYC targets')
dev.off()

