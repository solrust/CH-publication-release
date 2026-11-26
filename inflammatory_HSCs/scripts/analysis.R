library(ggplot2)

### load hscI and hscII gene sets
hsc_gsets <- readRDS(file = 'Zeng2026_HSCiM/hsc_gene_set.rds')
hscI.genes <- hsc_gsets$hscI.genes
hscII.genes <- hsc_gsets$hscII.genes
sig.genes <- union(hscI.genes, hscII.genes)


### upload RNA data
hsc.env <- new.env()
load("../Association_analysis/data/processed/hsc/hsc_axes.RData", env = hsc.env)

rna.env <- new.env()
load("../RNA_analysis/data/processed/scrna.RData", env = rna.env)

clone.env <- new.env()
load("../Association_analysis/data/processed/default/association_summary.RData", env = clone.env)

### analysis
cls <- rownames(hsc.env$prr$x)
sig.use <- intersect(sig.genes, colnames(rna.env$p2counts))
hsc.sig.expr <- as.matrix(rna.env$p2counts[cls, sig.use, drop = FALSE])
hsc.sig.z <- scale(hsc.sig.expr, center = TRUE, scale = TRUE)

## remove genes with zero variance
keep <- apply(hsc.sig.z, 2, function(x) all(is.finite(x)))
hsc.sig.z <- hsc.sig.z[, keep, drop = FALSE]
dim(hsc.sig.z)

hscI.use <- intersect(hscI.genes, colnames(hsc.sig.z)); hscII.use <- intersect(hscII.genes, colnames(hsc.sig.z))
hsc12.pca <- prcomp(hsc.sig.z, center = FALSE, scale. = FALSE)
summary(hsc12.pca)

### prepare data to train a gene-level classifier
pcs.use <- 1:100
gene.pc <- sweep(
  hsc12.pca$rotation[, pcs.use, drop = FALSE],
  MARGIN = 2,
  STATS = hsc12.pca$sdev[pcs.use],
  FUN = "*"
)
dim(gene.pc)
gene.class <- rep(NA_integer_, nrow(gene.pc))
names(gene.class) <- rownames(gene.pc)
gene.class[rownames(gene.pc) %in% hscI.genes] <- 0
gene.class[ rownames(gene.pc) %in% hscII.genes] <- 1
gene.pc.train <- gene.pc
gene.class.train <- gene.class


### find optimal linear separation of HSC-I and HSC-II gene sets
library(glmnet)
set.seed(123)

cvfit <- cv.glmnet(
  x = gene.pc.train,
  y = gene.class.train,
  family = "binomial",
  alpha = 0, # ridge
  standardize = T,
  nfolds = 10
)
cvfit
plot(cvfit)

fit.coef <- coef(cvfit, s = "lambda.1se")
beta <- as.numeric(fit.coef[-1, 1])
names(beta) <- colnames(gene.pc.train)
intercept <- as.numeric(fit.coef[1, 1])
gene.linear.score <- drop(gene.pc.train %*% beta)

gene.df <- data.frame(
  gene = rownames(gene.pc.train),
  score = gene.linear.score,
  class = factor(
    gene.class.train,
    levels = c(0, 1),
    labels = c("HSC-I", "HSC-II")
  )
)

pl.box <- ggplot(gene.df, aes(x = class, y = score, fill = class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15,size = 1) +
  theme_bw() +
  labs(x = NULL, y = "HSC-I/HSC-II discriminant score") +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))
pl.box

pdf("figures/HSCI_HSCII_gene_discriminant_boxplot.pdf", 5, 5)
print(pl.box)
dev.off()

roc.obj <- pROC::roc(
  response = gene.df$class,
  predictor = gene.df$score,
  levels = c("HSC-I", "HSC-II"),
  direction = "<",
  quiet = TRUE
)
pROC::auc(roc.obj)
pROC::ci.auc(roc.obj)

roc.df <- data.frame(FPR = 1 - roc.obj$specificities, TPR = roc.obj$sensitivities)
auc.val <- as.numeric(pROC::auc(roc.obj))

pl.roc <- ggplot(roc.df, aes(FPR, TPR)) +
  geom_line(linewidth = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() + theme_bw() +
  labs(x = "False positive rate", y = "True positive rate",
       title = sprintf("HSC-I vs HSC-II genes, AUC = %.2f", auc.val)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))
pl.roc

pdf("figures/HSCI_HSCII_gene_discriminant_ROC.pdf", 5, 5)
print(pl.roc)
dev.off()


### get  cell projection of optimal axis
cell.axis <- beta * hsc12.pca$sdev[pcs.use]
names(cell.axis) <- colnames(hsc12.pca$x[, pcs.use, drop = FALSE])
hscII.axis.score <- drop(hsc12.pca$x[, pcs.use, drop = FALSE] %*% cell.axis)
names(hscII.axis.score) <- rownames(hsc12.pca$x)

### Visualize the supervised HSC-I → HSC-II axis on your HSC embedding
scor <- hscII.axis.score[cls]
cols <- colorRampPalette(c("blue", "grey90", "red"))(10)[as.numeric(cut(rank(scor), breaks = 10, include.lowest = TRUE))]

pdf('figures/HSCI_HSCII_emb_score.pdf', width = 6, height = 6)
plot(rna.env$emb[cls, ], pch = 19, cex = 0.5, col = adjustcolor(cols, 0.5), xlab = '', ylab = '',
  xaxt = "n",  yaxt = "n", xlim = c(-2, 5), ylim = c(2, 10), cex.main = 2, font.main = 1,  main = "HSC-II / HSC-I axis")
dev.off()
