library(parallel)
library(Matrix)
library(MASS)

load('../RNA_analysis/data/processed/scrna.RData')
load('../RNA_clone_analysis/data/processed/clone_info.RData')
load('../DNA_analysis/data/processed/DNA_clones.RData')

### specify matrix of clone phenotypes:
clons <- rownames(clone.states2)
summary.size <- aggregate(clone.info[clons,]$size, by = list(clone.info[clons,]$id), sum)
rownames(summary.size) <- summary.size[,1]; summary.size <- summary.size[,-1,drop = F]
clone.pheno <- cbind('growth rate' = rate[clons],
                     'size' = log(clone.info[clons,]$size / summary.size[clone.info[clons,]$id,]),
                     clone.states2[clons, c(1, 4, 5, 8, 10)]
)
rownames(clone.pheno) <- clons
colnames(clone.pheno)[-c(1:2)] <-  c('HSPC', 'Lymphoid', 'DC', 'Myeloid', 'Erythroid') #setdiff(coarse_ann, 'pre-B')
#
pheno_names_new <- c('growth rate', 'clone size', 'HSC bias', 'HSPC bias', 'lymphoid bias', 'dendritic bias', 'pre-B bias', 'myeloid bias', 'erythroid bias')
colnames(clone.pheno) <-  pheno_names_new[c(1,2,3,5,6,8,9)]
head(clone.pheno)

### 1. Match HSC-iM scores to clone identities
shared.cells <- intersect(names(hscII.axis.score), rownames(cell.clone))
length(shared.cells)
hsc.clone.df <- data.frame(cell = shared.cells, cloneid = cell.clone[shared.cells, "cloneid"], HSCiM = hscII.axis.score[shared.cells], stringsAsFactors = FALSE)
head(hsc.clone.df)

### 2. Estimate clone-level HSC-iM score
clone.hscim <- aggregate(HSCiM ~ cloneid, data = hsc.clone.df, FUN = mean, na.rm = TRUE)
clone.nhsc <- table(hsc.clone.df$cloneid)
clone.hscim$n_HSC <- as.integer(clone.nhsc[clone.hscim$cloneid])
rownames(clone.hscim) <- clone.hscim$cloneid
#
head(clone.hscim)
summary(clone.hscim$n_HSC)
table(pmin(clone.hscim$n_HSC, 10))
hist(clone.hscim$n_HSC, breaks = 50, xlab = "Number of HSCs per clone", main = "")

minimum.HSCs <- 5
clone.hscim.use <- clone.hscim[clone.hscim$n_HSC >= minimum.HSCs, , drop = FALSE]
nrow(clone.hscim.use)

### 3. Match clone-level HSC-iM scores and clone phenotypes
shared.clones <- intersect(rownames(clone.hscim.use), rownames(clone.pheno))
length(shared.clones)

clone.analysis <- cbind(clone.hscim.use[shared.clones,c("HSCiM", "n_HSC"), drop = FALSE],
  clone.pheno[shared.clones, , drop = FALSE])
head(clone.analysis)

### 4. Correlation with clonal phenotypes
phenotypes <- colnames(clone.pheno)

hscim.pheno.cor <- do.call(
  rbind,
  lapply(
    phenotypes,
    function(pheno) {

      x <- clone.analysis$HSCiM
      y <- clone.analysis[[pheno]]

      keep <- is.finite(x) & is.finite(y)

      ft <- cor.test(
        x[keep],
        y[keep],
        method = "spearman",
        exact = FALSE
      )

      data.frame(
        phenotype = pheno,
        n_clones = sum(keep),
        rho = unname(ft$estimate),
        pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

hscim.pheno.cor$padj <- p.adjust(
  hscim.pheno.cor$pvalue,
  method = "BH"
)

hscim.pheno.cor <- hscim.pheno.cor[
  order(hscim.pheno.cor$pvalue),
]

hscim.pheno.cor

library(ggplot2)

hscim.pheno.cor$logP <- pmin(
  -log10(hscim.pheno.cor$pvalue),
  10
)

hscim.pheno.cor$program <- "HSC-II"

pl <- ggplot(
  hscim.pheno.cor,
  aes(
    y = phenotype,
    x = program,
    size = logP,
    fill = rho
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  scale_size_continuous(
    range = c(0, 20),
    limits = c(0, 10),
    breaks = c(1, 2, 5, 10)
  ) +
  theme_minimal() +
  labs(
    y = "clonal phenotype",
    x = NULL,
    size = "-log10(P)",
    fill = "Spearman rho"
  ) +
  theme(
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(
      size = 15,
      angle = 30,
      vjust = 1,
      hjust = 1
    )
  )

pl

pdf('figures/HSC-iM_clones.pdf', width = 5, height = 5)
pl
dev.off()
