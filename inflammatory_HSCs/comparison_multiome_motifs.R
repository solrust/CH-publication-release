dir.create("Zeng2026_HSCiM", showWarnings = FALSE)

ext5.file <- file.path(
  "Zeng2026_HSCiM",
  "Source_Data_Extended_Data_Fig5.xlsx"
)

download.file(
  url = paste0(
    "https://media.springernature.com/original/",
    "springer-static/esm/art%3A10.1038%2Fs41586-026-10522-7/",
    "MediaObjects/41586_2026_10522_MOESM14_ESM.xlsx"
  ),
  destfile = ext5.file,
  mode = "wb"
)

file.info(ext5.file)$size
library(readxl)
excel_sheets(ext5.file)
sheets <- readxl::excel_sheets(ext5.file)
sheets
grep("5g|Fig.*5g|ED.*5g", sheets, value = TRUE, ignore.case = TRUE)
pub.motif <- as.data.frame(readxl::read_excel(ext5.file, sheet = "ExtFig5g"))
dim(pub.motif)
head(pub.motif)

shared <- intersect(motif.res$motif, pub.motif$Motif_JASPAR)
length(shared)

cmp <- merge(
  motif.res,
  pub.motif,
  by.x = "motif",
  by.y = "Motif_JASPAR"
)
dim(cmp)
cor.test(cmp$beta, cmp$FoldEnrichment_HSCiMvsHSCI, method = "spearman")
qplot(cmp$beta, log(cmp$FoldEnrichment_HSCiMvsHSCI))


library(ggrepel)
cmp$logP <- -log10(pmax(cmp$pval, .Machine$double.xmin))
cmp$sig <- factor(ifelse(cmp$padj < 0.01, "FDR < 0.01",
                         ifelse(cmp$padj < 0.05, "FDR < 0.05", "NS")),
                  levels = c("NS", "FDR < 0.05", "FDR < 0.01"))

cmp$label <- ifelse(cmp$padj < 0.05, cmp$Motif_Name, NA)

pl <- ggplot(cmp, aes(beta, log(FoldEnrichment_HSCiMvsHSCI))) +
  geom_point(aes(color = sig), size = 3) +
  scale_color_manual(values = c("NS" = "grey75",
                                "FDR < 0.05" = "lightcoral",
                                "FDR < 0.01" = "red")) +
  ggrepel::geom_text_repel(aes(label = label), size = 4,
                           na.rm = TRUE, max.overlaps = 25) +
  theme_bw(base_size = 18) +
  labs(x = "Motif association with HSC-iM score",
       y = "log(published HSC-iM motif enrichment)",
       color = NULL)
pl

ct <- cor.test(cmp$beta, cmp$FoldEnrichment_HSCiMvsHSCI, method = "spearman")
lab.cor <- sprintf("Spearman = %.2f\nP = %.2g", unname(ct$estimate), ct$p.value)

pl <- ggplot(cmp, aes(beta, log(FoldEnrichment_HSCiMvsHSCI))) +
  geom_point(aes(color = sig), size = 1.5) +
  scale_color_manual(values = c("NS" = "grey75",
                                "FDR < 0.05" = "lightcoral",
                                "FDR < 0.01" = "red")) +
  ggrepel::geom_text_repel(aes(label = label), size = 4,
                           na.rm = TRUE, max.overlaps = 15) +
  annotate("text", x = -Inf, y = Inf, label = lab.cor,
           hjust = -0.1, vjust = 1.1, size = 5) +
  theme_bw(base_size = 18) +
  labs(x = "Motif association with HSC-iM score",
       y = "log(published HSC-iM motif enrichment)",
       color = NULL)

pl

pdf('figures/multiome_TF_motif_comparison.pdf')
pl
dev.off()
