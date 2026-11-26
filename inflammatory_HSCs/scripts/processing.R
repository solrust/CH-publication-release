## Download Zeng et al. Nature 2026 HSC-iM resources
resource.dir <- "Zeng2026_HSCiM"

dir.create(
  resource.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

download_if_missing <- function(url, destination) {
  if (!file.exists(destination)) {

    message("Downloading: ", basename(destination))

    download.file(
      url = url,
      destfile = destination,
      mode = "wb",
      method = "libcurl"
    )
  }
}

## 1. Supplementary Tables 1-22
supp.file <- file.path(resource.dir, "Zeng2026_Supplementary_Tables.xlsx")
download_if_missing(url = paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1038%2Fs41586-026-10522-7/",
  "MediaObjects/41586_2026_10522_MOESM3_ESM.xlsx"
),
destination = supp.file
)

library(readxl)
hsc12.de <- read_excel( supp.file, sheet = "ST9_Xenograft_HSCII_vs_HSCI_DE")
hsc12.de <- as.data.frame(hsc12.de)
dim(hsc12.de)
colnames(hsc12.de)
head(hsc12.de, 10)

gene.col <- 'Gene'; lfc.col <- 'DESeq_logFC'; padj.col <- 'DESeq_padj'
de <- hsc12.de[
  !is.na(hsc12.de[[gene.col]]) &
    !is.na(hsc12.de[[lfc.col]]) &
    !is.na(hsc12.de[[padj.col]]) &
    hsc12.de[[padj.col]] < 0.01,
]
dim(de)
n.program <- 200

## HSC-II-high genes
hscII.de <- de[de[[lfc.col]] > 0,]
hscII.de <- hscII.de[order(hscII.de[[lfc.col]], decreasing = TRUE),]
hscII.genes <- head(unique(hscII.de[[gene.col]]), n.program)

hscI.de <- de[de[[lfc.col]] < 0, ]
hscI.de <- hscI.de[order(hscI.de[[lfc.col]], decreasing = FALSE), ]
hscI.genes <- head(unique(hscI.de[[gene.col]]), n.program)

length(hscI.genes)
length(hscII.genes)
head(hscI.genes, 20)
head(hscII.genes, 20)

saveRDS(list(hscI.genes = hscI.genes, hscII.genes = hscII.genes), file = 'Zeng2026_HSCiM/hsc_gene_set.rds')


