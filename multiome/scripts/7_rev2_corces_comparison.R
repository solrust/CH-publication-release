library(Seurat)
library(Signac)
library(Matrix)
library(data.table)
library(edgeR)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(pheatmap)
set.seed(123)

dir.create("data/external/corces2016", recursive = TRUE, showWarnings = FALSE)

dir.create("figures/corces_atac_benchmark", recursive = TRUE, showWarnings = FALSE)

dir.create("data/processed/corces_atac_benchmark", recursive = TRUE, showWarnings = FALSE)

### upload multiome object
if (is.null(obj.combined)){
  obj.combined <- readRDS(
    "data/processed/multiome_seurat.rds"
  )
}
DefaultAssay(obj.combined) <- "ATAC"
multiome.fragments <- Fragments(obj.combined[["ATAC"]])
length(multiome.fragments)

### download the Corces count matrix and readable sample header
download_if_missing <- function(url, destination){
  if (!file.exists(destination)) {
    message("Downloading: ", basename(destination))
    status <- download.file(
      url = url,
      destfile = destination,
      mode = "wb",
      quiet = FALSE
    )

    if (!identical(status, 0L)) {
      stop("Download failed: ", url)
    }
  }
}

corces.file <- file.path("data/external/corces2016", "GSE74912_ATACseq_All_Counts.txt.gz")
corces.header.file <- file.path("data/external/corces2016", "Corces_HemeAML_ReHeader.txt")

download_if_missing(
  url = paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/",
    "GSE74912/suppl/GSE74912_ATACseq_All_Counts.txt.gz"
  ),
  destination = corces.file
)

download_if_missing(
  url = paste0(
    "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/",
    "Corces_HemeAML_ReHeader.txt"
  ),
  destination = corces.header.file
)

# Read the original column names without loading the full matrix
original.names <- names(fread(corces.file, nrows = 0, check.names = FALSE))
corces.header <- fread(corces.header.file, data.table = FALSE)
new.names <- corces.header$NewHeader

# Immunophenotypically purified HSPC types to use
reference.types <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP")

all.sample.names <- setdiff(new.names, c("chr", "start", "end"))

reference.samples <- grep(paste0("_(", paste(reference.types, collapse = "|"), ")$"),
  all.sample.names,
  value = TRUE
)
reference.samples
if (length(reference.samples) == 0) stop("No Corces HSPC samples were recognized from the new headers.")

columns.to.read <- c(match(c("chr", "start", "end"), new.names), match(reference.samples, new.names))

corces.data <- fread(
  corces.file,
  select = columns.to.read,
  check.names = FALSE
)

setnames(corces.data, new.names[columns.to.read])

corces.counts.hg19 <- as.matrix(
  corces.data[, reference.samples, with = FALSE]
)
storage.mode(corces.counts.hg19) <- "integer"

corces.metadata <- data.frame(
  sample = reference.samples,
  celltype = sub("^.*_", "", reference.samples),
  donor = sub("_.*$", "", reference.samples),
  stringsAsFactors = FALSE
)
table(corces.metadata$celltype)
corces.metadata

### Convert Corces BED coordinates from hg19 to hg38
corces.gr.hg19 <- GRanges(
  seqnames = corces.data$chr,
  ranges = IRanges(
    start = corces.data$start + 1L,
    end = corces.data$end
  )
)

chain.gz <- file.path("data/external/corces2016", "hg19ToHg38.over.chain.gz")

chain.file <- file.path("data/external/corces2016", "hg19ToHg38.over.chain")

download_if_missing(
  url = paste0(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/",
    "liftOver/hg19ToHg38.over.chain.gz"
  ),
  destination = chain.gz
)

# Base-R gunzip to avoid another package dependency
gunzip_file <- function(source, destination) {

  input <- gzfile(source, open = "rb")
  output <- file(destination, open = "wb")

  on.exit(close(input))
  on.exit(close(output), add = TRUE)

  repeat {

    buffer <- readBin(
      input,
      what = "raw",
      n = 1024 * 1024
    )

    if (length(buffer) == 0) {
      break
    }

    writeBin(buffer, output)
  }
}

if (!file.exists(chain.file)) {
  gunzip_file(chain.gz, chain.file)
}

hg19.to.hg38 <- import.chain(chain.file)

lifted <- liftOver(
  corces.gr.hg19,
  hg19.to.hg38
)

# Retain only unambiguous one-to-one mappings
one.to.one <- elementNROWS(lifted) == 1L

corces.gr.hg38 <- unlist(
  lifted[one.to.one],
  use.names = FALSE
)

corces.counts <- corces.counts.hg19[
  one.to.one,
  ,
  drop = FALSE
]

cat(
  "Original Corces peaks:", length(corces.gr.hg19), "\n",
  "One-to-one hg38 peaks:", length(corces.gr.hg38), "\n"
)


### Filter to autosomes, sensible widths and unique intervals.
#Excluding sex chromosomes avoids sex-associated accessibility entering the cell-identity signatures.
standard.autosomes <- paste0("chr", 1:22)

keep <- (
  as.character(seqnames(corces.gr.hg38)) %in% standard.autosomes &
    width(corces.gr.hg38) >= 100 &
    width(corces.gr.hg38) <= 2000
)

corces.gr.hg38 <- corces.gr.hg38[keep]
corces.counts <- corces.counts[keep, , drop = FALSE]

peak.ids <- paste0(
  seqnames(corces.gr.hg38),
  ":",
  start(corces.gr.hg38),
  "-",
  end(corces.gr.hg38)
)

not.duplicated <- !duplicated(peak.ids)

corces.gr.hg38 <- corces.gr.hg38[not.duplicated]
corces.counts <- corces.counts[not.duplicated, , drop = FALSE]
peak.ids <- peak.ids[not.duplicated]

names(corces.gr.hg38) <- peak.ids
rownames(corces.counts) <- peak.ids

### Optionally restrict to distal regions
use.distal.peaks <- FALSE

if (use.distal.peaks) {

  gene.ranges <- ensembldb::genes(
    EnsDb.Hsapiens.v86,
    return.type = "GRanges"
  )

  seqlevelsStyle(gene.ranges) <- "UCSC"

  promoter.ranges <- promoters(
    gene.ranges,
    upstream = 2000,
    downstream = 500
  )

  is.distal <- countOverlaps(
    corces.gr.hg38,
    promoter.ranges
  ) == 0

  corces.gr.hg38 <- corces.gr.hg38[is.distal]
  corces.counts <- corces.counts[is.distal, , drop = FALSE]
}

cat(
  "Corces hg38 regions after filtering:",
  length(corces.gr.hg38),
  "\n"
)

### normalize Corces samples and derive population-specific peaks
corces.dge <- DGEList(counts = corces.counts)
# Retain peaks detectably accessible in at least two reference samples
corces.cpm.initial <- cpm(corces.dge, log = FALSE)

keep.detected <- rowSums(
  corces.cpm.initial > 1
) >= 2

corces.counts <- corces.counts[
  keep.detected,
  ,
  drop = FALSE
]

corces.gr.hg38 <- corces.gr.hg38[
  keep.detected
]

corces.dge <- DGEList(
  counts = corces.counts
)

corces.dge <- calcNormFactors(
  corces.dge,
  method = "TMM"
)

corces.logcpm <- cpm(
  corces.dge,
  log = TRUE,
  prior.count = 2
)

reference.types.present <- reference.types[
  reference.types %in% corces.metadata$celltype
]

corces.mean.by.type <- vapply(
  reference.types.present,
  function(celltype) {

    rowMeans(
      corces.logcpm[
        ,
        corces.metadata$celltype == celltype,
        drop = FALSE
      ]
    )
  },
  FUN.VALUE = numeric(nrow(corces.logcpm))
)
rownames(corces.mean.by.type) <- rownames(corces.logcpm)

########## Select the top regions per immunophenotypic population:
n.signature.peaks <- 5000

signature.scores <- setNames(
  vector("list", length(reference.types.present)),
  reference.types.present
)

for (celltype in reference.types.present){

  other.types <- setdiff(
    reference.types.present,
    celltype
  )

  closest.other <- apply(
    corces.mean.by.type[
      ,
      other.types,
      drop = FALSE
    ],
    1,
    max
  )

  specificity <- (
    corces.mean.by.type[, celltype] -
      closest.other
  )

  specificity <- sort(
    specificity,
    decreasing = TRUE
  )

  # Require target accessibility to exceed every other reference type
  specificity <- specificity[
    specificity > 0
  ]

  signature.scores[[celltype]] <- head(
    specificity,
    n.signature.peaks
  )
}

sapply(signature.scores, length)


### build the combined reference region set:
signature.peak.ids <- lapply(
  signature.scores,
  names
)

reference.union.ids <- unique(
  unlist(
    signature.peak.ids,
    use.names = FALSE
  )
)

reference.union.gr <- corces.gr.hg38[
  match(
    reference.union.ids,
    names(corces.gr.hg38)
  )
]

names(reference.union.gr) <- reference.union.ids

length(reference.union.gr)

########################## Save the signature definition:
signature.table <- rbindlist(
  lapply(
    names(signature.scores),
    function(celltype) {

      scores <- signature.scores[[celltype]]
      ids <- names(scores)

      ranges <- reference.union.gr[
        match(ids, names(reference.union.gr))
      ]

      data.frame(
        signature = celltype,
        peak_id = ids,
        chromosome = as.character(seqnames(ranges)),
        start = start(ranges),
        end = end(ranges),
        specificity_score = as.numeric(scores),
        stringsAsFactors = FALSE
      )
    }
  )
)

fwrite(
  signature.table,
  file = paste0(
    "data/processed/corces_atac_benchmark/",
    "corces_hspc_hg38_signatures.tsv"
  ),
  sep = "\t"
)

saveRDS(
  list(
    signatures = signature.peak.ids,
    ranges = reference.union.gr,
    mean_reference_profiles = corces.mean.by.type,
    reference_metadata = corces.metadata
  ),
  file = paste0(
    "data/processed/corces_atac_benchmark/",
    "corces_hspc_hg38_signatures.rds"
  )
)

################## Recount your Multiome fragments over the Corces regions
multiome.reference.counts <- FeatureMatrix(
  fragments = multiome.fragments,
  features = reference.union.gr,
  cells = colnames(obj.combined),
  process_n = 2000,
  verbose = TRUE
)
rownames(multiome.reference.counts) <- names(reference.union.gr)
dim(multiome.reference.counts)

################### Construct sample-level ATAC pseudobulks
clust_labels <- c(
  `1` = "HSC",
  `2` = "LMPP-1",
  `3` = "MLP",
  `4` = "Pro-B",
  `5` = "pDC",
  `6` = "LMPP-2",
  `7` = "early GMP",
  `8` = "late GMP",
  `9` = "EMPP/MkP",
  `10` = "EryP"
)
# recreate metadata
multiome.metadata <- obj.combined[[]][colnames(multiome.reference.counts), , drop = FALSE]
raw.population <- as.character(multiome.metadata$projection)
mapped.population <- unname(clust_labels[raw.population])
# Preserve labels when projection already contains biological names
mapped.population[is.na(mapped.population)] <- raw.population[is.na(mapped.population)]
multiome.metadata$population <- mapped.population


### set the grouping variables. Start with one pseudobulk per dataset × population.
sample.column <- "dataset"
genotype.column <- NULL
group.columns <- c(sample.column, "population")

# Activate this when genotype varies within a dataset
if (!is.null(genotype.column) && genotype.column %in% colnames(multiome.metadata)){
  group.columns <- c(sample.column, genotype.column, "population")
}

group.data <- multiome.metadata[, group.columns, drop = FALSE]
complete.cells <- complete.cases(group.data)
multiome.reference.counts <- multiome.reference.counts[, complete.cells, drop = FALSE]
group.data <- group.data[complete.cells, ,drop = FALSE]

pseudobulk.group <- do.call(
  interaction,
  c(
    group.data,
    list(
      drop = TRUE,
      sep = "__"
    )
  )
)
minimum.cells.per.pseudobulk <- 30
group.cell.counts <- table(pseudobulk.group)
valid.groups <- names(group.cell.counts[group.cell.counts >= minimum.cells.per.pseudobulk])
keep.cells <- pseudobulk.group %in% valid.groups
multiome.reference.counts <- multiome.reference.counts[, keep.cells, drop = FALSE]
group.data <- group.data[keep.cells, , drop = FALSE]
pseudobulk.group <- droplevels(pseudobulk.group[keep.cells])

### Aggregate sparse counts:
pseudobulk.design <- sparseMatrix(
  i = seq_along(pseudobulk.group),
  j = as.integer(pseudobulk.group),
  x = 1,
  dims = c(
    length(pseudobulk.group),
    nlevels(pseudobulk.group)
  ),
  dimnames = list(
    NULL,
    levels(pseudobulk.group)
  )
)

multiome.pseudobulk.counts <- multiome.reference.counts %*% pseudobulk.design
first.cell.in.group <- match(levels(pseudobulk.group), as.character(pseudobulk.group))
pseudobulk.metadata <- group.data[first.cell.in.group, , drop = FALSE]
rownames(pseudobulk.metadata) <- levels(pseudobulk.group)
pseudobulk.metadata$n_cells <- as.integer(table(pseudobulk.group)[rownames(pseudobulk.metadata)])
pseudobulk.metadata

### normalize
multiome.dge <- DGEList(
  counts = as.matrix(
    multiome.pseudobulk.counts
  )
)

multiome.dge <- calcNormFactors(
  multiome.dge,
  method = "TMM"
)

multiome.logcpm <- cpm(
  multiome.dge,
  log = TRUE,
  prior.count = 2
)

#### Save cell numbers and pseudobulk metadata:
write.csv(
  pseudobulk.metadata,
  file = paste0(
    "data/processed/corces_atac_benchmark/",
    "multiome_atac_pseudobulk_metadata.csv"
  ),
  row.names = TRUE
)

######################## Calculate Corces signature scores in your populations ############
# Corces mean bulk profiles on the signature-union regions
corces.reference.profiles <- corces.mean.by.type[
  reference.union.ids,
  reference.types.present,
  drop = FALSE
]

# Average Multiome pseudobulk profiles by RNA-derived population
multiome.population.profiles <- vapply(
  population.order,
  function(population) {

    rowMeans(
      multiome.logcpm[
        ,
        pseudobulk.metadata$population == population,
        drop = FALSE
      ]
    )
  },
  FUN.VALUE = numeric(nrow(multiome.logcpm))
)

rownames(multiome.population.profiles) <- rownames(
  multiome.logcpm
)

shared.regions <- intersect(
  rownames(corces.reference.profiles),
  rownames(multiome.population.profiles)
)

# Exclude regions with no observed fragments in the Multiome pseudobulks
shared.regions <- shared.regions[
  Matrix::rowSums(
    multiome.pseudobulk.counts[
      shared.regions,
      ,
      drop = FALSE
    ]
  ) > 0
]

length(shared.regions)

corces.multiome.correlation <- cor(
  x = corces.reference.profiles[
    shared.regions,
    ,
    drop = FALSE
  ],
  y = multiome.population.profiles[
    shared.regions,
    ,
    drop = FALSE
  ],
  method = "pearson",
  use = "pairwise.complete.obs"
)

# Rows after transposition: your Multiome populations
# Columns: purified Corces HSPC populations
correlation.plot <- t(corces.multiome.correlation)
correlation.color.normalized <- apply(
  correlation.plot,
  2,
  function(x) {
    rng <- range(x, na.rm = TRUE)

    if (diff(rng) == 0) {
      rep(0.5, length(x))
    } else {
      (x - rng[1]) / diff(rng)
    }
  }
)
######
#correlation.colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
#correlation.breaks <- seq(-1,  1, length.out = length(correlation.colors) + 1)

pdf(paste0("figures/corces_atac_benchmark/", "corces_multiome_accessibility_correlation.pdf"),
    width = 5, height = 5)

## Restore dimensions and labels explicitly
correlation.color.normalized <- matrix(
  correlation.color.normalized,
  nrow = nrow(correlation.plot),
  ncol = ncol(correlation.plot),
  dimnames = dimnames(correlation.plot)
)
## Show raw correlations as text
display.correlations <- matrix(
  sprintf("%.2f", correlation.plot),
  nrow = nrow(correlation.plot),
  ncol = ncol(correlation.plot),
  dimnames = dimnames(correlation.plot)
)
## mark the best match in every reference column
for (j in seq_len(ncol(correlation.plot))) {

  best.row <- which.max(
    correlation.plot[, j]
  )

  display.correlations[best.row, j] <- paste0(
    display.correlations[best.row, j],
    " *"
  )
}

#plot.new()
ph <- pheatmap::pheatmap(
  correlation.color.normalized,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  #color = relative.colors,
  #breaks = seq(
  #  0,
  #  1,
  #  length.out = length(relative.colors) + 1
  #),
  legend_breaks = c(0, 0.5, 1),
  legend_labels = c(
    "Low",
    "Mid",
    "High"
  ),

  # Numbers remain raw Spearman correlations
  display_numbers = display.correlations,
  number_color = "grey25",

  border_color = NA,
  angle_col = 90,

  main = paste0(
    "Multiome versus purified HSPC ATAC\n",
    length(shared.regions), " reference regions\n"),
  silent = TRUE
)
#ph

library(grid)
library(gridExtra)

ph.with.axis.titles <- gridExtra::arrangeGrob(
  ph$gtable,
  left = grid::textGrob(
    "Multiome populations",
    rot = 90,
    gp = grid::gpar(
      fontsize = 16,
      fontface = "bold"
    )
  ),
  bottom = grid::textGrob(
    "HSPC immunophenotypes",
    gp = grid::gpar(
      fontsize = 16,
      fontface = "bold"
    )
  )
)

grid::grid.newpage()
grid::grid.draw(ph.with.axis.titles)
dev.off()
