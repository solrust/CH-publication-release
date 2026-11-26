# In vitro split-pool lineage tracing

Analysis of lineage-bias consistency between technical splits in the GS175 AAVS1-control and TET2 experiments. The workflow generates Fig. 5L and Supplementary Fig. 8H/I.

## Data

The **DNA clonal tracking** [Zenodo record](https://doi.org/10.5281/zenodo.24545120) contains the required inputs:

```text
clones/data/GS175/GS_175_larry in vitro 3.csv
clones/data/GS175/starcode/     # 112 barcode-count tables
```

The sample sheet identifies genotype, collection day, biological replicate, technical split and cell population. The analysis begins with processed Starcode counts; it does not convert FASTQs to barcode counts.

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("clones/publication/run_publication.R")
result <- clones_pub_run(
  output_root = file.path(getwd(), "clones/publication/reproduction"),
  validate = TRUE
)
```

Choose a new output directory. The correlation intervals use `confintr` 1.0.2, included in the shared package inventory.

## Outputs and interpretation

The output directory contains five numerical tables in `source_data/` and three PDFs in `figures/`. Fig. 5L uses day 5 and both biological replicates; S8H shows early-erythroid clone fractions for day 5, replicate 1; S8I summarizes day 10, replicate 1.

Lineage bias is `target / (target + bulk)` within each technical split. The analysis retains barcodes with at least 20 reads overall and requires at least four target-plus-bulk reads in each split for a comparison. Pearson correlations and 90% confidence intervals describe agreement between splits within a biological replicate, not between independent experiments. The original analysis is in `scripts/1_processing.R`.
