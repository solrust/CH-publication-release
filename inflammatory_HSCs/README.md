# Inflammatory HSC signatures

Comparison of inflammatory HSC transcriptional signatures with clonal phenotypes and chromatin accessibility. The workflows generate Supplementary Fig. 10B from scRNA-seq and Fig. 10C from Multiome data.

## Data

Use **RNA/LARRY**, **Multiome** and **DNA clonal tracking** from the [dataset index](../publication/datasets.tsv). The analyses use the shared RNA, HSC-axis, clone-growth and Multiome objects. `Zeng2026_HSCiM/` contains the signature gene set and published motif-reference workbook used for the comparison; retain their attribution to the original study.

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("inflammatory_HSCs/publication/run_publication.R")
result <- run_inflammatory_publication(
  panels = c("S10B", "S10C"),
  output_root = file.path(getwd(), "inflammatory_HSCs/publication/reproduction")
)
```

Choose a new output directory. The analysis loads large RNA and Multiome objects; allow enough memory for the uncompressed objects.

## Outputs and interpretation

The output directory contains ten numerical tables in `source_data/` and eight PDF components in `figures/`. Tables describe gene-signature discrimination, cell-level scores, clone associations and motif comparisons.

S10C uses SCT counts divided by total SCT counts per cell before standardization and the PCA-based discriminant. Substituting `SCT@scale.data` changes the signature coverage and does not reproduce the same analysis. S10B uses the RNA-based HSC representation. The original analysis code is in `scripts/`.
