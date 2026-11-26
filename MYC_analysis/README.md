# MYC-associated clonal phenotypes

MYC expression, MYC-target pathway activity and MYC motif accessibility in relation to clonal phenotypes. The analysis generates Supplementary Fig. 7G–H and 10D.

## Data

Use the **RNA/LARRY** and **Multiome** records in the [dataset index](../publication/datasets.tsv). Required inputs are:

```text
RNA_analysis/data/processed/scrna.RData
Association_analysis/data/processed/hsc/hsc_axes.RData
Association_analysis/data/processed/use/association_summary.RData
pipeline/resource/motif_associations.xlsx
```

The `use/` association results are intentional for this project. The motif workbook is read from inside the repository, not a sibling directory.

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("MYC_analysis/publication/run_publication.R")
result <- run_myc_publication(
  panels = c("S7G", "S7H", "S10D"),
  output_root = file.path(getwd(), "MYC_analysis/publication/reproduction")
)
```

Choose a new output directory. S7G loads the large RNA expression object; S7H and S10D use smaller saved results.

## Outputs

The run writes three numerical tables in `source_data/` and three PDFs in `figures/`. S7G reports HSC expression/program scores and coordinates; S7H reports MYC targets V1/V2 enrichment across clone phenotypes; S10D reports motif-accessibility associations.

S7H dot colour represents the normalized enrichment score. S10D uses the motif association effect size and P value. The original analyses are in `scripts/`.
