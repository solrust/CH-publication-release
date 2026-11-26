# Target-site editing measurements

Target-spanning transcript measurements for TET2 and DNMT3A editing in the LARRY RNA experiments, shown in Supplementary Fig. 5G.

## Data

The saved eight-run count table is supplied at `publication/source_data/FigS5G_target_site_editing_counts.tsv.gz`. Sample-to-genotype metadata are in `RNA_analysis/data/processed/hashtag_meta.rds`, supplied by the **RNA/LARRY** record in the [dataset index](../publication/datasets.tsv).

Recounting alignments requires the original indexed BAMs from `human_larry_exp123` and `human_larry_exp127`. Those files are not included in the processed-data archives and are not needed for the plot-only workflow below.

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("indels_calling/publication/run_publication.R")
run_indels_publication(
  mode = "plot",
  checkpoint_path = "indels_calling/publication/source_data/FigS5G_target_site_editing_counts.tsv.gz",
  output_root = "indels_calling/publication/reproduction",
  control_label = "WT"
)
```

Use a new output directory. The result is one numerical count table and a two-panel PDF under `source_data/` and `figures/`. `mode = "recount"` is available when the indexed BAM inputs have been restored.

## Interpretation

The control metadata use `WT`; the manuscript labels the same control `AAVS1`. The historical numerator counts distinct cell-barcode/UMI pairs from target-spanning alignments containing an insertion or deletion anywhere in the CIGAR string. It does not establish that the indel overlaps the cut site. The denominator counts distinct UMIs among target-spanning alignments; the table also retains a cell-aware denominator. These are editing-proxy molecule counts, not precise population-level editing-efficiency estimates. The original counting code is in `scripts/`.
