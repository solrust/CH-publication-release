# Transcriptional associations with clonal phenotypes

Associations between gene expression, pathway activity and clonal growth, size and lineage bias. The analyses also examine HSC transcriptional programs and their relationship to clonal behaviour.

## Data

Download **RNA/LARRY** and, for analyses using growth rates, **DNA clonal tracking** from the [dataset index](../publication/datasets.tsv).

The main figure analyses use `data/processed/default/association_summary.RData`. The separate `use/` results support downstream MYC and TF analyses and should not be substituted for `default/`. HSC projections are stored in `data/processed/hsc/hsc_axes.RData`.

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("Association_analysis/publication/run_publication.R")
run_association_publication(
  stages = "5C",
  dry_run = FALSE,
  work_root = file.path(getwd(), "reproduction_associations_5C")
)
```

Choose a new `work_root`. Omit `dry_run = FALSE` to list the planned steps without running them.

| Stage | Figures or analysis |
| --- | --- |
| `5C` | Gene-association counts |
| `4FH_S7F` | Gene and pathway associations |
| `4EG` | LSP1 and NPM1 expression versus clone phenotypes |
| `4I` | Comparison of pathway associations |
| `4JK` | HSC gene signatures and ROC curves |
| `4L` | Pathway-association heatmap |
| `4M` | HSC gene-program projections |
| `5B` | Mutation and AML signature associations |
| `S8ABC` | Clone-phenotype correlations and HSC projections |
| `5J_plot` | Fig. 5J from the saved correlation summary |

## Outputs and interpretation

Tables and PDFs are written under `<work_root>/Association_analysis/publication/source_data/` and `figures/`.

`publication/resources/Fig4H_accepted_pathway_ids.tsv` specifies the 16 pathways and their displayed order in Fig. 4H. The runner matches this selection to the saved enrichment results before plotting. Its audit and curation scripts are part of that calculation, not optional development files.

Fig. 4I/L pathway scores are recovered, but the exact functional-category assignments in the assembled artwork have not been recovered; reconstructed labels remain provisional. Fig. 5J normally uses saved means and standard deviations. Its optional `5J_refit` stage requires `allow_refit = TRUE` and repeats cell-level training/test splits, not held-out-clone validation. A standalone S8ABC run reuses the saved Fig. 4M projections.

Original model-fitting code is in `scripts/` and `R/`.
