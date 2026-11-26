# RNA-based clonal lineage tracing

Identification of LARRY-labelled clones and analysis of their size, lineage composition and mutation-associated transcriptional profiles.

## Data

The **RNA/LARRY** record provides `clone_info.RData`, the experiment-specific LARRY objects, mutation scores and clone-pattern results. Some analyses also use growth estimates from the **DNA clonal tracking** record. Download the records listed in the [dataset index](../publication/datasets.tsv) into the repository root.

Fig. 4B uses the saved Fig. 3A UMAP coordinates supplied with the code. It exports the highlighted clone cells rather than duplicating the complete embedding.

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("RNA_clone_analysis/publication/run_publication.R")
run_rna_clone_publication(
  stages = "4C",
  dry_run = FALSE,
  work_root = file.path(getwd(), "reproduction_clones_4C")
)
```

Choose a new `work_root`. Calling the function without arguments lists the analyses without executing them.

| Stage | Figures |
| --- | --- |
| `4B` | Example clones on the RNA UMAP |
| `4C` | Clone size and lineage-bias distributions |
| `S7BCD` | Barcode calling, sample assignment and clone counts |
| `5E_S8D` | Mutation prediction and clone scores |
| `5GH_S8E` | Lineage-bias means and variation |

## Outputs

Numerical tables and figure PDFs are written under `<work_root>/RNA_clone_analysis/publication/source_data/` and `figures/`. These commands use processed clone assignments and saved model results. The upstream barcode-calling code is in `scripts/` and `R/`.

Clone size and the number of cells visible on a shared embedding are distinct quantities; both are retained where relevant. Fig. 4C uses cell state 1 for HSC bias. The ROC and bootstrap analyses retain the cohorts and resampling definitions used in their source analyses.
