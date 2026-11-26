# Single-cell RNA analysis

Cell-state annotation, differential gene expression and pathway analysis of WT, DNMT3A-mutant and TET2-mutant hematopoietic cells.

## Data

Download the **RNA/LARRY** dataset listed in the [dataset index](../publication/datasets.tsv) and extract it into the repository root. The main inputs are `RNA_analysis/data/processed/scrna.RData`, `mutant_de.RData`, the Pagoda2 object, hashtag count matrices and gene sets. Small reference files are supplied with the code.

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("RNA_analysis/publication/run_publication.R")
run_rna_publication(
  stages = "3DE",
  dry_run = FALSE,
  work_root = file.path(getwd(), "reproduction_RNA_3DE")
)
```

Choose a new `work_root` for each run. Calling `run_rna_publication()` without arguments lists the available analyses without running them.

| Stage | Figures |
| --- | --- |
| `3A_3B_S6AB` | 3A–B and S6A–B: cell states, genotype distributions and abundance |
| `3DE` | 3D–E: HSC differential expression and selected pathways |
| `3F` | 3F: pathway-association heatmap |
| `S6C` | S6C: comparison of mutation-associated gene effects |
| `S5A` | S5A: hashtag count plots |
| `S5B` | S5B: cell quality control |
| `S5D` | S5D: gene-program scores |

## Outputs and figure definitions

Tables and PDFs are written under `<work_root>/RNA_analysis/publication/source_data/` and `figures/`. The scripts use saved expression and model results; they do not repeat raw-read processing or refit DESeq2.

Two small files in `publication/resources/` define the displayed figures:

| File | Purpose |
| --- | --- |
| `Fig3DE_published_pathway_ids.tsv` | Pathway IDs, labels, order and direction shown in Fig. 3D/E |
| `Fig3F_functional_categories.tsv` | Names assigned to the six pathway clusters in Fig. 3F |

These are figure inputs, not validation reports. Keep them with the scripts. Source Data retain the underlying statistical values separately from the display selections.

S5C and S6D are not included in this runner. S6A exports the sample distributions without reconstructing the manually added significance symbols. Label placement and point order may differ from the assembled figures. The original processing and modelling code is in `scripts/`.
