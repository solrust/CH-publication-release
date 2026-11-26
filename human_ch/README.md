# Comparison with human clonal hematopoiesis

Comparison of DNMT3A- and TET2-associated transcriptional changes in xenografted hematopoietic cells and human clonal hematopoiesis datasets. The analyses cover HSC-to-myeloid differentiation, pathway enrichment and agreement between cell populations.

## Data

The **RNA/LARRY** [Zenodo record](https://doi.org/10.5281/zenodo.24545548) supplies the processed inputs under `human_ch/data/`:

| Input | Analysis |
| --- | --- |
| `human_ch_seurat.rds` | Human HSC-to-GMP pseudotime, Fig. 3C |
| `datasets_DE_GSEA.rds` | Pathway comparisons, Fig. 3G |
| `integrated_seurat_metadata.rds` | Joint embedding and cell annotations, Fig. S5E |
| `integrated_seurat.rds` | Cell-population expression comparisons, Fig. S5F |

Fig. 3C also uses `RNA_analysis/data/processed/scrna.RData` and `hashtag_meta.rds` from the same record. External reference datasets are described in the article's Data Availability statement.

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("human_ch/publication/run_publication.R")

human_ch_export(panels = c("3G", "S5E"), validate_replay = TRUE)
hc_s5f_export()
hc_fig3c_export()
```

Run these commands in a fresh checkout with the data restored. They write to `human_ch/publication/source_data/` and `figures/` and refuse to replace existing outputs.

## Outputs and interpretation

The outputs contain cell coordinates, pseudotime distributions, pathway statistics and cell-population expression correlations, with corresponding PDFs. Fig. 3C retains the distinction between WT cells from CH samples and WT controls. Undefined pathway statistics remain missing rather than being replaced by zero.

These analyses start from processed objects, not raw sequencing reads. Fig. 3C recalculates the PCA-based ordering; S5F calculates population means from the saved integrated expression data. The PDFs reproduce analytical panels, not the complete manuscript layouts. The original integration and comparison analyses are in `scripts/` and `2_run_integrations_comparison.R`.
