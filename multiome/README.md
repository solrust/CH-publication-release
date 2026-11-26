# Single-cell multiome analysis

Joint analysis of RNA expression, chromatin accessibility, transcription factor motifs and LARRY clone identities. These workflows generate the computational panels of Fig. 6 and Supplementary Figs. 9–10.

## Data

Download **Multiome** and **RNA/LARRY** from the [dataset index](../publication/datasets.tsv). The main inputs are the processed Seurat object, cell and clone annotations, motif-association workbooks and saved pathway, pseudobulk and prediction results.

Extract both records into the repository root. The motif workbooks belong in `pipeline/resource/` inside that root. Small saved results supplied with the code are listed in [checkpoints.tsv](../publication/checkpoints.tsv).

## Run

Run from the repository root after completing the [environment setup](../publication/RELEASE_INSTALLATION.md).

```r
source("multiome/publication/run_publication.R")
result <- run_multiome_publication(
  mode = "checkpoint",
  output_root = "multiome/publication/reproduction"
)
```

Use a new output directory under `multiome/publication/`. Loading the Seurat object requires substantially more memory than its compressed file size.

## Outputs and scope

The run produces 15 numerical tables in `source_data/` and 17 PDF components in `figures/`, beneath the requested output directory.

`mode = "checkpoint"` regenerates 12 tables from processed inputs and reuses the supplied tables for Fig. 6F, 6H and S9A. It also uses saved motif-enrichment, pseudobulk and prediction results. `mode = "full"` additionally recalculates those three tables but still starts from processed data; fragment-based steps need the original fragment files.

Fig. 6A/E are schematics. S9D and S10E/G are outside this runner; S10B/C are covered by `inflammatory_HSCs` and S10D by `MYC_analysis`. The regenerated Fig. 6H coverage plot does not include the full gene-model track in the assembled artwork. Fig. 6J reports variability across repeated cell-level prediction splits. Original analyses are retained in `scripts/`.
