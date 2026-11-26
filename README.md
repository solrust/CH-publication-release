# Clonal hematopoiesis

Analysis code for our study of clonal behavior in human hematopoietic cells carrying DNMT3A or TET2 mutations. The analyses combine LARRY lineage tracing with single-cell RNA sequencing, chromatin accessibility profiling and transcription factor perturbations.

## Data

Processed datasets are listed below. See each record for its files and access information:

| Dataset | DOI |
| --- | --- |
| Single-cell RNA-seq and LARRY lineage tracing | [10.5281/zenodo.24545548](https://doi.org/10.5281/zenodo.24545548) |
| Single-cell multiome | [10.5281/zenodo.24546656](https://doi.org/10.5281/zenodo.24546656) |
| Transcription factor perturbations | [10.5281/zenodo.24545415](https://doi.org/10.5281/zenodo.24545415) |
| DNA barcode sequencing and clonal tracking | [10.5281/zenodo.24545120](https://doi.org/10.5281/zenodo.24545120) |

The [dataset index](publication/datasets.tsv) lists archive names and cross-project dependencies. See the article's Data Availability statement for sequencing data and external reference datasets.

## Getting started

Clone this repository and extract the relevant Zenodo ZIP archives into the repository root, preserving their directory structure. The analyses use R 4.2.3, Bioconductor 3.16 and Matrix 1.6-4.

See the [installation guide](publication/RELEASE_INSTALLATION.md) for data setup, R package versions and examples of how to run the analyses.

## Analyses

| Directory | Analysis |
| --- | --- |
| [RNA_analysis](RNA_analysis/README.md) | Cell states and differential gene expression |
| [RNA_clone_analysis](RNA_clone_analysis/README.md) | Clonal assignments and phenotypes |
| [Association_analysis](Association_analysis/README.md) | Gene expression associations with clonal phenotypes |
| [DNA_analysis](DNA_analysis/README.md) | Longitudinal clonal tracking and growth rates |
| [clones](clones/README.md) | In vitro split-pool lineage tracing |
| [human_ch](human_ch/README.md) | Comparison with human CH datasets |
| [multiome](multiome/README.md) | Gene expression, chromatin accessibility and motif analysis |
| [TF_overexpression](TF_overexpression/README.md) | Transcription factor perturbation experiments |
| [MYC_analysis](MYC_analysis/README.md) | MYC-associated clonal phenotypes |
| [inflammatory_HSCs](inflammatory_HSCs/README.md) | Inflammatory HSC signatures |
| [indels_calling](indels_calling/README.md) | Target-site editing analysis |

The `scripts/` directories contain the original analyses. The `publication/` directories contain code for generating figure data and plots from the deposited processed data and, where indicated, saved numerical results. The README in each project describes its inputs, commands and outputs. Figure-selection tables are kept in `publication/resources/` within the projects that use them. Saved numerical inputs distributed with the code are listed in [checkpoints.tsv](publication/checkpoints.tsv).

## Citation

Please cite the accompanying article and the relevant Zenodo dataset(s) when using this code or data.
