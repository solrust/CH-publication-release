# Data and analysis setup

## 1. Download the code and data

```bash
git clone https://github.com/solrust/CH-publication-release.git
cd CH-publication-release
```

Check out the code version cited in the article. Download the ZIP archives and SHA-256 manifests from the Zenodo records listed in the [README](../README.md). Each ZIP is a separate archive; the parts do not need to be joined.

For example, save the downloads in a directory named `zenodo_downloads` alongside the Git checkout, using the dataset names below. From the checkout root, extract the archives as follows:

```bash
(
  set -e
  for rel in \
    RNA_LARRY/RNA_LARRY_part01.zip \
    RNA_LARRY/RNA_LARRY_part02.zip \
    Multiome/Multiome_part01.zip \
    Multiome/Multiome_part02.zip \
    Multiome/Multiome_part03.zip \
    TF_perturbation/TF_perturbation_part01.zip \
    DNA_clonal_tracking/DNA_clonal_tracking_part01.zip
  do
    unzip -n -q "../zenodo_downloads/$rel" -d .
  done
)
```

Extract only the records needed for your analysis if disk space is limited. The `-n` option preserves files already present in Git. Four small files are present in both the code repository and data archives; the copies were checked to be identical when the datasets were prepared. Verify extracted files against the SHA-256 manifests provided with each dataset. Do not retain a conflicting copy if its checksum differs.

The extracted data should be under the corresponding analysis directories (for example, `RNA_analysis/data/processed/` and `multiome/data/processed/`). The motif association workbooks should be under `pipeline/resource/` within the checkout. The [dataset registry](datasets.tsv) lists the archive filenames and their DOIs.

### Which datasets are needed?

| Analysis | Zenodo datasets |
| --- | --- |
| RNA, comparative human CH | RNA/LARRY |
| Clone phenotypes and gene associations | RNA/LARRY; DNA clonal tracking for growth-rate analyses |
| Longitudinal DNA tracking | DNA clonal tracking, RNA/LARRY |
| In vitro split-pool experiment (`clones/`) | DNA clonal tracking |
| Multiome | Multiome, RNA/LARRY for RNA clone comparisons |
| TF perturbations | TF perturbation; RNA/LARRY for gene sets and clonal associations |
| MYC, inflammatory HSCs | RNA/LARRY, Multiome; DNA clonal tracking for growth-rate analyses |
| Target-site editing plots | RNA/LARRY and the numerical table included in Git |

## 2. R environment

The analyses were run with R 4.2.3, Bioconductor 3.16 and Matrix 1.6-4. Package versions are listed in [packages_validated.tsv](environment/packages_validated.tsv), with additional details in [R_environment_validated.txt](environment/R_environment_validated.txt). Install the listed compatible packages and their system dependencies in a dedicated R library. The package inventory records the original environment; it is not an automated installer.

In a fresh R session, set the working directory to the repository root and specify your R library before loading analysis packages:

```r
Sys.setenv(CH_R_LIB = "/absolute/path/to/R/library")
source("publication/tools/configure_environment.R")
ch_configure_environment()
```

The initialization function checks the R and Matrix versions and prepares the required directories. Set `CH_R_LIB` explicitly; no cluster-specific directory is required.

## 3. Run an analysis

Each project has one README in its top-level directory and a `publication/run_publication.R` script. The project README describes the available functions, required datasets and output locations. For example, the following regenerates the hashtag count plots in Supplementary Fig. 5A:

```r
source("RNA_analysis/publication/run_publication.R")
run_rna_publication(stages = "S5A", dry_run = FALSE,
                    work_root = file.path(getwd(), "reproduction_S5A"))
```

For the RNA, RNA-clone and association runners, omitting `dry_run = FALSE` displays the execution plan without running the analyses. Use a new output directory when repeating an analysis.

The following examples use inputs from the other three datasets:

```r
# Multiome and RNA/LARRY: MYC motif associations (Supplementary Fig. 10D)
source("MYC_analysis/publication/run_publication.R")
run_myc_publication(panels = "S10D",
                    output_root = file.path(getwd(), "MYC_analysis/publication/reproduction_S10D"))

# TF perturbation (Fig. 7C and D)
source("TF_overexpression/publication/run_publication.R")
run_tf_publication(panels = c("7C", "7D"),
                   output_root = file.path(getwd(), "TF_overexpression/publication/reproduction_7CD"))

# DNA split-pool lineage tracing (Fig. 5L and Supplementary Fig. 8H and I)
source("clones/publication/run_publication.R")
clones_pub_run(output_root = file.path(getwd(), "clones/publication/reproduction_GS175"),
               validate = TRUE)
```

These commands generate numerical tables and plots from the archived inputs. Individual project guides describe which results are recalculated and which use previously saved numerical results. Some analyses require substantial memory. Supplementary Fig. 5C and 6D are not part of the routine RNA runner; see the RNA project guide for details.

The original processing scripts are also provided for reference. Reprocessing sequencing reads, recounting ATAC fragments or BAM alignments, and repeating upstream model fits may require additional data, software and path configuration. These steps are distinct from regenerating the figures with the processed datasets.
