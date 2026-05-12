# Integrative Multi-Omic Data Analysis Practicals

This repository contains the practical material for the course **Integrative Multi-Omic Data Analysis for Biological and Biomedical Data**.

The practicals introduce basic preprocessing, quality control, exploratory analysis, and multi-omics integration using **MOFA2**.

## Repository structure

```text
.
├── day1/
│   ├── Day1_single_omics_practical.qmd
│   ├── RNAseq_count.tsv
│   ├── Methylation_Mvalues.tsv
│   ├── Genetics_mutationStatus.tsv
│   ├── drugScreen.tsv
│   └── metadata.tsv
│
├── day2/
│   ├── Day2_MOFA_practical.qmd
│   └── reactomeGS.rda
│
└── final_assignment/
    ├── AML_MOFA_exam_student_assignment.html
    ├── rna_count.tsv
    ├── gene_mutation.tsv
    ├── drug_auc.tsv
    └── metadata.tsv
```

## Recommended computational environment

The practicals are written as **Quarto (`.qmd`) documents** and are intended to be run in **RStudio**.

Please install the following software before the course:

1. **R** version 4.3 or newer  
2. **RStudio Desktop**  
3. **Quarto**  
4. The required R and Bioconductor packages listed below

## Installing required R packages

Open R or RStudio and run the following commands.

```r
# Install CRAN packages
install.packages(c(
  "tidyverse",
  "data.table",
  "matrixStats",
  "knitr",
  "ggpubr",
  "UpSetR",
  "survival",
  "BiocManager"
))

# Install Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "SummarizedExperiment",
  "MOFA2"
))
```

MOFA2 requires a Python backend called `mofapy2`. After installing `MOFA2`, run:

```r
library(MOFA2)
install_mofapy2(method = "basilisk")
```

To check whether MOFA2 is working, run:

```r
library(MOFA2)
library(DESeq2)
library(tidyverse)
```

If these commands run without error, the environment should be ready for the practicals.

## Troubleshooting

If package installation fails, please first check that:

- your R version is recent enough;
- `BiocManager` is installed;
- Quarto is installed and available in RStudio;
- the working directory is set to the repository folder;
- the input data files are in the expected subfolders.

For MOFA2-specific issues, try reinstalling the Python backend:

```r
library(MOFA2)
install_mofapy2(method = "basilisk", force = TRUE)
```

If this still doesn't work, a Docker image will be provided. Instructions will be shared during the first practical session. 

## Final assignment

The folder `final_assignment/` contains the dataset and instructions for the course assignment. The student-facing assignment is provided as an HTML file.

Please follow the submission instructions given during the course.


## Python dependencies and data for day 3

If you do not use the provided Docker container, you can reproduce the Python
environments locally with [`uv`](https://docs.astral.sh/uv/). Install `uv` once, for
example with `curl -LsSf https://astral.sh/uv/install.sh | sh` on Linux/macOS, or
follow the official installation instructions for your operating system. In general,
`uv init` creates a new Python project and `uv add <package>` adds dependencies, but
you do not need to do this for the course material: the dependencies are already
declared in the individual project folders such as `00_prep_python_uv`
`01_single_cell_processing`, `02_batch_integration`, and so on. To set up one of these
environments, open a terminal, change into the folder that contains the `pyproject
toml` file, and run `uv sync`. This will create the local `.venv` environment and
install the pinned dependencies from `uv.lock`. After that, run commands from the same
folder with `uv run ...`, or activate the environment manually with `source .venv/bin
activate`.

`04_cite_seq_application` and of course `XX_exam_data` require extra data that cannot be loaded like the tutorial data. If you are not using Dokcer, you may download under the following link:
https://heibox.uni-heidelberg.de/d/dcb5d295521b47108715/

## Day 3 repository structure

The `day3` directory is organized as a set of mostly
self-contained teaching modules. Each numbered folder
corresponds to one topic in the course, contains its own
notebooks and supporting files, and declares its Python
environment with a `pyproject.toml` and `uv.lock`. To run a
module locally, change into the corresponding folder and run
`uv sync`. The `data`, `models`, `results`, and
`prepare_data_scripts` folders contain the files needed for
the practicals or for preparing them. `XX_exam_data` contains 
the exam data preparation material in the same project 
structure. The `conda_envs` folder is included only as an 
example for users who are more familiar with Conda, since 
Conda environments are still widely used in scientific 
computing; the primary reproducible setup for this course 
uses `uv`.

```text
day3
├── 00_prep_python_uv
│   ├── notebooks
│   ├── pyproject.toml
│   ├── README.md
│   └── uv.lock
├── 01_single_cell_processing
│   ├── data
│   ├── notebooks
│   ├── pyproject.toml
│   ├── README.md
│   └── uv.lock
├── 02_batch_integration
│   ├── data
│   ├── main.py
│   ├── notebooks
│   ├── pyproject.toml
│   ├── README.md
│   └── uv.lock
├── 03_cite_seq
│   ├── data
│   ├── _extensions
│   ├── layouts
│   ├── models
│   ├── notebooks
│   ├── pyproject.toml
│   ├── README.md
│   ├── results
│   └── uv.lock
├── 04_cite_seq_application
│   ├── data
│   ├── notebooks
│   ├── prepare_data_scripts
│   ├── pyproject.toml
│   ├── README.md
│   └── uv.lock
├── 05_sc_atac_seq
│   ├── notebooks
│   ├── pyproject.toml
│   ├── README.md
│   └── uv.lock
├── conda_envs
│   ├── lecture__processing.yaml
│   └── lecture__r_interop.yaml
├── Dockerfile
└── XX_exam_data
    ├── data
    ├── notebooks
    ├── pyproject.toml
    ├── README.md
    └── uv.lock
```