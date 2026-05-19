# Integrative Multi-Omic Data Analysis Practicals

This repository contains the practical materials for the course **Integrative Multi-Omic Data Analysis for Biological and Biomedical Data**.

The course introduces practical workflows for preprocessing, quality control, exploratory analysis, and integration of multi-omic datasets. Days 1 and 2 focus on bulk multi-omics integration using a Chronic Lymphocytic Leukemia (CLL) cohort and **Multi-Omics Factor Analysis (MOFA)**. Day 3 focuses on single-cell multi-omics, and Day 4 introduces MOFAcell / knowledge-based interpretation of multi-cellular data.

---

## Repository overview

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
├── day3/
│   ├── 00_prep_python_uv/
│   ├── 01_single_cell_processing/
│   ├── 02_batch_integration/
│   ├── 03_cite_seq/
│   ├── 04_cite_seq_application/
│   ├── 05_sc_atac_seq/
│   ├── conda_envs/
│   ├── Dockerfile
│   └── XX_exam_data/
│
└── final_assignment/
    ├── AML_MOFA_exam_student_assignment.html
    ├── rna_count.tsv
    ├── gene_mutation.tsv
    ├── drug_auc.tsv
    └── metadata.tsv
```

---

## Quick start

For **Days 1 and 2**, please install the following before the course:

1. **R** version 4.3 or newer
2. **RStudio Desktop**
3. **Quarto**
4. Required R and Bioconductor packages listed below

The Day 1 and Day 2 practicals are written as **Quarto (`.qmd`) documents** and are intended to be run in **RStudio**.

---

## Installing software

### R

Download and install R from:

<https://cran.r-project.org/>

Please use **R 4.3 or newer**.

### RStudio Desktop

Download and install RStudio Desktop from:

<https://posit.co/download/rstudio-desktop/>


If you are not familiar with Quarto, we recommend that you read: https://quarto.org/. We will also give a brief introduction during the first practical session. 

---

## Installing R packages for Days 1 and 2

Open RStudio and run the following commands in the R console.

```r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN packages
install.packages(c(
  "tidyverse",
  "data.table",
  "matrixStats",
  "knitr",
  "ggpubr",
  "UpSetR",
  "psych",
  "survival",
  "survminer",
  "pheatmap",
  "RColorBrewer",
  "patchwork"
))

# Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "SummarizedExperiment",
  "MOFA2"
))
```

### MOFA2 backend

`MOFA2` requires a Python backend called `mofapy2`. In many cases, this can be automatically installed when the MOFA2 R package is installed.

If this fails, please check the official MOFA2 installation guide:

<https://biofam.github.io/MOFA2/installation.html>

---

## Test your setup before the course

Please run the following in RStudio:

```r
library(tidyverse)
library(DESeq2)
library(MOFA2)
library(UpSetR)

sessionInfo()
```

If all packages load without errors, your R environment is ready for Days 1 and 2.

You can also open `day1/Day1_single_omics_practical.qmd` in RStudio and render it with the command:

```r
quarto::quarto_render("day1/Day1_single_omics_practical.qmd")
```

Alternatively, click **Render** in RStudio.

---

## Day 1: single-omics preprocessing and exploratory analysis

Folder: `day1/`

Main file:

```text
day1/Day1_single_omics_practical.qmd
```

Day 1 covers:

- loading multi-omic data tables and sample metadata
- checking sample overlap across omics layers
- RNA-seq preprocessing with `DESeq2`
- methylation data quality control
- mutation matrix exploration
- drug response data preprocessing
- exploratory PCA for individual omics layers
- comparison of single-omics structure before integration

---

## Day 2: MOFA integration

Folder: `day2/`

Main file:

```text
day2/Day2_MOFA_practical.qmd
```

Day 2 builds on Day 1 and covers:

- preparing integration-ready matrices
- creating and training a MOFA model
- inspecting variance explained by latent factors
- associating factors with metadata
- interpreting feature weights
- pathway-level interpretation of MOFA factors
- linking latent factors to biological and clinical hypotheses

---

## Day 3: single-cell multi-omics

Folder: `day3/`

Day 3 uses Python-based single-cell analysis workflows. Because the Python environment is more complex, we recommend using the provided **Docker image** (will be provided prior to the day 3 practical session) for this part of the course.

If you prefer a local setup, the Day 3 modules use [`uv`](https://docs.astral.sh/uv/) for reproducible Python environments. Each module contains its own `pyproject.toml` and `uv.lock` file.

To set up one module locally:

```bash
cd day3/01_single_cell_processing
uv sync
uv run jupyter lab
```

or run Python scripts with:

```bash
uv run python main.py
```

Some Day 3 modules require additional data. If you are not using Docker, download the data from:

<https://heibox.uni-heidelberg.de/d/dcb5d295521b47108715/> (data will be uploaded later)

For the scATAC-seq module using `snapATAC2`, a Rust compiler may be required if no precompiled binary is available for your system:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

---

## Day 4: MOFAcell / knowledge-based interpretation

Day 4 focuses on MOFAcell and multi-cellular factor interpretation.

We will mainly use **Google Colab** to keep the setup lightweight. As a backup, a Conda environment will be provided. If you do not already have Conda, we recommend installing Miniconda:

<https://www.anaconda.com/docs/getting-started/miniconda/install/overview>

You can find everything we need for the day here: https://drive.google.com/drive/folders/1EI9_uoCQ7AIdq8X-c3r9RPDSApxY0NRY?usp=sharing

---

## Final assignment

Folder: `final_assignment/`

The student-facing assignment is provided as:

```text
final_assignment/AML_MOFA_exam_student_assignment.html
```

The main assignment uses a synthetic AML multi-omic dataset and asks you to apply concepts from the practicals, including preprocessing, exploratory analysis, and MOFA-based integration.

Please follow the submission instructions given during the course.

An optional/bonus assignment of analysing a single-cell COVID-19 CITE-seq dataset is also provided in the day3/XX_exam_data folder and detailed will be introduced in the Day 3 practical session. 

---


## Recommended preparation before the course

Before the first practical, please try to:

1. Install R, RStudio, and Quarto.
2. Install the required R packages.
3. Clone or download this repository.
4. Open `day1/Day1_single_omics_practical.qmd` in RStudio.
5. Check that the setup test above runs without errors.

Do not worry if everything is not perfect. We will reserve some time during the practical sessions for setup and troubleshooting.
