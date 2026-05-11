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


