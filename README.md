# DESeq2 Analysis for SMOC2

## Introduction
This repository contains a script that serves as a basic outline or blueprint for performing Differential Expression (DE) analysis using DESeq2, with a specific focus on the SMOC2 gene. It illustrates the typical workflow for DE analysis, including data preprocessing, normalization, visualization, and differential expression testing.

## Script Description
The `DESeq2_smoc2_DE_analysis.r` script performs several key tasks in DE analysis:
- Data loading and cleaning from `smoc2_rawcounts.csv` and `smoc2_metadata.csv`.
- Quality control checks using principal component analysis (PCA) and heatmap visualizations.
- Differential expression analysis of the SMOC2 gene across different conditions and genotypes.
- Generation of various plots including PCA plots, MA plots, and volcano plots to visualize the results.
- Extraction and listing of significantly differentially expressed genes.

## Installation
To run this script, ensure R is installed on your system along with the following R libraries:
- DESeq2
- RColorBrewer
- pheatmap
- tidyverse
- ggplot2

Install these packages using the R command:
```R
install.packages(c("DESeq2", "RColorBrewer", "pheatmap", "tidyverse", "ggplot2"))
```

## Usage
Follow these steps to use the script:

1. Clone or download the script to your local machine.
2. Place the data files smoc2_rawcounts.csv and smoc2_metadata.csv in '~/Documents/git/DESeq2/data/'.
3. Open the script in R or an R IDE (e.g., RStudio).
4. Run the script. Modify parameters or input files as needed.

## Data Files
The script utilizes:

-smoc2_metadata.csv: Metadata for the samples.
-smoc2_rawcounts.csv: Raw counts of RNAseq data.
Ensure these files are in the specified directory for the analysis.

## Acknowledgements
This analysis uses the DESeq2 package, a standard in analyzing count-based NGS data like RNAseq, for differential expression analysis and quality control checks.
