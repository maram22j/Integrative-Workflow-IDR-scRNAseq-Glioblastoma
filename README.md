# Integrative Workflow for the Analysis of Intrinsically Disordered Proteins from Single-Cell RNA Sequencing: A Glioblastoma Case Study
[![R](https://img.shields.io/badge/R-4.5+-276DC3.svg)]()
[![Seurat](https://img.shields.io/badge/Seurat-v5-00A087.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
---
## Overview

This repository provides an end-to-end computational workflow integrating single-cell RNA sequencing (scRNA-seq) data with protein intrinsic disorder annotations to investigate the role of intrinsically disordered proteins (IDPs) and intrinsically disordered regions (IDRs) in glioblastoma.

The workflow combines:
- single-cell transcriptomic analysis,
- malignant cell identification,
- Pseudobulk differential gene expression analysis,
- protein annotation,
- IDR characterization,
- and functional enrichment analysis.

The objective is to integrate transcriptomic and protein structural information to identify disorder-associated molecular signatures in glioblastoma.


## Biological Motivation

Glioblastoma (GBM) is the most aggressive primary brain tumor in adults and exhibits remarkable cellular heterogeneity.

While scRNA-seq enables the characterization of transcriptional heterogeneity at single-cell resolution, protein structural properties remain largely unexplored in these datasets.

Intrinsically disordered proteins (IDPs) participate in

- transcriptional regulation
- signal transduction
- phase separation
- protein interaction networks
- cancer progression

This project integrates transcriptomic and structural information to identify disorder-associated molecular signatures in malignant glioblastoma cells.

## Objectives

- Identify proteins containing intrinsically disordered regions (IDRs) encoded by genes expressed in the GSE131928 dataset.
- Integrate IDR information with single-cell gene expression data to investigate correlations between protein disorder and gene expression levels.
- Compare IDR profiles between adult and pediatric age group to identify potential differences associated with the disease state.
- Investigate the biological functions of genes encoding intrinsically disordered proteins and evaluate their potential roles in cancer through functional enrichment analysis.


## Repository Organization
```text
scRNAseq-IDR-Glioblastoma/
│
├── data/
│   ├── processed/
│   └── metadata/
│
├── scripts/
│   ├── 01_Single_cell_RNA_seq_with_CNA_filtration.R
│   ├── 02_Pseudobulk_differential_expression_analysis.R
│   ├── 03_Gene_filtration.R
│   ├── 04_IDR_conversion_adult.R
│   ├── 05_IDR_conversion_pediatric.R
│   ├── 06_
│   ├── 07_
│   ├── 08_
│   ├── 09_
│   ├── 10_
│   └── 11_
│
├── results/
│   ├── single_cell_analysis/
│   │   ├── PCA_results/
│   │   └── Plots/
│   │
│   ├── differential_expression/
│   │   ├── Files/
│   │   │   ├── DE_pseudobulk_patient_limma_DEGs.csv
│   │   │   └── DE_pseudobulk_patient_limma_full.csv
│   │   │
│   │   └── plots/
│   │       ├── Boxplot_limma_DEGs_logFC_0.5
│   │       ├── Density_limma_DEGs_logFC_0.5
│   │       ├── PCA_plot_DE_limma_logFC_0.5
│   │       ├── Volcano_plot_DE_limma_logFC_0.5
│   │       └── Heatmap_limma_DEGs
│   │
│   ├── gene_filtration/
│   │   ├── DE_adult.csv
│   │   └── DE_pediatric.csv
│   │
│   ├── IDR_annotation/
│   │   ├── adult/
│   │   │   ├── AIUPred_files/
│   │   │   │   ├── Binding_prediction/
│   │   │   │   ├── disorder_prediction/
│   │   │   │   ├── linker_protein/
│   │   │   │   └── redox_state/
│   │   │   │
│   │   │   ├── FASTA_files/
│   │   │   └── AIUPred_adult_results.csv
│   │   │
│   │   └── pediatric/
│   │       ├── AIUPred_files/
│   │       │   ├── Binding_prediction/
│   │       │   ├── disorder_prediction/
│   │       │   ├── linker_protein/
│   │       │   └── redox_state/
│   │       │
│   │       ├── FASTA_files/
│   │       └── AIUPred_pediatric_results.csv
│   │
│   └── figures/
│       └── workflows/
│
├── docs/
├── images/
├── README.md
└── LICENSE
```
## Data

| Item | Description |
|------|-------------|
| Dataset | GSE131928 |
| Disease | Glioblastoma |
| Species | Homo sapiens |
| Platform | Smart-seq2 |
| Data Type | Single-cell RNA sequencing of adult and paediatric IDH-wildtype Glioblastomas|
## Methodology

### 1. Data preprocessing

- Import expression matrix
- Metadata integration
- Quality control
- Identification of malignant cells using CNA-based approach
- Feature selection
- Scaling
  

### 2.Downstream analysis
- Cell Clustering
- Linear Dimensionality Reduction (PCA)
- Non-linear Dimensionality Reduction (UMAP/ tSNE)
- Pseudobulk differential expression analysis

### 3. Gene filtration

Differentially expressed genes are separated according to patient age groups:

- Adult glioblastoma
- Pediatric glioblastoma

### 4. Protein Annotation

Gene symbols are converted into UniProt identifiers using Ensembl BioMart.

Protein information retrieved:

- UniProt accession IDs
- Protein length
- Protein sequences
  
The resulting UniProt identifiers are then used to retrieve the corresponding protein sequences, which are required for downstream intrinsic disorder prediction.

### 5. Intrinsic Disorder Annotation

Protein intrinsic disorder regions (IDRs) are annotated using AIUPred.

For each protein, the following information is obtained:
- Protein length
- mean_score
- disordered_residues
- disorder_percentage

The percentage of intrinsically disordered regions (IDR percentage) is calculated as:

**IDR percentage = (Total IDR length / Protein length) × 100**

This provides the proportion of each protein sequence predicted to correspond to intrinsically disordered regions.

For proteins with multiple predicted IDR regions, the lengths of all identified IDR regions are summed to obtain the total IDR length before calculating the percentage.
## Installation

Clone the repository:

```bash
git clone https://github.com/maram22j/scRNAseq-IDR-Glioblastoma.git
cd scRNAseq-IDR-Glioblastoma
```

## Software

| Package | Purpose |
|----------|---------|
| Seurat | Single-cell analysis |
| EnsDb.Hsapiens.v86 | Provides human gene annotations. |
| patchwork | Combines ggplot2 plots. |
|Matrix |Provides efficient dense and sparse matrix classes and operations for large-scale data analysis. |
| Limma | Linear modeling and differential expression analysis for transcriptomic data. |
| tibble | data frame format for data manipulation in R. |
| ggrepel | Improves ggplot2 visualizations by preventing overlapping text labels. |
| reshape2 | Provides tools for reshaping and transforming data between wide and long formats. |
|pheatmap| Creates clustered heatmaps for visualizing gene expression patterns across samples or conditions.|
| dplyr | Data manipulation |
| biomaRt | Gene annotation |
| UniProt.ws | Protein annotation |
| protti | Protein utilities |
| ggplot2 | Visualization |
| HGNChelper | Corrects and updates gene symbols. |
|httr| Enables HTTP requests to retrieve data from web APIs.|
|AIUPred| Predicts protein intrinsic disorder and disorder-related functional properties from protein sequences.|
|Biostrings| Provides tools for efficient manipulation and analysis of biological sequences.|



## Outputs
The workflow produces: 

- Quality control plots
- UMAP projections
- Differential expression and gene filtration tables
- IDR annotation tables
- Protein FASTA files
- AIUPred prediction files
- IDR percentage calculations
- Downstream analysis plots and results
  
## Citation

If you use this workflow in your research, please cite this repository and the original GSE131928 dataset.

---

## License

MIT License

---

## Contact

**Maram Nhaili**

Industrial Biology Engineering Student

National Institute of Applied Science and Technology (INSAT)

Tunisia

Email: maram.nhaili@insat.ucar.tn

