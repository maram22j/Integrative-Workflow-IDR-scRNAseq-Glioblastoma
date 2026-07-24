# Integrative Workflow for the Analysis of Intrinsically Disordered Proteins from Single-Cell RNA Sequencing: A Glioblastoma Case Study
[![R](https://img.shields.io/badge/R-4.5+-276DC3.svg)]()
[![Seurat](https://img.shields.io/badge/Seurat-v5-00A087.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
---
## Overview

This repository presents an end-to-end computational workflow for integrating single-cell transcriptomic profiles with protein intrinsic disorder annotations in glioblastoma.

The workflow combines single-cell RNA sequencing (scRNA-seq), protein structural annotation, statistical analyses, and functional enrichment to investigate the contribution of intrinsically disordered proteins (IDPs) and intrinsically disordered regions (IDRs) to glioblastoma biology.


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
- Compare IDR profiles between healthy and malignant cells to identify potential differences associated with the disease state.
- Investigate the biological functions of genes encoding intrinsically disordered proteins and evaluate their potential roles in cancer through functional enrichment analysis.


## Repository Organization
```text
scRNAseq-IDR-Glioblastoma
│
├── data
│   ├── processed
│   └── metadata
│
├── scripts
│   ├── 01_Single_cell_RNA_seq_with_CNA_filtration.R
│   ├── 02_Gene_filtration.R
│   ├── 03_Conversion_adult_to_IDR_percentage.R
│   ├── 04_Conversion_pediatric_to_IDR_percentage.R
│   ├── 05_
│   ├── 06_
│   ├── 07_
│   ├── 08_
│   ├── 09_
│   └── 10_
│
├── results
│   ├── single_cell_analysis
│   │   ├── quality_control/
│   │   ├── variable_features/
│   │   ├── PCA/
│   │   ├── UMAP/
│   │   ├── tSNE/
│   │   └── malignant_cells/
│   ├── differential_expression
│   │   ├── adult_vs_pediatric.csv
│   ├── gene_filtration
│   │   ├── adult_filtered_genes.csv
│   │   └── pediatric_filtered_genes.csv
|   ├── IDR_annotation
│   │   ├── IDR_protein_for_adult.csv
|   |   ├── IDR_protein_for_pediatric.csv
|   └── figures
│       ├── workflows
├── docs
├── images
├── README.md
└── LICENSE
|
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
- Cell filtering using CNA method and malignant cell identification
- Feature selection
- Scaling

### 2.Downstream analysis
- Linear Dimensionality Reduction (PCA)
- Cell Clustering
- Cell annotation
- Non-linear Dimensionality Reduction (UMAP/ tSNE)
- Differential expression analysis

### 3. Gene filtration

Gene filtration based on age_group feature (adult and pediatric).

### 4. Protein Annotation

Gene symbols are converted into UniProt identifiers.

Protein characteristics include

- protein length

### 5. Intrinsic Disorder Annotation

Protein disorder annotations are retrieved from MobiDB.

For each protein:
- number of IDRs
- total IDR length
- percentage disorder
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
| SingleR | cell annotation |
| EnsDb.Hsapiens.v86 | Provides human gene annotations. |
| patchwork | Combines ggplot2 plots. |
| SingleCellExperiment | Stores single-cell RNA-seq data. |
| celldex | Provides reference datasets for cell annotation. |
| EnhancedVolcano | Generates volcano plots for DE analysis. |
| dplyr | Data manipulation |
| biomaRt | Gene annotation |
| UniProt.ws | Protein annotation |
| protti | Protein utilities |
| ggplot2 | Visualization |
| HGNChelper | Corrects and updates gene symbols. |

## Outputs

The workflow produces

- Quality control graphs
- UMAP projections
- Cluster annotations
- Gene filtration tables
- UniProt annotation tables
- IDR annotation tables
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

