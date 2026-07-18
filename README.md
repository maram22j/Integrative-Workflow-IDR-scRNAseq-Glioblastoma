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

- Perform preprocessing and quality control of scRNA-seq data
- Identify malignant glioblastoma cells
- Map expressed genes to UniProt proteins
- Retrieve IDR annotations from MobiDB
- Quantify protein disorder characteristics
- Integrate IDR information with gene expression

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
│   ├── 02_Differential_expression_analysis.R
│   ├── 03_Gene_filtration.R
│   ├── 04_Conversion_adult_to_IDR_percentage.R
│   ├── 05_Conversion_pediatric_to_IDR_percentage.R
│   ├── 06_mobidb_annotation.R
│   ├── 07_expression_integration.R
│   ├── 08_statistical_analysis.R
│   ├── 09_functional_enrichment.R
│   └── 10_visualization.R
│
├── results
│   ├── figures
│   ├── tables
│   └── supplementary
│
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

