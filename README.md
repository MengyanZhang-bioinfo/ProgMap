# ProgMap
# Analysis Code for: [ProgMap defines clinically actionable cancer progression trajectories]

![Language](https://img.shields.io/badge/Language-R-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
## Overview

This repository contains the source code and scripts used for the data analysis and visualization presented in the manuscript **"[ProgMap defines clinically actionable cancer progression trajectories]"**. The study integrates single-cell RNA sequencing (scRNA-seq), spatial transcriptomics, and bulk RNA-seq data to elucidate the molecular drivers of tumor progression and immune microenvironment remodeling.

## 📂 Repository Structure

The scripts are organized into four analytical modules mirroring the study's workflow:

### 1. Single-Cell & Spatial Transcriptomics Analysis
Scripts for processing high-resolution omics data to define cell identities and spatial architecture.

| Script Name | Description | Key Methods/Packages |
| :--- | :--- | :--- |
| **`Scrna-seq annotation.R`** | Quality control, dimensionality reduction, clustering, and cell type annotation of scRNA-seq data. | `Seurat`, `SingleR` |
| **`Infercnv.R`** | Inference of large-scale copy number variations (CNVs) to distinguish malignant epithelial cells from non-malignant cells. | `infercnv` |
| **`Pseudotime analysis.R`** | Trajectory inference to reconstruct the developmental lineage of malignant cells and map CPSig dynamics. | `Monocle2` / `Monocle3` |
| **`SCT.R`** | Processing of spatial transcriptomics data, including normalization, clustering, and spatial mapping of gene signatures (e.g., via AUCell). | `Seurat`, `AUCell` |

### 2. Bulk Data Processing & Molecular Stratification
Scripts for bridging single-cell insights to large-scale bulk cohorts.

| Script Name | Description | Key Methods/Packages |
| :--- | :--- | :--- |
| **`Pre-treatment.R`** | Standardization and preprocessing of raw transcriptomic data from bulk cohorts (e.g., TCGA, GEO). | `limma`, `sva` |
| **`Deconvolution.R`** | Estimation of cellular abundance in bulk tissues using scRNA-seq derived signatures. | `CIBERSORTx`, `MuSiC` |
| **`Clustering.R`** | Unsupervised hierarchical clustering of samples based on CPSig-related metrics (e.g., MECor values) to define distinct molecular subtypes. | `ConsensusClusterPlus` |

### 3. Functional & Genomic Characterization
Scripts for exploring biological pathways and genomic instability.

| Script Name | Description | Key Methods/Packages |
| :--- | :--- | :--- |
| **`Functional enrichment analysis.R`** | Functional annotation of CPSig and differentially expressed genes (GO, KEGG, GSEA). | `clusterProfiler` |
| **`EMT.R`** | Quantification of Epithelial-Mesenchymal Transition (EMT) activity scores across samples. | `GSVA` / `ssGSEA` |
| **`Mutation analysis.R`** | Analysis of somatic mutation landscapes and calculation of Tumor Mutational Burden (TMB). | `maftools` |

### 4. Immune Landscape & Clinical Translation
Scripts for evaluating the tumor microenvironment (TME) and clinical relevance.

| Script Name | Description | Key Methods/Packages |
| :--- | :--- | :--- |
| **`ESTIMATE.R`** | Calculation of ImmuneScore and StromalScore. | `ESTIMATE` |
| **`T cell exhaustion+Check point.R`** | Evaluation of T-cell exhaustion markers and immune checkpoint expression profiles. | `ggplot2`, `ggpubr` |
| **`TIDE.R`** | Prediction of immunotherapy response and tumor immune dysfunction/exclusion scores. | `TIDE` algorithm |
| **`Survival analysis.R`** | Kaplan-Meier survival curves and Cox proportional hazards regression to assess prognostic value. | `survival`, `survminer` |

### 5. Plots
| Script Name | Description | Key Methods/Packages |
| :--- | :--- | :--- |
| **`UpSet_plot.R`** |The intersection and specificity of the CPSig gene in pan-cancer | `UpSetR` package | 

## 🛠️ Prerequisites & Installation

To reproduce the analysis, ensure **R (version >= 4.0.0)** is installed. Key dependencies include:

```r
# Core Bioinformatics Packages
install.packages(c("Seurat", "Monocle", "survival", "survminer", "estimate", "maftools"))

# Tidyverse & Visualization
install.packages(c("tidyverse", "ggplot2", "pheatmap", "ggpubr"))

# Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "ComplexHeatmap", "clusterProfiler", "infercnv", "GSVA"))
