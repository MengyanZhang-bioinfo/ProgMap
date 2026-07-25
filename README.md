# ProgMap

![Language](https://img.shields.io/badge/Language-R%20%7C%20Python-blue.svg)
![Python](https://img.shields.io/badge/Python-3.9--3.11-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
[![Linux CI](https://github.com/MengyanZhang-bioinfo/ProgMap/actions/workflows/linux-ci.yml/badge.svg)](https://github.com/MengyanZhang-bioinfo/ProgMap/actions/workflows/linux-ci.yml)

ProgMap is a command-line Python package for three-class molecular-stage modeling and progression-feature attribution from paired gene-expression and DNA-methylation matrices. The repository also includes optional R scripts for downstream analysis and visualization.

## ProgMap Python package

The Python implementation performs fold-specific data preprocessing, MECor feature construction, nested epoch selection, three-class deep-learning classification, enhanced integrated-gradient attribution, and statistical identification of cancer progression signatures.

### Quick installation

```bash
git clone https://github.com/MengyanZhang-bioinfo/ProgMap.git
cd ProgMap/progmap-python
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-linux-cpu.txt
python -m pip install --no-deps .
progmap --version
```

The reference environment uses Python 3.11 and TensorFlow/Keras 2.14.0. A Conda environment, CPU and GPU Dockerfiles, a Slurm script, and pinned Linux dependencies are included in [`progmap-python/`](progmap-python/).

### Input data

The input root contains one directory per dataset. Each dataset directory contains six gene-by-sample CSV matrices:

| Class | Expression | DNA methylation |
|---|---|---|
| Normal | `en.csv` | `mn.csv` |
| Stage I | `e1.csv` or `exp1.csv` | `m1.csv` or `me1.csv` |
| Stage II/III | `e2.csv` or `exp2.csv` | `m2.csv` or `me2.csv` |

```text
/data/progmap_inputs/
└── example_dataset/
    ├── en.csv
    ├── e1.csv
    ├── e2.csv
    ├── mn.csv
    ├── m1.csv
    └── m2.csv
```

The first column contains gene identifiers and the remaining columns contain samples. Expression and methylation matrices are aligned by common gene and sample identifiers. Automatic discovery includes every subdirectory containing all six required matrices and ignores incomplete subdirectories.

### One-command analysis

Run all complete dataset directories with the default settings:

```bash
progmap --data-root /data/progmap_inputs --output /results/progmap
```

Run selected datasets with custom settings:

```bash
progmap \
  --data-root /data/progmap_inputs \
  --output /results/custom_run \
  --cancers dataset_A,dataset_B \
  --folds 3 \
  --inner-folds 3 \
  --epochs 1000 \
  --learning-rate 0.001 \
  --early-stopping nested \
  --early-stopping-monitor val_loss \
  --warmup-epochs 20 \
  --patience 50 \
  --seed 42 \
  --correlation-method spearman \
  --test wilcoxon \
  --top-n 100
```

### Default workflow

- Stratified three-fold outer cross-validation estimates held-out performance.
- Stratified three-fold inner cross-validation selects the training duration independently within each outer-training fold.
- Inner models use class-balanced validation loss, a maximum of 1,000 epochs, 20 warm-up epochs, patience of 50, and a minimum improvement of `1e-4`.
- The median inner-fold best epoch is used to train a new model on the complete outer-training fold.
- Expression and methylation data are independently imputed and min-max scaled to `[0, 1]` within each training partition.
- Per-gene Pearson correlation is calculated within each training partition before constructing MECor features. Spearman correlation is available with `--correlation-method spearman`.
- The model uses Dense(2048)-Dense(128) hidden layers, dropout, an input skip connection, and a three-class softmax output.
- Adam optimization uses an initial learning rate of `1e-3`, exponential decay, batch size 16, and seed 42.
- Enhanced integrated gradients quantify held-out feature attribution.
- A one-sided Welch t-test with Benjamini-Hochberg FDR correction at 0.01 identifies significant genes by default.

Outer test samples are excluded from imputation, scaling, correlation calculation, MECor construction, epoch selection, and model fitting.

### Configurable options

The command-line interface supports:

- maximum epochs, learning rate, batch size, dropout, and random seed;
- outer and inner cross-validation fold counts;
- nested early stopping, one internal holdout, or fixed-epoch training;
- validation-loss or validation-AUC monitoring, patience, warm-up, and minimum improvement;
- Pearson or Spearman expression-methylation correlation;
- t-test, Wilcoxon rank-sum, or permutation feature testing;
- multiple-testing method and adjusted significance threshold;
- all significant genes, the top `N` ranked genes, or all ranked genes;
- optional feature-effect bootstrap confidence intervals;
- CPU/GPU selection, thread count, model saving, attribution saving, and progress display.

Use `progmap --help` for the complete parameter list.

### Main outputs

Each dataset-specific result directory contains:

- `cross_validated_predictions.csv` and `fold_metrics.csv`;
- `pooled_out_of_fold_metrics.json`;
- `epoch_selection_all_folds.csv`;
- `features_by_class_all.csv`;
- `features_ranked_genes.csv`;
- `features_selected.csv` and `progression_genes.csv`;
- raw class-specific attribution arrays;
- one saved model, preprocessor, correlation table, training history, and sample-role table per outer fold.

`progression_genes.csv` includes the selected gene, associated stage class, target and background median absolute attribution, raw P value, adjusted P value, attribution effect, significance flag, statistical test, and correction method.

### Synthetic example

A small synthetic `DEMO_DATASET` is included for testing the complete workflow. It follows the required six-file layout but contains no patient data.

```bash
cd progmap-python
progmap \
  --data-root examples/data/datasets \
  --cancers DEMO_DATASET \
  --output demo_results \
  --folds 2 \
  --inner-folds 2 \
  --epochs 2 \
  --patience 0 \
  --warmup-epochs 0 \
  --ig-steps 4 \
  --ig-baselines 1 \
  --top-n 5 \
  --device cpu
```

See the [complete Python package documentation](progmap-python/README.md) for Conda, Docker, Slurm, server execution, all default values, output structure, and additional examples.

---

## Additional analysis scripts

### Overview

The R scripts provide optional workflows for single-cell RNA sequencing, spatial transcriptomics, bulk transcriptomics, molecular characterization, clinical analyses, and visualization. They can be adapted to compatible user datasets.

### Repository structure

#### 1. Single-cell and spatial transcriptomics

| Script | Description | Main packages or methods |
|---|---|---|
| `Scrna-seq_annotation.R` | Quality control, dimensionality reduction, clustering, and cell-type annotation. | Seurat, SingleR |
| `InferCNV.R` | Large-scale copy-number inference for malignant-cell identification. | infercnv |
| `Pseudotime analysis.R` | Trajectory inference and analysis of CPSig dynamics. | Monocle |
| `SCT.R` | Spatial-transcriptomics normalization, clustering, and signature mapping. | Seurat, AUCell |

#### 2. Bulk data processing and molecular stratification

| Script | Description | Main packages or methods |
|---|---|---|
| `Deconvolution.R` | Estimation of cellular abundance in bulk tissues. | CIBERSORTx, MuSiC |

#### 3. Functional and genomic characterization

| Script | Description | Main packages or methods |
|---|---|---|
| `Functional enrichment analysis.R` | GO, KEGG, and GSEA functional analyses. | clusterProfiler |
| `EMT.R` | Quantification of epithelial-mesenchymal-transition activity. | GSVA, ssGSEA |
| `Mutation analysis.R` | Somatic-mutation landscapes and tumor mutational burden. | maftools |

#### 4. Immune landscape and clinical translation

| Script | Description | Main packages or methods |
|---|---|---|
| `ESTIMATE.R` | ImmuneScore and StromalScore calculation. | ESTIMATE |
| `T cell exhaustion+check.R` | T-cell exhaustion and immune-checkpoint analyses. | ggplot2, ggpubr |
| `TIDE.R` | Immunotherapy-response and immune dysfunction/exclusion analyses. | TIDE |
| `Survival analysis.R` | Kaplan-Meier and Cox regression analyses. | survival, survminer |

#### 5. Visualization

| Script | Description | Main packages or methods |
|---|---|---|
| `UpSet_plot.R` | Pan-cancer CPSig intersection and specificity visualization. | UpSetR |

### R prerequisites

R version 4.0.0 or later is recommended. Representative dependencies include:

```r
install.packages(c(
  "Seurat", "Monocle", "survival", "survminer", "estimate",
  "maftools", "tidyverse", "ggplot2", "pheatmap", "ggpubr"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "limma", "ComplexHeatmap", "clusterProfiler", "infercnv", "GSVA"
))
```

## License

The ProgMap Python package is distributed under the MIT License. See [`progmap-python/LICENSE`](progmap-python/LICENSE).
