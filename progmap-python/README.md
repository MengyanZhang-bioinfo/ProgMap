# ProgMap

Command-line software for fold-specific MECor feature construction, three-class stage classification, and progression-associated feature ranking from paired gene-expression and DNA-methylation matrices.

## Requirements

- Linux x86_64
- Python 3.9-3.11
- TensorFlow/Keras 2.14.0
- NumPy 1.26.4
- pandas 2.2.3
- SciPy 1.13.1
- scikit-learn 1.5.2
- statsmodels 0.14.4
- joblib 1.4.2

## Installation

### `venv`

```bash
git clone https://github.com/MengyanZhang-bioinfo/ProgMap.git
cd ProgMap/progmap-python
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-linux-cpu.txt
python -m pip install --no-deps .
progmap --help
```

### Conda

```bash
conda env create -f environment-linux.yml
conda activate progmap-linux
python -m pip install --no-deps .
progmap --help
```

## Input files

`--data-root` must contain one directory per cancer type. Each cancer directory must contain six gene-by-sample CSV matrices. The first column contains gene identifiers, and the header contains sample identifiers.

| Sample group | Expression | Methylation |
|---|---|---|
| Normal | `en.csv` | `mn.csv` |
| Stage I | `exp1.csv` or `e1.csv` | `me1.csv` or `m1.csv` |
| Stage II/III | `exp2.csv` or `e2.csv` | `me2.csv` or `m2.csv` |

Example:

```text
/data/PANCANCER/
├── BRCA/
│   ├── en.csv
│   ├── e1.csv
│   ├── e2.csv
│   ├── mn.csv
│   ├── m1.csv
│   └── m2.csv
└── LUAD/
    └── ...
```

Filenames are case-sensitive. A directory named `GEO` is excluded.

## Usage

Validate the input files without training:

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/input_check \
  --dry-run \
  --device cpu
```

Analyze all cancer directories with the default t-test and retain all ranked genes:

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap \
  --cancers all \
  --test ttest \
  --top-n all \
  --device auto \
  --threads 8
```

Select the top 100 genes with the Wilcoxon rank-sum test:

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap_wilcoxon \
  --cancers all \
  --test wilcoxon \
  --top-n 100 \
  --device auto \
  --threads 8
```

Use a permutation test and 1,000 feature-effect bootstrap iterations:

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap_permutation \
  --cancers BRCA,LUAD \
  --test permutation \
  --permutations 1000 \
  --bootstrap-iterations 1000 \
  --top-n 100 \
  --device auto
```

Run `progmap --help` for all command-line options.

## Statistical options

- `--test ttest`: one-sided Welch's two-sample t-test.
- `--test wilcoxon`: one-sided Wilcoxon rank-sum/Mann-Whitney test.
- `--test permutation`: one-sided permutation test of the mean difference.
- `--correction`: multiple-testing correction (`bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, or `fdr_tsbky`).
- `--alpha`: adjusted significance threshold.
- `--top-n all`: retain all ranked genes.
- `--top-n N`: retain the top `N` unique genes.
- `--bootstrap-iterations N`: calculate bootstrap confidence intervals for feature effects; `0` disables this option.

## Cross-validation and preprocessing

The outer cross-validation loop is stratified. Within each outer fold, the software fits imputation, expression min-max scaling, methylation min-max scaling, per-gene Pearson correlations, MECor construction, and the neural network using only the training portion of that fold. Early stopping uses an internal split of the outer-training data. The held-out outer-fold samples are transformed with training-fold parameters and are used only for out-of-fold evaluation and attribution.

The neural network has fixed hidden layers of 2,048 and 128 units with a skip connection to the three-class output layer.

## Output files

The output root contains `run_config.json` and `run_summary.json`. Each cancer-specific directory contains:

- `data_summary.json`
- `cross_validated_predictions.csv`
- `fold_metrics.csv`
- `features_by_class_all.csv`
- `features_ranked_genes.csv`
- `features_selected.csv`
- `summary.json`
- one `fold_N` directory per outer fold

Each `fold_N` directory contains the fitted preprocessor, training-fold correlation values, sample roles, training history, and the saved model unless `--no-save-models` is specified.

## Server execution

Run in the background:

```bash
nohup progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap \
  --cancers all \
  --test ttest \
  --top-n all \
  --device auto \
  --threads 8 \
  > progmap.log 2>&1 &
```

Run the supplied shell wrapper:

```bash
bash examples/run_all.sh /data/PANCANCER /results/progmap
```

Submit the supplied Slurm script:

```bash
sbatch server/slurm_run.sh /data/PANCANCER /results/progmap
```

`--device gpu` stops with an error if TensorFlow cannot detect a GPU. `--device auto` uses an available GPU and otherwise runs on CPU.

## Docker

CPU image:

```bash
docker build -t progmap:0.2.0 .
docker run --rm \
  -v /absolute/path/PANCANCER:/data/PANCANCER:ro \
  -v /absolute/path/results:/results \
  progmap:0.2.0 \
  --data-root /data/PANCANCER \
  --output /results/run1 \
  --device cpu \
  --threads 8
```

GPU image:

```bash
docker build -f Dockerfile.gpu -t progmap:0.2.0-gpu .
docker run --rm --gpus all \
  -v /absolute/path/PANCANCER:/data/PANCANCER:ro \
  -v /absolute/path/results:/results \
  progmap:0.2.0-gpu \
  --data-root /data/PANCANCER \
  --output /results/run1 \
  --device gpu
```

GPU execution requires a compatible NVIDIA driver and NVIDIA Container Toolkit.

## Tests

```bash
python -m pytest
python -m build
```

## License

See [LICENSE](LICENSE).
