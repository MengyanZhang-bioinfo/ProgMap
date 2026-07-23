# ProgMap Python Tool for Linux Servers

ProgMap integrates expression and methylation matrices from normal tissue, Stage I, and Stage II/III disease into fold-specific MECor features and trains a fixed 2048-128 skip-connected neural network for three-class prediction and progression-associated feature discovery. This release is designed for headless Linux workstations, servers, containers, and Slurm clusters.

## Reference environment

- Ubuntu 22.04/24.04 x86_64 or a compatible Linux distribution
- Python 3.9-3.11; Python 3.11 recommended
- TensorFlow/Keras 2.14.0
- NumPy 1.26.4, pandas 2.2.3, SciPy 1.13.1
- scikit-learn 1.5.2, statsmodels 0.14.4, joblib 1.4.2

The archived models record TensorFlow/Keras 2.14.0 and the original bytecode records Python 3.9. This release therefore pins TensorFlow 2.14.0 and supports Python 3.9-3.11. At least 16 GB RAM is recommended; 32 GB is safer for approximately 13,000 input genes.

## Install

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

Alternatively:

```bash
conda env create -f environment-linux.yml
conda activate progmap-linux
```

## Input

Each cancer directory under `--data-root` must contain six gene-by-sample CSV files. The first column contains gene identifiers and the header contains sample identifiers.

| Group | Expression | Methylation |
|---|---|---|
| Normal | `en.csv` | `mn.csv` |
| Stage I | `exp1.csv` or `e1.csv` | `me1.csv` or `m1.csv` |
| Stage II/III | `exp2.csv` or `e2.csv` | `me2.csv` or `m2.csv` |

Linux filenames are case-sensitive. A directory named `GEO` is always excluded.

## Run

Validate inputs without training:

```bash
progmap --data-root /data/PANCANCER --output /results/input_check --dry-run --device cpu
```

Run all cancers and retain all ranked genes:

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

Select the top 100 genes with a rank-sum test:

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

Use `nohup` for a standalone server or submit `server/slurm_run.sh` on Slurm. `--device gpu` fails immediately if TensorFlow cannot see a GPU; `--device auto` uses a visible GPU and otherwise runs on CPU.

## Docker

CPU:

```bash
docker build -t progmap:0.2.0 .
docker run --rm \
  -v /absolute/path/PANCANCER:/data/PANCANCER:ro \
  -v /absolute/path/results:/results \
  progmap:0.2.0 \
  --data-root /data/PANCANCER --output /results/run1 --device cpu --threads 8
```

GPU, with a compatible NVIDIA driver and NVIDIA Container Toolkit:

```bash
docker build -f Dockerfile.gpu -t progmap:0.2.0-gpu .
docker run --rm --gpus all \
  -v /absolute/path/PANCANCER:/data/PANCANCER:ro \
  -v /absolute/path/results:/results \
  progmap:0.2.0-gpu \
  --data-root /data/PANCANCER --output /results/run1 --device gpu
```

## Leakage control

Every outer fold fits imputation, expression scaling, methylation scaling, per-gene Pearson correlations, MECor construction, and the neural network using only the corresponding training data. Early stopping uses a separate internal split of that outer training fold. Each outer test sample is transformed with parameters learned without that sample and contributes one out-of-fold prediction.

Audit files include `fold_N/training_fold_correlations.csv`, `fold_N/sample_roles.csv`, `cross_validated_predictions.csv`, and `run_config.json`.

## Tests

```bash
python -m pytest
python -m build
```

The Ubuntu workflow in `.github/workflows/linux-ci.yml` installs the pinned environment, runs unit and end-to-end tests, exercises the CLI, and builds the wheel and source distribution. Data are not distributed with this repository.

