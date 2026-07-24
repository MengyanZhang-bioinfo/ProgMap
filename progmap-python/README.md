# ProgMap

ProgMap is a command-line package for three-class molecular-stage modeling and progression-feature attribution from paired gene-expression and DNA-methylation matrices.

## Installation

### Linux virtual environment

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

### Conda

```bash
conda env create -f environment-linux.yml
conda activate progmap-linux
progmap --version
```

The reference environment uses Python 3.11, TensorFlow/Keras 2.14.0, NumPy 1.26.4, pandas 2.2.3, SciPy 1.13.1, scikit-learn 1.5.2, statsmodels 0.14.4, joblib 1.4.2, and tqdm 4.67.1.

## Input data

`--data-root` contains one directory per cancer type. Each cancer directory contains six gene-by-sample CSV matrices. The first column contains gene identifiers, and the remaining columns are samples.

| Class | Expression | DNA methylation |
|---|---|---|
| Normal | `en.csv` | `mn.csv` |
| Stage I | `e1.csv` or `exp1.csv` | `m1.csv` or `me1.csv` |
| Stage II/III | `e2.csv` or `exp2.csv` | `m2.csv` or `me2.csv` |

```text
/data/PANCANCER/
├── BRCA/
│   ├── en.csv
│   ├── e1.csv
│   ├── e2.csv
│   ├── mn.csv
│   ├── m1.csv
│   └── m2.csv
└── another_cancer/
    └── ...
```

Filenames are case-sensitive. A directory named `GEO` is ignored. Expression and methylation matrices are aligned by common gene and sample identifiers before analysis.

## One-command analysis

Run every complete cancer directory with the manuscript defaults:

```bash
progmap --data-root /data/PANCANCER --output /results/progmap
```

The command reads the raw expression and methylation CSV files, performs fold-specific preprocessing, trains the models, calculates held-out attributions, tests the features, and writes all results under `--output`. Progress bars are enabled by default; use `--no-progress` or `--quiet` to disable them.

Validate input files without training:

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/input_check \
  --dry-run \
  --device cpu
```

## Default analysis settings

| Setting | Default |
|---|---|
| Classes | normal, stage I, stage II/III |
| Outer cross-validation | stratified 3-fold |
| Epoch selection | nested stratified 3-fold within each outer-training fold |
| Maximum epochs | 1,000 |
| Early-stopping monitor | class-balanced validation loss |
| Warm-up / patience / minimum change | 20 / 50 / `1e-4` |
| Selected outer-fold epoch | median best epoch across the inner folds |
| Model | Dense(2048)-Dropout-Dense(128)-Dropout with an input skip connection to the softmax output |
| Optimizer | Adam with exponentially decayed learning rate |
| Initial learning rate / decay | `1e-3` / 0.9 every 10,000 steps |
| Batch size / dropout | 16 / 0.1 |
| Training class weights | 0.25, 0.50, 0.25 |
| Random seed | 42 |
| Missing values | 10-nearest-neighbor imputation |
| Expression and methylation scaling | independently scaled to [0, 1] within each training partition |
| MECor correlation | per-gene Pearson correlation within each training partition |
| Attribution | enhanced integrated gradients, 64 steps and 3 baselines |
| Feature test | one-sided Welch t-test |
| Multiple-testing correction | Benjamini-Hochberg FDR at 0.01 |
| Selected output | all significant unique genes |

For each outer fold, all imputation, min-max scaling, expression-methylation correlation, MECor construction, and epoch selection exclude the outer test fold. Inner cross-validation selects an epoch count. A new model is then trained from scratch on the complete outer-training fold for that fixed number of epochs, and the outer test fold is evaluated once.

## Optional parameters

All requested training and feature-selection settings are available from the command line:

```bash
progmap \
  --data-root /data/PANCANCER \
  --output /results/custom_run \
  --cancers BRCA,COAD \
  --folds 5 \
  --inner-folds 3 \
  --epochs 600 \
  --learning-rate 0.0005 \
  --early-stopping nested \
  --early-stopping-monitor val_loss \
  --warmup-epochs 20 \
  --patience 40 \
  --min-delta 0.0001 \
  --seed 2026 \
  --correlation-method spearman \
  --test wilcoxon \
  --correction fdr_bh \
  --alpha 0.01 \
  --top-n 100
```

Important choices include:

- `--early-stopping nested`: nested cross-validation followed by full outer-training-fold refitting.
- `--early-stopping holdout`: one stratified internal holdout controlled by `--inner-validation-fraction`.
- `--early-stopping off`: no epoch selection; every outer model uses `--epochs`.
- `--early-stopping-monitor val_loss` or `val_auc`.
- `--correlation-method pearson` or `spearman`.
- `--test ttest`, `wilcoxon`, or `permutation`.
- `--top-n significant`: output every gene passing the adjusted-P threshold.
- `--top-n N`: output the top `N` unique ranked genes.
- `--top-n all`: output every ranked gene.
- `--bootstrap-iterations N`: add feature-effect bootstrap confidence intervals; `0` disables them.
- `--no-save-models` and `--no-save-attributions`: reduce output size.

Run `progmap --help` for every available option.

## BRCA-format example

The repository includes a small synthetic `BRCA_DEMO` dataset with the same six-file layout as a real BRCA analysis. It contains no patient data.

Run a short end-to-end example:

```bash
progmap \
  --data-root examples/data/PANCANCER \
  --cancers BRCA_DEMO \
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

To create another synthetic dataset:

```bash
python examples/create_synthetic_data.py \
  --output /tmp/PANCANCER \
  --cancer BRCA_DEMO \
  --genes 12 \
  --samples-per-class 9 \
  --seed 7
```

## Output files

```text
OUTPUT/
├── run_config.json
├── run_summary.json
└── BRCA/
    ├── data_summary.json
    ├── cross_validated_predictions.csv
    ├── fold_metrics.csv
    ├── pooled_out_of_fold_metrics.json
    ├── epoch_selection_all_folds.csv
    ├── features_by_class_all.csv
    ├── features_ranked_genes.csv
    ├── features_selected.csv
    ├── progression_genes.csv
    ├── attributions_class_0.npz
    ├── attributions_class_1.npz
    ├── attributions_class_2.npz
    └── fold_1/
        ├── model.keras
        ├── preprocessor.joblib
        ├── training_fold_correlations.csv
        ├── training_history.csv
        ├── epoch_selection.csv
        ├── sample_roles.csv
        └── inner_fold_1/
            ├── training_fold_correlations.csv
            ├── training_history.csv
            └── sample_roles.csv
```

`progression_genes.csv` and `features_selected.csv` contain the selected genes and include the best-associated class, class name, target and background median absolute attribution, raw P value, adjusted P value, attribution effect, significance flag, test, and correction method. `features_by_class_all.csv` retains class-specific statistics for every input gene. A separate `.keras` model is saved for every outer fold.

## Server execution

```bash
nohup progmap \
  --data-root /data/PANCANCER \
  --output /results/progmap \
  --device auto \
  --threads 8 \
  > progmap.log 2>&1 &
```

Shell wrapper:

```bash
bash examples/run_all.sh /data/PANCANCER /results/progmap
```

Slurm:

```bash
sbatch server/slurm_run.sh /data/PANCANCER /results/progmap
```

`--device gpu` stops if TensorFlow cannot detect a GPU. `--device auto` uses an available GPU and otherwise runs on CPU.

## Docker

```bash
docker build -t progmap:0.3.0 .
docker run --rm \
  -v /absolute/path/PANCANCER:/data/PANCANCER:ro \
  -v /absolute/path/results:/results \
  progmap:0.3.0 \
  --data-root /data/PANCANCER \
  --output /results/run1 \
  --device cpu
```

For an NVIDIA host, build `Dockerfile.gpu` and run the image with `--gpus all` and `--device gpu`.

## Tests

```bash
python -m pytest
python -m build
```

## License

See [LICENSE](LICENSE).
