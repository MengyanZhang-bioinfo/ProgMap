# ProgMap v0.1.0 test report

Date: 2026-07-25

## Reference environment

- Python 3.11.6
- TensorFlow/Keras 2.14.0
- NumPy 1.26.4
- pandas 2.2.3
- SciPy 1.13.1
- scikit-learn 1.5.2
- statsmodels 0.14.4
- joblib 1.4.2
- tqdm 4.67.1
- protobuf 4.25.9
- ml-dtypes 0.2.0

The dependency versions are pinned in the Linux requirements and Conda environment files.

## Automated tests

```bash
python -m pytest
```

Result: **10 passed**.

The tests cover input aliases and modality alignment, exclusion of incomplete dataset directories, training-partition-only imputation/scaling/correlation, Pearson and Spearman correlations, command-line defaults and overrides, the Dense(2048)-Dense(128) skip architecture, enhanced integrated gradients, t-test/Wilcoxon/permutation feature ranking, and an end-to-end nested-cross-validation pipeline.

## Installed-command tests

A generic synthetic dataset containing 27 paired samples and 12 genes was analyzed with the installed `progmap` command.

The nested test used two outer folds, two inner folds, two maximum epochs, fold-specific Pearson MECor features, t-tests, saved raw attributions, and saved outer-fold models. The result contained:

- out-of-fold predictions for all 27 samples;
- one saved `.keras` model per outer fold;
- four inner-fold epoch-selection histories;
- fold-specific preprocessors and correlations;
- attribution arrays for all three output classes; and
- selected-gene CSV files containing attribution, raw P value, adjusted P value, effect, and significance columns.

A second installed-command test used a custom seed, learning rate, fixed-epoch training, Spearman correlation, Wilcoxon testing, disabled raw-attribution output, and `--top-n 7`. All requested values were recorded in `run_config.json`, and exactly seven genes were written to `progression_genes.csv`.

## Packaging and Linux checks

- `python -m compileall` completed without syntax errors.
- `git diff --check` reported no whitespace errors.
- `python -m build` produced `progmap-0.1.0-py3-none-any.whl` and `progmap-0.1.0.tar.gz`.
- The source distribution includes the generic synthetic six-file example.
- The Ubuntu 22.04 GitHub Actions workflow installs the pinned environment, runs all tests, executes the installed command, and builds both distributions.

## Test boundary

Local CPU execution is verified. GPU execution requires a compatible Linux NVIDIA/TensorFlow environment and is guarded by an explicit `--device gpu` availability check.
