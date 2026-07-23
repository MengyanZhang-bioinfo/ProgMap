# ProgMap 0.2.0 server-release test report

Date: 2026-07-23

## Reference dependency set

- Python 3.11.6
- TensorFlow/Keras 2.14.0
- NumPy 1.26.4
- pandas 2.2.3
- SciPy 1.13.1
- scikit-learn 1.5.2
- statsmodels 0.14.4
- joblib 1.4.2
- protobuf 4.25.9
- ml-dtypes 0.2.0

The pinned dependencies in `requirements-linux-cpu.txt` installed successfully in a new Python 3.11 virtual environment, and `pip check` reported no broken requirements.

## Automated tests

Command:

```bash
python -m pytest
```

Result: **7 passed**.

Coverage includes input aliases and expression/methylation pairing, exclusion of the GEO directory, training-only imputation/scaling/correlation estimation, the fixed Dense(2048)-Dense(128)-skip architecture, enhanced integrated gradients, t-test/Wilcoxon/permutation feature ranking, and an end-to-end TensorFlow training pipeline.

## Installed-command smoke test

The built `progmap-0.2.0-py3-none-any.whl` was force-installed instead of the editable source tree. The installed console command reported `progmap 0.2.0` from `site-packages`.

A separate synthetic dataset containing 27 paired samples and 12 genes was generated. The installed command completed:

- input-only validation;
- two outer folds with one training epoch per fold;
- training-fold-only MECor preprocessing;
- held-out predictions for all 27 samples;
- Wilcoxon feature testing; and
- selection of the requested top five genes.

The resulting `run_config.json` correctly recorded Python and package versions, requested CPU execution, two CPU threads, zero visible GPUs, and effective CPU execution.

## Packaging and shell checks

- `python -m build` produced `progmap-0.2.0-py3-none-any.whl` and `progmap-0.2.0.tar.gz`.
- Git Bash `bash -n` accepted `examples/run_all.sh` and `server/slurm_run.sh`.
- The Ubuntu 22.04 workflow at `.github/workflows/linux-ci.yml` performs the same installation, test, CLI smoke-test, and build sequence after publication to GitHub.

## Test boundary

The local host used for the above verification was not a Linux kernel and did not provide Docker, WSL, or a GPU. Therefore the local result verifies the Linux-pinned dependency set, platform-neutral package, wheel, CLI, shell syntax, and complete analysis path, but it is not labeled as a native-Linux execution. Native Ubuntu status must be taken from the GitHub Actions run generated after the commit is pushed. GPU execution likewise requires a compatible NVIDIA runtime and is checked at startup with `--device gpu`.

