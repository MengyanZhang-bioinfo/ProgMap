# Changelog

## 0.3.0 - 2026-07-25

- Changed the article defaults to stratified three-fold outer cross-validation and three-fold nested epoch selection.
- Added configurable early-stopping strategy, monitor, patience, warm-up, minimum improvement, epoch limit, fold counts, learning rate, and seed.
- Added fold-specific Pearson or Spearman expression-methylation correlation.
- Added progress bars, pooled out-of-fold metrics, epoch-selection audit tables, and detailed progression-gene output.
- Changed the default selected output to all FDR-significant genes while retaining top-N and all-ranked options.
- Added a compact BRCA-format example dataset and expanded the English command-line documentation.

## 0.2.0 - 2026-07-23

- Added a Linux-first, headless server release.
- Added reproducible Ubuntu CPU requirements and a Conda environment.
- Added CPU and NVIDIA-container-compatible Dockerfiles.
- Added shell and Slurm launch scripts using POSIX paths.
- Added CPU thread controls and explicit GPU availability validation.
- Added runtime platform, Python, and device metadata to `run_config.json`.
- Added an Ubuntu 22.04 GitHub Actions workflow and a synthetic input generator.
- Removed Windows-only PowerShell launch and dependency files from the server package.

## 0.1.0 - 2026-07-22

- Initial packaged implementation of the fixed ProgMap 2048-128 skip model.
