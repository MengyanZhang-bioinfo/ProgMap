from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version
import json
import math
import platform
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from tqdm.auto import tqdm

from .attribution import (
    balanced_attribution_indices,
    integrated_gradients,
    select_median_baselines,
)
from .data import CLASS_NAMES, MultiOmicsData, load_cancer
from .model import (
    EARLY_STOPPING_MONITORS,
    EARLY_STOPPING_STRATEGIES,
    build_skip_model,
    configure_runtime,
    fit_fixed_epochs,
    fit_model,
    require_tensorflow,
    set_global_seed,
)
from .preprocessing import CORRELATION_METHODS, FoldPreprocessor
from .statistics import collapse_gene_ranking, rank_features


@dataclass(frozen=True)
class RunConfig:
    data_root: str
    output: str
    cancers: tuple[str, ...]
    folds: int = 3
    inner_folds: int = 3
    early_stopping: str = "nested"
    early_stopping_monitor: str = "val_loss"
    inner_validation_fraction: float = 0.15
    epochs: int = 1000
    patience: int = 50
    min_delta: float = 1e-4
    warmup_epochs: int = 20
    balance_validation_loss: bool = True
    batch_size: int = 16
    dropout: float = 0.1
    learning_rate: float = 1e-3
    decay_steps: int = 10_000
    decay_rate: float = 0.9
    class_weights: tuple[float, float, float] = (0.25, 0.5, 0.25)
    imputation: str = "knn"
    knn_neighbors: int = 10
    correlation_method: str = "pearson"
    ig_steps: int = 64
    ig_baselines: int = 3
    attribution_batch_size: int = 32
    max_attribution_samples_per_class: int | None = None
    test: str = "ttest"
    correction: str = "fdr_bh"
    alpha: float = 0.01
    permutations: int = 1000
    bootstrap_iterations: int = 0
    top_n: int | str = "significant"
    seed: int = 42
    save_attributions: bool = True
    save_models: bool = True
    overwrite: bool = False
    dry_run: bool = False
    verbose: int = 1
    progress: bool = True
    device: str = "auto"
    threads: int | None = None


def _versions() -> dict[str, str]:
    packages = [
        "progmap",
        "tensorflow",
        "numpy",
        "pandas",
        "scikit-learn",
        "scipy",
        "statsmodels",
        "tqdm",
    ]
    result = {}
    for package in packages:
        try:
            result[package] = version(package)
        except PackageNotFoundError:
            result[package] = "not-installed-as-distribution"
    result["python"] = platform.python_version()
    result["platform"] = platform.platform()
    result["executable"] = sys.executable
    return result


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def _prepare_output(path: Path, overwrite: bool) -> None:
    if path.exists() and any(path.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"Output folder is not empty: {path}. Use --overwrite to replace this cancer result."
            )
        resolved = path.resolve()
        if resolved.parent == resolved or len(resolved.parts) < 3:
            raise RuntimeError(f"Refusing to remove unsafe output path: {resolved}")
        shutil.rmtree(resolved)
    path.mkdir(parents=True, exist_ok=True)


def _inner_split(
    outer_train_indices: np.ndarray,
    labels: np.ndarray,
    fraction: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    if not 0.0 < fraction < 0.5:
        raise ValueError("inner_validation_fraction must be between 0 and 0.5")
    outer_labels = labels[outer_train_indices]
    class_count = len(np.unique(outer_labels))
    validation_size = max(class_count, int(math.ceil(len(outer_train_indices) * fraction)))
    validation_size = min(validation_size, len(outer_train_indices) - class_count)
    if validation_size < class_count:
        raise ValueError("Not enough outer-training samples for stratified internal validation")
    fit_indices, validation_indices = train_test_split(
        outer_train_indices,
        test_size=validation_size,
        random_state=seed,
        shuffle=True,
        stratify=outer_labels,
    )
    return np.asarray(fit_indices), np.asarray(validation_indices)


def _usable_folds(requested: int, labels: np.ndarray, context: str) -> int:
    counts = np.bincount(np.asarray(labels, dtype=np.int64), minlength=3)
    minimum = int(counts.min())
    folds = min(requested, minimum)
    if folds < 2:
        raise ValueError(
            f"{context} needs at least two samples in every class; class counts are {counts.tolist()}"
        )
    return folds


def _fold_metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predictions = np.argmax(probabilities, axis=1)
    one_hot = np.eye(3, dtype=np.float32)[labels]
    metrics = {
        "accuracy": float(accuracy_score(labels, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predictions)),
        "auc_micro": float(
            roc_auc_score(one_hot, probabilities, average="micro", multi_class="ovr")
        ),
        "auc_macro": float(
            roc_auc_score(one_hot, probabilities, average="macro", multi_class="ovr")
        ),
    }
    for label in (0, 1, 2):
        metrics[f"auc_class_{label}"] = float(
            roc_auc_score(one_hot[:, label], probabilities[:, label])
        )
    return metrics


def _save_correlations(
    path: Path,
    genes: np.ndarray,
    correlations: np.ndarray,
    method: str,
) -> None:
    pd.DataFrame(
        {
            "gene": genes,
            f"training_fold_{method}_correlation": correlations,
        }
    ).to_csv(path, index=False)


def _preprocessor(config: RunConfig) -> FoldPreprocessor:
    return FoldPreprocessor(
        imputation=config.imputation,
        knn_neighbors=config.knn_neighbors,
        correlation_method=config.correlation_method,
    )


def _model(config: RunConfig, n_genes: int):
    return build_skip_model(
        n_genes,
        dropout=config.dropout,
        learning_rate=config.learning_rate,
        decay_steps=config.decay_steps,
        decay_rate=config.decay_rate,
    )


def _class_weight(config: RunConfig) -> dict[int, float]:
    return {index: weight for index, weight in enumerate(config.class_weights)}


def _progress_description(config: RunConfig, value: str) -> str | None:
    return value if config.progress else None


def _save_sample_roles(
    path: Path,
    data: MultiOmicsData,
    roles: list[tuple[np.ndarray, str]],
) -> None:
    choices = [np.isin(np.arange(data.n_samples), indices) for indices, _ in roles]
    labels = [label for _, label in roles]
    pd.DataFrame(
        {
            "sample_id": data.sample_ids,
            "role": np.select(choices, labels, default="unused"),
            "label": data.labels,
            "class_name": [CLASS_NAMES[int(value)] for value in data.labels],
        }
    ).to_csv(path, index=False)


def _select_epochs_nested(
    config: RunConfig,
    data: MultiOmicsData,
    outer_train: np.ndarray,
    fold: int,
    fold_output: Path,
) -> tuple[int, pd.DataFrame]:
    outer_labels = data.labels[outer_train]
    inner_folds = _usable_folds(
        config.inner_folds,
        outer_labels,
        f"{data.cancer} outer fold {fold} inner cross-validation",
    )
    splitter = StratifiedKFold(
        n_splits=inner_folds,
        shuffle=True,
        random_state=config.seed + fold * 100,
    )
    rows: list[dict[str, Any]] = []
    iterator = splitter.split(outer_train, outer_labels)
    iterator = tqdm(
        iterator,
        total=inner_folds,
        desc=f"{data.cancer} fold {fold} inner CV",
        unit="fold",
        leave=False,
        disable=not config.progress,
        dynamic_ncols=True,
    )
    for inner_fold, (relative_fit, relative_validation) in enumerate(iterator, start=1):
        fit_indices = outer_train[relative_fit]
        validation_indices = outer_train[relative_validation]
        inner_seed = config.seed + fold * 100 + inner_fold
        set_global_seed(inner_seed)
        inner_output = fold_output / f"inner_fold_{inner_fold}"
        inner_output.mkdir(parents=True, exist_ok=True)

        processor = _preprocessor(config).fit(
            data.expression[fit_indices],
            data.methylation[fit_indices],
            data.genes,
        )
        x_fit = processor.transform(
            data.expression[fit_indices], data.methylation[fit_indices]
        )
        x_validation = processor.transform(
            data.expression[validation_indices], data.methylation[validation_indices]
        )
        y_fit = data.labels[fit_indices]
        y_validation = data.labels[validation_indices]

        _save_correlations(
            inner_output / "training_fold_correlations.csv",
            data.genes,
            processor.correlations,
            config.correlation_method,
        )
        _save_sample_roles(
            inner_output / "sample_roles.csv",
            data,
            [
                (fit_indices, "inner_fit"),
                (validation_indices, "inner_validation"),
            ],
        )
        model = _model(config, data.n_genes)
        result = fit_model(
            model,
            x_fit,
            y_fit,
            x_validation,
            y_validation,
            epochs=config.epochs,
            patience=config.patience,
            batch_size=config.batch_size,
            class_weight=_class_weight(config),
            monitor=config.early_stopping_monitor,
            min_delta=config.min_delta,
            warmup_epochs=config.warmup_epochs,
            balance_validation_loss=config.balance_validation_loss,
            progress_description=_progress_description(
                config, f"{data.cancer} O{fold}/I{inner_fold}"
            ),
            verbose=config.verbose,
        )
        pd.DataFrame(result.history).to_csv(
            inner_output / "training_history.csv", index=False
        )
        rows.append(
            {
                "outer_fold": fold,
                "inner_fold": inner_fold,
                "seed": inner_seed,
                "fit_samples": len(fit_indices),
                "validation_samples": len(validation_indices),
                "best_epoch": result.best_epoch,
                "stopped_epoch": result.stopped_epoch,
                "monitor": config.early_stopping_monitor,
                "best_monitor_value": result.best_monitor_value,
            }
        )
        require_tensorflow().keras.backend.clear_session()

    selection = pd.DataFrame(rows)
    selected_epoch = max(1, int(np.rint(np.median(selection["best_epoch"]))))
    return selected_epoch, selection


def _select_epochs_holdout(
    config: RunConfig,
    data: MultiOmicsData,
    outer_train: np.ndarray,
    fold: int,
    fold_output: Path,
) -> tuple[int, pd.DataFrame]:
    fold_seed = config.seed + fold * 100
    fit_indices, validation_indices = _inner_split(
        outer_train,
        data.labels,
        config.inner_validation_fraction,
        fold_seed,
    )
    set_global_seed(fold_seed)
    inner_output = fold_output / "holdout_selection"
    inner_output.mkdir(parents=True, exist_ok=True)
    processor = _preprocessor(config).fit(
        data.expression[fit_indices], data.methylation[fit_indices], data.genes
    )
    x_fit = processor.transform(data.expression[fit_indices], data.methylation[fit_indices])
    x_validation = processor.transform(
        data.expression[validation_indices], data.methylation[validation_indices]
    )
    y_fit = data.labels[fit_indices]
    y_validation = data.labels[validation_indices]
    _save_correlations(
        inner_output / "training_fold_correlations.csv",
        data.genes,
        processor.correlations,
        config.correlation_method,
    )
    _save_sample_roles(
        inner_output / "sample_roles.csv",
        data,
        [(fit_indices, "holdout_fit"), (validation_indices, "holdout_validation")],
    )
    model = _model(config, data.n_genes)
    result = fit_model(
        model,
        x_fit,
        y_fit,
        x_validation,
        y_validation,
        epochs=config.epochs,
        patience=config.patience,
        batch_size=config.batch_size,
        class_weight=_class_weight(config),
        monitor=config.early_stopping_monitor,
        min_delta=config.min_delta,
        warmup_epochs=config.warmup_epochs,
        balance_validation_loss=config.balance_validation_loss,
        progress_description=_progress_description(config, f"{data.cancer} O{fold} holdout"),
        verbose=config.verbose,
    )
    pd.DataFrame(result.history).to_csv(inner_output / "training_history.csv", index=False)
    require_tensorflow().keras.backend.clear_session()
    selection = pd.DataFrame(
        [
            {
                "outer_fold": fold,
                "inner_fold": 1,
                "seed": fold_seed,
                "fit_samples": len(fit_indices),
                "validation_samples": len(validation_indices),
                "best_epoch": result.best_epoch,
                "stopped_epoch": result.stopped_epoch,
                "monitor": config.early_stopping_monitor,
                "best_monitor_value": result.best_monitor_value,
            }
        ]
    )
    return result.best_epoch, selection


def run_cancer(config: RunConfig, cancer: str) -> dict[str, Any]:
    if config.early_stopping not in EARLY_STOPPING_STRATEGIES:
        raise ValueError(f"early_stopping must be one of {EARLY_STOPPING_STRATEGIES}")
    if config.early_stopping_monitor not in EARLY_STOPPING_MONITORS:
        raise ValueError(
            f"early_stopping_monitor must be one of {EARLY_STOPPING_MONITORS}"
        )
    if config.correlation_method not in CORRELATION_METHODS:
        raise ValueError(f"correlation_method must be one of {CORRELATION_METHODS}")

    data = load_cancer(config.data_root, cancer)
    cancer_output = Path(config.output) / cancer
    _prepare_output(cancer_output, config.overwrite)
    class_counts = data.class_counts()
    folds = _usable_folds(config.folds, data.labels, cancer)
    sparse_class_warning = min(class_counts.values()) < 10

    metadata = {
        "cancer": cancer,
        "samples": data.n_samples,
        "genes": data.n_genes,
        "class_counts": class_counts,
        "requested_folds": config.folds,
        "used_folds": folds,
        "sparse_class_warning": sparse_class_warning,
        "file_aliases": {
            "stage_I_expression": "exp1.csv or e1.csv",
            "stage_II_III_expression": "exp2.csv or e2.csv",
            "stage_I_methylation": "me1.csv or m1.csv",
            "stage_II_III_methylation": "me2.csv or m2.csv",
        },
    }
    _write_json(cancer_output / "data_summary.json", metadata)
    if config.dry_run:
        return metadata

    set_global_seed(config.seed)
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=config.seed)
    prediction_frames: list[pd.DataFrame] = []
    metric_rows: list[dict[str, Any]] = []
    epoch_frames: list[pd.DataFrame] = []
    attribution_values: dict[int, list[np.ndarray]] = {0: [], 1: [], 2: []}
    attribution_labels: dict[int, list[np.ndarray]] = {0: [], 1: [], 2: []}
    attribution_samples: dict[int, list[np.ndarray]] = {0: [], 1: [], 2: []}

    outer_iterator = splitter.split(data.expression, data.labels)
    outer_iterator = tqdm(
        outer_iterator,
        total=folds,
        desc=f"{cancer} outer CV",
        unit="fold",
        disable=not config.progress,
        dynamic_ncols=True,
    )
    for fold, (outer_train, outer_test) in enumerate(outer_iterator, start=1):
        outer_train = np.asarray(outer_train)
        outer_test = np.asarray(outer_test)
        fold_seed = config.seed + fold
        fold_output = cancer_output / f"fold_{fold}"
        fold_output.mkdir(parents=True, exist_ok=True)

        if config.early_stopping == "nested":
            selected_epoch, epoch_selection = _select_epochs_nested(
                config, data, outer_train, fold, fold_output
            )
        elif config.early_stopping == "holdout":
            selected_epoch, epoch_selection = _select_epochs_holdout(
                config, data, outer_train, fold, fold_output
            )
        else:
            selected_epoch = config.epochs
            epoch_selection = pd.DataFrame(
                [
                    {
                        "outer_fold": fold,
                        "inner_fold": 0,
                        "seed": fold_seed,
                        "fit_samples": len(outer_train),
                        "validation_samples": 0,
                        "best_epoch": selected_epoch,
                        "stopped_epoch": selected_epoch,
                        "monitor": "off",
                        "best_monitor_value": np.nan,
                    }
                ]
            )
        epoch_selection["selected_outer_training_epoch"] = selected_epoch
        epoch_selection.to_csv(fold_output / "epoch_selection.csv", index=False)
        epoch_frames.append(epoch_selection)

        set_global_seed(fold_seed + 10_000)
        processor = _preprocessor(config).fit(
            data.expression[outer_train], data.methylation[outer_train], data.genes
        )
        x_train = processor.transform(
            data.expression[outer_train], data.methylation[outer_train]
        )
        x_test = processor.transform(data.expression[outer_test], data.methylation[outer_test])
        y_train = data.labels[outer_train]
        y_test = data.labels[outer_test]

        processor.save(fold_output / "preprocessor.joblib")
        _save_correlations(
            fold_output / "training_fold_correlations.csv",
            data.genes,
            processor.correlations,
            config.correlation_method,
        )
        _save_sample_roles(
            fold_output / "sample_roles.csv",
            data,
            [(outer_train, "outer_train"), (outer_test, "outer_test")],
        )

        model = _model(config, data.n_genes)
        history = fit_fixed_epochs(
            model,
            x_train,
            y_train,
            epochs=selected_epoch,
            batch_size=config.batch_size,
            class_weight=_class_weight(config),
            progress_description=_progress_description(config, f"{cancer} O{fold} final"),
            verbose=config.verbose,
        )
        pd.DataFrame(history).to_csv(fold_output / "training_history.csv", index=False)
        if config.save_models:
            model.save(fold_output / "model.keras")

        probabilities = model.predict(x_test, batch_size=config.batch_size, verbose=0)
        fold_metrics = _fold_metrics(y_test, probabilities)
        fold_metrics["fold"] = fold
        fold_metrics["selected_epoch"] = selected_epoch
        fold_metrics["outer_train_samples"] = len(outer_train)
        fold_metrics["outer_test_samples"] = len(outer_test)
        metric_rows.append(fold_metrics)
        prediction_frames.append(
            pd.DataFrame(
                {
                    "sample_id": data.sample_ids[outer_test],
                    "fold": fold,
                    "true_label": y_test,
                    "predicted_label": np.argmax(probabilities, axis=1),
                    "probability_normal": probabilities[:, 0],
                    "probability_stage_I": probabilities[:, 1],
                    "probability_stage_II_III": probabilities[:, 2],
                }
            )
        )

        rng = np.random.default_rng(fold_seed)
        attribution_index = balanced_attribution_indices(
            y_test, config.max_attribution_samples_per_class, rng
        )
        x_attribution = x_test[attribution_index]
        y_attribution = y_test[attribution_index]
        sample_attribution = data.sample_ids[outer_test][attribution_index]
        target_iterator = tqdm(
            (0, 1, 2),
            total=3,
            desc=f"{cancer} fold {fold} attribution",
            unit="class",
            leave=False,
            disable=not config.progress,
            dynamic_ncols=True,
        )
        for target_class in target_iterator:
            baselines = select_median_baselines(
                x_train, y_train, target_class, config.ig_baselines
            )
            values = integrated_gradients(
                model,
                x_attribution,
                baselines,
                target_class,
                steps=config.ig_steps,
                batch_size=config.attribution_batch_size,
            )
            attribution_values[target_class].append(values)
            attribution_labels[target_class].append(y_attribution.copy())
            attribution_samples[target_class].append(sample_attribution.copy())

        require_tensorflow().keras.backend.clear_session()

    predictions = pd.concat(prediction_frames, ignore_index=True)
    metrics = pd.DataFrame(metric_rows).sort_values("fold")
    epoch_table = pd.concat(epoch_frames, ignore_index=True)
    predictions.to_csv(cancer_output / "cross_validated_predictions.csv", index=False)
    metrics.to_csv(cancer_output / "fold_metrics.csv", index=False)
    epoch_table.to_csv(cancer_output / "epoch_selection_all_folds.csv", index=False)
    probability_columns = [
        "probability_normal",
        "probability_stage_I",
        "probability_stage_II_III",
    ]
    pooled_metrics = _fold_metrics(
        predictions["true_label"].to_numpy(dtype=np.int64),
        predictions[probability_columns].to_numpy(dtype=np.float32),
    )
    _write_json(cancer_output / "pooled_out_of_fold_metrics.json", pooled_metrics)

    feature_frames = []
    for target_class in tqdm(
        (0, 1, 2),
        total=3,
        desc=f"{cancer} feature statistics",
        unit="class",
        disable=not config.progress,
        dynamic_ncols=True,
    ):
        values = np.concatenate(attribution_values[target_class], axis=0)
        labels = np.concatenate(attribution_labels[target_class])
        samples = np.concatenate(attribution_samples[target_class])
        ranked = rank_features(
            values,
            labels,
            data.genes,
            target_class,
            test=config.test,
            correction=config.correction,
            alpha=config.alpha,
            permutations=config.permutations,
            bootstrap_iterations=config.bootstrap_iterations,
            seed=config.seed,
            progress=config.progress,
            progress_description=f"{cancer} class {target_class}",
        )
        ranked.insert(1, "class_name", CLASS_NAMES[target_class])
        feature_frames.append(ranked)
        if config.save_attributions:
            np.savez_compressed(
                cancer_output / f"attributions_class_{target_class}.npz",
                attributions=values,
                true_labels=labels,
                sample_ids=samples,
                genes=data.genes,
            )

    all_features = pd.concat(feature_frames, ignore_index=True)
    all_features.to_csv(cancer_output / "features_by_class_all.csv", index=False)
    gene_ranking = collapse_gene_ranking(all_features)
    gene_ranking.to_csv(cancer_output / "features_ranked_genes.csv", index=False)
    if config.top_n == "significant":
        selected = gene_ranking.loc[gene_ranking["significant_in_best_class"]].copy()
        selection_mode = "all significant genes"
    elif config.top_n == "all":
        selected = gene_ranking.copy()
        selection_mode = "all ranked genes"
    else:
        selected = gene_ranking.head(int(config.top_n)).copy()
        selection_mode = f"top {int(config.top_n)} ranked genes"
    selected.to_csv(cancer_output / "features_selected.csv", index=False)
    selected.to_csv(cancer_output / "progression_genes.csv", index=False)

    summary = {
        **metadata,
        "mean_fold_metrics": {
            column: float(metrics[column].mean())
            for column in metrics.columns
            if column not in {"fold"} and pd.api.types.is_numeric_dtype(metrics[column])
        },
        "pooled_out_of_fold_metrics": pooled_metrics,
        "early_stopping_strategy": config.early_stopping,
        "early_stopping_monitor": config.early_stopping_monitor,
        "correlation_method": config.correlation_method,
        "feature_test": config.test,
        "multiple_testing_correction": config.correction,
        "alpha": config.alpha,
        "selection_mode": selection_mode,
        "selected_feature_count": int(len(selected)),
        "significant_unique_genes": int(gene_ranking["significant_in_best_class"].sum()),
    }
    _write_json(cancer_output / "summary.json", summary)
    return summary


def run(config: RunConfig) -> list[dict[str, Any]]:
    output = Path(config.output)
    output.mkdir(parents=True, exist_ok=True)
    runtime_device = (
        {"requested": config.device, "effective": "not-initialized-dry-run"}
        if config.dry_run
        else configure_runtime(config.device, config.threads)
    )
    run_metadata = {
        "config": asdict(config),
        "versions": _versions(),
        "runtime_device": runtime_device,
        "architecture": {
            "hidden_layers": [2048, 128],
            "skip_connection": "concatenate input MECor vector with second hidden layer",
            "output_classes": [CLASS_NAMES[i] for i in (0, 1, 2)],
        },
        "preprocessing": (
            "Raw expression and methylation matrices are independently imputed and min-max "
            "scaled to [0, 1] inside each applicable training partition before the selected "
            f"per-gene {config.correlation_method} correlation and MECor values are calculated."
        ),
        "leakage_control": (
            "Outer test samples are excluded from imputation, scaling, correlation estimation, "
            "MECor construction, epoch selection, and model fitting. With nested early stopping, "
            "inner folds select the median best epoch and a fresh model is refitted on the full "
            "outer training fold before one-time outer evaluation."
        ),
    }
    _write_json(output / "run_config.json", run_metadata)
    summaries = []
    cancer_iterator = tqdm(
        config.cancers,
        total=len(config.cancers),
        desc="ProgMap cancers",
        unit="cancer",
        disable=not config.progress,
        dynamic_ncols=True,
    )
    for cancer in cancer_iterator:
        summaries.append(run_cancer(config, cancer))
    pd.DataFrame(summaries).to_json(
        output / "run_summary.json", orient="records", indent=2, force_ascii=False
    )
    return summaries
