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

from .attribution import (
    balanced_attribution_indices,
    integrated_gradients,
    select_median_baselines,
)
from .data import CLASS_NAMES, load_cancer
from .model import (
    build_skip_model,
    configure_runtime,
    fit_model,
    require_tensorflow,
    set_global_seed,
)
from .preprocessing import FoldPreprocessor
from .statistics import collapse_gene_ranking, rank_features


@dataclass(frozen=True)
class RunConfig:
    data_root: str
    output: str
    cancers: tuple[str, ...]
    folds: int = 5
    inner_validation_fraction: float = 0.15
    epochs: int = 1000
    patience: int = 200
    batch_size: int = 16
    dropout: float = 0.1
    learning_rate: float = 1e-3
    decay_steps: int = 10_000
    decay_rate: float = 0.9
    class_weights: tuple[float, float, float] = (0.25, 0.5, 0.25)
    imputation: str = "knn"
    knn_neighbors: int = 10
    ig_steps: int = 64
    ig_baselines: int = 3
    attribution_batch_size: int = 32
    max_attribution_samples_per_class: int | None = None
    test: str = "ttest"
    correction: str = "fdr_bh"
    alpha: float = 0.01
    permutations: int = 1000
    bootstrap_iterations: int = 0
    top_n: int | None = None
    seed: int = 42
    save_attributions: bool = False
    save_models: bool = True
    overwrite: bool = False
    dry_run: bool = False
    verbose: int = 1
    device: str = "auto"
    threads: int | None = None


def _versions() -> dict[str, str]:
    packages = ["progmap", "tensorflow", "numpy", "pandas", "scikit-learn", "scipy", "statsmodels"]
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


def _fold_metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predictions = np.argmax(probabilities, axis=1)
    one_hot = np.eye(3, dtype=np.float32)[labels]
    metrics = {
        "accuracy": float(accuracy_score(labels, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predictions)),
        "auc_micro": float(roc_auc_score(one_hot, probabilities, average="micro", multi_class="ovr")),
        "auc_macro": float(roc_auc_score(one_hot, probabilities, average="macro", multi_class="ovr")),
    }
    for label in (0, 1, 2):
        metrics[f"auc_class_{label}"] = float(roc_auc_score(one_hot[:, label], probabilities[:, label]))
    return metrics


def _save_correlations(path: Path, genes: np.ndarray, correlations: np.ndarray) -> None:
    pd.DataFrame({"gene": genes, "training_fold_pearson_correlation": correlations}).to_csv(
        path, index=False
    )


def run_cancer(config: RunConfig, cancer: str) -> dict[str, Any]:
    data = load_cancer(config.data_root, cancer)
    cancer_output = Path(config.output) / cancer
    _prepare_output(cancer_output, config.overwrite)
    class_counts = data.class_counts()
    minimum_class_count = min(class_counts.values())
    folds = min(config.folds, minimum_class_count)
    if folds < 2:
        raise ValueError(f"{cancer} needs at least two samples in every class")

    metadata = {
        "cancer": cancer,
        "samples": data.n_samples,
        "genes": data.n_genes,
        "class_counts": class_counts,
        "requested_folds": config.folds,
        "used_folds": folds,
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
    attribution_values: dict[int, list[np.ndarray]] = {0: [], 1: [], 2: []}
    attribution_labels: dict[int, list[np.ndarray]] = {0: [], 1: [], 2: []}
    attribution_samples: dict[int, list[np.ndarray]] = {0: [], 1: [], 2: []}

    for fold, (outer_train, outer_test) in enumerate(
        splitter.split(data.expression, data.labels), start=1
    ):
        fold_seed = config.seed + fold
        set_global_seed(fold_seed)
        fold_output = cancer_output / f"fold_{fold}"
        fold_output.mkdir(parents=True, exist_ok=True)
        fit_indices, validation_indices = _inner_split(
            np.asarray(outer_train), data.labels, config.inner_validation_fraction, fold_seed
        )

        preprocessor = FoldPreprocessor(
            imputation=config.imputation, knn_neighbors=config.knn_neighbors
        ).fit(
            data.expression[fit_indices], data.methylation[fit_indices], data.genes
        )
        x_fit = preprocessor.transform(
            data.expression[fit_indices], data.methylation[fit_indices]
        )
        x_validation = preprocessor.transform(
            data.expression[validation_indices], data.methylation[validation_indices]
        )
        x_test = preprocessor.transform(
            data.expression[outer_test], data.methylation[outer_test]
        )
        y_fit = data.labels[fit_indices]
        y_validation = data.labels[validation_indices]
        y_test = data.labels[outer_test]

        preprocessor.save(fold_output / "preprocessor.joblib")
        _save_correlations(
            fold_output / "training_fold_correlations.csv", data.genes, preprocessor.correlations
        )
        pd.DataFrame(
            {
                "sample_id": data.sample_ids,
                "role": np.select(
                    [
                        np.isin(np.arange(data.n_samples), fit_indices),
                        np.isin(np.arange(data.n_samples), validation_indices),
                        np.isin(np.arange(data.n_samples), outer_test),
                    ],
                    ["fit", "internal_validation", "outer_test"],
                    default="unused",
                ),
                "label": data.labels,
            }
        ).to_csv(fold_output / "sample_roles.csv", index=False)

        model = build_skip_model(
            data.n_genes,
            dropout=config.dropout,
            learning_rate=config.learning_rate,
            decay_steps=config.decay_steps,
            decay_rate=config.decay_rate,
        )
        history = fit_model(
            model,
            x_fit,
            y_fit,
            x_validation,
            y_validation,
            epochs=config.epochs,
            patience=config.patience,
            batch_size=config.batch_size,
            class_weight={index: weight for index, weight in enumerate(config.class_weights)},
            verbose=config.verbose,
        )
        pd.DataFrame(history).to_csv(fold_output / "training_history.csv", index=False)
        if config.save_models:
            model.save(fold_output / "model.keras")

        probabilities = model.predict(x_test, batch_size=config.batch_size, verbose=0)
        fold_metrics = _fold_metrics(y_test, probabilities)
        fold_metrics["fold"] = fold
        fold_metrics["fit_samples"] = len(fit_indices)
        fold_metrics["internal_validation_samples"] = len(validation_indices)
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
        for target_class in (0, 1, 2):
            baselines = select_median_baselines(
                x_fit, y_fit, target_class, config.ig_baselines
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
    predictions.to_csv(cancer_output / "cross_validated_predictions.csv", index=False)
    metrics.to_csv(cancer_output / "fold_metrics.csv", index=False)

    feature_frames = []
    for target_class in (0, 1, 2):
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
    selected = gene_ranking if config.top_n is None else gene_ranking.head(config.top_n)
    selected.to_csv(cancer_output / "features_selected.csv", index=False)

    summary = {
        **metadata,
        "mean_metrics": {
            column: float(metrics[column].mean())
            for column in metrics.columns
            if column not in {"fold"} and pd.api.types.is_numeric_dtype(metrics[column])
        },
        "feature_test": config.test,
        "multiple_testing_correction": config.correction,
        "alpha": config.alpha,
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
        "leakage_control": (
            "Imputation, expression scaling, methylation scaling, Pearson correlation, "
            "MECor construction and model fitting are learned inside each outer training fold. "
            "Early stopping uses a separate split inside that outer training fold."
        ),
    }
    _write_json(output / "run_config.json", run_metadata)
    summaries = []
    for cancer in config.cancers:
        summaries.append(run_cancer(config, cancer))
    pd.DataFrame(summaries).to_json(
        output / "run_summary.json", orient="records", indent=2, force_ascii=False
    )
    return summaries
