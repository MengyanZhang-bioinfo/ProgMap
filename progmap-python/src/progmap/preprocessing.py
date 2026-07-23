from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import joblib
import numpy as np
from sklearn.impute import KNNImputer
from sklearn.preprocessing import MinMaxScaler


def _nan_medians(values: np.ndarray) -> np.ndarray:
    with np.errstate(all="ignore"):
        medians = np.nanmedian(values, axis=0)
    return np.where(np.isfinite(medians), medians, 0.0).astype(np.float32)


def _impute(values: np.ndarray, medians: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=np.float32).copy()
    rows, cols = np.where(~np.isfinite(result))
    if len(rows):
        result[rows, cols] = medians[cols]
    return result


def pearson_per_gene(expression: np.ndarray, methylation: np.ndarray) -> np.ndarray:
    if expression.shape != methylation.shape:
        raise ValueError("Expression and methylation matrices must have identical shapes")
    expression_centered = expression - expression.mean(axis=0, keepdims=True)
    methylation_centered = methylation - methylation.mean(axis=0, keepdims=True)
    numerator = np.sum(expression_centered * methylation_centered, axis=0)
    denominator = np.sqrt(
        np.sum(expression_centered**2, axis=0)
        * np.sum(methylation_centered**2, axis=0)
    )
    correlation = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=np.float32),
        where=denominator > np.finfo(np.float32).eps,
    )
    return np.clip(correlation, -1.0, 1.0).astype(np.float32)


@dataclass
class FoldPreprocessor:
    """Training-fold-only imputation, scaling, correlation, and MECor construction."""

    imputation: str = "knn"
    knn_neighbors: int = 10
    expression_medians: np.ndarray | None = None
    methylation_medians: np.ndarray | None = None
    expression_imputer: KNNImputer | None = None
    methylation_imputer: KNNImputer | None = None
    expression_scaler: MinMaxScaler | None = None
    methylation_scaler: MinMaxScaler | None = None
    correlations: np.ndarray | None = None
    genes: np.ndarray | None = None

    def fit(
        self,
        expression_train: np.ndarray,
        methylation_train: np.ndarray,
        genes: np.ndarray,
    ) -> "FoldPreprocessor":
        if expression_train.shape != methylation_train.shape:
            raise ValueError("Expression and methylation matrices must have identical shapes")
        if expression_train.shape[1] != len(genes):
            raise ValueError("Gene count does not match the number of matrix columns")
        if self.imputation not in {"knn", "median"}:
            raise ValueError("imputation must be 'knn' or 'median'")
        if self.knn_neighbors < 1:
            raise ValueError("knn_neighbors must be positive")
        self.genes = np.asarray(genes, dtype=object)
        self.expression_medians = _nan_medians(expression_train)
        self.methylation_medians = _nan_medians(methylation_train)
        if self.imputation == "knn":
            neighbors = min(self.knn_neighbors, len(expression_train))
            self.expression_imputer = KNNImputer(
                n_neighbors=neighbors, keep_empty_features=True
            )
            self.methylation_imputer = KNNImputer(
                n_neighbors=neighbors, keep_empty_features=True
            )
            expression = self.expression_imputer.fit_transform(expression_train).astype(
                np.float32
            )
            methylation = self.methylation_imputer.fit_transform(methylation_train).astype(
                np.float32
            )
        else:
            expression = _impute(expression_train, self.expression_medians)
            methylation = _impute(methylation_train, self.methylation_medians)
        self.expression_scaler = MinMaxScaler(feature_range=(0.0, 1.0), clip=True)
        self.methylation_scaler = MinMaxScaler(feature_range=(0.0, 1.0), clip=True)
        expression_scaled = self.expression_scaler.fit_transform(expression).astype(np.float32)
        methylation_scaled = self.methylation_scaler.fit_transform(methylation).astype(np.float32)
        self.correlations = pearson_per_gene(expression_scaled, methylation_scaled)
        return self

    def transform_modalities(
        self, expression: np.ndarray, methylation: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        if any(
            item is None
            for item in (
                self.expression_medians,
                self.methylation_medians,
                self.expression_scaler,
                self.methylation_scaler,
                self.correlations,
            )
        ):
            raise RuntimeError("FoldPreprocessor has not been fitted")
        if self.imputation == "knn":
            if self.expression_imputer is None or self.methylation_imputer is None:
                raise RuntimeError("KNN imputers have not been fitted")
            expression_clean = self.expression_imputer.transform(expression).astype(np.float32)
            methylation_clean = self.methylation_imputer.transform(methylation).astype(np.float32)
        else:
            expression_clean = _impute(expression, self.expression_medians)
            methylation_clean = _impute(methylation, self.methylation_medians)
        expression_scaled = self.expression_scaler.transform(expression_clean).astype(np.float32)
        methylation_scaled = self.methylation_scaler.transform(methylation_clean).astype(np.float32)
        return expression_scaled, methylation_scaled

    def transform(self, expression: np.ndarray, methylation: np.ndarray) -> np.ndarray:
        expression_scaled, methylation_scaled = self.transform_modalities(expression, methylation)
        return (
            expression_scaled * methylation_scaled * self.correlations.reshape(1, -1)
        ).astype(np.float32)

    def save(self, path: str | Path) -> None:
        joblib.dump(self, path)

    @classmethod
    def load(cls, path: str | Path) -> "FoldPreprocessor":
        loaded = joblib.load(path)
        if not isinstance(loaded, cls):
            raise TypeError(f"Unexpected preprocessor type in {path}: {type(loaded)!r}")
        return loaded
