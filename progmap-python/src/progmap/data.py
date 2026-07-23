from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


CLASS_NAMES = {0: "normal", 1: "stage_I", 2: "stage_II_III"}

FILE_CANDIDATES = {
    0: {"expression": ("en.csv",), "methylation": ("mn.csv",)},
    1: {"expression": ("exp1.csv", "e1.csv"), "methylation": ("me1.csv", "m1.csv")},
    2: {"expression": ("exp2.csv", "e2.csv"), "methylation": ("me2.csv", "m2.csv")},
}


@dataclass(frozen=True)
class MultiOmicsData:
    expression: np.ndarray
    methylation: np.ndarray
    labels: np.ndarray
    sample_ids: np.ndarray
    genes: np.ndarray
    cancer: str

    @property
    def n_samples(self) -> int:
        return int(self.labels.shape[0])

    @property
    def n_genes(self) -> int:
        return int(self.genes.shape[0])

    def class_counts(self) -> dict[str, int]:
        return {
            CLASS_NAMES[label]: int(np.sum(self.labels == label))
            for label in sorted(CLASS_NAMES)
        }


def _resolve_file(folder: Path, candidates: Iterable[str]) -> Path:
    for name in candidates:
        path = folder / name
        if path.exists():
            return path
    raise FileNotFoundError(
        f"None of the expected files exist in {folder}: {', '.join(candidates)}"
    )


def _read_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, index_col=0)
    if matrix.empty:
        raise ValueError(f"Empty matrix: {path}")
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    if matrix.index.has_duplicates:
        duplicates = matrix.index[matrix.index.duplicated()].unique().tolist()[:5]
        raise ValueError(f"Duplicate gene identifiers in {path}: {duplicates}")
    if matrix.columns.has_duplicates:
        duplicates = matrix.columns[matrix.columns.duplicated()].unique().tolist()[:5]
        raise ValueError(f"Duplicate sample identifiers in {path}: {duplicates}")
    return matrix.apply(pd.to_numeric, errors="coerce")


def discover_cancers(data_root: str | Path) -> list[str]:
    root = Path(data_root)
    if not root.is_dir():
        raise NotADirectoryError(root)
    cancers = []
    for folder in sorted(root.iterdir()):
        if not folder.is_dir() or folder.name.upper() == "GEO":
            continue
        try:
            for label in FILE_CANDIDATES:
                _resolve_file(folder, FILE_CANDIDATES[label]["expression"])
                _resolve_file(folder, FILE_CANDIDATES[label]["methylation"])
        except FileNotFoundError:
            continue
        cancers.append(folder.name)
    return cancers


def load_cancer(data_root: str | Path, cancer: str) -> MultiOmicsData:
    folder = Path(data_root) / cancer
    if not folder.is_dir():
        raise NotADirectoryError(folder)

    raw: dict[int, dict[str, pd.DataFrame]] = {}
    for label, modalities in FILE_CANDIDATES.items():
        raw[label] = {
            modality: _read_matrix(_resolve_file(folder, names))
            for modality, names in modalities.items()
        }

    first_index = raw[0]["expression"].index
    common_genes = set(first_index)
    for label in raw:
        common_genes &= set(raw[label]["expression"].index)
        common_genes &= set(raw[label]["methylation"].index)
    genes = [gene for gene in first_index if gene in common_genes]
    if not genes:
        raise ValueError(f"No genes are shared by all six matrices in {folder}")

    expression_parts: list[np.ndarray] = []
    methylation_parts: list[np.ndarray] = []
    label_parts: list[np.ndarray] = []
    sample_parts: list[np.ndarray] = []

    for label in sorted(raw):
        expression = raw[label]["expression"]
        methylation = raw[label]["methylation"]
        common_samples = [sample for sample in expression.columns if sample in methylation.columns]
        if not common_samples:
            raise ValueError(
                f"No paired expression/methylation samples for {cancer}, class {CLASS_NAMES[label]}"
            )
        expression_values = expression.loc[genes, common_samples].T.to_numpy(dtype=np.float32)
        methylation_values = methylation.loc[genes, common_samples].T.to_numpy(dtype=np.float32)
        expression_parts.append(expression_values)
        methylation_parts.append(methylation_values)
        label_parts.append(np.full(len(common_samples), label, dtype=np.int64))
        sample_parts.append(
            np.asarray([f"{CLASS_NAMES[label]}::{sample}" for sample in common_samples], dtype=object)
        )

    return MultiOmicsData(
        expression=np.concatenate(expression_parts, axis=0),
        methylation=np.concatenate(methylation_parts, axis=0),
        labels=np.concatenate(label_parts),
        sample_ids=np.concatenate(sample_parts),
        genes=np.asarray(genes, dtype=object),
        cancer=cancer,
    )

