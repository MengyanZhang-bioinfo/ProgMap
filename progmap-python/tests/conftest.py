from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def write_matrix(path: Path, genes: list[str], samples: list[str], values: np.ndarray) -> None:
    pd.DataFrame(values, index=genes, columns=samples).to_csv(path)


@pytest.fixture
def small_data_root(tmp_path: Path) -> Path:
    root = tmp_path / "datasets"
    cancer = root / "TEST"
    cancer.mkdir(parents=True)
    genes = [f"GENE_{index}" for index in range(8)]
    rng = np.random.default_rng(7)
    aliases = {
        0: ("en.csv", "mn.csv"),
        1: ("e1.csv", "m1.csv"),
        2: ("e2.csv", "m2.csv"),
    }
    for label, (expression_name, methylation_name) in aliases.items():
        samples = [f"S{label}_{index}" for index in range(9)]
        base = label * 0.8
        expression = rng.normal(base, 0.25, size=(len(genes), len(samples)))
        methylation = 0.5 * expression + rng.normal(base, 0.2, size=expression.shape)
        write_matrix(cancer / expression_name, genes, samples, expression)
        write_matrix(cancer / methylation_name, genes, samples, methylation)

    incomplete = root / "INCOMPLETE"
    incomplete.mkdir()
    (incomplete / "en.csv").write_text("gene,S\nG,1\n", encoding="utf-8")
    return root
