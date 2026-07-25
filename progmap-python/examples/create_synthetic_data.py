from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a small ProgMap-format test dataset")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--dataset",
        default="DEMO_DATASET",
        help="Name of the dataset directory to create",
    )
    parser.add_argument("--genes", type=int, default=8)
    parser.add_argument("--samples-per-class", type=int, default=9)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    if args.genes < 2 or args.samples_per_class < 4:
        parser.error("Use at least 2 genes and 4 samples per class")

    dataset_dir = Path(args.output).expanduser().resolve() / args.dataset
    dataset_dir.mkdir(parents=True, exist_ok=True)
    genes = [f"GENE_{index}" for index in range(args.genes)]
    rng = np.random.default_rng(args.seed)
    names = {
        0: ("en.csv", "mn.csv"),
        1: ("e1.csv", "m1.csv"),
        2: ("e2.csv", "m2.csv"),
    }
    for label, (expression_name, methylation_name) in names.items():
        samples = [
            f"{args.dataset}_C{label}_S{index}"
            for index in range(args.samples_per_class)
        ]
        base = label * 0.8
        expression = rng.normal(base, 0.25, size=(args.genes, len(samples)))
        methylation = 0.5 * expression + rng.normal(base, 0.2, size=expression.shape)
        pd.DataFrame(expression, index=genes, columns=samples).to_csv(
            dataset_dir / expression_name
        )
        pd.DataFrame(methylation, index=genes, columns=samples).to_csv(
            dataset_dir / methylation_name
        )
    print(dataset_dir.parent)


if __name__ == "__main__":
    main()
