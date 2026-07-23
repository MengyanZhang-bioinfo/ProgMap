from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

from . import __version__
from .data import discover_cancers
from .pipeline import RunConfig, run
from .statistics import CORRECTIONS, TESTS


def _top_n(value: str) -> int | None:
    if value.lower() == "all":
        return None
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("--top-n must be a positive integer or 'all'")
    return number


def _optional_positive_int(value: str) -> int | None:
    if value.lower() == "all":
        return None
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer or 'all'")
    return number


def _class_weights(value: str) -> tuple[float, float, float]:
    try:
        weights = tuple(float(part.strip()) for part in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("class weights must be three comma-separated numbers") from exc
    if len(weights) != 3 or any(weight <= 0 for weight in weights):
        raise argparse.ArgumentTypeError("class weights must contain three positive values")
    return weights


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="progmap",
        description=(
            "Fold-safe MECor (methylation × expression × training-fold correlation), "
            "fixed 2048-128 skip model, enhanced integrated gradients, and feature testing."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--data-root", required=True, help="PANCANCER root containing one folder per cancer")
    parser.add_argument("--output", required=True, help="Result directory")
    parser.add_argument(
        "--cancers",
        default="all",
        help="Comma-separated cancer folders, or 'all' (GEO is always excluded)",
    )
    parser.add_argument("--folds", type=int, default=5, help="Outer stratified CV folds; auto-reduced to the smallest class")
    parser.add_argument("--inner-validation-fraction", type=float, default=0.15)
    parser.add_argument("--epochs", type=int, default=1000)
    parser.add_argument("--patience", type=int, default=200, help="Early stopping patience on inner validation AUC; 0 disables")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--decay-steps", type=int, default=10_000)
    parser.add_argument("--decay-rate", type=float, default=0.9)
    parser.add_argument("--class-weights", type=_class_weights, default=(0.25, 0.5, 0.25), metavar="W0,W1,W2")
    parser.add_argument("--imputation", choices=("knn", "median"), default="knn")
    parser.add_argument("--knn-neighbors", type=int, default=10)
    parser.add_argument("--ig-steps", type=int, default=64)
    parser.add_argument("--ig-baselines", type=int, default=3)
    parser.add_argument("--attribution-batch-size", type=int, default=32)
    parser.add_argument(
        "--max-attribution-samples-per-class",
        type=_optional_positive_int,
        default=None,
        metavar="N|all",
        help="Randomly cap held-out attribution samples per true class in each fold",
    )
    parser.add_argument("--test", choices=TESTS, default="ttest", help="Feature significance test; wilcoxon is rank-sum/Mann-Whitney")
    parser.add_argument("--correction", choices=CORRECTIONS, default="fdr_bh")
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--bootstrap-iterations", type=int, default=0, help="0 disables feature-effect bootstrap CIs")
    parser.add_argument("--top-n", type=_top_n, default=None, metavar="N|all", help="Number of unique ranked genes in features_selected.csv")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--save-attributions", action="store_true")
    parser.add_argument("--no-save-models", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Validate and summarize inputs without TensorFlow training")
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "gpu"),
        default="auto",
        help="Execution device; 'gpu' fails early when TensorFlow cannot see a GPU",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="CPU threads for TensorFlow/OpenMP; omit to use the server runtime default",
    )
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.folds < 2:
        parser.error("--folds must be at least 2")
    if args.epochs < 1:
        parser.error("--epochs must be at least 1")
    if not 0.0 <= args.dropout < 1.0:
        parser.error("--dropout must be in [0, 1)")
    if not 0.0 < args.alpha < 1.0:
        parser.error("--alpha must be between 0 and 1")
    if args.knn_neighbors < 1:
        parser.error("--knn-neighbors must be positive")
    if args.threads is not None and args.threads < 1:
        parser.error("--threads must be a positive integer")
    if args.device == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    if args.threads is not None:
        thread_count = str(args.threads)
        os.environ["OMP_NUM_THREADS"] = thread_count
        os.environ["TF_NUM_INTRAOP_THREADS"] = thread_count
        os.environ["TF_NUM_INTEROP_THREADS"] = "1"

    available = discover_cancers(args.data_root)
    if not available:
        parser.error("No cancer folders with the six required matrices were found")
    if args.cancers.lower() == "all":
        cancers = tuple(available)
    else:
        requested = tuple(part.strip() for part in args.cancers.split(",") if part.strip())
        missing = sorted(set(requested) - set(available))
        if missing:
            parser.error(f"Cancer folders not found or incomplete: {', '.join(missing)}")
        cancers = requested

    config = RunConfig(
        data_root=str(Path(args.data_root).resolve()),
        output=str(Path(args.output).resolve()),
        cancers=cancers,
        folds=args.folds,
        inner_validation_fraction=args.inner_validation_fraction,
        epochs=args.epochs,
        patience=args.patience,
        batch_size=args.batch_size,
        dropout=args.dropout,
        learning_rate=args.learning_rate,
        decay_steps=args.decay_steps,
        decay_rate=args.decay_rate,
        class_weights=args.class_weights,
        imputation=args.imputation,
        knn_neighbors=args.knn_neighbors,
        ig_steps=args.ig_steps,
        ig_baselines=args.ig_baselines,
        attribution_batch_size=args.attribution_batch_size,
        max_attribution_samples_per_class=args.max_attribution_samples_per_class,
        test=args.test,
        correction=args.correction,
        alpha=args.alpha,
        permutations=args.permutations,
        bootstrap_iterations=args.bootstrap_iterations,
        top_n=args.top_n,
        seed=args.seed,
        save_attributions=args.save_attributions,
        save_models=not args.no_save_models,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
        verbose=0 if args.quiet else 1,
        device=args.device,
        threads=args.threads,
    )
    try:
        summaries = run(config)
    except Exception as exc:
        print(f"progmap failed: {exc}", file=sys.stderr)
        raise
    for summary in summaries:
        print(
            f"{summary['cancer']}: samples={summary['samples']}, genes={summary['genes']}, "
            f"folds={summary['used_folds']}"
        )
    print(f"Results written to {config.output}")


if __name__ == "__main__":
    main()
