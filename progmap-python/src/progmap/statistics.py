from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


TESTS = ("ttest", "wilcoxon", "permutation")
CORRECTIONS = (
    "bonferroni",
    "sidak",
    "holm-sidak",
    "holm",
    "simes-hochberg",
    "hommel",
    "fdr_bh",
    "fdr_by",
    "fdr_tsbh",
    "fdr_tsbky",
)


def _ttest(target: np.ndarray, background: np.ndarray) -> np.ndarray:
    result = stats.ttest_ind(
        target,
        background,
        axis=0,
        equal_var=False,
        nan_policy="omit",
        alternative="greater",
    )
    return np.asarray(result.pvalue, dtype=float)


def _wilcoxon_rank_sum(target: np.ndarray, background: np.ndarray) -> np.ndarray:
    pvalues = np.ones(target.shape[1], dtype=float)
    for index in range(target.shape[1]):
        try:
            pvalues[index] = stats.mannwhitneyu(
                target[:, index],
                background[:, index],
                alternative="greater",
                method="auto",
            ).pvalue
        except ValueError:
            pvalues[index] = 1.0
    return pvalues


def _permutation(
    target: np.ndarray,
    background: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    if iterations < 1:
        raise ValueError("Permutation testing requires at least one iteration")
    pooled = np.concatenate([target, background], axis=0)
    n_target = len(target)
    observed = target.mean(axis=0) - background.mean(axis=0)
    exceedances = np.zeros(target.shape[1], dtype=np.int64)
    for _ in range(iterations):
        order = rng.permutation(len(pooled))
        permuted_target = pooled[order[:n_target]]
        permuted_background = pooled[order[n_target:]]
        permuted_effect = permuted_target.mean(axis=0) - permuted_background.mean(axis=0)
        exceedances += permuted_effect >= observed
    return (exceedances + 1.0) / (iterations + 1.0)


def _bootstrap_intervals(
    target: np.ndarray,
    background: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
    confidence: float = 0.95,
    chunk_size: int = 512,
) -> tuple[np.ndarray, np.ndarray]:
    n_features = target.shape[1]
    low = np.empty(n_features, dtype=np.float32)
    high = np.empty(n_features, dtype=np.float32)
    alpha = (1.0 - confidence) / 2.0
    for start in range(0, n_features, chunk_size):
        stop = min(start + chunk_size, n_features)
        differences = np.empty((iterations, stop - start), dtype=np.float32)
        for iteration in range(iterations):
            target_rows = rng.integers(0, len(target), size=len(target))
            background_rows = rng.integers(0, len(background), size=len(background))
            differences[iteration] = (
                np.median(target[target_rows, start:stop], axis=0)
                - np.median(background[background_rows, start:stop], axis=0)
            )
        low[start:stop] = np.quantile(differences, alpha, axis=0)
        high[start:stop] = np.quantile(differences, 1.0 - alpha, axis=0)
    return low, high


def rank_features(
    attributions: np.ndarray,
    labels: np.ndarray,
    genes: np.ndarray,
    target_class: int,
    *,
    test: str = "ttest",
    correction: str = "fdr_bh",
    alpha: float = 0.01,
    permutations: int = 1000,
    bootstrap_iterations: int = 0,
    seed: int = 42,
) -> pd.DataFrame:
    if test not in TESTS:
        raise ValueError(f"Unknown test {test!r}; choose one of {TESTS}")
    if correction not in CORRECTIONS:
        raise ValueError(f"Unknown correction {correction!r}")
    values = np.abs(np.asarray(attributions, dtype=np.float32))
    labels = np.asarray(labels)
    target = values[labels == target_class]
    background = values[labels != target_class]
    if len(target) < 2 or len(background) < 2:
        raise ValueError(
            f"Feature testing for class {target_class} requires at least two target and two background samples"
        )
    rng = np.random.default_rng(seed + target_class)
    if test == "ttest":
        pvalue = _ttest(target, background)
    elif test == "wilcoxon":
        pvalue = _wilcoxon_rank_sum(target, background)
    else:
        pvalue = _permutation(target, background, permutations, rng)
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
    adjusted = multipletests(pvalue, alpha=alpha, method=correction)[1]
    target_median = np.median(target, axis=0)
    background_median = np.median(background, axis=0)
    result = pd.DataFrame(
        {
            "gene": np.asarray(genes, dtype=object),
            "target_class": int(target_class),
            "median_abs_attribution_target": target_median,
            "median_abs_attribution_background": background_median,
            "effect": target_median - background_median,
            "pvalue": pvalue,
            "pvalue_adjusted": adjusted,
            "significant": adjusted <= alpha,
            "test": test,
            "correction": correction,
        }
    )
    if bootstrap_iterations > 0:
        low, high = _bootstrap_intervals(
            target, background, bootstrap_iterations, rng
        )
        result["effect_ci_low"] = low
        result["effect_ci_high"] = high
    return result.sort_values(
        ["pvalue_adjusted", "effect", "gene"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def collapse_gene_ranking(class_results: pd.DataFrame) -> pd.DataFrame:
    ordered = class_results.sort_values(
        ["pvalue_adjusted", "effect"], ascending=[True, False], kind="mergesort"
    )
    best = ordered.drop_duplicates("gene", keep="first").copy()
    best = best.rename(
        columns={
            "target_class": "best_class",
            "pvalue_adjusted": "best_adjusted_pvalue",
            "effect": "best_effect",
            "significant": "significant_in_best_class",
        }
    )
    columns = [
        "gene",
        "best_class",
        "best_adjusted_pvalue",
        "best_effect",
        "significant_in_best_class",
        "test",
        "correction",
    ]
    return best[columns].sort_values(
        ["best_adjusted_pvalue", "best_effect", "gene"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

