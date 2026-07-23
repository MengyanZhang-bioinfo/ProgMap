import numpy as np
import pytest

from progmap.statistics import rank_features


@pytest.mark.parametrize("test", ["ttest", "wilcoxon", "permutation"])
def test_feature_tests_rank_signal_first(test):
    rng = np.random.default_rng(11)
    labels = np.repeat([0, 1, 2], 12)
    values = rng.normal(0, 0.05, size=(36, 5)).astype(np.float32)
    values[labels == 1, 0] += 4.0
    ranked = rank_features(
        values,
        labels,
        np.asarray(["SIGNAL", "B", "C", "D", "E"]),
        1,
        test=test,
        correction="fdr_bh",
        alpha=0.05,
        permutations=99,
        seed=3,
    )
    assert ranked.iloc[0]["gene"] == "SIGNAL"
    assert ranked.iloc[0]["effect"] > 1
