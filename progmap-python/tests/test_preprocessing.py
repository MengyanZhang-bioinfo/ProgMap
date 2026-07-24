import numpy as np

from progmap.preprocessing import FoldPreprocessor, pearson_per_gene, spearman_per_gene


def test_preprocessing_is_fit_only_on_training_data(tmp_path):
    expression_train = np.asarray(
        [[0, 4, 2], [1, 3, 2], [2, 2, 2], [3, 1, 2]], dtype=np.float32
    )
    methylation_train = np.asarray(
        [[0, 1, 3], [1, 2, 3], [2, 3, 3], [3, 4, 3]], dtype=np.float32
    )
    processor = FoldPreprocessor().fit(
        expression_train, methylation_train, np.asarray(["A", "B", "C"])
    )
    expected = pearson_per_gene(*processor.transform_modalities(expression_train, methylation_train))
    np.testing.assert_allclose(processor.correlations, expected)

    correlations_before = processor.correlations.copy()
    extreme_test = np.asarray([[1e9, -1e9, np.nan]], dtype=np.float32)
    expression_scaled, methylation_scaled = processor.transform_modalities(
        extreme_test, -extreme_test
    )
    assert np.min(expression_scaled) >= 0 and np.max(expression_scaled) <= 1
    assert np.min(methylation_scaled) >= 0 and np.max(methylation_scaled) <= 1
    np.testing.assert_array_equal(processor.correlations, correlations_before)

    path = tmp_path / "processor.joblib"
    processor.save(path)
    loaded = FoldPreprocessor.load(path)
    np.testing.assert_allclose(loaded.correlations, correlations_before)


def test_spearman_correlation_and_configurable_preprocessor():
    expression = np.asarray(
        [[1, 4], [2, 1], [3, 3], [4, 2]], dtype=np.float32
    )
    methylation = np.asarray(
        [[10, 1], [20, 4], [30, 3], [40, 2]], dtype=np.float32
    )
    expected = np.asarray([1.0, -0.8], dtype=np.float32)
    np.testing.assert_allclose(
        spearman_per_gene(expression, methylation), expected, atol=1e-6
    )
    processor = FoldPreprocessor(
        imputation="median", correlation_method="spearman"
    ).fit(expression, methylation, np.asarray(["A", "B"]))
    np.testing.assert_allclose(processor.correlations, expected, atol=1e-6)
