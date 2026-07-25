from progmap.data import discover_cancers, load_cancer


def test_aliases_pairing_and_incomplete_directory_exclusion(small_data_root):
    assert discover_cancers(small_data_root) == ["TEST"]
    data = load_cancer(small_data_root, "TEST")
    assert data.expression.shape == (27, 8)
    assert data.methylation.shape == (27, 8)
    assert data.class_counts() == {"normal": 9, "stage_I": 9, "stage_II_III": 9}
    assert data.sample_ids[0].startswith("normal::")
