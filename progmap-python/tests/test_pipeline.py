import json

import pandas as pd

from progmap.pipeline import RunConfig, run


def test_end_to_end_training(small_data_root, tmp_path):
    output = tmp_path / "results"
    config = RunConfig(
        data_root=str(small_data_root),
        output=str(output),
        cancers=("TEST",),
        folds=2,
        inner_validation_fraction=0.25,
        epochs=1,
        patience=0,
        batch_size=8,
        ig_steps=3,
        ig_baselines=1,
        max_attribution_samples_per_class=2,
        test="ttest",
        top_n=3,
        save_models=False,
        verbose=0,
        device="cpu",
        threads=2,
    )
    summary = run(config)[0]
    assert summary["selected_feature_count"] == 3
    predictions = pd.read_csv(output / "TEST" / "cross_validated_predictions.csv")
    assert len(predictions) == 27
    assert predictions["fold"].nunique() == 2
    selected = pd.read_csv(output / "TEST" / "features_selected.csv")
    assert len(selected) == 3
    with (output / "run_config.json").open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    assert "inside each outer training fold" in metadata["leakage_control"]
    assert metadata["runtime_device"]["effective"] == "cpu"
    assert metadata["config"]["threads"] == 2
