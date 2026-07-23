import numpy as np

from progmap.attribution import integrated_gradients
from progmap.model import build_skip_model


def test_fixed_architecture_and_integrated_gradients():
    model = build_skip_model(7)
    assert model.get_layer("dense_2048").units == 2048
    assert model.get_layer("dense_128").units == 128
    assert model.get_layer("stage_probabilities").units == 3
    assert model.output_shape == (None, 3)
    samples = np.ones((2, 7), dtype=np.float32)
    baselines = np.zeros((1, 7), dtype=np.float32)
    values = integrated_gradients(model, samples, baselines, 1, steps=3, batch_size=4)
    assert values.shape == samples.shape
    assert np.isfinite(values).all()
