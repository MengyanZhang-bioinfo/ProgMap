from __future__ import annotations

import numpy as np

from .model import require_tensorflow


def select_median_baselines(
    x_train: np.ndarray,
    y_train: np.ndarray,
    target_class: int,
    n_baselines: int,
) -> np.ndarray:
    candidates = x_train[y_train != target_class]
    if not len(candidates):
        candidates = x_train
    n_baselines = max(1, min(int(n_baselines), len(candidates)))
    center = np.median(candidates, axis=0)
    distances = np.sum(np.abs(candidates - center), axis=1)
    return candidates[np.argsort(distances)[:n_baselines]].astype(np.float32)


def _path_gradients(model, path: np.ndarray, target_class: int, batch_size: int) -> np.ndarray:
    tf = require_tensorflow()
    gradients = []
    for start in range(0, len(path), batch_size):
        batch = tf.convert_to_tensor(path[start : start + batch_size], dtype=tf.float32)
        with tf.GradientTape() as tape:
            tape.watch(batch)
            output = model(batch, training=False)[:, target_class]
        gradients.append(tape.gradient(output, batch).numpy())
    return np.concatenate(gradients, axis=0)


def integrated_gradients(
    model,
    samples: np.ndarray,
    baselines: np.ndarray,
    target_class: int,
    *,
    steps: int = 64,
    batch_size: int = 32,
) -> np.ndarray:
    """Enhanced IG using several training-derived median baselines and a linear path."""
    if steps < 2:
        raise ValueError("steps must be at least 2")
    alphas = np.linspace(0.0, 1.0, steps, dtype=np.float32).reshape(-1, 1)
    all_samples = []
    for sample in np.asarray(samples, dtype=np.float32):
        baseline_attributions = []
        for baseline in np.asarray(baselines, dtype=np.float32):
            delta = sample - baseline
            path = baseline.reshape(1, -1) + alphas * delta.reshape(1, -1)
            gradients = _path_gradients(model, path, target_class, batch_size)
            # Integral over a unit-length, uniformly sampled path.  Writing the
            # trapezoidal mean explicitly keeps compatibility with NumPy 1.26.
            average_gradient = np.mean(
                0.5 * (gradients[:-1] + gradients[1:]), axis=0
            )
            baseline_attributions.append(delta * average_gradient)
        all_samples.append(np.median(np.stack(baseline_attributions, axis=0), axis=0))
    return np.asarray(all_samples, dtype=np.float32)


def balanced_attribution_indices(
    labels: np.ndarray,
    max_per_class: int | None,
    rng: np.random.Generator,
) -> np.ndarray:
    indices = []
    for label in (0, 1, 2):
        available = np.where(labels == label)[0]
        if max_per_class is not None and len(available) > max_per_class:
            available = np.sort(rng.choice(available, size=max_per_class, replace=False))
        indices.extend(available.tolist())
    return np.asarray(sorted(indices), dtype=int)
