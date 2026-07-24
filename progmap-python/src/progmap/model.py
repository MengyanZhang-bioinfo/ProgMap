from __future__ import annotations

from dataclasses import dataclass
import os
import random
from typing import Any

import numpy as np


EARLY_STOPPING_STRATEGIES = ("nested", "holdout", "off")
EARLY_STOPPING_MONITORS = ("val_loss", "val_auc")


@dataclass(frozen=True)
class FitResult:
    history: dict[str, list[float]]
    best_epoch: int
    stopped_epoch: int
    best_monitor_value: float | None


def require_tensorflow():
    try:
        import tensorflow as tf
    except ImportError as exc:
        raise RuntimeError(
            "TensorFlow is required for model training. Install the package with `pip install progmap`."
        ) from exc
    return tf


def configure_runtime(requested_device: str, threads: int | None = None) -> dict[str, Any]:
    """Validate the requested TensorFlow device and enable safe GPU memory growth."""
    if requested_device == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    if threads is not None:
        os.environ["OMP_NUM_THREADS"] = str(threads)
        os.environ["TF_NUM_INTRAOP_THREADS"] = str(threads)
        os.environ["TF_NUM_INTEROP_THREADS"] = "1"
    tf = require_tensorflow()
    gpus = tf.config.list_physical_devices("GPU")
    if requested_device == "gpu" and not gpus:
        raise RuntimeError(
            "--device gpu was requested, but TensorFlow cannot see a GPU. "
            "Check the NVIDIA driver/container runtime or use --device cpu."
        )
    if requested_device == "cpu":
        try:
            tf.config.set_visible_devices([], "GPU")
        except RuntimeError as exc:
            raise RuntimeError(
                "TensorFlow was initialized before CPU-only mode could be applied. "
                "Start a new process and set --device cpu before importing TensorFlow."
            ) from exc
        gpus = []
    else:
        for gpu in gpus:
            try:
                tf.config.experimental.set_memory_growth(gpu, True)
            except RuntimeError:
                pass
    if threads is not None:
        try:
            tf.config.threading.set_intra_op_parallelism_threads(threads)
            tf.config.threading.set_inter_op_parallelism_threads(1)
        except RuntimeError:
            pass
    return {
        "requested": requested_device,
        "gpu_count": len(gpus),
        "gpu_devices": [device.name for device in gpus],
        "effective": "gpu" if gpus and requested_device != "cpu" else "cpu",
    }


def set_global_seed(seed: int, deterministic: bool = True) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    tf = require_tensorflow()
    tf.keras.utils.set_random_seed(seed)
    if deterministic:
        try:
            tf.config.experimental.enable_op_determinism()
        except Exception:
            pass


def balanced_sample_weights(labels: np.ndarray) -> np.ndarray:
    """Give every represented class equal total weight."""
    labels = np.asarray(labels, dtype=np.int64)
    counts = np.bincount(labels, minlength=3)
    represented = counts > 0
    if not represented.all():
        raise ValueError("Balanced validation loss requires all three classes")
    weights = len(labels) / (represented.sum() * counts.astype(np.float64))
    return weights[labels].astype(np.float32)


def _monitor_summary(
    history: dict[str, list[float]],
    monitor: str,
    warmup_epochs: int,
) -> tuple[int, float | None]:
    values = np.asarray(history.get(monitor, []), dtype=float)
    if not len(values):
        return len(next(iter(history.values()), [])), None
    start = min(max(int(warmup_epochs), 0), max(len(values) - 1, 0))
    eligible = values[start:]
    finite = np.isfinite(eligible)
    if not finite.any():
        return len(values), None
    if monitor == "val_loss":
        local_index = int(np.nanargmin(np.where(finite, eligible, np.nan)))
    else:
        local_index = int(np.nanargmax(np.where(finite, eligible, np.nan)))
    index = start + local_index
    return index + 1, float(values[index])


def _epoch_progress_callback(description: str, epochs: int):
    tf = require_tensorflow()
    from tqdm.auto import tqdm

    bar = tqdm(total=epochs, desc=description, unit="epoch", leave=False, dynamic_ncols=True)

    class EpochProgress(tf.keras.callbacks.Callback):
        def on_epoch_end(self, epoch, logs=None):
            logs = logs or {}
            bar.update(1)
            shown = {
                key: f"{float(logs[key]):.4g}"
                for key in ("loss", "val_loss", "val_auc")
                if key in logs
            }
            if shown:
                bar.set_postfix(shown, refresh=False)

        def on_train_end(self, logs=None):
            bar.close()

    return EpochProgress()


def build_skip_model(
    n_features: int,
    dropout: float = 0.1,
    learning_rate: float = 1e-3,
    decay_steps: int = 10_000,
    decay_rate: float = 0.9,
):
    """Fixed final ProgMap architecture: 2048 -> 128 -> input skip -> 3 classes."""
    tf = require_tensorflow()
    keras = tf.keras
    constraint = keras.constraints.MaxNorm(3.0)
    inputs = keras.Input(shape=(n_features,), name="mecor_input")
    x = keras.layers.Dense(
        2048,
        activation="relu",
        kernel_constraint=constraint,
        bias_constraint=constraint,
        name="dense_2048",
    )(inputs)
    x = keras.layers.BatchNormalization(name="batch_norm_2048")(x)
    x = keras.layers.Dropout(dropout, name="dropout_2048")(x)
    x = keras.layers.Dense(
        128,
        activation="relu",
        kernel_constraint=constraint,
        bias_constraint=constraint,
        name="dense_128",
    )(x)
    x = keras.layers.BatchNormalization(name="batch_norm_128")(x)
    x = keras.layers.Dropout(dropout, name="dropout_128")(x)
    skip = keras.layers.Concatenate(name="input_skip")([inputs, x])
    outputs = keras.layers.Dense(3, activation="softmax", name="stage_probabilities")(skip)
    model = keras.Model(inputs=inputs, outputs=outputs, name="ProgMapSkip2048x128")

    schedule = keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate=learning_rate,
        decay_steps=decay_steps,
        decay_rate=decay_rate,
    )
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=schedule),
        loss=keras.losses.CategoricalCrossentropy(),
        metrics=[
            keras.metrics.CategoricalAccuracy(name="accuracy"),
            keras.metrics.AUC(name="auc", multi_label=True, num_labels=3),
        ],
        weighted_metrics=[],
    )
    return model


def fit_model(
    model,
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_validation: np.ndarray,
    y_validation: np.ndarray,
    *,
    epochs: int = 1000,
    patience: int = 50,
    batch_size: int = 16,
    class_weight: dict[int, float] | None = None,
    monitor: str = "val_loss",
    min_delta: float = 1e-4,
    warmup_epochs: int = 20,
    balance_validation_loss: bool = True,
    progress_description: str | None = None,
    verbose: int = 1,
) -> FitResult:
    tf = require_tensorflow()
    if monitor not in EARLY_STOPPING_MONITORS:
        raise ValueError(f"monitor must be one of {EARLY_STOPPING_MONITORS}")
    y_train_one_hot = tf.keras.utils.to_categorical(y_train, num_classes=3)
    y_validation_one_hot = tf.keras.utils.to_categorical(y_validation, num_classes=3)
    callbacks: list[Any] = []
    if progress_description:
        callbacks.append(_epoch_progress_callback(progress_description, epochs))
    if patience > 0:
        callbacks.append(
            tf.keras.callbacks.EarlyStopping(
                monitor=monitor,
                mode="min" if monitor == "val_loss" else "max",
                patience=patience,
                min_delta=min_delta,
                restore_best_weights=True,
                start_from_epoch=warmup_epochs,
                verbose=0 if progress_description else verbose,
            )
        )
    validation_data: tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]
    if balance_validation_loss:
        validation_data = (
            x_validation,
            y_validation_one_hot,
            balanced_sample_weights(y_validation),
        )
    else:
        validation_data = (x_validation, y_validation_one_hot)
    history = model.fit(
        x_train,
        y_train_one_hot,
        validation_data=validation_data,
        epochs=epochs,
        batch_size=batch_size,
        callbacks=callbacks,
        class_weight=class_weight,
        shuffle=True,
        verbose=0 if progress_description else verbose,
    )
    serializable = {
        key: [float(value) for value in values]
        for key, values in history.history.items()
    }
    best_epoch, best_value = _monitor_summary(serializable, monitor, warmup_epochs)
    return FitResult(
        history=serializable,
        best_epoch=best_epoch,
        stopped_epoch=len(next(iter(serializable.values()), [])),
        best_monitor_value=best_value,
    )


def fit_fixed_epochs(
    model,
    x_train: np.ndarray,
    y_train: np.ndarray,
    *,
    epochs: int,
    batch_size: int = 16,
    class_weight: dict[int, float] | None = None,
    progress_description: str | None = None,
    verbose: int = 1,
) -> dict[str, list[float]]:
    tf = require_tensorflow()
    y_train_one_hot = tf.keras.utils.to_categorical(y_train, num_classes=3)
    callbacks: list[Any] = []
    if progress_description:
        callbacks.append(_epoch_progress_callback(progress_description, epochs))
    history = model.fit(
        x_train,
        y_train_one_hot,
        epochs=epochs,
        batch_size=batch_size,
        callbacks=callbacks,
        class_weight=class_weight,
        shuffle=True,
        verbose=0 if progress_description else verbose,
    )
    return {
        key: [float(value) for value in values]
        for key, values in history.history.items()
    }
