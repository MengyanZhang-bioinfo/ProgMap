import pytest

from progmap.cli import build_parser


def test_public_version_and_generic_help(capsys):
    parser = build_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--version"])
    assert exc_info.value.code == 0
    assert capsys.readouterr().out.strip() == "progmap v0.1.0"

    help_text = parser.format_help()
    assert "Input root containing one directory per dataset" in help_text
    assert "Comma-separated dataset directory names" in help_text
    assert "Nested CV, one holdout split, or fixed epochs" in help_text


def test_defaults_and_configurable_options():
    parser = build_parser()
    defaults = parser.parse_args(["--data-root", "/data/progmap_inputs"])
    assert defaults.output == "progmap_results"
    assert defaults.folds == 3
    assert defaults.inner_folds == 3
    assert defaults.early_stopping == "nested"
    assert defaults.early_stopping_monitor == "val_loss"
    assert defaults.epochs == 1000
    assert defaults.patience == 50
    assert defaults.learning_rate == 1e-3
    assert defaults.seed == 42
    assert defaults.correlation_method == "pearson"
    assert defaults.test == "ttest"
    assert defaults.top_n == "significant"

    custom = parser.parse_args(
        [
            "--data-root",
            "/data/progmap_inputs",
            "--output",
            "/results/custom",
            "--folds",
            "4",
            "--inner-folds",
            "2",
            "--early-stopping",
            "holdout",
            "--early-stopping-monitor",
            "val_auc",
            "--epochs",
            "250",
            "--patience",
            "25",
            "--learning-rate",
            "0.0005",
            "--seed",
            "7",
            "--correlation-method",
            "spearman",
            "--test",
            "wilcoxon",
            "--top-n",
            "100",
        ]
    )
    assert custom.output == "/results/custom"
    assert custom.folds == 4
    assert custom.inner_folds == 2
    assert custom.early_stopping == "holdout"
    assert custom.early_stopping_monitor == "val_auc"
    assert custom.epochs == 250
    assert custom.patience == 25
    assert custom.learning_rate == 0.0005
    assert custom.seed == 7
    assert custom.correlation_method == "spearman"
    assert custom.test == "wilcoxon"
    assert custom.top_n == 100
