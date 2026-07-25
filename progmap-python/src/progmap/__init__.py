"""ProgMap: leakage-safe multi-omics stage classification and feature discovery."""

from .data import MultiOmicsData, load_cancer
from .preprocessing import FoldPreprocessor
from .statistics import rank_features

__all__ = ["MultiOmicsData", "FoldPreprocessor", "load_cancer", "rank_features"]
__version__ = "0.1.0"
