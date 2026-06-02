"""Benchmark framework for comparing taxonomy assignment methods.

This module provides tools to:
1. Simulate taxonomy assignments with different threshold strategies
2. Compute accuracy metrics (precision, recall, F1) per rank
3. Compare: fixed thresholds, confidence-based, consensus voting, symmetric coverage
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Callable
import numpy as np
import polars as pl


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""
    method_name: str
    rank: str
    precision: float
    recall: float
    f1: float
    accuracy: float
    true_positives: int
    false_positives: int
    false_negatives: int
    true_negatives: int


def evaluate_assignment(
    predicted: pl.DataFrame,
    ground_truth: pl.DataFrame,
    rank: str,
    id_col: str = "SequenceID",
    taxid_col: str = "rank_taxid",
) -> BenchmarkResult:
    """Evaluate taxonomy assignment accuracy at a given rank.

    Args:
        predicted: DataFrame with predicted assignments.
        ground_truth: DataFrame with true assignments.
        rank: Taxonomic rank to evaluate.
        id_col: Column name for sequence IDs.
        taxid_col: Column name for taxids.

    Returns:
        BenchmarkResult with precision, recall, F1, etc.
    """
    merged = predicted.join(ground_truth, on=id_col, how="full", suffix="_true")

    # Count outcomes
    tp = merged.filter(
        (pl.col(taxid_col) == pl.col(f"{taxid_col}_true"))
        & pl.col(taxid_col).is_not_null()
    ).height

    fp = merged.filter(
        (pl.col(taxid_col) != pl.col(f"{taxid_col}_true"))
        & pl.col(taxid_col).is_not_null()
    ).height

    fn = merged.filter(
        pl.col(taxid_col).is_null() & pl.col(f"{taxid_col}_true").is_not_null()
    ).height

    tn = merged.filter(
        pl.col(taxid_col).is_null() & pl.col(f"{taxid_col}_true").is_null()
    ).height

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    accuracy = (tp + tn) / (tp + fp + fn + tn) if (tp + fp + fn + tn) > 0 else 0.0

    return BenchmarkResult(
        method_name="",
        rank=rank,
        precision=round(precision, 4),
        recall=round(recall, 4),
        f1=round(f1, 4),
        accuracy=round(accuracy, 4),
        true_positives=tp,
        false_positives=fp,
        false_negatives=fn,
        true_negatives=tn,
    )


def grid_search_thresholds(
    scores: np.ndarray,
    labels: np.ndarray,
    thresholds: np.ndarray,
    metric: str = "f1",
) -> tuple[float, float]:
    """Find optimal threshold by grid search.

    Args:
        scores: Array of scores (e.g., tani values).
        labels: Binary array (1 = same rank, 0 = different rank).
        thresholds: Array of threshold candidates.
        metric: Metric to optimize ("f1", "precision", "recall", "accuracy").

    Returns:
        (best_threshold, best_score)
    """
    best_thresh = thresholds[0]
    best_score = 0.0

    for thresh in thresholds:
        pred = (scores >= thresh).astype(int)
        tp = np.sum((pred == 1) & (labels == 1))
        fp = np.sum((pred == 1) & (labels == 0))
        fn = np.sum((pred == 0) & (labels == 1))
        tn = np.sum((pred == 0) & (labels == 0))

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        accuracy = (tp + tn) / (tp + fp + fn + tn) if (tp + fp + fn + tn) > 0 else 0

        score = {"f1": f1, "precision": precision, "recall": recall, "accuracy": accuracy}[metric]
        if score > best_score:
            best_score = score
            best_thresh = thresh

    return float(best_thresh), float(best_score)


def compare_methods(
    methods: Dict[str, Callable[[], pl.DataFrame]],
    ground_truth: pl.DataFrame,
    ranks: List[str],
) -> pl.DataFrame:
    """Compare multiple assignment methods against ground truth.

    Args:
        methods: Dict of method_name → function that returns predictions.
        ground_truth: DataFrame with true assignments.
        ranks: List of ranks to evaluate.

    Returns:
        DataFrame with columns [method, rank, precision, recall, f1, accuracy].
    """
    rows = []
    for method_name, predict_fn in methods.items():
        predicted = predict_fn()
        for rank in ranks:
            result = evaluate_assignment(predicted, ground_truth, rank)
            result.method_name = method_name
            rows.append({
                "method": method_name,
                "rank": rank,
                "precision": result.precision,
                "recall": result.recall,
                "f1": result.f1,
                "accuracy": result.accuracy,
                "tp": result.true_positives,
                "fp": result.false_positives,
                "fn": result.false_negatives,
            })

    return pl.DataFrame(rows)
