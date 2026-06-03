"""Benchmark framework for comparing taxonomy assignment methods.

This module provides tools to:
1. Simulate taxonomy assignments with different threshold strategies
2. Compute accuracy metrics (precision, recall, F1) per rank
3. Compare: fixed thresholds, confidence-based, consensus voting, symmetric coverage
4. Stratify results by fragment length
5. Run leave-one-out benchmarks
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Callable, Tuple
from pathlib import Path
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
) -> Tuple[float, float]:
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


# ---------------------------------------------------------------------------
# Fragment-length stratification
# ---------------------------------------------------------------------------

DEFAULT_LENGTH_BINS = [
    (0, 1000, "<1kb"),
    (1000, 5000, "1-5kb"),
    (5000, 10000, "5-10kb"),
    (10000, float("inf"), ">10kb"),
]


def stratify_by_length(
    df: pl.DataFrame,
    length_col: str = "Seqlen",
    bins: Optional[List[Tuple[float, float, str]]] = None,
) -> pl.DataFrame:
    """Add a length_bin column to a DataFrame based on sequence lengths.

    Args:
        df: DataFrame with a length column.
        length_col: Name of the length column.
        bins: List of (min, max, label) tuples. Default: <1kb, 1-5kb, 5-10kb, >10kb.

    Returns:
        DataFrame with added 'length_bin' column.
    """
    bins = bins or DEFAULT_LENGTH_BINS

    expr = pl.lit("unknown")
    for min_len, max_len, label in bins:
        if max_len == float("inf"):
            expr = pl.when(pl.col(length_col) >= min_len).then(pl.lit(label)).otherwise(expr)
        else:
            expr = pl.when(
                (pl.col(length_col) >= min_len) & (pl.col(length_col) < max_len)
            ).then(pl.lit(label)).otherwise(expr)

    return df.with_columns(expr.alias("length_bin"))


def evaluate_stratified(
    predicted: pl.DataFrame,
    ground_truth: pl.DataFrame,
    ranks: List[str],
    length_col: str = "Seqlen",
    bins: Optional[List[Tuple[float, float, str]]] = None,
) -> pl.DataFrame:
    """Evaluate assignments stratified by fragment length.

    Args:
        predicted: DataFrame with predictions and length column.
        ground_truth: DataFrame with ground truth and length column.
        ranks: List of ranks to evaluate.
        length_col: Column name for sequence lengths.
        bins: Length bin definitions.

    Returns:
        DataFrame with columns [rank, length_bin, precision, recall, f1, accuracy, n].
    """
    pred_binned = stratify_by_length(predicted, length_col, bins)
    gt_binned = stratify_by_length(ground_truth, length_col, bins)

    rows = []
    all_bins = pred_binned["length_bin"].unique().to_list()

    for rank in ranks:
        for bin_label in all_bins:
            pred_bin = pred_binned.filter(pl.col("length_bin") == bin_label)
            gt_bin = gt_binned.filter(pl.col("length_bin") == bin_label)

            if pred_bin.is_empty() and gt_bin.is_empty():
                continue

            result = evaluate_assignment(pred_bin, gt_bin, rank)
            rows.append({
                "rank": rank,
                "length_bin": bin_label,
                "precision": result.precision,
                "recall": result.recall,
                "f1": result.f1,
                "accuracy": result.accuracy,
                "n": pred_bin.height + gt_bin.height,
            })

    return pl.DataFrame(rows)


# ---------------------------------------------------------------------------
# Leave-one-out benchmark runner
# ---------------------------------------------------------------------------

def run_leaveoneout_benchmark(
    results_dir: Path,
    ground_truth_file: Path,
    ranks: List[str] = None,
    length_col: str = "Seqlen",
) -> pl.DataFrame:
    """Run a leave-one-out benchmark from existing result files.

    Args:
        results_dir: Directory containing _leaveout_*.tsv result files.
        ground_truth_file: Path to ground truth TSV with columns [SequenceID, rank_taxid, ...].
        ranks: Ranks to evaluate. Default: ["species", "genus", "family"].
        length_col: Column for sequence length stratification.

    Returns:
        DataFrame with benchmark results per rank and length bin.
    """
    ranks = ranks or ["species", "genus", "family"]

    result_files = list(results_dir.glob("*_leaveout_*.tsv"))
    if not result_files:
        raise FileNotFoundError(f"No leave-out result files found in {results_dir}")

    ground_truth = pl.read_csv(ground_truth_file, separator="\t")

    all_results = []
    for rf in result_files:
        predicted = pl.read_csv(rf, separator="\t")
        stratified = evaluate_stratified(predicted, ground_truth, ranks, length_col)
        stratified = stratified.with_columns(pl.lit(rf.name).alias("file"))
        all_results.append(stratified)

    return pl.concat(all_results) if all_results else pl.DataFrame()


# ---------------------------------------------------------------------------
# Comparison report generator
# ---------------------------------------------------------------------------

def generate_comparison_report(
    results: pl.DataFrame,
    output_path: Optional[Path] = None,
) -> str:
    """Generate a markdown comparison report from benchmark results.

    Args:
        results: DataFrame from compare_methods() or evaluate_stratified().
        output_path: If provided, write report to this file.

    Returns:
        Markdown-formatted report string.
    """
    lines = ["# SKADI Benchmark Report\n"]

    # Overall results
    if "method" in results.columns:
        lines.append("## Method Comparison\n")
        lines.append("| Method | Rank | Precision | Recall | F1 | Accuracy |")
        lines.append("|--------|------|-----------|--------|-----|----------|")
        for row in results.iter_rows(named=True):
            lines.append(
                f"| {row['method']} | {row['rank']} | "
                f"{row['precision']:.3f} | {row['recall']:.3f} | "
                f"{row['f1']:.3f} | {row['accuracy']:.3f} |"
            )
    else:
        lines.append("## Stratified Results\n")
        lines.append("| Rank | Length Bin | Precision | Recall | F1 | Accuracy | N |")
        lines.append("|------|------------|-----------|--------|-----|----------|---|")
        for row in results.iter_rows(named=True):
            lines.append(
                f"| {row['rank']} | {row['length_bin']} | "
                f"{row['precision']:.3f} | {row['recall']:.3f} | "
                f"{row['f1']:.3f} | {row['accuracy']:.3f} | {row.get('n', '-')} |"
            )

    report = "\n".join(lines) + "\n"

    if output_path:
        output_path.write_text(report)

    return report
