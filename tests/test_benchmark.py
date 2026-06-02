"""Tests for benchmark framework."""

import numpy as np
import polars as pl
import pytest

from skadi.benchmark_thresholds import (
    evaluate_assignment,
    grid_search_thresholds,
    BenchmarkResult,
)


class TestEvaluateAssignment:
    def test_perfect_prediction(self):
        predicted = pl.DataFrame({
            "SequenceID": ["seq1", "seq2", "seq3"],
            "rank_taxid": [100, 200, 300],
        })
        ground_truth = pl.DataFrame({
            "SequenceID": ["seq1", "seq2", "seq3"],
            "rank_taxid": [100, 200, 300],
        })
        result = evaluate_assignment(predicted, ground_truth, "species")
        assert result.precision == 1.0
        assert result.recall == 1.0
        assert result.f1 == 1.0
        assert result.accuracy == 1.0
        assert result.true_positives == 3

    def test_all_wrong(self):
        predicted = pl.DataFrame({
            "SequenceID": ["seq1", "seq2"],
            "rank_taxid": [100, 200],
        })
        ground_truth = pl.DataFrame({
            "SequenceID": ["seq1", "seq2"],
            "rank_taxid": [200, 100],
        })
        result = evaluate_assignment(predicted, ground_truth, "species")
        assert result.precision == 0.0
        assert result.recall == 0.0
        assert result.f1 == 0.0
        assert result.false_positives == 2

    def test_partial_assignment(self):
        predicted = pl.DataFrame({
            "SequenceID": ["seq1", "seq2"],
            "rank_taxid": [100, None],
        })
        ground_truth = pl.DataFrame({
            "SequenceID": ["seq1", "seq2"],
            "rank_taxid": [100, 200],
        })
        result = evaluate_assignment(predicted, ground_truth, "species")
        assert result.precision == 1.0  # No false positives
        assert result.recall == 0.5     # One false negative
        assert result.true_positives == 1
        assert result.false_negatives == 1


class TestGridSearchThresholds:
    def test_finds_optimal_threshold(self):
        # Scores: higher = more likely same rank
        scores = np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3])
        labels = np.array([1, 1, 1, 0, 0, 0, 0])  # First 3 are same rank
        thresholds = np.arange(0.3, 0.9, 0.1)

        best_thresh, best_f1 = grid_search_thresholds(
            scores, labels, thresholds, metric="f1"
        )
        # Optimal should be between 0.6 and 0.7
        assert 0.5 <= best_thresh <= 0.7
        assert best_f1 > 0.8

    def test_perfect_separation(self):
        scores = np.array([0.9, 0.8, 0.1, 0.05])
        labels = np.array([1, 1, 0, 0])
        thresholds = np.arange(0.05, 0.9, 0.05)

        best_thresh, best_f1 = grid_search_thresholds(scores, labels, thresholds)
        assert best_f1 == 1.0
        assert 0.1 <= best_thresh <= 0.8
