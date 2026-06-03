"""Tests for benchmark stratification and report generation."""

import numpy as np
import polars as pl
import pytest

from skadi.benchmark_thresholds import (
    stratify_by_length,
    evaluate_stratified,
    generate_comparison_report,
    DEFAULT_LENGTH_BINS,
)


class TestStratifyByLength:
    def test_bins_correctly(self):
        df = pl.DataFrame({
            "SequenceID": ["s1", "s2", "s3", "s4"],
            "Seqlen": [500, 3000, 7000, 15000],
        })
        result = stratify_by_length(df, length_col="Seqlen")
        bins = result["length_bin"].to_list()
        assert bins[0] == "<1kb"
        assert bins[1] == "1-5kb"
        assert bins[2] == "5-10kb"
        assert bins[3] == ">10kb"

    def test_custom_bins(self):
        df = pl.DataFrame({
            "SequenceID": ["s1", "s2"],
            "Seqlen": [100, 500],
        })
        custom_bins = [(0, 200, "small"), (200, float("inf"), "large")]
        result = stratify_by_length(df, length_col="Seqlen", bins=custom_bins)
        bins = result["length_bin"].to_list()
        assert bins[0] == "small"
        assert bins[1] == "large"

    def test_empty_dataframe(self):
        df = pl.DataFrame({"SequenceID": [], "Seqlen": []})
        result = stratify_by_length(df, length_col="Seqlen")
        assert result.is_empty()


class TestEvaluateStratified:
    def test_stratified_evaluation(self):
        predicted = pl.DataFrame({
            "SequenceID": ["seq1", "seq2", "seq3"],
            "rank_taxid": [100, 200, 300],
            "Seqlen": [500, 3000, 15000],
        })
        ground_truth = pl.DataFrame({
            "SequenceID": ["seq1", "seq2", "seq3"],
            "rank_taxid": [100, 200, 400],
            "Seqlen": [500, 3000, 15000],
        })

        result = evaluate_stratified(
            predicted, ground_truth, ranks=["species"], length_col="Seqlen"
        )

        assert len(result) == 3  # 3 bins
        # seq1: <1kb, correct (TP)
        # seq2: 1-5kb, correct (TP)
        # seq3: >10kb, wrong (FP)
        row_1kb = result.filter(pl.col("length_bin") == "<1kb").to_dicts()[0]
        assert row_1kb["precision"] == 1.0
        assert row_1kb["recall"] == 1.0

        row_10kb = result.filter(pl.col("length_bin") == ">10kb").to_dicts()[0]
        assert row_10kb["precision"] == 0.0
        assert row_10kb["recall"] == 0.0


class TestGenerateComparisonReport:
    def test_method_comparison_report(self):
        results = pl.DataFrame({
            "method": ["cascade", "cascade", "consensus", "consensus"],
            "rank": ["species", "genus", "species", "genus"],
            "precision": [0.9, 0.8, 0.95, 0.85],
            "recall": [0.85, 0.75, 0.9, 0.8],
            "f1": [0.87, 0.77, 0.92, 0.82],
            "accuracy": [0.95, 0.9, 0.97, 0.93],
        })
        report = generate_comparison_report(results)
        assert "# SKADI Benchmark Report" in report
        assert "cascade" in report
        assert "consensus" in report
        assert "species" in report

    def test_stratified_report(self):
        results = pl.DataFrame({
            "rank": ["species", "species"],
            "length_bin": ["<1kb", ">10kb"],
            "precision": [0.9, 0.7],
            "recall": [0.85, 0.65],
            "f1": [0.87, 0.67],
            "accuracy": [0.95, 0.85],
            "n": [100, 50],
        })
        report = generate_comparison_report(results)
        assert "<1kb" in report
        assert ">10kb" in report
