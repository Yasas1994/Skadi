"""Tests for consensus taxonomy assignment."""

import numpy as np
import polars as pl
import pytest

from skadi.consensus import (
    _score_to_confidence,
    build_consensus_assignment,
    DEFAULT_WEIGHTS,
    RANKS,
)


class TestScoreToConfidence:
    def test_at_threshold(self):
        # At score == threshold, confidence should be ~0.5
        conf = _score_to_confidence(0.5, threshold=0.5, steepness=10.0)
        assert abs(conf - 0.5) < 0.01

    def test_above_threshold(self):
        conf = _score_to_confidence(0.8, threshold=0.5, steepness=10.0)
        assert conf > 0.9

    def test_below_threshold(self):
        conf = _score_to_confidence(0.2, threshold=0.5, steepness=10.0)
        assert conf < 0.1

    def test_exactly_zero(self):
        conf = _score_to_confidence(0.0, threshold=0.5)
        assert conf < 0.01


class TestBuildConsensus:
    def test_empty_inputs(self):
        result = build_consensus_assignment(
            ani_df=None,
            aai_df=None,
            api_df=None,
            thresholds={"species": 0.81, "genus": 0.49},
            taxdb=None,
        )
        assert result.is_empty()

    def test_ani_only_wins_when_strong(self):
        """When ANI has high confidence at species, it should win."""
        ani_df = pl.DataFrame({
            "query": ["seq1"],
            "taxid": [1],
            "tani": [0.95],
            "ani": [0.95],
            "qcov": [1.0],
        })

        # Mock taxdb that returns taxid itself for any rank
        class MockTaxon:
            def __init__(self, taxid, taxdb):
                self.taxid = taxid
                self.rank_taxid_dictionary = {"species": taxid, "genus": taxid}

        import taxopy
        original_taxon = taxopy.Taxon
        taxopy.Taxon = MockTaxon

        try:
            result = build_consensus_assignment(
                ani_df=ani_df,
                aai_df=None,
                api_df=None,
                thresholds={"species": 0.81, "genus": 0.49},
                taxdb="mock",
                min_confidence=0.5,
            )
            assert not result.is_empty()
            # With tani=0.95 and threshold=0.81, confidence is very high
            # The best rank should be species since it has highest confidence
            assert result["confidence"][0] > 0.5
        finally:
            taxopy.Taxon = original_taxon

    def test_weights_sum_to_one_per_rank(self):
        for rank in RANKS:
            total = sum(DEFAULT_WEIGHTS[m].get(rank, 0) for m in DEFAULT_WEIGHTS)
            assert abs(total - 1.0) < 0.01, f"Weights for {rank} sum to {total}"
