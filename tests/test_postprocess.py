"""Tests for postprocess.py consensus and cascade integration."""

import polars as pl
import pytest
from click.testing import CliRunner

from scripts.postprocess import _build_rank_lookup, _build_lineage_lookup, _apply_taxon_columns


class TestTaxonLookups:
    def test_build_rank_lookup_empty(self):
        result = _build_rank_lookup(None, [], "species")
        assert result == {}

    def test_build_lineage_lookup_empty(self):
        result = _build_lineage_lookup(None, [])
        assert result == {}


class TestPostprocessCLI:
    def test_help_includes_method_option(self):
        from scripts.postprocess import main
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "--method" in result.output
        assert "cascade" in result.output
        assert "consensus" in result.output

    def test_help_includes_min_confidence(self):
        from scripts.postprocess import main
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "--min-confidence" in result.output


class TestConsensusIntegration:
    def test_consensus_output_format(self):
        """Test that consensus output has expected columns."""
        from skadi.consensus import build_consensus_assignment

        # Mock taxdb
        class MockTaxon:
            def __init__(self, taxid, taxdb=None):
                self.taxid = taxid
                self.rank_taxid_dictionary = {
                    "species": taxid,
                    "genus": taxid,
                    "family": taxid,
                    "order": taxid,
                    "class": taxid,
                    "phylum": taxid,
                    "kingdom": taxid,
                }

        import taxopy
        original_taxon = taxopy.Taxon
        taxopy.Taxon = MockTaxon

        try:
            ani_df = pl.DataFrame({
                "query": ["seq1", "seq1", "seq2"],
                "taxid": [1, 2, 3],
                "tani": [0.95, 0.85, 0.92],
                "ani": [0.95, 0.85, 0.92],
            })

            result = build_consensus_assignment(
                ani_df=ani_df,
                aai_df=None,
                api_df=None,
                thresholds={"species": 0.81, "genus": 0.49},
                taxdb="mock",
                min_confidence=0.5,
            )

            assert "SequenceID" in result.columns
            assert "rank" in result.columns
            assert "rank_taxid" in result.columns
            assert "confidence" in result.columns
            assert "methods" in result.columns
            assert "taxlineage" in result.columns
            assert "Score" in result.columns
            assert "level" in result.columns
            assert len(result) == 2  # seq1 and seq2 (best hit each)
        finally:
            taxopy.Taxon = original_taxon

    def test_threshold_resolution_legacy_format(self):
        """Test that legacy threshold keys like 'tanis' are resolved."""
        from skadi.consensus import _resolve_threshold

        thresholds = {"tanis": 0.81, "tanig": 0.49, "species": 0.75}
        assert _resolve_threshold(thresholds, "ani", "species") == 0.81
        assert _resolve_threshold(thresholds, "ani", "genus") == 0.49

    def test_threshold_resolution_modern_format(self):
        """Test modern threshold keys like 'tani_species'."""
        from skadi.consensus import _resolve_threshold

        thresholds = {"tani_species": 0.85, "species": 0.75}
        assert _resolve_threshold(thresholds, "ani", "species") == 0.85

    def test_threshold_resolution_fallback(self):
        """Test fallback when no matching key exists."""
        from skadi.consensus import _resolve_threshold

        thresholds = {"other": 0.5}
        assert _resolve_threshold(thresholds, "ani", "species") == 0.3  # default
