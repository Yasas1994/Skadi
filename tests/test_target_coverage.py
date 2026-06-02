"""Tests for target coverage and symmetric coverage computation."""

import numpy as np
import polars as pl
import pytest

from skadi.utils import compute_cov, compute_cov_bidirectional


class TestComputeCovBidirectional:
    def test_perfect_symmetry(self):
        """When query and target have identical coverage, scov = qcov = tcov."""
        qstart = pl.Series("qstart", [1, 100])
        qend = pl.Series("qend", [50, 200])
        qlen = pl.Series("qlen", [200])
        tstart = pl.Series("tstart", [1, 100])
        tend = pl.Series("tend", [50, 200])
        tlen = pl.Series("tlen", [200])

        qcov, tcov, scov = compute_cov_bidirectional(
            qstart, qend, qlen, tstart, tend, tlen
        )
        assert qcov == tcov
        assert scov == qcov

    def test_asymmetric_penalizes_small_query(self):
        """Small query aligning to large target should have low scov."""
        # Query: 100bp aligning fully
        qstart = pl.Series("qstart", [1])
        qend = pl.Series("qend", [100])
        qlen = pl.Series("qlen", [100])
        # Target: 10000bp, only 100bp covered
        tstart = pl.Series("tstart", [5000])
        tend = pl.Series("tend", [5100])
        tlen = pl.Series("tlen", [10000])

        qcov, tcov, scov = compute_cov_bidirectional(
            qstart, qend, qlen, tstart, tend, tlen
        )
        assert qcov == pytest.approx(0.99, abs=0.02)  # Near-full query coverage
        assert tcov == 0.01  # Only 1% of target covered
        assert scov < 0.15  # Geometric mean penalizes asymmetry
        assert scov == pytest.approx(0.1, abs=0.05)

    def test_zero_coverage(self):
        """Empty alignment should return zeros."""
        qstart = pl.Series("qstart", [])
        qend = pl.Series("qend", [])
        qlen = pl.Series("qlen", [100])
        tstart = pl.Series("tstart", [])
        tend = pl.Series("tend", [])
        tlen = pl.Series("tlen", [100])

        qcov, tcov, scov = compute_cov_bidirectional(
            qstart, qend, qlen, tstart, tend, tlen
        )
        assert qcov == 0.0
        assert tcov == 0.0
        assert scov == 0.0
