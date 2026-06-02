import numpy as np
import pytest
import polars as pl

from skadi.utils import (
    compute_cov,
    merge_intervals,
    merge_intervals_with_cutoff,
)


class TestMergeIntervals:
    def test_empty(self):
        arr = np.empty((0, 2), dtype=int)
        result = merge_intervals(arr)
        assert result.shape == (0, 2)

    def test_single_interval(self):
        arr = np.array([[1, 5]])
        result = merge_intervals(arr)
        np.testing.assert_array_equal(result, [[1, 5]])

    def test_overlapping(self):
        arr = np.array([[1, 5], [3, 7], [8, 10]])
        result = merge_intervals(arr)
        np.testing.assert_array_equal(result, [[1, 7], [8, 10]])

    def test_touching_merge(self):
        arr = np.array([[1, 5], [5, 10]])
        result = merge_intervals(arr, merge_touches=True)
        np.testing.assert_array_equal(result, [[1, 10]])

    def test_touching_no_merge(self):
        arr = np.array([[1, 5], [5, 10]])
        result = merge_intervals(arr, merge_touches=False)
        np.testing.assert_array_equal(result, [[1, 5], [5, 10]])

    def test_unsorted_input(self):
        arr = np.array([[5, 10], [1, 3], [2, 6]])
        result = merge_intervals(arr)
        np.testing.assert_array_equal(result, [[1, 10]])

    def test_reversed_intervals_swap(self):
        arr = np.array([[5, 1], [10, 3]])
        result = merge_intervals(arr, normalize="swap")
        np.testing.assert_array_equal(result, [[1, 10]])

    def test_reversed_intervals_error(self):
        arr = np.array([[5, 1]])
        with pytest.raises(ValueError, match="start > end"):
            merge_intervals(arr, normalize="error")

    def test_nan_rejection(self):
        arr = np.array([[1.0, np.nan]])
        with pytest.raises(ValueError, match="NaN"):
            merge_intervals(arr)


class TestMergeIntervalsWithCutoff:
    def test_no_merge_above_cutoff(self):
        intervals = np.array([[1, 3], [5, 7]])
        result = merge_intervals_with_cutoff(intervals, cutoff=0)
        np.testing.assert_array_equal(result, [[1, 3], [5, 7]])

    def test_merge_with_cutoff(self):
        intervals = np.array([[1, 3], [4, 7]])
        result = merge_intervals_with_cutoff(intervals, cutoff=1)
        np.testing.assert_array_equal(result, [[1, 7]])

    def test_empty(self):
        intervals = np.array([]).reshape(0, 2)
        result = merge_intervals_with_cutoff(intervals)
        assert result.shape == (0, 2)


class TestComputeCov:
    def test_full_coverage(self):
        qstart = pl.Series("qstart", [1, 100])
        qend = pl.Series("qend", [50, 200])
        qlen = pl.Series("qlen", [200])
        result = compute_cov([qstart, qend, qlen])
        # Merged intervals: [1,50] and [100,200] = total covered 149
        # 149/200 = 0.745, rounded to 0.74
        assert result == 0.74

    def test_no_overlap(self):
        qstart = pl.Series("qstart", [1, 100])
        qend = pl.Series("qend", [50, 200])
        qlen = pl.Series("qlen", [200])
        result = compute_cov([qstart, qend, qlen])
        assert isinstance(result, float)
