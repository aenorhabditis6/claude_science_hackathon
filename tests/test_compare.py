"""Sanity checks for the distance math in compare.py, the statistical core of ShiftScope.

CPU-only, no network, no API key. Run with `pytest tests/` or `python tests/test_compare.py`.
These pin down the properties the whole tool relies on: the energy distance is zero for
identical clouds, symmetric, large and significant for separated clouds, and non-significant
when two samples come from the same distribution.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from shiftscope import compare  # noqa: E402


def test_edistance_zero_for_identical_clouds():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(200, 10))
    assert abs(compare.edistance(A, A)) < 1e-6


def test_edistance_symmetric():
    rng = np.random.default_rng(1)
    A = rng.normal(size=(150, 8))
    B = rng.normal(loc=1.0, size=(150, 8))
    assert abs(compare.edistance(A, B) - compare.edistance(B, A)) < 1e-6


def test_separated_clouds_are_large_and_significant():
    rng = np.random.default_rng(2)
    A = rng.normal(loc=0.0, size=(120, 6))
    B = rng.normal(loc=5.0, size=(120, 6))  # well separated
    e, p = compare.etest(A, B, n_perm=200, seed=0)
    assert e > 0
    assert p <= 1.0 / (200 + 1) + 1e-9  # no permutation beats a huge real separation


def test_same_distribution_is_not_significant():
    rng = np.random.default_rng(3)
    A = rng.normal(size=(120, 6))
    B = rng.normal(size=(120, 6))  # drawn from the same distribution
    _, p = compare.etest(A, B, n_perm=200, seed=0)
    assert p > 0.05


def test_mmd_nonnegative_and_zero_for_identical():
    rng = np.random.default_rng(4)
    A = rng.normal(size=(100, 5))
    assert compare.mmd(A, A) >= -1e-9
    assert abs(compare.mmd(A, A)) < 1e-6


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok  ", name)
    print("all passed")
