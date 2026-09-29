"""Tests for the sanity helper: one-to-one (Hungarian) matched accuracy."""
import numpy as np
import pytest

import sanity


def test_hungarian_perfect_prediction():
    true = np.array([1, 1, 2, 2, 3, 3, 0])
    pred = np.array([7, 7, 5, 5, 9, 9, 0])       # relabelled, plus noise==noise
    res = sanity.hungarian_match_accuracy(pred, true, np.ones(7, bool))
    assert res["accuracy"] == 1.0
    assert res["n_predicted_labels"] == 4
    assert res["n_true_labels"] == 4


def test_hungarian_penalises_split():
    # one true cell split across two predicted labels: at most the majority
    # mass of the true cell can be recovered, never both halves
    true = np.array([1] * 10 + [2] * 10)
    pred = np.array([4] * 6 + [5] * 4 + [6] * 10)
    # best injective map: 4->1 (6 mol), 6->2 (10 mol); label 5 unmatched
    res = sanity.hungarian_match_accuracy(pred, true, np.ones(20, bool))
    assert res["accuracy"] == pytest.approx(16 / 20)


def test_hungarian_never_exceeds_many_to_one_upper_bound():
    # merging two true cells into one predicted label: at most one is credited
    true = np.array([1] * 10 + [2] * 10)
    pred = np.array([4] * 20)
    res = sanity.hungarian_match_accuracy(pred, true, np.ones(20, bool))
    assert res["accuracy"] == pytest.approx(10 / 20)


def test_hungarian_noise_semantics():
    # a cell prediction on background molecules is never correct, and noise
    # predictions are only correct on true background
    true = np.array([1] * 4 + [0] * 4)
    pred = np.array([3] * 4 + [3] * 4)      # cell 3 covers cell 1 and bg
    res = sanity.hungarian_match_accuracy(pred, true, np.ones(8, bool))
    assert res["accuracy"] == pytest.approx(4 / 8)


def test_hungarian_respects_mask_and_rejects_empty():
    true = np.array([1, 1, 2, 2])
    pred = np.array([1, 1, 2, 2])
    mask = np.array([True, False, True, False])
    res = sanity.hungarian_match_accuracy(pred, true, mask)
    assert res["accuracy"] == 1.0 and res["n_evaluated"] == 2
    with pytest.raises(ValueError, match="empty"):
        sanity.hungarian_match_accuracy(pred, true, np.zeros(4, bool))
