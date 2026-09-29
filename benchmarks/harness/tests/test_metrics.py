"""Hand-made unit tests for the pure metric functions.

Cases: perfect segmentation, one split, one merge, all noise, permuted
labels, plus real-vs-baseline comparisons and a concrete comparison with
st-recoverability's one-to-one ``matched_accuracy`` definition.
"""
import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

import metrics as m


# ---------------------------------------------------------------------------
# hand-made label cases
# ---------------------------------------------------------------------------

def base_truth():
    """10+8+10+10 molecules in true cells 1..4 plus 5 noise molecules."""
    truth = np.array([1] * 10 + [2] * 8 + [3] * 10 + [4] * 10 + [0] * 5,
                     dtype=np.int64)
    return truth


def interior_mask(truth, with_noise=True):
    mask = np.ones(len(truth), dtype=bool)
    if not with_noise:
        mask = truth > 0
    return mask


def test_perfect_segmentation():
    truth = base_truth()
    pred = truth.copy()
    mask = interior_mask(truth)
    assert m.matched_accuracy(pred, truth, mask) == pytest.approx(1.0)
    assert m.accuracy_1to1(pred, truth, mask) == pytest.approx(1.0)
    assert m.ari(pred, truth, mask) == pytest.approx(1.0)
    assert m.ari_assigned(pred, truth, mask) == pytest.approx(1.0)
    assert m.ami(pred, truth, mask) == pytest.approx(1.0)
    assert m.ami_assigned(pred, truth, mask) == pytest.approx(1.0)
    # only 5 noise molecules -> below NOISE_MIN_COUNT, skipped by default
    assert np.isnan(m.noise_precision(pred, truth, mask))
    assert np.isnan(m.noise_recall(pred, truth, mask))
    assert m.noise_precision(pred, truth, mask, min_count=0) == pytest.approx(1.0)
    assert m.noise_recall(pred, truth, mask, min_count=0) == pytest.approx(1.0)
    assert m.cell_count_ratio(pred, truth, mask) == pytest.approx(1.0)
    assert m.over_segmentation_rate(pred, truth, mask) == pytest.approx(0.0)
    assert m.under_segmentation_rate(pred, truth, mask) == pytest.approx(0.0)
    assert m.recovery_rate(pred, truth, mask) == pytest.approx(1.0)
    assert m.median_matched_jaccard(pred, truth, mask) == pytest.approx(1.0)
    assert m.assigned_fraction(pred, mask) == pytest.approx(38 / 43)
    assert m.assigned_fraction(truth, mask) == pytest.approx(38 / 43)
    assert m.oracle_gap(1.0, 0.99) == pytest.approx(-0.01)
    assert np.isnan(m.oracle_gap(1.0, None))


def test_one_split():
    """True cell 1 split into two predicted parts (6 and 4 molecules)."""
    truth = base_truth()
    pred = truth.copy()
    pred[:6] = 10          # part A of true cell 1
    pred[6:10] = 11        # part B of true cell 1
    mask = interior_mask(truth)
    # many-to-one mapping: both parts map back to true cell 1 -> still perfect
    assert m.matched_accuracy(pred, truth, mask) == pytest.approx(1.0)
    # one-to-one: only one part can be matched (6 of 10 molecules of cell 1)
    assert m.accuracy_1to1(pred, truth, mask) == pytest.approx(39 / 43)
    # cell 1 is over-split: largest part holds 6/10 < 80%
    assert m.over_segmentation_rate(pred, truth, mask) == pytest.approx(0.25)
    # each part is internally pure -> no under-segmentation
    assert m.under_segmentation_rate(pred, truth, mask) == pytest.approx(0.0)
    # best Jaccard for cell 1 = 6/10 >= 0.5 (part is contained in the cell)
    assert m.recovery_rate(pred, truth, mask) == pytest.approx(1.0)
    jac = m.median_matched_jaccard(pred, truth, mask)
    assert jac == pytest.approx(np.median([6 / 10, 1.0, 1.0, 1.0]))
    assert m.cell_count_ratio(pred, truth, mask) == pytest.approx(5 / 4)
    assert m.ari(pred, truth, mask) < 1.0


def test_one_merge():
    """Predicted cell A merges true cells 1 (10) and 2 (8)."""
    truth = base_truth()
    pred = truth.copy()
    pred[truth != 0] = np.where(truth[truth != 0] <= 2,
                                7, truth[truth != 0])
    mask = interior_mask(truth)
    # A maps to true cell 1 (majority 10 vs 8); true-cell-2 molecules are wrong
    expected = (10 + 10 + 10 + 5) / len(truth)
    assert m.matched_accuracy(pred, truth, mask) == pytest.approx(expected)
    assert m.under_segmentation_rate(pred, truth, mask) == pytest.approx(1 / 3)
    assert m.over_segmentation_rate(pred, truth, mask) == pytest.approx(0.0)
    # Jaccards of the merged cell: 10/18 and 8/18
    assert m.recovery_rate(pred, truth, mask) == pytest.approx(0.75)
    jac = m.median_matched_jaccard(pred, truth, mask)
    assert jac == pytest.approx(np.median([10 / 18, 8 / 18, 1.0, 1.0]))
    assert m.cell_count_ratio(pred, truth, mask) == pytest.approx(3 / 4)
    # noise stays noise (but only 5 noise molecules -> NaN by default)
    assert np.isnan(m.noise_precision(pred, truth, mask))
    assert m.noise_precision(pred, truth, mask, min_count=0) == pytest.approx(1.0)
    assert m.noise_recall(pred, truth, mask, min_count=0) == pytest.approx(1.0)
    # one-to-one: cell 7 takes true cell 1 (10), cells 3+4 stay, noise stays
    assert m.accuracy_1to1(pred, truth, mask) == pytest.approx(35 / 43)


def test_all_noise():
    truth = base_truth()
    pred = np.zeros_like(truth)
    mask = interior_mask(truth)
    assert m.matched_accuracy(pred, truth, mask) == pytest.approx(5 / 43)
    assert m.accuracy_1to1(pred, truth, mask) == pytest.approx(5 / 43)
    # tiny noise counts -> skipped by default, defined with min_count=0
    assert np.isnan(m.noise_precision(pred, truth, mask))
    assert m.noise_precision(pred, truth, mask, min_count=0) == pytest.approx(5 / 43)
    assert m.noise_recall(pred, truth, mask, min_count=0) == pytest.approx(1.0)
    assert m.ari(pred, truth, mask) == pytest.approx(0.0)
    # no molecule assigned on both sides -> assigned-only ARI is undefined
    assert np.isnan(m.ari_assigned(pred, truth, mask))
    assert m.cell_count_ratio(pred, truth, mask) == pytest.approx(0.0)
    # no predicted cells at all -> under-segmentation undefined
    assert np.isnan(m.under_segmentation_rate(pred, truth, mask))
    # every true cell's only "part" is noise (100%) -> not counted as a split
    assert m.over_segmentation_rate(pred, truth, mask) == pytest.approx(0.0)
    assert m.recovery_rate(pred, truth, mask) == pytest.approx(0.0)
    assert m.median_matched_jaccard(pred, truth, mask) == pytest.approx(0.0)


def test_noise_min_count_threshold():
    """noise precision/recall are NaN below NOISE_MIN_COUNT, defined above."""
    n_noise = m.NOISE_MIN_COUNT
    truth = np.array([0] * n_noise + [1] * 100, dtype=np.int64)
    pred = truth.copy()
    assert m.noise_precision(pred, truth) == pytest.approx(1.0)
    assert m.noise_recall(pred, truth) == pytest.approx(1.0)
    # one molecule short on the predicted-noise side -> precision undefined,
    # recall still defined (49 of 50 true-noise molecules called noise)
    pred2 = truth.copy()
    pred2[0] = 1
    assert np.isnan(m.noise_precision(pred2, truth))
    assert m.noise_recall(pred2, truth) == pytest.approx(49 / 50)
    # no true noise at all (st-recoverability style data) -> both skipped
    clean_t = np.array([1] * 100, dtype=np.int64)
    assert np.isnan(m.noise_recall(clean_t, clean_t))
    assert np.isnan(m.noise_precision(clean_t, clean_t))


def test_majority_noise_cell_gets_no_noise_credit():
    """A real predicted cell whose majority is true noise earns no credit
    for its noise molecules."""
    truth = np.array([1] * 12 + [0] * 20, dtype=np.int64)
    # cell 5 = 2 molecules of cell 1 + 8 noise molecules -> majority noise
    pred = np.array([5, 5] + [1] * 10 + [0] * 12 + [5] * 8, dtype=np.int64)
    # cell 5 gets no credit; pred noise 12/12 correct; cell 1: 10/10 correct
    expected = (10 + 12) / 32
    assert m.matched_accuracy(pred, truth) == pytest.approx(expected)
    # one-to-one agrees: cell 5's noise molecules are not credited
    assert m.accuracy_1to1(pred, truth) == pytest.approx(expected)


def test_all_noise_precision_vacuous():
    truth = np.array([1, 1, 2, 2], dtype=np.int64)
    pred = np.array([1, 1, 2, 2], dtype=np.int64)
    # no predicted noise and no true noise -> vacuously perfect (min_count=0);
    # with the default minimum the counts are too small -> skipped
    assert m.noise_precision(pred, truth, min_count=0) == 1.0
    assert m.noise_recall(pred, truth, min_count=0) == 1.0
    assert np.isnan(m.noise_precision(pred, truth))
    assert np.isnan(m.noise_recall(pred, truth))


def test_permuted_labels():
    truth = base_truth()
    mapping = {0: 0, 1: 7, 2: 3, 3: 9, 4: 2}
    pred = np.array([mapping[int(t)] for t in truth], dtype=np.int64)
    mask = interior_mask(truth)
    assert m.matched_accuracy(pred, truth, mask) == pytest.approx(1.0)
    assert m.ari(pred, truth, mask) == pytest.approx(1.0)
    assert m.ami(pred, truth, mask) == pytest.approx(1.0)
    assert m.cell_count_ratio(pred, truth, mask) == pytest.approx(1.0)
    assert m.over_segmentation_rate(pred, truth, mask) == pytest.approx(0.0)
    assert m.under_segmentation_rate(pred, truth, mask) == pytest.approx(0.0)
    assert m.recovery_rate(pred, truth, mask) == pytest.approx(1.0)
    assert m.median_matched_jaccard(pred, truth, mask) == pytest.approx(1.0)


def test_interior_mask_excludes_edges():
    truth = np.array([1, 1, 2, 2, 0], dtype=np.int64)
    pred = np.array([1, 2, 2, 2, 1], dtype=np.int64)  # molecule 0 and noise wrong
    mask = np.array([False, True, True, True, False])
    # evaluated on molecules 1..3: pred (2,2,2) vs truth (1,2,2)
    assert m.matched_accuracy(pred, truth, mask) == pytest.approx(2 / 3)


def test_empty_mask_is_nan():
    truth = np.array([1, 1], dtype=np.int64)
    pred = np.array([1, 1], dtype=np.int64)
    mask = np.array([False, False])
    assert np.isnan(m.matched_accuracy(pred, truth, mask))
    assert np.isnan(m.ari(pred, truth, mask))


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        m.matched_accuracy(np.array([1]), np.array([1, 2]))


def test_sim_metrics_bundle():
    truth = base_truth()
    out = m.sim_metrics(truth, truth, None, oracle_accuracy=0.99)
    assert out["accuracy_1to1"] == pytest.approx(1.0)
    assert out["matched_accuracy"] == pytest.approx(1.0)
    assert out["oracle_gap"] == pytest.approx(-0.01)
    assert set(out) == {"accuracy_1to1", "matched_accuracy", "ari",
                        "ari_assigned", "ami", "ami_assigned",
                        "assigned_agreement", "assigned_fraction_pred",
                        "assigned_fraction_truth", "noise_precision",
                        "noise_recall", "cell_count_ratio",
                        "over_segmentation_rate", "under_segmentation_rate",
                        "recovery_rate", "median_matched_jaccard",
                        "oracle_gap"}


def test_ari_assigned_ignores_unassigned_cluster():
    """ari_assigned excludes molecules unassigned on either side; the old
    all-molecule ARI keeps label 0 as one giant cluster."""
    truth = np.array([1] * 40 + [2] * 40 + [0] * 20, dtype=np.int64)
    # prediction agrees on all truly assigned molecules but calls 10 noise
    # molecules cell 1
    pred = truth.copy()
    pred[80:90] = 1
    assert m.ari_assigned(pred, truth) == pytest.approx(1.0)
    assert m.ari(pred, truth) < 1.0
    # status agreement and per-side assigned fractions are reported
    assert m.assigned_agreement(pred, truth) == pytest.approx(0.9)
    assert m.assigned_fraction(pred) == pytest.approx(0.9)
    assert m.assigned_fraction(truth) == pytest.approx(0.8)
    # a prediction that merges everything assigned into one cell: the
    # all-molecule ARI is inflated/deflated by the shared noise cluster,
    # the assigned-only ARI collapses to 0
    merged = truth.copy()
    merged[truth > 0] = 7
    assert m.ari_assigned(merged, truth) == pytest.approx(0.0)


def test_accuracy_1to1_split_penalised():
    truth = base_truth()
    pred = truth.copy()
    pred[:6] = 10
    pred[6:10] = 11
    mask = interior_mask(truth)
    # one-to-one: only one part of the split cell can be matched
    assert m.accuracy_1to1(pred, truth, mask) == pytest.approx(39 / 43)
    # while many-to-one still scores a perfect 1.0
    assert m.matched_accuracy(pred, truth, mask) == pytest.approx(1.0)
    # the oracle gap now uses the one-to-one accuracy
    out = m.sim_metrics(pred, truth, mask, oracle_accuracy=1.0)
    assert out["oracle_gap"] == pytest.approx(1.0 - 39 / 43)


# ---------------------------------------------------------------------------
# difference vs st-recoverability's one-to-one matched_accuracy
# ---------------------------------------------------------------------------

def st_recoverability_matched_accuracy(method_label, true_cell, interior):
    """Reimplementation of st-recoverability src/headroom_common.py.

    Optimal ONE-TO-ONE matching (Hungarian) of method cells to true cells;
    unassigned method labels participate like any other label and can match at
    most one true cell (their docstring: background/unassigned counts as
    errors). Differences from our ``matched_accuracy``: one-to-one vs
    many-to-one majority mapping, and no noise-to-noise special case.
    """
    ml = np.asarray(method_label)[interior]
    tl = np.asarray(true_cell)[interior]
    umeth, mi = np.unique(ml, return_inverse=True)
    utrue, ti = np.unique(tl, return_inverse=True)
    C = np.zeros((len(umeth), len(utrue)), dtype=np.int64)
    np.add.at(C, (mi, ti), 1)
    r, c = linear_sum_assignment(-C)
    return float(C[r, c].sum() / len(ml))


def test_matched_accuracy_vs_st_recoverability_on_split():
    """On a pure split the two definitions diverge by design.

    Ours (many-to-one): both parts map back to their true cell -> 1.0.
    Theirs (one-to-one): only one part can be matched -> penalized.
    """
    truth = base_truth()
    pred = truth.copy()
    pred[:6] = 10
    pred[6:10] = 11
    mask = interior_mask(truth, with_noise=False)  # truth > 0
    ours = m.matched_accuracy(pred, truth, mask)
    theirs = st_recoverability_matched_accuracy(pred, truth, mask)
    assert ours == pytest.approx(1.0)
    assert theirs == pytest.approx(34 / 38)
    assert ours > theirs


def test_matched_accuracy_vs_st_recoverability_on_merge_agree():
    """On a clean one-to-one situation (merge) both definitions agree."""
    truth = np.array([1] * 10 + [2] * 8 + [3] * 10, dtype=np.int64)
    pred = np.array([7] * 18 + [3] * 10, dtype=np.int64)  # cells 1+2 merged
    mask = np.ones(len(truth), dtype=bool)
    ours = m.matched_accuracy(pred, truth, mask)
    theirs = st_recoverability_matched_accuracy(pred, truth, mask)
    assert ours == pytest.approx(theirs)
    # merged cell -> true 1 (majority) scores its 10 molecules; cell 3 keeps 10
    assert ours == pytest.approx(20 / 28)


# ---------------------------------------------------------------------------
# real vs baseline segmentation
# ---------------------------------------------------------------------------

def test_real_perfect_agreement():
    a = np.array([1, 1, 2, 2, 3, 0, 0], dtype=np.int64)
    out = m.real_pair_metrics(a, a)
    assert out["molecule_ari"] == pytest.approx(1.0)
    assert out["ari_assigned"] == pytest.approx(1.0)
    assert out["assigned_agreement"] == pytest.approx(1.0)
    assert out["assigned_fraction_candidate"] == pytest.approx(5 / 7)
    assert out["assigned_fraction_reference"] == pytest.approx(5 / 7)
    assert out["frac_cells_matched"] == pytest.approx(1.0)
    assert out["median_jaccard"] == pytest.approx(1.0)
    assert out["cell_count_ratio"] == pytest.approx(1.0)
    assert out["median_mpc_rel_change"] == pytest.approx(0.0)
    assert out["n_cells_reference"] == 3
    assert out["n_cells_candidate"] == 3


def test_real_partial_disagreement():
    ref = np.array([1, 1, 1, 1, 2, 2, 2, 2, 0, 0], dtype=np.int64)
    new = np.array([1, 1, 1, 5, 2, 2, 2, 0, 0, 5], dtype=np.int64)
    out = m.real_pair_metrics(new, ref)
    # assigned status differs on molecules 3 and 7
    assert out["assigned_agreement"] == pytest.approx(0.8)
    # molecule 9 becomes assigned -> candidate has 8/10, reference 8/10
    assert out["assigned_fraction_candidate"] == pytest.approx(0.8)
    assert out["assigned_fraction_reference"] == pytest.approx(0.8)
    # ARI over the 8 molecules assigned in both: only molecule 3 differs
    assert out["ari_assigned"] == pytest.approx(
        m.ari_assigned(new, ref))
    assert out["ari_assigned"] < 1.0
    # noise precision/recall are skipped (2 predicted-noise molecules)
    assert np.isnan(out["noise_precision"])
    assert np.isnan(out["noise_recall"])
    # reference cell 1: overlap 3, J = 3/(4+3-3) = 0.75 -> matched (>= 0.5)
    # reference cell 2: overlap 3, J = 3/(4+3-3) = 0.75
    assert out["frac_cells_matched"] == pytest.approx(1.0)
    assert out["median_jaccard"] == pytest.approx(0.75)
    assert out["cell_count_ratio"] == pytest.approx(3 / 2)  # cells 1,2,5 vs 1,2
    # median mpc: ref [4,4] -> 4; new counts [3,3,2] -> 3; change = -0.25
    assert out["median_mpc_rel_change"] == pytest.approx(-0.25)
    assert out["molecule_ari"] < 1.0


def test_real_unmatched_reference_cells():
    ref = np.array([1, 1, 2, 2, 3, 3], dtype=np.int64)
    new = np.array([1, 1, 4, 4, 0, 0], dtype=np.int64)
    stats = m.cell_match_stats(ref, new)
    # cell 3 becomes unassigned in `new` -> shares nothing with any new cell
    assert stats["frac_cells_matched"] == pytest.approx(2 / 3)
    assert stats["median_jaccard"] == pytest.approx(1.0)  # median of [1, 1, 0]


def test_real_cell_count_ratio_direction():
    ref = np.array([1, 1, 2, 2, 3, 3], dtype=np.int64)
    new = np.array([1, 1, 2, 2], dtype=np.int64)
    # candidate lost one cell
    assert m.cell_count_ratio_between(new, ref) == pytest.approx(2 / 3)
    assert m.cell_count_ratio_between(ref, new) == pytest.approx(3 / 2)


def test_real_shape_mismatch_raises():
    with pytest.raises(ValueError):
        m.molecule_ari(np.array([1]), np.array([1, 2]))
    with pytest.raises(ValueError):
        m.cell_match_stats(np.array([1]), np.array([1, 2]))


def test_median_mpc_all_unassigned():
    zeros = np.zeros(4, dtype=np.int64)
    assert np.isnan(m.median_mpc_rel_change(zeros, zeros))
    assert np.isnan(m.cell_count_ratio_between(zeros, zeros))
