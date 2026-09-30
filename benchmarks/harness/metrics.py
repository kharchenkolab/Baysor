"""Pure metric functions for the Baysor benchmark harness.

All functions take int64 label vectors (0 = unassigned / noise) as numpy
arrays and return plain floats. Label 0 always counts as its own label in
cluster metrics. Nothing here reads or writes files.

Primary agreement metrics ignore the giant unassigned cluster: ``ari_assigned``
/ ``ami_assigned`` are computed only over molecules assigned in *both*
sides, and the assigned status of each side is reported separately
(``assigned_agreement``, ``assigned_fraction_*``, ``noise_precision`` /
``noise_recall``). The primary sim accuracy is the one-to-one (Hungarian)
``accuracy_1to1``; the many-to-one ``matched_accuracy`` is kept as a
secondary metric.

Sim-vs-truth metrics operate on molecules with ``interior == True`` (the mask
defaults to all molecules). Real-vs-baseline metrics compare two assignments
of the *same* input molecules, so both vectors must be aligned by
``mol_index`` beforehand.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import (adjusted_rand_score,
                             adjusted_mutual_info_score)

__all__ = [
    "accuracy_1to1", "matched_accuracy", "ari", "ami", "ari_assigned",
    "ami_assigned", "assigned_agreement", "assigned_fraction",
    "noise_precision", "noise_recall",
    "cell_count_ratio", "over_segmentation_rate", "under_segmentation_rate",
    "recovery_rate", "median_matched_jaccard", "oracle_gap", "sim_metrics",
    "molecule_ari", "cell_match_stats",
    "cell_count_ratio_between", "median_mpc_rel_change", "real_pair_metrics",
    "NOISE_MIN_COUNT", "JACCARD_THRESHOLD",
]

JACCARD_THRESHOLD = 0.5
OVERSPLIT_FRACTION = 0.8   # largest true->pred part must hold >= 80% of the cell
NOISE_MIN_COUNT = 50       # below this, noise precision/recall are undefined


# ---------------------------------------------------------------------------
# small internal helpers
# ---------------------------------------------------------------------------

def _apply_mask(pred, truth, mask):
    pred = np.asarray(pred)
    truth = np.asarray(truth)
    if pred.shape != truth.shape:
        raise ValueError(f"shape mismatch: pred {pred.shape} vs truth {truth.shape}")
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        pred, truth = pred[mask], truth[mask]
    return pred.astype(np.int64), truth.astype(np.int64)


def _contingency(pred: np.ndarray, truth: np.ndarray):
    """Contingency between two label vectors.

    Returns (p_labels, t_labels, row_idx, col_idx, counts) where row_idx/col_idx
    index into p_labels/t_labels respectively. Unique labels come out sorted.
    """
    if len(pred) == 0:
        z = np.zeros(0, dtype=np.int64)
        return z, z.copy(), z.copy(), z.copy(), z.copy()
    p_labels, pi = np.unique(pred, return_inverse=True)
    t_labels, ti = np.unique(truth, return_inverse=True)
    key = pi.astype(np.int64) * len(t_labels) + ti
    uk, counts = np.unique(key, return_counts=True)
    return p_labels, t_labels, uk // len(t_labels), uk % len(t_labels), counts


def _nan() -> float:
    return float("nan")


# ---------------------------------------------------------------------------
# sim vs truth
# ---------------------------------------------------------------------------

def matched_accuracy(pred, truth, mask=None) -> float:
    """Fraction of molecules whose predicted cell maps to their own true cell.

    Each predicted cell (including 0 = noise) is mapped to the true cell it
    shares most molecules with; ties break towards the smallest true label.
    Special rules for noise:

    * predicted noise is only ever correct against true noise (noise maps to
      noise, never to a real cell);
    * a real predicted cell whose majority is true noise earns **no credit**
      for its noise molecules (its winner is the noise label, which never
      counts for a real cell).

    Many-to-one: several predicted cells may map to the same true cell. The
    one-to-one variant is :func:`accuracy_1to1`.
    """
    pred, truth = _apply_mask(pred, truth, mask)
    n = len(pred)
    if n == 0:
        return _nan()
    p_labels, t_labels, ri, ti, counts = _contingency(pred, truth)
    if len(p_labels) == 0:
        return _nan()
    # winner per predicted label: highest count, tie -> smallest true label
    order = np.lexsort((t_labels[ti], -counts, ri))
    rows = ri[order]
    first = order[np.concatenate(([True], rows[1:] != rows[:-1]))]
    winner = np.full(len(p_labels), -1, dtype=np.int64)
    winner[ri[first]] = ti[first]
    correct = 0
    for r, c, cnt in zip(ri, ti, counts):
        if winner[r] < 0:
            continue
        if p_labels[r] == 0:
            # noise is correct only against true noise
            if t_labels[c] == 0:
                correct += int(cnt)
        elif t_labels[winner[r]] == 0:
            # majority is true noise -> this real cell earns no credit for
            # noise molecules (c == winner[r] is the only possible hit here)
            continue
        elif c == winner[r]:
            correct += int(cnt)
    return correct / n


def accuracy_1to1(pred, truth, mask=None) -> float:
    """One-to-one (Hungarian) assignment accuracy — the PRIMARY sim metric.

    Builds the predicted x true overlap matrix over all labels (0 included)
    and matches predicted cells to true cells with
    ``scipy.optimize.linear_sum_assignment`` to maximise credited molecules.
    A pair earns credit only when the two sides agree on noise status:
    predicted noise matches true noise (and nothing else), and a real
    predicted cell never earns credit from true-noise molecules. Unmatched
    labels and mismatched-noise pairs contribute 0. A pure split therefore
    scores below 1.0 (only one part can be matched).
    """
    pred, truth = _apply_mask(pred, truth, mask)
    n = len(pred)
    if n == 0:
        return _nan()
    p_labels, t_labels, ri, ti, counts = _contingency(pred, truth)
    overlap = np.zeros((len(p_labels), len(t_labels)), dtype=np.float64)
    np.add.at(overlap, (ri, ti), counts)
    noise_status_agrees = ((p_labels > 0)[:, None] == (t_labels > 0)[None, :])
    credit = np.where(noise_status_agrees, overlap, 0.0)
    r, c = linear_sum_assignment(-credit)
    return float(credit[r, c].sum() / n)


def ari(pred, truth, mask=None) -> float:
    """Adjusted Rand index; label 0 is its own label."""
    pred, truth = _apply_mask(pred, truth, mask)
    if len(pred) == 0:
        return _nan()
    return float(adjusted_rand_score(truth, pred))


AMI_MAX_LABELS = 3000


def ami(pred, truth, mask=None) -> float:
    """Adjusted mutual information; label 0 is its own label.

    NaN when either side has more than ``AMI_MAX_LABELS`` distinct labels:
    sklearn's expected-mutual-information term scales quadratically with the
    number of clusters and takes hours on full-tier crops. AMI is
    informational only; no gate uses it.
    """
    pred, truth = _apply_mask(pred, truth, mask)
    if len(pred) == 0:
        return _nan()
    if max(len(np.unique(pred)), len(np.unique(truth))) > AMI_MAX_LABELS:
        return _nan()
    return float(adjusted_mutual_info_score(truth, pred))


def ari_assigned(pred, truth, mask=None) -> float:
    """ARI over molecules assigned in BOTH sides (primary agreement metric).

    Unassigned molecules are excluded instead of forming one giant cluster;
    NaN when no molecule is assigned on both sides.
    """
    pred, truth = _apply_mask(pred, truth, mask)
    both = (pred > 0) & (truth > 0)
    if not both.any():
        return _nan()
    return float(adjusted_rand_score(truth[both], pred[both]))


def ami_assigned(pred, truth, mask=None) -> float:
    """AMI over molecules assigned in BOTH sides (NaN when none)."""
    pred, truth = _apply_mask(pred, truth, mask)
    both = (pred > 0) & (truth > 0)
    if not both.any():
        return _nan()
    return ami(pred[both], truth[both])


def assigned_fraction(a, mask=None) -> float:
    """Fraction of molecules with a label > 0 in one label vector."""
    a = np.asarray(a)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if a.shape != mask.shape:
            raise ValueError(f"shape mismatch: {a.shape} vs {mask.shape}")
        a = a[mask]
    if len(a) == 0:
        return _nan()
    return float((a > 0).mean())


def noise_precision(pred, truth, mask=None,
                    min_count: int = NOISE_MIN_COUNT) -> float:
    """P(true noise | predicted noise).

    NaN (skipped) when fewer than ``min_count`` molecules are predicted noise;
    vacuously 1.0 when ``min_count == 0`` and nothing is predicted noise.
    """
    pred, truth = _apply_mask(pred, truth, mask)
    pn = pred == 0
    n_pred = int(pn.sum())
    if n_pred < min_count:
        return _nan()
    if n_pred == 0:
        return 1.0
    return float((pn & (truth == 0)).sum() / n_pred)


def noise_recall(pred, truth, mask=None,
                 min_count: int = NOISE_MIN_COUNT) -> float:
    """P(predicted noise | true noise).

    NaN (skipped) when fewer than ``min_count`` molecules are true noise;
    vacuously 1.0 when ``min_count == 0`` and there is no true noise.
    """
    pred, truth = _apply_mask(pred, truth, mask)
    tn = truth == 0
    n_true = int(tn.sum())
    if n_true < min_count:
        return _nan()
    if n_true == 0:
        return 1.0
    return float((tn & (pred == 0)).sum() / n_true)


def cell_count_ratio(pred, truth, mask=None) -> float:
    """#predicted cells / #true cells (labels > 0)."""
    pred, truth = _apply_mask(pred, truth, mask)
    n_pred = int(len(np.unique(pred[pred > 0])))
    n_true = int(len(np.unique(truth[truth > 0])))
    if n_true == 0:
        return _nan()
    return n_pred / n_true


def over_segmentation_rate(pred, truth, mask=None) -> float:
    """Fraction of true cells split so that their largest predicted part holds < 80%."""
    pred, truth = _apply_mask(pred, truth, mask)
    _, t_labels, ri, ti, counts = _contingency(pred, truth)
    if len(t_labels) == 0:
        return _nan()
    totals = np.zeros(len(t_labels), dtype=np.int64)
    np.add.at(totals, ti, counts)
    best = np.zeros(len(t_labels), dtype=np.int64)
    np.maximum.at(best, ti, counts)
    frac = np.divide(best, totals, out=np.zeros(len(t_labels), dtype=float), where=totals > 0)
    real_cells = t_labels > 0
    if not real_cells.any():
        return _nan()
    return float((frac[real_cells] < OVERSPLIT_FRACTION).mean())


def under_segmentation_rate(pred, truth, mask=None) -> float:
    """Fraction of predicted cells whose largest true source holds < 80%."""
    pred, truth = _apply_mask(pred, truth, mask)
    p_labels, _, ri, ti, counts = _contingency(pred, truth)
    if len(p_labels) == 0:
        return _nan()
    totals = np.zeros(len(p_labels), dtype=np.int64)
    np.add.at(totals, ri, counts)
    best = np.zeros(len(p_labels), dtype=np.int64)
    np.maximum.at(best, ri, counts)
    frac = np.divide(best, totals, out=np.zeros(len(p_labels), dtype=float), where=totals > 0)
    real_cells = p_labels > 0
    if not real_cells.any():
        return _nan()
    return float((frac[real_cells] < OVERSPLIT_FRACTION).mean())


def _best_jaccards(pred, truth) -> tuple[np.ndarray, np.ndarray]:
    """For each true label > 0, the best Jaccard against predicted labels > 0.

    Ties break towards the smallest predicted label. True cells without any
    overlapping predicted cell get 0.
    """
    p_labels, t_labels, ri, ti, counts = _contingency(pred, truth)
    t_mask = t_labels > 0
    p_sizes = np.zeros(len(p_labels), dtype=np.int64)
    t_sizes = np.zeros(len(t_labels), dtype=np.int64)
    np.add.at(p_sizes, ri, counts)
    np.add.at(t_sizes, ti, counts)
    best = np.zeros(len(t_labels), dtype=np.int64)
    argmax_p = np.full(len(t_labels), -1, dtype=np.int64)
    sel = (p_labels[ri] > 0) & t_mask[ti]
    if sel.any():
        sri, sti, sc = ri[sel], ti[sel], counts[sel]
        # primary: true label asc, then count desc, then pred label asc
        order = np.lexsort((sri, -sc, sti))
        ts = sti[order]
        first = order[np.concatenate(([True], ts[1:] != ts[:-1]))]
        wt = sti[first]
        best[wt] = sc[first]
        argmax_p[wt] = sri[first]
    inter = best.astype(float)
    union = t_sizes.astype(float)
    has = argmax_p >= 0
    union[has] += p_sizes[argmax_p[has]]
    union[has] -= inter[has]
    jac = np.divide(inter, union, out=np.zeros(len(t_labels), dtype=float),
                    where=union > 0)
    return t_labels, jac


def recovery_rate(pred, truth, mask=None,
                  threshold: float = JACCARD_THRESHOLD) -> float:
    """Fraction of true cells recovered with molecule-set Jaccard >= threshold."""
    pred, truth = _apply_mask(pred, truth, mask)
    if len(pred) == 0:
        return _nan()
    t_labels, jac = _best_jaccards(pred, truth)
    real_cells = t_labels > 0
    if not real_cells.any():
        return _nan()
    return float((jac[real_cells] >= threshold).mean())


def median_matched_jaccard(pred, truth, mask=None) -> float:
    """Median, over true cells, of the best Jaccard against a predicted cell."""
    pred, truth = _apply_mask(pred, truth, mask)
    if len(pred) == 0:
        return _nan()
    t_labels, jac = _best_jaccards(pred, truth)
    real_cells = t_labels > 0
    if not real_cells.any():
        return _nan()
    return float(np.median(jac[real_cells]))


def oracle_gap(accuracy_1to1: float, oracle_accuracy) -> float:
    """Gap between the one-to-one accuracy and the generator's oracle
    accuracy (the oracle is identity/one-to-one), NaN when not provided.
    """
    if oracle_accuracy is None:
        return _nan()
    try:
        oa = float(oracle_accuracy)
    except (TypeError, ValueError):
        return _nan()
    if np.isnan(oa) or np.isnan(accuracy_1to1):
        return _nan()
    return oa - accuracy_1to1


def sim_metrics(pred, truth, interior=None, oracle_accuracy=None) -> dict:
    """All sim-vs-truth metrics in one dict (interior mask optional).

    Primary: ``accuracy_1to1`` (Hungarian), ``ari_assigned``,
    ``recovery_rate``, ``cell_count_ratio``. Everything else is secondary /
    informational. ``oracle_gap`` is computed from ``accuracy_1to1``.
    """
    mask = None if interior is None else np.asarray(interior, dtype=bool)
    acc1 = accuracy_1to1(pred, truth, mask)
    return {
        "accuracy_1to1": acc1,
        "matched_accuracy": matched_accuracy(pred, truth, mask),
        "ari": ari(pred, truth, mask),
        "ari_assigned": ari_assigned(pred, truth, mask),
        "ami": ami(pred, truth, mask),
        "ami_assigned": ami_assigned(pred, truth, mask),
        "assigned_agreement": assigned_agreement(pred, truth, mask),
        "assigned_fraction_pred": assigned_fraction(pred, mask),
        "assigned_fraction_truth": assigned_fraction(truth, mask),
        "noise_precision": noise_precision(pred, truth, mask),
        "noise_recall": noise_recall(pred, truth, mask),
        "cell_count_ratio": cell_count_ratio(pred, truth, mask),
        "over_segmentation_rate": over_segmentation_rate(pred, truth, mask),
        "under_segmentation_rate": under_segmentation_rate(pred, truth, mask),
        "recovery_rate": recovery_rate(pred, truth, mask),
        "median_matched_jaccard": median_matched_jaccard(pred, truth, mask),
        "oracle_gap": oracle_gap(acc1, oracle_accuracy),
    }


# ---------------------------------------------------------------------------
# real vs baseline segmentation
# ---------------------------------------------------------------------------

def molecule_ari(a, b) -> float:
    """Molecule-level ARI between two assignments (label 0 = unassigned)."""
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    if len(a) == 0:
        return _nan()
    return float(adjusted_rand_score(a, b))


def assigned_agreement(a, b, mask=None) -> float:
    """Fraction of molecules that agree on assigned-vs-unassigned status."""
    if mask is not None:
        a, b = _apply_mask(a, b, mask)
    else:
        a = np.asarray(a, dtype=np.int64)
        b = np.asarray(b, dtype=np.int64)
        if a.shape != b.shape:
            raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    if len(a) == 0:
        return _nan()
    return float(((a > 0) == (b > 0)).mean())


def cell_match_stats(reference, other,
                     threshold: float = JACCARD_THRESHOLD) -> dict:
    """Cell matching by molecule-set Jaccard, reference -> other.

    Returns ``frac_cells_matched`` (fraction of reference cells with best Jaccard
    >= threshold) and ``median_jaccard`` (median of the best Jaccards).
    """
    reference = np.asarray(reference, dtype=np.int64)
    other = np.asarray(other, dtype=np.int64)
    if reference.shape != other.shape:
        raise ValueError(f"shape mismatch: {reference.shape} vs {other.shape}")
    r_lab, o_lab, ri, oi, counts = _contingency(reference, other)
    n_ref_cells = int((r_lab > 0).sum())
    if n_ref_cells == 0:
        return {"frac_cells_matched": _nan(), "median_jaccard": _nan(),
                "n_reference_cells": 0, "n_other_cells": int((o_lab > 0).sum())}
    r_sizes = np.zeros(len(r_lab), dtype=np.int64)
    o_sizes = np.zeros(len(o_lab), dtype=np.int64)
    np.add.at(r_sizes, ri, counts)
    np.add.at(o_sizes, oi, counts)
    sel = (r_lab[ri] > 0) & (o_lab[oi] > 0)
    best = np.zeros(len(r_lab), dtype=np.int64)
    arg = np.full(len(r_lab), -1, dtype=np.int64)
    if sel.any():
        sri, soi, sc = ri[sel], oi[sel], counts[sel]
        # primary: reference label asc, then count desc, then other label asc
        order = np.lexsort((soi, -sc, sri))
        rs = sri[order]
        first = order[np.concatenate(([True], rs[1:] != rs[:-1]))]
        wr = sri[first]
        best[wr] = sc[first]
        arg[wr] = soi[first]
    jac = np.zeros(len(r_lab), dtype=float)
    has = (r_lab > 0) & (arg >= 0)
    inter = best[has].astype(float)
    union = r_sizes[has].astype(float) + o_sizes[arg[has]].astype(float) - inter
    jac[has] = np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)
    ref_jac = jac[r_lab > 0]
    return {
        "frac_cells_matched": float((ref_jac >= threshold).mean()),
        "median_jaccard": float(np.median(ref_jac)),
        "n_reference_cells": n_ref_cells,
        "n_other_cells": int((o_lab > 0).sum()),
    }


def cell_count_ratio_between(new, reference) -> float:
    """#cells(new) / #cells(reference), labels > 0."""
    new = np.asarray(new, dtype=np.int64)
    reference = np.asarray(reference, dtype=np.int64)
    n_new = int(len(np.unique(new[new > 0])))
    n_ref = int(len(np.unique(reference[reference > 0])))
    if n_ref == 0:
        return _nan()
    return n_new / n_ref


def _median_mpc(x: np.ndarray) -> float:
    assigned = x[x > 0]
    if len(assigned) == 0:
        return float("nan")
    _, counts = np.unique(assigned, return_counts=True)
    return float(np.median(counts))


def median_mpc_rel_change(new, reference) -> float:
    """Relative change of the median molecules-per-cell: (new - ref) / ref."""
    new_med = _median_mpc(np.asarray(new, dtype=np.int64))
    ref_med = _median_mpc(np.asarray(reference, dtype=np.int64))
    if np.isnan(new_med) or np.isnan(ref_med) or ref_med == 0:
        return _nan()
    return (new_med - ref_med) / ref_med


def real_pair_metrics(candidate, reference) -> dict:
    """All real-segmentation comparison metrics for one (candidate, reference) pair.

    ``reference`` is the baseline (or lower-index replicate); ``candidate`` is
    the new assignment. Primary: ``ari_assigned`` (ARI over molecules
    assigned in both), ``frac_cells_matched``, ``cell_count_ratio``. The
    assigned status of each side and the noise precision/recall (reference
    plays the role of truth) are reported separately; ``molecule_ari`` keeps
    the all-molecule value for reference only.
    """
    candidate = np.asarray(candidate, dtype=np.int64)
    reference = np.asarray(reference, dtype=np.int64)
    stats = cell_match_stats(reference, candidate)
    return {
        "molecule_ari": molecule_ari(candidate, reference),
        "ari_assigned": ari_assigned(candidate, reference),
        "assigned_agreement": assigned_agreement(candidate, reference),
        "assigned_fraction_candidate": assigned_fraction(candidate),
        "assigned_fraction_reference": assigned_fraction(reference),
        "noise_precision": noise_precision(candidate, reference),
        "noise_recall": noise_recall(candidate, reference),
        "frac_cells_matched": stats["frac_cells_matched"],
        "median_jaccard": stats["median_jaccard"],
        "cell_count_ratio": cell_count_ratio_between(candidate, reference),
        "median_mpc_rel_change": median_mpc_rel_change(candidate, reference),
        "n_cells_reference": float(stats["n_reference_cells"]),
        "n_cells_candidate": float(stats["n_other_cells"]),
    }
