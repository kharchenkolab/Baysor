"""Tests for the border/dilate segmentation degradation."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from degrade import border_reassign, cell_codes, dilate_cells, load_assignment_table, write_assignment


def two_cell_grid(step: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Two adjacent blocks of cells: cells A (x in [0,10)) and B (x in [10,20))."""
    xs, ys = np.meshgrid(np.arange(0, 20, step), np.arange(0, 5, step))
    xy = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    labels = np.where(xy[:, 0] < 10, "A", "B").astype(object)
    return xy, labels


def test_border_reassign_moves_only_border_molecules():
    xy, labels = two_cell_grid()
    out, stats = border_reassign(xy, labels, fraction=0.1)
    n_expected = int(round(0.1 * len(labels)))
    assert stats["n_assigned"] == len(labels)
    assert stats["n_reassigned"] == n_expected
    # Reassigned molecules must have changed label (targets are foreign cells).
    changed = out != labels
    assert changed.sum() == n_expected
    # The chosen molecules are the ones closest to a foreign molecule: they
    # live next to the A|B seam at x in {9, 10}.
    assert np.all(np.isin(xy[changed][:, 0], [9.0, 10.0]))
    # Molecules deep inside either block are untouched.
    deep = (xy[:, 0] <= 5) | (xy[:, 0] >= 15)
    assert np.all(out[deep] == labels[deep])


def test_border_reassign_fraction_prefix_is_monotone():
    xy, labels = two_cell_grid()
    out10, stats10 = border_reassign(xy, labels, fraction=0.1)
    out30, stats30 = border_reassign(xy, labels, fraction=0.3)
    assert stats10["n_reassigned"] < stats30["n_reassigned"]
    # The 30% variant is a strict superset of the 10% variant: every molecule
    # changed by 10% is also changed by 30%, to the same target cell.
    changed10 = out10 != labels
    changed30 = out30 != labels
    assert np.all(changed10 <= changed30)
    assert np.array_equal(out10[changed10], out30[changed10])


def test_border_reassign_targets_nearest_foreign_cell():
    # One molecule of A sits right next to a B molecule; it must go to B.
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [50.0, 50.0], [51.0, 50.0]])
    labels = np.array(["A", "A", "B", "B"], dtype=object)
    out, stats = border_reassign(xy, labels, fraction=0.5)
    assert stats["n_reassigned"] == 2
    # The two seam-adjacent molecules swap to the other cell.
    assert out[0] == "B" or out[1] == "B"


def test_border_reassign_background_untouched_and_invalid_fraction():
    xy, labels = two_cell_grid()
    labels = labels.copy()
    labels[0] = ""  # background molecule
    out, stats = border_reassign(xy, labels, fraction=0.1)
    assert out[0] == ""
    assert stats["n_unassigned"] == 1
    with pytest.raises(ValueError):
        border_reassign(xy, labels, fraction=0.0)
    with pytest.raises(ValueError):
        border_reassign(xy, labels, fraction=1.5)


def test_dilate_absorbs_background_within_distance():
    # A single cell (5,5) plus a background molecule 1um away and one far away.
    xy = np.array([[5.0, 5.0], [6.0, 5.0], [5.0, 6.0], [4.0, 4.0],
                   [30.0, 30.0], [31.0, 30.0]])
    labels = np.array(["A", "A", "A", "", "", ""], dtype=object)
    out, stats = dilate_cells(xy, labels, distance=2.0)
    # Background molecules within 2um of the hull are absorbed by A.
    assert out[3] == "A"
    # Background beyond the dilation distance stays background.
    assert out[4] == "" and out[5] == ""
    assert stats["n_absorbed_background"] == 1
    assert out[0] == "A" and out[1] == "A" and out[2] == "A"


def test_dilate_takes_border_molecules_from_neighbour():
    # Cells A and B share a border at x=10; A's molecules extend to x=9.5
    # (i.e. 0.5um past the seam into A... build explicit hulls instead):
    # A occupies x in [0,10], B occupies x in [10.5, 20]. A molecule of A at
    # x=9.9 is 0.6um from B's hull -> claimed by B under 2um dilation.
    ya = np.array([[0.0, 0.0], [0.0, 5.0], [10.0, 0.0], [10.0, 5.0], [9.9, 2.5]])
    yb = np.array([[10.5, 0.0], [10.5, 5.0], [20.0, 0.0], [20.0, 5.0]])
    xy = np.vstack([ya, yb])
    labels = np.array(["A"] * len(ya) + ["B"] * len(yb), dtype=object)
    out, stats = dilate_cells(xy, labels, distance=2.0)
    # The A molecule nearest the seam flips to B ...
    assert out[4] == "B"
    # ... and B's first band (0.5um inside B) is in turn claimed by A:
    # dilation churns the shared border symmetrically.
    assert out[len(ya)] == "A"
    # Deep interior molecules keep their label.
    assert out[0] == "A"
    assert out[-1] == "B"
    assert stats["n_reassigned_from_neighbours"] >= 1


def test_dilate_deterministic():
    rng = np.random.default_rng(7)
    xy = rng.uniform(0, 100, size=(500, 2))
    labels = np.where(xy[:, 0] < 50, "A", "B").astype(object)
    out1, _ = dilate_cells(xy, labels, distance=2.0)
    out2, _ = dilate_cells(xy, labels, distance=2.0)
    assert np.array_equal(out1, out2)


def test_cell_codes_roundtrip(tmp_path: Path):
    labels = np.array(["b", "", "a", "b", ""], dtype=object)
    codes, mapping = cell_codes(labels)
    assert mapping == {"a": 1, "b": 2}
    assert codes.tolist() == [2, 0, 1, 2, 0]
    out = tmp_path / "assignment.parquet"
    write_assignment(out, labels)
    frame = pd.read_parquet(out)
    assert frame["cell"].tolist() == codes.tolist()
    assert frame["cell_label"].tolist() == labels.astype(str).tolist()


def test_load_assignment_table(tmp_path: Path):
    src = tmp_path / "molecules.parquet"
    pd.DataFrame({
        "x": [1.0, 2.0, 3.0], "y": [3.0, 4.0, 5.0],
        "gene": ["g1", "g2", "g3"],
        "cell_vendor": ["c1", "", "UNASSIGNED"],
    }).to_parquet(src, index=False)
    frame = load_assignment_table(src, "cell_vendor")
    assert frame.columns.tolist() == ["x", "y", "cell_vendor"]
    assert len(frame) == 3
    # Loader-specific unassigned tokens become background, not a cell: the
    # Xenium "UNASSIGNED" pseudo-cell would otherwise dominate the hulls.
    assert frame["cell_vendor"].tolist() == ["c1", "", ""]


def test_unassigned_pseudo_cell_never_degraded(tmp_path: Path):
    src = tmp_path / "molecules.parquet"
    rng = np.random.default_rng(1)
    n = 200
    # background molecules spread over the whole tissue + two real cells
    xy_bg = rng.uniform(0, 50, size=(n, 2))
    xy_a = rng.uniform(0, 5, size=(n, 2))
    xy_b = rng.uniform(20, 25, size=(n, 2))
    xy = np.vstack([xy_a, xy_b, xy_bg])
    pd.DataFrame({
        "x": xy[:, 0], "y": xy[:, 1], "gene": "g",
        "cell_vendor": (["A"] * n + ["B"] * n + ["UNASSIGNED"] * n),
    }).to_parquet(src, index=False)
    frame = load_assignment_table(src, "cell_vendor")
    labels = frame["cell_vendor"].to_numpy()
    assert (labels == "").sum() == n
    # Background never forms a hull and is never a reassignment target.
    out, stats = dilate_cells(xy, labels, distance=2.0)
    assert not (out == "UNASSIGNED").any()
    assert stats["n_cells"] == 2
    out2, stats2 = border_reassign(xy, labels, fraction=0.1)
    # Only A/B molecules are candidates; targets are real cells.
    assert set(np.unique(out2[out2 != ""])) <= {"A", "B"}
    assert stats2["n_assigned"] == 2 * n
