"""Tests for cell-type transfer by molecule overlap."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from transfer import read_celltypes, read_target_cells, transfer_celltypes, write_celltypes


def test_majority_vote_and_tie_break():
    target = np.array(["t1", "t1", "t1", "t2", "t2", "t3", ""], dtype=object)
    baseline = np.array(["b1", "b1", "b2", "b2", "b3", "", ""], dtype=object)
    types = pd.Series({"b1": "typeA", "b2": "typeB", "b3": "typeC"}, name="celltype")
    out, stats = transfer_celltypes(target, baseline, types)
    assert out["t1"] == "typeA"          # 2 vs 1 majority
    # 1-1 tie -> lexicographically smallest type
    assert out["t2"] == "typeB"
    # t3 has no typed baseline molecules -> dropped from the annotation
    assert "t3" not in out.index
    assert stats["n_dropped_untyped_target_cells"] == 1
    assert stats["n_typed_target_cells"] == 2


def test_min_votes_filter():
    target = np.array(["t1", "t2"], dtype=object)
    baseline = np.array(["b1", ""], dtype=object)
    types = pd.Series({"b1": "typeA"}, name="celltype")
    out, stats = transfer_celltypes(target, baseline, types, min_votes=1)
    assert list(out.index) == ["t1"]
    out2, _ = transfer_celltypes(target, baseline, types, min_votes=2)
    assert len(out2) == 0


def test_unknown_baseline_cell_types_are_untyped():
    target = np.array(["t1", "t1"], dtype=object)
    baseline = np.array(["unknown_cell", "also_unknown"], dtype=object)
    types = pd.Series({"known_cell": "typeA"}, name="celltype")
    out, stats = transfer_celltypes(target, baseline, types)
    assert len(out) == 0
    assert stats["n_molecules_with_baseline_type"] == 0


def test_row_alignment_enforced():
    with pytest.raises(ValueError):
        transfer_celltypes(np.array(["t1"]), np.array(["A", "B"]), pd.Series({"a": "A"}))


def test_read_write_celltypes_roundtrip(tmp_path: Path):
    types = pd.Series({"c1": "ductal", "c2": "immune"}, name="celltype")
    path = tmp_path / "ct.parquet"
    write_celltypes(path, types)
    back = read_celltypes(path)
    assert back["c1"] == "ductal"
    assert back["c2"] == "immune"
    # CSV with nonstandard column names works via explicit columns.
    csv = tmp_path / "ct.csv"
    pd.DataFrame({"cid": ["c1"], "label": ["ductal"]}).to_csv(csv, index=False)
    back2 = read_celltypes(csv, cell_col="cid", type_col="label")
    assert back2["c1"] == "ductal"


def test_read_target_cells_cell_column_and_assignment(tmp_path: Path):
    molecules = tmp_path / "molecules.parquet"
    pd.DataFrame({
        "x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0],
        "gene": ["g", "g", "g"], "cell_vendor": ["v1", "", "v2"],
    }).to_parquet(molecules, index=False)
    labels = read_target_cells(molecules, assignment=None, cell_column="cell_vendor")
    assert labels.tolist() == ["v1", "", "v2"]

    # int assignment without labels -> stringified codes, 0 becomes ''
    assignment = tmp_path / "a.parquet"
    pd.DataFrame({"cell": [3, 0, 3]}).to_parquet(assignment, index=False)
    labels2 = read_target_cells(molecules, assignment=assignment, cell_column=None)
    assert labels2.tolist() == ["3", "", "3"]

    # cell_label sidecar wins over codes
    assignment2 = tmp_path / "a2.parquet"
    pd.DataFrame({"cell": [3, 0, 3], "cell_label": ["v1", "", "v2"]}).to_parquet(
        assignment2, index=False)
    labels3 = read_target_cells(molecules, assignment=assignment2, cell_column=None)
    assert labels3.tolist() == ["v1", "", "v2"]

    # row-order mismatch raises
    bad = tmp_path / "bad.parquet"
    pd.DataFrame({"cell": [1, 2]}).to_parquet(bad, index=False)
    with pytest.raises(ValueError):
        read_target_cells(molecules, assignment=bad, cell_column=None)

    # both or neither
    with pytest.raises(ValueError):
        read_target_cells(molecules, assignment=None, cell_column=None)
    with pytest.raises(ValueError):
        read_target_cells(molecules, assignment=assignment, cell_column="cell_vendor")
