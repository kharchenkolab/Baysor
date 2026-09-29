"""Tests for audit.py helper logic (metrics, top pairs, assignment)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit import (AUDIT_PARAM_DEFAULTS, UNASSIGNED_TOKENS, compute_metrics,
                   filter_rows, resolve_assignment, seed_id, sha256_file,
                   top_pairs, write_audit_molecules)


def make_pairs(rows):
    return pd.DataFrame(rows, columns=[
        "source", "target", "rate", "admixed_molecules", "excess", "coverage",
        "q_value", "detected", "n_exposed", "n_reference", "n_markers", "n_strict"])


def test_compute_metrics_total_is_sum_over_detected():
    pairs = make_pairs([
        ("A", "B", 0.1, 100.0, 50.0, 0.5, 0.001, True, 10, 5, 20, 8),
        ("B", "A", 0.2, 300.0, 60.0, 0.2, 0.001, True, 12, 6, 20, 4),
        ("A", "C", 0.0, 0.0, 0.0, 0.3, 0.5, False, 3, 9, 20, 2),
    ])
    m = compute_metrics(pairs, total_molecules=10_000.0)
    assert m["total_admixture_molecules"] == 400.0
    assert m["total_admixture_rate"] == pytest.approx(0.04)
    assert m["n_pairs_evaluated"] == 3
    assert m["n_pairs_detected"] == 2
    assert m["status"] == "ok"


def test_compute_metrics_degenerate_cases():
    empty = compute_metrics(make_pairs([]), 100.0)
    assert empty["total_admixture_rate"] == 0.0
    assert empty["status"] == "no_pairs_evaluated"
    none_detected = compute_metrics(
        make_pairs([("A", "B", 0.1, 10.0, 5.0, 0.5, 0.4, False, 1, 1, 5, 2)]), 100.0)
    assert none_detected["total_admixture_rate"] == 0.0
    assert none_detected["status"] == "no_detected_pairs"
    # Denominator is guarded against zero.
    zero = compute_metrics(
        make_pairs([("A", "B", 0.1, 10.0, 5.0, 0.5, 0.001, True, 1, 1, 5, 2)]), 0.0)
    assert np.isfinite(zero["total_admixture_rate"])


def test_top_pairs_orders_by_admixed_molecules_and_limits():
    pairs = make_pairs([
        ("A", "B", 0.1, 100.0, 50.0, 0.5, 0.001, True, 10, 5, 20, 8),
        ("B", "A", 0.2, 500.0, 60.0, 0.2, 0.001, True, 12, 6, 20, 4),
        ("A", "C", 0.0, 0.0, 0.0, 0.3, 0.5, False, 3, 9, 20, 2),
        ("C", "A", 0.05, 200.0, 10.0, 0.1, 0.001, True, 4, 4, 20, 3),
    ])
    top = top_pairs(pairs, 2)
    assert len(top) == 2
    assert top[0]["source"] == "B" and top[0]["target"] == "A"
    assert top[1]["source"] == "C"
    assert set(top[0]) >= {"source", "target", "rate", "admixed_molecules", "q_value"}
    assert top_pairs(make_pairs([]), 5) == []


def _args(molecules=None, cell_column=None, assignment=None, assignment_col="cell"):
    return argparse.Namespace(
        molecules=molecules, cell_column=cell_column,
        assignment=assignment, assignment_col=assignment_col)


def test_resolve_assignment_cell_column(tmp_path: Path):
    src = tmp_path / "m.parquet"
    pd.DataFrame({
        "x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0],
        "gene": ["g", "g", "g"],
        "cell_vendor": ["c1", "UNASSIGNED", "cell_0"],
    }).to_parquet(src, index=False)
    mol = pd.read_parquet(src)
    labels, info = resolve_assignment(_args(cell_column="cell_vendor"), mol)
    # All loader-specific unassigned tokens normalize to ''.
    assert labels.tolist() == ["c1", "", ""]
    assert info["assignment"]["kind"] == "cell_column"


def test_resolve_assignment_file_int_and_labels(tmp_path: Path):
    src = tmp_path / "m.parquet"
    pd.DataFrame({"x": [1.0, 2.0], "y": [1.0, 2.0], "gene": ["g", "g"]}).to_parquet(
        src, index=False)
    mol = pd.read_parquet(src)

    a = tmp_path / "a.parquet"
    pd.DataFrame({"cell": [7, 0]}).to_parquet(a, index=False)
    labels, info = resolve_assignment(_args(assignment=a), mol)
    assert labels.tolist() == ["7", ""]
    assert info["assignment"]["label_column"] is None

    b = tmp_path / "b.parquet"
    pd.DataFrame({"cell": [7, 0], "cell_label": ["v1", ""]}).to_parquet(b, index=False)
    labels2, info2 = resolve_assignment(_args(assignment=b), mol)
    assert labels2.tolist() == ["v1", ""]
    assert info2["assignment"]["label_column"] == "cell_label"

    # Row-count mismatch aborts.
    c = tmp_path / "c.parquet"
    pd.DataFrame({"cell": [1, 2, 3]}).to_parquet(c, index=False)
    with pytest.raises(SystemExit):
        resolve_assignment(_args(assignment=c), mol)


def test_filter_rows_drops_malformed():
    mol = pd.DataFrame({
        "x": [1.0, np.nan, 3.0],
        "y": [1.0, 2.0, np.inf],
        "gene": ["g", "g", ""],
    })
    labels = np.array(["a", "b", "c"], dtype=object)
    out_mol, out_labels, filters = filter_rows(mol, labels)
    assert len(out_mol) == 1
    assert out_labels.tolist() == ["a"]
    assert filters["malformed_rows_dropped"] == 2


def test_write_audit_molecules_excludes_unassigned(tmp_path: Path):
    mol = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0],
                        "gene": ["g1", "g2", "g3"]})
    labels = np.array(["c1", "", "c2"], dtype=object)
    out = tmp_path / "t.parquet"
    n = write_audit_molecules(mol, labels, out)
    assert n == 2
    table = pd.read_parquet(out)
    assert table["cell"].tolist() == ["c1", "c2"]
    assert table["gene"].tolist() == ["g1", "g3"]
    # Row order is preserved for the assigned rows.
    assert table["x"].tolist() == [1.0, 3.0]


def test_unassigned_tokens_match_loader():
    assert "UNASSIGNED" in UNASSIGNED_TOKENS
    assert "cell_0" in UNASSIGNED_TOKENS
    assert "" in UNASSIGNED_TOKENS


def test_sha256_and_seed_id(tmp_path: Path):
    p = tmp_path / "f.bin"
    p.write_bytes(b"abc")
    assert sha256_file(p) == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
    assert seed_id(1) == "1"
    assert seed_id("a/b") == "a_b"


def test_audit_param_defaults_documented():
    assert set(AUDIT_PARAM_DEFAULTS) == {
        "neighbor_k", "n_pool", "q_thresh", "min_excess",
        "min_target_cells", "min_reference_cells"}
