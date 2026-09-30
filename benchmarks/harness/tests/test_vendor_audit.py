"""Tests for vendor_audit.py (vendor segmentation scored on baseline typing)."""
from pathlib import Path

import numpy as np
import pandas as pd

import vendor_audit


def test_vendor_assignment_factorizes_labels(tmp_path):
    df = pd.DataFrame({
        "cell_vendor": ["V2", "", "V1", "V2", "0", None, "V1"],
    })
    mol = tmp_path / "molecules.parquet"
    df.to_parquet(mol, index=False)
    out = tmp_path / "vendor_assignment.parquet"
    stats = vendor_audit.vendor_assignment(mol, out)
    assert stats == {"n_cells": 3, "n_assigned": 5}
    got = pd.read_parquet(out)
    assert list(got.columns) == ["mol_index", "cell", "confidence"]
    # sorted unique non-empty labels -> 1..K ("0" sorts first); empty/None
    # -> 0 ("0" is reported by the validator as a data finding, factorized
    # here like any label so the audit sees exactly what cell_vendor holds)
    assert got["cell"].tolist() == [3, 0, 2, 3, 1, 0, 2]
    assert got["mol_index"].tolist() == list(range(7))
    assert got["confidence"].isna().all()


def test_audit_dataset_requires_typing(tmp_path):
    root = tmp_path
    (root / "baselines" / "mybase" / "real_x").mkdir(parents=True)
    res = vendor_audit.audit_dataset("real_x", root, "mybase", 1,
                                     tmp_path)
    assert res["status"] == "failed"
    assert "celltypes.parquet" in res["reason"]


def test_audit_dataset_missing_dataset(tmp_path):
    base = tmp_path / "baselines" / "mybase" / "real_y"
    base.mkdir(parents=True)
    (base / "celltypes.parquet").write_bytes(b"x")
    (base / "fixed_pairs.json").write_text("[]")
    res = vendor_audit.audit_dataset("real_y", tmp_path, "mybase", 1, tmp_path)
    assert res["status"] == "failed"
    assert "not found" in res["reason"]
