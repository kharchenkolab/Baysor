"""Tests for segmentation/dataset row alignment (align.py)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from align import align_segmentation


def make_dataset(tmp_path: Path) -> Path:
    path = tmp_path / "molecules.parquet"
    pd.DataFrame({
        "x": [1.0, 2.0, 3.0, 4.0],
        "y": [1.0, 2.0, 3.0, 4.0],
        "gene": ["g1", "g2", "g3", "g4"],
        "cell_vendor": ["v1", "", "v2", "v1"],
    }).to_parquet(path, index=False)
    return path


def test_align_happy_path_with_coordinates(tmp_path: Path):
    ds = make_dataset(tmp_path)
    seg = tmp_path / "seg.parquet"
    pd.DataFrame({
        "x": [1.0, 2.0, 3.0, 4.0],
        "y": [1.0, 2.0, 3.0, 4.0],
        "gene": ["g1", "g2", "g3", "g4"],
        "cell": ["cell_1", "0", "cell_2", "cell_1"],
    }).to_parquet(seg, index=False)
    labels, stats = align_segmentation(ds, seg, "cell")
    # string '0' is an unassigned token
    assert labels.tolist() == ["cell_1", "", "cell_2", "cell_1"]
    assert stats["n_cells"] == 2
    assert stats["n_unassigned"] == 1


def test_align_extra_unassigned_labels(tmp_path: Path):
    ds = make_dataset(tmp_path)
    seg = tmp_path / "seg.parquet"
    pd.DataFrame({"cell": ["a", "noise", "b", "b"]}).to_parquet(seg, index=False)
    labels, stats = align_segmentation(ds, seg, "cell", extra_unassigned={"noise"})
    assert labels.tolist() == ["a", "", "b", "b"]
    assert stats["n_unassigned"] == 1


def test_align_rejects_row_count_mismatch(tmp_path: Path):
    ds = make_dataset(tmp_path)
    seg = tmp_path / "seg.parquet"
    pd.DataFrame({"cell": ["a", "b"]}).to_parquet(seg, index=False)
    with pytest.raises(ValueError, match="same row order"):
        align_segmentation(ds, seg, "cell")


def test_align_rejects_reordered_coordinates(tmp_path: Path):
    ds = make_dataset(tmp_path)
    seg = tmp_path / "seg.parquet"
    pd.DataFrame({
        "x": [2.0, 1.0, 3.0, 4.0],   # rows swapped
        "y": [1.0, 2.0, 3.0, 4.0],
        "gene": ["g1", "g2", "g3", "g4"],
        "cell": ["a", "b", "c", "d"],
    }).to_parquet(seg, index=False)
    with pytest.raises(ValueError, match="do not match"):
        align_segmentation(ds, seg, "cell")
