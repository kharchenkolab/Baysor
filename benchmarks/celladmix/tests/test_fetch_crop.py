"""Tests for the pancreas crop builder (window selection + contract output)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fetch_pancreas import build_crop, pick_window


def test_pick_window_densest_and_deterministic():
    # Two clusters of cells: 5 near (0,0), 20 near (100,0).
    cells = np.array([[float(i % 5), float(i // 5)] for i in range(5)]
                     + [[100.0 + (i % 5), float(i // 5)] for i in range(20)])
    x0, y0, x1, y1 = pick_window(cells, side=50.0, grid_step=10.0)
    assert (x1 - x0) == 50.0 and (y1 - y0) == 50.0
    # The window must cover the dense cluster (20 cells), not the sparse one.
    inside = ((cells[:, 0] >= x0) & (cells[:, 0] < x1)
              & (cells[:, 1] >= y0) & (cells[:, 1] < y1))
    assert inside.sum() == 20
    # Deterministic across calls and starting offsets.
    assert pick_window(cells, side=50.0, grid_step=10.0) == (x0, y0, x1, y1)


def test_pick_window_side_rounded_to_grid():
    cells = np.array([[0.0, 0.0], [30.0, 30.0]])
    x0, y0, x1, y1 = pick_window(cells, side=47.0, grid_step=10.0)
    assert (x1 - x0) == 50.0  # 47um rounds up to 5 grid steps


def test_build_crop_contract_output(tmp_path: Path):
    raw = tmp_path / "raw"
    raw.mkdir()
    # 12 cells in a dense block + one far away; 60 molecules, mixed genes.
    rng = np.random.default_rng(3)
    cells = []
    for i in range(12):
        cx, cy = (i % 4) * 10.0, (i // 4) * 10.0
        cells.append((f"cell{i + 1}", cx, cy))
    cells.append(("cell_far", 500.0, 500.0))
    pd.DataFrame(cells, columns=["cell_id", "x_centroid", "y_centroid"]).to_parquet(
        raw / "cells.parquet", index=False)

    n = 300
    tx = pd.DataFrame({
        "feature_name": rng.choice(
            ["GeneA", "GeneB", "NegControlProbe-1", "UnassignedCodeword-100"], n),
        "x_location": rng.uniform(0, 45, n),
        "y_location": rng.uniform(0, 45, n),
        "qv": rng.uniform(10, 40, n).astype(np.float32),
        "cell_id": rng.choice([f"cell{i + 1}" for i in range(12)] + [""], n),
    })
    tx.to_parquet(raw / "transcripts.parquet", index=False)

    out = tmp_path / "crop"
    meta = build_crop(
        raw / "transcripts.parquet", raw / "cells.parquet", out,
        crop_id="test_crop", side_um=50.0, max_molecules=10_000, qv_min=20.0)

    mol = pd.read_parquet(out / "molecules.parquet")
    # Contract columns and dtypes.
    assert {"x", "y", "gene", "qv", "cell_vendor"} <= set(mol.columns)
    assert mol["x"].dtype == np.float64 and mol["y"].dtype == np.float64
    assert mol["qv"].dtype == np.float32
    # Sorted by (y, x).
    assert mol[["y", "x"]].equals(mol[["y", "x"]].sort_values(["y", "x"]).reset_index(drop=True))
    # Control features and low-qv rows filtered.
    assert not mol["gene"].str.startswith(("NegControl", "UnassignedCodeword")).any()
    assert (mol["qv"] >= 20.0).all()
    # Window covers the dense block only.
    assert (mol["x"] < 50).all() and (mol["y"] < 50).all()

    # meta.json is complete and consistent with the parquet.
    saved = json.loads((out / "meta.json").read_text())
    assert saved["stats"]["n_molecules"] == len(mol)
    assert saved["stats"]["n_genes"] == mol["gene"].nunique()
    assert saved["tier"] == "quick"
    assert saved["filters"]["qv_min"] == 20.0
    assert saved["crop"]["bbox_um"] == [0.0, 0.0, 50.0, 50.0]


def test_build_crop_shrinks_to_budget(tmp_path: Path):
    raw = tmp_path / "raw"
    raw.mkdir()
    cells = pd.DataFrame({
        "cell_id": [f"c{i}" for i in range(4)],
        "x_centroid": [0.0, 50.0, 0.0, 50.0],
        "y_centroid": [0.0, 0.0, 50.0, 50.0],
    })
    cells.to_parquet(raw / "cells.parquet", index=False)
    n = 1000
    rng = np.random.default_rng(0)
    tx = pd.DataFrame({
        "feature_name": ["GeneA"] * n,
        "x_location": rng.uniform(0, 90, n),
        "y_location": rng.uniform(0, 90, n),
        "qv": np.full(n, 30.0, dtype=np.float32),
        "cell_id": ["c1"] * n,
    })
    tx.to_parquet(raw / "transcripts.parquet", index=False)

    out = tmp_path / "crop"
    meta = build_crop(
        raw / "transcripts.parquet", raw / "cells.parquet", out,
        crop_id="shrunk", side_um=90.0, max_molecules=200, qv_min=20.0,
        shrink_step_um=10.0, min_side_um=30.0)
    mol = pd.read_parquet(out / "molecules.parquet")
    assert len(mol) <= 200
    assert meta["stats"]["n_molecules"] == len(mol)
    w = meta["crop"]["bbox_um"]
    assert (w[2] - w[0]) == (w[3] - w[1]) < 90.0
