"""Tests for the st-recoverability wrapper (strec.py).

The pinned repo and the reference downloads live under
``$BAYSOR_BENCH_DATA/cache/sim``; tests that need them skip when unavailable
(e.g. on a fresh machine without network).
"""
import json

import numpy as np
import pandas as pd
import pytest

import common
import strec


def _repo_or_skip():
    try:
        return strec.ensure_repo()
    except Exception as exc:  # no network, no cache
        pytest.skip(f"st-recoverability clone unavailable: {exc}")


def test_ensure_repo_pinned():
    repo = _repo_or_skip()
    import subprocess
    head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    assert head == strec.REPO_COMMIT
    for mod in ("generator.py", "expression.py", "oracle.py", "config.py"):
        assert (repo / "src" / mod).exists()


def test_repo_is_not_vendored():
    """The external code must stay in the cache, never in the repo."""
    repo_root = common.Path(__file__).resolve().parents[3]
    for mod in ("generator.py", "oracle.py"):
        assert not (repo_root / "benchmarks" / "simulate" / mod).exists()


def test_import_modules():
    _repo_or_skip()
    generator, expression, oracle, config = strec._import_strec()
    assert hasattr(generator, "build_field")
    assert hasattr(expression, "build_disjoint_model")
    assert hasattr(oracle, "oracle_assign")
    assert config.MASTER_SEED == 20260618


def test_generate_disjoint_contract():
    _repo_or_skip()
    df, meta = strec.generate(dataset_id="t_strec", tier="quick",
                              packing=6000.0, sigma=2.0, model_kind="disjoint",
                              seed=991, mean_tx=60.0, n_target=100)
    # contract
    assert list(df.columns) == ["x", "y", "gene", "prior", "cell", "interior",
                                "celltype"]
    assert df[["y", "x"]].equals(df[["y", "x"]].sort_values(["y", "x"],
                               kind="stable").reset_index(drop=True))
    assert df["prior"].dtype == np.int32 and df["cell"].dtype == np.int32
    assert df["cell"].min() >= 1                # no background in this generator
    assert meta["kind"] == "sim" and meta["tier"] == "quick"

    # truth block
    t = meta["truth"]
    assert t["external"]["commit"] == strec.REPO_COMMIT
    assert t["params"]["packing_cells_per_mm2"] == 6000.0
    assert t["params"]["sigma_um"] == 2.0
    assert t["seed"]["transcripts"] == t["seed"]["field"] + 1
    for key in ("oracle_accuracy", "naive_accuracy"):
        assert 0.0 <= t[key] <= 1.0
    assert t["oracle_accuracy"] >= t["naive_accuracy"] - 0.05

    # units + scale
    assert t["units"]["units"] == "um (verified)"
    assert 1.0 <= meta["baysor"]["scale_um"] <= 50.0
    assert abs(meta["baysor"]["scale_um"] - t["units"]["r_mean_um"]) < 1e-9
    assert meta["baysor"]["prior"] == "column"

    # stats consistency: density equals the requested packing
    assert abs(meta["stats"]["true_cells_per_mm2"] - 6000.0) < 1.0


def test_generate_deterministic():
    _repo_or_skip()
    kw = dict(tier="quick", packing=6000.0, sigma=1.0, model_kind="disjoint",
              seed=992, mean_tx=50.0, n_target=50)
    df1, m1 = strec.generate(dataset_id="t", **kw)
    df2, m2 = strec.generate(dataset_id="t", **kw)
    pd.testing.assert_frame_equal(df1, df2)
    assert m1 == m2


def test_prior_is_three_micron_nucleus():
    _repo_or_skip()
    df, meta = strec.generate(dataset_id="t", tier="quick", packing=6000.0,
                              sigma=1.0, model_kind="disjoint", seed=993,
                              mean_tx=50.0, n_target=50)
    # rebuild the field to get the true centres
    generator, _e, _o, _c = strec._import_strec()
    field = generator.build_field(6000.0, 1.0, 993, n_target=50,
                                  model=strec.build_model("disjoint")[0])
    prior, _ = common.nucleus_prior(df[["x", "y"]].to_numpy(),
                                    field.centers, strec.NUCLEUS_RADIUS_UM)
    assert np.array_equal(prior, df["prior"].to_numpy())


def test_merfish_model_cached():
    _repo_or_skip()
    try:
        model, info = strec.build_model("merfish")
    except Exception as exc:  # no network for the first download
        pytest.skip(f"MERFISH reference unavailable: {exc}")
    assert model.n_genes == 155
    assert model.n_types > 5
    assert np.allclose(model.composition.sum(axis=1), 1.0)
    assert "MERFISH" in info["reference"]


def test_xenium_model_cached():
    _repo_or_skip()
    try:
        model, info = strec.build_model("xenium")
    except Exception as exc:
        pytest.skip(f"Xenium reference unavailable: {exc}")
    assert 250 <= model.n_genes <= 700          # Rep1 launch panel, ~313 genes
    assert model.n_types >= 5
    assert np.allclose(model.composition.sum(axis=1), 1.0)
    assert "Xenium" in info["reference"]
