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


# ---------------------------------------------------------------------------
# extensions: ambient noise, anisotropic geometry, 3D, imperfect prior
# ---------------------------------------------------------------------------

def test_ambient_background_molecules():
    _repo_or_skip()
    df, meta = strec.generate(dataset_id="t_bg", tier="quick", packing=6000.0,
                              sigma=2.0, model_kind="disjoint", seed=961,
                              mean_tx=60.0, n_target=100, bg_frac=0.2)
    frac = float((df["cell"] == 0).mean())
    assert abs(frac - 0.2) < 0.01
    bg = df["cell"] == 0
    assert (df.loc[bg, "celltype"] == "").all()
    assert set(df["gene"]) <= set(df.loc[~bg, "gene"])  # genes from the panel
    # interior of background = the same margin rule as the cells
    assert df["interior"].dtype == bool and df["interior"].any()
    # prior is computed over every molecule; labels only ever 0 or a cell id
    n_cells = meta["stats"]["n_true_cells"]
    assert set(np.unique(df["prior"])) <= set(range(0, n_cells + 1))
    assert meta["truth"]["params"]["ambient_bg_frac"] == 0.2
    assert "excluded" in meta["truth"]["accuracy_subset"]
    assert meta["truth"]["oracle_accuracy"] is not None
    # truth columns for cell molecules unchanged vs the same field w/o bg
    df0, _ = strec.generate(dataset_id="t_bg0", tier="quick", packing=6000.0,
                            sigma=2.0, model_kind="disjoint", seed=961,
                            mean_tx=60.0, n_target=100, bg_frac=0.0)
    cols = ["x", "y", "gene", "cell"]
    got = df.loc[df["cell"] > 0, cols].reset_index(drop=True)
    want = df0[cols].reset_index(drop=True)          # df0 has no background
    pd.testing.assert_frame_equal(got, want)


def test_aniso_geometry():
    _repo_or_skip()
    df, meta = strec.generate(dataset_id="t_aniso", tier="quick", packing=6000.0,
                              sigma=2.0, model_kind="disjoint", seed=962,
                              mean_tx=60.0, n_target=100, geometry="aniso")
    assert meta["truth"]["params"]["geometry"] == "aniso"
    assert "aniso" in meta["difficulty"]["notes"]
    # the anisotropic label image really differs from the Voronoi one
    generator, _e, _o, _c = strec._import_strec()
    model = strec.build_model("disjoint")[0]
    f_v = generator.build_field(6000.0, 2.0, 962, n_target=100, model=model)
    f_a = generator.build_field(6000.0, 2.0, 962, n_target=100, model=model,
                                geometry="aniso")
    assert f_a.geometry == "aniso"
    assert not np.array_equal(f_a.label, f_v.label)
    assert f_a.label.shape == f_v.label.shape
    assert meta["stats"]["n_true_cells"] == f_a.n_cells
    assert df["cell"].max() <= f_a.n_cells


def test_3d_variant():
    _repo_or_skip()
    slab = 6.0
    df, meta = strec.generate(dataset_id="t_z", tier="quick", packing=6000.0,
                              sigma=2.0, model_kind="disjoint", seed=963,
                              mean_tx=60.0, n_target=100, z_slab_um=slab)
    assert "z" in df.columns and df["z"].dtype == np.float64
    assert meta["crop"]["z_range_um"] == [0.0, slab]
    assert meta["stats"]["volume_um3"] > 0
    # observed z carries displacement (can leave the slab), true z is inside
    assert df["z"].std() > 0
    t = meta["truth"]
    assert t["oracle_accuracy"] is None            # 2D oracle would mislead
    assert "oracle_accuracy_note" in t
    assert 0.0 <= t["naive_accuracy"] <= 1.0       # naive is nearest in 3D
    assert t["params"]["z_slab_um"] == slab
    # prior still uses only nucleus labels
    n_cells = meta["stats"]["n_true_cells"]
    assert set(np.unique(df["prior"])) <= set(range(0, n_cells + 1))
    assert (df["prior"] > 0).any()


def test_imprior_keeps_truth():
    _repo_or_skip()
    kw = dict(tier="quick", packing=6000.0, sigma=2.0, model_kind="disjoint",
              seed=964, mean_tx=60.0, n_target=100)
    df0, m0 = strec.generate(dataset_id="t_i", **kw)
    df1, m1 = strec.generate(dataset_id="t_i",
                             prior_opts={"kind": "imperfect", "seed": 7311,
                                         "base": "t_base"}, **kw)
    for col in ("x", "y", "gene", "cell", "interior", "celltype"):
        assert (df0[col].to_numpy() == df1[col].to_numpy()).all(), col
    assert (df1["prior"] != df0["prior"]).any()
    info = m1["truth"]["prior"]
    assert info["n_missed"] == round(0.2 * 100)
    assert info["n_merge_sources"] == round(0.05 * 100)
    assert info["base"] == "t_base"
    assert "prior" not in m0["truth"]


# ---------------------------------------------------------------------------
# large panels: build_model_from_xenium_h5 on a local 10x-style fixture
# ---------------------------------------------------------------------------

def _write_fixture_h5(path, n_genes=40, n_cells=400, seed=11):
    """Tiny 10x-style cell_feature_matrix.h5 (CSC features x cells)."""
    import h5py
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    X = rng.poisson(0.6, size=(n_genes, n_cells)).astype(np.float32)
    X[rng.random(X.shape) < 0.7] = 0.0
    X = sp.csc_matrix(X)
    names = [f"GENE{i:03d}" for i in range(n_genes - 5)] + [f"CTRL{i}" for i in range(5)]
    ftype = (["Gene Expression"] * (n_genes - 5)) + ["Control Probe"] * 5
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("data", data=X.data)
        m.create_dataset("indices", data=X.indices)
        m.create_dataset("indptr", data=X.indptr)
        m.create_dataset("shape", data=np.array(X.shape, dtype=np.int64))
        feat = m.create_group("features")
        feat.create_dataset("feature_type",
                            data=np.array([s.encode() for s in ftype]))
        feat.create_dataset("name",
                            data=np.array([s.encode() for s in names]))
    return names


def test_build_model_from_xenium_h5_fixture(tmp_path):
    _repo_or_skip()
    h5 = tmp_path / "cell_feature_matrix.h5"
    names = _write_fixture_h5(h5)
    model = strec.build_model_from_xenium_h5(
        h5, n_genes=20, n_types=3, cell_subsample=300, seed=7,
        min_cells_per_type=5, name="fixture")
    assert model.n_genes == 20                     # panel subsetted
    assert set(model.gene_names) <= {n for n in names
                                     if n.startswith("GENE")}  # no controls
    assert 1 <= model.n_types <= 3
    assert np.allclose(model.composition.sum(axis=1), 1.0)
    assert abs(model.proportions.sum() - 1.0) < 1e-9
    assert (model.proportions > 0).all()
    # deterministic
    again = strec.build_model_from_xenium_h5(
        h5, n_genes=20, n_types=3, cell_subsample=300, seed=7,
        min_cells_per_type=5, name="fixture")
    assert again.gene_names == model.gene_names
    assert np.array_equal(again.composition, model.composition)
    # full panel is kept when n_genes >= available genes
    full = strec.build_model_from_xenium_h5(
        h5, n_genes=None, n_types=3, cell_subsample=300, seed=7,
        min_cells_per_type=5, name="fixture")
    assert full.n_genes == 35                     # Gene Expression rows only


def test_build_model_unknown_kind():
    with pytest.raises(KeyError, match="prime5k1000"):
        strec.build_model("nope")
    with pytest.raises(ValueError, match="100..10000"):
        strec.build_model("prime5k42")
