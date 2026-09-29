"""Tests for the sim.yaml manifest and generate_all.py dispatch/verification."""
import inspect

import pytest
import yaml

import common
import generate_all
import strec
import trivial

MANIFEST = generate_all.DEFAULT_MANIFEST


@pytest.fixture(scope="module")
def manifest():
    return generate_all.load_manifest(MANIFEST)


def test_manifest_structure(manifest):
    assert manifest["version"] == 1
    assert {"trivial", "strec"} <= set(manifest["generators"])
    ext = manifest["generators"]["strec"]["external"]
    assert ext["commit"] == strec.REPO_COMMIT
    assert ext["license"] == "MIT"
    assert ext["imported_modules"] == strec.REPO_MODULES


def test_ids_unique_and_generators_known(manifest):
    ids = [d["id"] for d in manifest["datasets"]]
    assert len(ids) == len(set(ids))
    for d in manifest["datasets"]:
        assert d["generator"] in manifest["generators"]
        assert d["tier"] in ("quick", "full")
        assert isinstance(d["seed"], int)


def test_trivial_params_match_scenario_signature(manifest):
    for d in manifest["datasets"]:
        if d["generator"] != "trivial":
            continue
        assert d["scenario"] in trivial.SCENARIOS
        sig = inspect.signature(trivial.SCENARIOS[d["scenario"]])
        for key in (d.get("params") or {}):
            assert key in sig.parameters, f"{d['id']}: unknown param {key}"


def test_strec_params(manifest):
    for d in manifest["datasets"]:
        if d["generator"] != "strec":
            continue
        p = d["params"]
        assert set(p) <= {"packing", "sigma", "model", "mean_tx", "n_target"}
        assert p["packing"] in (2575, 6000, 13625)
        assert p["sigma"] in (1.0, 2.0, 3.0)
        assert p["model"] in ("disjoint", "merfish", "xenium")


def test_design_covers_all_factors(manifest):
    """The st-recoverability grid must cover each factor and the required
    dense + sigma=2 + realistic combination, without being cartesian."""
    strec_ds = [d for d in manifest["datasets"] if d["generator"] == "strec"]
    packings = {d["params"]["packing"] for d in strec_ds}
    sigmas = {d["params"]["sigma"] for d in strec_ds}
    models = {d["params"]["model"] for d in strec_ds}
    assert packings == {2575, 6000, 13625}
    assert sigmas == {1.0, 2.0, 3.0}
    assert "disjoint" in models
    assert models & {"merfish", "xenium"}                 # >=1 realistic model
    assert any(d["params"]["packing"] == 13625 and d["params"]["sigma"] == 2.0
               and d["params"]["model"] in ("merfish", "xenium")
               for d in strec_ds)                          # required combo
    # not the full cartesian product
    assert len(strec_ds) < 3 * 3 * 3


def test_panel_ablation_and_single_full_tier(manifest):
    trivial_ds = [d for d in manifest["datasets"] if d["generator"] == "trivial"]
    circles = {d["params"]["n_genes"] for d in trivial_ds
               if d["scenario"] == "circles_gaps"}
    tiled = {d["params"]["n_genes"] for d in trivial_ds
             if d["scenario"] == "tiled_distinct"}
    assert {100, 1000, 5000} <= circles
    assert {100, 1000, 5000} <= tiled
    full = [d for d in trivial_ds if d["tier"] == "full"]
    assert len(full) == 1
    assert full[0]["scenario"] == "tiled_distinct"
    assert full[0]["params"]["n_genes"] == 1000
    # every scenario present
    assert {d["scenario"] for d in trivial_ds} == set(trivial.SCENARIOS)


def test_build_entry_trivial_roundtrip(tmp_path, manifest):
    """A tiny manifest entry builds, writes and hash-verifies byte-identically."""
    import generate_all as ga
    entry = {"id": "t_roundtrip", "generator": "trivial",
             "scenario": "circles_gaps", "tier": "quick", "seed": 5,
             "params": {"extent_um": 160.0, "n_genes": 50}}
    df, meta = ga.build(entry)
    assert meta["id"] == "t_roundtrip"
    out_root = tmp_path / "sim"
    common.write_dataset(out_root / entry["id"], df, meta)
    result = ga.verify_one(entry, out_root)
    assert result["ok"], result


def test_generate_all_rejects_duplicate_ids(tmp_path):
    dup = {"version": 1, "generators": {},
           "datasets": [{"id": "a"}, {"id": "a"}]}
    p = tmp_path / "dup.yaml"
    p.write_text(yaml.safe_dump(dup))
    with pytest.raises(ValueError, match="duplicate"):
        generate_all.load_manifest(p)
