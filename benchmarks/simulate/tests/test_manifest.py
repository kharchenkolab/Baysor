"""Tests for the sim.yaml manifest and generate_all.py dispatch/verification."""
import inspect
import re

import pytest
import yaml

import common
import generate_all
import strec
import trivial

MANIFEST = generate_all.DEFAULT_MANIFEST
_HEX64 = re.compile(r"^[0-9a-f]{64}$")


@pytest.fixture(scope="module")
def manifest():
    return generate_all.load_manifest(MANIFEST)


@pytest.fixture(scope="module")
def index(manifest):
    return generate_all.validate_manifest(manifest)


def test_manifest_structure(manifest):
    assert manifest["version"] == 1
    assert {"trivial", "strec", "noprior"} <= set(manifest["generators"])
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


def test_every_entry_records_sha256(manifest):
    """Expected file hashes are committed for every dataset (review item)."""
    for d in manifest["datasets"]:
        h = d.get("sha256")
        assert h is not None, f"{d['id']}: no sha256 in the manifest"
        assert set(h) == {"molecules.parquet", "meta.json"}
        assert all(_HEX64.match(v) for v in h.values()), d["id"]


def test_on_disk_datasets_match_manifest_hashes(manifest):
    """If a dataset is present on disk, its files must match the manifest.

    Skipped for datasets that have not been generated on this machine.
    """
    out_root = common.data_root() / "sim"
    checked = 0
    for d in manifest["datasets"]:
        path = out_root / d["id"]
        if not (path / "molecules.parquet").exists():
            continue
        for name in ("molecules.parquet", "meta.json"):
            assert common.sha256(path / name) == d["sha256"][name], \
                f"{d['id']}/{name}: on disk differs from the manifest"
        checked += 1
    if checked == 0:
        pytest.skip("no datasets generated on this machine")


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
        assert set(p) <= {"packing", "sigma", "model", "mean_tx", "n_target",
                          "geometry", "bg_frac", "z_slab_um", "prior_opts"}
        assert p["packing"] in (1000, 2575, 6000, 13625)
        assert p["sigma"] in (1.0, 2.0, 3.0)
        assert p["model"] in ("disjoint", "merfish", "xenium",
                              "prime5k1000", "prime5k5000")
        assert p.get("geometry", "voronoi") in ("voronoi", "aniso")


def test_design_covers_all_factors(manifest):
    """The st-recoverability grid must cover each factor and the required
    dense + sigma=2 + realistic combination, without being cartesian."""
    strec_ds = [d for d in manifest["datasets"] if d["generator"] == "strec"]
    packings = {d["params"]["packing"] for d in strec_ds}
    sigmas = {d["params"]["sigma"] for d in strec_ds}
    models = {d["params"]["model"] for d in strec_ds}
    assert 1000 in packings                    # truly sparse (< 2500 cells/mm2)
    assert packings == {1000, 2575, 6000, 13625}
    assert sigmas == {1.0, 2.0, 3.0}
    assert "disjoint" in models
    assert models & {"merfish", "xenium"}                 # >=1 realistic model
    assert models & {"prime5k1000", "prime5k5000"}        # large panels
    assert any(d["params"]["packing"] == 13625 and d["params"]["sigma"] == 2.0
               and d["params"]["model"] in ("merfish", "xenium")
               for d in strec_ds)                          # required combo
    # extensions
    assert any(d["params"].get("bg_frac") == 0.05 for d in strec_ds)
    assert any(d["params"].get("bg_frac") == 0.2 for d in strec_ds)
    assert any(d["params"].get("geometry") == "aniso" for d in strec_ds)
    assert any(d["params"].get("z_slab_um") for d in strec_ds)
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
    # new datasets stay in the quick tier; only the two documented full-tier
    # datasets exist
    full_ids = {d["id"] for d in manifest["datasets"] if d["tier"] == "full"}
    assert full_ids == {"sim_tiled_distinct_g1000_full", "strec_dense_s2_xenium_full"}


def test_noprior_coverage(manifest):
    """No-prior variants: every quick trivial dataset and >=5 strec datasets
    covering sparse/medium/dense, several sigmas and realistic models."""
    trivial_ids = {d["id"] for d in manifest["datasets"]
                   if d["generator"] == "trivial" and d["tier"] == "quick"
                   and not d["id"].endswith("_imprior")}
    strec_ids = {d["id"] for d in manifest["datasets"]
                 if d["generator"] == "strec"}
    noprior = {d["id"]: d for d in manifest["datasets"]
               if d["generator"] == "noprior"}
    for base in trivial_ids:
        assert f"{base}_noprior" in noprior, f"missing noprior of {base}"
    strec_nop = [d for d in noprior.values() if d["base"] in strec_ids]
    assert len(strec_nop) >= 5
    base_by_id = {d["id"]: d for d in manifest["datasets"]}
    densities = {base_by_id[d["base"]]["params"]["packing"] for d in strec_nop}
    assert {2575, 6000, 13625} <= densities
    sigmas = {base_by_id[d["base"]]["params"]["sigma"] for d in strec_nop}
    assert sigmas == {1.0, 2.0, 3.0}
    models = {base_by_id[d["base"]]["params"]["model"] for d in strec_nop}
    assert "disjoint" in models and models & {"merfish", "xenium"}
    # bases are never themselves noprior; tiers/seeds match (validate checks)
    for d in noprior.values():
        assert d["base"] != d["id"]


def test_imprior_coverage(manifest):
    imp = [d for d in manifest["datasets"]
           if (d.get("params") or {}).get("prior_opts")]
    assert 3 <= len(imp) <= 4
    for d in imp:
        opts = d["params"]["prior_opts"]
        assert opts["kind"] == "imperfect"
        assert opts["miss_frac"] == 0.2
        assert list(opts["shift_um"]) == [1.0, 2.0]
        assert opts["merge_frac"] == 0.05
        assert opts["seed"] != d["seed"]      # independent RNG for the defects


def test_validate_rejects_unknown_base(manifest):
    bad = dict(manifest)
    bad["datasets"] = [dict(d) for d in manifest["datasets"]]
    bad["datasets"].append({"id": "x_noprior", "generator": "noprior",
                            "base": "does_not_exist", "tier": "quick", "seed": 1})
    with pytest.raises((KeyError, ValueError)):
        generate_all.validate_manifest(bad)


def test_build_entry_trivial_roundtrip(tmp_path, manifest):
    """A tiny manifest entry builds, writes and verifies against its hashes;
    a tampered file on disk is detected."""
    import generate_all as ga
    entry = {"id": "t_roundtrip", "generator": "trivial",
             "scenario": "circles_gaps", "tier": "quick", "seed": 5,
             "params": {"extent_um": 160.0, "n_genes": 50}}
    df, meta = ga.build(entry)
    assert meta["id"] == "t_roundtrip"
    out_root = tmp_path / "sim"
    entry["sha256"] = common.write_dataset(out_root / entry["id"], df, meta)
    result = ga.verify_one(entry, out_root)
    assert result["ok"], result
    # tamper: on-disk parquet no longer matches
    with open(out_root / entry["id"] / "molecules.parquet", "ab") as fh:
        fh.write(b"corrupt")
    result = ga.verify_one(entry, out_root)
    assert not result["ok"]
    assert any("on-disk" in m for m in result["mismatch"])


def test_build_noprior_from_base(tmp_path):
    """noprior = same molecules without the prior column, baysor.prior none."""
    base = {"id": "t_base_np", "generator": "trivial",
            "scenario": "circles_gaps", "tier": "quick", "seed": 5,
            "params": {"extent_um": 160.0, "n_genes": 50}}
    derived = {"id": "t_base_np_noprior", "generator": "noprior",
               "base": "t_base_np", "tier": "quick", "seed": 5}
    index = {e["id"]: e for e in (base, derived)}
    df0, meta0 = generate_all.build(base, index)
    df1, meta1 = generate_all.build(derived, index)
    assert "prior" in df0.columns and "prior" not in df1.columns
    assert list(df0.drop(columns=["prior"]).columns) == list(df1.columns)
    for col in df1.columns:
        assert (df0[col].to_numpy() == df1[col].to_numpy()).all(), col
    assert meta1["id"] == "t_base_np_noprior"
    assert meta1["baysor"]["prior"] == "none"
    assert meta1["truth"]["prior_variant_of"] == "t_base_np"
    assert meta1["stats"] == meta0["stats"]      # same molecules
    assert meta0["baysor"]["prior"] == "column"


def test_generate_all_rejects_duplicate_ids(tmp_path):
    dup = {"version": 1, "generators": {},
           "datasets": [{"id": "a"}, {"id": "a"}]}
    p = tmp_path / "dup.yaml"
    p.write_text(yaml.safe_dump(dup))
    with pytest.raises(ValueError, match="duplicate"):
        generate_all.load_manifest(p)


def test_update_hashes_roundtrip(tmp_path):
    """--update-hashes writes sha256 lines that yaml re-parses, preserving
    comments and existing entries."""
    text = """\
# header comment
version: 1
generators:
  trivial: {note: x}
datasets:
  # --- section ---
  - id: a
    generator: trivial
    tier: quick
    seed: 1
    params: {n_genes: 10}
  - id: b
    generator: trivial
    tier: quick
    seed: 2
"""
    p = tmp_path / "m.yaml"
    p.write_text(text)
    h = {name: "ab" * 32 for name in ("molecules.parquet", "meta.json")}
    changed = generate_all.update_hashes(p, {"a": h, "b": h})
    assert set(changed) == {"a", "b"}
    back = generate_all.load_manifest(p)
    assert {d["id"]: d.get("sha256") for d in back["datasets"]} == {"a": h, "b": h}
    out = p.read_text()
    assert "# header comment" in out and "# --- section ---" in out
    # replacement is idempotent
    assert generate_all.update_hashes(p, {"a": h, "b": h}) == []
