"""Tests for suites.py: manifest loading/validation, run-id grouping,
estimation from a resources CSV and the committed real manifests."""
import textwrap
from pathlib import Path

import pytest

import common
import resources
import suites

REPO = common.repo_root()
MANIFEST = REPO / "benchmarks" / "datasets" / "suites.yaml"


def _write_manifest(tmp_path, body: str) -> Path:
    p = tmp_path / "suites.yaml"
    p.write_text(textwrap.dedent(body))
    return p


# ---------------------------------------------------------------------------
# committed manifests
# ---------------------------------------------------------------------------

def test_real_suites_resolve():
    regular = suites.resolve("regular", str(MANIFEST))
    release = suites.resolve("release", str(MANIFEST))
    assert [s.name for s in regular.steps] == ["exact", "noise"]
    assert [s.name for s in release.steps] == ["quick6", "full6", "quick1"]

    exact, noise = regular.steps
    assert (exact.threads, exact.replicates, exact.expect) == (1, 1, "identical")
    assert exact.baseline == "bugfixes-35e8a7e-t1"
    assert exact.celladmix is False and exact.no_ami is True
    assert (noise.threads, noise.replicates, noise.expect) == (6, 1, "same")
    assert noise.baseline == "bugfixes-35e8a7e"
    assert noise.celladmix is True and noise.no_ami is True
    assert noise.celltypes_from == "bugfixes-35e8a7e"
    # the coverage list: 12 trivial (6 scenarios x prior/noprior), 2 strec,
    # one 3D sim, 6 platforms at quick tier, one huge-panel crop, one admix
    assert len(noise.datasets.split(",")) == 23

    quick6, full6, quick1 = release.steps
    assert (quick6.threads, quick6.replicates, quick6.timeout) == (6, 3, 1800)
    assert (full6.threads, full6.replicates, full6.timeout) == (6, 3, 5400)
    assert full6.datasets == "full" and quick6.datasets == "quick"
    assert (quick1.threads, quick1.replicates) == (1, 1)
    assert quick1.expect == "identical"
    # quick + full = the complete inventory, quick1 = the complete quick tier
    assert {s.group for s in release.steps} == {"noise", "exact"}


def test_group_run_ids_keep_bare_id_for_identical():
    regular = suites.resolve("regular", str(MANIFEST))
    ids = suites.group_run_ids(regular, "b260101120000")
    assert ids == {"exact": "b260101120000",
                   "noise": "b260101120000-noise"}
    release = suites.resolve("release", str(MANIFEST))
    ids = suites.group_run_ids(release, "b260101120000")
    assert ids["exact"] == "b260101120000"        # <= 17 chars: determinism
    assert ids["noise"] == "b260101120000-noise"
    with pytest.raises(ValueError):
        suites.group_run_ids(regular, "bad run id!")


def test_unknown_suite_lists_available():
    with pytest.raises(ValueError, match="regular, release"):
        suites.resolve("nope", str(MANIFEST))


def test_committed_resources_cover_release_suite():
    """Every dataset of the release suite has measured resources, and the
    regular suite is a subset of release (no guessed or stray entries)."""
    table = resources.load_csv(REPO / "benchmarks" / "baselines" /
                               "bugfixes-35e8a7e" / "resources.csv")
    regular = suites.resolve("regular", str(MANIFEST))
    release = suites.resolve("release", str(MANIFEST))
    reg_ids = set().union(*[set(s.datasets.split(","))
                            for s in regular.steps])
    assert all(step.datasets in ("quick", "full")
               for step in release.steps)
    assert len(table) == 78            # the full inventory: quick + full
    assert reg_ids <= set(table)
    for ds_id in reg_ids:
        assert table[ds_id]["wall6_mean_s"] is not None


# ---------------------------------------------------------------------------
# custom manifests: validation, expansion, estimation
# ---------------------------------------------------------------------------

CUSTOM = """
suites:
  t:
    description: test suite
    budget_minutes: 5
    resources: nowhere/resources.csv
    steps:
      - name: exact
        datasets: [ds_a]
        threads: 1
        replicates: 1
        celladmix: false
        expect: identical
        baseline: b1
      - name: noise
        group: shared
        datasets: quick
        threads: 6
        replicates: 2
        timeout: 60
        celladmix: false
        expect: same
        baseline: b2
      - name: noise2
        group: shared
        datasets: ds_b
        threads: 6
        replicates: 2
        timeout: 60
        celladmix: false
        expect: same
        baseline: b2
"""


def test_custom_manifest_grouping(tmp_path):
    path = _write_manifest(tmp_path, CUSTOM)
    t = suites.resolve("t", str(path))
    assert [s.name for s in t.steps] == ["exact", "noise", "noise2"]
    # group 'exact' keeps the bare run-id (identical step), 'shared' is suffixed
    ids = suites.group_run_ids(t, "run1")
    assert ids == {"exact": "run1", "shared": "run1-shared"}


def test_group_expect_mismatch_rejected(tmp_path):
    bad = CUSTOM.replace("expect: same\n        baseline: b2",
                         "expect: improved\n        baseline: b2", 1)
    path = _write_manifest(tmp_path, bad)
    with pytest.raises(ValueError, match="share one expect/baseline"):
        suites.resolve("t", str(path))


def test_invalid_expect_rejected(tmp_path):
    path = _write_manifest(tmp_path, CUSTOM.replace("expect: identical",
                                                    "expect: whatever"))
    with pytest.raises(ValueError, match="expect"):
        suites.resolve("t", str(path))


def test_step_defaults(tmp_path):
    path = _write_manifest(tmp_path, """
        suites:
          d:
            steps:
              - name: only
                datasets: quick
                expect: same
                baseline: b
        """)
    (s,) = suites.resolve("d", str(path)).steps
    assert (s.threads, s.replicates, s.timeout) == (6, 1, 1800)
    assert s.celladmix is True and s.no_ami is False
    assert s.group == "only" and s.celltypes_from is None
    assert suites.step_spec(s) == "quick"


def test_estimate_arithmetic(tmp_path):
    from fixtures import make_sim_dataset
    root = tmp_path / "data"
    make_sim_dataset(root / "sim" / "ds_a", ds_id="ds_a")
    make_sim_dataset(root / "sim" / "ds_b", ds_id="ds_b", tier="full")
    csv = tmp_path / "resources.csv"
    csv.write_text(
        "dataset,molecules,genes,cpu6_mean_s,cpu6_sd_s,wall6_mean_s,"
        "peak_rss6_kb,wall1_s,peak_rss1_kb,audit_wall_s,cpu_s_per_1k_mol\n"
        "ds_a,100,5,30.0,1.0,10.0,1000,20.0,900,,0.3\n"
        "ds_b,200,5,80.0,2.0,40.0,2000,,,,0.4\n")
    path = _write_manifest(tmp_path, CUSTOM)
    t = suites.resolve("t", str(path))
    est = suites.estimate(t, root, csv)
    by_step = {e["step"]: e for e in est["steps"]}
    # exact: 1 thread x 1 rep -> wall1_s of ds_a
    assert by_step["exact"]["wall_s"] == pytest.approx(20.0)
    assert by_step["exact"]["cpu_s"] == pytest.approx(20.0)
    # noise: quick tier (ds_a) and ds_b, 6 threads x 2 reps
    assert by_step["noise"]["wall_s"] == pytest.approx(10.0 * 2)
    assert by_step["noise"]["cpu_s"] == pytest.approx(30.0 * 2)
    assert by_step["noise2"]["wall_s"] == pytest.approx(40.0 * 2)
    assert est["total"]["wall_s"] == pytest.approx(20 + 20 + 80)
    assert est["total"]["resources_found"] is True


def test_estimate_marks_missing_resources(tmp_path):
    from fixtures import make_sim_dataset
    root = tmp_path / "data"
    make_sim_dataset(root / "sim" / "ds_a", ds_id="ds_a")
    make_sim_dataset(root / "sim" / "ds_b", ds_id="ds_b", tier="full")
    csv = tmp_path / "resources.csv"
    csv.write_text("dataset,molecules,genes,cpu6_mean_s,cpu6_sd_s,"
                   "wall6_mean_s,peak_rss6_kb,wall1_s,peak_rss1_kb,"
                   "audit_wall_s,cpu_s_per_1k_mol\n")
    t = suites.resolve("t", str(_write_manifest(tmp_path, CUSTOM)))
    est = suites.estimate(t, root, csv)
    exact = est["steps"][0]
    assert exact["wall_s"] == 0.0
    assert exact["missing_resources"] == ["ds_a"]
    assert est["total"]["resources_found"] is False


def test_plan_text_lists_steps_and_totals(tmp_path):
    from fixtures import make_sim_dataset
    root = tmp_path / "data"
    make_sim_dataset(root / "sim" / "ds_a", ds_id="ds_a")
    make_sim_dataset(root / "sim" / "ds_b", ds_id="ds_b", tier="full")
    t = suites.resolve("t", str(_write_manifest(tmp_path, CUSTOM)))
    csv = tmp_path / "resources.csv"
    csv.write_text("dataset,molecules,genes,cpu6_mean_s,cpu6_sd_s,"
                   "wall6_mean_s,peak_rss6_kb,wall1_s,peak_rss1_kb,"
                   "audit_wall_s,cpu_s_per_1k_mol\n"
                   "ds_a,100,5,30.0,1.0,10.0,1000,20.0,900,,0.3\n"
                   "ds_b,200,5,80.0,2.0,40.0,2000,,,,0.4\n")
    selections = {"exact": ["ds_a"], "noise": ["ds_a"], "noise2": ["ds_b"]}
    text = suites.plan_text(t, "run1", selections, resources_csv=csv,
                            root=root)
    assert "step exact: group=exact run-id=run1" in text
    assert "step noise: group=shared run-id=run1-shared" in text
    assert "expect=identical baseline=b1" in text
    assert "datasets (1): ds_a" in text
    assert "estimated total" in text


def test_membership_tolerates_partial_root(tmp_path):
    from fixtures import make_sim_dataset
    root = tmp_path / "data"
    make_sim_dataset(root / "sim" / "ds_a", ds_id="ds_a")
    # 'quick' matches ds_a; ds_b does not exist here and is skipped silently
    path = _write_manifest(tmp_path, CUSTOM)
    out = suites.membership(root, manifest_path=str(path))
    assert out["t"] == {"ds_a"}
