"""Unit tests for the runner: command construction, output normalization,
dataset selection and /usr/bin/time parsing (no Baysor binary required)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import celladmix as camix
import common
import run as runner
from fixtures import make_sim_dataset, make_real_dataset

REPO = common.repo_root()

FAKE_HELP_RUN = """\
Usage: baysor run [OPTIONS] coordinates [prior_segmentation]
Options:
  -c,--config TEXT            TOML file with configuration
  -x,--x-column TEXT          Name of x column (default: x)
  -y,--y-column TEXT
  -z,--z-column TEXT
  -g,--gene-column TEXT
  -m,--min-molecules-per-cell INT
  -s,--scale FLOAT
  --scale-std TEXT
  --prior-segmentation-confidence FLOAT
  -o,--output TEXT            Output directory
  --output-style TEXT         Output bundle style: legacy or parquet
"""

FAKE_PROBE = {
    "path": "baysor", "sha256": "abc", "version_info": "Baysor",
    "flags": {"output-style": True, "prior-segmentation-confidence": True},
}


# ---------------------------------------------------------------------------
# command construction
# ---------------------------------------------------------------------------

def test_build_command_minimal(tmp_path):
    ds = make_sim_dataset(tmp_path / "data" / "sim" / "ds1",
                          baysor_overrides={"prior": "none", "config": None})
    ds = common.load_dataset(ds)
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               FAKE_PROBE, REPO)
    assert cmd[:4] == ["/bin/baysor", "run", str(ds.molecules_path), "-x"]
    assert "-s" in cmd and cmd[cmd.index("-s") + 1] == "5.0"
    assert cmd[cmd.index("--scale-std") + 1] == "25%"
    assert cmd[cmd.index("-m") + 1] == "10"
    assert "--output-style" in cmd and "--prior-segmentation-confidence" not in cmd
    assert cmd[-2:] == ["-o", str(tmp_path / "seg")]
    assert "-c" not in cmd


def test_build_command_prior_column(tmp_path):
    ds = make_sim_dataset(tmp_path / "data" / "sim" / "dsc",
                          baysor_overrides={"prior": "column"})
    ds = common.load_dataset(ds)
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               FAKE_PROBE, REPO)
    assert cmd[3] == ":prior"
    assert cmd[cmd.index("--prior-segmentation-confidence") + 1] == "0.5"


def test_build_command_prior_image(tmp_path):
    img_dir = make_sim_dataset(tmp_path / "data" / "sim" / "dsi",
                               baysor_overrides={"prior": "image:prior.tif"})
    (img_dir / "prior.tif").write_bytes(b"II*\x00fake")
    ds = common.load_dataset(img_dir)
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               FAKE_PROBE, REPO)
    assert cmd[3] == str(img_dir / "prior.tif")


def test_build_command_missing_prior_image(tmp_path):
    ds = make_sim_dataset(tmp_path / "data" / "sim" / "dsx",
                          baysor_overrides={"prior": "image:missing.tif"})
    ds = common.load_dataset(ds)
    with pytest.raises(FileNotFoundError):
        runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                             FAKE_PROBE, REPO)


def test_build_command_config_relative_to_repo_root(tmp_path):
    ds = make_sim_dataset(tmp_path / "data" / "sim" / "dsk",
                          baysor_overrides={"config": "configs/xenium.toml"})
    ds = common.load_dataset(ds)
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               FAKE_PROBE, REPO)
    assert cmd[cmd.index("-c") + 1] == str(REPO / "configs" / "xenium.toml")
    assert (REPO / "configs" / "xenium.toml").is_file()


def test_build_command_extra_args_and_no_output_style(tmp_path):
    ds = make_sim_dataset(tmp_path / "data" / "sim" / "dse",
                          baysor_overrides={"extra_args": ["--iters", "5"]})
    ds = common.load_dataset(ds)
    probe = dict(FAKE_PROBE, flags={"output-style": False})
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               probe, REPO)
    assert "--output-style" not in cmd
    assert cmd[cmd.index("--iters") + 1] == "5"
    # extra args come before -o
    assert cmd.index("--iters") < cmd.index("-o")


def test_build_command_scale_factor(tmp_path):
    ds = make_sim_dataset(tmp_path / "data" / "sim" / "dss")
    ds = common.load_dataset(ds)
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               FAKE_PROBE, REPO, scale_factor=0.5)
    assert cmd[cmd.index("-s") + 1] == "2.5"


def test_build_command_extra_args_verbatim_after_config(tmp_path):
    """Xenium-shaped meta (BENCH-REALX): extra_args carries the column flags
    plus --qv-column/--unassigned-prior-label because configs/xenium.toml
    maps vendor column names and sets unassigned_label=UNASSIGNED. Every
    extra_args token must pass through verbatim, after -c, and each option
    must appear exactly once (CLI11 rejects repeated scalar options)."""
    extra = ["-x", "x", "-y", "y", "-g", "gene", "--qv-column", "qv",
             "--unassigned-prior-label", "0"]
    ds = make_sim_dataset(
        tmp_path / "data" / "sim" / "dsx",
        baysor_overrides={"config": "configs/xenium.toml",
                          "prior": "column", "extra_args": extra})
    ds = common.load_dataset(ds)
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               FAKE_PROBE, REPO)
    # extra_args appears verbatim, in order, entirely after -c
    c_idx = cmd.index("-c")
    positions = []
    pos = c_idx + 2  # start scanning after '-c <path>'
    for tok in extra:
        assert tok in cmd, f"{tok} missing from command"
        idx = cmd.index(tok, pos)
        positions.append(idx)
        pos = idx + 1
    assert all(i > c_idx for i in positions), "extra_args must follow -c"
    # exactly one occurrence of each column flag, with the dataset's value
    assert cmd.count("-x") == 1 and cmd[cmd.index("-x") + 1] == "x"
    assert cmd.count("-y") == 1 and cmd[cmd.index("-y") + 1] == "y"
    assert cmd.count("-g") == 1 and cmd[cmd.index("-g") + 1] == "gene"
    assert cmd[cmd.index("--unassigned-prior-label") + 1] == "0"
    assert cmd[cmd.index("--qv-column") + 1] == "qv"


def test_build_command_builder_yields_to_extra_args(tmp_path):
    """When extra_args provides a builder flag (even with a different value),
    the builder omits its own occurrence: exactly one survives."""
    ds = make_sim_dataset(
        tmp_path / "data" / "sim" / "dsy",
        baysor_overrides={"extra_args": ["-g", "other_column", "--plot"]})
    ds = common.load_dataset(ds)
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               FAKE_PROBE, REPO)
    assert cmd.count("-g") == 1
    assert cmd[cmd.index("-g") + 1] == "other_column"
    assert "--plot" in cmd
    # non-conflicting builder flags are still emitted
    assert cmd.count("-x") == 1 and cmd.count("-y") == 1


def test_build_command_has_z_column(tmp_path):
    ds_dir = make_sim_dataset(tmp_path / "data" / "sim" / "dsz")
    df = pd.read_parquet(ds_dir / "molecules.parquet")
    df["z"] = 1.0
    df.to_parquet(ds_dir / "molecules.parquet", index=False)
    ds = common.load_dataset(ds_dir)
    cmd = runner.build_command(Path("/bin/baysor"), ds, tmp_path / "seg",
                               FAKE_PROBE, REPO)
    assert cmd[cmd.index("-z") + 1] == "z"


# ---------------------------------------------------------------------------
# segmentation output parsing / alignment
# ---------------------------------------------------------------------------

def _write_seg(seg_dir: Path, df: pd.DataFrame):
    seg_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(seg_dir / "molecules.parquet", index=False)


def test_normalize_assignment_positional(tmp_path):
    ds_dir = make_sim_dataset(tmp_path / "ds")
    src = pd.read_parquet(ds_dir / "molecules.parquet")
    seg = src[["gene", "x", "y"]].copy()
    seg.insert(0, "cell", ["0" if i % 4 == 0 else f"cell_{i % 4}"
                           for i in range(len(seg))])
    seg["assignment_confidence"] = np.linspace(0, 1, len(seg))
    _write_seg(tmp_path / "seg", seg)
    out = common.normalize_assignment(ds_dir, tmp_path / "seg")
    assert list(out.columns) == ["mol_index", "cell", "confidence"]
    assert out["mol_index"].tolist() == list(range(len(src)))
    assert out["cell"].dtype == np.int64
    assert out["cell"].iloc[0] == 0
    assert out["cell"].iloc[1] == 1
    np.testing.assert_allclose(out["confidence"], np.linspace(0, 1, len(seg)))


def test_normalize_assignment_row_shuffled(tmp_path):
    """Fallback alignment must re-map reordered segmentation rows correctly."""
    ds_dir = make_sim_dataset(tmp_path / "ds")
    src = pd.read_parquet(ds_dir / "molecules.parquet")
    rng = np.random.default_rng(3)
    perm = rng.permutation(len(src))
    seg = src.iloc[perm][["gene", "x", "y"]].reset_index(drop=True)
    # encode the original index in the cell label so we can verify alignment
    seg.insert(0, "cell", [f"cell_{i + 1}" for i in perm])
    seg["assignment_confidence"] = 0.5
    _write_seg(tmp_path / "seg", seg)
    out = common.normalize_assignment(ds_dir, tmp_path / "seg")
    np.testing.assert_array_equal(out["cell"].to_numpy(),
                                  np.arange(1, len(src) + 1))


def test_normalize_assignment_legacy_csv(tmp_path):
    ds_dir = make_sim_dataset(tmp_path / "ds")
    src = pd.read_parquet(ds_dir / "molecules.parquet")
    seg = src[["gene", "x", "y"]].copy()
    seg.insert(0, "cell", [f"cell_{i % 3 + 1}" for i in range(len(seg))])
    seg["assignment_confidence"] = 0.9
    seg_dir = tmp_path / "seg"
    seg_dir.mkdir()
    seg.to_csv(seg_dir / "segmentation.csv", index=False)
    out = common.normalize_assignment(ds_dir, seg_dir)
    assert out["cell"].min() == 1


def test_normalize_assignment_loader_dropped_rows(tmp_path):
    """Baysor's loader may drop molecules (min_molecules_per_gene etc.):
    they get cell 0 and NaN confidence, all other rows stay aligned."""
    ds_dir = make_sim_dataset(tmp_path / "ds")
    src = pd.read_parquet(ds_dir / "molecules.parquet")
    keep = np.arange(len(src)) % 3 != 0          # drop every third molecule
    seg = src.iloc[keep][["gene", "x", "y"]].reset_index(drop=True)
    seg.insert(0, "cell", [f"cell_{i % 3 + 1}" for i in range(len(seg))])
    seg["assignment_confidence"] = 0.7
    _write_seg(tmp_path / "seg", seg)
    out = common.normalize_assignment(ds_dir, tmp_path / "seg")
    assert len(out) == len(src)
    dropped = np.where(~keep)[0]
    assert (out.loc[dropped, "cell"] == 0).all()
    assert out.loc[dropped, "confidence"].isna().all()
    assert out.loc[keep, "cell"].min() == 1
    assert out.loc[keep, "confidence"].notna().all()


def test_normalize_assignment_foreign_rows_raise(tmp_path):
    ds_dir = make_sim_dataset(tmp_path / "ds")
    src = pd.read_parquet(ds_dir / "molecules.parquet").iloc[:5]
    seg = src[["gene", "x", "y"]].copy()
    seg["x"] = seg["x"] + 1000.0                # coordinates not in input
    seg.insert(0, "cell", "cell_1")
    seg["assignment_confidence"] = 1.0
    _write_seg(tmp_path / "seg", seg)
    with pytest.raises(RuntimeError):
        common.normalize_assignment(ds_dir, tmp_path / "seg")


def test_cells_to_int():
    vals = ["0", "cell_12", 3, "cell_3", None, "", "weird", "weird"]
    out = common.cells_to_int(vals)
    assert out.tolist() == [0, 12, 3, 3, 0, 0, 1, 1]


# ---------------------------------------------------------------------------
# dataset selection
# ---------------------------------------------------------------------------

def test_select_datasets_tiers_and_globs(tmp_path):
    make_sim_dataset(tmp_path / "sim" / "quick_a", tier="quick")
    make_sim_dataset(tmp_path / "sim" / "full_b", tier="full")
    make_real_dataset(tmp_path / "real" / "quick_r", tier="quick")
    assert len(common.select_datasets(tmp_path, "quick")) == 2
    assert len(common.select_datasets(tmp_path, "full")) == 1
    assert len(common.select_datasets(tmp_path, "all")) == 3
    assert len(common.select_datasets(tmp_path, "quick", kind="real")) == 1
    ids = [d.id for d in common.select_datasets(tmp_path, "full_b")]
    assert ids == ["full_b"]
    globbed = [d.id for d in common.select_datasets(tmp_path, "quick_*")]
    assert sorted(globbed) == ["quick_a", "quick_r"]
    multi = [d.id for d in common.select_datasets(tmp_path, "quick_a,quick_r")]
    assert sorted(multi) == ["quick_a", "quick_r"]
    with pytest.raises(ValueError):
        common.select_datasets(tmp_path, "no_such_dataset")
    with pytest.raises(ValueError):
        common.select_datasets(tmp_path, "")


# ---------------------------------------------------------------------------
# /usr/bin/time -v parsing & subprocess timeout
# ---------------------------------------------------------------------------

TIME_V_SAMPLE = """\
\tCommand being timed: "baysor run x.parquet"
\tUser time (seconds): 10.5
\tSystem time (seconds): 0.7
\tPercent of CPU this job got: 320%
\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02.50
\tMaximum resident set size (kbytes): 1234567
\tExit status: 0
"""


def test_parse_time_v():
    out = runner.parse_time_v(TIME_V_SAMPLE)
    assert out["peak_rss_kb"] == 1234567
    assert out["wall_s"] == pytest.approx(62.5)
    assert out["exit_status"] == 0
    assert out["user_s"] == pytest.approx(10.5)
    assert out["sys_s"] == pytest.approx(0.7)
    assert out["cpu_percent"] == 320


def test_parse_time_v_hour_format():
    out = runner.parse_time_v("Elapsed (wall clock) time "
                              "(h:mm:ss or m:ss): 1:02:03.50\n")
    assert out["wall_s"] == pytest.approx(3723.5)


def test_parse_time_v_empty():
    out = runner.parse_time_v("no timing here")
    assert out == {"wall_s": None, "peak_rss_kb": None, "exit_status": None,
                   "user_s": None, "sys_s": None, "cpu_percent": None}


def test_parse_time_v_without_cpu_lines():
    # older / unusual time output: no CPU%% or user/sys lines -> None fields
    out = runner.parse_time_v("\tElapsed (wall clock) time "
                              "(h:mm:ss or m:ss): 0:09.99\n")
    assert out["wall_s"] == pytest.approx(9.99)
    assert out["user_s"] is None and out["sys_s"] is None
    assert out["cpu_percent"] is None


def test_execute_timeout():
    res = runner.execute([sys.executable, "-c", "import time; time.sleep(30)"],
                         {}, timeout=1.0)
    assert res["timed_out"] is True
    assert res["wall_s"] < 10


def test_execute_exit_code():
    res = runner.execute([sys.executable, "-c", "raise SystemExit(7)"], {}, None)
    assert res["returncode"] == 7
    assert res["timed_out"] is False


# ---------------------------------------------------------------------------
# vendor labels
# ---------------------------------------------------------------------------

def test_vendor_labels(tmp_path):
    ds_dir = make_real_dataset(tmp_path / "ds", with_vendor=True)
    df = pd.read_parquet(ds_dir / "molecules.parquet")
    labels = runner.vendor_labels(df)
    assert labels.dtype == np.int64
    assert labels[df["cell_vendor"] == ""].max() == 0
    assert set(np.unique(labels)) - {0} == {1, 2, 3, 4, 5, 6}


# ---------------------------------------------------------------------------
# provenance, skip-existing reuse, celladmix aggregation
# ---------------------------------------------------------------------------

def _args_ns(**kw):
    import argparse
    base = dict(threads=6, scale_factor=1.0, timeout=0, no_celladmix=False,
                celltypes_from=None, run_id="r", replicates=1)
    base.update(kw)
    return argparse.Namespace(**base)


def test_reusable_run_checks_status_binary_threads_scale(tmp_path):
    probe = {"sha256": "abc"}
    prev = {"status": "ok", "binary_sha256": "abc", "threads": 6,
            "scale_factor": 1.0}
    ok, why = runner.reusable_run(prev, probe, _args_ns())
    assert ok and why == ""
    ok, why = runner.reusable_run(dict(prev, status="failed"), probe, _args_ns())
    assert not ok and "status" in why
    ok, why = runner.reusable_run(dict(prev, binary_sha256="def"), probe,
                                  _args_ns())
    assert not ok and "sha256" in why
    ok, why = runner.reusable_run(dict(prev, threads=1), probe, _args_ns())
    assert not ok and "threads" in why
    ok, why = runner.reusable_run(dict(prev, scale_factor=0.9), probe,
                                  _args_ns())
    assert not ok and "scale factor" in why
    # run.json from before provenance was recorded cannot be reused
    ok, why = runner.reusable_run({"status": "ok", "threads": 6}, probe,
                                  _args_ns())
    assert not ok and "sha256" in why


def test_metrics_record_input_and_rep_provenance(tmp_path):
    from fixtures import make_run
    root = tmp_path / "data"
    ds_dir = make_sim_dataset(root / "sim" / "prov")
    cells = pd.read_parquet(ds_dir / "molecules.parquet")["cell"].to_numpy(np.int64)
    make_run(root, "runprov", ds_dir, [cells] * 3, threads=6)
    m = common.read_json(root / "runs" / "runprov" / "prov" / "metrics.json")
    assert m["inputs"]["molecules_sha256"] == \
        common.sha256_file(ds_dir / "molecules.parquet")
    assert m["inputs"]["meta_sha256"] == common.sha256_file(ds_dir / "meta.json")
    for rec in m["reps"]:
        assert rec["threads"] == 6
        assert rec["binary_sha256"] == "0" * 64
        assert rec["assignment_sha256"]


def test_celladmix_aggregation_statuses(tmp_path):
    """Metrics derived from the adapter's normalized rep blocks: ok rates are
    averaged, unavailable audits keep the rate None (never 0.0)."""
    from fixtures import make_run, write_assignment, FAKE_BINARY
    root = tmp_path / "data"
    ds_dir = make_real_dataset(root / "real" / "agg")
    cells = pd.read_parquet(ds_dir / "molecules.parquet")["cell"].to_numpy(np.int64)

    # all ok: mean/sd over the rep rates
    m = make_run(root, "ag1", ds_dir, [cells] * 3,
                 celladmix_rates=[0.05, 0.06, 0.07])
    cm = m["real"]["celladmix"]
    assert cm["status"] == camix.STATUS_OK
    assert cm["mean_total"] == pytest.approx(0.06)
    assert cm["sd_total"] == pytest.approx(0.01)
    assert cm["per_rep_total"] == [0.05, 0.06, 0.07]

    # one replicate's audit is unavailable (small crop / no evaluated pairs)
    ds = common.load_dataset(ds_dir)
    dirs, records = {}, []
    for k in range(3):
        rep_dir = root / "runs" / "ag2" / ds.id / f"rep{k}"
        write_assignment(rep_dir, cells)
        records.append({
            "rep": k, "status": "ok", "exit_code": 0, "wall_s": 1.0,
            "peak_rss_kb": 1, "command": [], "command_str": "",
            "assignment": f"rep{k}/assignment.parquet",
            "assignment_sha256": common.sha256_file(rep_dir / "assignment.parquet"),
            "threads": 6, "binary_sha256": "0" * 64,
        })
        dirs[k] = rep_dir
    records[0]["celladmix"] = {"status": "ok", "total_admixture_rate": 0.05,
                               "admixture_capable": True}
    records[1]["celladmix"] = {"status": camix.STATUS_UNAVAILABLE,
                               "reason": "audit status 'no_detected_pairs'"}
    records[2]["celladmix"] = {"status": "ok", "total_admixture_rate": 0.07,
                               "admixture_capable": True}
    out = runner.compute_dataset_metrics(
        ds, "ag2", lambda r: dirs[r["rep"]], records, dict(FAKE_BINARY), 6,
        typing={"mode": "quick_cluster", "baseline": None})
    cm = out["real"]["celladmix"]
    assert cm["status"] == "partial"
    assert cm["mean_total"] == pytest.approx(0.06)   # None not averaged as 0
    assert cm["per_rep_total"] == [0.05, None, 0.07]
    assert cm["admixture_capable"] is True
    assert cm["typing"]["mode"] == "quick_cluster"
    assert cm["per_rep"][1]["status"] == camix.STATUS_UNAVAILABLE
    assert cm["per_rep"][1]["reason"]

    # no audit scored at all -> unavailable, rate None
    for rec in records:
        rec["celladmix"] = {"status": camix.STATUS_UNAVAILABLE, "reason": "x"}
    out2 = runner.compute_dataset_metrics(ds, "ag3",
                                          lambda r: dirs[r["rep"]], records,
                                          dict(FAKE_BINARY), 6)
    cm2 = out2["real"]["celladmix"]
    assert cm2["status"] == camix.STATUS_UNAVAILABLE
    assert cm2["mean_total"] is None


# ---------------------------------------------------------------------------
# /usr/bin/time -v CPU fields end to end (stub binary, failed replicate)
# ---------------------------------------------------------------------------

def test_run_replicate_records_cpu_fields(tmp_path):
    stub = tmp_path / "baysor-stub"
    stub.write_text("#!/bin/sh\nexit 3\n")
    stub.chmod(0o755)
    root = tmp_path / "data"
    ds_dir = make_sim_dataset(root / "sim" / "cpu")
    ds = common.load_dataset(ds_dir)
    probe = {"sha256": "0" * 64, "path": str(stub), "label": "t",
             "flags": {"output-style": False}}
    args = _args_ns(threads=2, run_id="rcpu", scale_factor=1.0, timeout=30,
                    no_celladmix=True, celltypes_from=None)
    rep_dir = root / "runs" / "rcpu" / "cpu" / "rep0"
    rep_dir.mkdir(parents=True)
    rec = runner.run_replicate(ds, 0, rep_dir, stub, probe, args,
                               common.repo_root(), root)
    assert rec["status"] == "failed"          # stub exits 3
    # the time -v block is parsed and recorded regardless of the exit status
    assert rec["cpu_user_s"] is not None and rec["cpu_user_s"] >= 0
    assert rec["cpu_sys_s"] is not None and rec["cpu_sys_s"] >= 0
    assert rec["cpu_percent"] is not None and rec["cpu_percent"] >= 0
    on_disk = common.read_json(rep_dir / "run.json")
    assert on_disk["cpu_user_s"] == rec["cpu_user_s"]
    assert on_disk["cpu_percent"] == rec["cpu_percent"]


# ---------------------------------------------------------------------------
# suite mode / dry-run / --no-ami
# ---------------------------------------------------------------------------

def _suite_files(tmp_path):
    root = tmp_path / "data"
    make_sim_dataset(root / "sim" / "sim_cpu", ds_id="sim_cpu")
    mp = tmp_path / "suites.yaml"
    mp.write_text(
        "suites:\n"
        "  t:\n"
        "    description: fixture suite\n"
        "    steps:\n"
        "      - name: exact\n"
        "        datasets: [sim_cpu]\n"
        "        threads: 1\n"
        "        replicates: 1\n"
        "        celladmix: false\n"
        "        no_ami: true\n"
        "        expect: identical\n"
        "        baseline: b1\n"
        "      - name: noise\n"
        "        group: noise\n"
        "        datasets: quick\n"
        "        threads: 6\n"
        "        replicates: 2\n"
        "        timeout: 60\n"
        "        expect: same\n"
        "        baseline: b2\n")
    return root, mp


def test_suite_dry_run_resolves_without_baysor(tmp_path, capsys):
    root, mp = _suite_files(tmp_path)
    rc = runner.main(["--suite", "t", "--manifest", str(mp),
                      "--run-id", "rr", "--dry-run", "--data-root", str(root)])
    assert rc == 0
    out = capsys.readouterr().out
    # the identical group keeps the bare run-id, the other group is suffixed
    assert "step exact: group=exact run-id=rr" in out
    assert "step noise: group=noise run-id=rr-noise" in out
    assert "threads=1 replicates=1" in out
    assert "threads=6 replicates=2" in out
    assert "timeout=60s" in out
    assert "ami=skipped" in out              # exact step carries no_ami
    assert "datasets (1): sim_cpu" in out
    assert not (root / "runs").exists()      # nothing was executed


def test_suite_rejects_cli_overrides(tmp_path):
    root, mp = _suite_files(tmp_path)
    with pytest.raises(SystemExit) as exc:
        runner.main(["--suite", "t", "--manifest", str(mp), "--run-id", "rr",
                     "--dry-run", "--data-root", str(root), "--threads", "2"])
    assert exc.value.code == 2
    with pytest.raises(SystemExit) as exc:
        runner.main(["--suite", "t", "--manifest", str(mp), "--run-id", "rr",
                     "--dry-run", "--data-root", str(root),
                     "--datasets", "quick"])
    assert exc.value.code == 2


def test_suite_requires_baysor_unless_dry_run(tmp_path):
    root, mp = _suite_files(tmp_path)
    with pytest.raises(SystemExit) as exc:
        runner.main(["--suite", "t", "--manifest", str(mp), "--run-id", "rr",
                     "--data-root", str(root)])
    assert exc.value.code == 2


def test_unknown_suite_is_setup_error(tmp_path, capsys):
    root, mp = _suite_files(tmp_path)
    with pytest.raises(SystemExit) as exc:
        runner.main(["--suite", "nope", "--manifest", str(mp),
                     "--run-id", "rr", "--dry-run", "--data-root", str(root)])
    assert exc.value.code == 2
    assert "unknown suite" in capsys.readouterr().err


def test_step_requires_suite(tmp_path):
    root, mp = _suite_files(tmp_path)
    with pytest.raises(SystemExit) as exc:
        runner.main(["--run-id", "rr", "--datasets", "quick", "--step",
                     "exact", "--dry-run", "--data-root", str(root)])
    assert exc.value.code == 2


def test_single_dry_run(tmp_path, capsys):
    root = tmp_path / "data"
    make_sim_dataset(root / "sim" / "solo", ds_id="solo")
    rc = runner.main(["--run-id", "r", "--datasets", "quick", "--dry-run",
                      "--data-root", str(root)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "plan: run-id=r threads=6 replicates=1" in out
    assert "datasets (1): solo" in out
    assert not (root / "runs").exists()


def test_no_ami_skips_ami_metric(monkeypatch):
    """--no-ami works by zeroing metrics.AMI_MAX_LABELS (guard = always NaN);
    AMI is informational, no gate reads it."""
    import metrics as mmod
    rng = np.random.default_rng(0)
    pred = rng.integers(0, 8, 100).astype(np.int64)
    truth = rng.integers(0, 8, 100).astype(np.int64)
    assert not np.isnan(mmod.ami(pred, truth))
    monkeypatch.setattr(mmod, "AMI_MAX_LABELS", 0)
    assert np.isnan(mmod.ami(pred, truth))


def test_suite_run_executes_all_groups_with_stub(tmp_path, capsys):
    """Full --suite execution path (no real Baysor): every step runs into
    its group folder, writes _suite.json/selection/metrics and the CPU
    fields, and the failing replicates surface in the exit code."""
    root, mp = _suite_files(tmp_path)
    stub = tmp_path / "baysor-stub"
    stub.write_text("#!/bin/sh\nexit 3\n")
    stub.chmod(0o755)
    rc = runner.main(["--suite", "t", "--manifest", str(mp), "--run-id", "rr",
                      "--baysor", str(stub), "--data-root", str(root)])
    assert rc == 1                              # every replicate failed
    out = capsys.readouterr().out
    assert "== suite t / step exact" in out
    assert "== suite t / step noise" in out

    exact = root / "runs" / "rr"
    noise = root / "runs" / "rr-noise"
    assert common.read_json(exact / "_suite.json")["group"] == "exact"
    assert common.read_json(noise / "_suite.json")["group"] == "noise"

    rec = common.read_json(exact / "sim_cpu" / "rep0" / "run.json")
    assert rec["status"] == "failed"
    assert rec["cpu_user_s"] is not None        # time -v recorded anyway
    assert rec["threads"] == 1

    m_exact = common.read_json(exact / "sim_cpu" / "metrics.json")
    assert m_exact["failures"]
    assert m_exact["metric_options"]["ami"].startswith("skipped")
    m_noise = common.read_json(noise / "sim_cpu" / "metrics.json")
    assert "metric_options" not in m_noise      # noise step computes AMI
    sel = common.read_json(noise / "_selection.json")
    assert sel["datasets"] == ["sim_cpu"]
    assert sel["invocations"][-1]["threads"] == 6
    assert sel["invocations"][-1]["replicates"] == 2
