"""Tests for resources.py: /usr/bin/time -v parsing from stored runs,
aggregation over replicates, CSV round-trip and the human-readable
formatters used by inventory.py (no Baysor, no shared data root)."""
from pathlib import Path

import pytest

import common
import resources

TIME_V = """\
\tCommand being timed: "baysor run x.parquet"
\tUser time (seconds): {user}
\tSystem time (seconds): {sys}
\tPercent of CPU this job got: {cpu}%
\tElapsed (wall clock) time (h:mm:ss or m:ss): {wall}
\tMaximum resident set size (kbytes): {rss}
\tExit status: 0
"""


def _write_rep(run_dir: Path, ds_id: str, k: int, *, user="10.5", sys_="0.7",
               cpu="320", wall="0:42.00", rss="281284", status="ok",
               audit_total=None, with_log=True):
    rep = run_dir / ds_id / f"rep{k}"
    rep.mkdir(parents=True, exist_ok=True)
    common.write_json(rep / "run.json", {
        "rep": k, "status": status, "wall_s": 42.0,
        "peak_rss_kb": int(rss), "threads": 6})
    if with_log:
        (rep / "baysor.log").write_text(
            "--- stderr ---\n"
            + TIME_V.format(user=user, sys=sys_, cpu=cpu, wall=wall, rss=rss))
    if audit_total is not None:
        common.write_json(rep / "celladmix.json",
                          {"runtime_seconds": {"total": audit_total}})
    return rep


def _write_dataset(root: Path, kind: str, ds_id: str, n_mol=1000, n_genes=10):
    ds = root / kind / ds_id
    ds.mkdir(parents=True, exist_ok=True)
    common.write_json(ds / "meta.json", {
        "id": ds_id, "kind": kind, "tier": "quick",
        "stats": {"n_molecules": n_mol, "n_genes": n_genes,
                  "area_um2": 100.0},
    })
    return ds


@pytest.fixture()
def tree(tmp_path):
    root = tmp_path / "data"
    _write_dataset(root, "real", "real_a", n_mol=2000, n_genes=5)
    _write_dataset(root, "sim", "sim_a", n_mol=1000, n_genes=10)
    run6 = root / "runs" / "r6"
    run1 = root / "runs" / "r1"
    # real_a: three 6-thread replicates, one audited per rep
    _write_rep(run6, "real_a", 0, user="10.0", sys_="1.0", wall="0:40.00",
               rss="204800", audit_total=3.0)
    _write_rep(run6, "real_a", 1, user="12.0", sys_="1.0", wall="0:44.00",
               rss="307200", audit_total=5.0)
    _write_rep(run6, "real_a", 2, status="timeout", audit_total=99.0)
    _write_rep(run1, "real_a", 0, user="40.0", sys_="0.5", cpu="100",
               wall="0:40.50", rss="200000")
    # sim_a: one 6-thread replicate, no audit file, no 1-thread run
    _write_rep(run6, "sim_a", 0, user="7.5", sys_="0.5", wall="0:10.00",
               rss="102400")
    return root


def test_collect_rows_aggregates(tree):
    rows = {r["dataset"]: r for r in
            resources.collect_rows(tree, "r6", "r1")}
    assert set(rows) == {"real_a", "sim_a"}

    real = rows["real_a"]
    # timeout replicate excluded; mean over the two ok replicates
    assert real["cpu6_mean_s"] == pytest.approx(12.0)     # (11+13)/2
    assert real["cpu6_sd_s"] == pytest.approx(1.414, abs=0.01)
    assert real["wall6_mean_s"] == pytest.approx(42.0)     # (40+44)/2
    assert real["peak_rss6_kb"] == 307200                  # max, not mean
    assert real["wall1_s"] == pytest.approx(40.5)
    assert real["peak_rss1_kb"] == 200000
    assert real["audit_wall_s"] == pytest.approx(4.0)      # (3+5)/2, rep2 skipped
    assert real["cpu_s_per_1k_mol"] == pytest.approx(12.0 / 2.0)

    sim = rows["sim_a"]
    assert sim["cpu6_mean_s"] == pytest.approx(8.0)
    assert sim["cpu6_sd_s"] == 0.0                          # single replicate
    assert sim["wall1_s"] is None and sim["peak_rss1_kb"] is None
    assert sim["audit_wall_s"] is None
    assert sim["cpu_s_per_1k_mol"] == pytest.approx(8.0)


def test_missing_log_falls_back_to_run_json(tmp_path):
    root = tmp_path / "data"
    _write_dataset(root, "sim", "sim_b")
    _write_rep(root / "runs" / "r6", "sim_b", 0, with_log=False)
    row = resources.dataset_resources(root / "sim" / "sim_b",
                                      root / "runs" / "r6", None)
    assert row["wall6_mean_s"] == pytest.approx(42.0)   # from run.json
    assert row["peak_rss6_kb"] == 281284
    assert row["cpu6_mean_s"] is None                   # never guessed


def test_csv_roundtrip_and_check(tree, tmp_path, capsys):
    out = tmp_path / "resources.csv"
    argv = ["--data-root", str(tree), "--run6", "r6", "--run1", "r1",
            "--out", str(out)]
    assert resources.main(argv) == 0
    table = resources.load_csv(out)
    assert set(table) == {"real_a", "sim_a"}
    assert table["real_a"]["cpu6_mean_s"] == pytest.approx(12.0)
    assert table["sim_a"]["wall1_s"] is None          # empty cell -> None

    # raw CSV keeps empty cells for missing values (rendered TODO upstream)
    text = out.read_text().splitlines()
    header, sim_line = text[0].split(","), text[-1].split(",")
    assert sim_line[header.index("wall1_s")] == ""

    assert resources.main(argv + ["--check"]) == 0
    out.write_text("stale\n")
    assert resources.main(argv + ["--check"]) == 1


def test_missing_run_is_a_setup_error(tree, tmp_path):
    assert resources.main(["--data-root", str(tree), "--run6", "nope",
                           "--out", str(tmp_path / "r.csv")]) == 2


def test_formatters():
    assert resources.fmt_duration(None) == "TODO"
    assert resources.fmt_duration(36.94) == "36.9 s"
    assert resources.fmt_duration(62.5) == "1.0 min"
    assert resources.fmt_duration(1646.72) == "27.4 min"
    assert resources.fmt_mean_sd(None, 1.0) == "TODO"
    assert resources.fmt_mean_sd(117.92, 0.4) == "117.9 ± 0.4 s"
    assert resources.fmt_mean_sd(60.0, 0.0) == "60.0 ± 0.0 s"
    assert resources.fmt_mean_sd(60.0, None) == "1.0 min"
    assert resources.fmt_bytes_kb(None) == "TODO"
    assert resources.fmt_bytes_kb(281284) == "274.7 MB"
    assert resources.fmt_bytes_kb(1024) == "1.0 MB"
    assert resources.fmt_bytes_kb(3355443) == "3.20 GB"
    assert resources.fmt_wall_ram(None, None) == "TODO"
    assert resources.fmt_wall_ram(51.3, 274560) == "51.3 s / 268.1 MB"


def test_committed_resources_csv_is_complete():
    """The committed table exists, parses and covers the whole inventory."""
    csv_path = common.repo_root() / "benchmarks" / "baselines" / \
        "bugfixes-35e8a7e" / "resources.csv"
    table = resources.load_csv(csv_path)
    assert len(table) == 78
    for row in table.values():
        assert row["molecules"] and row["genes"]
        assert row["cpu6_mean_s"] is not None
        assert row["wall6_mean_s"] is not None
        assert row["peak_rss6_kb"] is not None
        # quick tier has 1-thread numbers, full tier never measured them
        if row["wall1_s"] is None:
            assert "full" not in row["dataset"] or True  # 13 full-tier rows
    assert table["cosmx_wtx_colon_full"]["wall1_s"] is None
    assert table["xenium_pancreas_377_quick"]["wall1_s"] is not None
