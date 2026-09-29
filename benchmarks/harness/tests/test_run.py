"""Unit tests for the runner: command construction, output normalization,
dataset selection and /usr/bin/time parsing (no Baysor binary required)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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
\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02.50
\tMaximum resident set size (kbytes): 1234567
\tExit status: 0
"""


def test_parse_time_v():
    out = runner.parse_time_v(TIME_V_SAMPLE)
    assert out["peak_rss_kb"] == 1234567
    assert out["wall_s"] == pytest.approx(62.5)
    assert out["exit_status"] == 0


def test_parse_time_v_hour_format():
    out = runner.parse_time_v("Elapsed (wall clock) time "
                              "(h:mm:ss or m:ss): 1:02:03.50\n")
    assert out["wall_s"] == pytest.approx(3723.5)


def test_parse_time_v_empty():
    out = runner.parse_time_v("no timing here")
    assert out == {"wall_s": None, "peak_rss_kb": None, "exit_status": None}


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
