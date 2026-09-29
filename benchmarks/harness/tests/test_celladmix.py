"""Tests for the optional cellAdmix integration (fake audit scripts)."""
import json
import sys

import celladmix as camix


def _make_repo(tmp_path, script_body):
    script = tmp_path / "benchmarks" / "celladmix" / "audit.py"
    script.parent.mkdir(parents=True)
    script.write_text(script_body)
    return tmp_path


def test_absent_when_module_missing(tmp_path):
    out = camix.run_audit(tmp_path / "mols.parquet", tmp_path / "a.parquet",
                          tmp_path / "out.json", repo=tmp_path)
    assert out["status"] == camix.STATUS_ABSENT


def test_format_command(tmp_path):
    cmd = camix.format_command("python", tmp_path / "audit.py",
                               tmp_path / "m.parquet", tmp_path / "a.parquet",
                               tmp_path / "o.json")
    assert "--image" not in cmd
    cmd = camix.format_command("python", tmp_path / "audit.py",
                               tmp_path / "m.parquet", tmp_path / "a.parquet",
                               tmp_path / "o.json", image=tmp_path / "dapi.tif")
    assert cmd[-2:] == ["--image", str(tmp_path / "dapi.tif")]


def test_successful_audit(tmp_path):
    body = """
import argparse, json
p = argparse.ArgumentParser()
p.add_argument("--molecules"); p.add_argument("--assignment")
p.add_argument("--out"); p.add_argument("--image", default=None)
a = p.parse_args()
json.dump({"total_admixture_rate": 0.123,
           "pairs": [{"pair": "A|B", "rate": 0.4}]}, open(a.out, "w"))
"""
    repo = _make_repo(tmp_path, body)
    out = camix.run_audit(tmp_path / "m.parquet", tmp_path / "a.parquet",
                          tmp_path / "o.json", repo=repo)
    assert out["status"] == camix.STATUS_OK
    assert out["total_admixture_rate"] == 0.123
    assert out["pairs"][0]["rate"] == 0.4


def test_audit_nonzero_exit(tmp_path):
    repo = _make_repo(tmp_path, "import sys; print('boom', file=sys.stderr); sys.exit(3)")
    out = camix.run_audit(tmp_path / "m.parquet", tmp_path / "a.parquet",
                          tmp_path / "o.json", repo=repo)
    assert out["status"] == camix.STATUS_FAILED
    assert out["returncode"] == 3
    assert "boom" in out["stderr_tail"]


def test_audit_bad_output(tmp_path):
    repo = _make_repo(tmp_path, """
import argparse, json
p = argparse.ArgumentParser()
p.add_argument("--molecules"); p.add_argument("--assignment")
p.add_argument("--out"); p.add_argument("--image", default=None)
a = p.parse_args()
json.dump({"something": 1}, open(a.out, "w"))
""")
    out = camix.run_audit(tmp_path / "m.parquet", tmp_path / "a.parquet",
                          tmp_path / "o.json", repo=repo)
    assert out["status"] == camix.STATUS_FAILED
    assert "total_admixture_rate" in out["reason"]


def test_audit_unparsable_json(tmp_path):
    repo = _make_repo(tmp_path, f"""
import argparse
p = argparse.ArgumentParser()
p.add_argument("--molecules"); p.add_argument("--assignment")
p.add_argument("--out"); p.add_argument("--image", default=None)
a = p.parse_args()
open(a.out, "w").write("not json {{")
""")
    out = camix.run_audit(tmp_path / "m.parquet", tmp_path / "a.parquet",
                          tmp_path / "o.json", repo=repo)
    assert out["status"] == camix.STATUS_FAILED
