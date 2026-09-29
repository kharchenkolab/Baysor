"""Tests for the shared download helper (benchmarks/fetch/download.py).

The ``http_server`` fixture (tests/conftest.py) serves files with Range
support and scripted failure responses (429/5xx with Retry-After), so
retries, resume, hashing and atomicity are exercised end to end without
the network.

Run from the repository root:

    .deps/bench/bin/python -m pytest benchmarks/fetch/tests
"""

from __future__ import annotations

import sys
import time
from collections import deque
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import download as dl  # noqa: E402



# ---------------------------------------------------------------------------
# fetch_url
# ---------------------------------------------------------------------------


def test_fetch_url_basic(http_server):
    base, state = http_server
    payload = bytes(range(256)) * 512  # 128 KiB
    path = state.put("f.bin", payload)
    dst = Path(state.root).parent / "cache" / "f.bin"

    got = dl.fetch_url(base + path, dst, expected_size=len(payload))
    assert dst.read_bytes() == payload
    assert got == state.digest(payload)
    assert not dst.with_name("f.bin.part").exists()
    assert state.hits[path] == 1


def test_fetch_url_reuses_verified_copy_offline(http_server):
    base, state = http_server
    payload = b"hello benchmark" * 100
    path = state.put("reuse.bin", payload)
    dst = Path(state.root).parent / "cache" / "reuse.bin"
    spec = dict(expected_size=len(payload), expected_sha256=state.digest(payload))

    dl.fetch_url(base + path, dst, **spec)
    hits = state.hits[path]
    got = dl.fetch_url(base + path, dst, **spec)  # verified reuse
    assert got == state.digest(payload)
    assert state.hits[path] == hits  # no network on reuse


def test_fetch_url_corrupt_same_size_copy_is_refetched(http_server):
    base, state = http_server
    payload = b"A" * 4096
    path = state.put("corrupt.bin", payload)
    dst = Path(state.root).parent / "cache" / "corrupt.bin"
    dl.fetch_url(base + path, dst, expected_sha256=state.digest(payload))
    dst.write_bytes(b"B" * 4096)  # same size, different bytes
    got = dl.fetch_url(base + path, dst, expected_sha256=state.digest(payload))
    assert dst.read_bytes() == payload
    assert got == state.digest(payload)
    assert state.hits[path] == 2


def test_fetch_url_resumes_part_with_range(http_server):
    base, state = http_server
    payload = bytes(range(256)) * 256  # 64 KiB
    path = state.put("resume.bin", payload)
    dst = Path(state.root).parent / "cache" / "resume.bin"
    part = dst.with_name("resume.bin.part")
    part.parent.mkdir(parents=True, exist_ok=True)
    part.write_bytes(payload[:1000])  # interrupted download

    got = dl.fetch_url(
        base + path, dst,
        expected_size=len(payload), expected_sha256=state.digest(payload),
    )
    assert state.ranges[path] == ["bytes=1000-"]  # true byte-range resume
    assert dst.read_bytes() == payload
    assert got == state.digest(payload)
    assert not part.exists()


def test_fetch_url_retries_429_honouring_retry_after(http_server):
    base, state = http_server
    payload = b"retry me"
    path = state.put("rate.bin", payload)
    state.script[path] = deque([
        (429, {"Retry-After": "1"}),
        (429, {"Retry-After": "1"}),
    ])
    dst = Path(state.root).parent / "cache" / "rate.bin"

    t0 = time.monotonic()
    got = dl.fetch_url(base + path, dst, backoff_s=0.05)
    elapsed = time.monotonic() - t0
    assert dst.read_bytes() == payload
    assert got == state.digest(payload)
    assert state.hits[path] == 3  # two scripted 429s, then success
    assert elapsed >= 1.8  # Retry-After: 1 honoured twice


def test_fetch_url_retries_5xx(http_server):
    base, state = http_server
    payload = b"flaky server"
    path = state.put("flaky.bin", payload)
    state.script[path] = deque([(503, {}), (500, {})])
    dst = Path(state.root).parent / "cache" / "flaky.bin"

    dl.fetch_url(base + path, dst, backoff_s=0.01)
    assert dst.read_bytes() == payload
    assert state.hits[path] == 3


def test_fetch_url_404_is_permanent(http_server):
    base, state = http_server
    dst = Path(state.root).parent / "cache" / "missing.bin"
    with pytest.raises(dl.DownloadError):
        dl.fetch_url(base + "/missing.bin", dst, backoff_s=0.01)
    assert state.hits["/missing.bin"] == 1  # no retries on 4xx
    assert not dst.exists()


def test_fetch_url_wrong_expected_sha_fails_fast(http_server):
    base, state = http_server
    payload = b"data"
    path = state.put("badsha.bin", payload)
    dst = Path(state.root).parent / "cache" / "badsha.bin"
    with pytest.raises(dl.DownloadError, match="sha256 mismatch"):
        dl.fetch_url(base + path, dst, expected_sha256="00" * 32, backoff_s=0.01)
    assert state.hits[path] == 1  # fresh attempt: permanent, not retried
    assert not dst.exists()
    assert not dst.with_name("badsha.bin.part").exists()


def test_fetch_url_wrong_expected_size_fails(http_server):
    base, state = http_server
    payload = b"data"
    path = state.put("badsize.bin", payload)
    dst = Path(state.root).parent / "cache" / "badsize.bin"
    with pytest.raises(dl.DownloadError, match="size mismatch"):
        dl.fetch_url(base + path, dst, expected_size=len(payload) + 7,
                     backoff_s=0.01)
    assert not dst.exists()


def test_fetch_url_resume_heals_poisoned_part(http_server):
    """A corrupt .part prefix is detected by the resumed sha256 and healed."""
    base, state = http_server
    payload = b"C" * 4096
    path = state.put("poison.bin", payload)
    dst = Path(state.root).parent / "cache" / "poison.bin"
    part = dst.with_name("poison.bin.part")
    part.parent.mkdir(parents=True, exist_ok=True)
    part.write_bytes(b"D" * 2048)  # wrong content: resume -> sha mismatch ->
    # part discarded, fresh full download succeeds
    got = dl.fetch_url(
        base + path, dst,
        expected_size=len(payload), expected_sha256=state.digest(payload),
        backoff_s=0.01,
    )
    assert dst.read_bytes() == payload
    assert got == state.digest(payload)


# ---------------------------------------------------------------------------
# fetch_zip_members
# ---------------------------------------------------------------------------


def test_fetch_zip_members_streams_and_records(http_server):
    base, state = http_server
    members = {
        "outs/transcripts.parquet": bytes(range(256)) * 1000,
        "outs/cells.parquet": b"cell,data" * 500,
    }
    path = state.make_zip("bundle.zip", members)
    cache = Path(state.root).parent / "cache"
    dsts = {m: cache / m for m in members}

    out = dl.fetch_zip_members(base + path, dsts)
    for m, data in members.items():
        assert dsts[m].read_bytes() == data
        assert out[m] == state.digest(data)
        assert not dsts[m].with_name(dsts[m].name + ".part").exists()

    # recorded expectations turn the second fetch into a fully offline verify
    expected = {m: {"bytes": len(d), "sha256": out[m]}
                for m, d in members.items()}
    hits = state.hits[path]
    out2 = dl.fetch_zip_members(base + path, dsts, expected=expected)
    assert out2 == out
    assert state.hits[path] == hits

    # same-size corruption is caught and refetched
    dsts["outs/cells.parquet"].write_bytes(b"X" * len(members["outs/cells.parquet"]))
    dl.fetch_zip_members(base + path, dsts, expected=expected)
    assert dsts["outs/cells.parquet"].read_bytes() == members["outs/cells.parquet"]
    assert state.hits[path] > hits


def test_fetch_zip_members_retries_transient_status(http_server):
    base, state = http_server
    members = {"a.txt": b"alpha"}
    path = state.make_zip("retry.zip", members)
    state.script[path] = deque([(503, {})])
    cache = Path(state.root).parent / "cache"

    out = dl.fetch_zip_members(base + path, {"a.txt": cache / "a.txt"},
                               backoff_s=0.01)
    assert (cache / "a.txt").read_bytes() == members["a.txt"]
    assert out["a.txt"] == state.digest(members["a.txt"])
    assert state.hits[path] >= 2  # one 503 then success


def test_fetch_zip_members_missing_member_raises(http_server):
    base, state = http_server
    path = state.make_zip("small.zip", {"a.txt": b"a"})
    cache = Path(state.root).parent / "cache"
    with pytest.raises(KeyError):
        dl.fetch_zip_members(base + path, {"nope.txt": cache / "nope.txt"})


# ---------------------------------------------------------------------------
# header parsing units
# ---------------------------------------------------------------------------


def test_retry_after_seconds():
    assert dl._retry_after_seconds(None) is None
    assert dl._retry_after_seconds("7") == 7.0
    assert dl._retry_after_seconds("garbage") is None
    from datetime import datetime, timedelta, timezone
    from email.utils import format_datetime
    when = datetime.now(timezone.utc) + timedelta(seconds=30)
    got = dl._retry_after_seconds(format_datetime(when))
    assert got is not None and 25 <= got <= 31


def test_sha256_file(tmp_path):
    p = tmp_path / "x"
    p.write_bytes(b"abc")
    import hashlib
    assert dl.sha256_file(p) == hashlib.sha256(b"abc").hexdigest()
