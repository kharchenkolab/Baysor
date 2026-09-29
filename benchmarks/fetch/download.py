"""Shared download helpers for the Baysor benchmark suite (BENCH-*).

Every raw download of the suite should go through this module so that all
fetchers get the same robustness properties:

* **Retries with exponential backoff** for transient failures - connection
  errors, timeouts and HTTP 429/5xx - honouring the ``Retry-After`` header
  (seconds or HTTP-date, capped at :data:`MAX_RETRY_AFTER_S`).
* **HTTP ``Range`` resume** of partial ``<name>.part`` files: a download
  interrupted mid-file continues from the last byte written instead of
  restarting (``fetch_url``).
* **Atomic writes**: bytes go to ``<name>.part`` and are renamed onto the
  destination only after the size and sha256 checks pass.
* **sha256 per downloaded member**: computed while streaming, returned to the
  caller so it can be recorded in a manifest, and verified against the
  expected value on reuse (a same-size corrupt file is re-downloaded).
* **Size checks**: against the caller's ``expected_size`` and, when known,
  the response ``Content-Length``/``Content-Range``.

Two entry points cover the suite's needs:

``fetch_url(url, dst, ...)``
    direct file downloads (other-platforms fetchers need only this), and
``fetch_zip_members(url, {member: dst}, ...)``
    streaming single members out of a remote zip with :mod:`remotezip`
    (used by the Xenium fetcher; the zip itself is never downloaded whole).

Note on member resume: zip members are DEFLATE streams, so decompression can
only start at the beginning of the member; a partially written member is
therefore restarted (bounded by the same retry/backoff policy) instead of
resumed byte-wise - resuming the *compressed* range would need the full
compressed prefix anyway.  Direct files in ``fetch_url`` do a true byte-range
resume.

Example (manifest-driven fetcher)::

    from download import fetch_url, sha256_file

    expected = manifest.get("downloads", {}).get(name, {})
    digest = fetch_url(url, cache_dir / name,
                       expected_size=expected.get("bytes"),
                       expected_sha256=expected.get("sha256"))
    manifest.setdefault("downloads", {})[name] = {
        "bytes": (cache_dir / name).stat().st_size, "sha256": digest}

CLI::

    python benchmarks/fetch/download.py <url> <dst> [--expected-sha256 HEX]
"""

from __future__ import annotations

import argparse
import hashlib
import os
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Mapping

#: attempts after the first try (first attempt + this many retries)
DEFAULT_MAX_RETRIES = 5
#: base of the exponential backoff: sleep = backoff * 2**attempt
DEFAULT_BACKOFF_S = 2.0
#: upper bound applied to a server-provided Retry-After hint
MAX_RETRY_AFTER_S = 300.0
#: statuses that are retried (others: 4xx is permanent, handled as such)
RETRYABLE_STATUS = frozenset({429, 500, 502, 503, 504})
_STREAM_CHUNK = 8 << 20
_HASH_CHUNK = 1 << 20


class DownloadError(IOError):
    """Permanent download failure (HTTP 4xx, size/sha256 mismatch, ...)."""


class _Retryable(Exception):
    """Internal: transient failure worth another attempt."""

    def __init__(self, msg: str, retry_after: float | None = None):
        super().__init__(msg)
        self.retry_after = retry_after


# ---------------------------------------------------------------------------
# hashing
# ---------------------------------------------------------------------------


def sha256_file(path: Path | str, chunk: int = _HASH_CHUNK) -> str:
    """Streaming sha256 of a local file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# retry plumbing
# ---------------------------------------------------------------------------


def _retry_after_seconds(value: str | None) -> float | None:
    """Parse a Retry-After header (delta-seconds or HTTP-date)."""
    if value is None:
        return None
    try:
        return max(0.0, float(value))
    except ValueError:
        pass
    try:
        when = parsedate_to_datetime(value)
        if when.tzinfo is None:
            when = when.replace(tzinfo=timezone.utc)
        return max(0.0, (when - datetime.now(timezone.utc)).total_seconds())
    except Exception:
        return None


def _sleep_before_retry(
    attempt: int,
    backoff_s: float,
    retry_after: float | None = None,
    max_retry_after: float = MAX_RETRY_AFTER_S,
) -> float:
    """Sleep ``backoff * 2**attempt``, at least a valid Retry-After hint."""
    delay = backoff_s * (2 ** attempt)
    if retry_after is not None:
        delay = max(delay, min(retry_after, max_retry_after))
    if delay > 0:
        time.sleep(delay)
    return delay


def _retry_session(max_retries: int, backoff_s: float, user_agent: str):
    """A requests.Session whose adapter retries 429/5xx per request.

    Used for :mod:`remotezip`, which issues its own range GETs; the
    ``Retry-After`` header is honoured by urllib3 (``respect_retry_after_header``).
    """
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    try:
        retry = Retry(
            total=max_retries,
            connect=max_retries,
            read=max_retries,
            status=max_retries,
            backoff_factor=backoff_s,
            status_forcelist=sorted(RETRYABLE_STATUS),
            allowed_methods=frozenset({"GET", "HEAD"}),
            respect_retry_after_header=True,
            raise_on_status=False,
        )
    except TypeError:  # pragma: no cover - older urllib3
        retry = Retry(total=max_retries, backoff_factor=backoff_s,
                      status_forcelist=sorted(RETRYABLE_STATUS),
                      method_whitelist=frozenset({"GET", "HEAD"}),
                      raise_on_status=False)
    sess = requests.Session()
    adapter = HTTPAdapter(max_retries=retry)
    sess.mount("http://", adapter)
    sess.mount("https://", adapter)
    sess.headers["User-Agent"] = user_agent
    # byte-transparent transfers: Range offsets and sha256 only line up with
    # the stored file when the body is not content-encoded
    sess.headers["Accept-Encoding"] = "identity"
    return sess


# ---------------------------------------------------------------------------
# direct URL download with Range resume
# ---------------------------------------------------------------------------


def fetch_url(
    url: str,
    dst: Path | str,
    *,
    expected_size: int | None = None,
    expected_sha256: str | None = None,
    resume: bool = True,
    max_retries: int = DEFAULT_MAX_RETRIES,
    backoff_s: float = DEFAULT_BACKOFF_S,
    timeout: float = 120.0,
    session=None,
) -> str:
    """Download ``url`` to ``dst``; return the sha256 of the stored file.

    * existing ``dst`` is verified (size, then sha256 when an expected digest
      is given) and reused without touching the network on success;
    * bytes stream through ``dst.name + ".part"`` and are renamed onto ``dst``
      only after all checks pass;
    * transient failures (429/5xx with Retry-After, connection drops, short
      reads) are retried with exponential backoff; a partial ``.part`` file
      is resumed with an HTTP ``Range`` request unless ``resume=False``.
    """
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_name(dst.name + ".part")

    if dst.exists():
        size_ok = expected_size is None or dst.stat().st_size == expected_size
        if size_ok:
            digest = sha256_file(dst)
            if expected_sha256 is None or digest == expected_sha256.lower():
                return digest
        dst.unlink()  # wrong size or corrupt: fetch again

    # plain session: this function owns the retry loop so it can see the
    # Retry-After header itself; the adapter-level retries are only used for
    # remotezip (see fetch_zip_members).
    if session is None:
        import requests

        session = requests.Session()
        session.headers["User-Agent"] = "baysor-bench/1.0"
        session.headers["Accept-Encoding"] = "identity"
    sess = session
    attempt = 0
    while True:
        try:
            return _fetch_url_once(
                url, dst, part,
                expected_size=expected_size,
                expected_sha256=expected_sha256,
                resume=resume,
                timeout=timeout,
                session=sess,
            )
        except _Retryable as exc:
            if attempt >= max_retries:
                raise DownloadError(
                    f"download failed after {attempt + 1} attempts: {url}: {exc}"
                ) from exc
            _sleep_before_retry(attempt, backoff_s, exc.retry_after)
            attempt += 1


def _fetch_url_once(
    url: str,
    dst: Path,
    part: Path,
    *,
    expected_size: int | None,
    expected_sha256: str | None,
    resume: bool,
    timeout: float,
    session,
) -> str:
    import requests

    if not resume and part.exists():
        part.unlink()
    offset = part.stat().st_size if part.exists() else 0

    headers = {"Range": f"bytes={offset}-"} if offset else {}
    try:
        resp = session.get(url, stream=True, timeout=timeout, headers=headers)
    except (requests.RequestException, OSError) as exc:
        raise _Retryable(f"request failed: {exc}") from exc

    status = resp.status_code
    if status in RETRYABLE_STATUS or status >= 500:
        retry_after = _retry_after_seconds(resp.headers.get("Retry-After"))
        resp.close()
        raise _Retryable(f"HTTP {status}", retry_after=retry_after)
    if status == 416:
        resp.close()
        total_known = expected_size
        if total_known is not None and offset == total_known:
            return _finalize(part, dst, expected_size, expected_sha256,
                             hasher=None, resumed=offset > 0, label=url)
        # stale/oversized .part: discard and start over
        part.unlink(missing_ok=True)
        raise _Retryable("HTTP 416 with incomplete part file")
    if status not in (200, 206):
        resp.close()
        raise DownloadError(f"HTTP {status} for {url}")

    append = status == 206 and offset > 0
    if status == 200 and offset > 0:
        offset = 0  # server ignored the Range header: restart from scratch

    # expected total length from the response (not from the manifest: a
    # stale expected_size must fail fast in _finalize, not loop forever)
    total = None
    content_range = resp.headers.get("Content-Range")  # "bytes a-b/total"
    if content_range and "/" in content_range:
        tail = content_range.split("/")[-1]
        if tail.isdigit():
            total = int(tail)
    if total is None and resp.headers.get("Content-Length"):
        total = offset + int(resp.headers["Content-Length"])

    hasher = hashlib.sha256()
    try:
        with open(part, "ab" if append else "wb") as fh:
            if append:  # seed the digest with the bytes already on disk
                with open(part, "rb") as existing:
                    while True:
                        block = existing.read(_HASH_CHUNK)
                        if not block:
                            break
                        hasher.update(block)
                fh.seek(0, os.SEEK_END)
            try:
                for block in resp.iter_content(_STREAM_CHUNK):
                    if block:
                        fh.write(block)
                        hasher.update(block)
            except (requests.RequestException, OSError) as exc:
                raise _Retryable(f"stream interrupted: {exc}") from exc
            fh.flush()
            os.fsync(fh.fileno())
    finally:
        resp.close()

    size = part.stat().st_size
    if total is not None and size != total:
        raise _Retryable(f"short read: {size} != {total} bytes")
    return _finalize(part, dst, expected_size, expected_sha256,
                     hasher=hasher, resumed=append, label=url)


def _finalize(
    part: Path,
    dst: Path,
    expected_size: int | None,
    expected_sha256: str | None,
    hasher,
    *,
    resumed: bool,
    label: str = "",
) -> str:
    """Verify the completed ``.part`` file and rename it onto ``dst``.

    A size mismatch means the expectation (or the server content) is wrong -
    retrying cannot help, so it fails immediately.  A sha256 mismatch after a
    resumed attempt may be a poisoned ``.part`` prefix, so the part file is
    discarded and one fresh attempt follows; on a fresh attempt it is
    permanent (e.g. a stale manifest digest).
    """
    size = part.stat().st_size
    if expected_size is not None and size != expected_size:
        part.unlink(missing_ok=True)
        raise DownloadError(
            f"size mismatch for {label}: got {size}, want {expected_size} bytes"
        )
    digest = hasher.hexdigest() if hasher is not None else sha256_file(part)
    if expected_sha256 is not None and digest != expected_sha256.lower():
        part.unlink(missing_ok=True)
        if resumed:
            raise _Retryable(
                f"sha256 mismatch after resume for {label}: got {digest}, "
                f"want {expected_sha256}; restarting"
            )
        raise DownloadError(
            f"sha256 mismatch for {label}: got {digest}, want {expected_sha256}"
        )
    os.replace(part, dst)
    return digest


# ---------------------------------------------------------------------------
# remote zip members
# ---------------------------------------------------------------------------


def fetch_zip_members(
    url: str,
    members: Mapping[str, Path | str],
    *,
    expected: Mapping[str, Mapping] | None = None,
    max_retries: int = DEFAULT_MAX_RETRIES,
    backoff_s: float = DEFAULT_BACKOFF_S,
    timeout: float = 120.0,
) -> dict[str, str]:
    """Ensure the listed members of the remote zip at ``url`` are on disk.

    ``members`` maps member name -> destination path; ``expected`` may map
    member name -> ``{"bytes": int, "sha256": str}`` (typically from a
    manifest).  Members already present with the advertised size are kept:
    they are re-hashed and checked against ``expected["sha256"]`` when one is
    recorded, otherwise the computed digest is returned for the caller to
    record.  Missing members are streamed from the zip with range requests
    (the zip itself is never downloaded whole), verified, and renamed from a
    ``.part`` file.

    Returns ``{member: sha256}``.  Transient HTTP/network failures are
    retried with exponential backoff (429/5xx also honour Retry-After through
    the urllib3 retry adapter).
    """
    from remotezip import RemoteZip, RemoteZipError

    expected = expected or {}
    wanted = dict(members)  # de-duplicate, keep order
    out: dict[str, str] = {}

    # Fast offline path: members recorded with sizes and digests are verified
    # locally; anything unrecorded goes through the zip so its advertised size
    # is checked as well.
    pending: list[str] = []
    for m, dst in wanted.items():
        dst = Path(dst)
        exp = expected.get(m) or {}
        if "bytes" in exp and "sha256" in exp and dst.exists():
            if dst.stat().st_size != exp["bytes"]:
                dst.unlink()  # wrong size: refetch
                pending.append(m)
                continue
            digest = sha256_file(dst)
            if digest == str(exp["sha256"]).lower():
                out[m] = digest
                continue
            dst.unlink()  # same size but different bytes: refetch
        pending.append(m)
    if not pending:
        return out
    if not pending:
        return out

    retryable = (
        RemoteZipError,
        OSError,
        TimeoutError,
    )
    import requests

    retryable = retryable + (requests.RequestException,)

    attempt = 0
    while True:
        try:
            with RemoteZip(
                url,
                session=_retry_session(max_retries, backoff_s, "baysor-bench/1.0"),
            ) as rz:
                info = {i.filename: i.file_size for i in rz.infolist()}
                for m in pending:
                    if m not in info:
                        raise KeyError(f"member {m!r} not found in {url}")
                for m in pending:
                    dst = Path(wanted[m])
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    exp = expected.get(m) or {}
                    size = info[m]
                    if dst.exists():
                        if dst.stat().st_size != size:
                            dst.unlink()
                        else:
                            digest = sha256_file(dst)
                            if (not exp.get("sha256")
                                    or digest == str(exp["sha256"]).lower()):
                                out[m] = digest  # size-checked reuse
                                continue
                            dst.unlink()
                    out[m] = _stream_member(
                        rz, m, dst, size=size,
                        expected_sha256=exp.get("sha256"),
                        max_retries=max_retries, backoff_s=backoff_s,
                    )
            return out
        except KeyError:
            raise
        except DownloadError:
            raise
        except retryable as exc:
            if attempt >= max_retries:
                raise DownloadError(
                    f"zip member download failed after {attempt + 1} attempts: "
                    f"{url}: {exc}"
                ) from exc
            _sleep_before_retry(attempt, backoff_s)
            attempt += 1


def _stream_member(
    rz, member: str, dst: Path, *, size: int,
    expected_sha256: str | None,
    max_retries: int, backoff_s: float,
) -> str:
    """Stream one member to ``dst`` via ``.part`` with retries.

    DEFLATE members restart from the beginning on retry (see module
    docstring); the size and sha256 checks and the atomic rename still apply.
    """
    import requests

    from remotezip import RemoteZipError

    part = dst.with_name(dst.name + ".part")
    retryable = (RemoteZipError, OSError, TimeoutError,
                 requests.RequestException, _Retryable)
    last_exc: Exception | None = None
    for attempt in range(max_retries + 1):
        try:
            if part.exists():
                part.unlink()
            hasher = hashlib.sha256()
            with rz.open(member) as src, open(part, "wb") as fh:
                while True:
                    block = src.read(_STREAM_CHUNK)
                    if not block:
                        break
                    fh.write(block)
                    hasher.update(block)
                fh.flush()
                os.fsync(fh.fileno())
            got = part.stat().st_size
            if got != size:
                raise _Retryable(f"short read for {member}: {got} != {size}")
            digest = hasher.hexdigest()
            if expected_sha256 is not None and digest != expected_sha256.lower():
                part.unlink(missing_ok=True)
                raise DownloadError(
                    f"sha256 mismatch for {member!r}: got {digest}, "
                    f"want {expected_sha256}"
                )
            os.replace(part, dst)
            return digest
        except DownloadError:
            raise
        except retryable as exc:
            last_exc = exc
            if attempt >= max_retries:
                break
            _sleep_before_retry(attempt, backoff_s)
    part.unlink(missing_ok=True)
    raise DownloadError(f"failed to fetch member {member!r}: {last_exc}") from last_exc


# ---------------------------------------------------------------------------
# CLI (handy for quick checks and for other fetchers adopting this module)
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("url")
    ap.add_argument("dst", type=Path)
    ap.add_argument("--expected-size", type=int, default=None)
    ap.add_argument("--expected-sha256", default=None)
    ap.add_argument("--no-resume", action="store_true")
    ap.add_argument("--max-retries", type=int, default=DEFAULT_MAX_RETRIES)
    args = ap.parse_args(argv)
    digest = fetch_url(
        args.url, args.dst,
        expected_size=args.expected_size,
        expected_sha256=args.expected_sha256,
        resume=not args.no_resume,
        max_retries=args.max_retries,
    )
    print(f"{digest}  {args.dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
