"""Shared helpers for the BENCH-REALO fetch scripts (other.py).

Contains the non-trivial machinery used by the real-dataset builders:

* data-root / cache path resolution,
* HTTP downloads with size checks, including hosts behind proof-of-work
  challenges (Dryad's Anubis) and bot-walled hosts (pklab, fetched through a
  real headless Firefox via Playwright),
* deterministic dense-window selection for cropping,
* label-image crop / binary-mask connected-component labelling /
  polygon-to-label helpers,
* dataset-stat helpers (panel size / density classes per the contract).

Everything that is pure logic takes explicit inputs and is unit-tested in
``benchmarks/fetch/tests/test_other.py``.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import requests

USER_AGENT = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36"
)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

def repo_root() -> Path:
    """Root of the Baysor checkout (benchmarks/ is a direct child)."""
    return Path(__file__).resolve().parents[2]


def bench_data_root() -> Path:
    """Shared benchmark data root ($BAYSOR_BENCH_DATA or <repo>/.bench-data)."""
    env = os.environ.get("BAYSOR_BENCH_DATA")
    if env:
        return Path(env)
    return repo_root() / ".bench-data"


def cache_dir(group: str = "real_other") -> Path:
    p = bench_data_root() / "cache" / group
    p.mkdir(parents=True, exist_ok=True)
    return p


def real_dir() -> Path:
    p = bench_data_root() / "real"
    p.mkdir(parents=True, exist_ok=True)
    return p


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------

def http_session() -> requests.Session:
    s = requests.Session()
    s.headers["User-Agent"] = USER_AGENT
    return s


def ensure_file(
    url: str,
    dest: Path,
    *,
    session: requests.Session | None = None,
    expected_size: int | None = None,
    timeout: int = 300,
) -> Path:
    """Download ``url`` to ``dest`` unless it already exists with the right size.

    Downloads go to a ``.part`` file first so an interrupted transfer never
    leaves a truncated file behind.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        size = dest.stat().st_size
        if expected_size is None or size == expected_size:
            return dest
        if expected_size is not None and size != expected_size:
            raise RuntimeError(
                f"cached {dest} has size {size}, expected {expected_size}; delete to re-download"
            )
    sess = session or http_session()
    tmp = dest.with_suffix(dest.suffix + ".part")
    with sess.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                if chunk:
                    fh.write(chunk)
    if expected_size is not None and tmp.stat().st_size != expected_size:
        raise RuntimeError(f"download of {url} has wrong size: {tmp.stat().st_size}")
    tmp.replace(dest)
    return dest


def browser_download(url: str, dest: Path, timeout_ms: int = 600_000) -> Path:
    """Fetch a file through a real headless Firefox (Playwright).

    Used for hosts (pklab.med.harvard.edu) that serve an Incapsula JS
    challenge to non-browser clients.  Falls back to an in-context request
    once the challenge cookies exist.  Works for bodies larger than
    Playwright's ``APIResponse.body()`` string limit by using the download
    channel.
    """
    from playwright.sync_api import sync_playwright

    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    with sync_playwright() as p:
        browser = p.firefox.launch(headless=True)
        ctx = browser.new_context(accept_downloads=True)
        page = ctx.new_page()
        try:
            got: dict = {}

            def on_download(dl) -> None:
                got.setdefault("dl", dl)

            page.on("download", on_download)
            try:
                with page.expect_download(timeout=timeout_ms):
                    try:
                        page.goto(url, wait_until="commit", timeout=30_000)
                    except Exception:
                        # Navigation aborts once a download starts; expected.
                        pass
            except Exception:
                pass
            if "dl" in got:
                got["dl"].save_as(str(tmp))
            else:
                # Challenge page instead of a direct download: use the
                # now-authenticated request context.
                resp = ctx.request.get(url, timeout=timeout_ms)
                if resp.status != 200:
                    raise RuntimeError(f"browser fetch of {url} failed: HTTP {resp.status}")
                body = resp.body()
                if body[:200].lstrip().startswith(b"<"):
                    raise RuntimeError(f"browser fetch of {url} returned HTML, not the file")
                tmp.write_bytes(body)
        finally:
            browser.close()
    tmp.replace(dest)
    return dest


# ---------------------------------------------------------------------------
# Anubis proof-of-work (datadryad.org download protection)
# ---------------------------------------------------------------------------

def solve_anubis_pow(random_data: str, difficulty: int) -> tuple[int, str]:
    """Solve an Anubis ``fast`` challenge (pure function, unit-tested).

    Mirrors ``worker/sha256-purejs.mjs``: hash ``random_data + str(nonce)``
    with SHA-256 and require ``difficulty`` leading zero *bits* nibbles:
    ``floor(difficulty/2)`` zero bytes plus, for odd difficulties, a zero
    high nibble on the next byte.  Returns ``(nonce, hex_digest)``.
    """
    prefix_bytes = difficulty // 2
    odd = difficulty % 2
    nonce = 0
    while True:
        digest = hashlib.sha256((random_data + str(nonce)).encode("utf-8")).digest()
        if digest[:prefix_bytes] == b"\x00" * prefix_bytes:
            if not odd or (digest[prefix_bytes] >> 4) == 0:
                return nonce, digest.hex()
        nonce += 1


def anubis_authenticated_session(url: str, timeout: int = 120) -> requests.Session:
    """Return a session that passed the Anubis challenge guarding ``url``.

    The challenge page embeds ``<script id="anubis_challenge" type=
    "application/json">``; after solving it we call ``pass-challenge`` which
    sets the ``techaro.lol-anubis-auth`` cookie.
    """
    sess = http_session()
    t0 = time.time()
    r = sess.get(url, timeout=timeout)
    m = re.search(
        r'<script id="anubis_challenge" type="application/json">(.*?)</script>',
        r.text,
        re.S,
    )
    if m is None:
        if r.headers.get("content-type", "").startswith(("application/zip", "application/x-zip")):
            return sess  # no challenge (e.g. cookie still valid)
        raise RuntimeError(f"no Anubis challenge found at {url} (got {r.status_code})")
    payload = json.loads(m.group(1))
    difficulty = int(payload["rules"]["difficulty"])
    challenge = payload["challenge"]
    nonce, digest = solve_anubis_pow(challenge["randomData"], difficulty)
    resp = sess.get(
        f"{url.split('/stash/')[0]}/.within.website/x/cmd/anubis/api/pass-challenge",
        params={
            "id": challenge["id"],
            "response": digest,
            "nonce": nonce,
            "redir": url,
            "elapsedTime": int((time.time() - t0) * 1000),
        },
        timeout=timeout,
        allow_redirects=False,
    )
    if resp.status_code not in (302, 200) or "techaro.lol-anubis-auth" not in resp.cookies:
        raise RuntimeError(f"Anubis pass-challenge failed: HTTP {resp.status_code}")
    return sess


def dryad_zip_url(file_url: str, session: requests.Session | None = None) -> str:
    """Resolve a Dryad ``file_stream`` URL to its presigned S3 URL.

    The Dryad route is: Anubis challenge -> ``/stash/downloads/...`` 301 ->
    ``/downloads/...`` 302 -> presigned S3 URL (24 h validity).
    """
    sess = session or anubis_authenticated_session(file_url)
    r = sess.get(file_url.replace("/stash/downloads/", "/downloads/"), allow_redirects=False, timeout=120)
    if r.status_code not in (301, 302) or "location" not in r.headers:
        raise RuntimeError(f"could not resolve Dryad download URL: HTTP {r.status_code}")
    return r.headers["location"]


# ---------------------------------------------------------------------------
# Deterministic crop selection
# ---------------------------------------------------------------------------

def densest_window(
    x: np.ndarray,
    y: np.ndarray,
    target: int,
    cap: int,
    seed: int,
    *,
    max_rounds: int = 16,
) -> tuple[float, float, float, float]:
    """Pick a square window ``(x0, y0, x1, y1)`` with about ``target`` points.

    The side length starts at ``sqrt(target / global_density)`` and is
    re-fitted each round as ``side *= sqrt(target / n)`` (clipped to the
    frame and to a 0.5--2x change); every round the window position is
    re-chosen by maximising the point count over a coarse histogram grid
    (window spans 8 grid bins), with ties broken by a seeded RNG so the
    result is reproducible.  Among rounds with ``n <= cap`` the window whose
    count is closest to ``target`` is returned; if no round stayed under the
    cap, a final geometric shrink around the densest position is applied.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) == 0:
        raise ValueError("densest_window: empty input")
    x0f, x1f = float(x.min()), float(x.max())
    y0f, y1f = float(y.min()), float(y.max())
    frame_w = max(x1f - x0f, 1e-9)
    frame_h = max(y1f - y0f, 1e-9)
    density = len(x) / (frame_w * frame_h)
    side = math.sqrt(target / max(density, 1e-12))
    rng = np.random.default_rng(seed)

    def pick(s: float) -> tuple[tuple[float, float, float, float], int]:
        """Densest square window of side ``s`` (grid search + exact recount)."""
        s = min(s, frame_w, frame_h)
        bx = s / 8.0
        nx = max(int(math.ceil(frame_w / bx)), 1)
        ny = max(int(math.ceil(frame_h / bx)), 1)
        while nx * ny > 40_000_000:  # pathological guard
            bx *= 2
            nx = max(int(math.ceil(frame_w / bx)), 1)
            ny = max(int(math.ceil(frame_h / bx)), 1)
        hist, xedges, yedges = np.histogram2d(
            x, y, bins=[np.linspace(x0f, x1f + bx, nx + 1), np.linspace(y0f, y1f + bx, ny + 1)]
        )
        k = max(int(round(s / bx)), 1)
        integ = np.pad(hist, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
        sums = integ[k:, k:] - integ[:-k, k:] - integ[k:, :-k] + integ[:-k, :-k]
        if sums.size == 0:
            sums, k = hist, 1
        mx = sums.max()
        cand = np.argwhere(sums >= mx * 0.9999)
        ix, iy = cand[int(rng.integers(len(cand)))]
        wx0 = float(xedges[ix]) if ix < len(xedges) - 1 else x0f
        wy0 = float(yedges[iy]) if iy < len(yedges) - 1 else y0f
        wx0 = min(max(wx0, x0f), max(x1f - s, x0f))
        wy0 = min(max(wy0, y0f), max(y1f - s, y0f))
        wx1, wy1 = wx0 + s, wy0 + s
        m = (x >= wx0) & (x <= wx1) & (y >= wy0) & (y <= wy1)
        return (wx0, wy0, wx1, wy1), int(m.sum())

    feasible: list[tuple[tuple[float, float, float, float], int]] = []
    fallback: tuple[tuple[float, float, float, float], int] | None = None
    for _ in range(max_rounds):
        bbox, n = pick(side)
        if fallback is None or n > fallback[1]:
            fallback = (bbox, n)
        if n <= cap:
            feasible.append((bbox, n))
            if n >= target:
                break
            if n >= target * 0.75:
                break
        ratio = target / max(n, 1)
        side *= min(max(math.sqrt(ratio), 0.5), 2.0)
        side = min(side, frame_w, frame_h)
    if feasible:
        bbox, _ = min(feasible, key=lambda bn: abs(bn[1] - target))
        return bbox
    # Nothing fit under the cap: geometric shrink around the densest position.
    assert fallback is not None
    bbox, n = fallback
    side = bbox[2] - bbox[0]
    for _ in range(60):
        if n <= cap:
            break
        side *= min(math.sqrt(cap / n) * 0.99, 0.99)
        bbox, n = pick(side)
    return bbox


# ---------------------------------------------------------------------------
# Gene filters
# ---------------------------------------------------------------------------

def gene_mask(genes: pd.Series, patterns: Sequence[str]) -> pd.Series:
    """True where the gene survives (i.e. matches none of the fnmatch patterns)."""
    if not patterns:
        return pd.Series(True, index=genes.index)
    keep = pd.Series(True, index=genes.index)
    g = genes.astype(str)
    for pat in patterns:
        keep &= ~g.str.match(pat_to_regex(pat))
    return keep


def pat_to_regex(pat: str) -> str:
    """Translate an fnmatch-style pattern to a regex anchored at the start."""
    out = re.escape(pat)
    out = out.replace(r"\*", ".*").replace(r"\?", ".")
    return "^" + out + "$"


def strip_hybridization_suffix(genes: pd.Series) -> pd.Series:
    """Drop trailing ``_Hybridization<digits>`` suffixes (osmFISH loom naming)."""
    return genes.astype(str).str.replace(r"_Hybridization\d+$", "", regex=True)


# ---------------------------------------------------------------------------
# Label-image helpers
# ---------------------------------------------------------------------------

def cc_label_binary(mask: np.ndarray) -> np.ndarray:
    """Connected-component labels (uint16) of a binary mask; 0 stays 0."""
    from scipy import ndimage

    labels, _ = ndimage.label(mask > 0)
    if labels.max() > np.iinfo(np.uint16).max:
        raise RuntimeError("too many components for uint16 label image")
    return labels.astype(np.uint16)


def rasterize_polygons(
    polygons: dict[int, Sequence[Sequence[float]]],
    bbox: tuple[int, int, int, int],
    order: Sequence | None = None,
) -> np.ndarray:
    """Rasterise ``{label: [[x, y], ...]}`` into a label image cropped to
    ``bbox = (px0, py0, px1, py1)`` (pixel half-open bounds).

    Drawing order is ``order`` (or ascending label) so overlaps resolve
    deterministically.
    """
    from PIL import Image, ImageDraw

    px0, py0, px1, py1 = bbox
    w, h = int(px1 - px0), int(py1 - py0)
    if w <= 0 or h <= 0:
        raise ValueError("empty raster bbox")
    img = Image.new("I", (w, h), 0)
    dr = ImageDraw.Draw(img)
    labels = sorted(order) if order is not None else sorted(polygons)
    for label in labels:
        pts = polygons.get(label)
        if pts is None:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        if max(xs) < px0 or min(xs) > px1 or max(ys) < py0 or min(ys) > py1:
            continue
        shifted = [(float(p[0]) - px0, float(p[1]) - py0) for p in pts]
        dr.polygon(shifted, fill=int(label))
    arr = np.asarray(img, dtype=np.int32)
    if arr.min() < 0:
        raise RuntimeError("negative labels from rasterisation")
    return arr


def lookup_labels(
    polygons: dict[int, Sequence[Sequence[float]]],
    xs: np.ndarray,
    ys: np.ndarray,
) -> np.ndarray:
    """Label per query point via shapely STRtree ``contains`` (0 = outside)."""
    from shapely import points as shapely_points
    from shapely.geometry import Polygon
    from shapely.strtree import STRtree

    labels = sorted(polygons)
    geoms = []
    good_labels = []
    for lab in labels:
        pts = polygons[lab]
        if len(pts) < 3:
            continue
        poly = Polygon(pts)
        if poly.is_empty:
            continue
        geoms.append(poly)
        good_labels.append(lab)
    if not geoms:
        return np.zeros(len(xs), dtype=np.int32)
    tree = STRtree(geoms)
    # STRtree.query applies the predicate as predicate(input, tree_geometry).
    pts_idx, geom_idx = tree.query(
        shapely_points(np.asarray(xs, float), np.asarray(ys, float)), predicate="within"
    )
    out = np.zeros(len(xs), dtype=np.int32)
    # tree.query returns (input_index, tree_index)
    out[pts_idx] = [good_labels[g] for g in geom_idx]
    return out


# ---------------------------------------------------------------------------
# Stats / contract helpers
# ---------------------------------------------------------------------------

def gene_panel_class(n_genes: int) -> str:
    if n_genes < 50:
        return "tiny"
    if n_genes < 250:
        return "small"
    if n_genes < 700:
        return "medium"
    if n_genes < 2000:
        return "large"
    return "huge"


def cell_density_class(cells_per_mm2: float) -> str:
    if cells_per_mm2 < 2500:
        return "sparse"
    if cells_per_mm2 <= 7000:
        return "medium"
    return "dense"


def sort_molecules(df: pd.DataFrame) -> pd.DataFrame:
    """Contract ordering: sorted by (y, x)."""
    return df.sort_values(["y", "x"], kind="mergesort").reset_index(drop=True)


def write_molecules_parquet(df: pd.DataFrame, path: Path) -> None:
    """Write the contract molecule table (dictionary-encoded gene)."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(df, preserve_index=False)
    gene_idx = table.schema.get_field_index("gene")
    if gene_idx >= 0 and pa.types.is_string(table.schema.field(gene_idx).type):
        table = table.set_column(
            gene_idx, table.schema.field(gene_idx).with_type(pa.dictionary(pa.int32(), pa.string())),
            table.column(gene_idx).cast(pa.dictionary(pa.int32(), pa.string())),
        )
    pq.write_table(table, path)


def json_dump(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=False)
        fh.write("\n")
