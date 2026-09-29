"""Shared fixtures for the BENCH-REALX fetch tests.

``http_server`` serves files from a temp directory with HTTP Range support
and scripted failure responses (429/5xx + Retry-After), so download retry,
resume and hashing logic can be exercised end to end without the network.
"""

from __future__ import annotations

import hashlib
import threading
from collections import defaultdict, deque
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest


class ServerState:
    """Server behaviour + small helpers for the tests (see ``http_server``)."""

    def __init__(self, root: Path):
        self.root = root
        self.hits = defaultdict(int)       # path -> request count
        self.ranges = defaultdict(list)    # path -> [Range header, ...]
        self.script = defaultdict(deque)   # path -> deque[(status, {hdr: val})]

    def put(self, name: str, data: bytes) -> str:
        """Write a file into the served root; return its request path."""
        (self.root / name).write_bytes(data)
        return f"/{name}"

    def make_zip(self, name: str, members: dict[str, bytes]) -> str:
        """Create a zip under the served root; return its request path."""
        import zipfile

        zpath = self.root / name
        with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as zf:
            for member, data in members.items():
                zf.writestr(member, data)
        return f"/{name}"

    @staticmethod
    def digest(data: bytes) -> str:
        return hashlib.sha256(data).hexdigest()


@pytest.fixture
def http_server(tmp_path):
    state = ServerState(tmp_path / "www")
    state.root.mkdir()

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *args):  # noqa: N802
            pass

        def do_HEAD(self):  # noqa: N802
            self._serve(head=True)

        def do_GET(self):  # noqa: N802
            self._serve(head=False)

        def _serve(self, head: bool) -> None:
            path = self.path.split("?", 1)[0]
            state.hits[path] += 1
            state.ranges[path].append(self.headers.get("Range"))
            script = state.script[path]
            if script:
                status, hdrs = script.popleft()
                self.send_response(status)
                for key, val in hdrs.items():
                    self.send_header(key, val)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            fp = state.root / path.lstrip("/")
            if not fp.is_file():
                self.send_response(404)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            data = fp.read_bytes()
            start, end = 0, len(data) - 1
            rng = self.headers.get("Range")
            if rng:
                spec = rng.removeprefix("bytes=")
                if spec.startswith("-"):  # suffix range (remotezip EOCD probe)
                    start = max(0, len(data) - int(spec[1:]))
                else:
                    lo, _, hi = spec.partition("-")
                    start = int(lo)
                    if hi:
                        end = min(int(hi), len(data) - 1)
                if start >= len(data) or start > end:
                    self.send_response(416)
                    self.send_header("Content-Range", f"bytes */{len(data)}")
                    self.send_header("Content-Length", "0")
                    self.end_headers()
                    return
                self.send_response(206)
                self.send_header("Content-Range", f"bytes {start}-{end}/{len(data)}")
            else:
                self.send_response(200)
            self.send_header("Content-Length", str(end - start + 1))
            self.end_headers()
            if not head:
                self.wfile.write(data[start:end + 1])

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}", state
    server.shutdown()
    server.server_close()
