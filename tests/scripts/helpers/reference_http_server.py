#!/usr/bin/env python3
"""Serve a small reference fixture with controlled HTTP failures."""

from collections import Counter
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
import sys


payload = Path(sys.argv[1]).read_bytes()
attempts = Counter()


class ReferenceHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        mode = self.path.lstrip("/")
        attempts[mode] += 1
        failures = {"retry_once": 1, "retry_four": 4, "conflicting": 1}
        interrupted = attempts[mode] <= failures.get(mode, 0)
        etag = '"reference-fixture-v1"'
        if mode == "incorrect" or (mode == "conflicting" and interrupted):
            etag = '"different-reference"'

        body = payload
        if mode == "bad_checksum":
            body = b"X" + payload[1:]
        elif mode == "bad_size":
            body = payload[:-1]

        self.send_response(200)
        self.send_header("Content-Length", str(len(body)))
        if mode != "missing":
            self.send_header("ETag", etag)
        self.end_headers()
        self.wfile.write(body[: len(body) // 2] if interrupted else body)
        self.wfile.flush()
        self.close_connection = True

    def log_message(self, message, *args):
        pass


with HTTPServer(("127.0.0.1", 0), ReferenceHandler) as server:
    print(f"http://127.0.0.1:{server.server_port}", flush=True)
    server.serve_forever()
