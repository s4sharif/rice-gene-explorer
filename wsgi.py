"""Serve the existing read-only routes behind mandatory reviewer authentication."""

import base64
import binascii
import hmac
import io
import os
from http import HTTPStatus

from app import DB, Handler


class Response(Handler):
    # Adapt existing route responses to WSGI without opening a second HTTP server.
    def __init__(self, environ):
        self.path = environ.get("PATH_INFO", "/")
        if environ.get("QUERY_STRING"):
            self.path += "?" + environ["QUERY_STRING"]
        self.wfile = io.BytesIO()
        self.status = 200
        self.response_headers = []

    # Capture the status code produced by the shared application routes.
    def send_response(self, code, message=None):
        self.status = code

    # Capture application headers for the production server to transmit.
    def send_header(self, keyword, value):
        self.response_headers.append((keyword, str(value)))

    # Leave header transmission to the WSGI server after route dispatch completes.
    def end_headers(self):
        pass


# Require explicitly configured credentials before exposing any application content.
def credentials_ready():
    return (
        bool(os.environ.get("RICELENS_USERNAME"))
        and len(os.environ.get("RICELENS_PASSWORD", "")) >= 16
    )


# Compare reviewer credentials without logging them or accepting malformed authorization.
def authorized(environ):
    header = environ.get("HTTP_AUTHORIZATION", "")
    if not header.startswith("Basic ") or len(header) > 4096:
        return False
    try:
        decoded = base64.b64decode(header[6:], validate=True).decode("utf-8")
        username, password = decoded.split(":", 1)
    except (ValueError, UnicodeError, binascii.Error):
        return False
    user_ok = hmac.compare_digest(
        username.encode(), os.environ["RICELENS_USERNAME"].encode()
    )
    password_ok = hmac.compare_digest(
        password.encode(), os.environ["RICELENS_PASSWORD"].encode()
    )
    return user_ok & password_ok


# Return a bounded plain response for authentication, health checks and unsupported methods.
def plain_response(start_response, code, text, method, extra=()):
    body = text.encode()
    start_response(
        f"{code} {HTTPStatus(code).phrase}",
        [
            ("Content-Type", "text/plain; charset=utf-8"),
            ("Content-Length", str(len(body))),
            ("Cache-Control", "no-store"),
            ("X-Content-Type-Options", "nosniff"),
            ("X-Robots-Tag", "noindex, nofollow"),
            *extra,
        ],
    )
    return [b"" if method == "HEAD" else body]


# Protect HTML, assets, JSON and CSV equally while exposing only a minimal health response.
def application(environ, start_response):
    method = environ.get("REQUEST_METHOD", "GET")
    if not credentials_ready():
        return plain_response(
            start_response, 503, "Reviewer access is not configured.", method
        )
    if environ.get("PATH_INFO") == "/healthz" and method in {"GET", "HEAD"}:
        return plain_response(
            start_response,
            200 if DB.exists() else 503,
            "ready" if DB.exists() else "not ready",
            method,
        )
    if not authorized(environ):
        return plain_response(
            start_response,
            401,
            "Reviewer authentication required.",
            method,
            [("WWW-Authenticate", 'Basic realm="RiceLENS review", charset="UTF-8"')],
        )
    if method not in {"GET", "HEAD"}:
        return plain_response(
            start_response, 405, "Read-only service.", method, [("Allow", "GET, HEAD")]
        )
    response = Response(environ)
    response.do_GET()
    headers = response.response_headers + [
        ("X-Robots-Tag", "noindex, nofollow"),
        ("Referrer-Policy", "no-referrer"),
    ]
    start_response(f"{response.status} {HTTPStatus(response.status).phrase}", headers)
    return [b"" if method == "HEAD" else response.wfile.getvalue()]
