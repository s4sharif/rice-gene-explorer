"""Verify reviewer access covers every serving route and fails closed."""

import base64
import json
import os
import unittest
from unittest.mock import patch
from urllib.parse import urlsplit
from wsgiref.util import setup_testing_defaults
from wsgiref.validate import validator

from backend import wsgi

USERNAME = "reviewer"
PASSWORD = "test-only-strong-password"


# Exercise the hosted application as a WSGI server would without opening a socket.
def request(path, credentials=None, method="GET"):
    environ = {}
    setup_testing_defaults(environ)
    parsed = urlsplit(path)
    environ.update(
        PATH_INFO=parsed.path, QUERY_STRING=parsed.query, REQUEST_METHOD=method
    )
    if credentials is not None:
        environ["HTTP_AUTHORIZATION"] = (
            "Basic " + base64.b64encode(credentials.encode()).decode()
        )
    captured = {}

    # Capture response status and headers for protocol and privacy assertions.
    def start_response(status, headers, exc_info=None):
        captured.update(status=status, headers=dict(headers))

    result = validator(wsgi.application)(environ, start_response)
    try:
        captured["body"] = b"".join(result)
    finally:
        result.close()
    return captured


class HostingTests(unittest.TestCase):
    # Configure disposable reviewer credentials without storing real passwords in tests.
    def setUp(self):
        self.environment = patch.dict(
            os.environ, RICELENS_USERNAME=USERNAME, RICELENS_PASSWORD=PASSWORD
        )
        self.environment.start()
        self.addCleanup(self.environment.stop)

    # Require authentication for the page, static resources, data APIs and every CSV export.
    def test_all_data_routes_are_private(self):
        for path in [
            "/",
            "/variation",
            "/gwas",
            "/expression",
            "/app.js",
            "/style.css",
            "/api/search?q=GS3",
            "/api/status",
            "/api/region?gene=Os03g0407400",
            "/api/variant-suggestions?gene=Os03g0407400&variant_query=chr3",
            "/api/export?gene=Os03g0407400",
            "/api/expression/datasets",
            "/api/expression/sources?gene=Os03g0407400",
            "/api/expression/export?genes=Os03g0407400",
            "/data/rice.sqlite",
        ]:
            with self.subTest(path=path):
                response = request(path)
                self.assertTrue(response["status"].startswith("401"))
                self.assertIn("WWW-Authenticate", response["headers"])
                self.assertNotIn(b"Os03g0407400", response["body"])

    # Serve directly bookmarked tool pages behind the same reviewer access controls.
    def test_tool_page_routes(self):
        for path in ["/", "/variation", "/gwas", "/expression"]:
            with self.subTest(path=path):
                response = request(path, USERNAME + ":" + PASSWORD)
                self.assertTrue(response["status"].startswith("200"))
                self.assertIn(b'aria-label="Main navigation"', response["body"])
                self.assertIn(b'id="cell-expression-page"', response["body"])

    # Serve the existing search API only when both reviewer credentials match exactly.
    def test_authorized_search_and_wrong_password(self):
        good = request("/api/search?q=GS3", USERNAME + ":" + PASSWORD)
        self.assertTrue(good["status"].startswith("200"))
        self.assertGreater(len(json.loads(good["body"])["matches"]), 1)
        self.assertTrue(request("/", USERNAME + ":wrong")["status"].startswith("401"))

    # Refuse application access when credentials are missing or too short.
    def test_unconfigured_service_fails_closed(self):
        for password in ["", "short"]:
            with patch.dict(os.environ, RICELENS_PASSWORD=password):
                self.assertTrue(request("/")["status"].startswith("503"))
                self.assertTrue(request("/healthz")["status"].startswith("503"))

    # Keep health output minimal and preserve correct read-only HEAD and method behavior.
    def test_health_head_and_write_methods(self):
        self.assertEqual(request("/healthz")["body"], b"ready")
        head = request("/", USERNAME + ":" + PASSWORD, "HEAD")
        self.assertEqual(head["body"], b"")
        self.assertGreater(int(head["headers"]["Content-Length"]), 0)
        self.assertTrue(
            request("/", USERNAME + ":" + PASSWORD, "POST")["status"].startswith("405")
        )

    # Reject malformed authentication data without returning an application error or dataset.
    def test_malformed_authorization(self):
        environ = {"HTTP_AUTHORIZATION": "Basic invalid-base64!?"}
        self.assertFalse(wsgi.authorized(environ))


if __name__ == "__main__":
    unittest.main()
