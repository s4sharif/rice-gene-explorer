"""Bound the free review service's worker and thread count."""

import os

bind = "0.0.0.0:" + os.environ.get("PORT", "8000")
workers = 1
worker_class = "gthread"
threads = 2
timeout = 60
keepalive = 5
max_requests = 1000
max_requests_jitter = 50
accesslog = "-"
errorlog = "-"

control_socket_disable = True
