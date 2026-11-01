"""Locate local or restored hosting databases without moving source metadata."""

import os
from pathlib import Path

DATA_DIR = Path(
    os.environ.get("RICELENS_DATA_DIR", str(Path(__file__).resolve().parent / "data"))
).resolve()
