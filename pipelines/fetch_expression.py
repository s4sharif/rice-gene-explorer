"""Fetch the frozen public inputs for the baseline leaf expression import."""

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DESTINATION = ROOT / "data/raw/expression"


# Calculate a file checksum with bounded memory use.
def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


# Fetch each required source and reject changes to the recorded input snapshot.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--study", choices=["robertson", "wang", "li", "li_author"], default="robertson"
    )
    study = parser.parse_args().study
    suffix = "" if study == "robertson" else "_" + study
    destination = ROOT / "data/raw" / ("expression" if study == "robertson" else study)
    destination.mkdir(parents=True, exist_ok=True)
    if study == "robertson":
        (destination / "selection.json").write_bytes(
            (ROOT / "data/expression_selection.json").read_bytes()
        )
    elif study == "li":
        (destination / "selection.json").write_bytes(
            (ROOT / "data/expression_li_selection.json").read_bytes()
        )
    manifest = json.loads((ROOT / f"data/expression{suffix}_inputs.json").read_text())
    for source in manifest["sources"]:
        path = destination / source["filename"]
        if path.exists() and sha256_file(path) == source["sha256"]:
            print(f"Verified {path.name}", flush=True)
            continue
        if "url" not in source:
            raise ValueError(
                f"Place the supplied publisher supplement at {path} and rerun."
            )
        temporary = path.with_suffix(path.suffix + ".part")
        with (
            urllib.request.urlopen(source["url"], timeout=60) as incoming,
            temporary.open("wb") as output,
        ):
            while block := incoming.read(1024 * 1024):
                output.write(block)
        if sha256_file(temporary) != source["sha256"]:
            raise ValueError(f"Source changed: {path.name}; review before importing.")
        temporary.replace(path)
        print(f"Downloaded and verified {path.name}", flush=True)


if __name__ == "__main__":
    main()
