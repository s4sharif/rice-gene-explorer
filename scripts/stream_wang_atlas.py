"""Extract only RNA counts and cell metadata from the deposited Wang Seurat RDS.

The supported XDR serialization is read sequentially; large sparse vectors go to
big-endian disk files rather than Python lists. Unsupported encodings fail closed.
"""

import argparse
import csv
import gzip
import json
from pathlib import Path

from scripts.extract_rds_metadata import (
    MetadataFound,
    SlotReader,
    attributes,
    column_values,
)


class AtlasReader(SlotReader):
    # Track the RNA assay candidate and stream its sparse arrays into a dedicated output folder.
    def __init__(self, handle, destination, assay_index=0):
        super().__init__(handle)
        self.destination = Path(destination)
        self.destination.mkdir(parents=True, exist_ok=True)
        self.path = []
        self.assay_index = None
        self.target_index = assay_index
        self.assay_names = None
        self.matrix = None
        self.in_counts = False

    # Parse supported values while capturing only the selected original-count matrix and metadata.
    def item(self, capture=False):
        flags = self.integer() & 0xFFFFFFFF
        kind = flags & 255
        has_attr, has_tag = bool(flags & 512), bool(flags & 1024)
        attrs = None
        value = None
        if kind == 255:
            return self.references[(flags >> 8 or self.integer()) - 1]
        if kind in {0, 241, 242, 251, 253, 254}:
            return None
        if kind == 1:
            node = self.item(True)
            self.references.append(node)
            return node
        if kind in {2, 6, 3, 5, 17}:
            if has_attr:
                attrs = self.item(capture)
            tag = self.item(True) if has_tag else None
            name = tag.get("value") if tag else None
            is_metadata = name == "meta.data" and self.s4_depth == 1
            is_counts = (
                name == "counts"
                and self.path == ["assays"]
                and self.assay_index == self.target_index
            )
            if name:
                self.path.append(name)
            previous_mode = self.in_counts
            if is_counts:
                if self.matrix is not None:
                    raise ValueError("Multiple count slots in the selected assay.")
                self.in_counts = True
            car = self.item(capture or is_metadata or is_counts)
            self.in_counts = previous_mode
            if name:
                self.path.pop()
            if is_counts:
                self.matrix = car
                print("Streamed original count matrix", flush=True)
            if is_metadata and not capture:
                self.metadata = car
                raise MetadataFound
            cdr = self.item(capture)
            return {"value": (name, car, cdr), "attrs": attrs} if capture else None
        if kind in {7, 8, 9}:
            size = self.integer()
            value = self.read(size).decode("utf-8") if size >= 0 else None
        elif kind in {10, 13, 14, 15, 24}:
            size = self.length()
            width = {10: 4, 13: 4, 14: 8, 15: 16, 24: 1}[kind]
            if self.in_counts and self.path[-1:] in (["i"], ["p"], ["x"]):
                slot = self.path[-1]
                if (slot in {"i", "p"} and kind != 13) or (slot == "x" and kind != 14):
                    raise ValueError(
                        "Expected dgCMatrix integer indices and double counts."
                    )
                target = self.destination / (slot + ".bin")
                remaining = size * width
                with target.open("wb") as output:
                    while remaining:
                        chunk = min(remaining, 1024 * 1024)
                        output.write(self.read(chunk))
                        remaining -= chunk
                value = {
                    "file": target.name,
                    "length": size,
                    "dtype": ">i4" if kind == 13 else ">f8",
                }
            elif capture:
                import array
                import sys

                if kind not in {10, 13, 14} or size > 5000000:
                    raise ValueError(
                        f"Unexpected captured numeric vector: {self.path}, type={kind}, length={size}"
                    )
                vector = array.array("d" if kind == 14 else "i")
                vector.frombytes(self.read(size * width))
                if sys.byteorder == "little":
                    vector.byteswap()
                value = vector.tolist()
            else:
                self.skip(size * width)
        elif kind in {16, 19, 20}:
            size = self.length()
            if self.path == ["assays"] and not capture:
                for index in range(size):
                    self.assay_index = index
                    self.item(False)
                self.assay_index = None
            elif capture:
                value = [self.item(True) for _ in range(size)]
            else:
                for _ in range(size):
                    self.item(False)
        elif kind in {22, 23}:
            if capture:
                raise ValueError(
                    "External pointers and weak references cannot be expression data."
                )
            self.references.append(None)
            if kind == 22:
                self.item(False)
                self.item(False)
        elif kind == 4:
            if capture:
                raise ValueError("Environments cannot be expression data.")
            self.integer()
            self.references.append(None)
            for _ in range(4):
                self.item(False)
            return None
        elif kind == 25:
            self.s4_depth += 1
        elif kind == 238:
            info = self.item(True)
            state = self.item(capture)
            attrs = self.item(capture)
            if capture:
                raise ValueError(f"Unsupported captured ALTREP: {info} {state}")
            return None
        else:
            raise ValueError(f"Unsupported R serialization type {kind}.")
        if has_attr:
            attrs = self.item(
                capture or (self.path == ["assays"] and kind in {16, 19, 20})
            )
        if self.path == ["assays"] and attrs and kind in {16, 19, 20}:
            self.assay_names = column_values(attributes(attrs)["names"])
        if kind == 25:
            self.s4_depth -= 1
        return {"value": value, "attrs": attrs} if capture else None

    # Require an explicitly named RNA assay and a complete sparse matrix before exporting it.
    def export(self):
        metadata = self.extract()
        if not self.assay_names or self.assay_names[self.target_index] != "RNA":
            raise ValueError(f"Selected assay is not verified RNA: {self.assay_names}")
        if self.matrix is None:
            raise ValueError(
                "RNA counts were not found; normalized data are never substituted."
            )
        slots = attributes(self.matrix["attrs"])
        if column_values(slots["class"]) != ["dgCMatrix"]:
            raise ValueError("Only the original dgCMatrix RNA count slot is supported.")
        dim = slots["Dim"]["value"]
        features, cells = [column_values(n) for n in slots["Dimnames"]["value"]]
        if (
            dim != [len(features), len(cells)]
            or len(set(features)) != len(features)
            or len(set(cells)) != len(cells)
        ):
            raise ValueError("RNA matrix names and dimensions do not agree.")
        specification = {
            "shape": dim,
            "assay": "RNA",
            "slot": "counts",
            "assay_names": self.assay_names,
            "arrays": {name: slots[name]["value"] for name in ["i", "p", "x"]},
        }
        (self.destination / "matrix.json").write_text(
            json.dumps(specification, indent=2) + "\n"
        )
        with (self.destination / "features.tsv").open("w") as output:
            writer = csv.writer(output, delimiter="\t")
            writer.writerow(["feature_id"])
            writer.writerows([[name] for name in features])
        with (self.destination / "barcodes.tsv").open("w") as output:
            writer = csv.writer(output, delimiter="\t")
            writer.writerow(["barcode"])
            writer.writerows([[name] for name in cells])
        attrs = attributes(metadata["attrs"])
        names, barcodes = column_values(attrs["names"]), column_values(
            attrs["row.names"]
        )
        columns = [column_values(column) for column in metadata["value"]]
        if (
            set(barcodes) != set(cells)
            or len(barcodes) != len(cells)
            or any(len(c) != len(cells) for c in columns)
        ):
            raise ValueError(
                "Cell metadata and the RNA count columns must match exactly."
            )
        with (self.destination / "cell_metadata.tsv").open("w") as output:
            writer = csv.writer(output, delimiter="\t")
            writer.writerow(["barcode", *names])
            for index, barcode in enumerate(barcodes):
                writer.writerow([barcode, *[c[index] for c in columns]])
        print(
            json.dumps(
                {"shape": dim, "metadata_columns": names, "assays": self.assay_names}
            ),
            flush=True,
        )
        return specification


# Extract the verified RNA count slot and original cell metadata using bounded vector buffers.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--assay-index", type=int, default=0)
    args = parser.parse_args()
    with gzip.open(args.source, "rb") as handle:
        AtlasReader(handle, args.destination, args.assay_index).export()


if __name__ == "__main__":
    main()
