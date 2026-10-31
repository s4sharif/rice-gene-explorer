"""Extract a data-frame slot from XDR RDS without materializing assay matrices.

This deliberately narrow reader supports the serialization types needed before a
Seurat meta.data slot and rejects unsupported types rather than guessing.
"""

import argparse
import array
import csv
import gzip
import struct
import sys
from contextlib import ExitStack
from pathlib import Path


class MetadataFound(Exception):
    pass


class SlotReader:
    # Initialize the streaming reader with an explicit target slot name.
    def __init__(self, handle, target="meta.data"):
        self.handle = handle
        self.target = target
        self.references = []
        self.metadata = None
        self.s4_depth = 0

    # Read exactly the requested bytes and fail on a truncated source object.
    def read(self, size):
        value = self.handle.read(size)
        if len(value) != size:
            raise ValueError("Truncated RDS stream.")
        return value

    # Decode one big-endian signed integer from the XDR stream.
    def integer(self):
        return struct.unpack(">i", self.read(4))[0]

    # Skip a numeric payload in small chunks to keep memory use bounded.
    def skip(self, size):
        while size:
            take = min(size, 1024 * 1024)
            self.read(take)
            size -= take

    # Decode an ordinary or long-vector length without accepting negative sizes.
    def length(self):
        size = self.integer()
        if size == -1:
            size = (self.integer() << 32) + self.integer()
        if size < 0:
            raise ValueError("Invalid R vector length.")
        return size

    # Read one supported serialized value while discarding unneeded assay payloads.
    def item(self, capture=False):
        flags = self.integer() & 0xFFFFFFFF
        kind = flags & 255
        has_attr = bool(flags & 512)
        has_tag = bool(flags & 1024)
        attrs = None
        value = None
        if kind == 255:
            reference = flags >> 8 or self.integer()
            return self.references[reference - 1]
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
            is_target = name == self.target and self.s4_depth == 1
            car = self.item(capture or is_target)
            if is_target and not capture:
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
            if capture:
                if kind not in {10, 13, 14}:
                    raise ValueError("Unsupported captured numeric type.")
                vector = array.array("d" if kind == 14 else "i")
                vector.frombytes(self.read(size * width))
                if sys.byteorder == "little":
                    vector.byteswap()
                value = vector.tolist()
            else:
                self.skip(size * width)
        elif kind in {16, 19, 20}:
            size = self.length()
            if capture:
                value = [self.item(True) for _ in range(size)]
            else:
                for _ in range(size):
                    self.item(False)
        elif kind == 25:
            self.s4_depth += 1
        elif kind == 238:
            info = self.item(True)
            state = self.item(capture)
            attrs = self.item(capture)
            if capture:
                raise ValueError(
                    f"Captured ALTREP requires explicit conversion: {info} {state}"
                )
            return None
        else:
            raise ValueError(
                f"Unsupported R serialization type {kind}; use R to extract metadata."
            )
        if has_attr:
            attrs = self.item(capture)
        if kind == 25:
            self.s4_depth -= 1
        return {"value": value, "attrs": attrs} if capture else None

    # Validate the RDS header and stop as soon as the requested data frame is found.
    def extract(self):
        if self.read(2) != b"X\n":
            raise ValueError("Only XDR RDS streams are supported.")
        version = self.integer()
        self.integer()
        self.integer()
        if version == 3:
            self.read(self.length())
        elif version != 2:
            raise ValueError("Unsupported R serialization version.")
        try:
            self.item()
        except MetadataFound:
            return self.metadata
        raise ValueError("Requested metadata slot was not found.")


# Convert a captured R pairlist into an attribute-name mapping.
def attributes(node):
    result = {}
    while node:
        name, value, node = node["value"]
        result[name] = value
    return result


# Decode a captured string or factor column without changing its labels.
def column_values(node):
    attrs = attributes(node["attrs"])
    if "levels" in attrs:
        levels = [x["value"] for x in attrs["levels"]["value"]]
        return [levels[x - 1] if x > 0 else "" for x in node["value"]]
    values = node["value"]
    return [x["value"] if isinstance(x, dict) else x for x in values]


# Export the original cell barcodes and metadata columns as a plain TSV file.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--gzip-layers", type=int, choices=[1, 2], default=1)
    args = parser.parse_args()
    with ExitStack() as stack:
        handle = stack.enter_context(args.source.open("rb"))
        for _ in range(args.gzip_layers):
            handle = stack.enter_context(gzip.GzipFile(fileobj=handle))
        node = SlotReader(handle).extract()
    attrs = attributes(node["attrs"])
    names = column_values(attrs["names"])
    rownames = column_values(attrs["row.names"])
    columns = [column_values(col) for col in node["value"]]
    if len(rownames) != len(set(rownames)) or any(
        len(c) != len(rownames) for c in columns
    ):
        raise ValueError("Metadata dimensions or barcode uniqueness failed.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["barcode", *names])
        for index, barcode in enumerate(rownames):
            writer.writerow([barcode, *[c[index] for c in columns]])
    print(f"Exported {len(rownames)} unique annotated cells; columns: {names}")


if __name__ == "__main__":
    main()
