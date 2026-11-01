"""Read literal worksheet values from published XLSX supplements without Excel."""

import posixpath
import xml.etree.ElementTree as ET
from zipfile import ZipFile

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RELATIONSHIP = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


# Read cells by Excel column letter while rejecting formulas that need recalculation.
def worksheet_rows(path, sheet_name):
    with ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings = [
                "".join(item.itertext())
                for item in ET.fromstring(archive.read("xl/sharedStrings.xml")).findall(
                    "m:si", NS
                )
            ]
        relationships = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        }
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = {
            item.attrib["name"]: item.attrib[RELATIONSHIP]
            for item in workbook.findall("m:sheets/m:sheet", NS)
        }
        target = relationships[sheets[sheet_name]]
        target = (
            target.lstrip("/")
            if target.startswith("/")
            else posixpath.normpath("xl/" + target)
        )
        rows = []
        root = ET.fromstring(archive.read(target))
        for row in root.findall("m:sheetData/m:row", NS):
            result = {}
            for cell in row.findall("m:c", NS):
                if cell.find("m:f", NS) is not None:
                    raise ValueError("Source annotations must contain literal values.")
                value = cell.find("m:v", NS)
                inline = cell.find("m:is", NS)
                text = value.text if value is not None else ""
                if cell.attrib.get("t") == "s":
                    text = strings[int(text)]
                elif inline is not None:
                    text = "".join(inline.itertext())
                result[cell.attrib["r"].rstrip("0123456789")] = text or ""
            rows.append(result)
        return rows
