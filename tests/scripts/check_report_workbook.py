"""Check workbook headers and Dictionary without extracting large worksheets."""

import json
import re
import sys
import zipfile
from xml.etree import ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def cell_value(cell):
    if cell.get("t") == "s":
        return int(cell.findtext(NS + "v"))
    return "".join(cell.itertext())


def inspect_workbook(path):
    with zipfile.ZipFile(path) as archive:
        sheets = ET.fromstring(archive.read("xl/workbook.xml")).find(NS + "sheets")
        names = [sheet.get("name") for sheet in sheets]
        assert names[:2] == ["Contents", "Dictionary"]
        assert sheets[1].get("state", "visible") == "visible"
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {item.get("Id"): item.get("Target") for item in relationships}
        headers = {}
        dictionary = []
        contents = []
        for sheet in sheets:
            relation = sheet.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
            target = targets[relation]
            member = target.lstrip("/") if target.startswith("/") else "xl/" + target
            with archive.open(member) as stream:
                tail = b""
                while chunk := stream.read(1024 * 1024):
                    block = tail + chunk
                    assert not re.search(rb"<(?:\w+:)?(?:autoFilter|pane|cols)\b", block), member
                    tail = block[-128:]
            with archive.open(member) as stream:
                for _, element in ET.iterparse(stream, events=("end",)):
                    if element.tag != NS + "row":
                        continue
                    values = [cell_value(cell) for cell in element.findall(NS + "c")]
                    if element.get("r") == "1":
                        headers[sheet.get("name")] = values
                        if sheet.get("name") not in ("Contents", "Dictionary"):
                            break
                    elif sheet.get("name") == "Dictionary":
                        dictionary.append(values)
                    elif sheet.get("name") == "Contents":
                        contents.append(values)
                    element.clear()
        wanted = {value for row in list(headers.values()) + dictionary + contents for value in row if isinstance(value, int)}
        strings = {}
        with archive.open("xl/sharedStrings.xml") as stream:
            position = 0
            for _, element in ET.iterparse(stream, events=("end",)):
                if element.tag == NS + "si":
                    if position in wanted:
                        strings[position] = "".join(element.itertext())
                    position += 1
                    element.clear()
        resolve = lambda row: [strings[value] if isinstance(value, int) else value for value in row]
        headers = {name: resolve(values) for name, values in headers.items()}
        records = [dict(zip(headers["Dictionary"], resolve(row))) for row in dictionary]
        content_records = [dict(zip(headers["Contents"], resolve(row))) for row in contents]
        for name, columns in headers.items():
            assert set(columns) == {row["column"] for row in records if row["worksheet"] == name}, name
        assert not any(name.startswith("xl/tables/table") and name.endswith(".xml") for name in archive.namelist())
        return {"sheets": names, "dictionary": records, "contents": content_records}


if __name__ == "__main__":
    result = inspect_workbook(sys.argv[1])
    if "--json" in sys.argv[2:]:
        print(json.dumps(result))
    else:
        unknown = sorted({row["column"] for row in result["dictionary"] if row["definition"].startswith("Description not supplied")})
        print(f"Workbook passed: {len(result['sheets'])} sheets; visible Dictionary second; every column indexed; no styled tables, filters, panes or custom column widths.")
        if unknown:
            print("Fields without source definitions: " + ", ".join(unknown))
