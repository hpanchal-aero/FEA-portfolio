#!/usr/bin/env python3
"""Dump all non-empty cells of .xlsx files using only the standard library (zipfile + ElementTree). Read-only.

usage: dump_xlsx.py <file.xlsx> [more.xlsx ...]
Prints, per sheet, one line per row: 'row N: A=value | B=value | ...' (first 120 rows). Cached formula values are shown.
"""
import sys
import xml.etree.ElementTree as ET
import zipfile

M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def cell_value(c, shared):
    t = c.get("t")
    if t == "inlineStr":
        return "".join(x.text or "" for x in c.iter(M + "t"))
    v = c.find(M + "v")
    if v is None or v.text is None:
        return None
    if t == "s":
        return shared[int(v.text)]
    if t == "b":
        return "TRUE" if v.text == "1" else "FALSE"
    return v.text


for fpath in sys.argv[1:]:
    z = zipfile.ZipFile(fpath)
    names = z.namelist()
    shared = []
    if "xl/sharedStrings.xml" in names:
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(M + "si"):
            shared.append("".join(x.text or "" for x in si.iter(M + "t")))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    relmap = {r.get("Id"): r.get("Target") for r in rels}
    print("=" * 100)
    print(fpath)
    for sh in wb.find(M + "sheets"):
        target = relmap[sh.get(R + "id")].lstrip("/")
        sheet_path = target if target.startswith("xl/") else "xl/" + target
        ws = ET.fromstring(z.read(sheet_path))
        print(f"\n--- sheet '{sh.get('name')}'  ({sheet_path})")
        n = 0
        for row in ws.iter(M + "row"):
            cells = []
            for c in row.findall(M + "c"):
                val = cell_value(c, shared)
                if val not in (None, ""):
                    cells.append(f"{c.get('r')}={val}")
            if cells:
                n += 1
                if n <= 120:
                    print(f"row {row.get('r')}: " + " | ".join(cells))
        if n > 120:
            print(f"... {n - 120} more non-empty rows not shown")
