#!/usr/bin/env python3
"""Reference cruise condition of the uCRM-9 from the dataset's design-specification workbook (ISA atmosphere).

usage: flight_conditions.py [path/to/uCRM-9_Design_Specifications.xlsx]

Reads span, MAC, reference area, altitude, Mach, CL, fuel load, OEW, payload and fuel capacity from the workbook (not copied into the
repository), then computes ISA static conditions, Sutherland viscosity, freestream velocity, dynamic pressure, Re on the MAC,
half-model reference area and target lift, and compares the full-aircraft lift at CL with the cruise weight (consistency check).
Standard library only.
"""
import math
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
default = Path(__file__).resolve().parent.parent / "data/ucrm9/uCRM9/Data/uCRM-9_Design_Specifications.xlsx"
xlsx = Path(sys.argv[1]) if len(sys.argv) > 1 else default


def read_labels(path):
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(M + "si"):
            shared.append("".join(x.text or "" for x in si.iter(M + "t")))
    ws = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    out = {}
    for row in ws.iter(M + "row"):
        vals = {}
        for c in row.findall(M + "c"):
            v = c.find(M + "v")
            if v is None or v.text is None:
                continue
            col = "".join(ch for ch in c.get("r") if ch.isalpha())
            vals[col] = shared[int(v.text)] if c.get("t") == "s" else v.text
        if vals.get("A") and vals.get("B"):
            out[vals["A"].strip()] = vals["B"]
    return out


d = read_labels(xlsx)
need = ["Span", "MAC", "Wimpress reference area", "Altitude", "Mach", "CL", "Fuel load", "Fuel capacity (both wings)",
        "Design payload", "Operational empty weight (OEW)"]
missing = [k for k in need if k not in d]
if missing:
    sys.exit(f"ABORT: labels not found in {xlsx}: {missing}")
f = {k: float(d[k]) for k in need}

g0, R, gam = 9.80665, 287.053, 1.4
h = f["Altitude"] * 0.3048                                   # m
if h < 11000.0:
    sys.exit("ABORT: script implements the isothermal ISA layer (11-20 km) only")
T = 216.65
p = 22632.06 * math.exp(-g0 * (h - 11000.0) / (R * T))
rho = p / (R * T)
a = math.sqrt(gam * R * T)
U = f["Mach"] * a
q = 0.5 * rho * U * U
mu = 1.458e-6 * T ** 1.5 / (T + 110.4)                       # Sutherland
nu = mu / rho
Re = rho * U * f["MAC"] / mu
S_full, S_half = f["Wimpress reference area"], f["Wimpress reference area"] / 2.0
L_full = f["CL"] * q * S_full
fuel = f["Fuel load"] / 100.0 * f["Fuel capacity (both wings)"]
mass = f["Operational empty weight (OEW)"] + f["Design payload"] + fuel
W = mass * g0

print(f"workbook: {xlsx.name}")
print(f"input: altitude {f['Altitude']:.0f} ft ({h:.1f} m), Mach {f['Mach']}, CL {f['CL']}, fuel load {f['Fuel load']:.0f} %,"
      f" MAC {f['MAC']} m, span {f['Span']} m, Sref {S_full} m2")
print(f"ISA (isothermal layer): T {T:.2f} K, p {p:.1f} Pa, rho {rho:.5f} kg/m3, speed of sound {a:.2f} m/s")
print(f"freestream: U {U:.2f} m/s, dynamic pressure q {q:.1f} Pa, mu {mu:.4e} Pa s, nu {nu:.4e} m2/s")
print(f"Reynolds number on the MAC: {Re:.4e}")
print(f"half model: Sref/2 {S_half:.3f} m2, semispan {f['Span'] / 2:.3f} m, lift target at CL {f['CL']}: {f['CL'] * q * S_half:.1f} N")
print(f"consistency: full-aircraft lift CL*q*Sref = {L_full:.1f} N; cruise weight {mass:.1f} kg x g = {W:.1f} N "
      f"(OEW + payload + {f['Fuel load']:.0f} % fuel = {fuel:.1f} kg); lift/weight = {L_full / W:.4f}")
