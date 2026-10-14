#!/usr/bin/env python3
"""
Single source of truth for the rib layout / rib thickness variant of the Project 03 geometry pipeline.

Environment (read by the geometry, imprint and mesh scripts through configuration()):
  RIB_SPACING   uniform spacing in mm (default 200)
  RIB_STATIONS  comma-separated rib centre stations in mm, overriding RIB_SPACING; must start at 0, end at the
                span and keep at least 3 rib thicknesses between rib centres; requires LAYOUT_ID
  LAYOUT_ID     name of a custom layout: 'L' followed by letters/digits (e.g. L012)
  RIB_T         rib thickness in mm (default 1.5)

Tag appended to every output name: '' for the default (200 mm, 1.5 mm); '_rs<N>' for uniform N mm;
'_<LAYOUT_ID>' for a custom layout; plus '_rt<T>' (p for the decimal point) if RIB_T differs from 1.5.
"""
import os
import re

DEFAULT_SPACING = 200.0
DEFAULT_RIB_T = 1.5
MIN_GAP_FACTOR = 3.0


def _fmt(v):
    return ("%g" % v).replace(".", "p")


def configuration(span):
    """Return (rib_spacing, rib_stations, rib_t, tag). rib_spacing is 0.0 for a custom layout."""
    rib_t = float(os.environ.get("RIB_T", str(DEFAULT_RIB_T)))
    if not 0.1 <= rib_t <= 5.0:
        raise SystemExit("RIB_T=%g mm is outside the supported range 0.1-5 mm" % rib_t)
    env_st = os.environ.get("RIB_STATIONS", "").strip()
    if env_st:
        lid = os.environ.get("LAYOUT_ID", "")
        if not re.fullmatch(r"L[0-9A-Za-z]*", lid):
            raise SystemExit("RIB_STATIONS requires LAYOUT_ID = 'L' followed by letters/digits (e.g. L012)")
        stations = [float(v) for v in env_st.split(",")]
        if abs(stations[0]) > 1e-9 or abs(stations[-1] - span) > 1e-9:
            raise SystemExit("RIB_STATIONS must start at 0 and end at %g" % span)
        for a, b in zip(stations, stations[1:]):
            if b - a < MIN_GAP_FACTOR * rib_t:
                raise SystemExit("RIB_STATIONS: gap %g -> %g is below %g x RIB_T" % (a, b, MIN_GAP_FACTOR))
        spacing, tag = 0.0, "_" + lid
    else:
        spacing = float(os.environ.get("RIB_SPACING", str(DEFAULT_SPACING)))
        if abs(span / spacing - round(span / spacing)) > 1e-9:
            raise SystemExit("RIB_SPACING=%g must divide the %g mm span" % (spacing, span))
        stations = [i * spacing for i in range(int(round(span / spacing)) + 1)]
        tag = "" if spacing == DEFAULT_SPACING else "_rs%g" % spacing
    if rib_t != DEFAULT_RIB_T:
        tag += "_rt" + _fmt(rib_t)
    return spacing, stations, rib_t, tag
