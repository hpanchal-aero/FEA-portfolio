#!/usr/bin/env python3
"""Minimal fixed-format Nastran BDF reader for the uCRM-9 files (read-only).

Handles small-field (8-char) and large-field (16-char, name ends with '*') cards
with continuation lines. Free-format (comma) cards are skipped. Not a general parser.
"""
import re

_EXP = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$")


def num(s):
    """Nastran real: handles blanks and the short exponent form 1.234-3."""
    s = s.strip()
    if not s:
        return 0.0
    try:
        return float(s)
    except ValueError:
        m = _EXP.match(s)
        if m:
            return float(m.group(1) + "e" + m.group(2))
        raise ValueError(f"cannot parse numeric field {s!r}")


def ival(s):
    s = s.strip()
    return int(float(s)) if s else 0


def _chunks(line, width):
    body = line[8:72].ljust(64)
    return [body[i:i + width] for i in range(0, 64, width)]


def read_cards(path):
    """Return a list of (NAME, [field strings]) for every fixed-format card."""
    cards = []
    cur_name, cur_fields = None, None
    with open(path, errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip() or line.startswith("$"):
                continue
            c0 = line[0]
            if c0.isalpha():
                if cur_name is not None:
                    cards.append((cur_name, cur_fields))
                head = line[:9]
                if "," in head:
                    cur_name, cur_fields = None, None
                    continue
                name = line[:8].strip()
                large = name.endswith("*")
                cur_name = name.rstrip("*").upper()
                cur_fields = _chunks(line, 16 if large else 8)
            elif c0 in "*+ " and cur_name is not None:
                cur_fields.extend(_chunks(line, 16 if c0 == "*" else 8))
        if cur_name is not None:
            cards.append((cur_name, cur_fields))
    return cards
