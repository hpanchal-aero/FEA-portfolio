#!/usr/bin/python3
"""Small readers for CalculiX output and for the mesh include files of Project 05."""
import re

import numpy as np

_HEAD = re.compile(r"^\s*(.+?)\s+for set\s+(\S+)\s+and time\s+(\S+)\s*$")
_PARAMS = re.compile(r"L=(\S+)\s+W=(\S+)\s+h=(\S+)\s+m,\s+nx=(\d+)\s+ny=(\d+)\s+nz=(\d+)")


def read_params(mesh_inp):
    """Panel dimensions and divisions from the first comment line of a mesh include file."""
    with open(mesh_inp) as f:
        first = f.readline()
    m = _PARAMS.search(first)
    if not m:
        raise ValueError("cannot read panel parameters from the first line of %s" % mesh_inp)
    L, W, h = (float(m.group(i)) for i in (1, 2, 3))
    nx, ny, nz = (int(m.group(i)) for i in (4, 5, 6))
    return {"L": L, "W": W, "h": h, "nx": nx, "ny": ny, "nz": nz}


def read_nodes(mesh_inp):
    """Return {node id: (x, y, z)} from the *NODE block of a mesh include file."""
    nodes = {}
    on = False
    with open(mesh_inp) as f:
        for line in f:
            s = line.strip()
            if s.startswith("**"):
                continue
            if s.startswith("*"):
                on = s.upper().startswith("*NODE")
                continue
            if on and s:
                p = s.split(",")
                nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
    return nodes


def read_elements(mesh_inp):
    """Return (element ids, node ids as an (E, 20) array) from the *ELEMENT block (20-node hexahedra)."""
    ids, conn, buf = [], [], []
    on = False
    with open(mesh_inp) as f:
        for line in f:
            s = line.strip()
            if s.startswith("**"):
                continue
            if s.startswith("*"):
                on = s.upper().startswith("*ELEMENT")
                continue
            if on and s:
                buf.extend(int(v) for v in s.split(",") if v.strip())
                if len(buf) == 21:
                    ids.append(buf[0])
                    conn.append(buf[1:])
                    buf = []
    if buf:
        raise ValueError("incomplete element record at the end of %s" % mesh_inp)
    return np.array(ids, dtype=np.int64), np.array(conn, dtype=np.int64)


def read_nset(mesh_inp, name):
    """Return the node ids of a *NSET block (in file order)."""
    ids = []
    on = False
    with open(mesh_inp) as f:
        for line in f:
            s = line.strip()
            if s.startswith("**"):
                continue
            if s.startswith("*"):
                u = s.upper().replace(" ", "")
                on = u.startswith("*NSET,NSET=") and u.split("NSET=")[1].split(",")[0] == name.upper()
                continue
            if on and s:
                ids.extend(int(v) for v in s.split(",") if v.strip())
    return ids


def parse_dat(dat_path):
    """Parse a CalculiX .dat file into blocks.

    Each block: {"title", "set", "time", "rows"}; rows is a list of float lists
    (for nodal output the first number of a row is the node id).
    """
    blocks = []
    cur = None
    with open(dat_path) as f:
        for line in f:
            m = _HEAD.match(line)
            if m:
                cur = {"title": m.group(1).strip(), "set": m.group(2),
                       "time": float(m.group(3).replace("D", "E")), "rows": []}
                blocks.append(cur)
                continue
            s = line.strip()
            if not s or cur is None:
                continue
            try:
                vals = [float(v) for v in s.split()]
            except ValueError:
                continue
            cur["rows"].append(vals)
    return blocks
