#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""02_run_modal_virtual_elastic.py

CalculiX modal (eigenfrequency) analysis of the combined system:

    * elastic dam-break gate   (real material:  E=1e6 Pa, nu=0, rho=2500 kg/m^3)
    * fluid domain             (virtual elastic body: per-cell layered stiffness
                                 E_i = 0.1 * (min_vol/vol_i)^2.5, nu=0, rho=0)

The fluid domain is modelled as a *virtual elastic body* whose stiffness is
layered by cell volume (small cells near the structure get a larger Young's
modulus, large cells far away get a smaller one), exactly mirroring the APDL
reference ``ModeCalculation.apdl``::

    E_i = 0.1 * (min_vol / vol_i) ** 2.5

The gate and the fluid domain are tied at the FSI interface with a CalculiX
``*TIE`` constraint (the CalculiX equivalent of the APDL ``CEINTF`` command),
so the (massless) fluid follows the gate and only contributes stiffness.

The first ``--n-modes`` natural frequencies are extracted and compared against
paper table 3-2, right column ("gate + virtual-elastic-body system"):
6.11 / 35.52 / 63.53 / 90.42 / 159.44 Hz.  The left column (bare gate)
6.02 / 34.98 / 62.51 / 89.07 / 156.83 Hz is printed for reference.

Only the Python standard library is used (no third-party dependencies).

Usage:
    python3 02_run_modal_virtual_elastic.py \
        --fluid-dir /root/Workspace/precice_damfailure_new/fluid-openfoam \
        --out-dir   /root/Workspace/precice_damfailure_modal_ve

    # quick coarse validation (much smaller fluid mesh)
    python3 02_run_modal_virtual_elastic.py --fluid-dir ... --coarse-fluid
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys

# --------------------------------------------------------------------------- #
#  Geometry / material defaults (from 01_generate_mesh.py + CaseDescription.md)  #
# --------------------------------------------------------------------------- #

LX = 0.584      # tank full length in x [m]
LY = 0.365      # tank full length in y [m]
LZ = 0.012      # tank thickness in z [m] (quasi-2D)

FSI_X0 = 0.292  # gate left edge x
FSI_X1 = 0.304  # gate right edge x (12 mm wide)
FSI_Y = 0.080   # gate height (80 mm)

# x / y / z split coordinates (block corners), identical to 01_generate_mesh.py
XS = [0.0, 0.146, FSI_X0, FSI_X1, LX]
YS = [0.0, FSI_Y, 0.292, LY]
ZS = [0.0, LZ]

# Block table: (xi, yi, nx_key, ny_key, gx, gy).  gz = 1 always.
BLOCK_TABLE = [
    (0, 0, "nx_a", "ny_l", 5.0, 15.0),   # block 1
    (1, 0, "nx_b", "ny_l", 0.12, 15.0),  # block 2
    (3, 0, "nx_d", "ny_l", 8.0, 15.0),   # block 3
    (0, 1, "nx_a", "ny_m", 5.0, 10.0),   # block 4
    (1, 1, "nx_b", "ny_m", 0.12, 10.0),  # block 5
    (2, 1, "nx_c", "ny_m", 1.0, 10.0),   # block 6
    (3, 1, "nx_d", "ny_m", 8.0, 10.0),   # block 7
    (0, 2, "nx_a", "ny_u", 5.0, 0.15),   # block 8
    (1, 2, "nx_b", "ny_u", 0.12, 0.15),  # block 9
    (2, 2, "nx_c", "ny_u", 1.0, 0.15),   # block 10
    (3, 2, "nx_d", "ny_u", 8.0, 0.15),   # block 11
]

# Fallback counts of the actual generated fluid mesh (res ~115).  Used only if
# the blockMeshDict cannot be parsed.  Sum of nx*ny over blocks == 15632,
# times Nz == 187584 cells, matching checkMesh on the shipped case.
DEFAULT_COUNTS = {
    "nx_a": 29, "nx_b": 58, "nx_c": 12, "nx_d": 83,
    "ny_l": 32, "ny_m": 42, "ny_u": 14, "Nz": 12,
}

# Gate mesh resolution (matches 01_generate_mesh.build_solid_mesh -> 868/540).
GATE_NX, GATE_NY, GATE_NZ = 6, 30, 3

GATE_E_DEFAULT = 1.0e6
GATE_NU_DEFAULT = 0.0
GATE_RHO_DEFAULT = 2500.0

# Virtual-elastic-body layered stiffness (ModeCalculation.apdl).
FLUID_BASE_E = 0.1      # E_i = 0.1 * (min_vol/vol_i)^2.5
FLUID_EXP = 2.5         # power exponent
FLUID_NU = 0.0
FLUID_RHO = 0.0         # massless virtual body

# Paper table 3-2 (Hz).
PAPER_BARE = [6.02, 34.98, 62.51, 89.07, 156.83]
PAPER_VE = [6.11, 35.52, 63.53, 90.42, 159.44]

CCX_EXE = "ccx"
JOBNAME = "gate_fluid_modal"

FLUID_DIR_DEFAULT = "/root/Workspace/precice_damfailure_new/fluid-openfoam"
OUT_DIR_DEFAULT = "/root/Workspace/precice_damfailure_modal_ve"


# --------------------------------------------------------------------------- #
#  Small helpers                                                               #
# --------------------------------------------------------------------------- #

def fmt(x):
    """Compact float formatting for CalculiX ASCII files."""
    if isinstance(x, int):
        return str(x)
    return format(float(x), ".10g")


def _to_float(token):
    return float(token.replace("D", "E").replace("d", "e"))


def format_f20(x):
    """Fortran F20.10 field (coordinates in the APDL mode files)."""
    return format(float(x), "20.10f")


def format_e20(x):
    """Fortran E20.10 field (0.ddddddddddE+ee) for the APDL mode files.

    Python's ``%.10E`` normalises to ``6.0181352040E+00``, while the reference
    files use ``0.6018135204E+01``.  Returns the 20-character right-justified
    field including the leading sign/spaces.
    """
    x = float(x)
    if x == 0.0:
        return format("0.0000000000E+00", ">20")
    e = math.floor(math.log10(abs(x))) + 1
    m = x / (10.0 ** e)
    if abs(m) >= 1.0 - 1e-12:  # mantissa rounding to 1.0 must carry to exponent
        m /= 10.0
        e += 1
    return format("%.10fE%+03d" % (m, e), ">20")


def tet_signed_volume(a, b, c, d):
    """Signed volume/6 of a tetrahedron (a,b,c,d)."""
    bx, by, bz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    cx, cy, cz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    dx, dy, dz = d[0] - a[0], d[1] - a[1], d[2] - a[2]
    cross_x = by * cz - bz * cy
    cross_y = bz * cx - bx * cz
    cross_z = bx * cy - by * cx
    return (cross_x * dx + cross_y * dy + cross_z * dz) / 6.0


def hex_volume(nodes):
    """Exact volume of a hexahedron in C3D8 order (6-tet decomposition)."""
    v = nodes
    tets = [
        (0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
        (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6),
    ]
    total = 0.0
    for a, b, c, d in tets:
        total += tet_signed_volume(v[a], v[b], v[c], v[d])
    return total


def grading_positions(L, n, g):
    """Node positions along a length ``L`` divided into ``n`` cells.

    ``g`` is either a single blockMesh ``simpleGrading`` expansion ratio
    (float, last/first cell width) or a multi-segment tuple of
    (blockFraction, nDivFraction, expansionRatio) triples, e.g.
    ((0.5, 0.5, g1), (0.5, 0.5, g2)).

    Returns n+1 coordinates in [0, L].  Both forms reproduce blockMesh
    exactly:

    * single ratio: the cell widths form a geometric series with per-cell
      ratio q = g**(1/(n-1)), so node k sits at L * (1 - q**k)/(1 - q**n);
    * multi-segment: mirrors blockMesh ``lineDivide::lineDivide`` -- each
      segment gets ``label(nDivFraction*n + 0.5)`` cells (the segment with
      the largest nDivFraction absorbs the rounding remainder so the total
      is exactly ``n``), and within each segment the nodes follow a geometric
      series with that segment's own expansion ratio, concatenated end to end.
    """
    if n <= 1:
        return [0.0, L]
    if isinstance(g, tuple):
        if not g:
            return [L * i / n for i in range(n + 1)]
        # normalise fractions exactly like gradingDescriptors::normalise()
        bf = [float(s[0]) for s in g]
        nf = [float(s[1]) for s in g]
        er = [float(s[2]) for s in g]
        sum_bf = sum(bf)
        sum_nf = sum(nf)
        bf = [x / sum_bf for x in bf]
        nf = [x / sum_nf for x in nf]
        # per-segment cell counts (blockMesh lineDivide::lineDivide)
        nd = [int(x * n + 0.5) for x in nf]
        if sum(nd) != n:
            imax = max(range(len(nf)), key=lambda i: nf[i])
            nd[imax] += n - sum(nd)
        # normalised [0,1] positions (blockMesh divisions_)
        divs = [0.0] * (n + 1)
        sec_start = 0.0
        secn_start = 1
        for si in range(len(g)):
            secn_div = nd[si]
            secn_end = secn_start + secn_div
            if abs(er[si] - 1.0) < 1e-12:
                for i in range(secn_start, secn_end):
                    divs[i] = (sec_start
                               + bf[si] * (i - secn_start + 1) / secn_div)
            else:
                exp_fact = (er[si] ** (1.0 / (secn_div - 1))
                            if secn_div > 1 else 0.0)
                for i in range(secn_start, secn_end):
                    divs[i] = (sec_start
                               + bf[si]
                               * (1.0 - exp_fact ** (i - secn_start + 1))
                               / (1.0 - exp_fact ** secn_div))
            sec_start = divs[secn_end - 1]
            secn_start = secn_end
        divs[n] = 1.0
        return [L * d for d in divs]
    # single ratio (existing blockMesh reproduction)
    if abs(g - 1.0) < 1e-12:
        return [L * i / n for i in range(n + 1)]
    q = g ** (1.0 / (n - 1))
    qn = q ** n
    return [L * (1.0 - q ** i) / (1.0 - qn) for i in range(n + 1)]


# --------------------------------------------------------------------------- #
#  blockMeshDict parsing                                                       #
# --------------------------------------------------------------------------- #

def _strip_comments(text):
    out = []
    for line in text.splitlines():
        line = line.split("//", 1)[0]
        out.append(line)
    return "\n".join(out)


def _parse_list_entries(body):
    """Yield tuples of numeric tokens from a `( ... )` body (vertices)."""
    for m in re.finditer(r"\(\s*([^()]*?)\s*\)", body):
        toks = m.group(1).split()
        if toks:
            yield [float(t) for t in toks]


def _skip_ws(s, i):
    while i < len(s) and s[i].isspace():
        i += 1
    return i


def _match_paren_end(s, start):
    """Index one past the ``)`` matching the ``(`` at ``s[start]``."""
    depth = 0
    i = start
    n = len(s)
    while i < n:
        ch = s[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise ValueError("unbalanced parentheses in blockMeshDict")


def _split_top_level(s):
    """Split ``s`` on whitespace at parenthesis depth zero."""
    toks = []
    cur = []
    depth = 0
    for ch in s:
        if ch == "(":
            depth += 1
            cur.append(ch)
        elif ch == ")":
            depth -= 1
            cur.append(ch)
        elif ch.isspace():
            if depth == 0:
                if cur:
                    toks.append("".join(cur))
                    cur = []
            else:
                cur.append(ch)
        else:
            cur.append(ch)
    if cur:
        toks.append("".join(cur))
    return toks


def _parse_grading_token(tok):
    """Parse one simpleGrading token into a float or segment tuple.

    A plain expansion ratio returns a float; the multi-segment form
    ``((f1 n1 e1) (f2 n2 e2) ...)`` returns a tuple of
    (blockFraction, nDivFraction, expansionRatio) triples.
    """
    tok = tok.strip()
    if not tok.startswith("("):
        return float(tok)
    segs = []
    depth = 0
    start = None
    for i, ch in enumerate(tok):
        if ch == "(":
            depth += 1
            if depth == 2:
                start = i + 1
        elif ch == ")":
            if depth == 2 and start is not None:
                inner = tok[start:i].split()
                if len(inner) == 3:
                    segs.append((float(inner[0]), float(inner[1]),
                                 float(inner[2])))
                start = None
            depth -= 1
    if not segs:
        raise ValueError("bad simpleGrading token: %r" % tok)
    return tuple(segs)


def parse_block_mesh_dict(path):
    """Parse a blockMeshDict into a list of fluid block descriptors.

    Each block is returned as
        {xi, yi, nx, ny, nz, gx, gy}
    where (xi, yi) indexes into XS/YS and gx/gy are the simpleGrading ratios.
    gx/gy are either floats (single ratio) or tuples of (blockFraction,
    nDivFraction, expansionRatio) triples (multi-segment grading).
    Raises ValueError when the file is missing or malformed.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = _strip_comments(fh.read())

    vm = re.search(r"\bvertices\s*\(([\s\S]*?)\)\s*;", text)
    if not vm:
        raise ValueError("blockMeshDict: cannot find vertices section")
    verts = []
    for toks in _parse_list_entries(vm.group(1)):
        if len(toks) == 3:
            verts.append(tuple(toks))
    if not verts:
        raise ValueError("blockMeshDict: empty vertices section")

    bm = re.search(r"\bblocks\s*\(([\s\S]*?)\)\s*;", text)
    if not bm:
        raise ValueError("blockMeshDict: cannot find blocks section")
    body = bm.group(1)

    blocks = []
    # block line: hex (v0 v1 ... v7) (nx ny nz) simpleGrading (gx gy gz).
    # The grading group may contain nested parens (multi-segment grading),
    # so match each group with balanced-parenthesis scanning.
    for hm in re.finditer(r"\bhex\b", body, re.IGNORECASE):
        try:
            i = _skip_ws(body, hm.end())
            if i >= len(body) or body[i] != "(":
                continue
            j = _match_paren_end(body, i)
            vid = [int(t) for t in body[i + 1:j - 1].split()]
            i = _skip_ws(body, j)
            if i >= len(body) or body[i] != "(":
                continue
            j = _match_paren_end(body, i)
            nxyz = [int(t) for t in body[i + 1:j - 1].split()]
            i = _skip_ws(body, j)
            mm = re.match(r"simpleGrading\s*\(", body[i:], re.IGNORECASE)
            if not mm:
                continue
            i = i + mm.end() - 1          # position of '(' after simpleGrading
            j = _match_paren_end(body, i)
            gtoks = _split_top_level(body[i + 1:j - 1])
            if len(gtoks) != 3:
                continue
            gxyz = [_parse_grading_token(t) for t in gtoks]
        except ValueError:
            continue
        if len(vid) != 8 or len(nxyz) != 3:
            continue
        xs = [verts[i][0] for i in vid]
        ys = [verts[i][1] for i in vid]
        x0, x1 = min(xs), max(xs)
        y0, y1 = min(ys), max(ys)
        xi = min(range(len(XS)), key=lambda i: abs(XS[i] - x0))
        yi = min(range(len(YS)), key=lambda i: abs(YS[i] - y0))
        blocks.append({
            "xi": xi, "yi": yi,
            "nx": nxyz[0], "ny": nxyz[1], "nz": nxyz[2],
            "gx": gxyz[0], "gy": gxyz[1],
        })

    if not blocks:
        raise ValueError("blockMeshDict: no hex blocks parsed")
    return blocks


def blocks_from_defaults():
    """Fallback: rebuild blocks from the hard-coded geometry + known counts."""
    counts = dict(DEFAULT_COUNTS)
    blocks = []
    for (xi, yi, nxk, nyk, gx, gy) in BLOCK_TABLE:
        blocks.append({
            "xi": xi, "yi": yi,
            "nx": counts[nxk], "ny": counts[nyk], "nz": counts["Nz"],
            "gx": gx, "gy": gy,
        })
    return blocks


def coarsen_blocks(blocks, factor):
    """Reduce each block resolution by ``factor`` (ceil), >=1 cell each."""
    out = []
    for b in blocks:
        nb = dict(b)
        nb["nx"] = max(1, int(math.ceil(b["nx"] / factor)))
        nb["ny"] = max(1, int(math.ceil(b["ny"] / factor)))
        nb["nz"] = max(1, int(math.ceil(b["nz"] / factor)))
        out.append(nb)
    return out


# --------------------------------------------------------------------------- #
#  Fluid mesh reconstruction (structured hex + blockMesh simpleGrading)        #
# --------------------------------------------------------------------------- #

def build_fluid_mesh(blocks):
    """Rebuild the fluid C3D8 mesh from block descriptors.

    Returns
        coords     : list of (x, y, z), index 0 == local node id - 1
        elements   : list of 8-tuples of local node ids (1-based), C3D8 order
        meta       : list of (bi, cx, cy, cz) per element
        fsi_face   : dict mapping (bi, cx, cy, cz) marker -> list of (elem_id,
                     face_label) for the three FSI surfaces
    """
    node_key = {}

    def get_node(x, y, z):
        key = (round(x, 9), round(y, 9), round(z, 9))
        nid = node_key.get(key)
        if nid is None:
            nid = len(coords) + 1
            node_key[key] = nid
            coords.append((x, y, z))
        return nid

    coords = []
    elements = []
    meta = []

    for bi, b in enumerate(blocks):
        xi, yi = b["xi"], b["yi"]
        nx, ny, nz = b["nx"], b["ny"], b["nz"]
        gx, gy = b["gx"], b["gy"]
        x0, x1 = XS[xi], XS[xi + 1]
        y0, y1 = YS[yi], YS[yi + 1]
        z0, z1 = ZS[0], ZS[1]
        Lx, Ly, Lz = x1 - x0, y1 - y0, z1 - z0
        xs = grading_positions(Lx, nx, gx)
        ys = grading_positions(Ly, ny, gy)
        zs = grading_positions(Lz, nz, 1.0)

        grid = {}
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    grid[(i, j, k)] = get_node(
                        x0 + xs[i], y0 + ys[j], z0 + zs[k])

        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    n1 = grid[(i, j, k)]
                    n2 = grid[(i + 1, j, k)]
                    n3 = grid[(i + 1, j + 1, k)]
                    n4 = grid[(i, j + 1, k)]
                    n5 = grid[(i, j, k + 1)]
                    n6 = grid[(i + 1, j, k + 1)]
                    n7 = grid[(i + 1, j + 1, k + 1)]
                    n8 = grid[(i, j + 1, k + 1)]
                    elements.append((n1, n2, n3, n4, n5, n6, n7, n8))
                    meta.append((bi, i, j, k))

    return coords, elements, meta


# Face labels (CalculiX C3D8 convention, with our node ordering
# n1=(i,j,k) .. n8=(i,j+1,k+1)):
#   Face 1: -z   Face 2: +z   Face 3: -y
#   Face 4: +x   Face 5: +y   Face 6: -x
FACE_PX, FACE_NX = 4, 6
FACE_PY, FACE_NY = 5, 3
FACE_PZ, FACE_NZ = 2, 1


def collect_fluid_sets(blocks, meta, elements):
    """Collect fluid node sets (local ids) and FSI element-face map.

    Returns (fsi_nodes, wall_nodes, sym_nodes, fsi_face_elements)
      - fsi_nodes / wall_nodes / sym_nodes : set of local node ids
      - fsi_face_elements : list of (elem_id, face_label) for the FSI surface
    """
    fsi_nodes = set()
    wall_nodes = set()
    sym_nodes = set()
    fsi_faces = []

    for eid, (bi, cx, cy, cz) in enumerate(meta, start=1):
        b = blocks[bi]
        xi, yi = b["xi"], b["yi"]
        nx, ny, nz = b["nx"], b["ny"], b["nz"]
        conn = elements[eid - 1]
        n1, n2, n3, n4, n5, n6, n7, n8 = conn

        # --- FSI surfaces (fluid side of the gate) -----------------------
        # gate left face  x=0.292 : block (xi=1, yi=0) +x face
        # gate right face x=0.304 : block (xi=3, yi=0) -x face
        # gate top face   y=0.080 : block (xi=2, yi=1) -y face
        if xi == 1 and yi == 0 and cx == nx - 1:
            fsi_faces.append((eid, FACE_PX))
            fsi_nodes.update((n2, n3, n7, n6))
        if xi == 3 and yi == 0 and cx == 0:
            fsi_faces.append((eid, FACE_NX))
            fsi_nodes.update((n1, n4, n8, n5))
        if xi == 2 and yi == 1 and cy == 0:
            fsi_faces.append((eid, FACE_NY))
            fsi_nodes.update((n1, n2, n6, n5))

        # --- outer walls (full constraint) --------------------------------
        if yi == 2 and cy == ny - 1:              # TOP  (y=0.365)
            wall_nodes.update((n3, n4, n8, n7))
        if xi == 0 and cx == 0:                   # LEFT (x=0)
            wall_nodes.update((n1, n4, n8, n5))
        if xi == 3 and cx == nx - 1:              # RIGHT (x=0.584)
            wall_nodes.update((n2, n3, n7, n6))
        if yi == 0 and cy == 0:                   # BOTTOM_FLUID (y=0)
            wall_nodes.update((n1, n2, n6, n5))

        # --- symmetry planes (constrain uz only) ---------------------------
        if cz == 0:                                # BACK  (z=0)
            sym_nodes.update((n1, n2, n3, n4))
        if cz == nz - 1:                           # FRONT (z=0.012)
            sym_nodes.update((n5, n6, n7, n8))

    # FSI nodes are exclusively tied to the gate: remove them from the SPC
    # sets so they are not also used as dependent SPC nodes (CalculiX forbids
    # dependent nodes of a *TIE to be used in other SPCs).
    wall_nodes -= fsi_nodes
    sym_nodes -= fsi_nodes

    return fsi_nodes, wall_nodes, sym_nodes, fsi_faces


# --------------------------------------------------------------------------- #
#  Gate mesh (identical to 01_generate_mesh.build_solid_mesh / nsets)            #
# --------------------------------------------------------------------------- #

def gate_node_id(i, j, k, nx, ny):
    return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i + 1


def build_gate_mesh(nx, ny, nz):
    coords = [None] * ((nx + 1) * (ny + 1) * (nz + 1))
    dx = (FSI_X1 - FSI_X0) / nx
    dy = FSI_Y / ny
    dz = LZ / nz
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid = gate_node_id(i, j, k, nx, ny)
                coords[nid - 1] = (FSI_X0 + i * dx, j * dy, k * dz)
    elements = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                elements.append((
                    gate_node_id(i, j, k, nx, ny),
                    gate_node_id(i + 1, j, k, nx, ny),
                    gate_node_id(i + 1, j + 1, k, nx, ny),
                    gate_node_id(i, j + 1, k, nx, ny),
                    gate_node_id(i, j, k + 1, nx, ny),
                    gate_node_id(i + 1, j, k + 1, nx, ny),
                    gate_node_id(i + 1, j + 1, k + 1, nx, ny),
                    gate_node_id(i, j + 1, k + 1, nx, ny),
                ))
    return coords, elements


def build_gate_nsets(nx, ny, nz):
    nbottom, nzmin, nzmax = [], [], []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid = gate_node_id(i, j, k, nx, ny)
                if j == 0:
                    nbottom.append(nid)
                if k == 0:
                    nzmin.append(nid)
                if k == nz:
                    nzmax.append(nid)
    return {"Nbottom": sorted(nbottom), "Nzmin": sorted(nzmin),
            "Nzmax": sorted(nzmax)}


def gate_interface_faces(nx, ny, nz):
    """(elem_id, face_label) for the gate interface surface (master)."""
    faces = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                eid = k * (nx * ny) + j * nx + i + 1
                if i == 0:            # x = 0.292  (-x face)
                    faces.append((eid, FACE_NX))
                if i == nx - 1:       # x = 0.304  (+x face)
                    faces.append((eid, FACE_PX))
                if j == ny - 1:       # y = 0.080  (+y face)
                    faces.append((eid, FACE_PY))
    return faces


# --------------------------------------------------------------------------- #
#  Combined .inp writer                                                        #
# --------------------------------------------------------------------------- #

def write_nset(fh, name, node_ids):
    fh.write("*NSET, NSET=%s\n" % name)
    ids = sorted(node_ids)
    for i in range(0, len(ids), 8):
        fh.write(", ".join(str(n) for n in ids[i:i + 8]) + ",\n")


def write_combined_inp(path, gate_coords, gate_elements, gate_nsets,
                       gate_iface_faces, fluid_coords, fluid_elements,
                       fluid_meta, fluid_sets, fluid_E, n_modes,
                       gate_e, gate_nu, gate_rho, export_mode=False):
    """Write the combined modal-analysis deck.

    fluid_E[i] is the Young's modulus of fluid element i (0-based).  Elements
    sharing the same E value are grouped into a single material + ELSET.
    """
    fsi_nodes, wall_nodes, sym_nodes, fsi_faces = fluid_sets
    n_gate_nodes = len(gate_coords)
    n_gate_elem = len(gate_elements)

    # group fluid elements by E value
    groups = {}
    for eid in range(len(fluid_elements)):
        key = fmt(fluid_E[eid])
        groups.setdefault(key, []).append(eid)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("** Gate + virtual-elastic-body modal analysis (C3D8)\n")
        fh.write("** generated by 02_run_modal_virtual_elastic.py\n")

        # ---- nodes (gate first, then fluid) -----------------------------
        fh.write("*NODE\n")
        for nid, (x, y, z) in enumerate(gate_coords, start=1):
            fh.write("%d, %s, %s, %s\n" % (nid, fmt(x), fmt(y), fmt(z)))
        for idx, (x, y, z) in enumerate(fluid_coords):
            nid = n_gate_nodes + idx + 1
            fh.write("%d, %s, %s, %s\n" % (nid, fmt(x), fmt(y), fmt(z)))

        # ---- gate elements ----------------------------------------------
        fh.write("*ELEMENT, TYPE=C3D8, ELSET=EGATE\n")
        for eid, conn in enumerate(gate_elements, start=1):
            fh.write("%d, %s\n" % (eid, ", ".join(str(n) for n in conn)))

        # ---- fluid elements, grouped by material ------------------------
        felem_base = n_gate_elem
        for gidx, (ekey, eids) in enumerate(sorted(groups.items()), start=1):
            fh.write("*ELEMENT, TYPE=C3D8, ELSET=VE_E%d\n" % gidx)
            for local_eid in eids:
                conn = fluid_elements[local_eid]
                global_conn = [n_gate_nodes + n for n in conn]
                gid = felem_base + local_eid + 1
                fh.write("%d, %s\n" % (gid, ", ".join(str(n) for n in global_conn)))

        # ---- node sets (BCs) --------------------------------------------
        write_nset(fh, "Nbottom", gate_nsets["Nbottom"])
        write_nset(fh, "Nzmin", gate_nsets["Nzmin"])
        write_nset(fh, "Nzmax", gate_nsets["Nzmax"])

        write_nset(fh, "FLUID_WALL", [n_gate_nodes + n for n in wall_nodes])
        write_nset(fh, "FLUID_SYM", [n_gate_nodes + n for n in sym_nodes])

        # ---- materials ---------------------------------------------------
        fh.write("*MATERIAL, NAME=GATE\n")
        fh.write("*ELASTIC\n")
        fh.write("%s, %s\n" % (fmt(gate_e), fmt(gate_nu)))
        fh.write("*DENSITY\n")
        fh.write("%s\n" % fmt(gate_rho))

        for gidx, (ekey, eids) in enumerate(sorted(groups.items()), start=1):
            fh.write("*MATERIAL, NAME=VE_M%d\n" % gidx)
            fh.write("*ELASTIC\n")
            fh.write("%s, %s\n" % (ekey, fmt(FLUID_NU)))
            fh.write("*DENSITY\n")
            fh.write("%s\n" % fmt(FLUID_RHO))

        # ---- solid sections ----------------------------------------------
        fh.write("*SOLID SECTION, ELSET=EGATE, MATERIAL=GATE\n")
        for gidx, (ekey, eids) in enumerate(sorted(groups.items()), start=1):
            fh.write("*SOLID SECTION, ELSET=VE_E%d, MATERIAL=VE_M%d\n"
                     % (gidx, gidx))

        # ---- surfaces + tie ----------------------------------------------
        fh.write("*SURFACE, NAME=Sgate, TYPE=ELEMENT\n")
        for (eid, fl) in gate_iface_faces:
            fh.write("%d, S%d\n" % (eid, fl))
        fh.write("*SURFACE, NAME=Sfluid_fsi, TYPE=NODE\n")
        fsi_global = sorted(n_gate_nodes + n for n in fsi_nodes)
        # nodal surfaces accept one node per line (trailing comma)
        for n in fsi_global:
            fh.write("%d,\n" % n)

        fh.write("*TIE, NAME=gate_fluid, POSITION TOLERANCE=0.001\n")
        fh.write("Sfluid_fsi, Sgate\n")

        # ---- frequency step ----------------------------------------------
        fh.write("*STEP\n")
        fh.write("*FREQUENCY\n")
        fh.write("%d\n" % n_modes)
        if export_mode:
            fh.write("*NODE FILE\n")
            fh.write("U,\n")
        fh.write("*BOUNDARY\n")
        fh.write("Nbottom, 1, 3, 0.0\n")
        fh.write("Nzmin, 3, 3, 0.0\n")
        fh.write("Nzmax, 3, 3, 0.0\n")
        fh.write("FLUID_WALL, 1, 3, 0.0\n")
        fh.write("FLUID_SYM, 3, 3, 0.0\n")
        fh.write("*END STEP\n")

    return len(groups)


# --------------------------------------------------------------------------- #
#  ccx runner + eigenvalue parsing (mirrors run_modal_analysis.py)             #
# --------------------------------------------------------------------------- #

def run_ccx(out_dir, jobname):
    exe = shutil.which(CCX_EXE)
    if exe is None:
        raise RuntimeError("ccx executable not found (%s)" % CCX_EXE)
    proc = subprocess.run(
        [exe, "-i", jobname], cwd=out_dir,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    return proc


def parse_eigenfrequencies(dat_path, n_modes):
    with open(dat_path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    freqs = []
    in_table = False
    for line in lines:
        compact = line.replace(" ", "")
        if not in_table:
            if "MODENO" in compact and "EIGENVALUE" in compact and "FREQUENCY" in compact:
                in_table = True
            continue
        s = line.strip()
        if not s:
            continue
        if "PARTICIPATION" in s or "EFFECTIVE" in s or "MODAL" in s:
            break
        parts = s.split()
        if not parts or not parts[0].isdigit():
            continue
        try:
            if len(parts) == 5:
                hz = _to_float(parts[3])
            elif len(parts) >= 2:
                hz = _to_float(parts[1])
            else:
                continue
        except ValueError:
            continue
        freqs.append(hz)
        if len(freqs) >= n_modes:
            break
    return freqs


# --------------------------------------------------------------------------- #
#  .frd parsing + FDM mode export                                              #
# --------------------------------------------------------------------------- #

# CalculiX .frd node-data records use fixed-width fields:
#     " -1" + I10 node-id + three E12.5 values (36 chars).
# Negative values run together without a separating space (e.g.
# " -1      1858-6.58233E-03-2.72648E-03 0.00000E+00"), so plain split()
# cannot parse them; slice on the fixed columns instead.
_FRD_NODE_ID = slice(3, 13)
_FRD_V1 = slice(13, 25)
_FRD_V2 = slice(25, 37)
_FRD_V3 = slice(37, 49)


def _frd_float(tok):
    return float(tok.strip().replace("D", "E").replace("d", "e"))


def parse_frd(path):
    """Parse a CalculiX .frd result file.

    Returns (nodes, modes):
        nodes : dict node_id -> (x, y, z)   from the 2C coordinate block
        modes : list of [freq, {node_id: (ux, uy, uz)}]  from the 1P blocks
    """
    nodes = {}
    modes = []
    in_nodes = False
    in_disp = False
    current = None

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            raw = line.rstrip("\r\n")
            s = raw.strip()
            if not s:
                continue
            tag = s.split()[0]

            if tag == "2C":
                in_nodes, in_disp, current = True, False, None
                continue
            if tag == "100CL":
                current = [_to_float(s.split()[2]), {}]
                modes.append(current)
                in_nodes, in_disp = False, False
                continue
            if tag == "-4":
                in_disp, in_nodes = True, False
                continue
            if tag == "-3":
                in_nodes, in_disp, current = False, False, None
                continue
            if tag in ("-5", "-2"):
                continue
            if tag != "-1":
                continue
            if not (in_nodes or (in_disp and current is not None)):
                continue
            nid = int(raw[_FRD_NODE_ID])
            v1 = _frd_float(raw[_FRD_V1])
            v2 = _frd_float(raw[_FRD_V2])
            v3 = _frd_float(raw[_FRD_V3])
            if in_nodes:
                nodes[nid] = (v1, v2, v3)
            else:
                current[1][nid] = (v1, v2, v3)

    return nodes, modes


def _write_csv_lines(lines, path):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        for ln in lines:
            fh.write(ln + "\r\n")


def _csv_row(fields):
    return ",".join(fields) + ","


def export_mode_files(mode_dir, fluid_coords, n_gate_nodes, modes):
    """Write FDM mode files into ``mode_dir``.

    Only the fluid (virtual-elastic-body) nodes are exported; the gate nodes
    (global ids 1..n_gate_nodes) are skipped.  ``fluid_coords`` is in local
    order, so the fluid node with local index ``i`` has global id
    ``n_gate_nodes + i + 1`` in the .frd file.
    """
    os.makedirs(mode_dir, exist_ok=True)
    n_fluid = len(fluid_coords)
    n_modes = len(modes)
    g0 = n_gate_nodes + 1

    coor_lines = [
        _csv_row((format_f20(0.0), format_f20(n_fluid), format_f20(n_modes))),
    ]
    for (x, y, z) in fluid_coords:
        coor_lines.append(
            _csv_row((format_f20(x), format_f20(y), format_f20(z))))
    _write_csv_lines(coor_lines, os.path.join(mode_dir, "FluidNodeCoor.csv"))

    for m, (freq, disp) in enumerate(modes, start=1):
        lines = [
            _csv_row((format_e20(freq), format_e20(n_fluid),
                      format_e20(n_modes))),
        ]
        for i in range(n_fluid):
            ux, uy, uz = disp.get(g0 + i, (0.0, 0.0, 0.0))
            lines.append(
                _csv_row((format_e20(ux), format_e20(uy), format_e20(uz))))
        _write_csv_lines(
            lines, os.path.join(mode_dir, "FluidNodeDisp%d.csv" % m))

    return n_fluid, n_modes


def print_comparison(freqs_ccx, ref_ve, ref_bare):
    print("")
    print("[FDM 02] 模态频率对比")
    print("阶数 | CalculiX组合 (Hz) | 论文右列VE (Hz) | 论文左列裸板 (Hz) | 误差 vs VE (%)")
    print("-----+-------------------+-----------------+-------------------+----------------")
    n = min(len(freqs_ccx), len(ref_ve))
    for i in range(n):
        fc = freqs_ccx[i]
        rv = ref_ve[i]
        rb = ref_bare[i] if i < len(ref_bare) else float("nan")
        err = (fc - rv) / rv * 100.0
        print("%4d | %17.4f | %15.4f | %17.4f | %+13.2f" % (i + 1, fc, rv, rb, err))
    if len(freqs_ccx) > n:
        print("      （额外算出的阶数）")
        for i in range(n, len(freqs_ccx)):
            print("%4d | %17.4f |                 |                   |" % (i + 1, freqs_ccx[i]))
    print("")
    return n


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

def parse_args(argv):
    p = argparse.ArgumentParser(
        description="CalculiX modal analysis of gate (real material) + fluid "
                    "domain (virtual elastic body with layered stiffness)."
    )
    p.add_argument("--fluid-dir", default=FLUID_DIR_DEFAULT,
                   help="OpenFOAM fluid case dir (reads system/blockMeshDict; "
                        "default %s)" % FLUID_DIR_DEFAULT)
    p.add_argument("--out-dir", default=OUT_DIR_DEFAULT,
                   help="output dir (default %s)" % OUT_DIR_DEFAULT)
    p.add_argument("--n-modes", type=int, default=10,
                   help="number of eigenmodes (default 10)")
    p.add_argument("--coarse-fluid", action="store_true",
                   help="coarsen the fluid mesh for a quick validation run")
    p.add_argument("--coarse-factor", type=int, default=4,
                   help="coarsening factor for --coarse-fluid (default 4)")
    p.add_argument("--gate-e", type=float, default=GATE_E_DEFAULT,
                   help="gate Young's modulus Pa (default 1e6)")
    p.add_argument("--gate-nu", type=float, default=GATE_NU_DEFAULT,
                   help="gate Poisson ratio (default 0)")
    p.add_argument("--gate-rho", type=float, default=GATE_RHO_DEFAULT,
                   help="gate density kg/m3 (default 2500)")
    p.add_argument("--run-ccx", action="store_true", default=True,
                   help="run ccx (default True)")
    p.add_argument("--no-run-ccx", dest="run_ccx", action="store_false",
                   help="only generate the .inp, do not run ccx")
    p.add_argument("--export-mode", action="store_true",
                   help="export FDM mode files (FluidNodeCoor.csv + "
                        "FluidNodeDisp1..N.csv) from the ccx .frd result")
    p.add_argument("--mode-out-dir", default=None,
                   help="output dir for FDM mode files (default <out-dir>/mode)")
    return p.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    print("[FDM-PIPE] ===== 步骤 2/5: 模态分析与导出 =====")
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # ---- fluid blocks ---------------------------------------------------
    bmd_path = os.path.join(args.fluid_dir, "system", "blockMeshDict")
    if os.path.isfile(bmd_path):
        blocks = parse_block_mesh_dict(bmd_path)
        src = "parsed blockMeshDict (%s)" % bmd_path
    else:
        blocks = blocks_from_defaults()
        src = "hard-coded default counts (blockMeshDict not found)"

    n_cells_full = sum(b["nx"] * b["ny"] for b in blocks) * blocks[0]["nz"]
    print("[FDM 02] 流体块来源   : %s" % src)
    print("[FDM 02] 全分辨率单元 : %d" % n_cells_full)

    if args.coarse_fluid:
        blocks = coarsen_blocks(blocks, args.coarse_factor)
        n_cells = sum(b["nx"] * b["ny"] for b in blocks) * blocks[0]["nz"]
        print("[FDM 02] coarse-fluid : factor=%d -> %d cells" % (args.coarse_factor, n_cells))
    else:
        n_cells = n_cells_full

    # ---- fluid mesh + layered stiffness --------------------------------
    print("[FDM 02] --- 重建流体网格 ...")
    fluid_coords, fluid_elements, fluid_meta = build_fluid_mesh(blocks)
    print("[FDM 02] 流体节点数   : %d" % len(fluid_coords))
    print("[FDM 02] 流体单元数   : %d (C3D8)" % len(fluid_elements))

    # volumes + layered stiffness
    volumes = []
    for conn in fluid_elements:
        pts = [fluid_coords[n - 1] for n in conn]
        volumes.append(hex_volume(pts))
    neg = sum(1 for v in volumes if v <= 0.0)
    if neg:
        print("[FDM WARN] %d fluid element(s) with non-positive volume" % neg)
    min_vol = min(volumes)
    fluid_E = [FLUID_BASE_E * (min_vol / v) ** FLUID_EXP for v in volumes]
    print("[FDM 02] 流体 min_vol  : %.6e" % min_vol)
    print("[FDM 02] 流体 E 范围   : [%.6g, %.6g] Pa" % (min(fluid_E), max(fluid_E)))

    fluid_sets = collect_fluid_sets(blocks, fluid_meta, fluid_elements)
    fsi_nodes, wall_nodes, sym_nodes, fsi_faces = fluid_sets
    print("[FDM 02] 流体节点集   : FSI=%d, wall=%d, sym=%d" %
          (len(fsi_nodes), len(wall_nodes), len(sym_nodes)))
    print("[FDM 02] FSI 界面面片 : %d (fluid)" % len(fsi_faces))

    # ---- gate mesh ------------------------------------------------------
    gate_coords, gate_elements = build_gate_mesh(GATE_NX, GATE_NY, GATE_NZ)
    gate_nsets = build_gate_nsets(GATE_NX, GATE_NY, GATE_NZ)
    gate_iface_faces = gate_interface_faces(GATE_NX, GATE_NY, GATE_NZ)
    print("[FDM 02] 挡板节点/单元 : %d / %d (C3D8)" % (len(gate_coords), len(gate_elements)))
    print("[FDM 02] 挡板界面面片 : %d (gate)" % len(gate_iface_faces))
    print("[FDM 02] 材料参数     : gate E=%.6g, nu=%.6g, rho=%.6g" %
          (args.gate_e, args.gate_nu, args.gate_rho))

    # ---- write combined inp --------------------------------------------
    inp_path = os.path.join(out_dir, JOBNAME + ".inp")
    n_mat = write_combined_inp(
        inp_path, gate_coords, gate_elements, gate_nsets, gate_iface_faces,
        fluid_coords, fluid_elements, fluid_meta, fluid_sets, fluid_E,
        args.n_modes, args.gate_e, args.gate_nu, args.gate_rho,
        export_mode=args.export_mode)
    size_mb = os.path.getsize(inp_path) / 1e6
    print("[FDM 02] 已生成组合 .inp : %s (%.1f MB, %d 个流体材料)" % (inp_path, size_mb, n_mat))

    if not args.run_ccx:
        print("[FDM 02] --no-run-ccx：仅生成 inp，跳过 ccx 求解。")
        if args.export_mode:
            print("[FDM WARN] --export-mode 需要 ccx 生成 .frd 结果，当前已跳过。")
        return 0

    # ---- run ccx --------------------------------------------------------
    proc = run_ccx(out_dir, JOBNAME)
    dat_path = os.path.join(out_dir, JOBNAME + ".dat")
    if proc.returncode != 0:
        print("[FDM ERR] ccx 运行失败 (rc=%d)" % proc.returncode)
        print(proc.stdout[-4000:])
        print(proc.stderr[-4000:])
        if os.path.isfile(dat_path):
            print(open(dat_path, "r", encoding="utf-8", errors="replace").read()[-3000:])
        return 1
    print("[FDM 02] ccx 运行成功 (rc=0)，输出 %s" % dat_path)

    # ---- parse + compare -------------------------------------------------
    freqs_ccx = parse_eigenfrequencies(dat_path, args.n_modes)
    if len(freqs_ccx) < 1:
        print("[FDM ERR] 未能从 %s 解析出特征值。" % dat_path)
        return 1
    print("[FDM 02] 解析到 %d 阶频率：" % len(freqs_ccx))
    for i, f in enumerate(freqs_ccx, start=1):
        print("[FDM 02]   第 %2d 阶 : %.6f Hz" % (i, f))

    print_comparison(freqs_ccx, PAPER_VE, PAPER_BARE)

    # ---- export FDM mode files (FluidNodeCoor.csv + FluidNodeDisp*.csv) ---
    if args.export_mode:
        mode_out_dir = os.path.abspath(
            args.mode_out_dir or os.path.join(out_dir, "mode"))
        frd_path = os.path.join(out_dir, JOBNAME + ".frd")
        if not os.path.isfile(frd_path):
            print("[FDM ERR] 未找到 .frd 结果文件 %s，无法导出模态。" % frd_path)
            return 1
        nodes, modes = parse_frd(frd_path)
        if len(modes) < 1:
            print("[FDM ERR] .frd 中未解析到模态块（请确认已加 *NODE FILE 卡片）。")
            return 1
        if len(modes) < len(freqs_ccx):
            print("[FDM WARN] .frd 解析到 %d 阶模态，但 .dat 有 %d 阶。" %
                  (len(modes), len(freqs_ccx)))
        n_fluid_export, n_modes_export = export_mode_files(
            mode_out_dir, fluid_coords, len(gate_coords), modes)
        print("[FDM 02] 模态导出完成 : %s" % mode_out_dir)
        print("[FDM 02]   FluidNodeCoor.csv      : %d 行 (头行 + %d 流体节点)" %
              (n_fluid_export + 1, n_fluid_export))
        print("[FDM 02]   FluidNodeDisp1..%d.csv : 每文件 %d 行, 共 %d 阶模态" %
              (n_modes_export, n_fluid_export + 1, n_modes_export))
        print("[FDM 02] 提示：导出流体节点数 %d 需与 FDM 使用的 OpenFOAM 网格点数一致，"
              "请核对。" % n_fluid_export)

    f1 = freqs_ccx[0]
    ref1 = PAPER_VE[0]
    rel = abs(f1 - ref1) / abs(ref1)
    print("[FDM 02] 一阶频率: %.4f Hz  模态数: %d  基准误差: %.2f%%"
          % (f1, len(freqs_ccx), rel * 100.0))
    if rel < 0.10:
        print("[FDM 02] PASS: 一阶频率 %.4f Hz 与论文右列 %.2f Hz 相对误差 %.2f%% (< 10%%)" %
              (f1, ref1, rel * 100.0))
        print("[FDM 02]      （虚拟弹性体方法有效：组合系统频率略高于裸板 %.2f Hz）" % PAPER_BARE[0])
        return 0
    print("[FDM WARN] 一阶频率 %.4f Hz 与论文右列 %.2f Hz 相对误差 %.2f%% (>= 10%%)" %
          (f1, ref1, rel * 100.0))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
