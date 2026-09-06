#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""01_generate_mesh.py

Generate the mesh for a partitioned (preCICE) dam-break FSI case:

  * OpenFOAM fluid mesh  (blockMeshDict -> blockMesh -> checkMesh)
  * CalculiX solid mesh  (dam_gate.inp, C3D8, APDL sparse default 868/540)
  * ParaView artefacts   (case.foam, fluid VTK, dam_gate.vtk) via --view

This script ONLY builds and validates the mesh geometry.  It does NOT run any
CFD or FEM solver (no interFoam, no ccx).

Geometry (hard-coded from the Fluent DamFailure.cas.h5 measurements):
    fluid tank  : 0.584 x 0.365 x 0.012 m
    water column: 0.146 x 0.292 m (bottom-left)
    elastic gate: 0.012 x 0.080 x 0.012 m, standing at x in [0.292, 0.304]

Usage:
    python3 01_generate_mesh.py --res 240 --out-dir /root/Workspace/precice_damfailure_new --view
"""

import argparse
import glob
import math
import os
import re
import shutil
import subprocess
import sys

from display_width import pad_right

# --- optional mesh rendering (numpy/matplotlib) ---------------------------
try:
    import numpy as _np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as _plt
    _HAS_RENDER = True
except Exception:
    _np = None
    _plt = None
    _HAS_RENDER = False

# --------------------------------------------------------------------------- #
#  Hard-coded geometry (metres)                                                #
# --------------------------------------------------------------------------- #

LX = 0.584      # tank full length in x
LY = 0.365      # tank full length in y
LZ = 0.012      # tank thickness in z (quasi-2D)

WATER_X = 0.146  # initial water column width
WATER_Y = 0.292  # initial water column height

FSI_X0 = 0.292   # gate left edge x
FSI_X1 = 0.304   # gate right edge x  (12 mm wide)
FSI_Y = 0.080    # gate height (80 mm)

# x / y / z split coordinates
XS = [0.0, WATER_X, FSI_X0, FSI_X1, LX]           # 5 points
YS = [0.0, FSI_Y, WATER_Y, LY]                    # 4 points
ZS = [0.0, LZ]                                    # 2 points

# OpenFOAM bashrc (source before blockMesh/checkMesh/foamToVTK)
OPENFOAM_BASHRC_DEFAULT = "/usr/lib/openfoam/openfoam2412/etc/bashrc"

# --------------------------------------------------------------------------- #
#  Small helpers                                                               #
# --------------------------------------------------------------------------- #

def g(x):
    """Compact float formatting for OpenFOAM / CalculiX ASCII files."""
    if isinstance(x, int):
        return str(x)
    s = format(float(x), ".10g")
    return s


def vid(xi, yi, zi):
    """Global vertex index for the 5x4x2 grid of corner points.

    Vertex ordering (matches the task spec):
        0-4   : y=0,       z=0        (x = 0, .146, .292, .304, .584)
        5-9   : y=0.080,   z=0
        10-14 : y=0.292,   z=0
        15-19 : y=0.365,   z=0
        20-39 : same as 0-19 but z=0.012
    """
    return zi * 20 + yi * 5 + xi


def build_vertices():
    verts = []
    for zi, z in enumerate(ZS):
        for yi, y in enumerate(YS):
            for xi, x in enumerate(XS):
                verts.append((x, y, z))
    return verts


# --------------------------------------------------------------------------- #
#  Density control                                                             #
# --------------------------------------------------------------------------- #

LEVELS = {"coarse": 120, "medium": 180, "fine": 240}

# Max allowed cell-size ratio in any block direction (OpenFOAM quality gate).
# simpleGrading g is the geometric-series ratio of a direction: adjacent cells
# grow by q = g**(1/(n-1)) and the largest cell is g times the smallest.  Every
# grading is clamped to [1/MAX_EXPANSION, MAX_EXPANSION], so BOTH the
# max/min cell ratio (g) and the adjacent-cell ratio (q <= max_exp**(1/(n-1)))
# stay <= MAX_EXPANSION.  When a target first-layer thickness would need a
# larger g, the first layer is relaxed instead (quality beats exact delta).
MAX_EXPANSION = 4.0   # max adjacent-cell size ratio (OpenFOAM quality gate)


def resolve_density(N, z_layers=12):
    """Single-parameter mesh density (see task section 4.3)."""
    nx_a = max(4, int(round(N * 0.146 / 0.584)))              # 保持
    nx_b = max(4, int(round(N * 0.146 / 0.584 * 2.0)))        # 1.5 → 2.0（挡板前加密）
    nx_c = max(2, int(round(N * 0.012 / 0.584 * 5.0)))        # 1.0 → 5.0（挡板上方 x 向加密：2→12列）
    nx_d = max(4, int(round(N * 0.280 / 0.584 * 1.5)))        # 1.2 → 1.5（挡板后加密）
    Ny = max(4, int(round(N * 0.365 / 0.584)))
    ny_l = max(4, int(round(Ny * 0.080 / 0.365 * 2.0)))       # 1.5 → 2.0（挡板高区加密）
    ny_m = max(4, int(round(Ny * 0.212 / 0.365)))             # 保持
    ny_u = max(4, int(round(Ny * 0.073 / 0.365)))             # 3.0 → 1.0 还原
    Nz = z_layers
    return {
        "nx_a": nx_a, "nx_b": nx_b, "nx_c": nx_c, "nx_d": nx_d,
        "ny_l": ny_l, "ny_m": ny_m, "ny_u": ny_u, "Nz": Nz,
    }


# --------------------------------------------------------------------------- #
#  Boundary-layer grading (uniform first-layer thickness on every wall)        #
# --------------------------------------------------------------------------- #

def first_layer_g(L, n, g):
    """First (start) cell thickness of a block of length L, n cells, grading g.

    g > 1 : start cells are denser (smaller).  g < 1 : end cells denser.
    Replicates OpenFOAM simpleGrading: cell widths follow a geometric series
    with ratio q = g^(1/(n-1)).
    """
    if n <= 1:
        return L
    if abs(g - 1) < 1e-12:
        return L / n
    q = g ** (1.0 / (n - 1))
    if g > 1:
        return L * (q - 1) / (q ** n - 1)
    return L * (1 - q) / (1 - q ** n)


def last_layer_g(L, n, g):
    """Last (end) cell thickness of a block."""
    if n <= 1:
        return L
    q = g ** (1.0 / (n - 1))
    return first_layer_g(L, n, g) * q ** (n - 1)


def _exp_bound(max_exp):
    """Effective grading cap for one direction: max_exp with a tiny safety
    margin so a measured cell-size ratio can never exceed max_exp through
    floating-point round-off."""
    return max_exp / (1.0 + 1e-9)


def _clamp_gradient(g, max_exp):
    """Clamp grading g so the adjacent-cell ratio q = g**(1/(n-1)) <= max_exp.

    A direction graded with g has geometric-series ratio q = g**(1/(n-1)) and
    its largest cell is g times its smallest, so capping g at [1/max_exp,
    max_exp] keeps both the max/min cell ratio and the adjacent-cell ratio
    (q <= max_exp**(1/(n-1)) <= max_exp for n >= 2) within the quality gate.
    """
    hi = _exp_bound(max_exp)
    lo = 1.0 / hi
    if g > hi:
        return hi
    if g < lo:
        return lo
    return g


def grading_start(L, n, delta, max_exp=MAX_EXPANSION):
    """Grading ratio g >= 1 making the first (start) cell thickness == delta.

    g is capped at max_exp (quality gate).  When the exact-delta solution
    would exceed the cap, the cap wins and the first-layer thickness is
    relaxed (delta becomes a soft target).
    """
    if n <= 1 or delta >= L / n - 1e-12:
        return 1.0
    lo, hi = 1.0, _exp_bound(max_exp)
    for _ in range(300):
        g = math.sqrt(lo * hi)
        if first_layer_g(L, n, g) > delta:
            lo = g
        else:
            hi = g
    return math.sqrt(lo * hi)


def grading_end(L, n, delta, max_exp=MAX_EXPANSION):
    """Grading ratio g <= 1 making the last (end) cell thickness == delta.

    g is floored at 1/max_exp; when the exact-delta solution would violate the
    quality gate, the target first-layer thickness is relaxed instead.
    """
    if n <= 1 or delta >= L / n - 1e-12:
        return 1.0
    lo, hi = 1.0 / _exp_bound(max_exp), 1.0
    for _ in range(300):
        g = math.sqrt(lo * hi)
        if last_layer_g(L, n, g) > delta:
            hi = g
        else:
            lo = g
    return math.sqrt(lo * hi)


def actual_first_layer(L, n, g):
    """First (start) cell thickness actually produced by a (possibly clamped)
    grading g -- the value to report when the target delta was relaxed."""
    return first_layer_g(L, n, g)


def actual_last_layer(L, n, g):
    """Last (end) cell thickness actually produced by a (possibly clamped)
    grading g -- the value to report when the target delta was relaxed."""
    return last_layer_g(L, n, g)


def multi_seg_divs(n):
    """Cell counts (d0, d1) of a 50/50 two-section split, replicating blockMesh.

    blockMesh rounds each section as int(0.5*n + 0.5) and folds the remainder
    into the first of the equally-weighted sections.
    """
    d0 = int(0.5 * n + 0.5)
    d1 = int(0.5 * n + 0.5)
    if d0 + d1 != n:
        d0 += n - (d0 + d1)
    return d0, d1


def graded_divisions(g, n):
    """Normalized (0..1) edge-point positions for a single graded section."""
    divs = [0.0] * (n + 1)
    divs[n] = 1.0
    if n <= 1:
        return divs
    if abs(g - 1) < 1e-12:
        for i in range(1, n):
            divs[i] = i / n
        return divs
    exp = g ** (1.0 / (n - 1))
    for i in range(1, n):
        divs[i] = (1.0 - exp ** i) / (1.0 - exp ** n)
    return divs


def multi_divisions(n, d0, d1, g1, g2):
    """Normalized edge-point positions for a 50/50 two-section grading."""
    divs = [0.0] * (n + 1)
    divs[n] = 1.0
    if d0 <= 1 or abs(g1 - 1) < 1e-12:
        for i in range(1, d0 + 1):
            divs[i] = 0.5 * i / d0
    else:
        e1 = g1 ** (1.0 / (d0 - 1))
        for i in range(1, d0 + 1):
            divs[i] = 0.5 * (1.0 - e1 ** i) / (1.0 - e1 ** d0)
    if d1 <= 1 or abs(g2 - 1) < 1e-12:
        for i in range(d0 + 1, n + 1):
            divs[i] = 0.5 + 0.5 * (i - d0) / d1
    else:
        e2 = g2 ** (1.0 / (d1 - 1))
        for i in range(d0 + 1, n + 1):
            divs[i] = 0.5 + 0.5 * (1.0 - e2 ** (i - d0)) / (1.0 - e2 ** d1)
    return divs


def direction_divisions(gv, n):
    """Normalized edge-point positions from a grading value (float or tuple)."""
    if isinstance(gv, tuple):
        g1, g2 = gv
        d0, d1 = multi_seg_divs(n)
        return multi_divisions(n, d0, d1, g1, g2)
    return graded_divisions(gv, n)


def wall_edges(counts):
    """(name, length, nCells) for every wall-normal block edge."""
    return [
        ("LEFT",       XS[1] - XS[0], counts["nx_a"]),
        ("BOTTOM",     YS[1] - YS[0], counts["ny_l"]),
        ("TOP",        YS[3] - YS[2], counts["ny_u"]),
        ("gate_left",  XS[2] - XS[1], counts["nx_b"]),
        ("gate_right", XS[4] - XS[3], counts["nx_d"]),
        ("gate_top",   YS[2] - YS[1], counts["ny_m"]),
        ("RIGHT",      XS[4] - XS[3], counts["nx_d"]),
    ]


def resolve_bl_thickness(counts, bl_thickness_mm):
    """Target first-layer thickness in metres (uniform on every wall)."""
    if bl_thickness_mm is not None:
        return bl_thickness_mm * 1e-3
    min_ln = min(L / n for (_name, L, n) in wall_edges(counts))
    return 0.5 * min_ln


def grade_spec(spec, L, n, delta, max_exp=MAX_EXPANSION):
    """Grading ratio (float) or (g1, g2) tuple for one block direction.

    Every returned g is clamped so the adjacent-cell ratio stays <= max_exp
    (see _clamp_gradient); this also caps the per-direction max/min cell ratio.
    """
    if spec == "start":
        return _clamp_gradient(grading_start(L, n, delta, max_exp), max_exp)
    if spec == "end":
        return _clamp_gradient(grading_end(L, n, delta, max_exp), max_exp)
    if spec == "uniform":
        return 1.0
    if spec == "multi":
        d0, d1 = multi_seg_divs(n)
        g1 = _clamp_gradient(grading_start(0.5 * L, d0, delta, max_exp), max_exp)
        g2 = _clamp_gradient(grading_end(0.5 * L, d1, delta, max_exp), max_exp)
        # Rebalance both sections onto a common extreme (small-end) cell
        # thickness so the combined direction's max/min ratio stays within
        # max_exp even when the sections have unequal cell counts.
        f1 = actual_first_layer(0.5 * L, d0, g1)
        f2 = actual_last_layer(0.5 * L, d1, g2)
        d_multi = max(f1, f2)
        return (
            _clamp_gradient(grading_start(0.5 * L, d0, d_multi, max_exp), max_exp),
            _clamp_gradient(grading_end(0.5 * L, d1, d_multi, max_exp), max_exp),
        )
    raise ValueError("unknown grading spec: %r" % (spec,))


# Block table: (xi, yi, nx_key, ny_key, gx_spec, gy_spec).
# gx_spec/gy_spec : "start" = dense at the low-index (start) end,
#                   "end"   = dense at the high-index (end) end,
#                   "uniform" = no grading,
#                   "multi" (x only) = dense at BOTH ends (two 50/50 sections).
BLOCK_TABLE = [
    (0, 0, "nx_a", "ny_l", "start", "start"),   # block 1: LEFT + BOTTOM
    (1, 0, "nx_b", "ny_l", "end", "start"),     # block 2: 挡板左(x=0.292) + BOTTOM
    (3, 0, "nx_d", "ny_l", "multi", "start"),   # block 3: 挡板右+RIGHT + BOTTOM
    (0, 1, "nx_a", "ny_m", "start", "start"),   # block 4: LEFT + 挡板顶
    (1, 1, "nx_b", "ny_m", "end", "start"),     # block 5: 挡板左 + 挡板顶
    (2, 1, "nx_c", "ny_m", "uniform", "start"), # block 6: (no x wall) + 挡板顶
    (3, 1, "nx_d", "ny_m", "multi", "start"),   # block 7: 挡板右+RIGHT + 挡板顶
    (0, 2, "nx_a", "ny_u", "start", "end"),     # block 8: LEFT + TOP
    (1, 2, "nx_b", "ny_u", "end", "end"),       # block 9: 挡板左 + TOP
    (2, 2, "nx_c", "ny_u", "uniform", "end"),   # block 10: (no x wall) + TOP
    (3, 2, "nx_d", "ny_u", "multi", "end"),     # block 11: 挡板右+RIGHT + TOP
]


def hex_corners(xi, yi):
    """8 vertices of a block (blockMesh hex convention: bottom CCW + top)."""
    return [
        vid(xi, yi, 0), vid(xi + 1, yi, 0), vid(xi + 1, yi + 1, 0), vid(xi, yi + 1, 0),
        vid(xi, yi, 1), vid(xi + 1, yi, 1), vid(xi + 1, yi + 1, 1), vid(xi, yi + 1, 1),
    ]


def resolved_blocks(counts, delta, max_exp=MAX_EXPANSION):
    blocks = []
    for (xi, yi, nxk, nyk, gxs, gys) in BLOCK_TABLE:
        nx = counts[nxk]
        ny = counts[nyk]
        Lx = XS[xi + 1] - XS[xi]
        Ly = YS[yi + 1] - YS[yi]
        blocks.append({
            "xi": xi, "yi": yi,
            "nx": nx, "ny": ny, "Nz": counts["Nz"],
            "gx": grade_spec(gxs, Lx, nx, delta, max_exp),
            "gy": grade_spec(gys, Ly, ny, delta, max_exp),
        })
    return blocks


def expected_cell_count(counts):
    total = 0
    for (xi, yi, nxk, nyk, gxs, gys) in BLOCK_TABLE:
        total += counts[nxk] * counts[nyk]
    return total * counts["Nz"]


# --------------------------------------------------------------------------- #
#  OpenFOAM case files                                                         #
# --------------------------------------------------------------------------- #

FOAM_HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2412                                  |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
"""


def fmt_gx(gx):
    """Format a grading value (float) or two-section tuple as blockMesh text."""
    if isinstance(gx, tuple):
        g1, g2 = gx
        return "((0.5 0.5 %s) (0.5 0.5 %s))" % (g(g1), g(g2))
    return g(gx)


def write_block_mesh_dict(case_dir, counts, delta, max_exp=MAX_EXPANSION):
    verts = build_vertices()
    blocks = resolved_blocks(counts, delta, max_exp)

    lines = []
    lines.append(FOAM_HEADER)
    lines.append("FoamFile\n{\n    version     2.0;\n    format      ascii;")
    lines.append("    class       dictionary;\n    location    \"system\";")
    lines.append("    object      blockMeshDict;\n}\n")
    lines.append("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")

    lines.append("scale 1;\n")
    lines.append("vertices\n(")
    for (x, y, z) in verts:
        lines.append("    (%s %s %s)" % (g(x), g(y), g(z)))
    lines.append(");\n")

    lines.append("blocks\n(")
    for b in blocks:
        c = hex_corners(b["xi"], b["yi"])
        hexline = " ".join(str(v) for v in c)
        lines.append(
            "    hex (%s) (%d %d %d) simpleGrading (%s %s 1)"
            % (hexline, b["nx"], b["ny"], b["Nz"], fmt_gx(b["gx"]), g(b["gy"]))
        )
    lines.append(");\n")

    lines.append("edges\n(\n);\n")

    lines.append("boundary\n(")

    # FRONT : z = 0.012 (top faces of all blocks)
    lines.append("    FRONT\n    {\n        type symmetryPlane;\n        faces\n        (")
    for b in blocks:
        xi, yi = b["xi"], b["yi"]
        f = [vid(xi, yi, 1), vid(xi + 1, yi, 1), vid(xi + 1, yi + 1, 1), vid(xi, yi + 1, 1)]
        lines.append("            (%s)" % " ".join(str(v) for v in f))
    lines.append("        );\n    }\n")

    # BACK : z = 0 (bottom faces of all blocks)
    lines.append("    BACK\n    {\n        type symmetryPlane;\n        faces\n        (")
    for b in blocks:
        xi, yi = b["xi"], b["yi"]
        f = [vid(xi, yi, 0), vid(xi + 1, yi, 0), vid(xi + 1, yi + 1, 0), vid(xi, yi + 1, 0)]
        lines.append("            (%s)" % " ".join(str(v) for v in f))
    lines.append("        );\n    }\n")

    # TOP : y = 0.365 (blocks 8, 9, 10, 11)
    lines.append("    TOP\n    {\n        type wall;\n        faces\n        (")
    for xi in (0, 1, 2, 3):
        f = [vid(xi, 3, 0), vid(xi + 1, 3, 0), vid(xi + 1, 3, 1), vid(xi, 3, 1)]
        lines.append("            (%s)" % " ".join(str(v) for v in f))
    lines.append("        );\n    }\n")

    # LEFT : x = 0 (blocks 1, 4, 8)
    lines.append("    LEFT\n    {\n        type wall;\n        faces\n        (")
    for yi in (0, 1, 2):
        f = [vid(0, yi, 0), vid(0, yi + 1, 0), vid(0, yi + 1, 1), vid(0, yi, 1)]
        lines.append("            (%s)" % " ".join(str(v) for v in f))
    lines.append("        );\n    }\n")

    # RIGHT : x = 0.584 (blocks 3, 7, 11)
    lines.append("    RIGHT\n    {\n        type wall;\n        faces\n        (")
    for yi in (0, 1, 2):
        f = [vid(4, yi, 0), vid(4, yi, 1), vid(4, yi + 1, 1), vid(4, yi + 1, 0)]
        lines.append("            (%s)" % " ".join(str(v) for v in f))
    lines.append("        );\n    }\n")

    # BOTTOM_FLUID : y = 0 (blocks 1, 2, 3)
    lines.append("    BOTTOM_FLUID\n    {\n        type wall;\n        faces\n        (")
    for xi in (0, 1, 3):
        f = [vid(xi, 0, 0), vid(xi + 1, 0, 0), vid(xi + 1, 0, 1), vid(xi, 0, 1)]
        lines.append("            (%s)" % " ".join(str(v) for v in f))
    lines.append("        );\n    }\n")

    # FSI_FLUID : three gate surfaces
    lines.append("    FSI_FLUID\n    {\n        type wall;\n        faces\n        (")
    fsi_faces = [
        [2, 7, 27, 22],   # x = 0.292, y in [0, 0.080]
        [3, 8, 28, 23],   # x = 0.304, y in [0, 0.080]
        [7, 8, 28, 27],   # y = 0.080, x in [0.292, 0.304]
    ]
    for f in fsi_faces:
        lines.append("            (%s)" % " ".join(str(v) for v in f))
    lines.append("        );\n    }\n")

    lines.append(");\n")
    lines.append("mergePatchPairs\n(\n);\n")
    lines.append("// ************************************************************************* //\n")

    path = os.path.join(case_dir, "system", "blockMeshDict")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def write_alpha_water(case_dir):
    content = FOAM_HEADER + """FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alpha.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    TOP
    {
        type            zeroGradient;
    }
    LEFT
    {
        type            zeroGradient;
    }
    RIGHT
    {
        type            zeroGradient;
    }
    BOTTOM_FLUID
    {
        type            zeroGradient;
    }
    FSI_FLUID
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
"""
    path = os.path.join(case_dir, "0", "alpha.water")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)


def write_fv_schemes(case_dir):
    content = FOAM_HEADER + """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""
    path = os.path.join(case_dir, "system", "fvSchemes")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)


def write_fv_solution(case_dir):
    content = FOAM_HEADER + """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    pcorr
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-5;
        relTol          0;
    }

    p_rgh
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0.05;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0;
    }
}

PIMPLE
{
    nOuterCorrectors        1;
    nCorrectors             2;
    nNonOrthogonalCorrectors 0;
}

// ************************************************************************* //
"""
    path = os.path.join(case_dir, "system", "fvSolution")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)


def write_control_dict(case_dir):
    content = FOAM_HEADER + """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     interFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         0.5;

deltaT          0.001;

writeControl    timeStep;

writeInterval   50;

purgeWrite      0;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

// ************************************************************************* //
"""
    path = os.path.join(case_dir, "system", "controlDict")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)


def write_case_foam(case_dir):
    path = os.path.join(case_dir, "case.foam")
    with open(path, "w", encoding="utf-8") as fh:
        # ParaView OpenFOAM reader marker; content is ignored.
        fh.write("OpenFOAM case marker for ParaView (built by 01_generate_mesh.py)\n")


# --------------------------------------------------------------------------- #
#  Mesh rendering (optional, numpy/matplotlib)                                 #
# --------------------------------------------------------------------------- #

def find_plot_directory(fluid_dir):
    """Locate the OpenFOAM polyMesh points file (case or case/constant)."""
    for cand in (os.path.join(fluid_dir, "constant", "polyMesh", "points"),
                 os.path.join(fluid_dir, "polyMesh", "points")):
        if os.path.isfile(cand):
            return cand
    return None


def render_mesh_plot(log, fluid_dir, counts, delta, max_exp=MAX_EXPANSION):
    """Render x-y mesh cross-section and x cell-size distribution PNGs.

    Uses the system python (numpy+matplotlib).  Writes two files beside the
    fluid case: mesh_plot_xy.png and mesh_dx_dist.png.  Never raises -- on any
    failure logs a warning and returns False.
    """
    if not _HAS_RENDER:
        log.warn("mesh render skipped (numpy/matplotlib not available)")
        return False
    pts = find_plot_directory(fluid_dir)
    if pts is None:
        log.warn("mesh render skipped (no polyMesh/points found)")
        return False
    try:
        xs, ys = [], []
        with open(pts, encoding="utf-8", errors="replace") as fh:
            fh.readline()          # header: count
            for line in fh:
                line_s = line.replace("(", " ").replace(")", " ").replace("|", " ")
                parts = line_s.split()
                if len(parts) >= 3:
                    try:
                        xs.append(float(parts[0])); ys.append(float(parts[1]))
                    except ValueError:
                        pass
        xs = _np.array(xs); ys = _np.array(ys)
        xl = _np.sort(_np.unique(_np.round(xs, 6)))
        yl = _np.sort(_np.unique(_np.round(ys, 6)))
        dx = _np.diff(xl)

        # ---- plot 1: x-y grid cross section ----
        fig, ax = _plt.subplots(figsize=(16, 7))
        for x in xl:
            ax.vlines(x, yl.min(), yl.max(), color="#2b6cb0", lw=0.2, alpha=0.45)
        for y in yl:
            ax.hlines(y, xl.min(), xl.max(), color="#2b6cb0", lw=0.2, alpha=0.45)
        for lab, a, b in [("water", 0.0, 0.146), ("pre-gate", 0.146, 0.292),
                          ("gate", 0.292, 0.304), ("post-gate", 0.304, 0.584)]:
            ax.axvspan(a, b, color="grey", alpha=0.12)
            ax.text((a + b) / 2, yl.max(), lab, ha="center", fontsize=9)
        ax.axvline(0.292, color="red", lw=1.0, alpha=0.8)
        ax.axvline(0.304, color="red", lw=1.0, alpha=0.8)
        ax.set_title("mesh res cross-section  (expansion<=%.1f)  --  x-y slice" % max_exp, fontsize=13)
        ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
        ax.set_xlim(-0.005, 0.59); ax.set_ylim(-0.005, 0.375)
        ax.text(0.02, 0.02, "points=%d  Max dx=%.1fmm" % (len(xs), dx.max() * 1000), fontsize=9)
        p1 = os.path.join(fluid_dir, "mesh_plot_xy.png")
        _plt.tight_layout(); fig.savefig(p1, dpi=140); _plt.close(fig)

        # ---- plot 2: x cell-size distribution ----
        fig2, ax2 = _plt.subplots(figsize=(14, 4))
        xc = (xl[:-1] + xl[1:]) / 2
        ax2.plot(xc, dx * 1000, "b-", lw=1.2)
        ax2.axvspan(0.292, 0.304, color="grey", alpha=0.2)
        ax2.text(0.298, 1.0, "gate", ha="center", fontsize=9)
        ax2.set_xlabel("x (m)"); ax2.set_ylabel("x cell size (mm)")
        ax2.set_title("x cell-size distribution (max expansion <= %.1f)" % max_exp, fontsize=12)
        ax2.grid(alpha=0.3)
        p2 = os.path.join(fluid_dir, "mesh_dx_dist.png")
        _plt.tight_layout(); fig2.savefig(p2, dpi=140); _plt.close(fig2)

        log.log("wrote mesh_plot_xy.png  (points=%d, Max dx=%.1fmm)" % (len(xs), dx.max() * 1000))
        log.log("wrote mesh_dx_dist.png")
        return True
    except Exception as exc:
        log.warn("mesh render failed (non-fatal): %s" % exc)
        return False


# --------------------------------------------------------------------------- #
#  CalculiX solid mesh (C3D8)                                                  #
# --------------------------------------------------------------------------- #

def node_id(i, j, k, nx, ny):
    """1-based CalculiX node id.  x fastest, then y, then z."""
    return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i + 1


def build_solid_mesh(nx, ny, nz):
    """Return (coords, elements) for the structured C3D8 gate mesh.

    coords[node_id - 1] = (x, y, z)
    elements[elem_id - 1] = (n1..n8) in CalculiX C3D8 order.
    """
    n_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    coords = [None] * n_nodes

    dx = (FSI_X1 - FSI_X0) / nx
    dy = FSI_Y / ny
    dz = LZ / nz

    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid = node_id(i, j, k, nx, ny)
                coords[nid - 1] = (FSI_X0 + i * dx, j * dy, k * dz)

    elements = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                n1 = node_id(i, j, k, nx, ny)
                n2 = node_id(i + 1, j, k, nx, ny)
                n3 = node_id(i + 1, j + 1, k, nx, ny)
                n4 = node_id(i, j + 1, k, nx, ny)
                n5 = node_id(i, j, k + 1, nx, ny)
                n6 = node_id(i + 1, j, k + 1, nx, ny)
                n7 = node_id(i + 1, j + 1, k + 1, nx, ny)
                n8 = node_id(i, j + 1, k + 1, nx, ny)
                elements.append((n1, n2, n3, n4, n5, n6, n7, n8))
    return coords, elements


def build_solid_nsets(nx, ny, nz):
    """Return dict of nset name -> sorted list of node ids (1-based)."""
    nbottom = []
    nzmin = []
    nzmax = []
    interface = set()

    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid = node_id(i, j, k, nx, ny)
                if j == 0:
                    nbottom.append(nid)
                if k == 0:
                    nzmin.append(nid)
                if k == nz:
                    nzmax.append(nid)
                if i == 0 or i == nx or j == ny:
                    interface.add(nid)

    return {
        "Nbottom": sorted(nbottom),
        "Nzmin": sorted(nzmin),
        "Nzmax": sorted(nzmax),
        "Ninterface": sorted(interface),
    }


def write_solid_inp(path, coords, elements, nsets, e_mod, nu, rho):
    lines = []
    lines.append("** Dam-break gate mesh (C3D8), generated by 01_generate_mesh.py")
    lines.append("*NODE")
    for nid, (x, y, z) in enumerate(coords, start=1):
        lines.append("%d, %s, %s, %s" % (nid, g(x), g(y), g(z)))

    lines.append("*ELEMENT, TYPE=C3D8, ELSET=Eall")
    for eid, conn in enumerate(elements, start=1):
        lines.append("%d, %s" % (eid, ", ".join(str(n) for n in conn)))

    for name, nodes in nsets.items():
        lines.append("*NSET, NSET=%s" % name)
        for idx in range(0, len(nodes), 8):
            chunk = nodes[idx:idx + 8]
            lines.append(", ".join(str(n) for n in chunk) + ",")

    lines.append("*MATERIAL, NAME=GATE")
    lines.append("*ELASTIC")
    lines.append("%s, %s" % (g(e_mod), g(nu)))
    lines.append("*DENSITY")
    lines.append(g(rho))
    lines.append("*SOLID SECTION, ELSET=Eall, MATERIAL=GATE")
    lines.append("*STEP, INC=1000000")
    lines.append("*DYNAMIC, ALPHA=0.0, DIRECT")
    lines.append("0.0005, 0.5")
    lines.append("*RESTART,WRITE,FREQUENCY=1")
    lines.append("*BOUNDARY")
    lines.append("Nbottom, 1, 3, 0.0")
    lines.append("Nzmin, 3, 3, 0.0")
    lines.append("Nzmax, 3, 3, 0.0")
    lines.append("*CLOAD")
    lines.append("Ninterface, 1, 0.0")
    lines.append("Ninterface, 2, 0.0")
    lines.append("Ninterface, 3, 0.0")
    lines.append("*END STEP")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# --- solid self-check (pure-Python volume via cross product) -----------------

def tet_signed_volume(a, b, c, d):
    """Signed volume/6 of a tetrahedron (a,b,c,d)."""
    bx, by, bz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    cx, cy, cz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    dx, dy, dz = d[0] - a[0], d[1] - a[1], d[2] - a[2]
    # (b x c) . d
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


def check_solid_mesh(path, nx, ny, nz):
    """Parse the generated .inp and verify node/element/nset counts + volumes."""
    expected_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    expected_elements = nx * ny * nz
    expected_nsets = {
        "Nbottom": (nx + 1) * (nz + 1),
        "Nzmin": (nx + 1) * (ny + 1),
        "Nzmax": (nx + 1) * (ny + 1),
        "Ninterface": None,  # computed below
    }

    text = open(path, "r", encoding="utf-8").read()
    lines = text.splitlines()

    nodes = {}
    elements = {}
    nsets = {}

    section = None
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("*"):
            section = line.split(",")[0].strip()
            if section == "*NODE":
                section = "NODE"
            elif section == "*ELEMENT":
                section = "ELEMENT"
            elif section == "*NSET":
                name = line.split("NSET=", 1)[1].split(",")[0].strip()
                section = ("NSET", name)
            else:
                section = None
            continue

        if section == "NODE":
            parts = [p.strip() for p in line.split(",")]
            nid = int(parts[0])
            nodes[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
        elif section == "ELEMENT":
            parts = [int(p.strip()) for p in line.split(",")]
            elements[parts[0]] = tuple(parts[1:9])
        elif isinstance(section, tuple) and section[0] == "NSET":
            name = section[1]
            ids = [int(p.strip()) for p in line.rstrip(",").split(",") if p.strip()]
            nsets.setdefault(name, []).extend(ids)

    errors = []

    if len(nodes) != expected_nodes:
        errors.append("node count %d != %d" % (len(nodes), expected_nodes))
    if len(elements) != expected_elements:
        errors.append("element count %d != %d" % (len(elements), expected_elements))

    # element volume check (all positive)
    min_vol = None
    neg = 0
    for eid, conn in elements.items():
        pts = [nodes[n] for n in conn]
        vol = hex_volume(pts)
        if min_vol is None or vol < min_vol:
            min_vol = vol
        if vol <= 0.0:
            neg += 1
    if neg:
        errors.append("%d element(s) with non-positive volume" % neg)

    # nset counts
    for name, want in expected_nsets.items():
        got = len(nsets.get(name, []))
        if name == "Ninterface":
            # union of i==0, i==nx, j==ny faces, deduped
            iface = set()
            for k in range(nz + 1):
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        if i == 0 or i == nx or j == ny:
                            iface.add(node_id(i, j, k, nx, ny))
            want = len(iface)
        if got != want:
            errors.append("nset %s count %d != %d" % (name, got, want))

    if errors:
        raise RuntimeError("solid mesh self-check FAILED: " + "; ".join(errors))

    return {
        "nodes": len(nodes),
        "elements": len(elements),
        "min_volume": min_vol,
        "nsets": {name: len(nsets.get(name, [])) for name in expected_nsets},
        "nsets_expected": expected_nsets,
    }


# --------------------------------------------------------------------------- #
#  VTK writers                                                                 #
# --------------------------------------------------------------------------- #

def write_solid_vtk(path, coords, elements):
    n = len(coords)
    lines = []
    lines.append("# vtk DataFile Version 3.0")
    lines.append("dam_gate solid mesh (C3D8)")
    lines.append("ASCII")
    lines.append("DATASET UNSTRUCTURED_GRID")
    lines.append("POINTS %d double" % n)
    for (x, y, z) in coords:
        lines.append("%s %s %s" % (g(x), g(y), g(z)))
    lines.append("CELLS %d %d" % (len(elements), len(elements) * 9))
    for conn in elements:
        lines.append("8 " + " ".join(str(n - 1) for n in conn))  # 0-based
    lines.append("CELL_TYPES %d" % len(elements))
    for _ in elements:
        lines.append("12")  # VTK_HEXAHEDRON
    lines.append("POINT_DATA %d" % n)
    lines.append("VECTORS DISPLACEMENT double")
    for _ in coords:
        lines.append("0.0 0.0 0.0")
    lines.append("CELL_DATA %d" % len(elements))
    lines.append("SCALARS Volume double 1")
    lines.append("LOOKUP_TABLE default")
    for conn in elements:
        pts = [coords[c - 1] for c in conn]
        lines.append("%s" % g(hex_volume(pts)))

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def write_fluid_vtk_direct(path, counts, delta, max_exp=MAX_EXPANSION):
    """Fallback fluid VTK: rebuild graded hex cells straight from block geometry."""
    blocks = resolved_blocks(counts, delta, max_exp)

    points = []
    cells = []
    for b in blocks:
        xi, yi = b["xi"], b["yi"]
        nx, ny, nz = b["nx"], b["ny"], b["Nz"]
        Lx = XS[xi + 1] - XS[xi]
        Ly = YS[yi + 1] - YS[yi]
        xdiv = direction_divisions(b["gx"], nx)
        ydiv = direction_divisions(b["gy"], ny)
        zdiv = [k / nz for k in range(nz + 1)]

        x0, y0 = XS[xi], YS[yi]
        base = len(points)
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    points.append((
                        x0 + xdiv[i] * Lx,
                        y0 + ydiv[j] * Ly,
                        zdiv[k] * LZ,
                    ))

        def pidx(i, j, k):
            return base + k * (nx + 1) * (ny + 1) + j * (nx + 1) + i

        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    cells.append([
                        pidx(i, j, k), pidx(i + 1, j, k),
                        pidx(i + 1, j + 1, k), pidx(i, j + 1, k),
                        pidx(i, j, k + 1), pidx(i + 1, j, k + 1),
                        pidx(i + 1, j + 1, k + 1), pidx(i, j + 1, k + 1),
                    ])

    lines = []
    lines.append("# vtk DataFile Version 3.0")
    lines.append("fluid mesh (hex)")
    lines.append("ASCII")
    lines.append("DATASET UNSTRUCTURED_GRID")
    lines.append("POINTS %d double" % len(points))
    for (x, y, z) in points:
        lines.append("%s %s %s" % (g(x), g(y), g(z)))
    lines.append("CELLS %d %d" % (len(cells), len(cells) * 9))
    for c in cells:
        lines.append("8 " + " ".join(str(v) for v in c))
    lines.append("CELL_TYPES %d" % len(cells))
    for _ in cells:
        lines.append("12")
    lines.append("CELL_DATA %d" % len(cells))
    lines.append("SCALARS cellId double 1")
    lines.append("LOOKUP_TABLE default")
    for idx in range(len(cells)):
        lines.append(str(idx))

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
#  OpenFOAM tool runner                                                         #
# --------------------------------------------------------------------------- #

def locate_bashrc():
    exe = shutil.which("blockMesh")
    if exe:
        proj = os.path.dirname(os.path.dirname(os.path.dirname(exe)))
        rc = os.path.join(proj, "etc", "bashrc")
        if os.path.isfile(rc):
            return rc
    if os.path.isfile(OPENFOAM_BASHRC_DEFAULT):
        return OPENFOAM_BASHRC_DEFAULT
    return None


def run_foam(cmd, cwd):
    rc = locate_bashrc()
    if rc:
        shell = 'source "%s" >/dev/null 2>&1 && %s' % (rc, cmd)
    else:
        shell = cmd
    proc = subprocess.run(
        ["bash", "-c", shell], cwd=cwd,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    return proc


# --------------------------------------------------------------------------- #
#  Logging                                                                     #
# --------------------------------------------------------------------------- #

class Logger:
    def __init__(self, path, prefix="[FDM 01]"):
        self.path = path
        self.prefix = prefix
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.fh = open(path, "w", encoding="utf-8")

    def _emit(self, msg, prefix):
        line = ("%s %s" % (prefix, msg)) if msg else prefix
        self.fh.write(line + "\n")
        self.fh.flush()
        print(line)

    def log(self, msg=""):
        self._emit(msg, self.prefix)

    def raw(self, msg=""):
        self.fh.write(msg + "\n")
        self.fh.flush()
        print(msg)

    def error(self, msg):
        self._emit(msg, "[FDM ERR]")

    def warn(self, msg):
        self._emit(msg, "[FDM WARN]")

    def close(self):
        self.fh.close()


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

def parse_block_cell_sizes(blockmesh_stdout):
    """Parse blockMesh verbose output -> list of {'i':(beg,end), 'j':(...), ...}."""
    blocks = []
    cur = None
    for line in blockmesh_stdout.splitlines():
        s = line.strip()
        m = re.match(r"Block\s+(\d+)\s+cell\s+size", s)
        if m:
            cur = {"i": None, "j": None, "k": None}
            blocks.append(cur)
            continue
        if cur is not None:
            m2 = re.match(
                r"([ijk])\s*:\s*([0-9eE+.\-]+)\s*\.\.\s*([0-9eE+.\-]+)", s
            )
            if m2:
                cur[m2.group(1)] = (float(m2.group(2)), float(m2.group(3)))
    return blocks


# Wall -> (block index in blockMesh order, direction, 'beg'/'end').
# blockMesh orders blocks exactly as written (BLOCK_TABLE order).
WALL_FIRST_CELL = [
    ("LEFT",       0, "i", "beg"),   # block 1  x-start (x=0)
    ("gate_left",  1, "i", "end"),   # block 2  x-end   (x=0.292)
    ("gate_right", 2, "i", "beg"),   # block 3  x-start (x=0.304, multi seg 1)
    ("RIGHT",      2, "i", "end"),   # block 3  x-end   (x=0.584, multi seg 2)
    ("BOTTOM",     0, "j", "beg"),   # block 1  y-start (y=0)
    ("gate_top",   5, "j", "beg"),   # block 6  y-start (y=0.080)
    ("TOP",        7, "j", "end"),   # block 8  y-end   (y=0.365)
]


def _wall_expected_first_layer(blocks, bi, direction, end):
    """First-layer thickness (m) a wall should have under the (clamped) grading.

    Lets the boundary-layer report distinguish an intentional relaxation
    (expansion limit) from a real blockMesh/grading mismatch.
    """
    b = blocks[bi]
    if direction == "i":
        L = XS[b["xi"] + 1] - XS[b["xi"]]
        n = b["nx"]
        gv = b["gx"]
    else:
        L = YS[b["yi"] + 1] - YS[b["yi"]]
        n = b["ny"]
        gv = b["gy"]
    if isinstance(gv, tuple):
        d0, d1 = multi_seg_divs(n)
        if end == "beg":
            return actual_first_layer(0.5 * L, d0, gv[0])
        return actual_last_layer(0.5 * L, d1, gv[1])
    if end == "beg":
        return actual_first_layer(L, n, gv)
    return actual_last_layer(L, n, gv)


def report_boundary_layers(log, counts, delta, blockmesh_stdout,
                           max_exp=MAX_EXPANSION):
    """Log target vs actual first-layer thickness on every wall.

    Walls whose grading hit the expansion cap no longer match the exact delta
    target; they are reported with the actual value plus a
    '(relaxed: expansion limited)' note instead of a hard MISMATCH.  A wall
    only fails when blockMesh deviates from the (clamped) grading itself, so
    the uniform-target check stays soft while the per-wall spec check stays
    strict.
    """
    parsed = parse_block_cell_sizes(blockmesh_stdout)
    blocks = resolved_blocks(counts, delta, max_exp)
    log.log("--- boundary-layer report ---")
    log.raw("expansion ratio cap   : %.2f (adjacent-cell / max-min ratio)"
            % max_exp)
    log.raw("target first-layer thickness : %.4f mm (%.6g m)"
            % (delta * 1e3, delta))
    ok = True
    relaxed = 0
    for name, bi, direction, end in WALL_FIRST_CELL:
        if bi >= len(parsed) or parsed[bi].get(direction) is None:
            log.raw("  %s : block %d %s data missing"
                    % (pad_right(name, 10), bi, direction))
            ok = False
            continue
        beg, fin = parsed[bi][direction]
        actual = beg if end == "beg" else fin
        expected = _wall_expected_first_layer(blocks, bi, direction, end)
        rel = (actual - delta) / delta
        rel_spec = abs(actual - expected) / expected if expected > 0 else 0.0
        if abs(rel) < 0.20:
            status = "OK"
        elif rel_spec >= 0.20:
            ok = False
            status = "MISMATCH"
        else:
            relaxed += 1
            status = "(relaxed: expansion limited)"
        log.raw("  %s : %.4f mm  (%+.1f%%)  %s"
                % (pad_right(name, 10), actual * 1e3, rel * 100, status))
    if relaxed:
        log.warn("boundary-layer: %d wall(s) relaxed first-layer thickness "
                 "(expansion ratio capped at %.2f); delta not met exactly, "
                 "mesh quality preserved." % (relaxed, max_exp))
    log.raw("boundary-layer uniformity : %s" % ("PASS" if ok else "FAIL"))
    return ok


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Generate and validate the preCICE dam-break FSI mesh "
                    "(fluid blockMesh + CalculiX C3D8 gate)."
    )
    p.add_argument("--res", type=int, default=None,
                   help="x-direction reference cell count (default 240)")
    p.add_argument("--level", choices=list(LEVELS), default=None,
                   help="preset: coarse=120, medium=180, fine=240")
    p.add_argument("--out-dir", default="/root/Workspace/precice_damfailure_new",
                   help="output case root directory")
    p.add_argument("--view", action="store_true",
                   help="generate ParaView artefacts (case.foam, VTK files)")
    p.add_argument("--render", action="store_true",
                   help="render fluid mesh to PNGs (mesh_plot_xy.png, mesh_dx_dist.png)")
    p.add_argument("--gate-nx", type=int, default=None,
                   help="solid mesh x cells (default 6)")
    p.add_argument("--gate-ny", type=int, default=None,
                   help="solid mesh y cells (default 30)")
    p.add_argument("--gate-nz", type=int, default=None,
                   help="solid mesh z cells (default 3)")
    p.add_argument("--gate-e", type=float, default=1.0e6,
                   help="gate Young's modulus Pa (default 1e6)")
    p.add_argument("--gate-nu", type=float, default=0.0,
                   help="gate Poisson ratio (default 0)")
    p.add_argument("--gate-rho", type=float, default=2500.0,
                   help="gate density kg/m3 (default 2500)")
    p.add_argument("--z-layers", type=int, default=12,
                   help="number of cell layers in z (APDL uses 12, 1mm/layer)")
    p.add_argument("--bl-thickness", type=float, default=None,
                   help="uniform first-layer thickness in mm on every wall "
                        "(default: 0.5 x smallest wall-normal cell size)")
    p.add_argument("--max-expansion", type=float, default=MAX_EXPANSION,
                   help="max cell-size ratio (adjacent-cell and max/min) "
                        "allowed in any block direction; first-layer thickness "
                        "is relaxed when needed to respect it "
                        "(default %.2f)" % MAX_EXPANSION)
    return p.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    if args.res is not None:
        N = args.res
    elif args.level is not None:
        N = LEVELS[args.level]
    else:
        N = 240

    gate_nx = args.gate_nx if args.gate_nx is not None else 6
    gate_ny = args.gate_ny if args.gate_ny is not None else 30
    gate_nz = args.gate_nz if args.gate_nz is not None else 3

    out_dir = os.path.abspath(args.out_dir)
    fluid_dir = os.path.join(out_dir, "fluid-openfoam")
    solid_dir = os.path.join(out_dir, "solid-calculix")

    log = Logger(os.path.join(out_dir, "build_log.txt"))
    log.raw("[FDM-PIPE] ===== 步骤 1/5: 网格生成 =====")
    log.log("res = %d, level = %s, out_dir = %s" % (N, args.level, out_dir))

    max_exp = args.max_expansion
    if max_exp < 1.0:
        log.error("--max-expansion must be >= 1.0 (got %r)" % max_exp)
        return 1

    rc = locate_bashrc()
    if rc is None:
        log.error("OpenFOAM bashrc not found; blockMesh/checkMesh unavailable.")
        return 1
    log.log("OpenFOAM bashrc: %s" % rc)

    counts = resolve_density(N, args.z_layers)
    delta = resolve_bl_thickness(counts, args.bl_thickness)
    log.log("fluid counts: %s" % counts)
    log.log("expected fluid cells: %d" % expected_cell_count(counts))
    log.log("target first-layer thickness: %.4f mm" % (delta * 1e3))
    log.log("max expansion ratio: %.2f (adjacent-cell / max-min)" % max_exp)

    # ---- fluid case -----------------------------------------------------
    os.makedirs(os.path.join(fluid_dir, "system"), exist_ok=True)
    os.makedirs(os.path.join(fluid_dir, "0"), exist_ok=True)

    write_block_mesh_dict(fluid_dir, counts, delta, max_exp)
    write_alpha_water(fluid_dir)
    write_fv_schemes(fluid_dir)
    write_fv_solution(fluid_dir)
    write_control_dict(fluid_dir)

    log.log("--- blockMesh ---")
    proc = run_foam("blockMesh", fluid_dir)
    bm_out = proc.stdout
    log.raw(bm_out.rstrip())
    if proc.stderr.strip():
        log.raw("stderr: " + proc.stderr.rstrip())
    if proc.returncode != 0:
        log.error("blockMesh failed (rc=%d)" % proc.returncode)
        return 1

    report_boundary_layers(log, counts, delta, bm_out, max_exp)

    log.log("--- checkMesh ---")
    proc = run_foam("checkMesh", fluid_dir)
    check_out = proc.stdout + "\n" + proc.stderr
    log.raw(check_out.rstrip())
    if proc.returncode != 0:
        log.error("checkMesh failed (rc=%d)" % proc.returncode)
        return 1

    mesh_ok = "Mesh OK" in check_out
    log.log("checkMesh mesh_ok = %s" % mesh_ok)
    if not mesh_ok:
        log.error("checkMesh did not report 'Mesh OK'.")
        return 1

    # ---- solid case -----------------------------------------------------
    os.makedirs(solid_dir, exist_ok=True)
    coords, elements = build_solid_mesh(gate_nx, gate_ny, gate_nz)
    nsets = build_solid_nsets(gate_nx, gate_ny, gate_nz)
    inp_path = os.path.join(solid_dir, "dam_gate.inp")
    write_solid_inp(inp_path, coords, elements, nsets,
                    args.gate_e, args.gate_nu, args.gate_rho)

    log.log("--- solid mesh self-check ---")
    solid_info = check_solid_mesh(inp_path, gate_nx, gate_ny, gate_nz)
    log.log("solid nodes        = %d" % solid_info["nodes"])
    log.log("solid elements     = %d" % solid_info["elements"])
    log.log("solid min volume   = %.6e" % solid_info["min_volume"])
    log.log("solid nsets        = %s" % solid_info["nsets"])
    log.log("solid self-check PASSED")

    # ---- ParaView artefacts --------------------------------------------
    if args.view:
        write_case_foam(fluid_dir)
        log.log("wrote fluid-openfoam/case.foam")

        vtk_solid = os.path.join(solid_dir, "dam_gate.vtk")
        write_solid_vtk(vtk_solid, coords, elements)
        log.log("wrote solid-calculix/dam_gate.vtk")

        log.log("--- foamToVTK (fluid) ---")
        proc = run_foam("foamToVTK", fluid_dir)
        if proc.returncode == 0:
            log.raw(proc.stdout.rstrip())
        else:
            log.raw("foamToVTK stderr: " + proc.stderr.rstrip())
        vtk_files = []
        for pat in ("*.vtk", "*.vtu", "*.vtp"):
            vtk_files.extend(glob.glob(os.path.join(fluid_dir, "VTK", "**", pat), recursive=True))
        vtk_files = sorted(vtk_files)
        if vtk_files:
            log.log("fluid VTK files: %s" % vtk_files)
        else:
            log.log("foamToVTK produced no VTK; falling back to direct writer")
            vtk_dir = os.path.join(fluid_dir, "VTK")
            os.makedirs(vtk_dir, exist_ok=True)
            direct_path = os.path.join(vtk_dir, "fluid_dam_gate.vtk")
            write_fluid_vtk_direct(direct_path, counts, delta, max_exp)
            log.log("wrote %s" % direct_path)

    # ---- mesh render (optional) -------------------------------------------
    if args.render:
        log.log("--- mesh render ---")
        render_mesh_plot(log, fluid_dir, counts, delta, max_exp)

    # ---- final report ---------------------------------------------------
    log.log("=== summary ===")
    log.log("fluid case        : %s" % fluid_dir)
    log.log("solid inp         : %s" % inp_path)
    log.log("checkMesh         : Mesh OK")
    log.log("expected cells    : %d" % expected_cell_count(counts))

    cells_s = points_s = ""
    for line in check_out.splitlines():
        s = line.strip()
        if s.startswith("cells:"):
            cells_s = s.split(":", 1)[1].strip()
            log.log("actual cells      : " + s)
        if s.startswith("points:"):
            points_s = s.split(":", 1)[1].strip()
            log.log("actual points     : " + s)
        if s.startswith("faces:"):
            log.log("actual faces      : " + s)
        if "Max skewness" in s or "non-orthogonality Max" in s or "Mesh OK" in s:
            log.log("checkMesh metric  : " + s)

    log.log("cells: %s  points: %s  Mesh OK" % (cells_s, points_s))

    log.log("DONE")
    log.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
