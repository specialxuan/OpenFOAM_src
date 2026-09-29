#!/usr/bin/env python3
"""Topology-preserving node smoothing for an OpenFOAM ``polyMesh`` case.

Purpose
-------
The AGARD 445.6 far-field ``polyMesh`` (388,465 points / 372,760 hex cells)
contains thin cells at the wing tip and in the wake sheet whose smallest
dimension is ~1e-6 m.  Those slivers force a very small acoustic timestep on a
compressible solver.  This tool moves **point coordinates only** - the point
count, the point ordering, the faces/owner/neighbour topology and the patch
definitions are preserved byte-for-byte - in an attempt to raise the smallest
cell dimension while keeping every cell volume positive.

The tool is standalone (numpy only) and deterministic: the same input and
options always produce the same output.  It is re-runnable and never
overwrites the input mesh; results go to ``<out-dir>/constant/polyMesh``.

Algorithm
---------
A weighed Laplacian / Taubin smoothing of the free (interior) nodes:

* ``laplacian`` : ``x_i <- x_i + relax * (mean_j x_j - x_i)``
* ``taubin``    : alternating ``relax`` / ``-mu`` Laplacian passes (volume
  preserving; suppresses the shrinkage of plain Laplacian)
* ``expand``    : Laplacian plus a short-edge repulsion term that pushes nodes
  apart when a neighbour is closer than ``--target-edge``

Every iteration is guarded: the full cell-volume / aspect-ratio field is
recomputed, and a step that creates a non-positive cell volume (or worsens the
worst aspect ratio beyond ``--ar-slack``) is rejected and retried with a
halved step.  The best state seen (by smallest cell dimension, then worst
aspect ratio) is written out.

Fixed nodes
-----------
* ``FSI`` patch nodes are ALWAYS fixed (the wing surface must not move).
* All other boundary-patch nodes are fixed by default.  ``--free-boundary``
  releases them; with ``--slide-boundary`` they are constrained to their own
  patch plane (x/z planes, y=0 and y=7 planes, z=+-5 planes) so the domain
  shape cannot change, while the half-cylinder ``IN`` patch (if present) is
  kept on its cylinder.

CLI
---
    smooth_farfield_polymesh.py --input <polyMesh dir> --out-dir <case dir> \
        [--iters N] [--relax L] [--mode laplacian|taubin|expand] \
        [--target-edge H] [--repel K] [--region all|bbox] \
        [--bbox X0 X1 Y0 Y1 Z0 Z1] [--free-boundary] [--slide-boundary] \
        [--ar-slack F] [--objective sliver|minedge|maxar] \
        [--report PATH] [--help]
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

_SKIP_HEADER = ("FoamFile", "version", "format", "class", "location", "object",
                "{", "}", "(", ")", "//", "/*", "*", "\\", "|")


def _read_list_head(path: Path) -> tuple[int, list[str], int]:
    lines = path.read_text(encoding="ascii").splitlines()
    i = 0
    n = None
    while i < len(lines):
        s = lines[i].strip()
        if s and s[0].isdigit():
            n = int(s)
            i += 1
            break
        i += 1
    if n is None:
        raise SystemExit(f"[smooth] ERROR: could not read list size from {path}")
    return n, lines, i


def read_points(path: Path) -> np.ndarray:
    n, lines, i = _read_list_head(path)
    pts = np.empty((n, 3), dtype=np.float64)
    j = 0
    while j < n:
        s = lines[i].strip()
        i += 1
        if not s or s in ("(", ")"):
            continue
        s = s.strip("()")
        f = s.split()
        pts[j, 0] = float(f[0]); pts[j, 1] = float(f[1]); pts[j, 2] = float(f[2])
        j += 1
    return pts


def read_faces(path: Path) -> np.ndarray:
    n, lines, i = _read_list_head(path)
    faces = np.empty((n, 4), dtype=np.int64)
    j = 0
    while j < n:
        s = lines[i].strip()
        i += 1
        if not s or s == "(":
            continue
        idx = s.replace(")", "").split("(")[-1].split()
        if len(idx) != 4:
            raise SystemExit(f"[smooth] ERROR: non-quad face at line {i}: {s!r}")
        faces[j] = [int(idx[0]), int(idx[1]), int(idx[2]), int(idx[3])]
        j += 1
    return faces


def read_labels(path: Path) -> np.ndarray:
    n, lines, i = _read_list_head(path)
    out = np.empty(n, dtype=np.int64)
    j = 0
    while j < n:
        s = lines[i].strip()
        i += 1
        if not s or s in ("(", ")"):
            continue
        out[j] = int(s)
        j += 1
    return out


def read_boundary(path: Path) -> list[tuple[str, str, int, int]]:
    txt = path.read_text(encoding="ascii")
    body = txt[txt.index("(") + 1:]
    patches: list[tuple[str, str, int, int]] = []
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", body):
        blk = m.group(2)
        patches.append((
            m.group(1),
            re.search(r"type\s+(\w+)", blk).group(1),
            int(re.search(r"nFaces\s+(\d+)", blk).group(1)),
            int(re.search(r"startFace\s+(\d+)", blk).group(1)),
        ))
    return patches


def _foam_header(class_name: str, object_name: str) -> str:
    return (
        "FoamFile\n{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        f"    class       {class_name};\n"
        '    location    "constant/polyMesh";\n'
        f"    object      {object_name};\n"
        "}\n\n"
    )


_CHUNK = 20000


def write_points(path: Path, pts: np.ndarray) -> None:
    with open(path, "w", encoding="ascii") as fh:
        fh.write(_foam_header("vectorField", "points"))
        fh.write(f"{pts.shape[0]}\n(\n")
        for start in range(0, pts.shape[0], _CHUNK):
            chunk = pts[start:start + _CHUNK]
            fh.write("\n".join("(%.12g %.12g %.12g)" % (x, y, z)
                               for x, y, z in chunk.tolist()))
            fh.write("\n")
        fh.write(")\n")


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def face_geometry(pos: np.ndarray, faces: np.ndarray):
    """Return (face vertices, face centres, area normals)."""
    fv = pos[faces]
    fc = fv.mean(axis=1)
    nxt = np.roll(fv, -1, axis=1)
    nx = np.sum((fv[:, :, 1] - nxt[:, :, 1]) * (fv[:, :, 2] + nxt[:, :, 2]), axis=1)
    ny = np.sum((fv[:, :, 2] - nxt[:, :, 2]) * (fv[:, :, 0] + nxt[:, :, 0]), axis=1)
    nz = np.sum((fv[:, :, 0] - nxt[:, :, 0]) * (fv[:, :, 1] + nxt[:, :, 1]), axis=1)
    return fv, fc, 0.5 * np.stack([nx, ny, nz], axis=1)


def cell_metrics(pos, faces, owner, neigh, ncells):
    """Volume / min-edge / aspect-ratio using OpenFOAM's own definitions.

    Aspect ratio is ``primitiveMeshTools::cellAspectRatio`` (v2412), so the
    numbers are directly comparable with ``checkMesh``.
    """
    fv, fc, an = face_geometry(pos, faces)
    farea = np.linalg.norm(an, axis=1)
    csum = np.zeros((ncells, 3))
    wsum = np.zeros(ncells)
    np.add.at(csum, owner, fc * farea[:, None])
    np.add.at(wsum, owner, farea)
    np.add.at(csum, neigh, fc[:neigh.size] * farea[:neigh.size, None])
    np.add.at(wsum, neigh, farea[:neigh.size])
    ctr0 = csum / wsum[:, None]
    # OpenFOAM cell::centre: volume-weighted average of the cell-centre-to-face
    # pyramids (see meshes/meshShapes/cell/cell.C).  Required for a
    # checkMesh-comparable non-orthogonality angle.
    fcells = np.concatenate([owner, neigh])
    fids = np.concatenate([np.arange(faces.shape[0]), np.arange(neigh.size)])
    fsign = np.concatenate([np.ones(faces.shape[0]), -np.ones(neigh.size)])
    pyr_vol = np.abs((1.0 / 3.0) * np.einsum(
        "ij,ij->i", an[fids] * fsign[:, None], ctr0[fcells] - fc[fids]))
    pyr_ctr = 0.75 * fc[fids] + 0.25 * ctr0[fcells]
    sum_v = np.bincount(fcells, weights=pyr_vol, minlength=ncells)
    sum_vc = np.stack(
        [np.bincount(fcells, weights=pyr_vol * pyr_ctr[:, d], minlength=ncells)
         for d in range(3)], axis=1)
    ccen = sum_vc / np.maximum(sum_v, 1e-300)[:, None]
    vol = np.zeros(ncells)
    np.add.at(vol, owner, (1.0 / 3.0) * np.einsum("ij,ij->i", an, fc - ccen[owner]))
    np.add.at(vol, neigh, (1.0 / 3.0) * np.einsum(
        "ij,ij->i", -an[:neigh.size], fc[:neigh.size] - ccen[neigh]))
    e = np.stack([np.linalg.norm(fv[:, k] - fv[:, (k + 1) % 4], axis=1)
                  for k in range(4)], axis=1)
    fmin = e.min(axis=1)
    fmax = e.max(axis=1)
    cmin = np.full(ncells, np.inf)
    cmax = np.zeros(ncells)
    np.minimum.at(cmin, owner, fmin)
    np.maximum.at(cmax, owner, fmax)
    np.minimum.at(cmin, neigh, fmin[:neigh.size])
    np.maximum.at(cmax, neigh, fmax[:neigh.size])
    summag = np.zeros((ncells, 3))
    np.add.at(summag, owner, np.abs(an))
    np.add.at(summag, neigh, np.abs(an[:neigh.size]))
    v = np.maximum(vol, 1e-300)
    ar = np.maximum(summag.max(axis=1) / (summag.min(axis=1) + 1e-300),
                    (1.0 / 6.0) * summag.sum(axis=1) / v ** (2.0 / 3.0))
    ni = neigh.size
    d = ccen[neigh] - ccen[owner[:ni]]
    anorm = np.linalg.norm(an[:ni], axis=1)
    cosang = np.einsum("ij,ij->i", an[:ni], d) / (
        anorm * np.linalg.norm(d, axis=1) + 1e-300)
    nonorth = float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))).max())
    return vol, cmin, ar, ccen, nonorth


def build_adjacency(faces: np.ndarray, npoints: int):
    """CSR adjacency from the unique undirected face edges."""
    f = faces
    e = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 3]], f[:, [3, 0]]])
    a = np.minimum(e[:, 0], e[:, 1]).astype(np.int64)
    b = np.maximum(e[:, 0], e[:, 1]).astype(np.int64)
    key = np.unique(a * np.int64(npoints) + b)
    a = (key // npoints).astype(np.int64)
    b = (key % npoints).astype(np.int64)
    src = np.concatenate([a, b])
    dst = np.concatenate([b, a])
    order = np.argsort(src, kind="stable")
    src = src[order]
    dst = dst[order]
    deg = np.bincount(src, minlength=npoints)
    starts = np.zeros(npoints + 1, dtype=np.int64)
    np.cumsum(deg, out=starts[1:])
    return src, dst, starts, deg


# ---------------------------------------------------------------------------
# Boundary classification
# ---------------------------------------------------------------------------


def patch_node_masks(faces, patches, npoints):
    """Return (node_patch_id, node_patch_name) with FSI first."""
    node_patch = np.full(npoints, -1, dtype=np.int64)
    names = [p[0] for p in patches]
    order = sorted(range(len(patches)), key=lambda i: 0 if names[i] == "FSI" else 1)
    for pi in order:
        name, typ, nf, sf = patches[pi]
        v = np.unique(faces[sf:sf + nf])
        node_patch[v[node_patch[v] < 0]] = pi
    return node_patch, names


def _plane_of_patch(name, typ):
    """Return (axis, value) for planar patches, else None."""
    if name == "FRONT":
        return (1, 0.0)
    if name == "BACK":
        return (1, 7.0)
    if name == "OUT":
        return (0, 6.0)
    if name == "TOP":
        return (2, 5.0)
    if name == "BOTTOM":
        return (2, -5.0)
    return None


# ---------------------------------------------------------------------------
# Smoother
# ---------------------------------------------------------------------------


def _neighbor_stats(cur, src, dst, npoints, weights):
    nb = np.zeros((npoints, 3))
    for d in range(3):
        nb[:, d] = np.bincount(src, weights=weights * cur[dst, d], minlength=npoints)
    w = np.bincount(src, weights=weights, minlength=npoints)
    return nb, w


def _sliver_score(vol, ar, vol_thresh, ar_thresh):
    """Number of degenerate cells (small volume or high aspect ratio)."""
    return int(np.count_nonzero(vol < vol_thresh) + np.count_nonzero(ar > ar_thresh))


def smooth(
    pos, faces, owner, neigh, patches, ncells,
    *, iters, relax, mode, target_edge, repel, region_mask,
    free_boundary, slide_boundary, ar_slack, nonorth_slack, objective,
    report_every, log,
):
    npoints = pos.shape[0]
    node_patch, names = patch_node_masks(faces, patches, npoints)

    fixed = np.zeros(npoints, dtype=bool)
    # FSI is always fixed; other patches fixed unless --free-boundary.
    for pi, (name, typ, nf, sf) in enumerate(patches):
        if name == "FSI":
            fixed[np.unique(faces[sf:sf + nf])] = True
        elif not free_boundary:
            fixed[np.unique(faces[sf:sf + nf])] = True

    free = (~fixed) & region_mask
    fidx = np.nonzero(free)[0]

    # Node masks / constraints for sliding boundary nodes.
    slide_axis = np.full(npoints, -1, dtype=np.int64)
    slide_value = np.zeros(npoints)
    cyl_mask = np.zeros(npoints, dtype=bool)
    if free_boundary and slide_boundary:
        for pi, (name, typ, nf, sf) in enumerate(patches):
            v = np.unique(faces[sf:sf + nf])
            plane = _plane_of_patch(name, typ)
            if plane is not None:
                slide_axis[v] = plane[0]
                slide_value[v] = plane[1]
            elif name == "IN":
                cyl_mask[v] = True
    # IN half-cylinder centre/axis (x-z cylinder centred at origin, axis y).
    cyl_r = np.zeros(npoints)
    if cyl_mask.any():
        cyl_r[cyl_mask] = np.hypot(pos[cyl_mask, 0], pos[cyl_mask, 2])

    src, dst, starts, deg = build_adjacency(faces, npoints)

    vol0, cmin0, ar0, _, nonorth0 = cell_metrics(pos, faces, owner, neigh, ncells)
    log(f"initial: minvol {vol0.min():.6e}  minedge {cmin0.min():.6e}  "
        f"maxAR {ar0.max():.6e}  maxNonOrth {nonorth0:.4f} deg  "
        f"negvol {int(np.count_nonzero(vol0 <= 0))}")

    cur = pos.copy()
    ar_cur = float(ar0.max())
    no_cur = float(nonorth0)
    best = {
        "pos": cur.copy(), "minedge": float(cmin0.min()), "ar": float(ar0.max()),
        "nonorth": float(nonorth0),
        "vol": float(vol0.min()),
        "sliver": _sliver_score(vol0, ar0, 1e-12, 1e5),
        "degen_vol": int(np.count_nonzero(vol0 < 1e-12)),
        "degen_ar": int(np.count_nonzero(ar0 > 1e5)),
        "negvol": int(np.count_nonzero(vol0 <= 0.0)),
        "iter": -1,
    }
    lam = relax
    mu = -relax * 1.05

    mean_len = np.bincount(src, weights=np.linalg.norm(
        cur[dst] - cur[np.repeat(np.arange(npoints), deg)], axis=1),
        minlength=npoints) / np.maximum(deg, 1)
    step_cap = 0.5 * mean_len
    stopped = False

    for it in range(iters):
        dvec = cur[dst] - cur[np.repeat(np.arange(npoints), deg)]
        dlen = np.linalg.norm(dvec, axis=1) + 1e-300
        w = np.ones_like(dlen)
        nb, wsum = _neighbor_stats(cur, src, dst, npoints, w)
        avg = np.zeros_like(cur)
        nz = wsum > 0
        avg[nz] = nb[nz] / wsum[nz, None]
        avg[~nz] = cur[~nz]

        if mode == "taubin":
            this_lam = mu if (it % 2 == 1) else lam
        else:
            this_lam = lam

        delta = this_lam * (avg - cur)

        if mode == "expand":
            # short-edge repulsion: push away from neighbours closer than target
            nprev = np.repeat(np.arange(npoints), deg)
            shorter = dlen < target_edge
            if shorter.any():
                fr = np.zeros_like(cur)
                direction = dvec[shorter] / dlen[shorter, None]
                mag = (target_edge - dlen[shorter]) * repel
                np.add.at(fr, nprev[shorter], direction * mag[:, None])
                dg = np.maximum(deg, 1)[:, None]
                np.divide(fr, dg, out=fr, where=dg > 0)
                delta = delta + fr

        disp = np.linalg.norm(delta, axis=1)
        over = disp > step_cap
        if over.any():
            delta[over] *= (step_cap[over] / disp[over])[:, None]

        step = np.zeros_like(cur)
        step[fidx] = delta[fidx]

        accepted = None
        for shrink in (1.0, 0.5, 0.25, 0.1, 0.05, 0.01):
            trial = cur.copy()
            trial[fidx] = cur[fidx] + shrink * step[fidx]
            # keep released boundary nodes on their planes / cylinder
            if free_boundary and slide_boundary:
                sl = slide_axis >= 0
                if sl.any():
                    trial[sl, slide_axis[sl]] = slide_value[sl]
                if cyl_mask.any():
                    xz = trial[cyl_mask][:, [0, 2]]
                    rr = np.hypot(xz[:, 0], xz[:, 1]) + 1e-300
                    scale = cyl_r[cyl_mask] / rr
                    trial[cyl_mask, 0] = xz[:, 0] * scale
                    trial[cyl_mask, 2] = xz[:, 1] * scale
            vol, cmin, ar, _, nonorth = cell_metrics(trial, faces, owner, neigh, ncells)
            # GUARD: positive volumes, no worsening of the worst aspect ratio
            # or of the maximum face non-orthogonality angle.
            if (vol.min() > 0.0
                    and ar.max() <= ar_cur * (1.0 + ar_slack) + 1e-30
                    and nonorth <= no_cur * (1.0 + nonorth_slack) + 1e-9):
                accepted = (trial, vol, cmin, ar, nonorth)
                break
        if accepted is None:
            stopped = True
            log(f"iter {it}: no acceptable step (volume/AR/non-orth guard); stopping")
            break
        cur, vol, cmin, ar, nonorth = accepted
        ar_cur = float(ar.max())
        no_cur = float(nonorth)

        # update local step cap from the current mesh
        step_cap = 0.5 * (np.bincount(src, weights=dlen, minlength=npoints)
                          / np.maximum(deg, 1))

        sliver = _sliver_score(vol, ar, 1e-12, 1e5)
        if objective == "minedge":
            better = cmin.min() > best["minedge"]
        elif objective == "maxar":
            better = ar.max() < best["ar"] - 1e-12
        else:
            better = (sliver < best["sliver"]) or (
                sliver == best["sliver"] and cmin.min() > best["minedge"])
        if better:
            best.update(pos=cur.copy(), minedge=float(cmin.min()),
                        ar=float(ar.max()), nonorth=float(nonorth),
                        vol=float(vol.min()), sliver=sliver,
                        degen_vol=int(np.count_nonzero(vol < 1e-12)),
                        degen_ar=int(np.count_nonzero(ar > 1e5)),
                        negvol=int(np.count_nonzero(vol <= 0.0)), iter=it)

        if report_every and (it % report_every == 0 or it == iters - 1):
            log(f"iter {it:4d} [{mode}]: minvol {vol.min():.6e}  "
                f"minedge {cmin.min():.6e}  maxAR {ar.max():.6e}  "
                f"maxNonOrth {nonorth:.3f}  slivers {sliver}  "
                f"(best {best['sliver']} @ {best['iter']})")

    log(f"best state: iter {best['iter']}  minvol {best['vol']:.6e}  "
        f"minedge {best['minedge']:.6e}  maxAR {best['ar']:.6e}  "
        f"maxNonOrth {best['nonorth']:.4f}  "
        f"slivers {best['sliver']} (vol<1e-12 {best['degen_vol']}, "
        f"AR>1e5 {best['degen_ar']})  stopped_early={stopped}")
    return best["pos"], {
        "initial_minvol": float(vol0.min()),
        "initial_minedge": float(cmin0.min()),
        "initial_maxar": float(ar0.max()),
        "initial_maxnonorth": float(nonorth0),
        "initial_slivers": _sliver_score(vol0, ar0, 1e-12, 1e5),
        "best_minvol": best["vol"],
        "best_minedge": best["minedge"],
        "best_maxar": best["ar"],
        "best_maxnonorth": best["nonorth"],
        "best_slivers": best["sliver"],
        "best_degen_vol": best["degen_vol"],
        "best_degen_ar": best["degen_ar"],
        "best_negvol": best["negvol"],
        "best_iter": best["iter"],
        "n_free": int(fidx.size),
        "n_fixed": int(fixed.sum()),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run(args) -> int:
    t0 = time.time()
    inp = Path(args.input)
    out = Path(args.out_dir)
    for name in ("points", "faces", "owner", "neighbour", "boundary"):
        if not (inp / name).is_file():
            raise SystemExit(f"[smooth] ERROR: missing {inp / name}")

    log_lines: list[str] = []

    def log(msg: str) -> None:
        line = f"[smooth] {msg}"
        print(line, flush=True)
        log_lines.append(line)

    pos = read_points(inp / "points")
    faces = read_faces(inp / "faces")
    owner = read_labels(inp / "owner")
    neigh = read_labels(inp / "neighbour")
    patches = read_boundary(inp / "boundary")
    ncells = int(owner.max()) + 1
    npoints = pos.shape[0]
    log(f"input: {npoints} points  {ncells} cells  {faces.shape[0]} faces  "
        f"{len(patches)} patches")

    # region mask
    region_mask = np.ones(npoints, dtype=bool)
    if args.region == "bbox":
        x0, x1, y0, y1, z0, z1 = args.bbox
        region_mask = ((pos[:, 0] >= x0) & (pos[:, 0] <= x1) &
                       (pos[:, 1] >= y0) & (pos[:, 1] <= y1) &
                       (pos[:, 2] >= z0) & (pos[:, 2] <= z1))
        log(f"region bbox x[{x0},{x1}] y[{y0},{y1}] z[{z0},{z1}]: "
            f"{int(region_mask.sum())} nodes in region")

    best_pos, stats = smooth(
        pos, faces, owner, neigh, patches, ncells,
        iters=args.iters, relax=args.relax, mode=args.mode,
        target_edge=args.target_edge, repel=args.repel,
        region_mask=region_mask, free_boundary=args.free_boundary,
        slide_boundary=args.slide_boundary, ar_slack=args.ar_slack,
        nonorth_slack=args.nonorth_slack, objective=args.objective,
        report_every=args.report_every, log=log,
    )

    # sanity: topology unchanged
    assert best_pos.shape == pos.shape
    max_move = float(np.linalg.norm(best_pos - pos, axis=1).max())
    log(f"max node displacement: {max_move:.6e} m  "
        f"({int(np.count_nonzero(np.any(best_pos != pos, axis=1)))} nodes moved)")

    mesh_dir = out / "constant" / "polyMesh"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    write_points(mesh_dir / "points", best_pos)
    for name in ("faces", "owner", "neighbour", "boundary"):
        shutil.copyfile(inp / name, mesh_dir / name)
    log(f"wrote {mesh_dir}")

    if args.report:
        rp = Path(args.report)
        rp.parent.mkdir(parents=True, exist_ok=True)
        with open(rp, "w", encoding="ascii") as fh:
            fh.write("mode=%s iters=%d relax=%g target_edge=%g repel=%g "
                     "region=%s free_boundary=%s slide_boundary=%s\n" % (
                         args.mode, args.iters, args.relax, args.target_edge,
                         args.repel, args.region, args.free_boundary,
                         args.slide_boundary))
            for k, v in stats.items():
                fh.write(f"{k}={v}\n")
            fh.write(f"max_node_displacement={max_move:.8e}\n")
            for line in log_lines:
                fh.write(line + "\n")
        log(f"report: {rp}")

    log(f"done in {time.time() - t0:.0f}s")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="smooth_farfield_polymesh.py",
        description="Topology-preserving node smoothing for an OpenFOAM polyMesh. "
                    "Moves point coordinates only; topology/patch definitions are copied verbatim.")
    p.add_argument("--input", required=True, type=Path,
                   help="Input polyMesh directory (containing points/faces/owner/neighbour/boundary)")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Output case directory; writes constant/polyMesh/*")
    p.add_argument("--iters", type=int, default=100, help="Number of smoothing iterations")
    p.add_argument("--relax", type=float, default=0.3, help="Laplacian relaxation (0<relax<1)")
    p.add_argument("--mode", choices=["laplacian", "taubin", "expand"], default="laplacian",
                   help="Smoothing kernel (default laplacian)")
    p.add_argument("--target-edge", type=float, default=1e-4,
                   help="[expand] desired minimum edge length (m)")
    p.add_argument("--repel", type=float, default=0.1,
                   help="[expand] repulsion gain")
    p.add_argument("--region", choices=["all", "bbox"], default="all",
                   help="Which free nodes may move (default all)")
    p.add_argument("--bbox", type=float, nargs=6, metavar=("X0", "X1", "Y0", "Y1", "Z0", "Z1"),
                   default=None, help="Region bounding box for --region bbox")
    p.add_argument("--free-boundary", action="store_true",
                   help="Also release non-FSI boundary-patch nodes (FSI stays fixed)")
    p.add_argument("--slide-boundary", action="store_true",
                   help="With --free-boundary: constrain released boundary nodes to their plane/cylinder")
    p.add_argument("--ar-slack", type=float, default=0.0,
                   help="Allowed fractional worsening of the worst aspect ratio per step")
    p.add_argument("--nonorth-slack", type=float, default=0.0,
                   help="Allowed fractional worsening of max face non-orthogonality (deg)")
    p.add_argument("--objective", choices=["sliver", "minedge", "maxar"], default="sliver",
                   help="Which accepted state to write: fewest degenerate cells (default), "
                        "largest minimum edge, or smallest worst aspect ratio")
    p.add_argument("--report-every", type=int, default=20, help="Log cadence (0 = final only)")
    p.add_argument("--report", type=Path, default=None, help="Write a text report to this path")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not (0.0 < args.relax < 1.0):
        raise SystemExit("[smooth] ERROR: --relax must be in (0,1)")
    if args.region == "bbox" and args.bbox is None:
        raise SystemExit("[smooth] ERROR: --region bbox requires --bbox X0 X1 Y0 Y1 Z0 Z1")
    if args.slide_boundary and not args.free_boundary:
        raise SystemExit("[smooth] ERROR: --slide-boundary requires --free-boundary")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
