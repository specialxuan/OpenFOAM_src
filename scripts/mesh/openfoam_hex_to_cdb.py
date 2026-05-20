#!/usr/bin/env python3
"""
Convert OpenFOAM snappyHexMesh polyMesh to ANSYS CDB for CalculiX modal analysis.

Uses the proven order_hex_nodes() from damfailure_mesh_independence for C3D8
ordering, and skips polyhedral cells (all their nodes are covered by hex cells).
"""

import argparse, re, math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


def parse_points(path: Path) -> List[Tuple[float, float, float]]:
    content = path.read_text()
    m = re.search(r"\(\s*\n(.*?)\n\s*\)", content, re.DOTALL)
    pts = []
    for line in m.group(1).strip().split("\n"):
        line = line.strip().strip("()")
        parts = line.split()
        if len(parts) >= 3:
            pts.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return pts


def parse_faces(path: Path) -> List[List[int]]:
    lines = path.read_text().strip().split("\n")
    ds = 0
    for i, l in enumerate(lines):
        if l.strip() == "(": ds = i + 1; break
    faces = []
    for i in range(ds, len(lines)):
        l = lines[i].strip()
        if l.startswith(")"): break
        m = re.match(r"(\d+)\s*\(([^)]*)\)", l)
        if m: faces.append([int(x) for x in m.group(2).split()])
    return faces


def parse_ilist(path: Path) -> List[int]:
    c = path.read_text()
    m = re.search(r"\(\s*\n(.*?)\n\s*\)", c, re.DOTALL)
    vals = []
    if m:
        for line in m.group(1).strip().split("\n"):
            for v in re.findall(r"\d+", line): vals.append(int(v))
    return vals


def parse_boundary(path: Path) -> List[Tuple[str, int, int]]:
    text = path.read_text()
    m = re.search(r"\(\s*\n(.*?)\n\s*\)", text, re.DOTALL)
    patches = []
    for pm in re.finditer(r"(\w+)\s*\{(.+?)\}", m.group(1), re.DOTALL):
        nf = re.search(r"nFaces\s+(\d+)", pm.group(2))
        sf = re.search(r"startFace\s+(\d+)", pm.group(2))
        if nf and sf: patches.append((pm.group(1), int(nf.group(1)), int(sf.group(1))))
    return patches


# ---- Taken from damfailure_mesh_independence/openfoam_mesh_to_cdb.py ----

def order_hex_nodes(node_ids: Sequence[int], points: Sequence[Tuple[float, float, float]]) -> List[int]:
    """Geometric C3D8 corner ordering (proven in damfailure workflow)."""
    if len(node_ids) != 8:
        raise ValueError(f"Expected 8-node hex, got {len(node_ids)}")
    coords = [points[n] for n in node_ids]
    xs = [c[0] for c in coords]; ys = [c[1] for c in coords]; zs = [c[2] for c in coords]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin, zmax = min(zs), max(zs)
    dx = max(xmax - xmin, 1.0e-16)
    dy = max(ymax - ymin, 1.0e-16)
    dz = max(zmax - zmin, 1.0e-16)
    norm = {nid: ((points[nid][0]-xmin)/dx, (points[nid][1]-ymin)/dy, (points[nid][2]-zmin)/dz)
            for nid in node_ids}
    corners = [
        (0.0,0.0,0.0), (1.0,0.0,0.0), (1.0,1.0,0.0), (0.0,1.0,0.0),
        (0.0,0.0,1.0), (1.0,0.0,1.0), (1.0,1.0,1.0), (0.0,1.0,1.0),
    ]
    unused = set(node_ids)
    ordered = []
    for cx, cy, cz in corners:
        best = None; best_d2 = 1e300
        for nid in unused:
            nx, ny, nz = norm[nid]
            d2 = (nx-cx)**2 + (ny-cy)**2 + (nz-cz)**2
            if d2 < best_d2: best_d2 = d2; best = nid
        if best is None: raise ValueError("Failed to order hex nodes")
        ordered.append(best); unused.remove(best)
    return ordered


def hex_volume(coords: Sequence[Tuple[float, float, float]]) -> float:
    if len(coords) != 8: return 0.0
    def tv(a, b, c, d):
        ad = (a[0]-d[0], a[1]-d[1], a[2]-d[2])
        bd = (b[0]-d[0], b[1]-d[1], b[2]-d[2])
        cd = (c[0]-d[0], c[1]-d[1], c[2]-d[2])
        cr = (bd[1]*cd[2]-bd[2]*cd[1], bd[2]*cd[0]-bd[0]*cd[2], bd[0]*cd[1]-bd[1]*cd[0])
        return abs(ad[0]*cr[0] + ad[1]*cr[1] + ad[2]*cr[2]) / 6.0
    return (tv(coords[0],coords[1],coords[3],coords[4]) +
            tv(coords[1],coords[2],coords[3],coords[6]) +
            tv(coords[1],coords[3],coords[4],coords[6]) +
            tv(coords[1],coords[4],coords[5],coords[6]) +
            tv(coords[3],coords[4],coords[6],coords[7]))


def write_cmblock(fp, name: str, ids: Sequence[int]):
    fp.write(f"CMBLOCK,{name},NODE,{len(ids)}\n(8i10)\n")
    for i in range(0, len(ids), 8):
        fp.write("".join(f"{v:10d}" for v in ids[i:i+8]) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--out-cdb", required=True)
    parser.add_argument("--fsi-patch", default="AGARD_WING")
    args = parser.parse_args()

    case = Path(args.case)
    poly = case / "constant" / "polyMesh"

    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_ilist(poly / "owner")
    neighbour = parse_ilist(poly / "neighbour")
    patches = parse_boundary(poly / "boundary")

    n_cells = max(max(owner), max(neighbour)) + 1
    cell_faces = [[] for _ in range(n_cells)]
    for fi, c in enumerate(owner): cell_faces[c].append(fi)
    for fi, c in enumerate(neighbour): cell_faces[c].append(fi)

    # Collect hex cells only - must have exactly 6 faces, all 4-node quads, good quality
    elements = []
    volumes = []
    skipped_quality = 0
    for cid in range(n_cells):
        node_set = set()
        n_faces = len(cell_faces[cid])
        all_quads = True
        for fi in cell_faces[cid]:
            node_set.update(faces[fi])
            if len(faces[fi]) != 4:
                all_quads = False
        if len(node_set) != 8 or n_faces != 6 or not all_quads:
            continue
        ordered = order_hex_nodes(sorted(node_set), points)
        conn = [n + 1 for n in ordered]
        coords = [points[n] for n in ordered]
        vol = hex_volume(coords)
        # Quality check: compute bounding box volume ratio
        xs = [c[0] for c in coords]; ys = [c[1] for c in coords]; zs = [c[2] for c in coords]
        bb_vol = (max(xs)-min(xs)) * (max(ys)-min(ys)) * (max(zs)-min(zs))
        if bb_vol > 0 and vol / bb_vol < 1e-6:  # extremely flat cell
            skipped_quality += 1
            continue
        elements.append(conn)
        volumes.append(vol)

    # Patch node sets
    patch_nodes = {"FLUIDDOMAIN": sorted(set().union(*(set(e) for e in elements)))}
    tie_master_nodes = set()
    tie_slave_nodes = set()

    # Compute cell volumes for coarse/fine classification
    cell_vols = {}
    for cid in range(n_cells):
        ns = set()
        for fi in cell_faces[cid]: ns.update(faces[fi])
        if len(ns) == 8:
            xs = [points[n][0] for n in ns]
            ys = [points[n][1] for n in ns]
            zs = [points[n][2] for n in ns]
            cell_vols[cid] = (max(xs)-min(xs))*(max(ys)-min(ys))*(max(zs)-min(zs))
    median_vol = sorted(cell_vols.values())[len(cell_vols)//2] if cell_vols else 1.0

    for name, nf, sf in patches:
        ids = set()
        for fi in range(sf, sf + nf):
            if fi < len(faces):
                for nid in faces[fi]: ids.add(nid + 1)
        upper = name.upper()
        if upper == "OLDINTERNALFACES":
            # Split into TIE_MASTER (coarse) and TIE_SLAVE (fine) based on owner cell size
            for fi in range(sf, sf + nf):
                if fi < len(owner):
                    c = owner[fi]
                    if cell_vols.get(c, 0) > median_vol:
                        for nid in faces[fi]: tie_master_nodes.add(nid + 1)
                    else:
                        for nid in faces[fi]: tie_slave_nodes.add(nid + 1)
        else:
            patch_nodes[upper] = sorted(ids)

    # Map patch names
    fsi_name = args.fsi_patch.upper()
    patch_nodes["FSI_FLUID"] = patch_nodes.pop(fsi_name, [])
    mapping = {"INLET":"LEFT","OUTLET":"RIGHT","FARFIELD":"TOP","ROOTSYMMETRY":"FRONT"}
    for old, new in mapping.items():
        if old in patch_nodes and new not in patch_nodes:
            patch_nodes[new] = patch_nodes.pop(old)
    for d in ["FRONT","BACK","TOP","BOTTOM_FLUID","LEFT","RIGHT"]:
        patch_nodes.setdefault(d, [])

    # Add TIE node sets
    patch_nodes["TIE_MASTER"] = sorted(tie_master_nodes) if tie_master_nodes else []
    patch_nodes["TIE_SLAVE"] = sorted(tie_slave_nodes) if tie_slave_nodes else []

    # Write CDB
    out = Path(args.out_cdb)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write("/BATCH\n/NOPR\n")
        f.write(f"/COM, AGARD Hex Fluid CDB: {len(points)} pts, {len(elements)} elems\n")
        f.write("/PREP7\nET,1,185,,0,,0,,0,0\n")
        f.write("MP,EX,  1,           1\nMP,PRXY,  1,         0.5\nMP,DENS,  1,      0.0025\n")
        f.write("/NOLIST\nNBLOCK,6,SOLID\n(3i8,6e16.9)\n")
        for nid, (x, y, z) in enumerate(points, start=1):
            f.write(f"{nid:8d}{0:8d}{0:8d}{x:16.9E}{y:16.9E}{z:16.9E}\n")
        f.write("N ,R5.3,LOC,     -1\nSHPP,WARN\nEBLOCK,19,SOLID\n(19i10)\n")
        for eid, conn in enumerate(elements, start=1):
            vals = [1,1,0,0,0,0,0,0,8,0,eid] + conn
            f.write("".join(f"{v:10d}" for v in vals) + "\n")
        f.write("-1\nEN ,R5.5,ATTR,     -1\n")
        for name in sorted(patch_nodes):
            if patch_nodes[name]:
                write_cmblock(f, name, patch_nodes[name])

    # Write volumes
    vpath = out.with_suffix(".volumes.txt")
    with vpath.open("w") as vf:
        for eid, vol in enumerate(volumes, start=1):
            vf.write(f"{eid} {vol:.15e}\n")

    print(f"[OK] {out}")
    print(f"  Nodes: {len(points)}, Elements: {len(elements)} (hex only, poly skipped)")
    print(f"  TIE_MASTER: {len(tie_master_nodes)} nodes, TIE_SLAVE: {len(tie_slave_nodes)} nodes")
    for name in sorted(patch_nodes):
        n = len(patch_nodes[name])
        if n > 0: print(f"    {name}: {n} nodes")
    print(f"  Volumes: {vpath}")


if __name__ == "__main__":
    main()
