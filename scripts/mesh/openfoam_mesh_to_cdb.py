#!/usr/bin/env python3
"""
Convert OpenFOAM polyMesh (HEX mesh from blockMesh) into a minimal ANSYS CDB-like
input file that can be consumed by build_calculix_modal_cases.py.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


INT_RE = re.compile(r"[+-]?\d+")
FLOAT_RE = re.compile(r"[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[EeDd][+-]?\d+)?")
PATCH_RE = re.compile(
    r"^\s*([A-Za-z0-9_]+)\s*\{[^{}]*?nFaces\s+(\d+)\s*;[^{}]*?startFace\s+(\d+)\s*;",
    re.MULTILINE | re.DOTALL,
)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def parse_points(path: Path) -> List[Tuple[float, float, float]]:
    lines = read_text(path).splitlines()
    count = None
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if count is None and s.isdigit():
            count = int(s)
            i += 1
            continue
        if count is not None and s.startswith("("):
            i += 1
            break
        i += 1
    if count is None:
        raise ValueError(f"Cannot parse points count from {path}")

    pts: List[Tuple[float, float, float]] = []
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith(")"):
            break
        m = re.match(r"^\(\s*([^\s]+)\s+([^\s]+)\s+([^\s]+)\s*\)$", s)
        if m:
            x = float(m.group(1).replace("D", "E").replace("d", "e"))
            y = float(m.group(2).replace("D", "E").replace("d", "e"))
            z = float(m.group(3).replace("D", "E").replace("d", "e"))
            pts.append((x, y, z))
        i += 1

    if len(pts) != count:
        raise ValueError(f"{path}: expected {count} points, got {len(pts)}")
    return pts


def parse_label_list(path: Path) -> List[int]:
    lines = read_text(path).splitlines()
    count = None
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if count is None and s.isdigit():
            count = int(s)
            i += 1
            continue
        if count is not None and s.startswith("("):
            i += 1
            break
        i += 1
    if count is None:
        raise ValueError(f"Cannot parse list count from {path}")

    vals: List[int] = []
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith(")"):
            break
        vals.extend(int(v) for v in INT_RE.findall(s))
        i += 1

    if len(vals) != count:
        raise ValueError(f"{path}: expected {count} values, got {len(vals)}")
    return vals


def parse_faces(path: Path) -> List[List[int]]:
    lines = read_text(path).splitlines()
    count = None
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if count is None and s.isdigit():
            count = int(s)
            i += 1
            continue
        if count is not None and s.startswith("("):
            i += 1
            break
        i += 1
    if count is None:
        raise ValueError(f"Cannot parse faces count from {path}")

    faces: List[List[int]] = []
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith(")"):
            break
        m = re.match(r"^\s*(\d+)\s*\(([^)]*)\)\s*$", s)
        if m:
            n_pts = int(m.group(1))
            node_ids = [int(v) for v in INT_RE.findall(m.group(2))]
            if len(node_ids) != n_pts:
                raise ValueError(f"{path}: malformed face line: {s}")
            faces.append(node_ids)
        i += 1

    if len(faces) != count:
        raise ValueError(f"{path}: expected {count} faces, got {len(faces)}")
    return faces


def parse_boundary(path: Path) -> List[Tuple[str, int, int]]:
    text = read_text(path)
    patches = []
    for name, nfaces, start in PATCH_RE.findall(text):
        patches.append((name.upper(), int(nfaces), int(start)))
    if not patches:
        raise ValueError(f"Cannot parse boundary patches from {path}")
    return patches


def order_hex_nodes(node_ids: Sequence[int], points: Sequence[Tuple[float, float, float]]) -> List[int]:
    if len(node_ids) != 8:
        raise ValueError(f"Expected 8-node hex, got {len(node_ids)}")

    coords = [points[n] for n in node_ids]
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin, zmax = min(zs), max(zs)

    dx = max(xmax - xmin, 1.0e-16)
    dy = max(ymax - ymin, 1.0e-16)
    dz = max(zmax - zmin, 1.0e-16)

    norm = {
        nid: ((points[nid][0] - xmin) / dx, (points[nid][1] - ymin) / dy, (points[nid][2] - zmin) / dz)
        for nid in node_ids
    }

    corners = [
        (0.0, 0.0, 0.0),  # n1
        (1.0, 0.0, 0.0),  # n2
        (1.0, 1.0, 0.0),  # n3
        (0.0, 1.0, 0.0),  # n4
        (0.0, 0.0, 1.0),  # n5
        (1.0, 0.0, 1.0),  # n6
        (1.0, 1.0, 1.0),  # n7
        (0.0, 1.0, 1.0),  # n8
    ]

    unused = set(node_ids)
    ordered: List[int] = []
    for cx, cy, cz in corners:
        best = None
        best_d2 = 1.0e300
        for nid in unused:
            nx, ny, nz = norm[nid]
            d2 = (nx - cx) ** 2 + (ny - cy) ** 2 + (nz - cz) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best = nid
        if best is None:
            raise ValueError("Failed to order hex nodes")
        ordered.append(best)
        unused.remove(best)
    return ordered


def write_cmblock(fp, name: str, ids: Sequence[int]) -> None:
    fp.write(f"CMBLOCK,{name},NODE,{len(ids)}\n")
    fp.write("(8i10)\n")
    for i in range(0, len(ids), 8):
        chunk = ids[i : i + 8]
        fp.write("".join(f"{v:10d}" for v in chunk) + "\n")


def convert_case(case_dir: Path, out_file: Path) -> None:
    poly = case_dir / "constant" / "polyMesh"
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_label_list(poly / "owner")
    neighbour = parse_label_list(poly / "neighbour")
    patches = parse_boundary(poly / "boundary")

    if len(owner) != len(faces):
        raise ValueError(f"{case_dir}: owner length != faces length")
    if len(neighbour) > len(faces):
        raise ValueError(f"{case_dir}: neighbour length > faces length")

    n_cells = max(max(owner), max(neighbour) if neighbour else 0) + 1
    cell_faces: List[List[int]] = [[] for _ in range(n_cells)]
    for fi, c in enumerate(owner):
        cell_faces[c].append(fi)
    for fi, c in enumerate(neighbour):
        cell_faces[c].append(fi)

    cell_conn: List[List[int]] = []
    fluid_nodes = set()
    for cid in range(n_cells):
        node_set = set()
        for fi in cell_faces[cid]:
            node_set.update(faces[fi])
        if len(node_set) != 8:
            raise ValueError(f"{case_dir}: cell {cid} has {len(node_set)} unique nodes (expect 8)")
        ordered = order_hex_nodes(sorted(node_set), points)
        cell_conn.append([nid + 1 for nid in ordered])  # 1-based ids for CDB
        fluid_nodes.update(nid + 1 for nid in node_set)

    patch_nodes: Dict[str, List[int]] = {}
    for name, nfaces, start in patches:
        ids = set()
        end = start + nfaces
        for fi in range(start, end):
            for nid in faces[fi]:
                ids.add(nid + 1)
        patch_nodes[name] = sorted(ids)

    patch_nodes["FLUIDDOMAIN"] = sorted(fluid_nodes)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", encoding="utf-8") as fp:
        fp.write("/BATCH\n")
        fp.write("/NOPR\n")
        fp.write("/PREP7\n")
        fp.write("ET,1,185,,0,,0,,0,0\n")
        fp.write("MP,EX,  1,           1\n")
        fp.write("MP,PRXY,  1,         0.5\n")
        fp.write("MP,DENS,  1,      0.0025\n")
        fp.write("/NOLIST\n")
        fp.write("NBLOCK,6,SOLID\n")
        fp.write("(3i8,6e16.9)\n")
        for nid0, (x, y, z) in enumerate(points, start=1):
            fp.write(f"{nid0:8d}{0:8d}{0:8d}{x:16.9E}{y:16.9E}{z:16.9E}\n")
        fp.write("N ,R5.3,LOC,     -1\n")
        fp.write("SHPP,WARN\n")
        fp.write("EBLOCK,19,SOLID\n")
        fp.write("(19i10)\n")
        for eid, conn in enumerate(cell_conn, start=1):
            vals = [1, 1, 0, 0, 0, 0, 0, 0, 8, 0, eid, *conn]
            fp.write("".join(f"{v:10d}" for v in vals) + "\n")
        fp.write("-1\n")
        fp.write("EN ,R5.5,ATTR,     -1\n")

        for name in sorted(patch_nodes.keys()):
            write_cmblock(fp, name, patch_nodes[name])


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert OpenFOAM polyMesh case to CDB-like fluid input")
    parser.add_argument("--case", action="append", required=True, help="OpenFOAM case directory")
    parser.add_argument("--out-dir", required=True, help="Output directory for generated CDB files")
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    for case in args.case:
        cdir = Path(case).resolve()
        name = cdir.name
        out_file = out_dir / f"{name}_fluid.in"
        convert_case(cdir, out_file)
        print(f"[OK] {cdir} -> {out_file}")


if __name__ == "__main__":
    main()
