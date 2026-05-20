#!/usr/bin/env python3
"""
Convert OpenFOAM tetrahedral polyMesh to a minimal ANSYS CDB-like fluid input.

Handles tetra cells by writing degenerate C3D8 elements (4 unique nodes
with quad-triplication to make 8 connectivity entries — standard ANSYS
SOLID185 degenerate tetra).

Node sets written:
  FLUIDDOMAIN - all fluid nodes
  FSI_FLUID   - nodes on the FSI interface patch (e.g., AGARD_WING)
  FRONT, BACK, TOP, BOTTOM_FLUID, LEFT, RIGHT - boundary sets

Each boundary patch in the OpenFOAM case is mapped to one of these sets
based on patch name matching rules.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple


INT_RE = re.compile(r"[+-]?\d+")
FLOAT_RE = re.compile(r"[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[EeDd][+-]?\d+)?")


def parse_points(path: Path) -> List[Tuple[float, float, float]]:
    content = path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"\(\s*\n(.*?)\n\s*\)", content, re.DOTALL)
    if not m:
        raise ValueError(f"Cannot parse points from {path}")
    pts = []
    for line in m.group(1).strip().split("\n"):
        line = line.strip().strip("()")
        parts = line.split()
        if len(parts) >= 3:
            pts.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return pts


def parse_faces(path: Path) -> List[List[int]]:
    lines = path.read_text(encoding="utf-8", errors="ignore").strip().split("\n")
    # Find count line
    count = 0
    data_start = 0
    for i, line in enumerate(lines):
        s = line.strip()
        if s.isdigit() and count == 0:
            count = int(s)
        if s == "(":
            data_start = i + 1
            break
    if count == 0:
        raise ValueError(f"Cannot parse face count from {path}")

    faces = []
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line.startswith(")"):
            break
        m = re.match(r"(\d+)\s*\(([^)]*)\)", line)
        if m:
            nodes = [int(x) for x in m.group(2).split()]
            faces.append(nodes)
    if len(faces) != count:
        raise ValueError(f"{path}: expected {count} faces, got {len(faces)}")
    return faces


def parse_label_list(path: Path) -> List[int]:
    content = path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"\(\s*\n(.*?)\n\s*\)", content, re.DOTALL)
    if not m:
        raise ValueError(f"Cannot parse label list from {path}")
    vals = []
    for line in m.group(1).strip().split("\n"):
        for v in INT_RE.findall(line):
            vals.append(int(v))
    return vals


def parse_boundary(path: Path) -> List[Tuple[str, int, int]]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"\(\s*\n(.*?)\n\s*\)", text, re.DOTALL)
    if not m:
        raise ValueError(f"Cannot parse boundary from {path}")
    patch_text = m.group(1)
    patches = []
    patch_re = re.compile(r"(\w+)\s*\{(.*?)\}", re.DOTALL)
    for pm in patch_re.finditer(patch_text):
        name = pm.group(1)
        body = pm.group(2)
        nf_match = re.search(r"nFaces\s+(\d+)", body)
        sf_match = re.search(r"startFace\s+(\d+)", body)
        if nf_match and sf_match:
            patches.append((name, int(nf_match.group(1)), int(sf_match.group(1))))
    return patches


def order_tet_nodes_for_degenerate_hex(
    node_ids: Sequence[int],
    points: Sequence[Tuple[float, float, float]],
) -> Tuple[int, int, int, int, int, int, int, int]:
    """
    Convert a 4-node tetra to an 8-node degenerate C3D8 connectivity.
    Standard ANSYS degenerated tetra pattern (SOLID185):
      face 1->2->3, node 4 is the apex.
    Connectivity: n1, n2, n3, n3, n4, n4, n4, n4
    """
    if len(node_ids) != 4:
        raise ValueError(f"Expected 4-node tetra, got {len(node_ids)}")

    coords = [(n, points[n]) for n in node_ids]

    # Find the apex: the node with the largest distance to the centroid of the other 3
    cx = sum(p[0] for _, p in coords) / 4.0
    cy = sum(p[1] for _, p in coords) / 4.0
    cz = sum(p[2] for _, p in coords) / 4.0

    def dist_sq(n, p):
        return (p[0] - cx) ** 2 + (p[1] - cy) ** 2 + (p[2] - cz) ** 2

    sorted_nodes = sorted(coords, key=lambda x: dist_sq(x[0], x[1]), reverse=True)
    n4 = sorted_nodes[0][0]  # apex
    base = [n for n, _ in sorted_nodes[1:4]]  # base triangle

    # Degenerate C3D8 for tetra: n1,n2,n3,n3, n4,n4,n4,n4
    return (base[0], base[1], base[2], base[2], n4, n4, n4, n4)


def tet_volume(
    a: Tuple[float, float, float],
    b: Tuple[float, float, float],
    c: Tuple[float, float, float],
    d: Tuple[float, float, float],
) -> float:
    """Absolute volume of a tetrahedron."""
    ad = (a[0] - d[0], a[1] - d[1], a[2] - d[2])
    bd = (b[0] - d[0], b[1] - d[1], b[2] - d[2])
    cd = (c[0] - d[0], c[1] - d[1], c[2] - d[2])
    cross = (
        bd[1] * cd[2] - bd[2] * cd[1],
        bd[2] * cd[0] - bd[0] * cd[2],
        bd[0] * cd[1] - bd[1] * cd[0],
    )
    return abs(ad[0] * cross[0] + ad[1] * cross[1] + ad[2] * cross[2]) / 6.0


def write_cmblock(fp, name: str, ids: Sequence[int]) -> None:
    fp.write(f"CMBLOCK,{name},NODE,{len(ids)}\n")
    fp.write("(8i10)\n")
    for i in range(0, len(ids), 8):
        chunk = ids[i : i + 8]
        fp.write("".join(f"{v:10d}" for v in chunk) + "\n")


def map_patch_to_set(patch_name: str) -> str:
    """Map OpenFOAM patch names to CalculiX/ANSYS node set names."""
    upper = patch_name.upper()
    mapping = {
        "AGARD_WING": "FSI_FLUID",
        "ROOTSYMMETRY": "FRONT",  # root plane - constrain UZ
        "FARFIELD": "TOP",         # farfield - fully fixed
        "INLET": "LEFT",
        "OUTLET": "RIGHT",
    }
    return mapping.get(upper, patch_name.upper())


def convert_case(case_dir: Path, out_file: Path, fsi_patch: str) -> None:
    poly = case_dir / "constant" / "polyMesh"
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_label_list(poly / "owner")
    neighbour = parse_label_list(poly / "neighbour")
    patches = parse_boundary(poly / "boundary")

    n_cells = max(max(owner), max(neighbour) if neighbour else 0) + 1

    # Build cell->faces mapping
    cell_faces: List[List[int]] = [[] for _ in range(n_cells)]
    for fi, c in enumerate(owner):
        cell_faces[c].append(fi)
    for fi, c in enumerate(neighbour):
        cell_faces[c].append(fi)

    # Build cell connectivity and track volumes
    cell_conn: List[List[int]] = []
    cell_vols: List[float] = []
    fluid_nodes: Set[int] = set()

    for cid in range(n_cells):
        node_set = set()
        for fi in cell_faces[cid]:
            node_set.update(faces[fi])
        nodes_list = sorted(node_set)

        if len(nodes_list) == 4:
            # Tetra cell: convert to degenerate C3D8
            conn = order_tet_nodes_for_degenerate_hex(nodes_list, points)
            # Compute actual tetra volume
            a = points[nodes_list[0]]
            b = points[nodes_list[1]]
            c = points[nodes_list[2]]
            d = points[nodes_list[3]]
            vol = tet_volume(a, b, c, d)
        elif len(nodes_list) == 8:
            # Hex cell: standard C3D8
            conn = tuple(nodes_list)  # keep original ordering
            vol = 0.0  # placeholder
            print(f"  Warning: cell {cid} is hex (8 nodes), using standard ordering")
        elif len(nodes_list) == 6:
            print(f"  Warning: cell {cid} is prism (6 nodes), skipping")
            continue
        else:
            print(f"  Warning: cell {cid} has {len(nodes_list)} nodes, skipping")
            continue

        # Convert to 1-based IDs for CDB
        conn_1b = [n + 1 for n in conn]
        cell_conn.append(conn_1b)
        cell_vols.append(vol)
        fluid_nodes.update(n + 1 for n in nodes_list)

    # Patch node sets
    patch_nodes: Dict[str, List[int]] = {}
    for name, nfaces, start in patches:
        ids = set()
        end = start + nfaces
        for fi in range(start, end):
            if fi < len(faces):
                for nid in faces[fi]:
                    ids.add(nid + 1)
        set_name = map_patch_to_set(name)
        patch_nodes[set_name] = sorted(ids)

    patch_nodes["FLUIDDOMAIN"] = sorted(fluid_nodes)

    # Assign remaining standard sets
    for default_set in ["FRONT", "BACK", "TOP", "BOTTOM_FLUID", "LEFT", "RIGHT"]:
        if default_set not in patch_nodes:
            patch_nodes[default_set] = []

    # If we have rootSymmetry mapped to FRONT, also set BOTTOM_FLUID = FRONT
    # for this case (symmetry plane is also a fixing plane)
    if "FRONT" in patch_nodes and patch_nodes["FRONT"] and not patch_nodes.get("BOTTOM_FLUID"):
        pass  # keep default empty

    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", encoding="utf-8") as f:
        f.write("/BATCH\n")
        f.write("/NOPR\n")
        f.write(f"/COM, AGARD 445.6 Fluid CDB - tetra mesh ({len(cell_conn)} elements)\n")
        f.write("/PREP7\n")
        f.write("ET,1,185,,0,,0,,0,0\n")
        f.write("MP,EX,  1,           1\n")
        f.write("MP,PRXY,  1,         0.5\n")
        f.write("MP,DENS,  1,      0.0025\n")
        f.write("/NOLIST\n")
        f.write("NBLOCK,6,SOLID\n")
        f.write("(3i8,6e16.9)\n")
        for nid0, (x, y, z) in enumerate(points, start=1):
            f.write(f"{nid0:8d}{0:8d}{0:8d}{x:16.9E}{y:16.9E}{z:16.9E}\n")
        f.write("N ,R5.3,LOC,     -1\n")
        f.write("SHPP,WARN\n")
        f.write("EBLOCK,19,SOLID\n")
        f.write("(19i10)\n")
        for eid, (conn, vol) in enumerate(zip(cell_conn, cell_vols), start=1):
            # EBLOCK format: mat=1, type=1, ..., num_nodes=4 (tet via degenerate hex), eid, conn
            vals = [1, 1, 0, 0, 0, 0, 0, 0, 4, 0, eid] + list(conn)
            f.write("".join(f"{v:10d}" for v in vals) + "\n")
        f.write("-1\n")
        f.write("EN ,R5.5,ATTR,     -1\n")

        for name in sorted(patch_nodes.keys()):
            if patch_nodes[name]:
                write_cmblock(f, name, patch_nodes[name])

    # Write volumes file for the build script
    vol_path = out_file.with_suffix(".volumes.txt")
    with vol_path.open("w") as vf:
        for eid, vol in enumerate(cell_vols, start=1):
            vf.write(f"{eid} {vol:.15e}\n")

    print(f"[OK] {out_file}")
    print(f"  Total fluid nodes: {len(fluid_nodes)} (OpenFOAM points: {len(points)})")
    print(f"  Fluid elements: {len(cell_conn)} (degenerate C3D8 from tetra)")
    print(f"  Node sets: {list(patch_nodes.keys())}")
    for sname in sorted(patch_nodes.keys()):
        n = len(patch_nodes[sname])
        if n > 0:
            print(f"    {sname}: {n} nodes")
    print(f"  Volumes saved to: {vol_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert OpenFOAM tetra mesh to CDB-like fluid input"
    )
    parser.add_argument("--case", required=True, help="OpenFOAM case directory")
    parser.add_argument("--out-cdb", required=True, help="Output CDB file path")
    parser.add_argument(
        "--fsi-patch",
        default="AGARD_WING",
        help="Patch name mapped to FSI_FLUID (default: AGARD_WING)",
    )
    args = parser.parse_args()

    case_dir = Path(args.case).resolve()
    out_file = Path(args.out_cdb).resolve()
    convert_case(case_dir, out_file, args.fsi_patch)


if __name__ == "__main__":
    main()
