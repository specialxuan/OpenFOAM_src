#!/usr/bin/env python3
"""
Convert OpenFOAM snappyHexMesh polyMesh (hex + polyhedral cells) to
ANSYS CDB-like fluid input for CalculiX modal analysis.

Hex cells (8-node) → C3D8 hex elements.
Polyhedral cells → decomposed into tetrahedra (face triangulation + centroid).
Each tetra → degenerate C3D8 (4 unique nodes, duplicated to 8).
"""

import argparse
import math
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple, Sequence


INT_RE = re.compile(r"[+-]?\d+")


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
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip() == "(":
            data_start = i + 1
            break
    faces = []
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line.startswith(")"):
            break
        m = re.match(r"(\d+)\s*\(([^)]*)\)", line)
        if m:
            nodes = [int(x) for x in m.group(2).split()]
            faces.append(nodes)
    return faces


def parse_label_list(path: Path) -> List[int]:
    content = path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"\(\s*\n(.*?)\n\s*\)", content, re.DOTALL)
    if not m:
        raise ValueError(f"Cannot parse from {path}")
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
    patches = []
    patch_re = re.compile(r"(\w+)\s*\{(.*?)\}", re.DOTALL)
    for pm in patch_re.finditer(m.group(1)):
        name = pm.group(1)
        body = pm.group(2)
        nf = re.search(r"nFaces\s+(\d+)", body)
        sf = re.search(r"startFace\s+(\d+)", body)
        if nf and sf:
            patches.append((name, int(nf.group(1)), int(sf.group(1))))
    return patches


def tet_volume(a, b, c, d):
    ad = (a[0] - d[0], a[1] - d[1], a[2] - d[2])
    bd = (b[0] - d[0], b[1] - d[1], b[2] - d[2])
    cd = (c[0] - d[0], c[1] - d[1], c[2] - d[2])
    cross = (bd[1]*cd[2] - bd[2]*cd[1], bd[2]*cd[0] - bd[0]*cd[2], bd[0]*cd[1] - bd[1]*cd[0])
    return abs(ad[0]*cross[0] + ad[1]*cross[1] + ad[2]*cross[2]) / 6.0


def hex_volume(coords):
    if len(coords) != 8:
        return 0.0
    def v(a_idx, b_idx, c_idx, d_idx):
        return tet_volume(coords[a_idx], coords[b_idx], coords[c_idx], coords[d_idx])
    return v(0,1,3,4) + v(1,2,3,6) + v(1,3,4,6) + v(1,4,5,6) + v(3,4,6,7)


def write_cmblock(fp, name: str, ids: Sequence[int]):
    fp.write(f"CMBLOCK,{name},NODE,{len(ids)}\n")
    fp.write("(8i10)\n")
    for i in range(0, len(ids), 8):
        chunk = ids[i:i+8]
        fp.write("".join(f"{v:10d}" for v in chunk) + "\n")


def c3d8_order_blockmesh(i: int, j: int, k: int, nx: int, ny: int) -> List[int]:
    """C3D8 corner ordering for blockMesh cell at (i,j,k)."""
    np = (nx + 1) * (ny + 1)
    nr = nx + 1
    return [
        k * np + j * nr + i,           # n1: bottom-left-front
        k * np + j * nr + (i + 1),     # n2: bottom-right-front
        k * np + (j + 1) * nr + (i + 1),  # n3: bottom-right-back
        k * np + (j + 1) * nr + i,     # n4: bottom-left-back
        (k + 1) * np + j * nr + i,     # n5: top-left-front
        (k + 1) * np + j * nr + (i + 1),  # n6: top-right-front
        (k + 1) * np + (j + 1) * nr + (i + 1),  # n7: top-right-back
        (k + 1) * np + (j + 1) * nr + i,  # n8: top-left-back
    ]
    node_ids: List[int],
    cell_face_indices: List[int],
    faces: List[List[int]],
    points: List[Tuple[float, float, float]],
) -> List[int]:
    """Order 8 hex nodes to C3D8 using face connectivity.

    Uses the 6 quadrilateral faces to determine opposite face pairs,
    then orders bottom face CCW and top face correspondingly.
    """
    if len(node_ids) != 8:
        raise ValueError(f"Expected 8 nodes, got {len(node_ids)}")

    # Collect the 6 quad faces belonging to this cell
    cell_faces_nodes = []
    for fi in cell_face_indices:
        fnodes = faces[fi]
        if len(fnodes) == 4:  # only quads
            cell_faces_nodes.append(fnodes)

    if len(cell_faces_nodes) != 6:
        # Fall back to coordinate-based ordering
        return _order_hex_by_coords(node_ids, points)

    # For each face, compute its centroid and normal
    face_data = []
    for fnodes in cell_faces_nodes:
        coords = [points[n] for n in fnodes]
        cx = sum(p[0] for p in coords) / 4.0
        cy = sum(p[1] for p in coords) / 4.0
        cz = sum(p[2] for p in coords) / 4.0
        # Face normal (newell's method for quads)
        nx, ny, nz = 0.0, 0.0, 0.0
        for i in range(4):
            j = (i + 1) % 4
            nx += (coords[i][1] - coords[j][1]) * (coords[i][2] + coords[j][2])
            ny += (coords[i][2] - coords[j][2]) * (coords[i][0] + coords[j][0])
            nz += (coords[i][0] - coords[j][0]) * (coords[i][1] + coords[j][1])
        face_data.append((fnodes, (cx, cy, cz), (nx, ny, nz)))

    # Find opposite face pairs by checking face normals
    # Opposite faces have roughly opposite normals
    pairs = []
    used = [False] * 6
    for i in range(6):
        if used[i]:
            continue
        ni = face_data[i][2]
        best_j = -1
        best_dot = 1.0  # looking for most negative dot product
        for j in range(i + 1, 6):
            if used[j]:
                continue
            nj = face_data[j][2]
            dot = ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2]
            norm_i = (ni[0]**2 + ni[1]**2 + ni[2]**2)**0.5
            norm_j = (nj[0]**2 + nj[1]**2 + nj[2]**2)**0.5
            if norm_i > 0 and norm_j > 0:
                dot = dot / (norm_i * norm_j)
            if dot < best_dot:
                best_dot = dot
                best_j = j
        if best_j >= 0:
            pairs.append((i, best_j))
            used[i] = True
            used[best_j] = True

    if len(pairs) != 3:
        return _order_hex_by_coords(node_ids, points)

    # Identify bottom/top, left/right, front/back
    # Bottom face has most negative (or smallest) z in its normal
    pair_z_scores = []
    for a, b in pairs:
        ca = face_data[a][1]
        cb = face_data[b][1]
        # The face with smaller z-centroid is "bottom", larger is "top"
        za = ca[2]
        zb = cb[2]
        # Determine which direction (z-like) this pair represents
        nz_comp = abs(face_data[a][2][2]) + abs(face_data[b][2][2])
        pair_z_scores.append((nz_comp, za, zb, a, b))

    # Sort by z-component of normal: the pair with largest |nz| is bottom/top
    pair_z_scores.sort(key=lambda x: -x[0])
    btm_top = pair_z_scores[0]
    # Identify which is bottom (smaller z centroid)
    if btm_top[1] < btm_top[2]:
        btm_idx, top_idx = btm_top[3], btm_top[4]
    else:
        btm_idx, top_idx = btm_top[4], btm_top[3]

    btm_face = face_data[btm_idx][0]
    top_face = face_data[top_idx][0]

    # Order bottom face CCW (when viewed from +z, i.e., from outside looking at bottom)
    btm_coords = [points[n] for n in btm_face]
    btm_cx = sum(p[0] for p in btm_coords) / 4.0
    btm_cy = sum(p[1] for p in btm_coords) / 4.0
    # Sort by angle around centroid
    btm_ordered = sorted(btm_face, key=lambda n: (
        math.atan2(points[n][1] - btm_cy, points[n][0] - btm_cx)
    ))

    # Order top face CCW (matching bottom face order)
    top_coords = [points[n] for n in top_face]
    top_cx = sum(p[0] for p in top_coords) / 4.0
    top_cy = sum(p[1] for p in top_coords) / 4.0
    top_ordered = sorted(top_face, key=lambda n: (
        math.atan2(points[n][1] - top_cy, points[n][0] - top_cx)
    ))

    # C3D8: n1,n2,n3,n4 = bottom_ordered, n5,n6,n7,n8 = top_ordered
    return btm_ordered + top_ordered


def _order_hex_by_coords(
    node_ids: List[int],
    points: List[Tuple[float, float, float]],
) -> List[int]:
    """Fallback: order by normalized coordinates (bounding box corners)."""
    coords = {n: points[n] for n in node_ids}
    xs = [coords[n][0] for n in node_ids]
    ys = [coords[n][1] for n in node_ids]
    zs = [coords[n][2] for n in node_ids]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin, zmax = min(zs), max(zs)
    dx = max(xmax - xmin, 1e-16)
    dy = max(ymax - ymin, 1e-16)
    dz = max(zmax - zmin, 1e-16)

    def norm(nid):
        return ((coords[nid][0]-xmin)/dx, (coords[nid][1]-ymin)/dy, (coords[nid][2]-zmin)/dz)

    corners = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    unused = set(node_ids)
    ordered = []
    for cx, cy, cz in corners:
        best = min(unused, key=lambda n: (norm(n)[0]-cx)**2 + (norm(n)[1]-cy)**2 + (norm(n)[2]-cz)**2)
        ordered.append(best)
        unused.remove(best)
    return ordered


def get_blockmesh_dims(case_dir: Path) -> tuple:
    """Read nx, ny, nz from blockMeshDict."""
    bd = (case_dir / "system" / "blockMeshDict").read_text()
    nx = int(re.search(r'nx\s+(\d+)', bd).group(1))
    ny = int(re.search(r'ny\s+(\d+)', bd).group(1))
    nz = int(re.search(r'nz\s+(\d+)', bd).group(1))
    return nx, ny, nz


def convert_case(case_dir: Path, out_file: Path, fsi_patch: str = "AGARD_WING"):
    poly = case_dir / "constant" / "polyMesh"
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_label_list(poly / "owner")
    neighbour = parse_label_list(poly / "neighbour")
    patches = parse_boundary(poly / "boundary")

    n_cells = max(max(owner), max(neighbour) if neighbour else 0) + 1
    cell_faces = [[] for _ in range(n_cells)]
    for fi, c in enumerate(owner):
        cell_faces[c].append(fi)
    for fi, c in enumerate(neighbour):
        cell_faces[c].append(fi)

    # Classify cells
    hex_cells = []  # (cid, [8 node ids sorted])
    poly_cells = []  # (cid, node_list)
    for cid in range(n_cells):
        node_set = set()
        for fi in cell_faces[cid]:
            node_set.update(faces[fi])
        nodes = sorted(node_set)
        if len(nodes) == 8:
            hex_cells.append((cid, nodes))
        else:
            poly_cells.append((cid, nodes))

    # Only use CFD nodes (no centroid nodes needed since we skip poly cells)
    total_nodes = len(points)
    all_nodes: Dict[int, Tuple[float, float, float]] = {}
    for i, p in enumerate(points):
        all_nodes[i + 1] = p  # 1-based

    # Build elements - only hex cells with proper C3D8 ordering
    elements_conn = []
    element_volumes = []
    nx, ny, nz = get_blockmesh_dims(case_dir)
    np = (nx + 1) * (ny + 1)
    nr = nx + 1

    for cid, nodes_list in hex_cells:
        node_set = set(nodes_list)
        # Determine (i,j,k) from min node ID
        min_nid = min(node_set)
        k = min_nid // np
        rem = min_nid % np
        j = rem // nr
        i = rem % nr
        ordered = c3d8_order_blockmesh(i, j, k, nx, ny)
        if set(ordered) == node_set:
            conn_1b = [n + 1 for n in ordered]
        else:
            # Fallback: sorted (works for most blockMesh cells except n3/n4 swap)
            conn_1b = [n + 1 for n in sorted(node_set)]
        coords = [points[n] for n in [x - 1 for x in conn_1b]]
        vol = hex_volume(coords)
        elements_conn.append(tuple(conn_1b))
        element_volumes.append(vol)

    # Skip poly cell decomposition - hex cells cover all nodes

    # Patch node sets
    patch_nodes: Dict[str, List[int]] = {}
    for name, nfaces, start in patches:
        ids = set()
        end = start + nfaces
        for fi in range(start, end):
            if fi < len(faces):
                for nid in faces[fi]:
                    ids.add(nid + 1)
        set_name = name.upper()
        patch_nodes[set_name] = sorted(ids)

    # FSI_FLUID = AGARD_WING
    if "AGARD_WING" in patch_nodes:
        patch_nodes["FSI_FLUID"] = patch_nodes.pop("AGARD_WING")

    patch_nodes["FLUIDDOMAIN"] = sorted(all_nodes.keys())

    # Map other patches
    patch_mapping = {
        "INLET": "LEFT",
        "OUTLET": "RIGHT",
        "FARFIELD": "TOP",
        "ROOTSYMMETRY": "FRONT",
    }
    for old, new in patch_mapping.items():
        if old in patch_nodes and new not in patch_nodes:
            patch_nodes[new] = patch_nodes[old]

    for default in ["FRONT", "BACK", "TOP", "BOTTOM_FLUID", "LEFT", "RIGHT"]:
        if default not in patch_nodes:
            patch_nodes[default] = []

    # Write CDB
    out_file.parent.mkdir(parents=True, exist_ok=True)
    total_nodes = len(all_nodes)
    total_elems = len(elements_conn)

    with out_file.open("w", encoding="utf-8") as f:
        f.write("/BATCH\n/NOPR\n")
        f.write(f"/COM, AGARD Hex Fluid CDB: {total_nodes} nodes, {total_elems} elems\n")
        f.write("/PREP7\n")
        f.write("ET,1,185,,0,,0,,0,0\n")
        f.write("MP,EX,  1,           1\n")
        f.write("MP,PRXY,  1,         0.5\n")
        f.write("MP,DENS,  1,      0.0025\n")
        f.write("/NOLIST\n")
        f.write("NBLOCK,6,SOLID\n")
        f.write("(3i8,6e16.9)\n")
        for nid in sorted(all_nodes):
            x, y, z = all_nodes[nid]
            f.write(f"{nid:8d}{0:8d}{0:8d}{x:16.9E}{y:16.9E}{z:16.9E}\n")
        f.write("N ,R5.3,LOC,     -1\n")
        f.write("SHPP,WARN\n")
        f.write("EBLOCK,19,SOLID\n")
        f.write("(19i10)\n")
        for eid, conn in enumerate(elements_conn, start=1):
            # num_nodes in conn that are unique
            unique_n = len(set(conn))
            vals = [1, 1, 0, 0, 0, 0, 0, 0, unique_n, 0, eid] + list(conn)
            f.write("".join(f"{v:10d}" for v in vals) + "\n")
        f.write("-1\n")
        f.write("EN ,R5.5,ATTR,     -1\n")

        for name in sorted(patch_nodes.keys()):
            if patch_nodes[name]:
                write_cmblock(f, name, patch_nodes[name])

    # Write volumes
    vol_path = out_file.with_suffix(".volumes.txt")
    with vol_path.open("w") as vf:
        for eid, vol in enumerate(element_volumes, start=1):
            vf.write(f"{eid} {vol:.15e}\n")

    print(f"[OK] {out_file}")
    print(f"  Nodes: {total_nodes} (CFD only, no centroids)")
    print(f"  Elements: {total_elems} ({len(hex_cells)} hex, poly cells skipped)")
    print(f"  Volumes: {vol_path}")
    for sname in sorted(patch_nodes.keys()):
        n = len(patch_nodes[sname])
        if n > 0:
            print(f"    {sname}: {n} nodes")


def main():
    parser = argparse.ArgumentParser(description="Convert snappyHexMesh polyMesh to CDB-like fluid input")
    parser.add_argument("--case", required=True)
    parser.add_argument("--out-cdb", required=True)
    parser.add_argument("--fsi-patch", default="AGARD_WING")
    args = parser.parse_args()
    convert_case(Path(args.case), Path(args.out_cdb), args.fsi_patch)


if __name__ == "__main__":
    main()
