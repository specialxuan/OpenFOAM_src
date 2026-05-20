#!/usr/bin/env python3
"""
Split snappyHexMesh at refinement boundaries: remove polyhedral cells,
expose interfaces as TIE_MASTER / TIE_SLAVE boundary patches.

After this step, the mesh contains ONLY hex cells, with two new boundary
patches at the (former) refinement interface. These patches will be
bonded via CalculiX *TIE constraints.
"""

import argparse, re, shutil
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict


def parse_points(path: Path) -> List[Tuple[float,float,float]]:
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


def parse_boundary(path: Path) -> List[Tuple[str, int, int, str]]:
    """Return list of (name, nFaces, startFace, type)."""
    text = path.read_text()
    m = re.search(r"\(\s*\n(.*?)\n\s*\)", text, re.DOTALL)
    if not m:
        return []
    patches = []
    patch_re = re.compile(r"(\w+)\s*\{(.+?)\}", re.DOTALL)
    for pm in patch_re.finditer(m.group(1)):
        name = pm.group(1)
        body = pm.group(2)
        nf = re.search(r"nFaces\s+(\d+)", body)
        sf = re.search(r"startFace\s+(\d+)", body)
        typ = "patch"
        tm = re.search(r"type\s+(\w+)", body)
        if tm: typ = tm.group(1)
        if nf and sf:
            patches.append((name, int(nf.group(1)), int(sf.group(1)), typ))
    return patches


def write_foam_header(fp):
    fp.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
    fp.write("| =========                 |                                                 |\n")
    fp.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
    fp.write("|  \\\\    /   O peration     | Version:  v2412                                 |\n")
    fp.write("|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n")
    fp.write("|    \\\\/     M anipulation  |                                                 |\n")
    fp.write("\\*---------------------------------------------------------------------------*/\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--backup", action="store_true", help="Backup original polyMesh before modifying")
    args = parser.parse_args()

    case = Path(args.case)
    poly = case / "constant" / "polyMesh"

    # ---- Parse original mesh ----
    points = parse_points(poly / "points")
    faces = parse_faces(poly / "faces")
    owner = parse_ilist(poly / "owner")
    neighbour = parse_ilist(poly / "neighbour")
    orig_patches = parse_boundary(poly / "boundary")

    n_cells = max(max(owner), max(neighbour)) + 1
    n_internal = len(neighbour)
    n_total_faces = len(owner)

    # Cell → faces map
    cf = [[] for _ in range(n_cells)]
    for fi, c in enumerate(owner): cf[c].append(fi)
    for fi, c in enumerate(neighbour): cf[c].append(fi)

    # Classify cells
    hex_cells: Set[int] = set()
    poly_cells: Set[int] = set()
    for cid in range(n_cells):
        ns = set()
        for fi in cf[cid]: ns.update(faces[fi])
        if len(ns) == 8: hex_cells.add(cid)
        else: poly_cells.add(cid)

    print(f"Total cells: {n_cells}")
    print(f"Hex cells: {len(hex_cells)}")
    print(f"Poly cells: {len(poly_cells)}")

    if len(poly_cells) == 0:
        print("No poly cells found. Nothing to split.")
        return

    # ---- Find interface faces (hex<->poly internal faces) ----
    interface_hex_to_poly: List[Tuple[int, int, int]] = []  # (face_idx, hex_cid, poly_cid)
    for fi in range(n_internal):
        c = owner[fi]
        n = neighbour[fi]
        if c in hex_cells and n in poly_cells:
            interface_hex_to_poly.append((fi, c, n))
        elif n in hex_cells and c in poly_cells:
            interface_hex_to_poly.append((fi, n, c))

    print(f"Interface faces (hex<->poly): {len(interface_hex_to_poly)}")

    # ---- Determine fine vs coarse hex cells ----
    # Use cell volume as proxy: smaller cells = finer (near wing)
    cell_sizes = {}
    for cid in hex_cells:
        ns = set()
        for fi in cf[cid]: ns.update(faces[fi])
        xs = [points[n][0] for n in ns]
        ys = [points[n][1] for n in ns]
        zs = [points[n][2] for n in ns]
        vol = (max(xs)-min(xs)) * (max(ys)-min(ys)) * (max(zs)-min(zs))
        cell_sizes[cid] = vol

    sizes = list(cell_sizes.values())
    sizes.sort()
    median_size = sizes[len(sizes)//2]
    print(f"Median hex cell volume: {median_size:.6e}")

    # Cells smaller than median are "fine", larger are "coarse"
    # (Fine cells are near the wing from refinement)
    # For TIE: master = coarser surface, slave = finer surface (standard convention)

    # ---- Assign interface faces to TIE patches ----
    tie_coarse_faces = []  # faces on the coarse hex side
    tie_fine_faces = []    # faces on the fine hex side

    for fi, hex_cid, poly_cid in interface_hex_to_poly:
        if cell_sizes.get(hex_cid, 0) > median_size:
            tie_coarse_faces.append(fi)
        else:
            tie_fine_faces.append(fi)

    print(f"TIE_COARSE faces: {len(tie_coarse_faces)}")
    print(f"TIE_FINE faces: {len(tie_fine_faces)}")

    # ---- Rebuild mesh: remove poly cells ----
    # face_map: old_face_idx -> new_face_idx
    # Poly cell faces are REMOVED (they would be orphaned)
    # Hex-hex internal faces are kept (reordered)
    # Original boundary faces are kept
    # TIE faces become new boundary patches

    # Identify all faces that need to be removed:
    # 1. Faces between poly-poly (internal to poly region)
    # 2. Faces between poly-hex (become TIE boundary faces)
    # We KEEP hex-hex faces and original boundary faces

    # Actually, we need to handle this carefully:
    # - Poly cell faces between POLY-POLY: REMOVE entirely
    # - Poly cell faces between POLY-HEX: become new TIE boundary
    # - Hex cell faces between HEX-HEX: KEEP as internal → may need reordering
    # - Original boundary faces: KEEP

    # Build set of faces to remove (poly-poly internal + extra poly faces)
    poly_internal_faces = set()
    for fi in range(n_internal):
        c = owner[fi]
        n = neighbour[fi]
        if c in poly_cells and n in poly_cells:
            poly_internal_faces.add(fi)

    # TIE faces (hex-poly) become new boundary faces
    tie_face_set = set(fi for fi, _, _ in interface_hex_to_poly)

    # Faces to remove: poly internal + poly cell boundary faces that are NOT TIE faces
    # (poly cells may have faces on existing boundary patches too)
    faces_to_remove = set(poly_internal_faces)
    # Also remove any boundary faces that belong to poly cells
    # (these are poly cell faces on external boundaries - orphaned after poly cell removal)
    # After removing poly cells, any face owned by a poly cell is orphaned

    # faces owned by poly cells
    poly_owned_faces = set()
    for fi in range(n_total_faces):
        if owner[fi] in poly_cells:
            poly_owned_faces.add(fi)

    # TIE faces: hex-owned but neighbor is poly → keep as TIE boundary
    # poly-owned faces: remove (they'll be orphaned)

    # Build new face order:
    # 1. Hex-hex internal faces (owner and neighbour both hex)
    # 2. TIE boundary faces (one side is hex, other side was poly)
    # 3. Original boundary faces that are NOT poly-owned
    new_faces_order = []
    face_remap = {}  # old_fi -> new_fi

    # Category 1: hex-hex internal faces
    new_internal = []
    for fi in range(n_internal):
        c = owner[fi]
        n = neighbour[fi]
        if c in hex_cells and n in hex_cells:
            new_internal.append(fi)
    new_n_internal = len(new_internal)

    # Category 2: TIE boundary faces (hex-poly → new boundary)
    # Include ALL hex-poly faces, regardless of which cell owns them.
    # If poly owns the face, reassign owner to the hex cell.
    tie_coarse_new = []
    tie_fine_new = []
    reassigned_owners = {}  # old_fi -> new_owner_cid

    for fi in tie_coarse_faces:
        hex_cid = owner[fi] if owner[fi] in hex_cells else neighbour[fi]
        tie_coarse_new.append(fi)
        if owner[fi] not in hex_cells:
            reassigned_owners[fi] = hex_cid

    for fi in tie_fine_faces:
        hex_cid = owner[fi] if owner[fi] in hex_cells else neighbour[fi]
        tie_fine_new.append(fi)
        if owner[fi] not in hex_cells:
            reassigned_owners[fi] = hex_cid

    # Category 3: Original boundary faces not owned by poly cells
    orig_boundary_new = {}
    boundary_base = n_internal
    for name, nf, sf, typ in orig_patches:
        kept = []
        for fi in range(sf, sf + nf):
            if fi not in poly_owned_faces and fi not in poly_internal_faces:
                kept.append(fi)
        if kept:
            orig_boundary_new[name] = (kept, typ)

    # Build new face list
    new_face_list = new_internal + tie_coarse_new + tie_fine_new
    face_start = len(new_face_list)
    for kept, _ in orig_boundary_new.values():
        new_face_list.extend(kept)

    # Build face remap
    for new_fi, old_fi in enumerate(new_face_list):
        face_remap[old_fi] = new_fi

    # ---- Rebuild owner and neighbour ----
    new_owner = []
    new_neighbour = []
    neighbor_for_internal = set(new_internal)  # only internal faces have neighbours

    for old_fi in new_face_list:
        if old_fi in reassigned_owners:
            new_owner.append(reassigned_owners[old_fi])
        else:
            new_owner.append(owner[old_fi])
        if old_fi in neighbor_for_internal:
            new_neighbour.append(neighbour[old_fi])
        else:
            # Boundary face → no neighbour
            # But we still need neighbour list to match count for faces in TIE
            # Actually, in O.F., neighbour only exists for internal faces
            pass

    # In OpenFOAM, neighbour list only has entries for internal faces
    new_neighbour = [neighbour[fi] for fi in new_internal]

    # ---- Write new boundary file ----
    tie_coarse_start = new_n_internal
    tie_fine_start = tie_coarse_start + len(tie_coarse_new)
    boundary_start = tie_fine_start + len(tie_fine_new)

    new_patches = [
        ("TIE_COARSE", len(tie_coarse_new), tie_coarse_start, "patch"),
        ("TIE_FINE", len(tie_fine_new), tie_fine_start, "patch"),
    ]
    for name, (kept, typ) in orig_boundary_new.items():
        sf = boundary_start
        new_patches.append((name, len(kept), sf, typ))
        boundary_start += len(kept)

    # ---- Backup original if requested ----
    if args.backup:
        backup_dir = poly.parent / "polyMesh.before_split"
        if not backup_dir.exists():
            shutil.copytree(poly, backup_dir)
            print(f"Backup: {backup_dir}")

    # ---- Write new mesh ----
    # points: unchanged
    # faces: reordered
    # owner: rebuilt
    # neighbour: rebuilt
    # boundary: rebuilt

    # Write faces
    with (poly / "faces").open("w") as f:
        write_foam_header(f)
        f.write("FoamFile { version 2.0; format ascii; class faceList; location \"constant/polyMesh\"; object faces; }\n")
        f.write(f"{len(new_face_list)}\n(\n")
        for old_fi in new_face_list:
            fnodes = faces[old_fi]
            f.write(f"{len(fnodes)}(" + " ".join(str(n) for n in fnodes) + ")\n")
        f.write(")\n")

    # Write owner
    with (poly / "owner").open("w") as f:
        write_foam_header(f)
        f.write("FoamFile { version 2.0; format ascii; class labelList; location \"constant/polyMesh\"; object owner; }\n")
        f.write(f"{len(new_owner)}\n(\n")
        for v in new_owner: f.write(f"{v}\n")
        f.write(")\n")

    # Write neighbour
    with (poly / "neighbour").open("w") as f:
        write_foam_header(f)
        f.write("FoamFile { version 2.0; format ascii; class labelList; location \"constant/polyMesh\"; object neighbour; }\n")
        f.write(f"{len(new_neighbour)}\n(\n")
        for v in new_neighbour: f.write(f"{v}\n")
        f.write(")\n")

    # Write boundary
    with (poly / "boundary").open("w") as f:
        write_foam_header(f)
        f.write("FoamFile { version 2.0; format ascii; class polyBoundaryMesh; location \"constant/polyMesh\"; object boundary; }\n")
        f.write(f"{len(new_patches)}\n(\n")
        for name, nf, sf, typ in new_patches:
            groups = ""
            if typ == "wall": groups = "        inGroups        1(wall);\n"
            elif typ == "symmetryPlane": groups = "        inGroups        1(symmetryPlane);\n"
            f.write(f"    {name}\n    {{\n        type            {typ};\n{groups}        nFaces          {nf};\n        startFace       {sf};\n    }}\n")
        f.write(")\n")

    # ---- Summary ----
    print(f"\n=== New mesh summary ===")
    print(f"  Internal faces: {new_n_internal}")
    print(f"  TIE_COARSE: {len(tie_coarse_new)} faces")
    print(f"  TIE_FINE: {len(tie_fine_new)} faces")
    print(f"  Original boundary patches:")
    for name, nf, sf, typ in new_patches[2:]:
        print(f"    {name}: {nf} faces")
    print(f"  Total faces: {len(new_face_list)}")
    print(f"  All cells: hex only ({len(hex_cells)} cells, poly removed)")
    print(f"\nRun checkMesh to verify.")


if __name__ == "__main__":
    main()
