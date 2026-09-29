#!/usr/bin/env python3
"""Build the body-fitted AGARD 445.6 fluid mesh (OpenFOAM polyMesh + CalculiX deck).

Replaces the ICEM-derived far-field mesh whose wing-tip cells degenerated to
~6 um.  Built with the gmsh 4.15 python API (OCC kernel for the CAD, the
built-in kernel for the 3D boundary-layer extrusion):

* domain: the EXISTING far field, reproduced exactly -- ``y in [0, 7]``;
  ``x <= 0`` a half-cylinder of radius 5 about the y axis; ``x >= 0`` a box
  ``x in [0, 6]``, ``|z| <= 5``.  Patches ``FRONT`` (y=0 symmetry plane),
  ``BACK`` (y=7), ``IN`` (cylinder), ``OUT`` (x=6), ``TOP``/``BOTTOM``
  (z=+/-5) and ``FSI`` (the wing wetted surface).
* wing: ``Wing.step`` (MILLIMETRES, half wing, root at span 0) read by OCC and
  permuted into the mesh frame: ``mesh_x = step_z/1000`` (chord),
  ``mesh_y = step_x/1000`` (span), ``mesh_z = step_y/1000`` (thickness).
* core: tetrahedra; boundary layer: triangular prisms (C3D6) extruded from the
  wing, y+~1 first layer, graded growth.

Symmetry clip
-------------
The wing root sits on the symmetry plane.  The swept wing's surface normal at
the root free edge has a y component, so the extruded BL collar leans across
y=0 by up to the BL thickness.  The collar nodes are therefore projected onto
y=0 after meshing and the mesh is relaxed with Netgen; this clips the BL at
the symmetry plane (the physically correct treatment) and keeps the domain at
exactly y in [0, 7] with a planar FRONT patch.

All size-field parameters are CLI options (see --help); the defaults are the
production values.  Deterministic: gmsh is deterministic for a fixed geometry
and field setup.

Outputs (under --out-dir):
  fluid_mesh.msh           gmsh MSH 2.2 (gmshToFoam input, physical groups)
  fluid_mesh.geo_unrolled  gmsh geometry snapshot
  fluid_mesh.inp           CalculiX solid deck (C3D4 + C3D6, NSETs)
  build.log                stdout of this run
  mesh_report.json         measured counts/metrics
  README.md                parameter + boundary documentation
  openfoam/                polyMesh + minimal system/ (if --openfoam)
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

STEP_DEFAULT = Path("/root/Workspace/AGARD/AGARD4456/Geometry/Wing.step")
R_CYL = 5.0
X_BOX = 6.0
Y_SPAN = 7.0
PATCH_TYPES = {
    "FRONT": "symmetryPlane",
    "FSI": "wall",
    "BACK": "patch",
    "IN": "patch",
    "OUT": "patch",
    "TOP": "patch",
    "BOTTOM": "patch",
}


class MeshError(RuntimeError):
    """A geometry or meshing contract was violated."""


def _log(msg: str) -> None:
    print(f"[agard-mesh] {msg}", flush=True)


def _bbox(gmsh, dim: int, tag: int) -> list[float]:
    """Bounding box of an OCC or discrete entity (no booleans are used)."""
    try:
        b = [float(v) for v in gmsh.model.getBoundingBox(dim, tag)]
        if all(np.isfinite(b)):
            return b
    except Exception:  # noqa: BLE001 - discrete BL entities have no geometry
        pass
    verts = [e[1] for e in gmsh.model.getBoundary(
        [(dim, tag)], recursive=True, oriented=False) if e[0] == 0]
    if not verts:
        raise MeshError(f"no bounding box for entity ({dim}, {tag})")
    c = np.array([gmsh.model.getValue(0, v, []) for v in verts])
    return [float(c[:, 0].min()), float(c[:, 1].min()), float(c[:, 2].min()),
            float(c[:, 0].max()), float(c[:, 1].max()), float(c[:, 2].max())]


def _vertex_bbox(gmsh, dim: int, tag: int) -> list[float]:
    return _bbox(gmsh, dim, tag)


def _at_y(gmsh, dim: int, tag: int, y: float, tol: float = 1e-6) -> bool:
    b = _vertex_bbox(gmsh, dim, tag)
    return abs(b[1] - y) < tol and abs(b[4] - y) < tol


def build_geometry(gmsh, step: Path, heal_tolerance: float = 0.0) -> dict:
    """Outer domain + transformed wing; return the tags needed downstream."""
    occ = gmsh.model.occ
    o = occ.addPoint(0.0, 0.0, 0.0)
    a = occ.addPoint(0.0, 0.0, R_CYL)
    b = occ.addPoint(-R_CYL, 0.0, 0.0)
    c = occ.addPoint(0.0, 0.0, -R_CYL)
    d = occ.addPoint(X_BOX, 0.0, -R_CYL)
    e = occ.addPoint(X_BOX, 0.0, R_CYL)
    arc1 = occ.addCircleArc(a, o, b)
    arc2 = occ.addCircleArc(b, o, c)
    line1 = occ.addLine(c, d)
    line2 = occ.addLine(d, e)
    line3 = occ.addLine(e, a)
    loop = occ.addCurveLoop([arc1, arc2, line1, line2, line3])
    face = occ.addPlaneSurface([loop])
    ext = occ.extrude([(2, face)], 0.0, Y_SPAN, 0.0)
    occ.synchronize()
    v_outer = next(tag for dim, tag in ext if dim == 3)
    outer_surfs = [e[1] for e in gmsh.model.getBoundary(
        [(3, v_outer)], oriented=False)]
    front_tag = next(t for t in outer_surfs if _at_y(gmsh, 2, t, 0.0))
    front_wire = gmsh.model.getBoundary(
        [(2, front_tag)], combined=True, oriented=True)

    imported = occ.importShapes(str(step))
    occ.synchronize()
    v_wing = next(e[1] for e in imported if e[0] == 3)
    transform = np.array(
        [[0, 0, 1.0 / 1000, 0],
         [1.0 / 1000, 0, 0, 0],
         [0, 1.0 / 1000, 0, 0],
         [0, 0, 0, 1]], dtype=float)
    occ.affineTransform([(3, v_wing)], transform.flatten().tolist())
    occ.synchronize()
    if heal_tolerance > 0.0:
        occ.healShapes([(3, v_wing)], heal_tolerance, fixDegenerated=True,
                       fixSmallEdges=True, fixSmallFaces=True, sewFaces=True,
                       makeSolids=True)
        occ.synchronize()
        candidates = [e[1] for e in gmsh.model.getEntities(3)
                      if e[1] != v_outer]
        v_wing = max(candidates, key=lambda t: occ.getMass(3, t))
        outer_surfs = [e[1] for e in gmsh.model.getBoundary(
            [(3, v_outer)], oriented=False)]
        front_tag = next(t for t in outer_surfs if _at_y(gmsh, 2, t, 0.0))
        front_wire = gmsh.model.getBoundary(
            [(2, front_tag)], combined=True, oriented=True)
    wing_bbox = list(gmsh.model.getBoundingBox(3, v_wing))
    if not (abs(wing_bbox[0]) < 1e-4 and abs(wing_bbox[3] - 1.177825) < 1e-3
            and abs(wing_bbox[1]) < 1e-6 and abs(wing_bbox[4] - 0.762) < 1e-3
            and abs(wing_bbox[2] + 0.01125) < 1e-3
            and abs(wing_bbox[5] - 0.01125) < 1e-3):
        raise MeshError(f"wing bbox after transform is wrong: {wing_bbox}")
    wing_surfs = [e[1] for e in gmsh.model.getBoundary(
        [(3, v_wing)], oriented=False)]
    root_face = next(t for t in wing_surfs if _at_y(gmsh, 2, t, 0.0))
    fsi = [t for t in wing_surfs if t != root_face]
    te_surfs = [t for t in fsi
                if (_bbox(gmsh, 2, t)[5] - _bbox(gmsh, 2, t)[2]) < 1e-3]
    root_wire = gmsh.model.getBoundary(
        [(2, root_face)], combined=True, oriented=True)
    return {"v_outer": v_outer, "v_wing": v_wing, "outer_surfs": outer_surfs,
            "front_tag": front_tag, "front_wire": front_wire,
            "root_face": root_face, "root_wire": root_wire, "fsi": fsi,
            "te_surfs": te_surfs, "wing_bbox": wing_bbox}


def extrude_boundary_layer(gmsh, fsi: list[int], layers: int, first: float,
                           ratio: float) -> dict:
    """Prism BL from the wing wetted surfaces; return top and free-edge collar."""
    heights = [first]
    for i in range(1, layers):
        heights.append(heights[-1] + first * ratio ** i)
    extbl = gmsh.model.geo.extrudeBoundaryLayer(
        [(2, t) for t in fsi], [1] * layers, heights, True)
    gmsh.model.geo.synchronize()
    top = [extbl[i - 1][1] for i in range(1, len(extbl)) if extbl[i][0] == 3]
    top_pairs = [extbl[i - 1] for i in range(1, len(extbl)) if extbl[i][0] == 3]
    top_set = set(top)
    side = [e[1] for e in extbl if e[0] == 2 and e[1] not in top_set]
    collar = [t for t in sorted(set(side))
              if len(gmsh.model.getAdjacencies(2, t)[0]) == 1]
    if len(top) != len(fsi):
        raise MeshError(f"BL top count {len(top)} != base count {len(fsi)}")
    if not collar:
        raise MeshError("no open BL collar found")
    return {"top": top, "top_pairs": top_pairs, "collar": collar,
            "thickness": heights[-1]}


def build_fluid_volume(gmsh, geom: dict, bl: dict,
                       te_surfaces: list[int]) -> dict:
    """Remove scaffolding and add the closed fluid volume; split root/TE collars."""
    occ = gmsh.model.occ
    for v in (geom["v_outer"], geom["v_wing"]):
        occ.remove([(3, v)], recursive=False)
    occ.remove([(2, geom["front_tag"])], recursive=False)
    occ.synchronize()
    gmsh.model.geo.synchronize()
    # The FRONT hole is the BL top's root free edge.  With the blunt TE
    # excluded from the BL the free boundary is a single root+TE loop, so the
    # root portion is closed with a straight seam line before making the hole.
    free = gmsh.model.getBoundary(bl["top_pairs"], combined=True,
                                  oriented=True)

    def endpoints(edge):
        return [x[1] for x in gmsh.model.getBoundary(
            [(1, edge)], recursive=True, oriented=False) if x[0] == 0]

    root_edges: list[tuple[int, list[int]]] = []
    te_edges: list[int] = []
    for _, tag in free:
        pts = endpoints(tag)
        coords = [gmsh.model.getValue(0, p, []) for p in pts]
        if (len(pts) == 2 and all(abs(c[1]) < 1e-6 for c in coords)):
            root_edges.append((tag, pts))
        else:
            te_edges.append(tag)

    def boundary_edges(surface):
        return {e[1] for e in gmsh.model.getBoundary(
            [(2, surface)], oriented=False)}

    te_excluded = bool(te_edges)
    if te_excluded:
        root_edge_tags = {tag for tag, _ in root_edges}
        root_collar = [t for t in bl["collar"]
                       if boundary_edges(t) & root_edge_tags]
        te_collar = [t for t in bl["collar"] if t not in root_collar]
        from collections import Counter
        degree: Counter = Counter()
        for _, pts in root_edges:
            degree[pts[0]] += 1
            degree[pts[1]] += 1
        seam = [p for p, n in degree.items() if n == 1]
        if len(seam) != 2:
            raise MeshError("cannot identify the root seam endpoints")
        line = gmsh.model.geo.addLine(seam[0], seam[1])
        try:
            cl_inner = gmsh.model.geo.addCurveLoop(
                [root_edges[0][0], line, root_edges[1][0]])
        except Exception:  # noqa: BLE001 - try the reversed orientation
            cl_inner = gmsh.model.geo.addCurveLoop(
                [root_edges[0][0], line, -root_edges[1][0]])
        extra = te_collar + list(te_surfaces)
    else:
        root_collar = list(bl["collar"])
        te_collar = []
        cl_inner = gmsh.model.geo.addCurveLoop([t[1] for t in free])
        extra = []
    cl_outer = gmsh.model.geo.addCurveLoop([t[1] for t in geom["front_wire"]])
    new_front = gmsh.model.geo.addPlaneSurface([cl_outer, cl_inner])
    gmsh.model.geo.synchronize()
    others = [t for t in geom["outer_surfs"] if t != geom["front_tag"]]
    loop = gmsh.model.geo.addSurfaceLoop(
        others + [new_front] + bl["top"] + extra)
    gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()
    return {"others": others, "new_front": new_front,
            "root_collar": root_collar, "te_collar": te_collar,
            "all_volumes": [e[1] for e in gmsh.model.getEntities(3)]}


def classify_outer(gmsh, others: list[int]) -> dict:
    """Group the six outer surfaces into named patches by their vertices."""
    groups: dict[str, list[int]] = {}
    for t in others:
        b = _vertex_bbox(gmsh, 2, t)
        if abs(b[1] - Y_SPAN) < 1e-6 and abs(b[4] - Y_SPAN) < 1e-6:
            name = "BACK"
        elif b[2] > R_CYL - 1e-3 and b[5] > R_CYL - 1e-3:
            name = "TOP"
        elif b[2] < -R_CYL + 1e-3 and b[5] < -R_CYL + 1e-3:
            name = "BOTTOM"
        elif b[0] > X_BOX - 1e-3 and b[3] > X_BOX - 1e-3:
            name = "OUT"
        elif b[3] < 1e-6:
            name = "IN"
        else:
            raise MeshError(f"unclassified outer surface {t}: {b}")
        groups.setdefault(name, []).append(t)
    return groups



def set_size_field(gmsh, wing_surfaces: list[int], args) -> None:
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", wing_surfaces)
    th = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(th, "InField", dist)
    gmsh.model.mesh.field.setNumber(th, "SizeMin", args.wing_size)
    gmsh.model.mesh.field.setNumber(th, "SizeMax", args.far_size)
    gmsh.model.mesh.field.setNumber(th, "DistMin", args.dist_min)
    gmsh.model.mesh.field.setNumber(th, "DistMax", args.dist_max)
    gmsh.model.mesh.field.setAsBackgroundMesh(th)
    gmsh.option.setNumber("Mesh.MeshSizeMin", args.mesh_size_min)
    gmsh.option.setNumber("Mesh.MeshSizeMax", args.far_size)
    gmsh.option.setNumber("Mesh.Smoothing", args.smoothing)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


def collar_inner_nodes(gmsh, surfaces: list[int]) -> tuple[set[int], dict]:
    """Collar nodes lying on y=0 before the symmetry projection, plus their positions."""
    tags: set[int] = set()
    for t in surfaces:
        tags.update(int(x) for x in gmsh.model.mesh.getNodes(2, t, True)[0])
    raw = gmsh.model.mesh.getNodes()
    coords = np.asarray(raw[1]).reshape(-1, 3)
    pos = {int(t): coords[i].copy() for i, t in enumerate(raw[0])}
    inner = {t for t in tags if abs(pos[t][1]) < 1e-9}
    return inner, pos


def restore_collar_first_layer(gmsh, surfaces: list[int], inner: set[int],
                               pre_pos: dict[int, np.ndarray]) -> int:
    """Re-space root-collar nodes along their post-snap direction so every BL chain keeps its pre-snap cumulative height."""
    tags: set[int] = set()
    for t in surfaces:
        tags.update(int(x) for x in gmsh.model.mesh.getNodes(2, t, True)[0])
    outer = sorted(tags - inner)
    if not outer or not inner:
        return 0
    inner_tags = sorted(inner)
    base_pts = np.array([pre_pos[t] for t in inner_tags])
    from scipy.spatial import cKDTree
    tree = cKDTree(base_pts)
    raw = gmsh.model.mesh.getNodes()
    coord = {int(t): np.asarray(raw[1]).reshape(-1, 3)[i]
             for i, t in enumerate(raw[0])}
    fixed = 0
    for t in outer:
        _, idx = tree.query(pre_pos[t], k=1)
        base = inner_tags[int(idx)]
        height = float(np.linalg.norm(pre_pos[t] - pre_pos[base]))
        delta = coord[t] - coord[base]
        length = float(np.linalg.norm(delta))
        if height > 1e-12 and length > 1e-12 and abs(length - height) > 1e-12:
            new_pos = coord[base] + delta / length * height
            gmsh.model.mesh.setNode(t, new_pos.tolist(), [])
            fixed += 1
    return fixed


def snap_collar_to_symmetry(gmsh, surfaces: list[int]) -> int:
    """Project the symmetry-plane faces onto y=0 (clip the BL at y=0)."""
    nodes = gmsh.model.mesh.getNodes()
    tag_to_row = {int(t): i for i, t in enumerate(nodes[0])}
    coords = np.asarray(nodes[1]).reshape(-1, 3)
    moved = 0
    seen: set[int] = set()
    for t in surfaces:
        for raw in gmsh.model.mesh.getNodes(2, t, True)[0]:
            tag = int(raw)
            if tag in seen:
                continue
            seen.add(tag)
            row = tag_to_row[tag]
            if abs(coords[row, 1]) > 0.0:
                gmsh.model.mesh.setNode(
                    tag, [float(coords[row, 0]), 0.0, float(coords[row, 2])], [])
                moved += 1
    return moved


def assign_physical_groups(gmsh, names: dict, new_front: int,
                           root_collar: list[int], trailing_collar: list[int],
                           fsi_surfaces: list[int],
                           all_volumes: list[int]) -> None:
    groups = {k: list(v) for k, v in names.items()}
    groups["FRONT"] = [new_front] + list(root_collar)
    groups["FSI"] = list(fsi_surfaces) + list(trailing_collar)
    for name, tags in groups.items():
        tag = gmsh.model.addPhysicalGroup(2, tags)
        gmsh.model.setPhysicalName(2, tag, name)
    vtag = gmsh.model.addPhysicalGroup(3, all_volumes)
    gmsh.model.setPhysicalName(3, vtag, "FLUID")


def _element_blocks(gmsh) -> dict:
    return {"tet": gmsh.model.mesh.getElementsByType(4),
            "prism": gmsh.model.mesh.getElementsByType(6),
            "pyramid": gmsh.model.mesh.getElementsByType(7)}


def _positive_tet(conn: tuple, xyz: dict) -> tuple:
    """Return the 4 nodes in an orientation with positive signed volume."""
    a, b, c, d = (xyz[t] for t in conn)
    if float(np.dot(np.cross(b - a, c - a), d - a)) < 0.0:
        return (conn[0], conn[2], conn[1], conn[3])
    return conn


def _split_pyramids(blocks: dict, xyz: dict) -> list[tuple]:
    """Split every quad-base pyramid into two C3D4 tets (CalculiX has no C3D5)."""
    enodes = blocks["pyramid"][1]
    out: list[tuple] = []
    for row in np.asarray(enodes).reshape(-1, 5):
        q0, q1, q2, q3, apex = (int(t) for t in row)
        out.append(_positive_tet((q0, q1, q2, apex), xyz))
        out.append(_positive_tet((q0, q2, q3, apex), xyz))
    return out


def write_calculix(gmsh, path: Path, bl_layers: int, first: float,
                   ratio: float) -> dict:
    """C3D4 + C3D6 deck with the sibling NSET names."""
    node_tags, node_xyz, _ = gmsh.model.mesh.getNodes()
    order = np.argsort(node_tags)
    node_tags = node_tags[order]
    node_xyz = np.asarray(node_xyz).reshape(-1, 3)[order]
    new_id = {int(t): i + 1 for i, t in enumerate(node_tags)}
    xyz = {int(t): node_xyz[i] for i, t in enumerate(node_tags)}

    blocks = _element_blocks(gmsh)
    pyramid_tets = _split_pyramids(blocks, xyz)

    nset: dict[str, set[int]] = {n: set() for n in PATCH_TYPES}
    for dim, ptag in gmsh.model.getPhysicalGroups(2):
        name = gmsh.model.getPhysicalName(2, ptag)
        for e in gmsh.model.getEntitiesForPhysicalGroup(2, ptag):
            tags = gmsh.model.mesh.getNodes(2, e, True)[0]
            nset.setdefault(name, set()).update(int(t) for t in tags)
    missing = [n for n in PATCH_TYPES if not nset.get(n)]
    if missing:
        raise MeshError(f"empty NSET surfaces: {missing}")

    header = [
        "** ---------------------------------------------------------------------------\n",
        "** AGARD 445.6 body-fitted fluid solid deck (C3D4 core + C3D6 boundary layer)\n",
        f"** BL: {bl_layers} layers, first {first:g} m, ratio {ratio:g}\n",
        "** generated by case/AGARD/Mesh/build_fluid_mesh.py\n",
        "** ---------------------------------------------------------------------------\n",
    ]
    tet_rows = [tuple(int(t) for t in row)
                for row in np.asarray(blocks["tet"][1]).reshape(-1, 4)]
    tet_rows += pyramid_tets
    prism_rows = [tuple(int(t) for t in row)
                  for row in np.asarray(blocks["prism"][1]).reshape(-1, 6)]
    with path.open("w", encoding="ascii") as fh:
        fh.writelines(header)
        fh.write("*NODE\n")
        for i in range(node_xyz.shape[0]):
            x, y, z = node_xyz[i]
            fh.write(f"{i + 1},{x:.9g},{y:.9g},{z:.9g}\n")
        elem_id = 0
        for etype, rows in (("C3D4", tet_rows), ("C3D6", prism_rows)):
            fh.write(f"*ELEMENT, TYPE={etype}, ELSET=FAR_PART_ELSET\n")
            for row in rows:
                elem_id += 1
                conn = [new_id[t] for t in row]
                fh.write(f"{elem_id}," + ",".join(str(v) for v in conn) + "\n")
        _write_nset(fh, "FAR_PART_NSET",
                    [new_id[int(t)] for t in node_tags])
        for name in PATCH_TYPES:
            _write_nset(fh, f"{name}_NSET",
                        sorted(new_id[t] for t in nset[name]))
    return {"nodes": int(node_xyz.shape[0]), "c3d4": len(tet_rows),
            "c3d6": len(prism_rows), "pyramids_split": len(pyramid_tets) // 2,
            "nset_sizes": {n: len(nset[n]) for n in PATCH_TYPES}}


def _write_nset(fh, name: str, values: list[int]) -> None:
    fh.write(f"*NSET, NSET={name}\n")
    for off in range(0, len(values), 16):
        fh.write(",".join(str(v) for v in values[off:off + 16]) + "\n")


def write_openfoam(out_dir: Path, msh: Path, log) -> dict:
    """gmshToFoam + patch retype + minimal system dicts."""
    case = out_dir / "openfoam"
    if case.exists():
        shutil.rmtree(case)
    (case / "system").mkdir(parents=True)
    shutil.copy(msh, case / "fluid_mesh.msh")
    (case / "system" / "controlDict").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }\n"
        "application rhoPimpleFoam;\nstartFrom startTime; startTime 0;\n"
        "stopAt endTime; endTime 1;\ndeltaT 1e-5;\n"
        "writeControl timeStep; writeInterval 1000;\nrunTimeModifiable true;\n")
    (case / "system" / "fvSchemes").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }\n"
        "ddtSchemes { default Euler; }\ngradSchemes { default Gauss linear; }\n"
        "divSchemes { default none; div(phi,U) Gauss upwind; }\n"
        "laplacianSchemes { default Gauss linear corrected; }\n"
        "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    (case / "system" / "fvSolution").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }\n"
        "solvers { p { solver GAMG; tolerance 1e-7; relTol 0.01; } }\n"
        "PIMPLE { nOuterCorrectors 1; nCorrectors 2; }\n")
    log_path = out_dir / "gmshToFoam.log"
    with log_path.open("w") as lf:
        rc = subprocess.run(["gmshToFoam", "fluid_mesh.msh"], cwd=case,
                            stdout=lf, stderr=subprocess.STDOUT).returncode
    log(f"gmshToFoam exit {rc} -> {log_path}")
    if rc != 0:
        raise MeshError("gmshToFoam failed; see gmshToFoam.log")
    boundary = case / "constant" / "polyMesh" / "boundary"
    text = boundary.read_text()
    for name, ptype in PATCH_TYPES.items():
        text = re.sub(rf"(\b{name}\s*\{{[^}}]*?type\s+)\w+(\s*;)",
                      rf"\1{ptype}\2", text, count=1)
    boundary.write_text(text)
    if "defaultFaces" in text:
        raise MeshError("gmshToFoam created defaultFaces; some boundary was "
                        "left without a physical group")
    return {"case": str(case), "gmshToFoam_log": str(log_path)}


def compare_coordinates(out_dir: Path, gmsh) -> float:
    """Max nearest-neighbour distance between OpenFOAM points and gmsh nodes."""
    pm = out_dir / "openfoam" / "constant" / "polyMesh" / "points"
    if not pm.is_file():
        return float("nan")
    lines = pm.read_text().splitlines()
    idx = next(i for i, ln in enumerate(lines) if ln.strip().isdigit())
    n = int(lines[idx].strip())
    pts = np.empty((n, 3))
    for i in range(n):
        s = lines[idx + 2 + i].strip().rstrip(")").lstrip("(")
        pts[i] = [float(v) for v in s.split()]
    node_xyz = gmsh.model.mesh.getNodes()[1]
    q = np.asarray(node_xyz).reshape(-1, 3)
    try:
        from scipy.spatial import cKDTree
        d, _ = cKDTree(q).query(pts, k=1)
        return float(np.max(d))
    except ImportError:
        return float("nan")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--step", type=Path, default=STEP_DEFAULT)
    p.add_argument("--heal-te", type=float, default=0.0,
                   help="OCC heal tolerance to sharpen the blunt TE (0 = off)")
    p.add_argument("--bl-layers", type=int, default=35,
                   help="boundary-layer prism layers (default 35)")
    p.add_argument("--bl-first", type=float, default=1.51e-5,
                   help="first-layer height m, y+~1 target (default 1.51e-5)")
    p.add_argument("--bl-ratio", type=float, default=1.15,
                   help="BL growth ratio (default 1.15)")
    p.add_argument("--wing-size", type=float, default=0.012,
                   help="surface element size at the wing m (default 0.012)")
    p.add_argument("--mesh-size-min", type=float, default=0.004,
                   help="global minimum element size m (default 0.004)")
    p.add_argument("--smoothing", type=int, default=10,
                   help="gmsh Mesh.Smoothing steps (default 10)")
    p.add_argument("--far-size", type=float, default=1.0,
                   help="far-field element size m (default 1.0)")
    p.add_argument("--dist-min", type=float, default=0.01,
                   help="distance where refinement starts m (default 0.01)")
    p.add_argument("--dist-max", type=float, default=2.0,
                   help="distance where far size is reached m (default 2.0)")
    p.add_argument("--algo3d", type=int, default=1,
                   help="gmsh Mesh.Algorithm3D (default 1 Delaunay)")
    p.add_argument("--te-bl", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="extrude BL from the thin blunt trailing-edge face "
                        "(default on; disable only for diagnostics)")
    p.add_argument("--clip", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="project the BL root collar onto y=0 (symmetry clip)")
    p.add_argument("--optimize", action=argparse.BooleanOptionalAction,
                   default=True, help="Netgen-relax after the symmetry clip")
    p.add_argument("--opt-passes", type=int, default=3,
                   help="number of Netgen+Laplace relaxation passes (default 3)")
    p.add_argument("--openfoam", action=argparse.BooleanOptionalAction,
                   default=True, help="export/convert the OpenFOAM polyMesh")
    p.add_argument("--check-mesh", action=argparse.BooleanOptionalAction,
                   default=True, help="run checkMesh (needs OpenFOAM env)")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    report: dict = {"params": {k: (str(v) if isinstance(v, Path) else v)
                               for k, v in vars(args).items()}}

    def log(msg: str) -> None:
        _log(msg)
        with (out / "build.log").open("a") as fh:
            fh.write(msg + "\n")

    (out / "build.log").write_text("")

    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", args.algo3d)
        gmsh.option.setNumber("Mesh.RandomFactor", 1e-9)
        gmsh.model.add("agard4456_fluid")

        geom = build_geometry(gmsh, args.step, args.heal_te)
        try:
            gmsh.write(str(out / "fluid_mesh.geo_unrolled"))
        except Exception as exc:  # noqa: BLE001 - snapshot is best-effort
            log(f"geo_unrolled snapshot skipped: {exc}")
        report["wing_bbox_m"] = geom["wing_bbox"]
        log(f"wing bbox {[round(v, 5) for v in geom['wing_bbox']]}")
        te = set() if args.te_bl else set(geom["te_surfs"])
        bl_base = [t for t in geom["fsi"] if t not in te]
        log(f"BL base {len(bl_base)} wetted surfaces; "
            f"trailing-edge faces {sorted(te)} left tetrahedral")

        bl = extrude_boundary_layer(
            gmsh, bl_base, args.bl_layers, args.bl_first, args.bl_ratio)
        for tag in bl["collar"]:
            gmsh.model.mesh.setRecombine(2, tag, 0)
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        report["bl"] = {"layers": args.bl_layers, "first_m": args.bl_first,
                        "ratio": args.bl_ratio, "thickness_m": bl["thickness"]}
        log(f"BL: {args.bl_layers} layers, first {args.bl_first:g}, ratio "
            f"{args.bl_ratio:g}, thickness {bl['thickness']*1e3:.3f} mm")

        vols = build_fluid_volume(
            gmsh, geom, bl, list(geom["te_surfs"]) if not args.te_bl else [])
        root_collar = vols["root_collar"]
        trailing_collar = vols["te_collar"]
        log(f"collars root={len(root_collar)} te={len(trailing_collar)}")
        names = classify_outer(gmsh, vols["others"])
        set_size_field(gmsh, geom["fsi"], args)
        assign_physical_groups(gmsh, names, vols["new_front"], root_collar,
                               trailing_collar, geom["fsi"],
                               vols["all_volumes"])
        gmsh.model.geo.synchronize()

        log("generating 3D mesh ...")
        tgen = time.time()
        gmsh.model.mesh.generate(3)
        report["generate_s"] = round(time.time() - tgen, 1)
        log(f"generated in {report['generate_s']} s")

        front_faces = [vols["new_front"]] + list(root_collar)
        clip_moved = 0
        if args.clip:
            inner, pre_pos = collar_inner_nodes(gmsh, root_collar)
            clip_moved = snap_collar_to_symmetry(gmsh, front_faces)
            log(f"symmetry clip: snapped {clip_moved} FRONT nodes to y=0")
            if args.optimize:
                for _ in range(max(1, args.opt_passes)):
                    gmsh.model.mesh.optimize("Netgen")
                    try:
                        gmsh.model.mesh.optimize("Laplace", niter=3)
                    except Exception:  # noqa: BLE001 - Laplace is best-effort
                        pass
                re_snapped = snap_collar_to_symmetry(gmsh, front_faces)
                log(f"relaxation ({args.opt_passes} passes) done; "
                    f"re-snapped {re_snapped} FRONT nodes")
            stretched = restore_collar_first_layer(
                gmsh, root_collar, inner, pre_pos)
            log(f"restored {stretched} root-collar nodes to their BL "
                "cumulative heights")
        report["collar_nodes_snapped"] = clip_moved

        msh = out / "fluid_mesh.msh"
        gmsh.write(str(msh))
        counts = write_calculix(gmsh, out / "fluid_mesh.inp",
                                args.bl_layers, args.bl_first, args.bl_ratio)
        report["mesh"] = counts
        log(f"mesh: nodes={counts['nodes']} C3D4={counts['c3d4']} "
            f"C3D6={counts['c3d6']}")

        if args.openfoam:
            report["openfoam"] = write_openfoam(out, msh, log)
            report["coord_max_diff_m"] = compare_coordinates(out, gmsh)
            log(f"OpenFOAM vs gmsh node max|diff| = "
                f"{report['coord_max_diff_m']:.3e} m")
            if args.check_mesh:
                cmlog = out / "checkMesh.log"
                with cmlog.open("w") as lf:
                    subprocess.run(["checkMesh"], cwd=out / "openfoam", stdout=lf,
                                   stderr=subprocess.STDOUT)
                report["checkMesh_log"] = str(cmlog)
                log(f"checkMesh -> {cmlog}")

        report["elapsed_s"] = round(time.time() - started, 1)
        (out / "mesh_report.json").write_text(json.dumps(report, indent=2) + "\n")
        write_readme(out, args, report)
        log(f"done in {report['elapsed_s']} s")
    finally:
        gmsh.finalize()
    return 0


def write_readme(out: Path, args, report: dict) -> None:
    r = report
    bl = r.get("bl", {})
    text = f"""# AGARD 445.6 body-fitted fluid mesh

Generated by `case/AGARD/Mesh/build_fluid_mesh.py`.

## Domain (exactly the existing far field)

`y in [0, 7]`; `x <= 0` half-cylinder r=5 about the y axis; `x >= 0` box
`x in [0, 6]`, `|z| <= 5`.  Patches: FRONT (y=0, symmetryPlane), BACK (y=7),
IN (cylinder), OUT (x=6), TOP/BOTTOM (z=+/-5), FSI (wing wall).

## Parameters

| parameter | value |
| --- | --- |
| BL layers | {args.bl_layers} |
| BL first layer | {args.bl_first:g} m |
| BL ratio | {args.bl_ratio:g} |
| BL total thickness | {bl.get('thickness_m', float('nan')):.5g} m |
| wing surface size | {args.wing_size:g} m |
| far-field size | {args.far_size:g} m |
| refinement transition | {args.dist_min:g} -> {args.dist_max:g} m |
| 3D algorithm | {args.algo3d} |
| symmetry-clipped collar nodes | {r.get('collar_nodes_snapped', '?')} |

## Measured

- nodes {r.get('mesh', {}).get('nodes', '?')},
  C3D4 {r.get('mesh', {}).get('c3d4', '?')},
  C3D6 {r.get('mesh', {}).get('c3d6', '?')}
- OpenFOAM vs gmsh node max|diff| = {r.get('coord_max_diff_m', float('nan')):.3e} m
- checkMesh log: `checkMesh.log`

The wing root sits on the symmetry plane; the extruded BL collar (which leans
across y=0 because the swept wing's root-edge normal has a y component) is
projected onto y=0 and Netgen-relaxed, so the domain stays exactly y in [0, 7]
and FRONT remains a planar symmetryPlane.
"""
    (out / "README.md").write_text(text)


if __name__ == "__main__":
    sys.exit(main())
