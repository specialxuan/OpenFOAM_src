#!/usr/bin/env python3
"""Convert a CalculiX/Abaqus-format mesh (mixed C3D8 hex + C3D6 wedge) into an
OpenFOAM ``polyMesh`` (``points``/``faces``/``owner``/``neighbour``/``boundary``).

The converter is written for the AGARD 445.6 soft-wing FDM-AMR workflow: the
far-field pseudo-elastic body is turned into the fluid ``polyMesh`` used by
``myRhoPimpleFoam``.  It is a standalone script, not part of any ``wmake``
build.

Design contract (verified against the input, not guessed):

* ``*NODE`` ids may contain gaps (the AGARD Fluid7 deck has 280,770 present ids
  spanning ``1..285460``).  The OpenFOAM point index is the rank of the node id
  in ascending present-id order, which matches the row order of the
  ``FluidNodeCoor.csv`` cross-check; connectivity that references an absent id
  aborts loudly.
* A deck may mix element types.  Supported: ``C3D8`` (8-node hex, 6 quad faces)
  and ``C3D6`` (6-node wedge, Abaqus S1..S5: 2 triangles + 3 quads).  Element
  ids must be globally unique; cells are re-ordered by ascending element id, so
  cell ``i`` is the ``i``-th element by id, not ``id - 1``.
* Faces are built from the per-type face table and each face is re-oriented so
  its Newell normal points outward from the creating cell (cell-centroid test).
* Faces are de-duplicated by a canonical (sorted) node tuple padded to four
  entries with ``-1``.  The first (lowest-id) cell that creates a face is its
  ``owner``; the second is its ``neighbour``.
* Internal faces come first, upper-triangular ordered by ``(owner, neighbour)``;
  boundary faces follow in contiguous per-patch ranges.
* Boundary faces (used by exactly one cell) are classified by which boundary
  ``*NSET`` contains **all** of their real vertices (3 or 4).  Any unclassified
  or multiply-classified boundary face aborts the conversion loudly.
* ``FAR_PART_NSET`` is the full node set and is never used as a patch.

Usage
-----

    python3 inp_to_polymesh.py --input Fluid7.inp \
        --out-dir /path/to/case

Add ``--self-check`` to re-read the files just written and compare counts
against the parsed input.  ``--help`` documents every option.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import zipfile
from array import array
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Element type codes.  A deck may mix types; every element records its code so
# the face builder can pick the right per-type face table and node count.
_TYPE_C3D8 = 1
_TYPE_C3D6 = 2
_ELEM_TYPE_NAME: dict[int, str] = {_TYPE_C3D8: "C3D8", _TYPE_C3D6: "C3D6"}
_ELEM_NNODES: dict[int, int] = {_TYPE_C3D8: 8, _TYPE_C3D6: 6}
_MAX_ELEM_NODES = 8
_MAX_FACE_NODES = 4

# Per-element-type face topology as ``(local 0-based node tuple, face size)``.
# C3D8 gives the six standard hex quads; C3D6 (wedge/prism) gives Abaqus S1..S5:
# two triangular caps plus three quad sides.  Winding is corrected per cell by
# the outward-centroid test, so only the node *set* of each face matters here.
_FACE_SPECS: dict[int, tuple[tuple[tuple[int, ...], int], ...]] = {
    _TYPE_C3D8: (
        ((0, 1, 2, 3), 4),  # (1,2,3,4)
        ((4, 7, 6, 5), 4),  # (5,8,7,6)
        ((0, 4, 5, 1), 4),  # (1,5,6,2)
        ((1, 5, 6, 2), 4),  # (2,6,7,3)
        ((2, 6, 7, 3), 4),  # (3,7,8,4)
        ((3, 7, 4, 0), 4),  # (4,8,5,1)
    ),
    _TYPE_C3D6: (
        ((0, 1, 2), 3),     # S1 (1,2,3)
        ((3, 4, 5), 3),     # S2 (4,5,6)
        ((0, 1, 4, 3), 4),  # S3 (1,2,5,4)
        ((1, 2, 5, 4), 4),  # S4 (2,3,6,5)
        ((2, 0, 3, 5), 4),  # S5 (3,1,4,6)
    ),
}

# Patch write order and types.  `FSI` is the wing surface on the fluid side.
_PATCH_ORDER: tuple[str, ...] = ("FSI", "FRONT", "BACK", "IN", "OUT", "TOP", "BOTTOM")
_PATCH_TYPE: dict[str, str] = {
    "FSI": "wall",
    "FRONT": "patch",
    "BACK": "patch",
    "IN": "patch",
    "OUT": "patch",
    "TOP": "patch",
    "BOTTOM": "patch",
}

# The far-field NSET is the full node set; it is a sanity check, not a patch.
_FULL_SET = "FAR_PART_NSET"

_CHUNK = 20000
_COORD_CSV_DEFAULT = "/root/Workspace/AGARD/Mode/Fluid8/mode_soft/FluidNodeCoor.csv"


# ---------------------------------------------------------------------------
# Input reading
# ---------------------------------------------------------------------------


def _iter_input_lines(path: Path):
    """Yield ASCII lines from a plain, ``.gz`` or single-member ``.zip`` deck."""
    suffix = path.suffix.lower()
    if suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            members = archive.namelist()
            if len(members) != 1:
                raise SystemExit(
                    f"[inp_to_polymesh] ERROR: {path} holds {len(members)} members; "
                    "a single-member zip is required"
                )
            with archive.open(members[0]) as raw:
                for raw_line in raw:
                    yield raw_line.decode("ascii")
    elif suffix == ".gz":
        with gzip.open(path, "rt", encoding="ascii") as handle:
            yield from handle
    else:
        with open(path, "r", encoding="ascii") as handle:
            yield from handle


def _nset_name(header: str) -> str | None:
    for field in header.split(",")[1:]:
        if "=" in field:
            key, value = field.split("=", 1)
            if key.strip().upper() == "NSET":
                return value.strip()
    return None


def _element_type(header: str) -> int | None:
    upper = header.upper()
    for name, code in (("C3D8", _TYPE_C3D8), ("C3D6", _TYPE_C3D6)):
        if name in upper:
            return code
    return None


def _present_id_index(present_ids: np.ndarray, query_ids: np.ndarray) -> np.ndarray:
    """Map node ids to their rank in ascending ``present_ids``.

    Tolerates gaps in the id space.  Aborts if any query id is absent.
    """
    idx = np.searchsorted(present_ids, query_ids)
    idx_clipped = np.minimum(idx, present_ids.size - 1)
    absent = present_ids[idx_clipped] != query_ids
    if np.any(absent):
        missing = np.unique(query_ids[absent])
        preview = ", ".join(str(int(v)) for v in missing[:10])
        raise SystemExit(
            "[inp_to_polymesh] ERROR: connectivity references "
            f"{missing.size} node id(s) absent from *NODE: [{preview}"
            f"{', ...' if missing.size > 10 else ''}]"
        )
    return idx


def parse_inp(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Stream the deck and return ``(points, conn, elem_type, nset_membership)``.

    ``points`` has shape ``(nNodes, 3)`` in ascending present-node-id order;
    ``conn`` has shape ``(nCells, 8)`` padded with ``-1`` and holds 0-based point
    indices (cell ``i`` is the ``i``-th element in ascending element-id order);
    ``elem_type`` has shape ``(nCells,)`` with a ``_TYPE_*`` code per cell;
    ``nset_membership`` maps each NSET name to a boolean array of length
    ``nNodes`` indexed by point position.
    """
    node_ids = array("q")
    node_xyz = array("d")
    elem_ids = array("q")
    elem_types = array("q")
    elem_conn = array("q")
    nsets: dict[str, array] = {}

    state = "other"
    nset_name: str | None = None
    elem_type: int | None = None

    for line in _iter_input_lines(path):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("*"):
            upper = stripped.upper()
            if upper.startswith("*NODE"):
                state = "node"
            elif upper.startswith("*ELEMENT"):
                code = _element_type(stripped)
                if code is None:
                    raise SystemExit(
                        "[inp_to_polymesh] ERROR: unsupported element type in "
                        f"header {stripped!r} (only C3D8 and C3D6 are handled)"
                    )
                state = "elem"
                elem_type = code
            elif upper.startswith("*NSET"):
                name = _nset_name(stripped)
                if name is None:
                    raise SystemExit(
                        f"[inp_to_polymesh] ERROR: malformed NSET header: {stripped!r}"
                    )
                state = "nset"
                nset_name = name
                nsets.setdefault(name, array("q"))
            else:
                state = "other"
            continue

        if state == "node":
            fields = stripped.split(",")
            node_ids.append(int(fields[0]))
            node_xyz.append(float(fields[1]))
            node_xyz.append(float(fields[2]))
            node_xyz.append(float(fields[3]))
        elif state == "elem":
            if elem_type is None:
                raise SystemExit(
                    "[inp_to_polymesh] ERROR: element data without a preceding "
                    "*ELEMENT header"
                )
            fields = stripped.split(",")
            elem_ids.append(int(fields[0]))
            elem_types.append(elem_type)
            nodes = [int(value) for value in fields[1:] if value.strip()]
            expected = _ELEM_NNODES[elem_type]
            if len(nodes) != expected:
                raise SystemExit(
                    "[inp_to_polymesh] ERROR: "
                    f"{_ELEM_TYPE_NAME[elem_type]} element "
                    f"{fields[0].strip()} has {len(nodes)} nodes, expected "
                    f"{expected}"
                )
            nodes.extend([-1] * (_MAX_ELEM_NODES - len(nodes)))
            elem_conn.extend(nodes)
        elif state == "nset" and nset_name is not None:
            nsets[nset_name].extend(
                int(token) for token in stripped.replace(",", " ").split()
            )

    nid = np.frombuffer(node_ids, dtype=np.int64)
    xyz = np.frombuffer(node_xyz, dtype=np.float64).reshape(-1, 3)
    eid = np.frombuffer(elem_ids, dtype=np.int64)
    etype = np.frombuffer(elem_types, dtype=np.int64)
    conn = np.frombuffer(elem_conn, dtype=np.int64).reshape(-1, _MAX_ELEM_NODES)

    n_nodes = nid.size
    if n_nodes == 0:
        raise SystemExit("[inp_to_polymesh] ERROR: no *NODE data found")
    if np.unique(nid).size != n_nodes:
        raise SystemExit(
            "[inp_to_polymesh] ERROR: duplicate node ids present "
            f"(n={n_nodes}, unique={np.unique(nid).size})"
        )
    if nid.min() < 1:
        raise SystemExit(
            f"[inp_to_polymesh] ERROR: non-positive node id {int(nid.min())}"
        )
    # Node ids may contain gaps (the AGARD Fluid7 deck has 4,690 of them).  The
    # OpenFOAM point index is the rank of the node id in ascending present-id
    # order, which matches the row order of FluidNodeCoor.csv.
    present_ids = np.sort(nid)
    points = np.empty((n_nodes, 3), dtype=np.float64)
    points[np.searchsorted(present_ids, nid)] = xyz
    n_gaps = int(present_ids[-1] - n_nodes)
    print(
        "[inp_to_polymesh] node ids present "
        f"n={n_nodes}, range [1..{int(present_ids[-1])}], gaps={n_gaps} "
        "-> point index by ascending present id"
    )

    n_cells = eid.size
    if n_cells == 0:
        raise SystemExit("[inp_to_polymesh] ERROR: no *ELEMENT data found")
    if np.unique(eid).size != n_cells:
        raise SystemExit(
            "[inp_to_polymesh] ERROR: duplicate element ids present "
            f"(n={n_cells}, unique={np.unique(eid).size})"
        )
    order = np.argsort(eid, kind="stable")
    conn_ordered = conn[order]
    etype_ordered = etype[order]

    valid = conn_ordered >= 0
    mapped = conn_ordered.copy()
    mapped[valid] = _present_id_index(present_ids, conn_ordered[valid])

    counts = {
        name: int(np.count_nonzero(etype_ordered == code))
        for code, name in _ELEM_TYPE_NAME.items()
    }
    counts_str = ", ".join(f"{name}={n}" for name, n in counts.items())
    print(
        "[inp_to_polymesh] element ids unique "
        f"[{int(eid.min())}..{int(eid.max())}] ({counts_str}); "
        "cell index by ascending id"
    )

    membership: dict[str, np.ndarray] = {}
    for name, ids in nsets.items():
        arr = np.frombuffer(ids, dtype=np.int64)
        mask = np.zeros(n_nodes, dtype=bool)
        if arr.size:
            unique = np.unique(arr)
            if unique.size != arr.size:
                print(
                    f"[inp_to_polymesh] note: NSET {name} has duplicate ids; "
                    "de-duplicated for membership"
                )
            if unique.min() < 1:
                raise SystemExit(
                    f"[inp_to_polymesh] ERROR: NSET {name} has a non-positive id"
                )
            mask[_present_id_index(present_ids, unique)] = True
        membership[name] = mask

    return points, mapped, etype_ordered, membership


# ---------------------------------------------------------------------------
# Face construction
# ---------------------------------------------------------------------------


def _openfoam_cell_volumes(
    points: np.ndarray, conn: np.ndarray, oriented_faces: list[np.ndarray]
) -> np.ndarray:
    """Replicate OpenFOAM ``cellModel::mag`` for one block of cells.

    ``oriented_faces`` holds one ``(nCells, k)`` point-index array per face slot;
    ``k`` may differ between slots so triangular and quad faces coexist.  The
    volume is the sum of the face pyramids from the cell centre.  The
    two-triangle divergence formula is NOT used because it underflows on the
    strongly warped thin cells of this mesh, where OpenFOAM's per-face area
    normal keeps the volume positive.
    """
    cell_centroid = points[conn].mean(axis=1)
    volume = np.zeros(conn.shape[0], dtype=np.float64)
    for face in oriented_faces:
        pts = points[face]  # (nCells, k, 3)
        centre_point = pts.mean(axis=1, keepdims=True)
        nxt = np.roll(pts, -1, axis=1)
        area_normal = 0.5 * np.cross(nxt - pts, centre_point - pts).sum(axis=1)
        triangle_centres = pts + nxt + centre_point
        tri_area = np.linalg.norm(np.cross(pts - centre_point, nxt - centre_point), axis=2)
        sum_area = tri_area.sum(axis=1)
        sum_area_centre = (tri_area[:, :, None] * triangle_centres).sum(axis=1)
        face_centre = np.where(
            (sum_area > 1e-300)[:, None],
            sum_area_centre / (3.0 * sum_area)[:, None],
            centre_point[:, 0, :],
        )
        volume += np.einsum("ij,ij->i", area_normal, face_centre - cell_centroid) / 3.0
    return volume


def _orient_faces(
    points: np.ndarray, conn: np.ndarray, spec: tuple[tuple[int, ...], int]
) -> np.ndarray:
    """Return the ``spec`` face of every cell, wound outward from the cell."""
    local_nodes, _size = spec
    face = conn[:, list(local_nodes)].copy()  # (nCells, k)
    pts = points[face]
    nxt = np.roll(pts, -1, axis=1)
    nx = np.sum((pts[:, :, 1] - nxt[:, :, 1]) * (pts[:, :, 2] + nxt[:, :, 2]), axis=1)
    ny = np.sum((pts[:, :, 2] - nxt[:, :, 2]) * (pts[:, :, 0] + nxt[:, :, 0]), axis=1)
    nz = np.sum((pts[:, :, 0] - nxt[:, :, 0]) * (pts[:, :, 1] + nxt[:, :, 1]), axis=1)
    cell_centroid = points[conn].mean(axis=1)
    face_centroid = pts.mean(axis=1)
    outward = (
        nx * (face_centroid[:, 0] - cell_centroid[:, 0])
        + ny * (face_centroid[:, 1] - cell_centroid[:, 1])
        + nz * (face_centroid[:, 2] - cell_centroid[:, 2])
    )
    flip = outward < 0.0
    if flip.any():
        face[flip] = face[flip][:, ::-1]
    return face


def build_faces(
    points: np.ndarray, conn: np.ndarray, etype: np.ndarray
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Return ``(uniq_sorted, uniq_size, owner_oriented, owner_cell,
    neighbour_cell, is_internal, cell_volume)``.

    ``uniq_sorted`` is the sorted node key of every unique face padded to
    ``_MAX_FACE_NODES`` with ``-1`` (for membership tests); ``uniq_size`` is its
    true vertex count; ``owner_oriented`` is that face wound outward from
    ``owner_cell``; ``neighbour_cell`` is ``-1`` for boundary faces.
    """
    n_cells = conn.shape[0]
    faces_blocks: list[np.ndarray] = []
    size_blocks: list[np.ndarray] = []
    cell_blocks: list[np.ndarray] = []
    cell_volume = np.empty(n_cells, dtype=np.float64)

    for code in (_TYPE_C3D8, _TYPE_C3D6):
        cells = np.nonzero(etype == code)[0]
        if cells.size == 0:
            continue
        conn_t = conn[cells][:, : _ELEM_NNODES[code]]
        specs = _FACE_SPECS[code]
        oriented: list[np.ndarray] = []
        padded = np.full(
            (cells.size, len(specs), _MAX_FACE_NODES), -1, dtype=np.int64
        )
        for slot, spec in enumerate(specs):
            size = spec[1]
            face = _orient_faces(points, conn_t, spec)
            padded[:, slot, :size] = face
            oriented.append(face)
        faces_blocks.append(padded.reshape(-1, _MAX_FACE_NODES))
        size_blocks.append(np.tile([spec[1] for spec in specs], cells.size))
        cell_blocks.append(np.repeat(cells, len(specs)))
        cell_volume[cells] = _openfoam_cell_volumes(points, conn_t, oriented)

    if np.any(cell_volume <= 0.0):
        bad = int(np.count_nonzero(cell_volume <= 0.0))
        raise SystemExit(
            f"[inp_to_polymesh] ERROR: {bad} non-positive cell volumes "
            f"(OpenFOAM pyramid definition, min={cell_volume.min():.6e})"
        )

    faces = np.concatenate(faces_blocks, axis=0)
    face_size = np.concatenate(size_blocks)
    cell_of = np.concatenate(cell_blocks)
    # Face de-duplication takes the first occurrence as owner, so the flattened
    # faces MUST be cell-major with ascending cell index.  The per-type blocks
    # above are already cell-major; this stable sort merges them (and keeps each
    # cell's slot order) even when the element-id order interleaves types.
    face_order = np.argsort(cell_of, kind="stable")
    faces = faces[face_order]
    face_size = face_size[face_order]
    cell_of = cell_of[face_order]

    key = np.sort(faces, axis=1)
    key_view = (
        np.ascontiguousarray(key.astype(np.int64))
        .view(
            [
                ("n0", "<i8"),
                ("n1", "<i8"),
                ("n2", "<i8"),
                ("n3", "<i8"),
            ]
        )
        .reshape(-1)
    )
    _, first_idx, inverse, counts = np.unique(
        key_view, return_index=True, return_inverse=True, return_counts=True
    )
    if counts.max() > 2:
        hist = {int(k): int(v) for k, v in zip(*np.unique(counts, return_counts=True))}
        raise SystemExit(
            "[inp_to_polymesh] ERROR: non-manifold faces with >2 cells: "
            f"occurrence histogram {hist}"
        )

    uniq_sorted = key[first_idx]
    uniq_size = face_size[first_idx]
    owner_oriented = faces[first_idx]
    owner_cell = cell_of[first_idx]
    is_internal = counts == 2

    order = np.argsort(inverse, kind="stable")
    starts = np.zeros(counts.size, dtype=np.int64)
    np.cumsum(counts[:-1], out=starts[1:])
    if not np.array_equal(order[starts], first_idx):
        raise SystemExit(
            "[inp_to_polymesh] ERROR: face owner bookkeeping mismatch "
            "(first occurrence is not the lowest-id cell)"
        )
    second_occ = order[np.minimum(starts + 1, order.size - 1)]
    neighbour_cell = np.full(counts.size, -1, dtype=np.int64)
    neighbour_cell[is_internal] = cell_of[second_occ[is_internal]]
    if np.any(owner_cell[is_internal] >= neighbour_cell[is_internal]):
        raise SystemExit(
            "[inp_to_polymesh] ERROR: internal face owner/neighbour are not "
            "strictly ordered by cell id"
        )

    repeated = np.zeros(key.shape[0], dtype=bool)
    for slot in range(_MAX_FACE_NODES - 1):
        both_real = slot >= (_MAX_FACE_NODES - face_size)
        repeated |= both_real & (key[:, slot] == key[:, slot + 1])
    degenerate = int(np.count_nonzero(repeated))
    if degenerate:
        raise SystemExit(
            f"[inp_to_polymesh] ERROR: {degenerate} faces have repeated vertices"
        )

    return (
        uniq_sorted,
        uniq_size,
        owner_oriented,
        owner_cell,
        neighbour_cell,
        is_internal,
        cell_volume,
    )


def classify_boundary(
    uniq_sorted: np.ndarray,
    uniq_size: np.ndarray,
    owner_cell: np.ndarray,
    is_internal: np.ndarray,
    membership: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Classify boundary faces; return ``(boundary_face_idx, patch_id)``.

    A face belongs to a patch when *every* one of its real vertices (3 for a
    triangle, 4 for a quad) is in that patch's NSET.  Padding vertices (``-1``)
    are ignored.
    """
    boundary_idx = np.nonzero(~is_internal)[0]
    nodes = uniq_sorted[boundary_idx]  # sorted point indices, padded with -1
    sizes = uniq_size[boundary_idx]
    # np.sort puts the -1 padding at the FRONT, so the real vertices occupy the
    # trailing `size` slots of each row.
    slot_is_real = (np.arange(_MAX_FACE_NODES)[None, :] >= (_MAX_FACE_NODES - sizes)[:, None])
    assigned = np.zeros(boundary_idx.size, dtype=np.int16)
    patch_id = np.full(boundary_idx.size, -1, dtype=np.int16)

    for pid, name in enumerate(_PATCH_ORDER):
        key = f"{name}_NSET"
        mask = membership.get(key)
        if mask is None:
            raise SystemExit(
                f"[inp_to_polymesh] ERROR: boundary NSET {key!r} is missing from input"
            )
        hit = np.ones(boundary_idx.size, dtype=bool)
        for slot in range(_MAX_FACE_NODES):
            real = slot_is_real[:, slot]
            node = np.where(real, nodes[:, slot], 0)
            hit &= (mask[node] | ~real)
        patch_id[hit] = pid
        assigned += hit

    if np.any(assigned != 1):
        unclassified = int(np.count_nonzero(assigned == 0))
        multiple = int(np.count_nonzero(assigned > 1))
        per_patch = [
            (name, int(np.count_nonzero(patch_id == pid)))
            for pid, name in enumerate(_PATCH_ORDER)
        ]
        raise SystemExit(
            "[inp_to_polymesh] ERROR: boundary face classification failed: "
            f"{unclassified} unclassified, {multiple} multiply classified "
            f"(of {boundary_idx.size}); per-patch so far: {per_patch}"
        )

    # A boundary face's owner must exist and be unique.
    if np.any(owner_cell[boundary_idx] < 0):
        raise SystemExit("[inp_to_polymesh] ERROR: boundary face without owner")
    return boundary_idx, patch_id


# ---------------------------------------------------------------------------
# polyMesh writing
# ---------------------------------------------------------------------------


def _foam_header(class_name: str, object_name: str) -> str:
    return (
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        f"    class       {class_name};\n"
        '    location    "constant/polyMesh";\n'
        f"    object      {object_name};\n"
        "}\n\n"
    )


def write_points(path: Path, points: np.ndarray) -> None:
    with open(path, "w", encoding="ascii") as handle:
        handle.write(_foam_header("vectorField", "points"))
        handle.write(f"{points.shape[0]}\n(\n")
        for start in range(0, points.shape[0], _CHUNK):
            chunk = points[start : start + _CHUNK]
            handle.write(
                "\n".join(
                    "(%.12g %.12g %.12g)" % (x, y, z) for x, y, z in chunk.tolist()
                )
            )
            handle.write("\n")
        handle.write(")\n")


def write_faces(path: Path, faces: np.ndarray, face_size: np.ndarray) -> None:
    with open(path, "w", encoding="ascii") as handle:
        handle.write(_foam_header("faceList", "faces"))
        handle.write(f"{faces.shape[0]}\n(\n")
        for start in range(0, faces.shape[0], _CHUNK):
            chunk = faces[start : start + _CHUNK]
            sizes = face_size[start : start + _CHUNK]
            rows = (
                "%d(%s)" % (int(size), " ".join(str(int(v)) for v in row[:size]))
                for row, size in zip(chunk.tolist(), sizes.tolist())
            )
            handle.write("\n".join(rows))
            handle.write("\n")
        handle.write(")\n")


def write_labels(path: Path, class_name: str, object_name: str, labels: np.ndarray) -> None:
    with open(path, "w", encoding="ascii") as handle:
        handle.write(_foam_header(class_name, object_name))
        handle.write(f"{labels.shape[0]}\n(\n")
        for start in range(0, labels.shape[0], _CHUNK):
            chunk = labels[start : start + _CHUNK]
            handle.write("\n".join(str(int(value)) for value in chunk.tolist()))
            handle.write("\n")
        handle.write(")\n")


def write_boundary(path: Path, patches: list[tuple[str, str, int, int]]) -> None:
    with open(path, "w", encoding="ascii") as handle:
        handle.write(_foam_header("polyBoundaryMesh", "boundary"))
        handle.write(f"{len(patches)}\n(\n")
        for name, patch_type, n_faces, start_face in patches:
            handle.write(
                f"    {name}\n"
                "    {\n"
                f"        type        {patch_type};\n"
                f"        nFaces      {n_faces};\n"
                f"        startFace   {start_face};\n"
                "    }\n"
            )
        handle.write(")\n")


def write_control_dict(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     | Version:  v2412                                 |\n"
        "|   \\\\  /    A nd           | Website:  https://openfoam.org                  |\n"
        "|    \\\\/     M anipulation  |                                                 |\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        "    class       dictionary;\n"
        "    object      controlDict;\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
        "\n"
        "application     foamRun;\n"
        "\n"
        "startFrom       startTime;\n"
        "startTime       0;\n"
        "stopAt          endTime;\n"
        "endTime         0;\n"
        "deltaT          1;\n"
        "\n"
        "writeControl    timeStep;\n"
        "writeInterval   1;\n"
        "purgeWrite      0;\n"
        "writeFormat     ascii;\n"
        "writePrecision  6;\n"
        "writeCompression off;\n"
        "timeFormat      general;\n"
        "timePrecision   6;\n"
        "runTimeModifiable true;\n"
        "\n"
        "// ************************************************************************* //\n",
        encoding="ascii",
    )


def _dictionary_header(object_name: str) -> str:
    return (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     | Version:  v2412                                 |\n"
        "|   \\\\  /    A nd           | Website:  https://openfoam.org                  |\n"
        "|    \\\\/     M anipulation  |                                                 |\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        "    class       dictionary;\n"
        f"    object      {object_name};\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
        "\n"
    )


def write_check_mesh_dictionaries(system_dir: Path) -> None:
    """Write minimal fvSchemes/fvSolution.  OpenFOAM v2412 ``checkMesh`` builds
    an ``fvMesh`` and aborts on a missing ``system/fvSchemes`` or
    ``system/fvSolution``, so they are required for the mesh to be checkable.
    """
    system_dir.mkdir(parents=True, exist_ok=True)
    (system_dir / "fvSchemes").write_text(
        _dictionary_header("fvSchemes")
        + "ddtSchemes\n{\n    default         steadyState;\n}\n\n"
        + "gradSchemes\n{\n    default         Gauss linear;\n}\n\n"
        + "divSchemes\n{\n    default         none;\n}\n\n"
        + "laplacianSchemes\n{\n    default         Gauss linear corrected;\n}\n\n"
        + "interpolationSchemes\n{\n    default         linear;\n}\n\n"
        + "snGradSchemes\n{\n    default         corrected;\n}\n\n"
        + "// ************************************************************************* //\n",
        encoding="ascii",
    )
    (system_dir / "fvSolution").write_text(
        _dictionary_header("fvSolution")
        + "solvers\n{\n}\n\n"
        + "SIMPLE\n{\n}\n\n"
        + "PIMPLE\n{\n}\n\n"
        + "relaxationFactors\n{\n}\n\n"
        + "// ************************************************************************* //\n",
        encoding="ascii",
    )


def read_foam_count(path: Path) -> int:
    """Return the leading list-size integer of a written OpenFOAM ascii file."""
    skip_prefixes = (
        "FoamFile",
        "version",
        "format",
        "class",
        "location",
        "object",
        "{",
        "}",
        "(",
        ")",
        "//",
        "/*",
        "*",
        "\\",
        "|",
    )
    with open(path, "r", encoding="ascii") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(skip_prefixes):
                continue
            try:
                return int(stripped)
            except ValueError:
                continue
    raise SystemExit(f"[inp_to_polymesh] ERROR: could not read list size from {path}")


# ---------------------------------------------------------------------------
# Coordinate cross-check
# ---------------------------------------------------------------------------


def cross_check_coordinates(points: np.ndarray, csv_path: Path) -> float:
    """Compare ``points`` against FluidNodeCoor.csv; return max abs difference."""
    if not csv_path.is_file():
        raise SystemExit(
            f"[inp_to_polymesh] ERROR: coordinate CSV not found: {csv_path} "
            "(pass --coord-csv PATH, or --skip-coord-check to disable)"
        )
    with open(csv_path, "r", encoding="ascii") as handle:
        header = handle.readline().split(",")
        if len(header) < 2:
            raise SystemExit(
                f"[inp_to_polymesh] ERROR: malformed coordinate CSV header: {header!r}"
            )
        declared = int(round(float(header[1])))
        if declared != points.shape[0]:
            raise SystemExit(
                "[inp_to_polymesh] ERROR: coordinate CSV declares "
                f"{declared} nodes but mesh has {points.shape[0]}"
            )
        coords = np.empty((declared, 3), dtype=np.float64)
        row = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split(",")
            coords[row, 0] = float(fields[0])
            coords[row, 1] = float(fields[1])
            coords[row, 2] = float(fields[2])
            row += 1
        if row != declared:
            raise SystemExit(
                f"[inp_to_polymesh] ERROR: coordinate CSV has {row} data lines, "
                f"expected {declared}"
            )
    max_abs = float(np.max(np.abs(coords - points)))
    if max_abs > 1e-9:
        raise SystemExit(
            "[inp_to_polymesh] ERROR: coordinate cross-check failed: "
            f"max|diff| = {max_abs:.6e} > 1e-9"
        )
    return max_abs


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _assemble_poly_mesh(
    owner_oriented: np.ndarray,
    owner_cell: np.ndarray,
    neighbour_cell: np.ndarray,
    is_internal: np.ndarray,
    uniq_size: np.ndarray,
    boundary_idx: np.ndarray,
    patch_id: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[tuple[str, str, int, int]]]:
    int_idx = np.nonzero(is_internal)[0]
    int_order = np.lexsort((neighbour_cell[int_idx], owner_cell[int_idx]))
    internal_faces = owner_oriented[int_idx][int_order]
    internal_size = uniq_size[int_idx][int_order]
    internal_owner = owner_cell[int_idx][int_order]
    internal_neighbour = neighbour_cell[int_idx][int_order]

    boundary_owner_of = owner_cell[boundary_idx]
    boundary_faces_of = owner_oriented[boundary_idx]
    boundary_size_of = uniq_size[boundary_idx]

    patch_faces: list[np.ndarray] = []
    patch_sizes: list[np.ndarray] = []
    patch_owners: list[np.ndarray] = []
    patch_meta: list[tuple[str, str, int, int]] = []
    cursor = internal_faces.shape[0]

    for pid, name in enumerate(_PATCH_ORDER):
        sel = np.nonzero(patch_id == pid)[0]
        faces_here = boundary_faces_of[sel]
        owners_here = boundary_owner_of[sel]
        order = np.lexsort(
            (
                faces_here[:, 3],
                faces_here[:, 2],
                faces_here[:, 1],
                faces_here[:, 0],
                owners_here,
            )
        )
        faces_here = faces_here[order]
        owners_here = owners_here[order]
        sizes_here = boundary_size_of[sel][order]
        patch_faces.append(faces_here)
        patch_sizes.append(sizes_here)
        patch_owners.append(owners_here)
        patch_meta.append((name, _PATCH_TYPE[name], int(faces_here.shape[0]), cursor))
        cursor += int(faces_here.shape[0])

    all_faces = np.concatenate([internal_faces, *patch_faces], axis=0)
    all_size = np.concatenate([internal_size, *patch_sizes], axis=0)
    all_owner = np.concatenate([internal_owner, *patch_owners], axis=0)
    return all_faces, all_size, all_owner, internal_neighbour, patch_meta


def _print_summary(
    points: np.ndarray,
    etype: np.ndarray,
    faces: np.ndarray,
    face_size: np.ndarray,
    neighbour: np.ndarray,
    patch_meta: list[tuple[str, str, int, int]],
    max_coord_diff: float | None,
    volume_min: float,
    volume_max: float,
) -> None:
    n_internal = neighbour.shape[0]
    n_tri = int(np.count_nonzero(face_size == 3))
    n_quad = int(np.count_nonzero(face_size == 4))
    by_type = ", ".join(
        f"{name}={int(np.count_nonzero(etype == code))}"
        for code, name in _ELEM_TYPE_NAME.items()
    )
    print(
        f"[inp_to_polymesh] input : {points.shape[0]} nodes, "
        f"{etype.size} cells ({by_type})"
    )
    print("[inp_to_polymesh] output: "
          f"{points.shape[0]} points, {faces.shape[0]} faces "
          f"({n_internal} internal, {faces.shape[0] - n_internal} boundary; "
          f"{n_quad} quads, {n_tri} triangles)")
    for name, patch_type, n_faces, start_face in patch_meta:
        print(
            f"[inp_to_polymesh]   patch {name:<7} type={patch_type:<5} "
            f"nFaces={n_faces:<7} startFace={start_face}"
        )
    print(
        "[inp_to_polymesh] cell volume min/max = "
        f"{volume_min:.6e} / {volume_max:.6e} (all positive)"
    )
    if max_coord_diff is not None:
        print(f"[inp_to_polymesh] coordinate cross-check max|diff| = {max_coord_diff:.3e}")


def _run_self_check(
    points: np.ndarray,
    faces: np.ndarray,
    neighbour: np.ndarray,
    mesh_dir: Path,
) -> bool:
    checks = [
        ("points", points.shape[0]),
        ("faces", faces.shape[0]),
        ("owner", faces.shape[0]),
        ("neighbour", neighbour.shape[0]),
    ]
    ok = True
    print("[inp_to_polymesh] self-check (written vs expected):")
    for name, expected in checks:
        written = read_foam_count(mesh_dir / name)
        status = "OK" if written == expected else "MISMATCH"
        if written != expected:
            ok = False
        print(f"[inp_to_polymesh]   {name:<10} written={written:<9} expected={expected:<9} {status}")
    return ok


def convert(
    input_path: Path,
    out_dir: Path,
    coord_csv: Path | None,
    self_check: bool,
) -> int:
    print(f"[inp_to_polymesh] reading {input_path}")
    points, conn, etype, membership = parse_inp(input_path)

    full = membership.get(_FULL_SET)
    if full is None or not bool(np.all(full)):
        raise SystemExit(
            f"[inp_to_polymesh] ERROR: {_FULL_SET} is missing or does not cover "
            "every node"
        )
    for name in _PATCH_ORDER:
        key = f"{name}_NSET"
        if key not in membership:
            raise SystemExit(f"[inp_to_polymesh] ERROR: required NSET {key!r} missing")

    (
        uniq_sorted,
        uniq_size,
        owner_oriented,
        owner_cell,
        neighbour_cell,
        is_internal,
        cell_volume,
    ) = build_faces(points, conn, etype)
    boundary_idx, patch_id = classify_boundary(
        uniq_sorted, uniq_size, owner_cell, is_internal, membership
    )
    faces, all_size, owner, neighbour, patch_meta = _assemble_poly_mesh(
        owner_oriented,
        owner_cell,
        neighbour_cell,
        is_internal,
        uniq_size,
        boundary_idx,
        patch_id,
    )

    mesh_dir = out_dir / "constant" / "polyMesh"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    write_points(mesh_dir / "points", points)
    write_faces(mesh_dir / "faces", faces, all_size)
    write_labels(mesh_dir / "owner", "labelList", "owner", owner)
    write_labels(mesh_dir / "neighbour", "labelList", "neighbour", neighbour)
    write_boundary(mesh_dir / "boundary", patch_meta)
    write_control_dict(out_dir / "system" / "controlDict")
    write_check_mesh_dictionaries(out_dir / "system")

    max_coord_diff: float | None = None
    if coord_csv is not None:
        max_coord_diff = cross_check_coordinates(points, coord_csv)

    _print_summary(
        points,
        etype,
        faces,
        all_size,
        neighbour,
        patch_meta,
        max_coord_diff,
        float(cell_volume.min()),
        float(cell_volume.max()),
    )
    print(f"[inp_to_polymesh] wrote {mesh_dir}")
    print(f"[inp_to_polymesh] wrote {out_dir / 'system' / 'controlDict'}")

    if self_check:
        if not _run_self_check(points, faces, neighbour, mesh_dir):
            raise SystemExit("[inp_to_polymesh] ERROR: self-check failed")
        print("[inp_to_polymesh] self-check PASSED")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="inp_to_polymesh.py",
        description=(
            "Convert a CalculiX/Abaqus mesh (mixed C3D8 hex + C3D6 wedge) into "
            "an OpenFOAM polyMesh and a minimal system/controlDict."
        ),
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Plain .inp, .inp.gz or single-member .inp .zip deck",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Case directory; writes constant/polyMesh/* and system/controlDict",
    )
    parser.add_argument(
        "--coord-csv",
        type=Path,
        default=Path(_COORD_CSV_DEFAULT),
        help=(
            "FluidNodeCoor.csv for the ascending-node-id coordinate cross-check "
            f"(default: {_COORD_CSV_DEFAULT})"
        ),
    )
    parser.add_argument(
        "--skip-coord-check",
        action="store_true",
        help="Disable the coordinate cross-check (not recommended)",
    )
    parser.add_argument(
        "--self-check",
        action="store_true",
        help="Re-read the written polyMesh files and compare counts with the input",
    )
    args = parser.parse_args(argv)

    coord_csv = None if args.skip_coord_check else args.coord_csv
    return convert(args.input, args.out_dir, coord_csv, args.self_check)


if __name__ == "__main__":
    sys.exit(main())
