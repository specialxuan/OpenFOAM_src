#!/usr/bin/env python3
"""Stage 01: assemble the AGARD merged mesh deck (far-field C3D8 + wing C3D10).

Merges the trimmed far-field Abaqus deck (hexahedral C3D8 ``FAR_PART_ELSET``
plus its eight retained NSETs) with the second-order C3D10 wing deck
(``WingMode2/mesh.inp`` layout).  Output is a mesh-only deck (nodes, element
blocks, NSETs) that stage 02 turns into a complete CalculiX modal deck.

The two input meshes ship next to this script as ``mesh/farfield_c3d8.zip``
and ``mesh/wing_c3d10.zip`` (single-member archives; see
``mesh/PROVENANCE.md``), so the pipeline is reproducible from the repository
alone.  Inputs are opened transparently by extension: plain ASCII, ``.gz`` or
``.zip``; a zip must hold exactly one deck member unless an explicit member
name is passed in code.  The committed far-field archive is the contract: it
must hold exactly 388,465 nodes, 372,760 C3D8 elements in ``FAR_PART_ELSET``
and exactly the eight NSETs ``FAR_PART_NSET``/``FRONT_NSET``/``BACK_NSET``/
``IN_NSET``/``OUT_NSET``/``TOP_NSET``/``FSI_NSET``/``BOTTOM_NSET``.  A foreign
element type, an unexpected element set or NSET name, or a node id outside
``1..388465`` is rejected, never tolerated.

Merge rule
----------
* far-field nodes keep ids ``1..far_nodes``;
* wing nodes are rank-renumbered by ascending original id into the contiguous
  range ``far_nodes+1 .. far_nodes+wing_nodes`` (the source ids are
  non-contiguous; the rank order is what the FDM CSV contract relies on);
* elements are renumbered contiguously in file order:
  C3D8 ``1..far_elems``, C3D10 ``far_elems+1 .. far_elems+wing_elems``;
* coordinate strings are copied verbatim, never reformatted.

The deck header is emitted by this script and names its real producer.  The
structural bytes from the first ``*NODE`` line to EOF are the interchange
contract shared with the signed-off Mode3 artifact (enforced by the tests);
real provenance (source paths, sha256, every measured assertion) lives in the
sidecar metadata JSON.

Module size note: ``allow: SIZE_OK`` -- this is one CLI stage whose streaming
parsers, geometry assertions and writer intentionally stay in one file so the
stage boundary (and its metadata) is auditable in a single place.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import sys
import time
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Iterator, TextIO, TypedDict

import numpy as np
import numpy.typing as npt
from scipy.spatial import cKDTree

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]
NodeText = tuple[str, str, str]

type JsonValue = str | int | float | bool | None | list["JsonValue"] | dict[str, "JsonValue"]
type Report = dict[str, JsonValue]

SCRIPT_DIR: Final = Path(__file__).resolve().parent
DEFAULT_FARFIELD_DECK: Final = SCRIPT_DIR / "mesh" / "farfield_c3d8.zip"
DEFAULT_WING_DECK: Final = SCRIPT_DIR / "mesh" / "wing_c3d10.zip"
DEFAULT_TIE_TOLERANCE: Final = 0.05
AGARD_ROOT: Final = Path("/root/Workspace/AGARD")

# C3D10 mid-edge pairing (nodes 5..10, 0-based) in CalculiX order.
TET_EDGES: Final = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))

HEADING: Final = "AGARD 445.6 merged base: C3D8 far-field + C3D10 wing"


@dataclass(frozen=True, slots=True)
class MergeContract:
    """Expected source counts and set names for one merge invocation."""

    far_nodes: int
    far_elems: int
    wing_nodes: int
    wing_elems: int
    node_offset: int
    elem_offset: int
    wing_volume: float
    wing_node_gaps: int
    fsi_nodes: int
    far_elset: str = "FAR_PART_ELSET"
    wing_elset: str = "WING"
    far_sets: tuple[str, ...] = (
        "FAR_PART_NSET", "FSI_NSET", "FRONT_NSET", "BACK_NSET",
        "IN_NSET", "OUT_NSET", "TOP_NSET", "BOTTOM_NSET")
    outer_sets: tuple[str, ...] = (
        "FRONT", "BACK", "IN", "OUT", "TOP", "BOTTOM")
    wing_sets: tuple[str, ...] = ("ROOT", "TIP", "SIDE", "WING")

    @property
    def total_nodes(self) -> int:
        return self.far_nodes + self.wing_nodes

    @property
    def total_elems(self) -> int:
        return self.far_elems + self.wing_elems


AGARD_CONTRACT: Final = MergeContract(
    far_nodes=388465, far_elems=372760,
    wing_nodes=977086, wing_elems=532373,
    node_offset=388465, elem_offset=372760,
    wing_volume=4.4542836953e-3, wing_node_gaps=12964, fsi_nodes=4045)


class FarfieldData(TypedDict):
    node_text: list[NodeText | None]
    elem_ids: IntArray
    elem_conn: IntArray
    sets: dict[str, IntArray]


class WingData(TypedDict):
    node_ids: IntArray
    node_text: list[NodeText]
    elem_conn: IntArray
    sets: dict[str, IntArray]


class MergeError(RuntimeError):
    """A source deck violates the merge contract."""


@dataclass(frozen=True, slots=True)
class FarfieldCheck:
    coords: FloatArray
    fixed_sets: dict[str, IntArray]
    report: Report


@dataclass(frozen=True, slots=True)
class WingCheck:
    old_to_new: IntArray
    sets_new: dict[str, IntArray]
    report: Report


@dataclass(frozen=True, slots=True)
class GeometryCheck:
    wing_coords: FloatArray
    report: Report


def _header(line: str) -> str:
    """Normalised (spaces removed, upper-case) keyword header."""
    return line.replace(" ", "").strip().upper()


def _display_path(path: Path) -> str:
    """Reference-style relative path under the AGARD root or this directory."""
    resolved = path.resolve()
    for root in (AGARD_ROOT, SCRIPT_DIR):
        try:
            return str(resolved.relative_to(root))
        except ValueError:
            continue
    return resolved.name


def sha256(path: Path) -> str:
    """Streaming sha256 of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _zip_member(path: Path, archive: zipfile.ZipFile,
                member: str | None) -> str:
    """The one deck member of an archive, or the explicitly requested name."""
    names = [name for name in archive.namelist() if not name.endswith("/")]
    if member is not None:
        if member not in names:
            raise MergeError(
                f"{path}: zip archive has no member {member!r} "
                f"(members: {', '.join(names) or 'none'})")
        return member
    if len(names) != 1:
        raise MergeError(
            f"{path}: zip archive must hold exactly one deck member, found "
            f"{len(names)}: {', '.join(names) or 'none'}")
    return names[0]


@contextmanager
def open_deck_text(path: Path, member: str | None = None) -> Iterator[TextIO]:
    """Open a deck by extension: plain ASCII, ``.gz`` or one-member ``.zip``.

    ``member`` names the archive entry explicitly; without it a zip must hold
    exactly one non-directory entry.  Malformed archives, non-ASCII bytes and
    missing members are reported as :class:`MergeError`, never a traceback.
    """
    if not path.is_file():
        raise MergeError(f"source deck not found: {path}")
    suffix = path.suffix.lower()
    try:
        if suffix == ".zip":
            with zipfile.ZipFile(path) as archive:
                text = archive.read(_zip_member(path, archive, member)).decode(
                    "ascii")
        elif suffix == ".gz":
            with gzip.open(path, "rt", encoding="ascii") as source:
                text = source.read()
        else:
            text = path.read_text(encoding="ascii")
    except (OSError, zipfile.BadZipFile, UnicodeDecodeError) as error:
        raise MergeError(f"cannot read source deck {path}: {error}") from error
    with io.StringIO(text) as source:
        yield source


def parse_farfield(path: Path, contract: MergeContract,
                   member: str | None = None) -> FarfieldData:
    """Stream the far-field deck: kept nodes, the C3D8 block, the kept NSETs."""
    node_text: list[NodeText | None] = [None] * (contract.far_nodes + 1)
    elem_ids: list[int] = []
    elem_conn: list[tuple[int, ...]] = []
    sets: dict[str, list[int]] = {name: [] for name in contract.far_sets}
    section: str | None = None
    with open_deck_text(path, member) as source:
        for line in source:
            if line.startswith("*"):
                head = _header(line)
                if head == "*NODE" or head.startswith("*NODE,"):
                    section = "node"
                elif head.startswith("*ELEMENT"):
                    tokens = {token.split("=", 1)[0]: token.split("=", 1)[1]
                              for token in head.split(",") if "=" in token}
                    if (tokens.get("TYPE") != "C3D8"
                            or tokens.get("ELSET") != contract.far_elset):
                        raise MergeError(
                            f"unexpected far-field element card {line.strip()!r} "
                            f"(expected TYPE=C3D8, ELSET={contract.far_elset})")
                    section = "elem"
                elif head.startswith("*NSET"):
                    if "NSET=" not in head:
                        raise MergeError(
                            f"far-field NSET card without a name: {line.strip()!r}")
                    name = head.split("NSET=")[1].split(",")[0]
                    if name not in sets:
                        raise MergeError(
                            f"unexpected far-field NSET {name!r} (expected: "
                            f"{', '.join(contract.far_sets)})")
                    section = f"nset:{name}"
                else:
                    section = None
                continue
            if section == "node":
                fields = line.split(",")
                node_id = int(fields[0])
                if not 1 <= node_id <= contract.far_nodes:
                    raise MergeError(
                        f"far-field node id {node_id} outside "
                        f"1..{contract.far_nodes}")
                node_text[node_id] = (fields[1].strip(), fields[2].strip(),
                                      fields[3].strip())
            elif section == "elem":
                fields = [int(value) for value in line.split(",")[:9]]
                elem_ids.append(fields[0])
                elem_conn.append(tuple(fields[1:9]))
            elif section is not None and section.startswith("nset:"):
                name = section.split(":", 1)[1]
                sets[name].extend(int(value) for value in
                                  line.strip().split(",") if value.strip())
    return {"node_text": node_text,
            "elem_ids": np.asarray(elem_ids, dtype=np.int64),
            "elem_conn": np.asarray(elem_conn, dtype=np.int64),
            "sets": {key: np.asarray(values, dtype=np.int64)
                     for key, values in sets.items()}}


def parse_wing(path: Path, contract: MergeContract,
               member: str | None = None) -> WingData:
    """Stream the C3D10 wing deck: all nodes, the C3D10 block, its NSETs."""
    node_ids: list[int] = []
    node_text: list[NodeText] = []
    elem_conn: list[tuple[int, ...]] = []
    sets: dict[str, list[int]] = {name: [] for name in contract.wing_sets}
    section: str | None = None
    with open_deck_text(path, member) as source:
        for line in source:
            if line.startswith("*"):
                head = _header(line)
                if head == "*NODE" or head.startswith("*NODE,"):
                    section = "node"
                elif head.startswith("*ELEMENT"):
                    section = "elem" if "TYPE=C3D10" in head else None
                elif head.startswith("*NSET"):
                    name = head.split("NSET=")[1].split(",")[0]
                    section = f"nset:{name}" if name in sets else None
                else:
                    section = None
                continue
            if section == "node":
                fields = line.split(",")
                node_ids.append(int(fields[0]))
                node_text.append((fields[1].strip(), fields[2].strip(),
                                  fields[3].strip()))
            elif section == "elem":
                fields = [int(value) for value in line.split(",")[:11]]
                elem_conn.append(tuple(fields[1:11]))
            elif section is not None and section.startswith("nset:"):
                name = section.split(":", 1)[1]
                sets[name].extend(int(value) for value in
                                  line.strip().split(",") if value.strip())
    return {"node_ids": np.asarray(node_ids, dtype=np.int64),
            "node_text": node_text,
            "elem_conn": np.asarray(elem_conn, dtype=np.int64),
            "sets": {key: np.asarray(values, dtype=np.int64)
                     for key, values in sets.items()}}


def quadratic_tet_volume(corners: FloatArray) -> FloatArray:
    """Volume of a C3D10 from its four corner nodes (mids at edge midpoints)."""
    a = corners[:, 1] - corners[:, 0]
    b = corners[:, 2] - corners[:, 0]
    c = corners[:, 3] - corners[:, 0]
    return np.abs(np.einsum("ij,ij->i", np.cross(a, b), c)) / 6.0


def write_nset(output: TextIO, name: str, values: IntArray) -> None:
    """One ``*NSET`` with 16 ids per continuation line."""
    output.write(f"*NSET, NSET={name}\n")
    text = [str(int(value)) for value in values]
    for offset in range(0, len(text), 16):
        output.write(",".join(text[offset:offset + 16]) + "\n")


def deck_header(farfield_deck: Path, wing_deck: Path,
                contract: MergeContract) -> str:
    """The merged-deck header, templated with the live contract and producer."""
    return (
        "** ---------------------------------------------------------------------------\n"
        "** AGARD 445.6 merged base deck (far-field C3D8 + C3D10 wing)\n"
        f"** far-field : {_display_path(farfield_deck)}  "
        f"nodes 1..{contract.far_nodes}\n"
        f"**             C3D8 ELSET={contract.far_elset}  "
        f"ids renumbered 1..{contract.far_elems}\n"
        f"** wing      : {_display_path(wing_deck)}  C3D10  "
        f"ids renumbered {contract.elem_offset + 1}..{contract.total_elems}\n"
        f"**             nodes rank-renumbered to "
        f"{contract.node_offset + 1}..{contract.total_nodes}\n"
        "** generated by case/AGARD/Mode/01_assemble_mesh.py -- "
        "mesh only, placeholder material\n"
        "** ---------------------------------------------------------------------------\n")


def check_farfield(far: FarfieldData, contract: MergeContract) -> FarfieldCheck:
    """Validate the far-field block; return coordinates, FIXED sets, report."""
    far_text = far["node_text"]
    if not all(far_text[node] is not None for node in range(1, contract.far_nodes + 1)):
        raise MergeError("far-field node block is missing ids inside "
                         f"1..{contract.far_nodes}")
    elem_ids = far["elem_ids"]
    if elem_ids.shape[0] != contract.far_elems:
        raise MergeError(f"far C3D8 count {elem_ids.shape[0]} != {contract.far_elems}")
    if not np.all(np.diff(elem_ids) > 0):
        raise MergeError("far element ids are not ascending")
    conn = far["elem_conn"]
    if conn.min() < 1 or conn.max() > contract.far_nodes:
        raise MergeError(f"far connectivity outside 1..{contract.far_nodes}: "
                         f"[{int(conn.min())},{int(conn.max())}]")
    if not np.array_equal(np.unique(conn), np.arange(1, contract.far_nodes + 1)):
        raise MergeError("far connectivity does not use every node 1..far_nodes")
    far_part = far["sets"]["FAR_PART_NSET"]
    if far_part.shape[0] != contract.far_nodes or not np.array_equal(
            np.sort(far_part), np.arange(1, contract.far_nodes + 1)):
        raise MergeError("FAR_PART_NSET is not exactly 1..far_nodes")
    fsi_ids = far["sets"]["FSI_NSET"]
    if fsi_ids.shape[0] != contract.fsi_nodes:
        raise MergeError(f"FSI_NSET size {fsi_ids.shape[0]} != {contract.fsi_nodes}")
    coords = np.zeros((contract.far_nodes + 1, 3), dtype=np.float64)
    for node_id in range(1, contract.far_nodes + 1):
        text = far_text[node_id]
        coords[node_id] = (float(text[0]), float(text[1]), float(text[2]))
    fsi_set = set(fsi_ids.tolist())
    fixed: dict[str, IntArray] = {}
    for name in contract.outer_sets:
        values = far["sets"][f"{name}_NSET"]
        keep = np.array(sorted(set(values.tolist()) - fsi_set), dtype=np.int64)
        if keep.shape[0] == 0:
            raise MergeError(f"FIXED_{name} is empty after removing FSI nodes")
        fixed[f"FIXED_{name}"] = keep
    report: Report = {
        "nodes": contract.far_nodes,
        "c3d8_elements": int(elem_ids.shape[0]),
        "fsi_nodes": int(fsi_ids.shape[0]),
        "outer_set_sizes": {name: int(far["sets"][f"{name}_NSET"].shape[0])
                            for name in contract.outer_sets},
        "fixed_set_sizes": {name: int(values.shape[0])
                            for name, values in fixed.items()},
        "element_id_min": int(elem_ids.min()),
        "element_id_max": int(elem_ids.max()),
    }
    return FarfieldCheck(coords, fixed, report)


def check_wing(wing: WingData, contract: MergeContract) -> WingCheck:
    """Validate the wing block; return the id map, renumbered sets, report."""
    old_ids = wing["node_ids"]
    if old_ids.shape[0] != contract.wing_nodes:
        raise MergeError(f"wing node count {old_ids.shape[0]} != {contract.wing_nodes}")
    if not np.all(np.diff(old_ids) > 0):
        raise MergeError("wing node ids are not ascending")
    gaps = int(old_ids.max() - old_ids.min() + 1 - old_ids.shape[0])
    if gaps != contract.wing_node_gaps:
        raise MergeError(f"wing node id gaps {gaps} != {contract.wing_node_gaps}")
    old_to_new = np.zeros(int(old_ids.max()) + 1, dtype=np.int64)
    old_to_new[old_ids] = np.arange(contract.node_offset + 1,
                                    contract.node_offset + 1 + contract.wing_nodes,
                                    dtype=np.int64)
    expected = np.arange(contract.node_offset + 1, contract.total_nodes + 1,
                         dtype=np.int64)
    if not np.array_equal(old_to_new[old_ids], expected):
        raise MergeError("wing node renumbering is not the contiguous rank mapping")
    conn_old = wing["elem_conn"]
    if conn_old.shape[0] != contract.wing_elems:
        raise MergeError(f"wing C3D10 count {conn_old.shape[0]} != {contract.wing_elems}")
    if conn_old.min() < 1 or conn_old.max() > int(old_ids.max()):
        raise MergeError("wing connectivity references an unknown node id")
    if old_to_new[conn_old].min() <= contract.far_nodes:
        raise MergeError("wing connectivity reuses far-field node ids (must be disjoint)")
    sizes = {name: int(wing["sets"][name].shape[0]) for name in contract.wing_sets}
    if sizes.get("WING") != contract.wing_nodes:
        raise MergeError(f"WING NSET size {sizes.get('WING')} != {contract.wing_nodes}")
    sets_new = {name: old_to_new[wing["sets"][name]]
                for name in contract.wing_sets}
    report: Report = {
        "nodes": contract.wing_nodes,
        "c3d10_elements": int(conn_old.shape[0]),
        "node_id_min": int(old_ids.min()),
        "node_id_max": int(old_ids.max()),
        "node_id_gaps": gaps,
        "set_sizes": sizes,
    }
    return WingCheck(old_to_new, sets_new, report)


def check_geometry(wing: WingData, contract: MergeContract) -> GeometryCheck:
    """Wing volume, C3D10 mid-edge positions, and the wing coordinate table."""
    wing_coords = np.array([(float(x), float(y), float(z))
                            for x, y, z in wing["node_text"]], dtype=np.float64)
    row_of_id = np.zeros(int(wing["node_ids"].max()) + 1, dtype=np.int64)
    row_of_id[wing["node_ids"]] = np.arange(contract.wing_nodes, dtype=np.int64)
    conn_old = wing["elem_conn"]
    corners = wing_coords[row_of_id[conn_old[:, :4]]]
    volume = float(quadratic_tet_volume(corners).sum())
    rel_volume = abs(volume - contract.wing_volume) / contract.wing_volume
    if rel_volume >= 1e-6:
        raise MergeError(f"wing volume {volume:.12e} != {contract.wing_volume:.12e} "
                         f"(rel {rel_volume:.2e})")
    # The source stores 10-significant-digit ASCII, so the achievable tolerance
    # is relative to the COORDINATE scale, not the edge length: the worst
    # edge-relative deviation sits on a sub-0.1 mm edge and is print rounding.
    mid_abs = mid_rel_coord = mid_rel_edge = 0.0
    for mid, (i, j) in zip(range(4, 10), TET_EDGES):
        midpoint = 0.5 * (corners[:, i] + corners[:, j])
        got = wing_coords[row_of_id[conn_old[:, mid]]]
        deviation = np.linalg.norm(got - midpoint, axis=1)
        length = np.linalg.norm(corners[:, j] - corners[:, i], axis=1)
        scale = np.maximum(1.0, np.linalg.norm(midpoint, axis=1))
        mid_abs = max(mid_abs, float(deviation.max()))
        mid_rel_coord = max(mid_rel_coord, float((deviation / scale).max()))
        mid_rel_edge = max(mid_rel_edge,
                           float((deviation / np.maximum(length, 1e-300)).max()))
    if mid_rel_coord > 1e-9:
        raise MergeError(f"C3D10 mid-edge deviation {mid_rel_coord:.3e} "
                         "(coordinate-relative) > 1e-9")
    report: Report = {
        "volume_m3": volume,
        "volume_rel_err": rel_volume,
        "mid_edge_abs_max_m": mid_abs,
        "mid_edge_rel_coord_max": mid_rel_coord,
        "mid_edge_rel_edge_max": mid_rel_edge,
    }
    return GeometryCheck(wing_coords, report)


def fsi_side_distances(coords: FloatArray, far: FarfieldData,
                       wing: WingData, wing_coords: FloatArray) -> FloatArray:
    """Nearest-neighbour distance of every FSI node to the wing SIDE set."""
    side_rows = np.searchsorted(wing["node_ids"], wing["sets"]["SIDE"])
    tree = cKDTree(wing_coords[side_rows])
    distance, _ = tree.query(coords[far["sets"]["FSI_NSET"]], k=1)
    return np.asarray(distance, dtype=np.float64)


def assemble(farfield_deck: Path, wing_deck: Path, output: Path, meta: Path,
             contract: MergeContract = AGARD_CONTRACT,
             tie_tolerance: float = DEFAULT_TIE_TOLERANCE,
             farfield_member: str | None = None,
             wing_member: str | None = None) -> Report:
    """Merge the two source decks, write deck + metadata; return the report."""
    started = time.time()
    for path in (farfield_deck, wing_deck):
        if not path.is_file():
            raise MergeError(f"source deck not found: {path}")

    far = parse_farfield(farfield_deck, contract, farfield_member)
    far_check = check_farfield(far, contract)
    wing = parse_wing(wing_deck, contract, wing_member)
    wing_check = check_wing(wing, contract)
    geometry = check_geometry(wing, contract)

    distance = fsi_side_distances(far_check.coords, far, wing, geometry.wing_coords)
    distance_report: Report = {
        "max_m": float(distance.max()),
        "mean_m": float(distance.mean()),
        "count": int(distance.shape[0]),
        "all_within_tolerance": bool((distance <= tie_tolerance).all()),
        "tie_tolerance": tie_tolerance,
    }
    if (distance > tie_tolerance).any():
        raise MergeError(f"{int((distance > tie_tolerance).sum())} FSI nodes "
                         f"farther than {tie_tolerance} m from the wing SIDE set")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="ascii") as dst:
        dst.write(deck_header(farfield_deck, wing_deck, contract))
        dst.write(f"*HEADING\n{HEADING}\n*NODE\n")
        far_text = far["node_text"]
        for node_id in range(1, contract.far_nodes + 1):
            text = far_text[node_id]
            dst.write(f"{node_id},{text[0]},{text[1]},{text[2]}\n")
        for rank, (x, y, z) in enumerate(wing["node_text"]):
            dst.write(f"{contract.node_offset + 1 + rank},{x},{y},{z}\n")
        dst.write(f"*ELEMENT, TYPE=C3D8, ELSET={contract.far_elset}\n")
        for row, element in enumerate(far["elem_conn"]):
            dst.write(f"{row + 1}," + ",".join(str(int(v)) for v in element) + "\n")
        dst.write(f"*ELEMENT, TYPE=C3D10, ELSET={contract.wing_elset}\n")
        for row, element in enumerate(wing["elem_conn"]):
            dst.write(f"{contract.elem_offset + row + 1}," +
                      ",".join(str(int(v)) for v in
                               wing_check.old_to_new[element]) + "\n")
        for name in contract.far_sets:
            write_nset(dst, name, far["sets"][name])
        for name in contract.outer_sets:
            write_nset(dst, f"FIXED_{name}", far_check.fixed_sets[f"FIXED_{name}"])
        for name in contract.wing_sets:
            write_nset(dst, name, wing_check.sets_new[name])

    report: Report = {
        "stage": "01_assemble_mesh",
        "sources": {"farfield": str(farfield_deck.resolve()),
                    "wing": str(wing_deck.resolve())},
        "sha256": {"farfield": sha256(farfield_deck), "wing": sha256(wing_deck)},
        "counts": {"nodes": contract.total_nodes,
                   "c3d8": contract.far_elems,
                   "c3d10": contract.wing_elems,
                   "elements": contract.total_elems},
        "node_numbering": {"far": f"1..{contract.far_nodes}",
                           "wing": f"{contract.node_offset + 1}..{contract.total_nodes}"},
        "element_numbering": {"c3d8": f"1..{contract.far_elems}",
                              "c3d10": f"{contract.elem_offset + 1}..{contract.total_elems}"},
        "far": far_check.report,
        "wing": {**wing_check.report, **geometry.report},
        "fsi_wing_distance": distance_report,
        "output": {"path": str(output.resolve()), "bytes": output.stat().st_size,
                   "sha256": sha256(output)},
        "elapsed_s": round(time.time() - started, 1),
    }
    if output.stat().st_size <= 100_000 and contract.total_nodes > 10_000:
        raise MergeError("merged deck looks truncated")
    meta.parent.mkdir(parents=True, exist_ok=True)
    meta.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--farfield-deck", type=Path, default=DEFAULT_FARFIELD_DECK,
                        help="far-field Abaqus deck with C3D8 elements")
    parser.add_argument("--wing-deck", type=Path, default=DEFAULT_WING_DECK,
                        help="second-order C3D10 wing deck")
    parser.add_argument("--output", type=Path, required=True,
                        help="merged mesh-only deck to write")
    parser.add_argument("--meta", type=Path, default=None,
                        help="metadata JSON (default: <output stem>_meta.json)")
    parser.add_argument("--tie-tolerance", type=float, default=DEFAULT_TIE_TOLERANCE,
                        help="FSI-to-SIDE distance gate in metres (default 0.05)")
    args = parser.parse_args()
    if args.tie_tolerance <= 0:
        parser.error("--tie-tolerance must be positive")
    meta = args.meta if args.meta is not None else args.output.with_name(
        args.output.stem + "_meta.json")
    try:
        report = assemble(args.farfield_deck, args.wing_deck, args.output, meta,
                          tie_tolerance=args.tie_tolerance)
    except MergeError as error:
        print(f"[AGARD ERR] {error}", file=sys.stderr)
        return 1
    counts = report["counts"]
    distance = report["fsi_wing_distance"]
    print(f"[AGARD 01] merged={args.output} nodes={counts['nodes']} "
          f"c3d8={counts['c3d8']} c3d10={counts['c3d10']} "
          f"fsi_max_dist={distance['max_m']:.6g}m meta={meta}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
