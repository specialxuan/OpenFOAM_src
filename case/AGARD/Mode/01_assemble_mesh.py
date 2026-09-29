#!/usr/bin/env python3
"""Stage 01: assemble the AGARD merged mesh deck (far-field C3D8 + wing C3D10).

Merges the trimmed far-field Abaqus deck (hexahedral C3D8 ``FAR_PART_ELSET``
plus its eight retained NSETs) with the second-order C3D10 wing deck
(``case/AGARD/Mode/mesh/wing_c3d10.zip`` member ``wing_c3d10.inp``).  Output is
a mesh-only deck (nodes, element blocks, NSETs) that stage 02 turns into a
complete CalculiX modal deck.

The two input meshes ship next to this script as ``mesh/farfield_c3d8.zip``
and ``mesh/wing_c3d10.zip`` (single-member archives; see
``mesh/PROVENANCE.md``), so the pipeline is reproducible from the repository
alone.  Inputs are opened transparently by extension: plain ASCII, ``.gz`` or
``.zip``; a zip must hold exactly one deck member unless an explicit member
name is passed in code.

The far field is selected by an explicit named *profile*
(``--farfield-profile`` ``legacy`` or ``fluid``, default ``legacy``).  The
legacy profile is the signed-off baseline: it pins exactly 388,465 nodes,
372,760 C3D8 elements in ``FAR_PART_ELSET`` and exactly the eight NSETs
``FAR_PART_NSET``/``FRONT_NSET``/``BACK_NSET``/``IN_NSET``/``OUT_NSET``/
``TOP_NSET``/``FSI_NSET``/``BOTTOM_NSET``.  The fluid profile is the actively
iterated ICEM family (Fluid7, Fluid8, ...): it requires an explicit
``--farfield-deck`` and *derives* the far-field node count, id gaps, per-type
element counts and ``FSI_NSET`` size from that deck instead of pinning one
member's numbers, recording every derived count and ignored NSET alias in the
metadata.  A derived deck still has to satisfy the family structure: ascending
positive node ids rank-renumbered into ``1..far_nodes``, the eight contract
NSETs present and non-empty, ``FAR_PART_NSET`` exactly the full node id set,
element connectivity inside and exactly covering that set, only ``C3D8``/
``C3D6`` blocks on ``FAR_PART_ELSET`` with ascending ids per block, and the
wing side pinned exactly as for ``legacy``.  Unknown NSET names are ignored
and recorded, so a new ICEM alias never needs a code change; a foreign element
type, an unexpected element set, or a structural violation is rejected, never
tolerated.

Merge rule
----------
* far-field nodes are rank-renumbered by ascending original id into the dense
  range ``1..far_nodes`` (the fluid source ids are gapped; for the gap-free
  legacy deck this rank map is the identity), and every connectivity and NSET
  reference follows the same map;
* wing nodes are rank-renumbered by ascending original id into the contiguous
  range ``far_nodes+1 .. far_nodes+wing_nodes`` (the source ids are
  non-contiguous; the rank order is what the FDM CSV contract relies on);
* elements are renumbered contiguously in writer order per type:
  C3D8 ``1..c3d8``, optional C3D6 ``c3d8+1..far_elems``,
  C3D10 ``far_elems+1 .. far_elems+wing_elems``;
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
from dataclasses import dataclass, replace
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
    """Expected source counts and set names for one merge invocation.

    ``derive_far=False`` pins the far-field counts (legacy): every pinned
    value is checked and must match.  ``derive_far=True`` (fluid) treats the
    far-field count fields as unset placeholders and resolves them from the
    deck; only the wing side, the set names and the declared element types
    stay pinned, so any further ICEM export of the same family merges with no
    code change.
    """

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
    # Named profile this contract belongs to; recorded in the metadata.
    profile: str = "legacy"
    # Far-field element types accepted, in writer order: ``("C3D8",)`` is the
    # legacy single-type contract, the fluid family adds the ``C3D6`` wedge
    # sheet as an optional second type.
    far_elem_types: tuple[str, ...] = ("C3D8",)
    # Explicit per-type element counts; empty means the single-type legacy
    # default ``{"C3D8": far_elems}``.  A derived contract replaces this with
    # the counts measured in the deck before anything is written.
    far_type_counts: tuple[tuple[str, int], ...] = ()
    # Far-field node ids are rank-renumbered to a dense ``1..far_nodes``.  The
    # accepted source envelope is ``1..(far_nodes + far_node_id_gaps)``: the
    # legacy deck has no gaps, a fluid deck scatters its gaps across the range.
    far_node_id_gaps: int = 0
    # NSET names skipped without error (comparison is case-insensitive because
    # headers are normalised to upper-case).  A pinned profile rejects any
    # other unknown set name; a derived profile ignores and records it.
    ignored_sets: tuple[str, ...] = ()
    # Derive the far-field counts/gaps/type split from the deck instead of
    # checking them against the fields above (the fluid family).
    derive_far: bool = False

    @property
    def total_nodes(self) -> int:
        return self.far_nodes + self.wing_nodes

    @property
    def total_elems(self) -> int:
        return self.far_elems + self.wing_elems

    @property
    def far_node_id_max(self) -> int | None:
        """Declared source-id envelope; ``None`` for a derived profile."""
        if self.derive_far:
            return None
        return self.far_nodes + self.far_node_id_gaps

    @property
    def expected_type_counts(self) -> dict[str, int]:
        if self.far_type_counts:
            return {name: count for name, count in self.far_type_counts}
        return {"C3D8": self.far_elems}


# Abaqus element type -> node count (element data lines are copied by width).
ELEM_NODE_COUNT: Final = {"C3D8": 8, "C3D6": 6}


LEGACY_CONTRACT: Final = MergeContract(
    far_nodes=388465, far_elems=372760,
    wing_nodes=977086, wing_elems=532373,
    node_offset=388465, elem_offset=372760,
    wing_volume=4.4542836953e-3, wing_node_gaps=12964, fsi_nodes=4045,
    profile="legacy")


# The far-field zero/unset placeholders are resolved from the deck by
# ``resolve_contract`` before any check, header or write; the wing side and
# the set/type names stay pinned so every fluid-family member is identified
# by structure, not by one member's counts.
FLUID_CONTRACT: Final = MergeContract(
    far_nodes=0, far_elems=0,
    wing_nodes=977086, wing_elems=532373,
    node_offset=0, elem_offset=0,
    wing_volume=4.4542836953e-3, wing_node_gaps=12964, fsi_nodes=0,
    profile="fluid",
    far_elem_types=("C3D8", "C3D6"),
    derive_far=True)


# Backwards-compatible alias: the default profile stays the legacy contract.
AGARD_CONTRACT: Final = LEGACY_CONTRACT
PROFILES: Final = {"legacy": LEGACY_CONTRACT, "fluid": FLUID_CONTRACT}


class FarfieldBlock(TypedDict):
    type: str
    ids: IntArray
    conn: IntArray


class FarfieldData(TypedDict):
    node_text: list[NodeText | None]
    blocks: list[FarfieldBlock]
    sets: dict[str, IntArray]
    node_id_gaps: int
    node_id_max: int
    ignored_sets: list[str]


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


def _nset_spelling(line: str) -> str | None:
    """Source-spelled name of the first ``NSET=`` token, or ``None``."""
    compact = line.replace(" ", "")
    index = compact.upper().find("NSET=")
    if index == -1:
        return None
    return compact[index + 5:].split(",")[0].strip() or None


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
    """Stream the far-field deck: nodes, element blocks and the kept NSETs.

    Node ids may be gapped (the fluid family): they must be positive and
    strictly ascending and are rank-renumbered to a dense ``1..far_nodes``,
    with every connectivity and NSET reference following the same map.  A
    pinned profile checks the source-id envelope and rejects any unknown NSET;
    a derived profile has no count envelope and ignores+records every unknown
    NSET alias so a new ICEM export needs no code change.
    """
    source_ids: list[int] = []
    source_text: list[NodeText] = []
    block_types: list[str] = []
    block_ids: list[list[int]] = []
    block_conn: list[list[tuple[int, ...]]] = []
    sets_src: dict[str, list[int]] = {name: [] for name in contract.far_sets}
    ignored = {name.upper() for name in contract.ignored_sets}
    ignored_names: set[str] = set()
    id_max = contract.far_node_id_max
    section: str | None = None
    active_type: str | None = None
    with open_deck_text(path, member) as source:
        for line in source:
            if not line.strip():
                continue
            if line.startswith("*"):
                head = _header(line)
                if head == "*NODE" or head.startswith("*NODE,"):
                    section = "node"
                    spelling = _nset_spelling(line)
                    if spelling is not None and spelling.upper() not in sets_src:
                        ignored_names.add(spelling)
                elif head.startswith("*ELEMENT"):
                    tokens = {token.split("=", 1)[0]: token.split("=", 1)[1]
                              for token in head.split(",") if "=" in token}
                    etype = tokens.get("TYPE", "")
                    if (etype not in contract.far_elem_types
                            or tokens.get("ELSET") != contract.far_elset):
                        raise MergeError(
                            f"unexpected far-field element card {line.strip()!r} "
                            f"(expected TYPE in "
                            f"{'/'.join(contract.far_elem_types)}, "
                            f"ELSET={contract.far_elset})")
                    section = "elem"
                    active_type = etype
                    block_types.append(etype)
                    block_ids.append([])
                    block_conn.append([])
                elif head.startswith("*NSET"):
                    if "NSET=" not in head:
                        raise MergeError(
                            f"far-field NSET card without a name: {line.strip()!r}")
                    name = head.split("NSET=")[1].split(",")[0]
                    if name in sets_src:
                        section = f"nset:{name}"
                    elif name in ignored or contract.derive_far:
                        ignored_names.add(_nset_spelling(line) or name)
                        section = "ignore"
                    else:
                        raise MergeError(
                            f"unexpected far-field NSET {name!r} (expected: "
                            f"{', '.join(contract.far_sets)})")
                else:
                    section = None
                continue
            if section == "node":
                fields = line.split(",")
                node_id = int(fields[0])
                if node_id < 1 or (id_max is not None and node_id > id_max):
                    limit = f"1..{id_max}" if id_max is not None else "the positive range"
                    raise MergeError(
                        f"far-field node id {node_id} outside {limit}")
                source_ids.append(node_id)
                source_text.append((fields[1].strip(), fields[2].strip(),
                                    fields[3].strip()))
            elif section == "elem":
                width = ELEM_NODE_COUNT[active_type or ""]
                fields = [int(value) for value in line.split(",")[:width + 1]]
                block_ids[-1].append(fields[0])
                block_conn[-1].append(tuple(fields[1:width + 1]))
            elif section == "ignore":
                continue
            elif section is not None and section.startswith("nset:"):
                name = section.split(":", 1)[1]
                sets_src[name].extend(int(value) for value in
                                      line.strip().split(",") if value.strip())
    ids = np.asarray(source_ids, dtype=np.int64)
    if ids.shape[0] == 0:
        raise MergeError("far-field deck declares no nodes")
    if not np.all(np.diff(ids) > 0):
        raise MergeError("far-field node ids are not strictly ascending")
    rank_of_source = np.zeros(int(ids.max()) + 1, dtype=np.int64)
    rank_of_source[ids] = np.arange(1, ids.shape[0] + 1, dtype=np.int64)
    node_text: list[NodeText | None] = [None] * (ids.shape[0] + 1)
    for rank, text in enumerate(source_text, start=1):
        node_text[rank] = text

    def remap(values: object) -> IntArray:
        array = np.asarray(values, dtype=np.int64)
        if array.size and (array.min() < 1
                           or array.max() >= rank_of_source.shape[0]):
            raise MergeError(
                f"far-field reference outside the declared node id range "
                f"1..{int(rank_of_source.shape[0]) - 1}")
        ranks = rank_of_source[array]
        if array.size and np.any(ranks == 0):
            missing = int(array[ranks == 0][0])
            raise MergeError(
                f"far-field reference {missing} names a gapped/absent node "
                f"id (declared range 1..{int(rank_of_source.shape[0]) - 1}, "
                f"no source node {missing})")
        return ranks

    blocks: list[FarfieldBlock] = [
        {"type": etype, "ids": np.asarray(eids, dtype=np.int64),
         "conn": remap(econn)}
        for etype, eids, econn in zip(block_types, block_ids, block_conn)]
    return {"node_text": node_text,
            "blocks": blocks,
            "sets": {key: remap(values) for key, values in sets_src.items()},
            "node_id_gaps": int(ids.max() - ids.min() + 1 - ids.shape[0]),
            "node_id_max": int(ids.max()),
            "ignored_sets": sorted(ignored_names)}


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
    """The merged-deck header, templated with the live contract and producer.

    The legacy profile's bytes are frozen; only a non-legacy profile gains a
    profile line and a per-type element description.
    """
    lines = [
        "** ---------------------------------------------------------------------------\n",
        "** AGARD 445.6 merged base deck (far-field C3D8 + C3D10 wing)\n",
    ]
    if contract.profile != "legacy":
        lines.append(f"** far-field profile: {contract.profile}\n")
    lines.append(f"** far-field : {_display_path(farfield_deck)}  "
                 f"nodes 1..{contract.far_nodes}\n")
    if contract.far_elem_types == ("C3D8",):
        lines.append(f"**             C3D8 ELSET={contract.far_elset}  "
                     f"ids renumbered 1..{contract.far_elems}\n")
    else:
        counts = contract.expected_type_counts
        types = "+".join(contract.far_elem_types)
        detail = " ".join(f"{name}={counts.get(name, 0)}"
                          for name in contract.far_elem_types)
        lines.append(f"**             {types} ELSET={contract.far_elset}  "
                     f"ids renumbered 1..{contract.far_elems}  ({detail})\n")
    lines.extend([
        f"** wing      : {_display_path(wing_deck)}  C3D10  "
        f"ids renumbered {contract.elem_offset + 1}..{contract.total_elems}\n",
        f"**             nodes rank-renumbered to "
        f"{contract.node_offset + 1}..{contract.total_nodes}\n",
        "** generated by case/AGARD/Mode/01_assemble_mesh.py -- "
        "mesh only, placeholder material\n",
        "** ---------------------------------------------------------------------------\n",
    ])
    return "".join(lines)


def resolve_contract(far: FarfieldData, contract: MergeContract) -> MergeContract:
    """Check a pinned contract or derive the far-field counts from the deck.

    A pinned (legacy) contract keeps every equality check it has today and is
    returned unchanged.  A derived (fluid) contract replaces its placeholder
    far-field counts with the deck's measured node count, id gaps, per-type
    element split and ``FSI_NSET`` size, after enforcing that every contract
    NSET is present and non-empty.
    """
    node_count = len(far["node_text"]) - 1
    type_counts: dict[str, int] = {}
    for block in far["blocks"]:
        type_counts[block["type"]] = (type_counts.get(block["type"], 0)
                                      + int(block["ids"].shape[0]))
    total_elems = sum(type_counts.values())
    fsi_count = int(far["sets"]["FSI_NSET"].shape[0])
    got = {name: type_counts.get(name, 0) for name in contract.far_elem_types}
    if not contract.derive_far:
        if node_count != contract.far_nodes:
            raise MergeError(f"far node count {node_count} != {contract.far_nodes}")
        if far["node_id_gaps"] != contract.far_node_id_gaps:
            raise MergeError(f"far node id gaps {far['node_id_gaps']} != "
                             f"{contract.far_node_id_gaps}")
        if total_elems != contract.far_elems:
            raise MergeError(f"far element count {total_elems} != {contract.far_elems}")
        expected = contract.expected_type_counts
        expected_selected = {name: expected.get(name, 0)
                             for name in contract.far_elem_types}
        if got != expected_selected:
            raise MergeError(f"far element type counts {got} != {expected_selected}")
        if fsi_count != contract.fsi_nodes:
            raise MergeError(f"FSI_NSET size {fsi_count} != {contract.fsi_nodes}")
        return contract
    empty = [name for name in contract.far_sets
             if far["sets"][name].shape[0] == 0]
    if empty:
        raise MergeError("far-field contract NSET(s) missing or empty: "
                         f"{', '.join(empty)}")
    return replace(contract,
                   far_nodes=node_count, far_elems=total_elems,
                   node_offset=node_count, elem_offset=total_elems,
                   fsi_nodes=fsi_count,
                   far_node_id_gaps=far["node_id_gaps"],
                   far_type_counts=tuple((name, got[name])
                                         for name in contract.far_elem_types))


def check_farfield(far: FarfieldData, contract: MergeContract) -> FarfieldCheck:
    """Validate the far-field block; return coordinates, FIXED sets, report."""
    far_text = far["node_text"]
    type_counts: dict[str, int] = {}
    used: list[IntArray] = []
    elem_min: int | None = None
    elem_max: int | None = None
    for block in far["blocks"]:
        ids = block["ids"]
        conn = block["conn"]
        etype = block["type"]
        count = int(ids.shape[0])
        type_counts[etype] = type_counts.get(etype, 0) + count
        if count:
            if not np.all(np.diff(ids) > 0):
                raise MergeError(f"far {etype} element ids are not ascending")
            low, high = int(ids.min()), int(ids.max())
            elem_min = low if elem_min is None else min(elem_min, low)
            elem_max = high if elem_max is None else max(elem_max, high)
        if conn.size:
            if int(conn.min()) < 1 or int(conn.max()) > contract.far_nodes:
                raise MergeError(f"far connectivity outside 1..{contract.far_nodes}: "
                                 f"[{int(conn.min())},{int(conn.max())}]")
            used.append(np.unique(conn))
    got = {name: type_counts.get(name, 0) for name in contract.far_elem_types}
    union = (np.unique(np.concatenate(used)) if used
             else np.empty(0, dtype=np.int64))
    if not np.array_equal(union, np.arange(1, contract.far_nodes + 1)):
        raise MergeError("far connectivity does not use every node 1..far_nodes")
    far_part = far["sets"]["FAR_PART_NSET"]
    if far_part.shape[0] != contract.far_nodes or not np.array_equal(
            np.sort(far_part), np.arange(1, contract.far_nodes + 1)):
        raise MergeError("FAR_PART_NSET is not exactly 1..far_nodes")
    fsi_ids = far["sets"]["FSI_NSET"]
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
        "node_id_gaps": far["node_id_gaps"],
        "node_id_max": far["node_id_max"],
        "elements": contract.far_elems,
        "element_type_counts": {name: int(got[name])
                                for name in contract.far_elem_types},
        "c3d8_elements": int(got.get("C3D8", 0)),
        "fsi_nodes": int(fsi_ids.shape[0]),
        "ignored_sets": list(far["ignored_sets"]),
        "outer_set_sizes": {name: int(far["sets"][f"{name}_NSET"].shape[0])
                            for name in contract.outer_sets},
        "fixed_set_sizes": {name: int(values.shape[0])
                            for name, values in fixed.items()},
        "element_id_min": elem_min,
        "element_id_max": elem_max,
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
    contract = resolve_contract(far, contract)
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
        element_id = 0
        for etype in contract.far_elem_types:
            rows = [block for block in far["blocks"] if block["type"] == etype]
            if not sum(int(block["conn"].shape[0]) for block in rows):
                continue
            dst.write(f"*ELEMENT, TYPE={etype}, ELSET={contract.far_elset}\n")
            for block in rows:
                for element in block["conn"]:
                    element_id += 1
                    dst.write(f"{element_id}," +
                              ",".join(str(int(v)) for v in element) + "\n")
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

    far_type_counts: dict[str, int] = {}
    for block in far["blocks"]:
        far_type_counts[block["type"]] = (
            far_type_counts.get(block["type"], 0) + int(block["ids"].shape[0]))
    c3d8 = far_type_counts.get("C3D8", 0)
    c3d6 = far_type_counts.get("C3D6", 0)
    counts: dict[str, JsonValue] = {
        "nodes": contract.total_nodes,
        "c3d8": c3d8,
        "c3d10": contract.wing_elems,
        "elements": contract.total_elems,
    }
    if c3d6:
        counts["c3d6"] = c3d6
    element_numbering: dict[str, JsonValue] = {"c3d8": f"1..{c3d8}"}
    if c3d6:
        element_numbering["c3d6"] = f"{c3d8 + 1}..{contract.far_elems}"
    element_numbering["c3d10"] = (f"{contract.elem_offset + 1}.."
                                  f"{contract.total_elems}")
    report: Report = {
        "stage": "01_assemble_mesh",
        "farfield_profile": contract.profile,
        "sources": {"farfield": str(farfield_deck.resolve()),
                    "wing": str(wing_deck.resolve())},
        "sha256": {"farfield": sha256(farfield_deck), "wing": sha256(wing_deck)},
        "counts": counts,
        "node_numbering": {"far": f"1..{contract.far_nodes}",
                           "wing": f"{contract.node_offset + 1}..{contract.total_nodes}"},
        "element_numbering": element_numbering,
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
    parser.add_argument("--farfield-deck", type=Path, default=None,
                        help="far-field Abaqus deck with C3D8 (and optional C3D6) "
                             "elements; required for --farfield-profile fluid, "
                             "defaults to mesh/farfield_c3d8.zip for legacy")
    parser.add_argument("--farfield-profile", choices=sorted(PROFILES),
                        default="legacy",
                        help="named far-field profile: legacy (pinned) or fluid "
                             "(derived from the deck) (default: legacy)")
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
    if args.farfield_deck is None:
        if args.farfield_profile != "legacy":
            parser.error(f"--farfield-profile {args.farfield_profile} requires "
                         "--farfield-deck (no committed fluid deck default)")
        args.farfield_deck = DEFAULT_FARFIELD_DECK
    meta = args.meta if args.meta is not None else args.output.with_name(
        args.output.stem + "_meta.json")
    try:
        report = assemble(args.farfield_deck, args.wing_deck, args.output, meta,
                          contract=PROFILES[args.farfield_profile],
                          tie_tolerance=args.tie_tolerance)
    except MergeError as error:
        print(f"[AGARD ERR] {error}", file=sys.stderr)
        return 1
    counts = report["counts"]
    distance = report["fsi_wing_distance"]
    print(f"[AGARD 01] merged={args.output} profile={report['farfield_profile']} "
          f"nodes={counts['nodes']} c3d8={counts['c3d8']} "
          f"c3d6={counts.get('c3d6', 0)} c3d10={counts['c3d10']} "
          f"fsi_max_dist={distance['max_m']:.6g}m meta={meta}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
