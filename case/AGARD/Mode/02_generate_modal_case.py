#!/usr/bin/env python3
"""Stage 02: generate the combined AGARD modal CalculiX deck.

Consumes the stage-01 merged mesh deck (far-field C3D8 + wing C3D10) and emits
the complete modal deck: element-type-aware FSI tie surfaces, the fixed
boundaries, the calibrated orthotropic wing material and the pseudo-elastic
far-field materials.

Element types
-------------
* far-field : C3D8 faces are emitted when all four face nodes lie in the
  far-field inner set (``FSI_NSET``); the six quad faces are CalculiX S1..S6.
* wing      : C3D10 faces are emitted when all three CORNER nodes lie in the
  wing ``SIDE`` set, using the CalculiX tet face table
  ``S1=(1,2,3) S2=(1,4,2) S3=(2,4,3) S4=(3,4,1)``; the mid-edge nodes belong
  to the face implicitly.

Wing material (values from ``../CaseDescription.md``)
-----------------------------------------------------
``soft`` -> ``soft_anchored`` (Model No.3) / ``hard`` -> ``hard_e1_0p75``
(Model No.1); ``*ELASTIC, TYPE=ENGINEERING CONSTANTS`` with the single
``*ORIENTATION`` data line ``1,1,0, 0,0,1`` (axis1 span -> E1, axis2 thickness
-> E2, axis3 chord -> E3).

Far-field rule (bounded; see ``README.md``)
-------------------------------------------
``E(V) = clip(E_max * (V_min/V)**p, E_max/floor, E_max)`` with
``E_max = factor * E_ref`` (default ``1e-6 * E1``), one material + elset +
solid section per volume-ranked layer, ``*DENSITY 0.`` and Poisson ``0.`` for
every layer.  A stiff grounded pseudo-elastic medium stiffens the wing
(measured +594.8 % at a uniform 1e5 Pa), so the band is capped and floored:
the emitted span can never collapse to the old 2.4e16 contrast.

Module size note: ``allow: SIZE_OK`` -- one CLI stage owning parsing, layering,
deck emission and read-back verification so the artifact contract is auditable
in a single place.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import shutil
import subprocess
import sys
import time
import zipfile
from array import array
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Iterator, TextIO

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]

type JsonValue = str | int | float | bool | None | list["JsonValue"] | dict[str, "JsonValue"]
type Report = dict[str, JsonValue]

SCRIPT_DIR: Final = Path(__file__).resolve().parent
DEFAULT_WING_SOURCE: Final = SCRIPT_DIR / "mesh" / "wing_c3d10.zip"
DEFAULT_TIE_TOLERANCE: Final = 0.05
DEFAULT_MODES: Final = 10
DEFAULT_FAR_FIELD_MODULUS: Final = 28.0
DEFAULT_LAYERS: Final = 20
DEFAULT_LAYER_P: Final = 0.25
DEFAULT_LAYER_FLOOR: Final = 1000.0
DEFAULT_EMAX_FACTOR: Final = 1e-6
FAR_FIELD_ELSET: Final = "FAR_PART_ELSET"
FAR_FIELD_NAME: Final = "PFLUID"
WING_ELSET: Final = "WING"
BOUNDARY_SETS: Final = ("FIXED_FRONT", "FIXED_BACK", "FIXED_IN", "FIXED_OUT",
                        "FIXED_TOP", "FIXED_BOTTOM")
FIXED_SETS: Final = BOUNDARY_SETS + ("ROOT",)
REF_MODULUS_KEYS: Final = {"E1": "E1", "E2": "E2", "E3": "E3"}

# CalculiX hexahedral face table (0-based indices).
HEX_FACES: Final = (("S1", (0, 1, 2, 3)), ("S2", (4, 5, 6, 7)),
                    ("S3", (0, 1, 5, 4)), ("S4", (1, 2, 6, 5)),
                    ("S5", (2, 3, 7, 6)), ("S6", (3, 0, 4, 7)))
# CalculiX C3D10 face table by CORNER nodes: S1=(1,2,3) S2=(1,4,2)
# S3=(2,4,3) S4=(3,4,1); mid-edge nodes belong to the face implicitly.
TET_FACES: Final = (("S1", (0, 1, 2)), ("S2", (0, 3, 1)),
                    ("S3", (1, 3, 2)), ("S4", (2, 3, 0)))

FORBIDDEN_CARD_TOKENS: Final = ("TYPE=ORTHO", "STORAGE=", "SOLVER=", "*MPC")
ALLOWED_CARD_FAMILIES: Final = frozenset({
    "*HEADING", "*NODE", "*ELEMENT", "*NSET", "*SURFACE", "*TIE", "*MATERIAL",
    "*ELASTIC", "*DENSITY", "*ORIENTATION", "*SOLID", "*ELSET", "*STEP",
    "*FREQUENCY", "*BOUNDARY", "*END"})

# 2x2x2 Gauss integration for the far-field brick volumes.
NATURAL: Final = np.array(
    [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
     [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=np.float64)
DERIVATIVES: Final = tuple(
    np.array([[0.125 * a * (1 + b * e) * (1 + c * z),
               0.125 * b * (1 + a * x) * (1 + c * z),
               0.125 * c * (1 + a * x) * (1 + b * e)] for a, b, c in NATURAL],
             dtype=np.float64)
    for x, e, z in NATURAL / np.sqrt(3.0))


@dataclass(frozen=True, slots=True)
class Material:
    """One calibrated orthotropic wing material card."""

    name: str
    E1: float
    E2: float
    E3: float
    nu12: float
    nu13: float
    nu23: float
    G12: float
    G13: float
    G23: float
    rho: float
    orientation: str = "1,1,0, 0,0,1"

    def stability_ratio(self) -> float:
        """``E3/E2``; must stay below ``1/nu23**2`` for a positive-definite law."""
        return self.E3 / self.E2

    def stability_limit(self) -> float:
        return 1.0 / self.nu23**2


MATERIALS: Final = {
    "soft": Material(name="soft_anchored", E1=3.518452897e9, E2=2.673062111e8,
                     E3=1.313278384e9, nu12=0.034, nu13=0.033, nu23=0.326,
                     G12=3.904665311e8, G13=4.201475033e8, G23=9.268418440e7,
                     rho=413.5),
    "hard": Material(name="hard_e1_0p75", E1=1.013605240e10, E2=7.295992132e8,
                     E3=1.738397081e9, nu12=0.034, nu13=0.033, nu23=0.326,
                     G12=2.341301081e8, G13=5.652122962e8, G23=8.280724206e8,
                     rho=495.5),
}

type WingType = str
type FarFieldMode = str


@dataclass(frozen=True, slots=True)
class CaseSpec:
    """User controls for one modal-deck generation."""

    wing_type: WingType = "soft"
    modes: int = DEFAULT_MODES
    tie_tolerance: float = DEFAULT_TIE_TOLERANCE
    far_field_mode: FarFieldMode = "layered"
    layer_emax_factor: float = DEFAULT_EMAX_FACTOR
    layer_ref_modulus: str = "E1"
    layer_p: float = DEFAULT_LAYER_P
    layer_floor: float = DEFAULT_LAYER_FLOOR
    layers: int = DEFAULT_LAYERS
    far_field_modulus: float = DEFAULT_FAR_FIELD_MODULUS

    def material(self) -> Material:
        """Resolve ``wing_type`` to its calibrated material."""
        try:
            return MATERIALS[self.wing_type]
        except KeyError as error:
            raise DeckError(f"unknown wing type {self.wing_type!r}: "
                            f"expected one of {tuple(MATERIALS)}") from error


@dataclass(frozen=True, slots=True)
class LayerRow:
    """One emitted volume-ranked far-field layer."""

    index: int
    count: int
    vmin: float
    vmax: float
    vmean_log: float
    E: float
    floored: bool


@dataclass(frozen=True, slots=True)
class LayeredConfig:
    """The bounded far-field band and its element assignment."""

    factor: float
    ref: str
    ref_modulus_Pa: float
    p: float
    floor_div: float
    layers: int
    E_max: float
    E_min: float
    contrast: float
    floored_layers: int
    table: tuple[LayerRow, ...]
    assignment: IntArray


@dataclass(frozen=True, slots=True)
class MergedMesh:
    """The stage-01 merged deck parsed into arrays."""

    coords: FloatArray
    far_conn: IntArray
    wing_conn: IntArray
    sets: dict[str, IntArray]
    total_nodes: int
    far_nodes: int
    far_elems: int


@dataclass(frozen=True, slots=True)
class TieSurfaces:
    """Sorted FSI tie faces: element rows plus CalculiX face flags."""

    fluid_rows: IntArray
    fluid_flags: IntArray
    wing_rows: IntArray
    wing_flags: IntArray

    def write(self, output: TextIO, tie_tolerance: float, far_elems: int) -> None:
        """Emit both ``*SURFACE`` blocks in element/face order."""
        output.write("*SURFACE, NAME=SURF_FSI_FLUID, TYPE=ELEMENT\n")
        for row, flag in zip(self.fluid_rows.tolist(), self.fluid_flags.tolist(),
                             strict=True):
            output.write(f"{row + 1},S{flag}\n")
        output.write("*SURFACE, NAME=SURF_FSI_WING, TYPE=ELEMENT\n")
        for row, flag in zip(self.wing_rows.tolist(), self.wing_flags.tolist(),
                             strict=True):
            output.write(f"{far_elems + row + 1},S{flag}\n")


@dataclass(frozen=True, slots=True)
class PureWing:
    """Requested pure-wing reference deck (built from the wing source deck)."""

    source: Path
    output: Path


@dataclass(frozen=True, slots=True)
class BuildRequest:
    """Everything one ``generate`` call needs."""

    merged_deck: Path
    output: Path
    meta: Path
    spec: CaseSpec
    pure: PureWing | None = None
    run_ccx: bool = False
    run_pure_ccx: bool = False
    ccx_command: str = "ccx"


@dataclass(frozen=True, slots=True)
class LayerCard:
    """A layer material card parsed back from a generated deck."""

    name: str
    elastic: tuple[float, ...]
    density: str
    ids: IntArray


class DeckError(RuntimeError):
    """A deck or control set violates the modal-case contract."""


def sha256(path: Path) -> str:
    """Streaming sha256 of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_spec(spec: CaseSpec) -> None:
    """Reject control sets that would emit an invalid or unbounded band."""
    if spec.wing_type not in MATERIALS:
        raise DeckError(f"--wing-type must be one of {tuple(MATERIALS)}, "
                        f"got {spec.wing_type!r}")
    if spec.modes < 1:
        raise DeckError("--modes must be positive")
    if not spec.tie_tolerance > 0.0:
        raise DeckError("--tie-tolerance must be positive")
    if spec.far_field_mode not in ("layered", "uniform"):
        raise DeckError(f"--far-field-mode must be layered or uniform, "
                        f"got {spec.far_field_mode!r}")
    if not spec.layer_emax_factor > 0.0:
        raise DeckError("--layer-emax-factor must be positive")
    if spec.layer_ref_modulus not in REF_MODULUS_KEYS:
        raise DeckError(f"--layer-ref-modulus must be one of "
                        f"{tuple(REF_MODULUS_KEYS)}, got {spec.layer_ref_modulus!r}")
    if not spec.layer_p > 0.0:
        raise DeckError("--layer-p must be positive")
    if not spec.layer_floor > 1.0:
        raise DeckError("--layer-floor must be greater than 1")
    if spec.layers < 1:
        raise DeckError("--layers must be at least 1")
    if spec.far_field_mode == "uniform" and not spec.far_field_modulus > 0.0:
        raise DeckError("--far-field-modulus must be positive")


def parse_merged(deck: Path) -> MergedMesh:
    """Stream the merged mesh deck into arrays, asserting its contiguity contract."""
    coords = array("d")
    far_conn = array("q")
    wing_conn = array("q")
    sets: dict[str, array] = {}
    current_set: str | None = None
    section: str | None = None
    node_count = 0
    far_rows = 0
    wing_rows = 0
    with deck.open("r", encoding="ascii") as source:
        for line in source:
            if line.startswith("*"):
                head = line.replace(" ", "").strip().upper()
                current_set = None
                if head == "*NODE" or head.startswith("*NODE,"):
                    section = "node"
                elif head.startswith("*ELEMENT"):
                    if "TYPE=C3D8" in head:
                        section = "far"
                    elif "TYPE=C3D10" in head:
                        section = "wing"
                    else:
                        section = None
                elif head.startswith("*NSET"):
                    name = head.split("NSET=")[1].split(",")[0]
                    current_set = name
                    sets[name] = array("q")
                    section = "nset"
                else:
                    section = None
                continue
            if section == "node":
                fields = line.split(",")
                node_count += 1
                if int(fields[0]) != node_count:
                    raise DeckError(f"node ids are not contiguous at line {node_count}")
                coords.append(float(fields[1]))
                coords.append(float(fields[2]))
                coords.append(float(fields[3]))
            elif section == "far":
                fields = line.split(",")
                far_rows += 1
                if int(fields[0]) != far_rows:
                    raise DeckError(f"far element ids are not contiguous at {far_rows}")
                far_conn.extend(int(value) for value in fields[1:9])
            elif section == "wing":
                fields = line.split(",")
                wing_rows += 1
                expected = len(far_conn) // 8 + wing_rows
                if int(fields[0]) != expected:
                    raise DeckError(f"wing element ids are not contiguous at {expected}")
                wing_conn.extend(int(value) for value in fields[1:11])
            elif section == "nset" and current_set is not None:
                sets[current_set].extend(
                    int(value) for value in line.strip().split(",") if value.strip())
    far_conn_np = np.asarray(far_conn, dtype=np.int64).reshape(-1, 8)
    wing_conn_np = np.asarray(wing_conn, dtype=np.int64).reshape(-1, 10)
    raw_coords = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
    total_nodes = int(raw_coords.shape[0])
    # Row 0 is padding: coordinates are indexed by node id, not by row.
    coords_np = np.zeros((total_nodes + 1, 3), dtype=np.float64)
    coords_np[1:] = raw_coords
    sets_np = {name: np.asarray(values, dtype=np.int64)
               for name, values in sets.items()}
    for required in ("FAR_PART_NSET", "FSI_NSET", "SIDE", "WING", *FIXED_SETS):
        if required not in sets_np:
            raise DeckError(f"merged deck is missing NSET={required}")
    far_part = sets_np["FAR_PART_NSET"]
    far_nodes = int(far_part.shape[0])
    if not np.array_equal(np.sort(far_part), np.arange(1, far_nodes + 1)):
        raise DeckError("FAR_PART_NSET is not exactly 1..far_nodes")
    wing_set = sets_np["WING"]
    if not np.array_equal(np.sort(wing_set),
                          np.arange(far_nodes + 1, total_nodes + 1)):
        raise DeckError("WING NSET is not the contiguous wing node range")
    if far_conn_np.min() < 1 or far_conn_np.max() > far_nodes:
        raise DeckError("far connectivity escapes the far-field node range")
    if wing_conn_np.min() <= far_nodes or wing_conn_np.max() > total_nodes:
        raise DeckError("wing connectivity escapes the wing node range")
    return MergedMesh(coords=coords_np, far_conn=far_conn_np,
                      wing_conn=wing_conn_np, sets=sets_np,
                      total_nodes=total_nodes, far_nodes=far_nodes,
                      far_elems=int(far_conn_np.shape[0]))


def extract_surfaces(mesh: MergedMesh) -> TieSurfaces:
    """Element-type-aware FSI tie faces on both sides of the interface."""
    in_fsi = np.zeros(mesh.total_nodes + 1, dtype=bool)
    in_fsi[mesh.sets["FSI_NSET"]] = True
    fluid_rows: list[IntArray] = []
    fluid_flags: list[IntArray] = []
    for flag, (_name, indices) in enumerate(HEX_FACES, start=1):
        nodes = mesh.far_conn[:, list(indices)]
        rows = np.nonzero(in_fsi[nodes].all(axis=1))[0]
        fluid_rows.append(rows)
        fluid_flags.append(np.full(rows.shape[0], flag, dtype=np.int64))
    fluid_row = np.concatenate(fluid_rows)
    fluid_flag = np.concatenate(fluid_flags)
    fluid_order = np.lexsort((fluid_flag, fluid_row))
    in_side = np.zeros(mesh.total_nodes + 1, dtype=bool)
    in_side[mesh.sets["SIDE"]] = True
    wing_rows: list[IntArray] = []
    wing_flags: list[IntArray] = []
    for flag, (_name, indices) in enumerate(TET_FACES, start=1):
        nodes = mesh.wing_conn[:, list(indices)]
        rows = np.nonzero(in_side[nodes].all(axis=1))[0]
        wing_rows.append(rows)
        wing_flags.append(np.full(rows.shape[0], flag, dtype=np.int64))
    wing_row = np.concatenate(wing_rows)
    wing_flag = np.concatenate(wing_flags)
    wing_order = np.lexsort((wing_flag, wing_row))
    surfaces = TieSurfaces(fluid_rows=fluid_row[fluid_order],
                           fluid_flags=fluid_flag[fluid_order],
                           wing_rows=wing_row[wing_order],
                           wing_flags=wing_flag[wing_order])
    if surfaces.fluid_rows.shape[0] == 0 or surfaces.wing_rows.shape[0] == 0:
        raise DeckError("no faces found on one side of the FSI tie")
    return surfaces


def hexa_volumes(coords: FloatArray, conn: IntArray) -> FloatArray:
    """2x2x2 Gauss volumes of eight-node bricks."""
    volumes = np.zeros(conn.shape[0], dtype=np.float64)
    points = coords[conn]
    for derivatives in DERIVATIVES:
        volumes += np.abs(np.linalg.det(
            np.einsum("ia,mij->maj", derivatives, points)))
    return volumes


def rank_layers(volumes: FloatArray, layers: int) -> IntArray:
    """Assign every element a layer index by ascending volume rank.

    Contiguous chunks of the stable volume sort; the first ``N % layers``
    chunks take one extra element, so no layer can be empty.
    """
    if layers < 1:
        raise DeckError(f"layers must be >= 1, got {layers}")
    order = np.argsort(volumes, kind="stable")
    assignment = np.empty(volumes.shape[0], dtype=np.int64)
    base, extra = divmod(volumes.shape[0], layers)
    start = 0
    for index in range(layers):
        stop = start + base + (1 if index < extra else 0)
        assignment[order[start:stop]] = index
        start = stop
    if start != volumes.shape[0]:
        raise DeckError("layer chunks do not cover every element")
    return assignment


def layer_table(volumes: FloatArray, assignment: IntArray, p: float,
                floor_div: float, emax: float) -> tuple[LayerRow, ...]:
    """Emitted layers: E is the pointwise law at the layer's smallest volume."""
    vmin_global = float(volumes.min())
    emin = emax / floor_div
    rows: list[LayerRow] = []
    for index in range(int(assignment.max()) + 1):
        values = volumes[assignment == index]
        if values.shape[0] == 0:
            raise DeckError(f"empty layer {index}")
        low = float(values.min())
        pointwise = emax * (vmin_global / low) ** p
        modulus = min(emax, max(emin, pointwise))
        rows.append(LayerRow(index=index, count=int(values.shape[0]), vmin=low,
                             vmax=float(values.max()),
                             vmean_log=float(np.exp(np.log(values).mean())),
                             E=float(modulus),
                             floored=bool(pointwise < emin)))
    moduli = [row.E for row in rows]
    if not all(later <= earlier for earlier, later in zip(moduli, moduli[1:])):
        raise DeckError(f"layer moduli are not non-increasing: {moduli}")
    if abs(max(moduli) - emax) > abs(emax) * 1e-12:
        raise DeckError(f"emitted max modulus {max(moduli)} != cap {emax}")
    if min(moduli) < emin * (1.0 - 1e-12):
        raise DeckError(f"emitted min modulus {min(moduli)} is below floor {emin}")
    return tuple(rows)


def layered_config(spec: CaseSpec, volumes: FloatArray) -> LayeredConfig:
    """Build the bounded layer band; asserts the cap and the span bound."""
    material = spec.material()
    ref_modulus = float(getattr(material, REF_MODULUS_KEYS[spec.layer_ref_modulus]))
    emax = spec.layer_emax_factor * ref_modulus
    emin = emax / spec.layer_floor
    if abs(emax - spec.layer_emax_factor * ref_modulus) > abs(emax) * 1e-12:
        raise DeckError("E_max != factor * E_ref")
    assignment = rank_layers(volumes, spec.layers)
    table = layer_table(volumes, assignment, spec.layer_p, spec.layer_floor, emax)
    moduli = [row.E for row in table]
    span = max(moduli) / min(moduli)
    if span > spec.layer_floor * (1.0 + 1e-12):
        raise DeckError(f"layer modulus span {span:.6g} exceeds floor "
                        f"{spec.layer_floor:.6g}")
    return LayeredConfig(factor=spec.layer_emax_factor, ref=spec.layer_ref_modulus,
                         ref_modulus_Pa=ref_modulus, p=spec.layer_p,
                         floor_div=spec.layer_floor, layers=spec.layers,
                         E_max=emax, E_min=emin, contrast=span,
                         floored_layers=sum(1 for row in table if row.floored),
                         table=table, assignment=assignment)


def wing_material_card(material: Material) -> str:
    """The calibrated orthotropic wing material + solid section, byte-stable."""
    elastic = ", ".join(f"{value:.10e}" for value in
                        (material.E1, material.E2, material.E3, material.nu12,
                         material.nu13, material.nu23, material.G12,
                         material.G13))
    return (
        f"*MATERIAL, NAME={material.name}\n"
        "*ELASTIC, TYPE=ENGINEERING CONSTANTS\n"
        f"{elastic}\n"
        f"{material.G23:.10e}, 0.\n"
        "*DENSITY\n"
        f"{material.rho:.10e},\n"
        "*ORIENTATION, NAME=WING_AXES\n"
        f"{material.orientation}\n"
        f"*SOLID SECTION, ELSET={WING_ELSET}, MATERIAL={material.name}, "
        "ORIENTATION=WING_AXES\n,\n"
    )


def uniform_farfield_card(modulus: float) -> str:
    """ONE uniform massless isotropic material over the whole far-field set."""
    if not modulus > 0.0:
        raise DeckError(f"far-field modulus must be positive, got {modulus!r}")
    return (
        f"*MATERIAL, NAME={FAR_FIELD_NAME}\n"
        "*ELASTIC\n"
        f"{modulus:.9e},0.\n"
        "*DENSITY\n"
        "0.\n"
        f"*SOLID SECTION, ELSET={FAR_FIELD_ELSET}, MATERIAL={FAR_FIELD_NAME}\n"
    )


def layered_slot_header(config: LayeredConfig) -> str:
    """Comment block describing the emitted bounded layered band."""
    lines = [
        f"** far-field: {config.layers} massless pseudo-elastic layers on "
        f"{FAR_FIELD_ELSET} (bounded layered rule)",
        f"** E_max={config.E_max:.9e} Pa = "
        f"{config.E_max / config.ref_modulus_Pa:g}*{config.ref}="
        f"{config.ref_modulus_Pa:.9e} Pa; E(V)=E_max*(V_min/V)^p with "
        f"p={config.p:g}, floor_div={config.floor_div:g} -> "
        f"E_min={config.E_min:.9e} Pa",
        f"** contrast={config.contrast:.6g} (bounded; the reference 2.4e16 "
        "band is abandoned)",
    ]
    for row in config.table:
        lines.append(
            f"** layer={row.index} count={row.count} "
            f"V=[{row.vmin:.9e},{row.vmax:.9e}] E={row.E:.9e} Pa"
            + (" [floored]" if row.floored else ""))
    return "\n".join(lines) + "\n"


def layered_cards(config: LayeredConfig) -> str:
    """One *MATERIAL + *ELSET + *SOLID SECTION per bounded layer."""
    lines: list[str] = []
    for row in config.table:
        name = f"PFLUID_L{row.index:02d}"
        ids = np.nonzero(config.assignment == row.index)[0] + 1
        if ids.shape[0] != row.count:
            raise DeckError(f"{name} id count {ids.shape[0]} != {row.count}")
        lines.append(f"*MATERIAL, NAME={name}")
        lines.append("*ELASTIC")
        lines.append(f"{row.E:.9e},0.")
        lines.append("*DENSITY")
        lines.append("0.")
        lines.append(f"*ELSET, ELSET={name}")
        for offset in range(0, ids.shape[0], 16):
            lines.append(",".join(str(int(value))
                                  for value in ids[offset:offset + 16]))
        lines.append(f"*SOLID SECTION, ELSET={name}, MATERIAL={name}")
    return "\n".join(lines) + "\n"


def material_slot(spec: CaseSpec, config: LayeredConfig | None) -> str:
    """The full material region, markers included, exactly as the reference."""
    material = spec.material()
    if spec.far_field_mode == "uniform":
        head = (
            "** === MATERIAL SLOT (rewritten by "
            "case/AGARD/Mode/02_generate_modal_case.py) ===\n"
            f"** far-field: ONE uniform massless pseudo-elastic medium "
            f"({FAR_FIELD_NAME}) on {FAR_FIELD_ELSET}\n"
            f"** E={spec.far_field_modulus:.9e} Pa, nu=0.0, rho=0.0 (uniform "
            "rule; the 20-layer reference band is abandoned)\n")
        return (head + wing_material_card(material)
                + uniform_farfield_card(spec.far_field_modulus)
                + "** === END MATERIAL SLOT ===\n")
    if config is None:
        raise DeckError("layered mode requires a layer configuration")
    return (
        "** === MATERIAL SLOT (rewritten by "
        "case/AGARD/Mode/02_generate_modal_case.py) ===\n"
        + layered_slot_header(config)
        + wing_material_card(material)
        + layered_cards(config)
        + "** === END MATERIAL SLOT ===\n")


def step_card(spec: CaseSpec) -> str:
    """``*STEP`` / ``*FREQUENCY`` / ``*NODE FILE`` / ``*BOUNDARY`` / ``*END STEP``."""
    lines = ["*STEP", "*FREQUENCY", f"{spec.modes},0.334,1.E6",
             "*NODE FILE", "U", "*BOUNDARY"]
    lines.extend(f"{name},1,3" for name in FIXED_SETS)
    lines.extend(["*END STEP", ""])
    return "\n".join(lines)


def write_combined(deck: Path, request: BuildRequest, mesh: MergedMesh,
                   surfaces: TieSurfaces,
                   config: LayeredConfig | None) -> None:
    """Copy the merged deck verbatim, then append surfaces, tie, materials, step."""
    request.output.parent.mkdir(parents=True, exist_ok=True)
    with deck.open("rb") as source, request.output.open("wb") as sink:
        shutil.copyfileobj(source, sink, 1 << 20)
    with request.output.open("a", encoding="ascii") as output:
        surfaces.write(output, request.spec.tie_tolerance, mesh.far_elems)
        output.write(
            f"*TIE, NAME=FSI_COUPLING, POSITION TOLERANCE="
            f"{request.spec.tie_tolerance:.9g}, ADJUST=NO\n"
            "SURF_FSI_FLUID,SURF_FSI_WING\n")
        output.write(material_slot(request.spec, config))
        output.write(step_card(request.spec))


@dataclass(frozen=True, slots=True)
class DeckFacts:
    """Read-back evidence collected from one generated deck."""

    material_names: tuple[str, ...]
    card_families: frozenset[str]
    boundary_cards: frozenset[str]
    frequency_lines: int
    tie_line: str | None
    slot_text: str
    forbidden_hits: tuple[str, ...]


def deck_facts(path: Path) -> DeckFacts:
    """Stream the generated deck once, collecting the read-back evidence."""
    material_names: list[str] = []
    card_families: set[str] = set()
    boundary_cards: set[str] = set()
    slot_lines: list[str] = []
    forbidden_hits: list[str] = []
    frequency_lines = 0
    tie_line: str | None = None
    in_slot = False
    in_boundary = False
    with path.open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith("** === MATERIAL SLOT"):
                in_slot = True
                continue
            if line.startswith("** === END MATERIAL SLOT"):
                in_slot = False
                continue
            if line.startswith("*") and not line.startswith("**"):
                family = line.split(",")[0].split()[0].upper()
                card_families.add(family)
                in_boundary = family == "*BOUNDARY"
                if family == "*FREQUENCY":
                    frequency_lines += 1
                elif family == "*TIE":
                    tie_line = line.strip()
                elif family == "*MATERIAL":
                    material_names.append(
                        line.split("NAME=", 1)[1].split(",")[0].strip())
            elif in_boundary and "," in line:
                name = line.split(",", 1)[0].strip()
                if name in FIXED_SETS:
                    boundary_cards.add(name)
            if in_slot:
                slot_lines.append(line)
            for token in FORBIDDEN_CARD_TOKENS:
                if token in line:
                    forbidden_hits.append(token)
    return DeckFacts(material_names=tuple(material_names),
                     card_families=frozenset(card_families),
                     boundary_cards=frozenset(boundary_cards),
                     frequency_lines=frequency_lines, tie_line=tie_line,
                     slot_text="".join(slot_lines),
                     forbidden_hits=tuple(dict.fromkeys(forbidden_hits)))


def parse_layered_slot(text: str) -> list[LayerCard]:
    """Parse the emitted layer cards: material, elastic E, density, ids."""
    lines = text.splitlines()
    cards: list[LayerCard] = []
    index = 0
    while index < len(lines):
        match = re.match(r"\*MATERIAL,\s*NAME=(\S+)", lines[index])
        if match and match.group(1).startswith("PFLUID_L"):
            name = match.group(1)
            if lines[index + 1].strip() != "*ELASTIC":
                raise DeckError(f"{name}: expected *ELASTIC")
            elastic = tuple(float(value) for value in
                            lines[index + 2].replace(" ", "").split(",") if value)
            if lines[index + 3].strip() != "*DENSITY":
                raise DeckError(f"{name}: expected *DENSITY")
            density = lines[index + 4].strip()
            elset = re.match(r"\*ELSET,\s*ELSET=(\S+)", lines[index + 5].strip())
            if elset is None or elset.group(1) != name:
                raise DeckError(f"{name}: expected *ELSET, ELSET={name}")
            ids: list[int] = []
            cursor = index + 6
            while cursor < len(lines) and not lines[cursor].startswith("*"):
                ids.extend(int(value) for value in
                           lines[cursor].replace(" ", "").split(",") if value)
                cursor += 1
            section = re.match(
                r"\*SOLID SECTION,\s*ELSET=(\S+),\s*MATERIAL=(\S+)",
                lines[cursor].strip())
            if section is None or section.group(1) != name or section.group(2) != name:
                raise DeckError(f"{name}: unexpected *SOLID SECTION")
            cards.append(LayerCard(name=name, elastic=elastic, density=density,
                                   ids=np.asarray(ids, dtype=np.int64)))
            index = cursor + 1
        else:
            index += 1
    return cards


def verify_deck(path: Path, spec: CaseSpec, config: LayeredConfig | None) -> Report:
    """Assert the written deck carries exactly the configured artifact."""
    material = spec.material()
    facts = deck_facts(path)
    if facts.frequency_lines != 1:
        raise DeckError(f"expected one *FREQUENCY card, got {facts.frequency_lines}")
    if facts.forbidden_hits:
        raise DeckError(f"forbidden card tokens emitted: {facts.forbidden_hits}")
    unknown = facts.card_families - ALLOWED_CARD_FAMILIES
    if unknown:
        raise DeckError(f"unexpected card families emitted: {sorted(unknown)}")
    expected_tie = (f"*TIE, NAME=FSI_COUPLING, POSITION TOLERANCE="
                    f"{spec.tie_tolerance:.9g}, ADJUST=NO")
    if facts.tie_line != expected_tie:
        raise DeckError(f"tie card mismatch: {facts.tie_line!r} != {expected_tie!r}")
    for name in FIXED_SETS:
        if name not in facts.boundary_cards:
            raise DeckError(f"missing *BOUNDARY entry for {name}")
    names = list(facts.material_names)
    if not names or names[0] != material.name:
        raise DeckError(f"first material must be {material.name}, got {names[:1]}")
    if material.stability_ratio() >= material.stability_limit():
        raise DeckError(f"wing material E3/E2={material.stability_ratio():.6g} "
                        f"violates 1/nu23^2={material.stability_limit():.6g}")

    if spec.far_field_mode == "uniform":
        if names != [material.name, FAR_FIELD_NAME]:
            raise DeckError(f"unexpected material list for uniform mode: {names}")
        pattern = (
            rf"\*MATERIAL, NAME={FAR_FIELD_NAME}\n\*ELASTIC\n([^\n]*)\n"
            rf"\*DENSITY\n([^\n]*)\n\*SOLID SECTION, ELSET={FAR_FIELD_ELSET}, "
            rf"MATERIAL={FAR_FIELD_NAME}\n")
        match = re.search(pattern, facts.slot_text)
        if match is None:
            raise DeckError("uniform far-field material card shape unexpected")
        elastic = [float(value) for value in match.group(1).split(",")]
        density = [float(value) for value in match.group(2).split(",")]
        if abs(elastic[0] - spec.far_field_modulus) > spec.far_field_modulus * 1e-9:
            raise DeckError(f"far-field E {elastic[0]} != {spec.far_field_modulus}")
        if elastic[1] != 0.0:
            raise DeckError(f"far-field Poisson must be 0.0, got {elastic[1]}")
        if density != [0.0]:
            raise DeckError(f"far-field density must be 0., got {density}")
        return {"mode": "uniform", "far_field_elastic": elastic,
                "far_field_density": density, "material_names": names}

    if config is None:
        raise DeckError("layered mode requires a layer configuration")
    expected_names = [f"PFLUID_L{index:02d}" for index in range(config.layers)]
    if names != [material.name, *expected_names]:
        raise DeckError(f"unexpected material list for layered mode: {names}")
    cards = parse_layered_slot(facts.slot_text)
    if len(cards) != config.layers:
        raise DeckError(f"parsed {len(cards)} layers != {config.layers}")
    moduli = [card.elastic[0] for card in cards]
    if any(card.elastic[1] != 0.0 for card in cards):
        raise DeckError("every layer must have Poisson 0.0")
    if moduli != sorted(moduli, reverse=True):
        raise DeckError(f"layer moduli are not non-increasing: {moduli}")
    if abs(max(moduli) - config.E_max) > config.E_max * 1e-12:
        raise DeckError(f"emitted E_max {max(moduli)} != {config.E_max}")
    if min(moduli) < config.E_min * (1.0 - 1e-12):
        raise DeckError(f"emitted E_min {min(moduli)} is below {config.E_min}")
    if max(moduli) / min(moduli) > config.floor_div * (1.0 + 1e-12):
        raise DeckError("layer span exceeds the configured floor")
    if any(card.density != "0." for card in cards):
        raise DeckError("every layer density must be exactly 0.")
    seen = np.zeros(int(sum(row.count for row in config.table)) + 1, dtype=np.int64)
    for card, row in zip(cards, config.table, strict=True):
        if card.ids.shape[0] != row.count:
            raise DeckError(f"{card.name} has {card.ids.shape[0]} ids != {row.count}")
        if card.ids.min() < 1 or card.ids.max() > seen.shape[0] - 1:
            raise DeckError(f"{card.name} id out of range")
        seen[card.ids] += 1
    if not np.all(seen[1:] == 1):
        raise DeckError("layer ids are not an exact partition of the far-field set")
    return {"mode": "layered", "layer_moduli": moduli,
            "layer_counts": [int(card.ids.shape[0]) for card in cards],
            "layer_densities": [card.density for card in cards],
            "E_max": config.E_max, "E_min": config.E_min,
            "contrast": max(moduli) / min(moduli),
            "floored_layers": config.floored_layers,
            "material_names": names}


@contextmanager
def open_wing_source(source: Path) -> Iterator[TextIO]:
    """Open the wing source deck: plain ASCII, ``.gz`` or one-member ``.zip``."""
    if not source.is_file():
        raise DeckError(f"wing source not found: {source}")
    try:
        if source.suffix.lower() == ".zip":
            with zipfile.ZipFile(source) as archive:
                names = [name for name in archive.namelist()
                         if not name.endswith("/")]
                if len(names) != 1:
                    raise DeckError(
                        f"{source}: zip archive must hold exactly one deck "
                        f"member, found {len(names)}: {', '.join(names) or 'none'}")
                text = archive.read(names[0]).decode("ascii")
        elif source.suffix.lower() == ".gz":
            with gzip.open(source, "rt", encoding="ascii") as handle:
                text = handle.read()
        else:
            text = source.read_text(encoding="ascii")
    except (OSError, zipfile.BadZipFile, UnicodeDecodeError) as error:
        raise DeckError(f"cannot read wing source {source}: {error}") from error
    with io.StringIO(text) as stream:
        yield stream


def build_pure_deck(material: Material, source: Path, output: Path,
                    modes: int) -> Report:
    """Copy the wing source up to its placeholder material, then calibrate it."""
    root_count = 0
    in_root = False
    copied = False
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_wing_source(source) as src, \
            output.open("w", encoding="ascii") as dst:
        for line in src:
            if line.startswith("*"):
                head = line.replace(" ", "").upper()
                in_root = head.startswith("*NSET") and "NSET=ROOT" in head
                if head.startswith("*MATERIAL,"):
                    dst.write(wing_material_card(material))
                    dst.write("\n".join(["*STEP", "*FREQUENCY",
                                         f"{modes},0.334,1.E6", "*NODE FILE", "U",
                                         "*BOUNDARY", "ROOT,1,3", "*END STEP", ""]))
                    copied = True
                    break
            elif in_root:
                root_count += sum(1 for value in line.strip().split(",") if value.strip())
            dst.write(line)
    if not copied:
        raise DeckError(f"no *MATERIAL card found in wing source {source}")
    if root_count == 0:
        raise DeckError(f"no ROOT NSET found in wing source {source}")
    return {"case": "pure", "material": material.name, "modes": modes,
            "root_nodes": root_count, "source": str(source.resolve()),
            "output": str(output.resolve())}


def run_solver(command: str, deck: Path) -> Report:
    """Run one CalculiX job in the deck's directory; explicit request only."""
    subprocess.run([command, "-i", deck.stem], cwd=deck.parent, check=True)
    return {"command": command, "deck": str(deck.resolve()), "returncode": 0}


def generate(request: BuildRequest) -> Report:
    """Build the combined modal deck (and optional pure deck) and verify it."""
    started = time.time()
    validate_spec(request.spec)
    if not request.merged_deck.is_file():
        raise DeckError(f"merged deck not found: {request.merged_deck}")
    if request.run_pure_ccx and request.pure is None:
        raise DeckError("--run-pure-ccx requires --pure-wing")
    mesh = parse_merged(request.merged_deck)
    surfaces = extract_surfaces(mesh)
    volumes = hexa_volumes(mesh.coords[:mesh.far_nodes + 1], mesh.far_conn)
    config = (layered_config(request.spec, volumes)
              if request.spec.far_field_mode == "layered" else None)
    write_combined(request.merged_deck, request, mesh, surfaces, config)
    verified = verify_deck(request.output, request.spec, config)
    report: Report = {
        "stage": "02_generate_modal_case",
        "case": "combined",
        "input": str(request.merged_deck.resolve()),
        "output": {"path": str(request.output.resolve()),
                   "bytes": request.output.stat().st_size,
                   "sha256": sha256(request.output)},
        "wing_type": request.spec.wing_type,
        "material": request.spec.material().name,
        "modes": request.spec.modes,
        "frequency_line": f"{request.spec.modes},0.334,1.E6",
        "tie_tolerance": request.spec.tie_tolerance,
        "far_field_mode": request.spec.far_field_mode,
        "surfaces": {"fluid_faces": int(surfaces.fluid_rows.shape[0]),
                     "wing_faces": int(surfaces.wing_rows.shape[0])},
        "verified": verified,
        "step_card": step_card(request.spec).splitlines(),
    }
    if config is not None:
        report["layered"] = {
            "factor": config.factor, "ref": config.ref,
            "ref_modulus_Pa": config.ref_modulus_Pa, "p": config.p,
            "floor_div": config.floor_div, "layers": config.layers,
            "E_max_Pa": config.E_max, "E_min_Pa": config.E_min,
            "contrast": config.contrast, "floored_layers": config.floored_layers,
            "table": [{"index": row.index, "count": row.count, "vmin": row.vmin,
                       "vmax": row.vmax, "E": row.E, "floored": row.floored}
                      for row in config.table],
        }
        if abs(config.E_max - config.factor * config.ref_modulus_Pa) > \
                config.E_max * 1e-12:
            raise DeckError("E_max != factor * E_ref")
        if config.contrast > config.floor_div * (1.0 + 1e-12):
            raise DeckError("layer span exceeds the configured floor")
    else:
        report["far_field_modulus_Pa"] = request.spec.far_field_modulus
    if request.pure is not None:
        report["pure"] = build_pure_deck(request.spec.material(),
                                         request.pure.source,
                                         request.pure.output,
                                         request.spec.modes)
        report["pure"]["sha256"] = sha256(request.pure.output)
    if request.run_ccx:
        report["ccx"] = run_solver(request.ccx_command, request.output)
    if request.run_pure_ccx and request.pure is not None:
        report["pure_ccx"] = run_solver(request.ccx_command, request.pure.output)
    report["elapsed_s"] = round(time.time() - started, 1)
    request.meta.parent.mkdir(parents=True, exist_ok=True)
    request.meta.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="stage-01 merged mesh deck")
    parser.add_argument("--output", type=Path, required=True,
                        help="combined modal deck to write")
    parser.add_argument("--meta", type=Path, default=None,
                        help="metadata JSON (default: <output stem>_meta.json)")
    parser.add_argument("--modes", type=int, default=DEFAULT_MODES)
    parser.add_argument("--tie-tolerance", type=float,
                        default=DEFAULT_TIE_TOLERANCE)
    parser.add_argument("--wing-type", choices=tuple(MATERIALS), default="soft")
    parser.add_argument("--far-field-mode", choices=("layered", "uniform"),
                        default="layered",
                        help="bounded layered band (default) or one uniform "
                             "massless medium")
    parser.add_argument("--far-field-modulus", type=float,
                        default=DEFAULT_FAR_FIELD_MODULUS,
                        help="uniform-mode E in Pa")
    parser.add_argument("--layer-emax-factor", type=float,
                        default=DEFAULT_EMAX_FACTOR,
                        help="E_max = factor * reference modulus (default 1e-6)")
    parser.add_argument("--layer-ref-modulus", choices=tuple(REF_MODULUS_KEYS),
                        default="E1",
                        help="wing modulus the E_max cap references")
    parser.add_argument("--layer-p", type=float, default=DEFAULT_LAYER_P,
                        help="grading exponent in E(V)=E_max*(V_min/V)^p")
    parser.add_argument("--layer-floor", type=float, default=DEFAULT_LAYER_FLOOR,
                        help="E_min = E_max / floor (band contrast cap)")
    parser.add_argument("--layers", type=int, default=DEFAULT_LAYERS)
    parser.add_argument("--pure-wing", action="store_true",
                        help="also write the standalone pure-wing reference deck")
    parser.add_argument("--pure-output", type=Path, default=None,
                        help="pure-wing deck path (default: <output dir>/pure_<type>.inp)")
    parser.add_argument("--wing-source", type=Path, default=DEFAULT_WING_SOURCE,
                        help="second-order wing source deck for --pure-wing "
                             "(plain, .gz or single-member .zip)")
    parser.add_argument("--run-ccx", action="store_true",
                        help="run CalculiX on the combined deck (explicit only)")
    parser.add_argument("--run-pure-ccx", action="store_true",
                        help="run CalculiX on the pure-wing deck (requires --pure-wing)")
    parser.add_argument("--ccx-command", default="ccx")
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    if args.modes < 1:
        parser.error("--modes must be positive")
    if args.layers < 1:
        parser.error("--layers must be at least 1")
    if args.tie_tolerance <= 0:
        parser.error("--tie-tolerance must be positive")
    if args.layer_emax_factor <= 0:
        parser.error("--layer-emax-factor must be positive")
    if args.layer_p <= 0:
        parser.error("--layer-p must be positive")
    if args.layer_floor <= 1:
        parser.error("--layer-floor must be greater than 1")
    if args.run_pure_ccx and not args.pure_wing:
        parser.error("--run-pure-ccx requires --pure-wing")
    meta = args.meta if args.meta is not None else args.output.with_name(
        args.output.stem + "_meta.json")
    pure = None
    if args.pure_wing:
        pure_output = args.pure_output if args.pure_output is not None else \
            args.output.with_name(f"pure_{args.wing_type}.inp")
        pure = PureWing(source=args.wing_source, output=pure_output)
    spec = CaseSpec(wing_type=args.wing_type, modes=args.modes,
                    tie_tolerance=args.tie_tolerance,
                    far_field_mode=args.far_field_mode,
                    layer_emax_factor=args.layer_emax_factor,
                    layer_ref_modulus=args.layer_ref_modulus,
                    layer_p=args.layer_p, layer_floor=args.layer_floor,
                    layers=args.layers,
                    far_field_modulus=args.far_field_modulus)
    request = BuildRequest(merged_deck=args.input, output=args.output, meta=meta,
                           spec=spec, pure=pure, run_ccx=args.run_ccx,
                           run_pure_ccx=args.run_pure_ccx,
                           ccx_command=args.ccx_command)
    try:
        report = generate(request)
    except DeckError as error:
        print(f"[AGARD ERR] {error}", file=sys.stderr)
        return 1
    surfaces = report["surfaces"]
    print(f"[AGARD 02] deck={args.output} modes={args.modes} "
          f"wing={spec.wing_type} far-field={spec.far_field_mode} "
          f"fluid_faces={surfaces['fluid_faces']} "
          f"wing_faces={surfaces['wing_faces']} meta={meta}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
