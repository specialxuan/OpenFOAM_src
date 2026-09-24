#!/usr/bin/env python3
"""Stage 03: verify the CalculiX modal results, then export the FDM CSV files.

Gate order (NO CSV is written unless every gate passes)
-------------------------------------------------------
1. the ``.dat`` must hold exactly ``--modes`` strictly positive, strictly
   ascending frequencies;
2. every frequency must sit inside the structural sanity window
   (``--min-frequency .. --max-frequency``, default 1..1e4 Hz);
3. when ``--pure-dat`` is given, the combined and pure-wing frequencies are
   compared mode by mode and ``max |relative difference|`` must stay within
   ``--max-rel`` (default 0.01); the full table is always printed;
4. when ``--pure-frd`` is given, every combined mode must match the
   corresponding pure-wing mode's displacement field with a cosine of at least
   ``--min-cosine`` (default 0.9) and the best match must be diagonal;
5. only then are the CSVs written with the unchanged FDM contract: ascending
   global node id order, the three-field header ``frequency,count,modes,``
   with trailing commas, ``FluidNodeCoor.csv``/``StructNodeCoor.csv``,
   ``FluidNodeDisp{n}.csv``/``StructNodeDisp{n}.csv`` and the legacy
   ``FluidPara.csv``.

Node set names follow the merged deck: the fluid side is ``FAR_PART_NSET``
(alias ``FAR_PART``) and the structural side is ``WING`` (alias
``WING_PART``); ``FSI_NSET``/``ROOT``/``TIP``/``SIDE`` belong to the solver
deck and are not exported.

Module size note: ``allow: SIZE_OK`` -- one CLI stage whose gates and the
legacy export contract share the FRD reader, so they stay in one file.
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]

FLUID_SET_ALIASES: Final = ("FAR_PART_NSET", "FAR_PART")
STRUCT_SET_ALIASES: Final = ("WING", "WING_PART")
DEFAULT_MIN_FREQUENCY: Final = 1.0
DEFAULT_MAX_FREQUENCY: Final = 1.0e4
DEFAULT_MAX_REL: Final = 0.01
DEFAULT_MIN_COSINE: Final = 0.9
FREQUENCY = re.compile(r"100CL\s+\d+\s+([-+0-9.Ee]+)")


class ExportError(RuntimeError):
    """A gate failed or a result file violates the export contract."""


class GateError(ExportError):
    """The verification gate refuses to export."""


@dataclass(frozen=True, slots=True)
class ModalBlock:
    """One FRD displacement block: node ids and their vectors."""

    nodes: IntArray
    vectors: FloatArray


def read_set_ids(deck: Path, name: str) -> list[int]:
    """Ascending unique node ids of ``NSET=name`` from the modal deck."""
    values: list[int] = []
    active = False
    target = name.upper()
    for line in deck.read_text(encoding="ascii").splitlines():
        header = line.replace(" ", "").upper()
        if header.startswith("*NSET"):
            declared = [part.split("=", 1)[1] for part in header.split(",")
                        if part.startswith("NSET=")]
            active = target in declared
        elif line.startswith("*"):
            active = False
        elif active:
            values.extend(int(value) for value in line.split(",") if value.strip())
    if not values:
        raise ExportError(f"missing or empty NSET={name}")
    return sorted(set(values))


def read_any_set(deck: Path, aliases: tuple[str, ...]) -> list[int]:
    """First present NSET among ``aliases`` (backward-compatible rename path)."""
    for name in aliases:
        try:
            return read_set_ids(deck, name)
        except ExportError:
            continue
    raise ExportError(f"missing or empty NSET among {aliases}")


def read_coordinates(deck: Path, wanted: set[int]) -> dict[int, tuple[float, float, float]]:
    """Coordinates of the wanted node ids, keyed by node id."""
    values: dict[int, tuple[float, float, float]] = {}
    active = False
    for line in deck.read_text(encoding="ascii").splitlines():
        if line.upper().startswith("*NODE"):
            active = True
            continue
        if line.startswith("*"):
            active = False
            continue
        if active:
            columns = line.split(",")
            try:
                node_id = int(columns[0])
            except ValueError:
                continue
            if node_id in wanted and len(columns) >= 4:
                values[node_id] = (float(columns[1]), float(columns[2]),
                                   float(columns[3]))
    return values


def dat_frequencies(path: Path, modes: int) -> list[float]:
    """Parse the eigenvalue table; require exactly ``modes`` frequencies."""
    values: list[float] = []
    in_table = False
    for line in path.read_bytes().decode("ascii", errors="replace").splitlines():
        if "MODE NO" in line and "FREQUENCY" in line:
            in_table = True
            continue
        if not in_table:
            continue
        columns = line.split()
        if len(columns) < 4 or not columns[0].isdigit():
            if values:
                break
            continue
        values.append(float(columns[3]))
        if len(values) == modes:
            return values
    raise GateError(f"expected {modes} frequencies in DAT result, found {len(values)}")


def gate_frequencies(dat: Path, modes: int, window: tuple[float, float]) -> list[float]:
    """Gate 1+2: mode count, positivity, strict ascent, sanity window."""
    if modes < 1:
        raise GateError("modes must be positive")
    if not dat.is_file():
        raise GateError(f"missing CalculiX DAT result: {dat}")
    frequencies = dat_frequencies(dat, modes)
    if len(frequencies) != modes:
        raise GateError(f"expected exactly {modes} frequencies, found {len(frequencies)}")
    if any(value <= 0.0 for value in frequencies):
        raise GateError(f"non-positive frequency in {dat}: {frequencies}")
    for earlier, later in zip(frequencies, frequencies[1:]):
        if later <= earlier:
            raise GateError(f"frequencies are not strictly ascending: {frequencies}")
    low, high = window
    for index, value in enumerate(frequencies):
        if not low <= value <= high:
            raise GateError(f"mode {index + 1} frequency {value:.6g} Hz outside the "
                            f"sanity window [{low:.6g}, {high:.6g}] Hz")
    return frequencies


def _fields(line: str) -> tuple[int, tuple[float, float, float]]:
    """FRD fixed-width record: node id at 3:13, components at 13:25:37:49."""
    return int(line[3:13]), (float(line[13:25]), float(line[25:37]),
                             float(line[37:49]))


def _write_group(output: Path, name: str, ids: list[int],
                 data: dict[int, tuple[float, float, float]],
                 frequency: float | None, modes: int) -> None:
    """One FDM CSV: three-field header with trailing commas, then 3 values/node."""
    with (output / name).open("w", encoding="ascii") as stream:
        header_frequency = 0.0 if frequency is None else frequency
        stream.write(f"{header_frequency:20.10E},{len(ids):20.10E},{modes:20.10E},\n")
        for node_id in ids:
            vector = data[node_id]
            stream.write(f"{vector[0]:20.10E},{vector[1]:20.10E},{vector[2]:20.10E},\n")


def export_modes(deck: Path, dat: Path, frd: Path, output: Path,
                 modes: int = 10) -> list[float]:
    """Export fluid/struct coordinates and modal displacements (unchanged contract)."""
    if modes < 1:
        raise ExportError("modes must be positive")
    if not dat.is_file():
        raise ExportError(f"missing CalculiX DAT result: {dat}")
    dat_values = dat_frequencies(dat, modes)
    fluid_ids = read_any_set(deck, FLUID_SET_ALIASES)
    struct_ids = read_any_set(deck, STRUCT_SET_ALIASES)
    wanted = set(fluid_ids + struct_ids)
    output.mkdir(parents=True, exist_ok=True)
    coordinates = read_coordinates(deck, wanted)
    phase, active, current_frequency = "none", False, 0.0
    frequencies: list[float] = []
    displacements: dict[int, tuple[float, float, float]] = {}
    with frd.open(encoding="ascii") as source:
        for line in source:
            record = line.lstrip()
            if record.startswith("2C"):
                phase, active = "coordinates", False
            elif record.startswith("100CL"):
                if phase == "coordinates":
                    _write_group(output, "FluidNodeCoor.csv", fluid_ids,
                                 coordinates, None, modes)
                    _write_group(output, "StructNodeCoor.csv", struct_ids,
                                 coordinates, None, modes)
                phase, active = "displacements", False
                match = FREQUENCY.search(record)
                if match is None:
                    raise ExportError("FRD modal frequency header is malformed")
                current_frequency = float(match.group(1))
                frequencies.append(current_frequency)
            elif record.startswith("-4  DISP"):
                active = True
            elif line.startswith(" -3"):
                if phase == "displacements" and active:
                    if len(displacements) != len(wanted):
                        raise ExportError("FRD modal displacement data is incomplete")
                    mode = len(frequencies)
                    _write_group(output, f"FluidNodeDisp{mode}.csv", fluid_ids,
                                 displacements, dat_values[mode - 1], modes)
                    _write_group(output, f"StructNodeDisp{mode}.csv", struct_ids,
                                 displacements, dat_values[mode - 1], modes)
                    displacements = {}
                active = False
            elif active and line.startswith(" -1"):
                node_id, vector = _fields(line)
                if node_id in wanted and phase == "displacements":
                    displacements[node_id] = vector
    if phase != "displacements" or active or len(frequencies) != modes:
        raise ExportError(f"expected {modes} FRD modes, found {len(frequencies)}")
    (output / "FluidPara.csv").write_text(
        "19,\n11,\n0.0,\n0.0,\n0.0,\n0.0,\n0.0,\n", encoding="ascii")
    return dat_values


def parse_frd_modes(path: Path) -> list[ModalBlock]:
    """Parse every ``-4 DISP`` block of an FRD file into node ids + vectors."""
    if not path.is_file():
        raise GateError(f"missing FRD result: {path}")
    modes: list[ModalBlock] = []
    with path.open("rb") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith(b" -4  DISP"):
                continue
            records: list[bytes] = []
            while True:
                position = handle.tell()
                record = handle.readline()
                if not record:
                    break
                if record.startswith(b" -1"):
                    records.append(record)
                elif record.startswith(b" -5"):
                    continue
                elif record.startswith(b" -3"):
                    break
                else:
                    handle.seek(position)
                    break
            if not records:
                continue
            width = len(records[0])
            if width not in (49, 50) or any(len(item) != width for item in records):
                raise GateError(f"non-uniform FRD record width in {path}")
            if width == 49:
                records = [item + b"\n" for item in records]
            raw = np.frombuffer(b"".join(records), dtype=np.uint8).reshape(-1, 50)
            nodes = raw[:, 3:13].view("S10").astype(np.int64).reshape(-1)
            components = [raw[:, 13 + 12 * k:25 + 12 * k].view("S12")
                          .astype(np.float64).reshape(-1) for k in range(3)]
            modes.append(ModalBlock(nodes=nodes,
                                    vectors=np.stack(components, axis=1)))
    return modes


def combined_wing_field(block: ModalBlock, struct_lookup: npt.NDArray[np.bool_],
                        struct_count: int, label: str) -> FloatArray:
    """Combined FRD block restricted to wing nodes, ascending node id order."""
    mask = struct_lookup[block.nodes]
    nodes = block.nodes[mask]
    if nodes.shape[0] != struct_count:
        raise GateError(f"{label}: {nodes.shape[0]} wing nodes != {struct_count}")
    order = np.argsort(nodes, kind="stable")
    return block.vectors[mask][order]


def pure_wing_field(block: ModalBlock, struct_count: int,
                    label: str) -> FloatArray:
    """Pure wing FRD block in ascending node id order.

    The merged deck renumbers the wing by ascending original id, so the k-th
    smallest combined wing id is the k-th smallest pure wing node: sorting both
    blocks by node id puts them in the same modal DOF order.
    """
    nodes = block.nodes
    if nodes.shape[0] != struct_count:
        raise GateError(f"{label}: {nodes.shape[0]} wing nodes != {struct_count}")
    order = np.argsort(nodes, kind="stable")
    return block.vectors[order]


def identity_gate(struct_ids: list[int], combined_frd: Path, pure_frd: Path,
                  modes: int, min_cosine: float) -> FloatArray:
    """Gate 4: combined modes must be wing modes of the pure-wing reference."""
    struct_count = len(struct_ids)
    struct_lookup = np.zeros(max(struct_ids) + 1, dtype=bool)
    struct_lookup[np.asarray(struct_ids, dtype=np.int64)] = True
    pure_blocks = parse_frd_modes(pure_frd)
    combined_blocks = parse_frd_modes(combined_frd)
    if len(pure_blocks) != modes:
        raise GateError(f"pure FRD holds {len(pure_blocks)} modes, expected {modes}")
    if len(combined_blocks) != modes:
        raise GateError(f"combined FRD holds {len(combined_blocks)} modes, "
                        f"expected {modes}")
    pure_fields = np.stack(
        [pure_wing_field(block, struct_count, pure_frd.name)
         for block in pure_blocks[:modes]])
    combined_fields = np.stack(
        [combined_wing_field(block, struct_lookup, struct_count, combined_frd.name)
         for block in combined_blocks[:modes]])
    flat_pure = pure_fields.reshape(modes, -1)
    flat_combined = combined_fields.reshape(modes, -1)
    norm_pure = np.linalg.norm(flat_pure, axis=1, keepdims=True)
    norm_combined = np.linalg.norm(flat_combined, axis=1, keepdims=True)
    if np.any(norm_pure == 0.0) or np.any(norm_combined == 0.0):
        raise GateError("zero displacement norm in an FRD mode")
    cosine = (flat_combined @ flat_pure.T) / (norm_combined * norm_pure.T)
    best = np.argmax(np.abs(cosine), axis=1)
    best_cosine = np.abs(cosine)[np.arange(modes), best]
    for index in range(modes):
        if int(best[index]) != index:
            raise GateError(f"mode {index + 1} best-matches pure mode "
                            f"{int(best[index]) + 1} (expected diagonal)")
        if float(best_cosine[index]) < min_cosine:
            raise GateError(f"mode {index + 1} field cosine "
                            f"{float(best_cosine[index]):.6f} < {min_cosine}")
    print("[AGARD 03] mode identity (cosine vs best pure-wing mode)")
    print("  mode  best_pure  cosine")
    for index in range(modes):
        print(f"  {index + 1:4d}  {int(best[index]) + 1:9d}  "
              f"{float(best_cosine[index]):.6f}")
    return cosine


def frequency_gate(combined: list[float], pure: list[float],
                   max_rel: float) -> float:
    """Gate 3: print the per-mode table and return the max |relative| delta."""
    print("[AGARD 03] mode  f_combined(Hz)    f_pure(Hz)       delta(Hz)    rel")
    relative: list[float] = []
    for index, (combined_value, pure_value) in enumerate(zip(combined, pure,
                                                             strict=True)):
        delta = combined_value - pure_value
        ratio = delta / pure_value
        relative.append(ratio)
        print(f"  {index + 1:4d}  {combined_value:14.6f}  {pure_value:14.6f}  "
              f"{delta:+13.6f}  {ratio * 100:+.4f} %")
    max_rel_value = max(abs(value) for value in relative)
    status = "PASS" if max_rel_value <= max_rel else "FAIL"
    print(f"[AGARD 03] max |rel delta| = {max_rel_value * 100:.4f} % "
          f"(limit {max_rel * 100:.4f} %) -> {status}")
    if max_rel_value > max_rel:
        raise GateError(f"combined-vs-pure max |relative difference| "
                        f"{max_rel_value:.6g} > {max_rel}")
    return max_rel_value


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="combined modal deck (node sets + coordinates)")
    parser.add_argument("--dat", type=Path, required=True,
                        help="combined CalculiX DAT result")
    parser.add_argument("--frd", type=Path, required=True,
                        help="combined CalculiX FRD result")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="FDM CSV directory")
    parser.add_argument("--modes", type=int, default=10)
    parser.add_argument("--pure-dat", type=Path, default=None,
                        help="pure-wing DAT reference for the frequency gate")
    parser.add_argument("--pure-frd", type=Path, default=None,
                        help="pure-wing FRD reference for the mode-identity gate")
    parser.add_argument("--max-rel", type=float, default=DEFAULT_MAX_REL,
                        help="max |relative frequency difference| (default 0.01)")
    parser.add_argument("--min-cosine", type=float, default=DEFAULT_MIN_COSINE,
                        help="min field cosine against the pure wing (default 0.9)")
    parser.add_argument("--min-frequency", type=float, default=DEFAULT_MIN_FREQUENCY,
                        help="lower sanity-window bound in Hz (default 1)")
    parser.add_argument("--max-frequency", type=float, default=DEFAULT_MAX_FREQUENCY,
                        help="upper sanity-window bound in Hz (default 1e4)")
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    if args.modes < 1:
        parser.error("--modes must be positive")
    if args.max_rel <= 0:
        parser.error("--max-rel must be positive")
    if not 0.0 < args.min_cosine <= 1.0:
        parser.error("--min-cosine must be in (0, 1]")
    if args.min_frequency < 0 or args.max_frequency <= args.min_frequency:
        parser.error("--min-frequency/--max-frequency sanity window is invalid")
    try:
        frequencies = gate_frequencies(args.dat, args.modes,
                                       (args.min_frequency, args.max_frequency))
        if not args.pure_dat and not args.pure_frd:
            print("[AGARD 03] no pure-wing reference supplied; frequency and "
                  "identity gates skipped")
            print(f"[AGARD 03] combined frequencies: {frequencies}")
        if args.pure_dat is not None:
            pure = gate_frequencies(args.pure_dat, args.modes,
                                    (args.min_frequency, args.max_frequency))
            frequency_gate(frequencies, pure, args.max_rel)
        if args.pure_frd is not None:
            struct_ids = read_any_set(args.input, STRUCT_SET_ALIASES)
            identity_gate(struct_ids, args.frd, args.pure_frd, args.modes,
                          args.min_cosine)
        exported = export_modes(args.input, args.dat, args.frd,
                                args.output_dir, args.modes)
    except ExportError as error:
        print(f"[AGARD ERR] {error}", file=sys.stderr)
        print("[AGARD ERR] verification failed; no CSV was written", file=sys.stderr)
        return 1
    print(f"[AGARD 03] gates passed; exported={len(exported)} "
          f"output={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
