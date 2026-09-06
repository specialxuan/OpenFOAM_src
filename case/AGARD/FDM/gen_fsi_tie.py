#!/usr/bin/env python3
"""Generate CalculiX element surfaces and the AGARD FSI tie constraint.

Run from /root/Workspace/tmp:
    python3 /root/OpenFOAM/user-v2412/src/case/AGARD/FDM/gen_fsi_tie.py

The fluid-side surface is identified by NSET=FSI. The wing-side surface is
identified by the union of NSET=SIDE and NSET=TIP. ROOT is intentionally
excluded because it is fixed rather than tied.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Final, TextIO

from agard_mesh import read_node_set

DEFAULT_INPUT: Final = Path(
    "/root/Workspace/Oldversion/AGARD4456/Mesh/Assembly/Assembly.inp"
)
DEFAULT_OUTPUT: Final = Path(
    "/root/Workspace/AGARD/Mode/fsi_tie.inp"
)
FACES: Final = (
    ("S1", (0, 1, 2, 3)),
    ("S2", (4, 5, 6, 7)),
    ("S3", (0, 1, 5, 4)),
    ("S4", (1, 2, 6, 5)),
    ("S5", (2, 3, 7, 6)),
    ("S6", (3, 0, 4, 7)),
)
FARFIELD_SETS: Final = ("IN", "OUT", "TOP", "BOTTOM", "FRONT", "BACK")


@dataclass(frozen=True, slots=True)
class Config:
    input_path: Path
    output_path: Path
    tolerance: float


def parse_args() -> Config:
    """Parse generation parameters."""
    parser = argparse.ArgumentParser(description="Generate AGARD FSI tie surfaces.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--tolerance", type=float, default=0.05)
    arguments = parser.parse_args()
    if arguments.tolerance <= 0:
        parser.error("--tolerance must be positive")
    return Config(arguments.input, arguments.output, arguments.tolerance)


def write_surface_faces(
    source_path: Path,
    output: TextIO,
    elset_name: str,
    surface_nodes: frozenset[int],
) -> int:
    """Write C3D8 faces whose four nodes belong to the requested node set."""
    count = 0
    in_elset = False
    with source_path.open(encoding="ascii") as source:
        for line in source:
            if line.startswith("*ELEMENT"):
                in_elset = f"ELSET={elset_name}" in line
                continue
            if line.startswith("*"):
                in_elset = False
                continue
            if not in_elset:
                continue
            fields = line.strip().split(",")
            element_id = fields[0]
            nodes = tuple(int(value) for value in fields[1:9])
            for face_name, indices in FACES:
                if all(nodes[index] in surface_nodes for index in indices):
                    output.write(f"{element_id},{face_name}\n")
                    count += 1
    return count


def write_node_set(output: TextIO, name: str, node_ids: list[int]) -> None:
    output.write(f"*NSET, NSET={name}\n")
    for offset in range(0, len(node_ids), 16):
        output.write(",".join(str(value) for value in node_ids[offset : offset + 16]))
        output.write("\n")


def main() -> None:
    """Generate fluid and wing surfaces followed by a tied contact."""
    config = parse_args()
    fluid_nodes = frozenset(read_node_set(config.input_path, "FSI").tolist())
    wing_nodes = frozenset(
        read_node_set(config.input_path, "SIDE").tolist()
        + read_node_set(config.input_path, "TIP").tolist()
    )

    with config.output_path.open("w", encoding="ascii") as output:
        output.write("** AGARD fluid-wing interface reconstructed from ICEM node sets\n")
        output.write("*SURFACE, NAME=SURF_FSI_FLUID, TYPE=ELEMENT\n")
        fluid_faces = write_surface_faces(
            config.input_path, output, "ETYPE1", fluid_nodes
        )
        output.write("*SURFACE, NAME=SURF_FSI_WING, TYPE=ELEMENT\n")
        wing_faces = write_surface_faces(
            config.input_path, output, "ETYPE2", wing_nodes
        )
        output.write(
            "*TIE, NAME=FSI_COUPLING, "
            f"POSITION TOLERANCE={config.tolerance:.9g}, ADJUST=NO\n"
        )
        output.write("SURF_FSI_FLUID,SURF_FSI_WING\n")
        for set_name in FARFIELD_SETS:
            original_nodes = read_node_set(config.input_path, set_name).tolist()
            fixed_nodes = sorted(set(original_nodes) - fluid_nodes)
            write_node_set(output, f"FIXED_{set_name}", fixed_nodes)

    if fluid_faces == 0 or wing_faces == 0:
        raise RuntimeError("no faces found on one side of the FSI tie")
    print(
        f"fluid nodes={len(fluid_nodes)}, faces={fluid_faces}; "
        f"wing nodes={len(wing_nodes)}, faces={wing_faces}"
    )
    print(f"written: {config.output_path}")


if __name__ == "__main__":
    main()
