#!/usr/bin/env python3
"""Generate small-cell-focused pseudo-fluid material layers.

Run from /root/Workspace/tmp:
    python3 /root/OpenFOAM/user-v2412/src/case/AGARD/FDM/gen_pseudofluid_sections.py
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import numpy as np
import numpy.typing as npt

from agard_mesh import read_coordinates, read_element_volumes, read_node_set

DEFAULT_INPUT: Final = Path(
    "/root/Workspace/Oldversion/AGARD4456/Mesh/Assembly/Assembly.inp"
)
DEFAULT_OUTPUT: Final = Path(
    "/root/Workspace/AGARD/Mode/fluid_stiffness.inp"
)

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]


@dataclass(frozen=True, slots=True)
class Config:
    input_path: Path
    output_path: Path
    layer_count: int
    maximum_modulus: float
    protected_modulus: float
    exponent: float


def parse_args() -> Config:
    """Parse generation parameters."""
    parser = argparse.ArgumentParser(
        description="Generate volume-layered pseudo-fluid materials."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--layers", type=int, default=20)
    parser.add_argument("--maximum-modulus", type=float, default=8.17e6)
    parser.add_argument("--protected-modulus", type=float, default=8.17e5)
    parser.add_argument("--exponent", type=float, default=2.5)
    arguments = parser.parse_args()
    if arguments.layers < 4:
        parser.error("--layers must be at least 4")
    if not 0 < arguments.protected_modulus <= arguments.maximum_modulus:
        parser.error("protected modulus must be positive and no greater than maximum")
    if arguments.exponent <= 0:
        parser.error("--exponent must be positive")
    return Config(
        input_path=arguments.input,
        output_path=arguments.output,
        layer_count=arguments.layers,
        maximum_modulus=arguments.maximum_modulus,
        protected_modulus=arguments.protected_modulus,
        exponent=arguments.exponent,
    )


def calculate_modulus(
    volume: float,
    minimum_volume: float,
    protected_limit: float,
    config: Config,
) -> float:
    """Map representative volume to protected or decaying stiffness."""
    if volume <= protected_limit:
        span = np.log(protected_limit / minimum_volume)
        score = np.log(protected_limit / volume) / span
        return config.protected_modulus + (
            config.maximum_modulus - config.protected_modulus
        ) * score**config.exponent
    return config.protected_modulus * (protected_limit / volume) ** config.exponent


def build_edges(
    volumes: FloatArray,
    interface_volumes: FloatArray,
    layer_count: int,
) -> tuple[FloatArray, float]:
    """Reserve 75% of log-volume layers for FSI-adjacent small cells."""
    small_layers = round(layer_count * 0.75)
    large_layers = layer_count - small_layers
    protected_limit = float(np.quantile(interface_volumes, 0.95))
    small_edges = np.geomspace(volumes.min(), protected_limit, small_layers + 1)
    large_edges = np.geomspace(protected_limit, volumes.max(), large_layers + 1)
    return np.concatenate((small_edges, large_edges[1:])), protected_limit


def write_layers(
    config: Config,
    element_ids: IntArray,
    volumes: FloatArray,
    interface_volumes: FloatArray,
) -> None:
    """Write 20 ELSET/material/section groups to a CalculiX include."""
    edges, protected_limit = build_edges(
        volumes, interface_volumes, config.layer_count
    )
    layer_ids = np.clip(
        np.digitize(volumes, edges[1:-1]), 0, config.layer_count - 1
    )
    minimum_volume = float(volumes.min())

    with config.output_path.open("w", encoding="ascii") as output:
        output.write("** Volume-layered pseudo-fluid materials\n")
        output.write(
            f"** Vlimit={protected_limit:.9e}; E={config.protected_modulus}.."
            f"{config.maximum_modulus} below Vlimit, then E~(Vlimit/V)^"
            f"{config.exponent}\n"
        )
        for layer_index in range(config.layer_count):
            mask = layer_ids == layer_index
            ids = element_ids[mask]
            layer_volumes = volumes[mask]
            if len(ids) == 0:
                continue
            representative_volume = float(np.exp(np.mean(np.log(layer_volumes))))
            stiffness_volume = (
                representative_volume
                if representative_volume <= protected_limit
                else float(layer_volumes.min())
            )
            modulus = (
                config.maximum_modulus
                if layer_index == 0
                else calculate_modulus(
                    stiffness_volume, minimum_volume, protected_limit, config
                )
            )
            name = f"PFLUID_L{layer_index:02d}"
            output.write(
                f"** layer={layer_index} count={len(ids)} "
                f"V=[{layer_volumes.min():.9e},{layer_volumes.max():.9e}] "
                f"E={modulus:.9e}\n"
            )
            output.write(f"*MATERIAL, NAME={name}\n*ELASTIC\n{modulus:.9e},0.\n")
            output.write("*DENSITY\n0.\n")
            output.write(f"*ELSET, ELSET={name}\n")
            for offset in range(0, len(ids), 16):
                output.write(",".join(str(value) for value in ids[offset : offset + 16]))
                output.write("\n")
            output.write(f"*SOLID SECTION, ELSET={name}, MATERIAL={name}\n")
            print(
                f"layer {layer_index:02d}: count={len(ids)}, "
                f"V={layer_volumes.min():.3e}..{layer_volumes.max():.3e}, "
                f"E={modulus:.3e}"
            )
    print(f"FSI protected volume limit: {protected_limit:.6e}")


def main() -> None:
    """Generate the pseudo-fluid stiffness include."""
    config = parse_args()
    coordinates = read_coordinates(config.input_path)
    interface_nodes = read_node_set(config.input_path, "FSI")
    element_ids, volumes, interface_volumes = read_element_volumes(
        config.input_path, coordinates, "ETYPE1", interface_nodes
    )
    if np.any(volumes <= 0) or len(interface_volumes) == 0:
        raise RuntimeError("invalid fluid or FSI-adjacent volume data")
    print(
        f"fluid elements={len(element_ids)}, FSI-adjacent={len(interface_volumes)}, "
        f"volume={volumes.min():.6e}..{volumes.max():.6e}"
    )
    write_layers(config, element_ids, volumes, interface_volumes)
    print(f"written: {config.output_path}")


if __name__ == "__main__":
    main()
