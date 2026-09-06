"""Read the converted AGARD CalculiX mesh and calculate C3D8 volumes."""

from __future__ import annotations

from pathlib import Path
from typing import Final

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]

NODES_PER_ELEMENT: Final = 8
ELEMENT_CHUNK_SIZE: Final = 200_000
NATURAL_COORDINATES: Final[FloatArray] = np.array(
    [
        [-1, -1, -1],
        [1, -1, -1],
        [1, 1, -1],
        [-1, 1, -1],
        [-1, -1, 1],
        [1, -1, 1],
        [1, 1, 1],
        [-1, 1, 1],
    ],
    dtype=np.float64,
)


def shape_derivatives(point: FloatArray) -> FloatArray:
    """Return C3D8 shape-function derivatives at one natural point."""
    xi, eta, zeta = point
    result = np.empty((NODES_PER_ELEMENT, 3), dtype=np.float64)
    for index, (a, b, c) in enumerate(NATURAL_COORDINATES):
        result[index, 0] = 0.125 * a * (1 + b * eta) * (1 + c * zeta)
        result[index, 1] = 0.125 * b * (1 + a * xi) * (1 + c * zeta)
        result[index, 2] = 0.125 * c * (1 + a * xi) * (1 + b * eta)
    return result


DERIVATIVES: Final[tuple[FloatArray, ...]] = tuple(
    shape_derivatives(point) for point in NATURAL_COORDINATES / np.sqrt(3.0)
)


def read_coordinates(path: Path) -> FloatArray:
    """Read the CalculiX node block into a node-id-indexed dense array."""
    maximum_node_id = 0
    with path.open(encoding="ascii") as source:
        in_nodes = False
        for line in source:
            if line.startswith("*NODE"):
                in_nodes = True
                continue
            if line.startswith("*"):
                if in_nodes:
                    break
                continue
            if in_nodes:
                maximum_node_id = max(maximum_node_id, int(line.split(",", 1)[0]))

    coordinates = np.zeros((maximum_node_id + 1, 3), dtype=np.float64)
    with path.open(encoding="ascii") as source:
        in_nodes = False
        for line in source:
            if line.startswith("*NODE"):
                in_nodes = True
                continue
            if line.startswith("*"):
                if in_nodes:
                    break
                continue
            if in_nodes:
                fields = line.split(",")
                coordinates[int(fields[0])] = tuple(float(value) for value in fields[1:4])
    return coordinates


def read_node_set(path: Path, set_name: str) -> IntArray:
    """Read one named CalculiX node set."""
    node_ids: list[int] = []
    with path.open(encoding="ascii") as source:
        in_set = False
        for line in source:
            if line.startswith("*NSET"):
                in_set = f"NSET={set_name}" in line
                continue
            if line.startswith("*"):
                in_set = False
                continue
            if in_set:
                node_ids.extend(int(value) for value in line.strip().split(",") if value)
    return np.asarray(node_ids, dtype=np.int64)


def calculate_volumes(coordinates: FloatArray, connectivity: IntArray) -> FloatArray:
    """Integrate positive C3D8 Jacobian determinants at 2x2x2 Gauss points."""
    element_coordinates = coordinates[connectivity]
    volumes = np.zeros(connectivity.shape[0], dtype=np.float64)
    for derivatives in DERIVATIVES:
        jacobians = np.einsum("ia,mij->maj", derivatives, element_coordinates)
        volumes += np.abs(np.linalg.det(jacobians))
    return volumes


def read_element_volumes(
    path: Path,
    coordinates: FloatArray,
    elset_name: str,
    interface_nodes: IntArray,
) -> tuple[IntArray, FloatArray, FloatArray]:
    """Calculate element and interface-adjacent volumes in bounded chunks."""
    id_parts: list[IntArray] = []
    volume_parts: list[FloatArray] = []
    interface_parts: list[FloatArray] = []
    element_ids: list[int] = []
    connectivity: list[list[int]] = []

    def flush() -> None:
        if not element_ids:
            return
        connections = np.asarray(connectivity, dtype=np.int64)
        chunk_volumes = calculate_volumes(coordinates, connections)
        interface_mask = np.count_nonzero(np.isin(connections, interface_nodes), axis=1) >= 4
        id_parts.append(np.asarray(element_ids, dtype=np.int64))
        volume_parts.append(chunk_volumes)
        interface_parts.append(chunk_volumes[interface_mask])
        element_ids.clear()
        connectivity.clear()

    with path.open(encoding="ascii") as source:
        in_target = False
        for line in source:
            if line.startswith("*ELEMENT"):
                in_target = f"ELSET={elset_name}" in line
                continue
            if line.startswith("*"):
                in_target = False
                continue
            if in_target:
                fields = line.split(",")
                element_ids.append(int(fields[0]))
                connectivity.append([int(value) for value in fields[1:9]])
                if len(element_ids) == ELEMENT_CHUNK_SIZE:
                    flush()
    flush()
    return (
        np.concatenate(id_parts),
        np.concatenate(volume_parts),
        np.concatenate(interface_parts),
    )
