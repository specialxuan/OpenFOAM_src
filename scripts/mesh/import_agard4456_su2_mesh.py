#!/usr/bin/env python3
"""Import the public AGARD 445.6 SU2 tetrahedral mesh as OpenFOAM polyMesh."""

from __future__ import annotations

import argparse
import math
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


Point = Tuple[float, float, float]
Face = Tuple[int, int, int]

PATCH_MAP = {
    "1": ("farfield", "patch", None),
    "2": ("rootSymmetry", "symmetryPlane", None),
    "3": ("AGARD_WING", "wall", "wall"),
}


def parse_su2(path: Path) -> Tuple[List[Tuple[int, int, int, int]], List[Point], Dict[str, List[Face]]]:
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        line = f.readline().strip()
        if line != "NDIME= 3":
            raise ValueError(f"{path}: expected 'NDIME= 3', got {line!r}")

        line = f.readline().strip()
        if not line.startswith("NELEM="):
            raise ValueError(f"{path}: missing NELEM")
        n_elem = int(line.split("=", 1)[1])

        cells: List[Tuple[int, int, int, int]] = []
        for _ in range(n_elem):
            parts = f.readline().split()
            if not parts:
                raise ValueError(f"{path}: unexpected EOF while reading elements")
            elem_type = int(parts[0])
            if elem_type != 10:
                raise ValueError(f"{path}: only SU2 tetra type 10 is supported, got {elem_type}")
            cells.append(tuple(int(v) for v in parts[1:5]))  # type: ignore[arg-type]

        line = f.readline().strip()
        if not line.startswith("NPOIN="):
            raise ValueError(f"{path}: missing NPOIN")
        n_points = int(line.split("=", 1)[1])

        points: List[Point] = []
        for _ in range(n_points):
            parts = f.readline().split()
            if len(parts) < 3:
                raise ValueError(f"{path}: malformed point line")
            points.append((float(parts[0]), float(parts[1]), float(parts[2])))

        line = f.readline().strip()
        if not line.startswith("NMARK="):
            raise ValueError(f"{path}: missing NMARK")
        n_mark = int(line.split("=", 1)[1])

        markers: Dict[str, List[Face]] = {}
        for _ in range(n_mark):
            tag_line = f.readline().strip()
            elem_line = f.readline().strip()
            if not tag_line.startswith("MARKER_TAG=") or not elem_line.startswith("MARKER_ELEMS="):
                raise ValueError(f"{path}: malformed marker header")
            tag = tag_line.split("=", 1)[1].strip()
            n_marker_faces = int(elem_line.split("=", 1)[1])
            faces: List[Face] = []
            for _ in range(n_marker_faces):
                parts = f.readline().split()
                face_type = int(parts[0])
                if face_type != 5:
                    raise ValueError(f"{path}: only SU2 triangle marker type 5 is supported, got {face_type}")
                faces.append(tuple(int(v) for v in parts[1:4]))  # type: ignore[arg-type]
            markers[tag] = faces

    return cells, points, markers


def sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Point, b: Point) -> Point:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def center(ids: Sequence[int], points: Sequence[Point]) -> Point:
    inv = 1.0 / len(ids)
    return (
        sum(points[i][0] for i in ids) * inv,
        sum(points[i][1] for i in ids) * inv,
        sum(points[i][2] for i in ids) * inv,
    )


def orient_face_outward(face: Face, cell_center: Point, points: Sequence[Point]) -> Face:
    p0, p1, p2 = (points[i] for i in face)
    normal = cross(sub(p1, p0), sub(p2, p0))
    f_center = center(face, points)
    if dot(normal, sub(f_center, cell_center)) < 0.0:
        return (face[0], face[2], face[1])
    return face


def build_faces(
    cells: Sequence[Tuple[int, int, int, int]],
    points: Sequence[Point],
    markers: Dict[str, List[Face]],
) -> Tuple[List[Face], List[int], List[int], Dict[str, List[Tuple[Face, int]]], Dict[str, int]]:
    face_owner: Dict[Tuple[int, int, int], int] = {}
    face_oriented: Dict[Tuple[int, int, int], Face] = {}
    face_neighbour: Dict[Tuple[int, int, int], int] = {}

    for cell_i, tet in enumerate(cells):
        c_center = center(tet, points)
        raw_faces = [
            (tet[0], tet[1], tet[2]),
            (tet[0], tet[3], tet[1]),
            (tet[0], tet[2], tet[3]),
            (tet[1], tet[3], tet[2]),
        ]
        for raw in raw_faces:
            key = tuple(sorted(raw))
            if key in face_owner:
                face_neighbour[key] = cell_i
            else:
                face_owner[key] = cell_i
                face_oriented[key] = orient_face_outward(raw, c_center, points)

    marker_by_key: Dict[Tuple[int, int, int], str] = {}
    marker_counts: Dict[str, int] = {}
    for tag, faces in markers.items():
        marker_counts[tag] = len(faces)
        for face in faces:
            marker_by_key[tuple(sorted(face))] = tag

    internal_keys = sorted(
        (key for key in face_owner if key in face_neighbour),
        key=lambda key: (face_owner[key], face_neighbour[key]),
    )
    boundary_keys = [key for key in face_owner if key not in face_neighbour]

    patch_keys: Dict[str, List[Tuple[Face, int]]] = defaultdict(list)
    for key in boundary_keys:
        tag = marker_by_key.get(key, "unassigned")
        patch_keys[tag].append((face_oriented[key], face_owner[key]))

    ordered_faces: List[Face] = []
    owners: List[int] = []
    neighbours: List[int] = []

    for key in internal_keys:
        ordered_faces.append(face_oriented[key])
        owners.append(face_owner[key])
        neighbours.append(face_neighbour[key])

    return ordered_faces, owners, neighbours, patch_keys, marker_counts


def write_header(fp, cls: str, obj: str) -> None:
    fp.write("FoamFile\n")
    fp.write("{\n")
    fp.write("    version     2.0;\n")
    fp.write("    format      ascii;\n")
    fp.write(f"    class       {cls};\n")
    fp.write('    location    "constant/polyMesh";\n')
    fp.write(f"    object      {obj};\n")
    fp.write("}\n")
    fp.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")


def write_points(path: Path, points: Sequence[Point]) -> None:
    with path.open("w", encoding="utf-8") as fp:
        write_header(fp, "vectorField", "points")
        fp.write(f"{len(points)}\n(\n")
        for x, y, z in points:
            fp.write(f"({x:.12g} {y:.12g} {z:.12g})\n")
        fp.write(")\n")


def write_faces(path: Path, faces: Sequence[Face]) -> None:
    with path.open("w", encoding="utf-8") as fp:
        write_header(fp, "faceList", "faces")
        fp.write(f"{len(faces)}\n(\n")
        for face in faces:
            fp.write(f"3({face[0]} {face[1]} {face[2]})\n")
        fp.write(")\n")


def write_label_list(path: Path, obj: str, values: Iterable[int]) -> None:
    vals = list(values)
    with path.open("w", encoding="utf-8") as fp:
        write_header(fp, "labelList", obj)
        fp.write(f"{len(vals)}\n(\n")
        for v in vals:
            fp.write(f"{v}\n")
        fp.write(")\n")


def write_boundary(path: Path, patch_entries: Sequence[Tuple[str, str, str | None, int, int]]) -> None:
    with path.open("w", encoding="utf-8") as fp:
        write_header(fp, "polyBoundaryMesh", "boundary")
        fp.write(f"{len(patch_entries)}\n(\n")
        for name, patch_type, group, n_faces, start_face in patch_entries:
            fp.write(f"    {name}\n")
            fp.write("    {\n")
            fp.write(f"        type            {patch_type};\n")
            if group:
                fp.write(f"        inGroups        1({group});\n")
            fp.write(f"        nFaces          {n_faces};\n")
            fp.write(f"        startFace       {start_face};\n")
            fp.write("    }\n")
        fp.write(")\n")


def write_control_dict(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
application     myRhoPimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.005;
deltaT          1e-5;
writeControl    runTime;
writeInterval   0.001;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
""",
        encoding="utf-8",
    )


def write_fv_schemes(path: Path) -> None:
    path.write_text(
        """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}

ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}
""",
        encoding="utf-8",
    )


def write_fv_solution(path: Path) -> None:
    path.write_text(
        """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
}

PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors 1;
    nNonOrthogonalCorrectors 0;
}
""",
        encoding="utf-8",
    )


def write_source_notes(
    path: Path,
    su2_source: Path,
    bdf_source: Path | None,
    dat_source: Path | None,
    csm_source: Path | None,
    marker_counts: Dict[str, int],
    patch_counts: Sequence[Tuple[str, int]],
) -> None:
    lines = [
        "# AGARD 445.6 Imported Mesh Notes",
        "",
        "This case starts from the public AGARD445 resources in CODE-Lab-IASTATE/CRM_Flutter_Prediction_Framework.",
        "",
        "Public source:",
        "- https://github.com/CODE-Lab-IASTATE/CRM_Flutter_Prediction_Framework/tree/main/AGARD445",
        "",
        "Copied source files are stored in `reference_public_agard445/` so the case does not depend on `/tmp`.",
        "",
        "Imported flow mesh:",
        f"- Source SU2 mesh: `{su2_source}`",
        "- Converted to OpenFOAM `constant/polyMesh` with `scripts/import_agard4456_su2_mesh.py`.",
        "- SU2 marker mapping:",
        "  - `3` -> `AGARD_WING` (`wall`)",
        "  - `2` -> `rootSymmetry` (`symmetryPlane`)",
        "  - `1` -> `farfield` (`patch`)",
        "",
        "Source structural/material references copied for the later CalculiX step:",
    ]
    for label, src in (("MYSTRAN BDF", bdf_source), ("MYSTRAN DAT", dat_source), ("CAPS CSM", csm_source)):
        if src:
            lines.append(f"- {label}: `{src}`")
    lines.extend(
        [
            "",
            "Known material/property values from the public MYSTRAN deck:",
            "- `MAT1`: E = `5.34E+8`, G = `8.31E+7`, nu = `0.31`, rho = `147`.",
            "- `PSHELL` skin thickness = `1.65E-2`.",
            "- `PSHELL` root/rib-spar thickness = `5.00E-5`.",
            "",
            "Marker face counts in the source mesh:",
        ]
    )
    for tag in sorted(marker_counts):
        lines.append(f"- marker `{tag}`: {marker_counts[tag]}")
    lines.extend(["", "OpenFOAM patch face counts:"])
    for name, count in patch_counts:
        lines.append(f"- `{name}`: {count}")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def copy_reference_files(out_dir: Path, sources: Sequence[Path]) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    copied: List[Path] = []
    for src in sources:
        if src and src.exists():
            dst = out_dir / src.name
            shutil.copy2(src, dst)
            copied.append(dst)
    return copied


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--su2", required=True, help="Input AGARD 445.6 SU2 mesh")
    parser.add_argument("--case-dir", required=True, help="Output OpenFOAM case directory")
    parser.add_argument("--bdf", help="Optional MYSTRAN BDF source to copy")
    parser.add_argument("--dat", help="Optional MYSTRAN DAT source to copy")
    parser.add_argument("--csm", help="Optional CAPS CSM source to copy")
    args = parser.parse_args()

    su2 = Path(args.su2).resolve()
    case_dir = Path(args.case_dir).resolve()
    poly = case_dir / "constant" / "polyMesh"
    poly.mkdir(parents=True, exist_ok=True)

    cells, points, markers = parse_su2(su2)
    internal_faces, owners, neighbours, patch_keys, marker_counts = build_faces(cells, points, markers)

    faces = list(internal_faces)
    owner_values = list(owners)
    patch_entries = []
    patch_counts = []
    start_face = len(faces)

    for tag in ("3", "2", "1", "unassigned"):
        entries = patch_keys.get(tag, [])
        if not entries:
            continue
        if tag in PATCH_MAP:
            name, patch_type, group = PATCH_MAP[tag]
        else:
            name, patch_type, group = "unassignedBoundary", "patch", None
        patch_entries.append((name, patch_type, group, len(entries), start_face))
        patch_counts.append((name, len(entries)))
        for face, owner in entries:
            faces.append(face)
            owner_values.append(owner)
        start_face += len(entries)

    write_points(poly / "points", points)
    write_faces(poly / "faces", faces)
    write_label_list(poly / "owner", "owner", owner_values)
    write_label_list(poly / "neighbour", "neighbour", neighbours)
    write_boundary(poly / "boundary", patch_entries)
    system_dir = case_dir / "system"
    write_control_dict(system_dir / "controlDict")
    write_fv_schemes(system_dir / "fvSchemes")
    write_fv_solution(system_dir / "fvSolution")

    ref_dir = case_dir / "reference_public_agard445"
    bdf = Path(args.bdf).resolve() if args.bdf else None
    dat = Path(args.dat).resolve() if args.dat else None
    csm = Path(args.csm).resolve() if args.csm else None
    copy_reference_files(ref_dir, [su2, *(p for p in (bdf, dat, csm) if p)])
    write_source_notes(case_dir / "README_agard4456_imported_mesh.md", su2, bdf, dat, csm, marker_counts, patch_counts)

    print(f"[OK] case: {case_dir}")
    print(f"[OK] points: {len(points)}")
    print(f"[OK] tetra cells: {len(cells)}")
    print(f"[OK] internal faces: {len(neighbours)}")
    print(f"[OK] total faces: {len(faces)}")
    for name, count in patch_counts:
        print(f"[OK] patch {name}: {count}")


if __name__ == "__main__":
    main()
