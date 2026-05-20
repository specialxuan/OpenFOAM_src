#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BASELINE_ZIP="/root/Workspace/case_damfailure.zip"
WORK_ROOT="${SCRIPT_DIR}"
MESH_LEVEL="all"

LX="0.584"
LY="0.365"
LZ="0.012"
WATER_X="0.146"
WATER_Y="0.292"
FSI_X="0.376"

usage() {
    cat <<'EOF'
Usage:
  generate_blockMesh_cases.sh [options]
  generate_blockMesh_cases.sh [baseline_zip] [work_root]   # legacy positional form

Options:
  --baseline-zip <path>       Baseline zip (default: /root/Workspace/case_damfailure.zip)
  --work-root <path>          Output root (default: script directory)
  --mesh-level <level>        coarse | medium | fine | all (default: all)
  -h, --help                  Show help

Notes:
1. blockMesh is split using physical water limits (x=0.146, y=0.292) so initial water is consistent across mesh levels.
2. FRONT/BACK are generated as symmetryPlane patches.
EOF
}

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --baseline-zip)
            BASELINE_ZIP="$2"
            shift 2
            ;;
        --work-root)
            WORK_ROOT="$2"
            shift 2
            ;;
        --mesh-level)
            MESH_LEVEL="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

if [[ ${#POSITIONAL[@]} -ge 1 ]]; then
    BASELINE_ZIP="${POSITIONAL[0]}"
fi
if [[ ${#POSITIONAL[@]} -ge 2 ]]; then
    WORK_ROOT="${POSITIONAL[1]}"
fi
if [[ ${#POSITIONAL[@]} -gt 2 ]]; then
    echo "Error: too many positional arguments." >&2
    usage
    exit 1
fi

case "${MESH_LEVEL}" in
    coarse|medium|fine|all) ;;
    *)
        echo "Error: --mesh-level must be coarse|medium|fine|all, got: ${MESH_LEVEL}" >&2
        exit 1
        ;;
esac

if [[ ! -f "${BASELINE_ZIP}" ]]; then
    echo "Error: baseline zip not found: ${BASELINE_ZIP}" >&2
    exit 1
fi

TEMPLATE_DIR="${WORK_ROOT}/template_case"

extract_template() {
    rm -rf "${TEMPLATE_DIR}"
    mkdir -p "${TEMPLATE_DIR}"

    python3 - "${BASELINE_ZIP}" "${TEMPLATE_DIR}" <<'PY'
import sys
import zipfile
from pathlib import Path

zip_path = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
prefixes = ("0/", "constant/", "system/", "mode/")

with zipfile.ZipFile(zip_path, "r") as zf:
    for info in zf.infolist():
        name = info.filename
        rel = name
        if not rel.startswith(prefixes) and "/" in rel:
            rel = rel.split("/", 1)[1]
        if not rel.startswith(prefixes):
            continue
        if rel.startswith("constant/polyMesh/"):
            # Do not copy old mesh; regenerate with blockMesh.
            continue
        if not rel or rel.endswith("/"):
            (out_dir / rel).mkdir(parents=True, exist_ok=True)
            continue
        target = out_dir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        with zf.open(info) as src, target.open("wb") as dst:
            dst.write(src.read())
PY
}

write_alpha0() {
    local case_dir="$1"
    cat > "${case_dir}/0/alpha.water" <<'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2412                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alpha.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;

boundaryField
{
    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
EOF
}

write_set_fields_dict() {
    local case_dir="$1"
    cat > "${case_dir}/system/setFieldsDict" <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2412                                  |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      setFieldsDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defaultFieldValues
(
    volScalarFieldValue alpha.water 0
);

regions
(
    boxToCell
    {
        box (0 0 0) (${WATER_X} ${WATER_Y} ${LZ});
        fieldValues
        (
            volScalarFieldValue alpha.water 1
        );
    }
);
EOF
}

patch_front_back_symmetry_fields() {
    local case_dir="$1"
    python3 - "${case_dir}" <<'PY'
import re
import sys
from pathlib import Path

case_dir = Path(sys.argv[1])

replacement_map = {
    "U": ('"(FRONT|TOP|RIGHT|LEFT|BACK|BOTTOM_FLUID)"', '"(TOP|RIGHT|LEFT|BOTTOM_FLUID)"'),
    "p_rgh": ('"(FRONT|TOP|RIGHT|LEFT|BACK|BOTTOM_FLUID|FSI_FLUID)"', '"(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"'),
    "p": ('"(FRONT|TOP|RIGHT|LEFT|BACK|BOTTOM_FLUID|FSI_FLUID)"', '"(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"'),
    "k": ('"(FRONT|TOP|RIGHT|LEFT|BACK|BOTTOM_FLUID|FSI_FLUID)"', '"(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"'),
    "omega": ('"(FRONT|TOP|RIGHT|LEFT|BACK|BOTTOM_FLUID|FSI_FLUID)"', '"(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"'),
    "nut": ('"(FRONT|TOP|RIGHT|LEFT|BACK|BOTTOM_FLUID|FSI_FLUID)"', '"(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"'),
}

inject_block = """boundaryField
{
    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
"""

for rel, (old, new) in replacement_map.items():
    path = case_dir / "0" / rel
    text = path.read_text(encoding="utf-8")
    if old in text:
        text = text.replace(old, new, 1)
    elif not re.search(r"^\s*FRONT\s*\{", text, flags=re.MULTILINE):
        raise RuntimeError(f"{path}: cannot find expected patch group {old}")

    if not re.search(r"^\s*FRONT\s*\{", text, flags=re.MULTILINE):
        marker = "boundaryField\n{"
        if marker not in text:
            raise RuntimeError(f"{path}: cannot find boundaryField block")
        text = text.replace(marker, inject_block, 1)

    path.write_text(text, encoding="utf-8")
PY
}

write_block_mesh_dict() {
    local case_dir="$1"
    local nx="$2"
    local ny="$3"
    local nz="$4"
    local nx_a nx_b nx_c ny_a ny_b

    read -r nx_a nx_b nx_c ny_a ny_b <<< "$(python3 - "${nx}" "${ny}" "${LX}" "${WATER_X}" "${FSI_X}" "${LY}" "${WATER_Y}" <<'PY'
import sys

nx = int(sys.argv[1])
ny = int(sys.argv[2])
lx = float(sys.argv[3])
xw = float(sys.argv[4])
xf = float(sys.argv[5])
ly = float(sys.argv[6])
yw = float(sys.argv[7])

if not (0.0 < xw < xf < lx):
    raise SystemExit("Invalid x splits: require 0 < WATER_X < FSI_X < LX")
if not (0.0 < yw < ly):
    raise SystemExit("Invalid y split: require 0 < WATER_Y < LY")

def split_count(total, seg_lens):
    vals = [max(1, int(round(total * l / sum(seg_lens)))) for l in seg_lens]
    while sum(vals) < total:
        idx = max(range(len(vals)), key=lambda i: seg_lens[i] / vals[i])
        vals[idx] += 1
    while sum(vals) > total:
        candidates = [i for i, v in enumerate(vals) if v > 1]
        if not candidates:
            break
        idx = max(candidates, key=lambda i: vals[i])
        vals[idx] -= 1
    if sum(vals) != total or any(v < 1 for v in vals):
        raise SystemExit(f"Failed to split counts: total={total}, vals={vals}")
    return vals

sx = split_count(nx, [xw, xf - xw, lx - xf])
sy = split_count(ny, [yw, ly - yw])
print(sx[0], sx[1], sx[2], sy[0], sy[1])
PY
)"

    cat > "${case_dir}/system/blockMeshDict" <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2412                                  |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    (0         0         0)
    (${WATER_X} 0         0)
    (${FSI_X}   0         0)
    (${LX}      0         0)
    (0         ${WATER_Y} 0)
    (${WATER_X} ${WATER_Y} 0)
    (${FSI_X}   ${WATER_Y} 0)
    (${LX}      ${WATER_Y} 0)
    (0         ${LY}      0)
    (${WATER_X} ${LY}      0)
    (${FSI_X}   ${LY}      0)
    (${LX}      ${LY}      0)

    (0         0         ${LZ})
    (${WATER_X} 0         ${LZ})
    (${FSI_X}   0         ${LZ})
    (${LX}      0         ${LZ})
    (0         ${WATER_Y} ${LZ})
    (${WATER_X} ${WATER_Y} ${LZ})
    (${FSI_X}   ${WATER_Y} ${LZ})
    (${LX}      ${WATER_Y} ${LZ})
    (0         ${LY}      ${LZ})
    (${WATER_X} ${LY}      ${LZ})
    (${FSI_X}   ${LY}      ${LZ})
    (${LX}      ${LY}      ${LZ})
);

blocks
(
    hex (0 1 5 4 12 13 17 16) (${nx_a} ${ny_a} ${nz}) simpleGrading (1 1 1)
    hex (1 2 6 5 13 14 18 17) (${nx_b} ${ny_a} ${nz}) simpleGrading (1 1 1)
    hex (2 3 7 6 14 15 19 18) (${nx_c} ${ny_a} ${nz}) simpleGrading (1 1 1)
    hex (4 5 9 8 16 17 21 20) (${nx_a} ${ny_b} ${nz}) simpleGrading (1 1 1)
    hex (5 6 10 9 17 18 22 21) (${nx_b} ${ny_b} ${nz}) simpleGrading (1 1 1)
    hex (6 7 11 10 18 19 23 22) (${nx_c} ${ny_b} ${nz}) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    FRONT
    {
        type symmetryPlane;
        faces
        (
            (0 1 5 4)
            (1 2 6 5)
            (2 3 7 6)
            (4 5 9 8)
            (5 6 10 9)
            (6 7 11 10)
        );
    }
    BACK
    {
        type symmetryPlane;
        faces
        (
            (12 13 17 16)
            (13 14 18 17)
            (14 15 19 18)
            (16 17 21 20)
            (17 18 22 21)
            (18 19 23 22)
        );
    }
    TOP
    {
        type wall;
        faces
        (
            (8 9 21 20)
            (9 10 22 21)
            (10 11 23 22)
        );
    }
    LEFT
    {
        type wall;
        faces
        (
            (0 4 16 12)
            (4 8 20 16)
        );
    }
    RIGHT
    {
        type wall;
        faces
        (
            (3 7 19 15)
            (7 11 23 19)
        );
    }
    BOTTOM_FLUID
    {
        type wall;
        faces
        (
            (0 1 13 12)
            (1 2 14 13)
        );
    }
    FSI_FLUID
    {
        type wall;
        faces
        (
            (2 3 15 14)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
EOF
}

verify_front_back_boundary() {
    local case_dir="$1"
    python3 - "${case_dir}/constant/polyMesh/boundary" <<'PY'
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text(encoding="utf-8")
for name in ("FRONT", "BACK"):
    m = re.search(rf"{name}\s*\{{[^{{}}]*type\s+([A-Za-z]+)\s*;", text, flags=re.MULTILINE)
    if not m:
        raise SystemExit(f"Patch {name} not found in boundary file")
    ptype = m.group(1)
    if ptype != "symmetryPlane":
        raise SystemExit(f"Patch {name} type is {ptype}, expected symmetryPlane")
PY
}

make_case() {
    local name="$1"
    local nx="$2"
    local ny="$3"
    local nz="$4"
    local case_dir="${WORK_ROOT}/${name}"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"
    cp -a "${TEMPLATE_DIR}/." "${case_dir}/"

    write_alpha0 "${case_dir}"
    write_set_fields_dict "${case_dir}"
    write_block_mesh_dict "${case_dir}" "${nx}" "${ny}" "${nz}"
    patch_front_back_symmetry_fields "${case_dir}"

    blockMesh -case "${case_dir}" > "${case_dir}/log.blockMesh" 2>&1
    verify_front_back_boundary "${case_dir}"
    setFields -case "${case_dir}" > "${case_dir}/log.setFields" 2>&1
    checkMesh -case "${case_dir}" > "${case_dir}/log.checkMesh" 2>&1

    local cells_line
    cells_line="$(grep -m1 "cells:" "${case_dir}/log.checkMesh" | sed 's/^[[:space:]]*//')"
    echo "[${name}] ${cells_line}"
}

mkdir -p "${WORK_ROOT}"
extract_template

case "${MESH_LEVEL}" in
    coarse)
        make_case "coarse" 120 75 6
        ;;
    medium)
        make_case "medium" 180 110 8
        ;;
    fine)
        make_case "fine" 240 150 10
        ;;
    all)
        make_case "coarse" 120 75 6
        make_case "medium" 180 110 8
        make_case "fine" 240 150 10
        ;;
esac

echo "Mesh-independence cases created under: ${WORK_ROOT}"
if [[ "${MESH_LEVEL}" == "all" ]]; then
    echo "Generated cases: coarse, medium, fine"
else
    echo "Generated case: ${MESH_LEVEL}"
fi
