#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BASELINE_ZIP="/root/Workspace/case_damfailure.zip"
WORK_ROOT="${SCRIPT_DIR}"
MESH_LEVEL=""
STAGE="all"
SOLVER="myInterFoam"
PARALLEL="1"
FSI_END_TIME=""

SOLID_CDB=""
for candidate in \
    "${SCRIPT_DIR}/DamFailure.in" \
    "/root/Workspace/modeNew/DamFailure.in"
do
    if [[ -f "${candidate}" ]]; then
        SOLID_CDB="${candidate}"
        break
    fi
done

N_MODES="10"
FLUID_STIFF_DIGITS="12"
CCX_CMD="ccx_preCICE"
CEINTF_TOL="0.35"
INTERP_K="1"

CDB_OUT_DIR="${WORK_ROOT}/cdb_inputs"
MODAL_OUT_DIR="${WORK_ROOT}/calculix_modal_inputs"

usage() {
    cat <<'EOF'
Usage:
  run_damfailure_pipeline.sh --mesh-level <coarse|medium|fine> [options]

Options:
  --mesh-level <level>        coarse | medium | fine (required)
  --stage <name>              all | mesh | modal | export | fsi (default: all)
  --baseline-zip <path>       Baseline zip for mesh generation
  --work-root <path>          Root folder containing coarse/medium/fine cases
  --solid-cdb <path>          Structural CDB (default: auto-detect)
  --solver <name>             FSI solver command (default: myInterFoam)
  --parallel <n>              MPI ranks for FSI stage (default: 1)
  --fsi-end-time <value>      Temporary endTime override for current FSI run only
  --n-modes <int>             Number of modes for modal solve (default: 10)
  --fluid-stiff-digits <int>  Pseudo-fluid E grouping digits (default: 12)
  --ceintf-tol <float>        Interface tolerance for modal coupling (default: 0.35)
  --interp-k <int>            Nearest nodes used in interface interpolation (default: 1)
  --ccx-cmd <cmd>             CalculiX command (default: ccx_preCICE)
  -h, --help                  Show this help

Examples:
  ./run_damfailure_pipeline.sh --mesh-level coarse --stage all
  ./run_damfailure_pipeline.sh --mesh-level fine --stage modal --ceintf-tol 0.35 --interp-k 1
  ./run_damfailure_pipeline.sh --mesh-level medium --stage fsi --parallel 4
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mesh-level)
            MESH_LEVEL="$2"
            shift 2
            ;;
        --stage)
            STAGE="$2"
            shift 2
            ;;
        --baseline-zip)
            BASELINE_ZIP="$2"
            shift 2
            ;;
        --work-root)
            WORK_ROOT="$2"
            shift 2
            ;;
        --solid-cdb)
            SOLID_CDB="$2"
            shift 2
            ;;
        --solver)
            SOLVER="$2"
            shift 2
            ;;
        --parallel)
            PARALLEL="$2"
            shift 2
            ;;
        --fsi-end-time)
            FSI_END_TIME="$2"
            shift 2
            ;;
        --n-modes)
            N_MODES="$2"
            shift 2
            ;;
        --fluid-stiff-digits)
            FLUID_STIFF_DIGITS="$2"
            shift 2
            ;;
        --ceintf-tol)
            CEINTF_TOL="$2"
            shift 2
            ;;
        --interp-k)
            INTERP_K="$2"
            shift 2
            ;;
        --ccx-cmd)
            CCX_CMD="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

case "${MESH_LEVEL}" in
    coarse|medium|fine) ;;
    *)
        echo "Error: --mesh-level is required and must be coarse|medium|fine." >&2
        exit 1
        ;;
esac

case "${STAGE}" in
    all|mesh|modal|export|fsi) ;;
    *)
        echo "Error: --stage must be all|mesh|modal|export|fsi." >&2
        exit 1
        ;;
esac

if ! [[ "${PARALLEL}" =~ ^[0-9]+$ ]] || [[ "${PARALLEL}" -lt 1 ]]; then
    echo "Error: --parallel must be a positive integer." >&2
    exit 1
fi

if [[ -n "${FSI_END_TIME}" ]]; then
    if ! [[ "${FSI_END_TIME}" =~ ^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$ ]]; then
        echo "Error: --fsi-end-time must be numeric." >&2
        exit 1
    fi
fi

if [[ -z "${SOLID_CDB}" || ! -f "${SOLID_CDB}" ]]; then
    echo "Error: solid CDB not found. Use --solid-cdb <path>." >&2
    exit 1
fi

CASE_DIR="${WORK_ROOT}/${MESH_LEVEL}"
FLUID_CDB="${CDB_OUT_DIR}/${MESH_LEVEL}_fluid.in"
MODAL_BASE="${MODAL_OUT_DIR}/${MESH_LEVEL}_fluid"
MODAL_DAT="${MODAL_BASE}.dat"
MODAL_CSV="${MODAL_OUT_DIR}/${MESH_LEVEL}_modal_summary.csv"
MODAL_JSON="${MODAL_OUT_DIR}/${MESH_LEVEL}_modal_summary.json"

MESH_SCRIPT="${SCRIPT_DIR}/generate_blockMesh_cases.sh"
CDB_SCRIPT="${SCRIPT_DIR}/../mesh/openfoam_mesh_to_cdb.py"
MODAL_SCRIPT="${SCRIPT_DIR}/../modal/run_calculix_modal_cases.sh"
EXPORT_SCRIPT="${SCRIPT_DIR}/../modal/export_modal_info.py"

for req in "${MESH_SCRIPT}" "${CDB_SCRIPT}" "${MODAL_SCRIPT}" "${EXPORT_SCRIPT}"; do
    if [[ ! -f "${req}" ]]; then
        echo "Error: missing required script: ${req}" >&2
        exit 1
    fi
done

run_mesh_stage() {
    echo "[STAGE] mesh: generating ${MESH_LEVEL} case..."
    "${MESH_SCRIPT}" --baseline-zip "${BASELINE_ZIP}" --work-root "${WORK_ROOT}" --mesh-level "${MESH_LEVEL}"
}

run_modal_stage() {
    if [[ ! -d "${CASE_DIR}" ]]; then
        echo "Error: case directory not found: ${CASE_DIR}. Run --stage mesh first." >&2
        exit 1
    fi

    echo "[STAGE] modal: converting OpenFOAM mesh to fluid CDB..."
    python3 "${SCRIPT_DIR}/../mesh/openfoam_mesh_to_cdb.py" --case "${CASE_DIR}" --out-dir "${CDB_OUT_DIR}"

    if [[ ! -f "${FLUID_CDB}" ]]; then
        echo "Error: fluid CDB was not generated: ${FLUID_CDB}" >&2
        exit 1
    fi

    echo "[STAGE] modal: building and solving CalculiX modal case..."
    "${SCRIPT_DIR}/../modal/run_calculix_modal_cases.sh" \
        --solid-cdb "${SOLID_CDB}" \
        --out-dir "${MODAL_OUT_DIR}" \
        --n-modes "${N_MODES}" \
        --fluid-stiff-digits "${FLUID_STIFF_DIGITS}" \
        --ceintf-tol "${CEINTF_TOL}" \
        --interp-k "${INTERP_K}" \
        --ccx-cmd "${CCX_CMD}" \
        --run-ccx \
        "${FLUID_CDB}"

    if [[ ! -f "${MODAL_DAT}" ]]; then
        echo "Error: modal output dat not found: ${MODAL_DAT}" >&2
        exit 1
    fi
}

run_export_stage() {
    if [[ ! -f "${MODAL_DAT}" ]]; then
        echo "Error: modal dat not found: ${MODAL_DAT}. Run --stage modal first." >&2
        exit 1
    fi

    echo "[STAGE] export: exporting modal summary..."
    python3 "${SCRIPT_DIR}/../modal/export_modal_info.py" --dat "${MODAL_DAT}" --csv "${MODAL_CSV}" --json "${MODAL_JSON}"
}

run_fsi_stage() {
    if [[ ! -d "${CASE_DIR}" ]]; then
        echo "Error: case directory not found: ${CASE_DIR}. Run --stage mesh first." >&2
        exit 1
    fi

    echo "[STAGE] fsi: running ${SOLVER} in ${CASE_DIR}..."
    local control_dict="${CASE_DIR}/system/controlDict"
    local control_backup=""
    local solver_rc=0
    if [[ -n "${FSI_END_TIME}" ]]; then
        control_backup="${control_dict}.pipeline.bak"
        cp "${control_dict}" "${control_backup}"
        python3 - "${control_dict}" "${FSI_END_TIME}" <<'PY'
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
end_time = sys.argv[2]
text = path.read_text(encoding="utf-8")
new_text, count = re.subn(r"(^\s*endTime\s+)([^;]+)(;)", rf"\g<1>{end_time}\g<3>", text, flags=re.MULTILINE)
if count != 1:
    raise SystemExit(f"Failed to update endTime in {path}")
path.write_text(new_text, encoding="utf-8")
PY
        echo "[STAGE] fsi: temporary endTime override = ${FSI_END_TIME}"
    fi

    set +e
    pushd "${CASE_DIR}" >/dev/null
    rm -f "log.${SOLVER}" "log.${SOLVER}_parallel" "log.decomposePar" "log.reconstructPar"

    if [[ "${PARALLEL}" -gt 1 ]]; then
        decomposePar -force > "log.decomposePar" 2>&1 \
            && mpirun -np "${PARALLEL}" "${SOLVER}" -parallel > "log.${SOLVER}_parallel" 2>&1 \
            && reconstructPar -latestTime > "log.reconstructPar" 2>&1
    else
        "${SOLVER}" > "log.${SOLVER}" 2>&1
    fi

    solver_rc=$?
    popd >/dev/null
    set -e

    if [[ -n "${control_backup}" && -f "${control_backup}" ]]; then
        mv "${control_backup}" "${control_dict}"
    fi
    if [[ "${solver_rc}" -ne 0 ]]; then
        return "${solver_rc}"
    fi
}

case "${STAGE}" in
    mesh)
        run_mesh_stage
        ;;
    modal)
        run_modal_stage
        ;;
    export)
        run_export_stage
        ;;
    fsi)
        run_fsi_stage
        ;;
    all)
        run_mesh_stage
        run_modal_stage
        run_export_stage
        run_fsi_stage
        ;;
esac

echo "[DONE] stage=${STAGE}, mesh-level=${MESH_LEVEL}"
