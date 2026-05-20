#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/build_calculix_modal_cases.py"

if [[ ! -f "${PY_SCRIPT}" ]]; then
    echo "Error: missing ${PY_SCRIPT}" >&2
    exit 1
fi

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
if [[ -z "${SOLID_CDB}" ]]; then
    SOLID_CDB="${SCRIPT_DIR}/DamFailure.in"
fi

OUT_DIR="${SCRIPT_DIR}/calculix_modal_inputs"
N_MODES=10
FLUID_STIFF_DIGITS=12
RUN_CCX=0
CCX_CMD="ccx_preCICE"
CEINTF_TOL="0.35"
INTERP_K="1"

usage() {
    cat <<'EOF'
Usage:
  run_calculix_modal_cases.sh [options] <fluid_cdb1> [fluid_cdb2 ...]

Options:
  --solid-cdb <path>   Structural reference CDB (default: auto-detect DamFailure.in)
  --out-dir <path>     Output folder for generated .inp files
  --n-modes <int>      Number of modes (default: 10)
  --fluid-stiff-digits Scientific-notation rounding digits for pseudo-fluid E grouping (default: 12)
  --ceintf-tol <float> Interface match tolerance passed to build script (default: 0.35)
  --interp-k <int>     Number of nearest solid nodes for interface interpolation (default: 1)
  --ccx-cmd <cmd>      CalculiX command passed to build script (default: ccx_preCICE)
  --run-ccx            Run ccx after generating .inp
  -h, --help           Show this help

Notes:
1. You can pass one fluid CDB or multiple fluid CDB files.
2. The same structural mesh from --solid-cdb is reused for all fluid meshes.
3. If no fluid CDB is provided, the same file as --solid-cdb is used.
EOF
}

FLUID_CDBS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --solid-cdb)
            SOLID_CDB="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
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
        --run-ccx)
            RUN_CCX=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            FLUID_CDBS+=("$1")
            shift
            ;;
    esac
done

if [[ ${#FLUID_CDBS[@]} -eq 0 ]]; then
    FLUID_CDBS=("${SOLID_CDB}")
fi

CMD=(
    python3 "${PY_SCRIPT}"
    --solid-cdb "${SOLID_CDB}"
    --out-dir "${OUT_DIR}"
    --n-modes "${N_MODES}"
    --fluid-stiff-digits "${FLUID_STIFF_DIGITS}"
    --ceintf-tol "${CEINTF_TOL}"
    --interp-k "${INTERP_K}"
    --ccx-cmd "${CCX_CMD}"
)
for f in "${FLUID_CDBS[@]}"; do
    CMD+=(--fluid-cdb "${f}")
done
if [[ "${RUN_CCX}" -eq 1 ]]; then
    CMD+=(--run-ccx)
fi

echo "Running: ${CMD[*]}"
"${CMD[@]}"
