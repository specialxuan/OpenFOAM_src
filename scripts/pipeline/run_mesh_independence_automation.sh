#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORK_ROOT="${SCRIPT_DIR}"
BASELINE_ZIP="/root/Workspace/case_damfailure.zip"
SOLID_CDB="/root/Workspace/modeNew/DamFailure.in"
LEVELS="coarse,medium,fine"
STAGE="check"
PARALLEL="1"
SOLVER="myInterFoam"
CCX_CMD="ccx_preCICE"
N_MODES="10"
FSI_END_TIME=""
RUN_CCX=0
ALLOW_STALE_MODE=0
DRY_RUN=0

usage() {
    cat <<'EOF'
Usage:
  run_mesh_independence_automation.sh [options]

Stages:
  check      Verify inputs, dependencies, case/mode consistency, and existing outputs.
  mesh       Generate OpenFOAM mesh cases from the baseline zip.
  modal      Convert OpenFOAM meshes to CDB, build CalculiX inputs, optionally run ccx, export modal summaries.
  fsi        Run FSI cases after mode consistency checks.
  reconstruct Reconstruct latest time for parallel FSI runs.
  analyze    Summarize mesh, modal, FSI log, and reconstruction status.
  all        mesh -> modal -> fsi -> reconstruct -> analyze.

Options:
  --levels <csv>            Mesh levels, default: coarse,medium,fine.
  --stage <name>            Stage listed above, default: check.
  --baseline-zip <path>     Baseline OpenFOAM case archive.
  --work-root <path>        Root containing coarse/medium/fine cases.
  --solid-cdb <path>        Structural CDB, default: /root/Workspace/modeNew/DamFailure.in.
  --parallel <n>            MPI ranks for FSI, default: 1.
  --solver <cmd>            FSI solver, default: myInterFoam.
  --fsi-end-time <value>    Temporary endTime override for FSI runs.
  --ccx-cmd <cmd>           CalculiX command, default: ccx_preCICE.
  --run-ccx                 Actually run CalculiX in the modal stage.
  --n-modes <n>             Modal count for generated CalculiX input, default: 10.
  --allow-stale-mode        Allow FSI even if mode/FluidNodeCoor.csv node count differs from mesh points.
  --dry-run                 Print commands without executing stages.
  -h, --help                Show help.

Important:
  The FSI stage requires mode/FluidNodeCoor.csv to match the OpenFOAM mesh point
  count unless --allow-stale-mode is set. Use --stage modal --run-ccx to
  regenerate CalculiX FRD results and export matching fastDynamicFvMesh mode CSVs.
EOF
}

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

run_cmd() {
    log "RUN: $*"
    if [[ "${DRY_RUN}" -eq 0 ]]; then
        "$@"
    fi
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

split_levels() {
    python3 - "$LEVELS" <<'PY'
import sys
levels = [x.strip() for x in sys.argv[1].split(",") if x.strip()]
valid = {"coarse", "medium", "fine"}
bad = [x for x in levels if x not in valid]
if bad:
    raise SystemExit("invalid level(s): " + ",".join(bad))
print("\n".join(levels))
PY
}

foam_points_count() {
    local case_dir="$1"
    local points_file="${case_dir}/constant/polyMesh/points"
    [[ -f "${points_file}" ]] || return 1
    awk '/^[[:space:]]*[0-9]+[[:space:]]*$/{print $1; exit}' "${points_file}"
}

mode_nodes_count() {
    local case_dir="$1"
    local coor="${case_dir}/mode/FluidNodeCoor.csv"
    [[ -f "${coor}" ]] || return 1
    python3 - "${coor}" <<'PY'
import re
import sys
from pathlib import Path
line = Path(sys.argv[1]).read_text(errors="ignore").splitlines()[0]
vals = [float(x) for x in re.findall(r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line)]
if len(vals) < 3:
    raise SystemExit(1)
print(int(round(vals[1])))
PY
}

latest_time_dir() {
    local case_dir="$1"
    find "${case_dir}" -maxdepth 1 -type d -printf '%f\n' \
        | awk '/^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/' \
        | sort -g | tail -n 1
}

require_file() {
    [[ -f "$1" ]] || die "missing file: $1"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "missing command: $1"
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --levels) LEVELS="$2"; shift 2 ;;
            --stage) STAGE="$2"; shift 2 ;;
            --baseline-zip) BASELINE_ZIP="$2"; shift 2 ;;
            --work-root) WORK_ROOT="$2"; shift 2 ;;
            --solid-cdb) SOLID_CDB="$2"; shift 2 ;;
            --parallel) PARALLEL="$2"; shift 2 ;;
            --solver) SOLVER="$2"; shift 2 ;;
            --fsi-end-time) FSI_END_TIME="$2"; shift 2 ;;
            --ccx-cmd) CCX_CMD="$2"; shift 2 ;;
            --run-ccx) RUN_CCX=1; shift ;;
            --n-modes) N_MODES="$2"; shift 2 ;;
            --allow-stale-mode) ALLOW_STALE_MODE=1; shift ;;
            --dry-run) DRY_RUN=1; shift ;;
            -h|--help) usage; exit 0 ;;
            *) die "unknown argument: $1" ;;
        esac
    done

    case "${STAGE}" in
        check|mesh|modal|fsi|reconstruct|analyze|all) ;;
        *) die "--stage must be check|mesh|modal|fsi|reconstruct|analyze|all" ;;
    esac
    [[ "${PARALLEL}" =~ ^[1-9][0-9]*$ ]] || die "--parallel must be a positive integer"
}

check_environment() {
    require_file "${BASELINE_ZIP}"
    require_file "${SOLID_CDB}"
    require_file "${WORK_ROOT}/generate_blockMesh_cases.sh" || require_file "${SCRIPT_DIR}/generate_blockMesh_cases.sh"
    require_file "${WORK_ROOT}/openfoam_mesh_to_cdb.py" || require_file "${SCRIPT_DIR}/../mesh/openfoam_mesh_to_cdb.py"
    require_file "${WORK_ROOT}/run_calculix_modal_cases.sh" || require_file "${SCRIPT_DIR}/../modal/run_calculix_modal_cases.sh"
    require_file "${WORK_ROOT}/export_modal_info.py" || require_file "${SCRIPT_DIR}/../modal/export_modal_info.py"
    require_file "${WORK_ROOT}/export_calculix_modes_to_fastdynamic.py" || require_file "${SCRIPT_DIR}/../modal/export_calculix_modes_to_fastdynamic.py"
    require_cmd python3
    require_cmd unzip
    require_cmd blockMesh
    require_cmd setFields
    require_cmd checkMesh
    require_cmd "${SOLVER}"
    if [[ "${RUN_CCX}" -eq 1 ]]; then
        require_cmd "${CCX_CMD}"
    fi
    if [[ "${PARALLEL}" -gt 1 ]]; then
        require_cmd decomposePar
        require_cmd mpirun
        require_cmd reconstructPar
    fi
}

check_mode_consistency() {
    local level="$1"
    local case_dir="${WORK_ROOT}/${level}"
    local mesh_n mode_n
    mesh_n="$(foam_points_count "${case_dir}")" || die "${level}: missing OpenFOAM mesh points"
    mode_n="$(mode_nodes_count "${case_dir}")" || die "${level}: missing mode/FluidNodeCoor.csv"

    if [[ "${mesh_n}" != "${mode_n}" ]]; then
        local msg="${level}: mesh points=${mesh_n}, mode nodes=${mode_n}; FSI mode files do not match this mesh"
        if [[ "${ALLOW_STALE_MODE}" -eq 1 ]]; then
            echo "WARNING: ${msg}" >&2
        else
            die "${msg}. Generate matching FluidNodeCoor/FluidNodeDisp*.csv before FSI, or pass --allow-stale-mode for a diagnostic-only run."
        fi
    else
        log "${level}: mode nodes match mesh points (${mesh_n})"
    fi
}

stage_mesh() {
    local levels=()
    mapfile -t levels < <(split_levels)
    for level in "${levels[@]}"; do
        run_cmd "${SCRIPT_DIR}/generate_blockMesh_cases.sh" \
            --baseline-zip "${BASELINE_ZIP}" \
            --work-root "${WORK_ROOT}" \
            --mesh-level "${level}"
    done
}

stage_modal() {
    mkdir -p "${WORK_ROOT}/cdb_inputs" "${WORK_ROOT}/calculix_modal_inputs"
    local dats=()
    local levels=()
    mapfile -t levels < <(split_levels)
    for level in "${levels[@]}"; do
        local case_dir="${WORK_ROOT}/${level}"
        local fluid_cdb="${WORK_ROOT}/cdb_inputs/${level}_fluid.in"
        local modal_dat="${WORK_ROOT}/calculix_modal_inputs/${level}_fluid.dat"
        local modal_inp="${WORK_ROOT}/calculix_modal_inputs/${level}_fluid.inp"
        local modal_frd="${WORK_ROOT}/calculix_modal_inputs/${level}_fluid.frd"
        run_cmd python3 "${SCRIPT_DIR}/../mesh/openfoam_mesh_to_cdb.py" --case "${case_dir}" --out-dir "${WORK_ROOT}/cdb_inputs"
        local ccx_args=()
        if [[ "${RUN_CCX}" -eq 1 ]]; then
            ccx_args+=(--run-ccx)
        fi
        run_cmd "${SCRIPT_DIR}/../modal/run_calculix_modal_cases.sh" \
            --solid-cdb "${SOLID_CDB}" \
            --out-dir "${WORK_ROOT}/calculix_modal_inputs" \
            --n-modes "${N_MODES}" \
            --ccx-cmd "${CCX_CMD}" \
            "${ccx_args[@]}" \
            "${fluid_cdb}"
        if [[ -f "${modal_dat}" ]]; then
            dats+=("${modal_dat}")
            run_cmd python3 "${SCRIPT_DIR}/../modal/export_modal_info.py" \
                --dat "${modal_dat}" \
                --csv "${WORK_ROOT}/calculix_modal_inputs/${level}_modal_summary.csv" \
                --json "${WORK_ROOT}/calculix_modal_inputs/${level}_modal_summary.json"
            if [[ "${RUN_CCX}" -eq 1 && -f "${modal_inp}" && -f "${modal_frd}" ]]; then
                run_cmd python3 "${SCRIPT_DIR}/../modal/export_calculix_modes_to_fastdynamic.py" \
                    --inp "${modal_inp}" \
                    --dat "${modal_dat}" \
                    --frd "${modal_frd}" \
                    --out-mode-dir "${case_dir}/mode" \
                    --n-modes "${N_MODES}" \
                    --fluid-para-source "${case_dir}/mode/FluidPara.csv"
            elif [[ "${RUN_CCX}" -eq 0 ]]; then
                log "SKIP fastDynamicFvMesh mode export for ${level}: --run-ccx was not set"
            fi
        else
            log "SKIP modal export for ${level}: ${modal_dat} not found"
        fi
    done

    if [[ "${#dats[@]}" -gt 0 ]]; then
        local args=()
        for d in "${dats[@]}"; do
            args+=(--dat "${d}")
        done
        run_cmd python3 "${SCRIPT_DIR}/../modal/export_modal_info.py" \
            "${args[@]}" \
            --csv "${WORK_ROOT}/calculix_modal_inputs/mesh_independence_modal_summary.csv" \
            --json "${WORK_ROOT}/calculix_modal_inputs/mesh_independence_modal_summary.json"
    fi
}

stage_fsi() {
    local levels=()
    mapfile -t levels < <(split_levels)
    for level in "${levels[@]}"; do
        check_mode_consistency "${level}"
        local case_dir="${WORK_ROOT}/${level}"
        local control="${case_dir}/system/controlDict"
        local bak=""
        if [[ -n "${FSI_END_TIME}" ]]; then
            bak="${control}.automation.bak"
            run_cmd cp "${control}" "${bak}"
            run_cmd foamDictionary "${control}" -entry endTime -set "${FSI_END_TIME}"
        fi

        if [[ "${PARALLEL}" -gt 1 ]]; then
            if [[ "${DRY_RUN}" -eq 1 ]]; then
                log "RUN: decomposePar -case ${case_dir} -force > ${case_dir}/log.decomposePar 2>&1"
                log "RUN: mpirun --allow-run-as-root -np ${PARALLEL} ${SOLVER} -case ${case_dir} -parallel > ${case_dir}/log.${SOLVER}_parallel 2>&1"
            else
                decomposePar -case "${case_dir}" -force > "${case_dir}/log.decomposePar" 2>&1
                mpirun --allow-run-as-root -np "${PARALLEL}" "${SOLVER}" -case "${case_dir}" -parallel > "${case_dir}/log.${SOLVER}_parallel" 2>&1 < /dev/null
            fi
        else
            if [[ "${DRY_RUN}" -eq 1 ]]; then
                log "RUN: ${SOLVER} -case ${case_dir} > ${case_dir}/log.${SOLVER} 2>&1"
            else
                "${SOLVER}" -case "${case_dir}" > "${case_dir}/log.${SOLVER}" 2>&1 < /dev/null
            fi
        fi

        if [[ -n "${bak}" && -f "${bak}" && "${DRY_RUN}" -eq 0 ]]; then
            mv "${bak}" "${control}"
        fi
    done
}

stage_reconstruct() {
    if [[ "${PARALLEL}" -le 1 ]]; then
        log "parallel=1; reconstruction is not required"
        return
    fi
    local levels=()
    mapfile -t levels < <(split_levels)
    for level in "${levels[@]}"; do
        local case_dir="${WORK_ROOT}/${level}"
        if [[ "${DRY_RUN}" -eq 1 ]]; then
            log "RUN: reconstructPar -case ${case_dir} -latestTime > ${case_dir}/log.reconstructPar 2>&1"
        else
            reconstructPar -case "${case_dir}" -latestTime > "${case_dir}/log.reconstructPar" 2>&1
        fi
    done
}

stage_analyze() {
    local out_dir="${WORK_ROOT}/analysis"
    local out_csv="${out_dir}/mesh_independence_summary.csv"
    mkdir -p "${out_dir}"
    python3 - "${WORK_ROOT}" "${LEVELS}" "${out_csv}" <<'PY'
import csv
import re
import sys
from pathlib import Path

work = Path(sys.argv[1])
levels = [x.strip() for x in sys.argv[2].split(",") if x.strip()]
out = Path(sys.argv[3])

def foam_count(path: Path):
    if not path.exists():
        return ""
    for line in path.read_text(errors="ignore").splitlines():
        if re.fullmatch(r"\s*\d+\s*", line):
            return int(line)
    return ""

def mode_count(path: Path):
    if not path.exists():
        return ""
    line = path.read_text(errors="ignore").splitlines()[0]
    vals = re.findall(r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line)
    return int(round(float(vals[1]))) if len(vals) >= 2 else ""

def latest_time(case: Path):
    vals = []
    for p in case.iterdir() if case.exists() else []:
        if p.is_dir():
            try:
                vals.append(float(p.name))
            except ValueError:
                pass
    return max(vals) if vals else ""

def checkmesh_cells(case: Path):
    log = case / "log.checkMesh"
    if not log.exists():
        return ""
    m = re.search(r"^\s*cells:\s+(\d+)", log.read_text(errors="ignore"), re.M)
    return int(m.group(1)) if m else ""

def parse_log(case: Path):
    logs = sorted(case.glob("log.myInterFoam*"))
    if not logs:
        return "", "", "", ""
    text = logs[-1].read_text(errors="ignore")
    m = list(re.finditer(r"ExecutionTime =\s*([0-9.eE+-]+) s\s+ClockTime =\s*([0-9.eE+-]+) s", text))
    exe = m[-1].group(1) if m else ""
    clk = m[-1].group(2) if m else ""
    end_ok = "yes" if re.search(r"^End\s*$", text, re.M) else "no"
    iters = re.findall(r"FSI/PIMPLE Iteration\s+(\d+)", text)
    last_iter = iters[-1] if iters else ""
    return exe, clk, end_ok, last_iter

rows = []
for level in levels:
    case = work / level
    meta = work / "calculix_modal_inputs" / f"{level}_fluid.json"
    modal_csv = work / "calculix_modal_inputs" / f"{level}_modal_summary.csv"
    first_freq = ""
    if modal_csv.exists():
        with modal_csv.open(newline="") as f:
            r = csv.DictReader(f)
            first = next(r, None)
            if first:
                first_freq = first.get("frequency_cycles_per_time", "")
    exe, clk, end_ok, last_iter = parse_log(case)
    rows.append({
        "level": level,
        "cells": checkmesh_cells(case),
        "points": foam_count(case / "constant/polyMesh/points"),
        "mode_nodes": mode_count(case / "mode/FluidNodeCoor.csv"),
        "latest_time": latest_time(case),
        "modal_first_frequency": first_freq,
        "fsi_execution_time_s": exe,
        "fsi_clock_time_s": clk,
        "fsi_end_ok": end_ok,
        "fsi_last_outer_iter": last_iter,
        "modal_meta_exists": "yes" if meta.exists() else "no",
    })

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
    writer.writeheader()
    writer.writerows(rows)
print(out)
PY
    log "analysis CSV: ${out_csv}"
}

stage_check() {
    check_environment
    while IFS= read -r level; do
        local case_dir="${WORK_ROOT}/${level}"
        if [[ -d "${case_dir}" ]]; then
            local mesh_n mode_n
            mesh_n="$(foam_points_count "${case_dir}" || true)"
            mode_n="$(mode_nodes_count "${case_dir}" || true)"
            log "${level}: mesh_points=${mesh_n:-missing}, mode_nodes=${mode_n:-missing}, latest_time=$(latest_time_dir "${case_dir}" || true)"
        else
            log "${level}: case directory missing"
        fi
    done < <(split_levels)
}

parse_args "$@"
check_environment

case "${STAGE}" in
    check) stage_check ;;
    mesh) stage_mesh ;;
    modal) stage_modal ;;
    fsi) stage_fsi ;;
    reconstruct) stage_reconstruct ;;
    analyze) stage_analyze ;;
    all)
        stage_mesh
        stage_modal
        stage_fsi
        stage_reconstruct
        stage_analyze
        ;;
esac

log "DONE stage=${STAGE} levels=${LEVELS}"
