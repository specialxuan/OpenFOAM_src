#!/usr/bin/env bash
set -euo pipefail

workspace_root="${WORKSPACE_ROOT:-/root/Workspace}"
src_root="${SRC_ROOT:-/root/OpenFOAM/user-v2412/src}"
run_root="${RUN_ROOT:-/tmp/fastDynamicFvMesh_smoke_$(date +%Y%m%d_%H%M%S)}"
solver_override="${SMOKE_SOLVER:-}"
end_time="${SMOKE_END_TIME:-}"
dry_run=0

usage()
{
    cat <<'USAGE'
Usage: scripts/smoke_fastDynamicFvMesh.sh [--dry-run] [--run-root DIR] [--solver NAME] [--end-time VALUE]

Copies the baseline cases from /root/Workspace into a temporary run directory,
builds fastDynamicFvMesh, and runs a short solver smoke test when possible.
By default, each case uses its own system/controlDict application entry.

Environment overrides:
  WORKSPACE_ROOT   default: /root/Workspace
  SRC_ROOT         default: /root/OpenFOAM/user-v2412/src
  RUN_ROOT         default: /tmp/fastDynamicFvMesh_smoke_<timestamp>
  SMOKE_SOLVER     optional solver override for all cases
  SMOKE_END_TIME   optional controlDict endTime override
USAGE
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --dry-run)
            dry_run=1
            shift
            ;;
        --run-root)
            run_root="$2"
            shift 2
            ;;
        --solver)
            solver_override="$2"
            shift 2
            ;;
        --end-time)
            end_time="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

archives=(
    "$workspace_root/case_damfailure.zip"
    "$workspace_root/case_transient.zip"
)

echo "Smoke run root: $run_root"
echo "Source root: $src_root"
if [ -n "$solver_override" ]; then
    echo "Solver override: $solver_override"
else
    echo "Solver: use each case system/controlDict application"
fi

for archive in "${archives[@]}"; do
    if [ ! -f "$archive" ]; then
        echo "Missing baseline archive: $archive" >&2
        exit 1
    fi
done

if [ "$dry_run" -eq 1 ]; then
    echo "Dry run only. Baseline archives exist:"
    printf '  %s\n' "${archives[@]}"
    exit 0
fi

mkdir -p "$run_root/logs"

echo "Building fastDynamicFvMesh..."
(
    cd "$src_root"
    wmake libso fastDynamicFvMesh
) 2>&1 | tee "$run_root/logs/wmake_fastDynamicFvMesh.log"

for archive in "${archives[@]}"; do
    case_name="$(basename "$archive" .zip)"
    case_dir="$run_root/$case_name"

    echo "Preparing $case_name..."
    mkdir -p "$case_dir"
    if command -v unzip >/dev/null 2>&1; then
        unzip -q "$archive" -d "$case_dir"
    elif command -v python3 >/dev/null 2>&1; then
        python3 -m zipfile -e "$archive" "$case_dir"
    else
        echo "Neither unzip nor python3 is available to extract $archive" >&2
        exit 1
    fi

    control_dict="$(find "$case_dir" -maxdepth 4 -type f -path '*/system/controlDict' | head -n 1)"
    if [ -n "$control_dict" ]; then
        actual_case="${control_dict%/system/controlDict}"
    else
        actual_case="$case_dir"
    fi

    if [ -n "$end_time" ] && [ -f "$actual_case/system/controlDict" ]; then
        foamDictionary "$actual_case/system/controlDict" -entry endTime -set "$end_time" >/dev/null
    fi

    actual_solver="$solver_override"
    if [ -z "$actual_solver" ]; then
        if [ ! -f "$actual_case/system/controlDict" ]; then
            echo "No system/controlDict found for $case_name; pass --solver explicitly." >&2
            exit 1
        fi

        actual_solver="$(foamDictionary "$actual_case/system/controlDict" -entry application 2>/dev/null | awk '{print $2}' | tr -d ';')"
        if [ -z "$actual_solver" ]; then
            echo "Could not read application from $actual_case/system/controlDict" >&2
            exit 1
        fi
    fi

    echo "Running $actual_solver in $actual_case..."
    (
        cd "$actual_case"
        "$actual_solver"
    ) >"$run_root/logs/${case_name}_${actual_solver}.log" 2>&1 || {
        echo "Solver smoke failed for $case_name. See $run_root/logs/${case_name}_${actual_solver}.log" >&2
        exit 1
    }

    echo "Completed $case_name"
done

echo "Smoke run completed. Logs: $run_root/logs"
