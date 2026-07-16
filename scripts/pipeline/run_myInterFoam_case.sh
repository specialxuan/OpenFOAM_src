#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./run_myInterFoam_case.sh [-f|-r|-c|-s] [NPROCS]

Modes:
  -f  Fresh start (default): clean generated results, run checkMesh,
      decomposePar, solve in parallel, then reconstruct mesh and fields.
  -r  Restart from startTime with existing decomposition: clean generated
      result times/logs, keep processor*/constant and processor*/0, then solve.
  -c  Continue from latest written time with existing decomposition.
  -s  Solve with existing decomposition; do not change startFrom, then
      reconstruct mesh and fields. Useful for manual/diagnostic continuation.

Arguments:
  NPROCS   MPI process count (default: value from system/decomposeParDict, or 8)

The script is case-local: copy it into an OpenFOAM case directory and run it
from there. It resolves OpenFOAM commands from the loaded OpenFOAM environment.
EOF
}

MODE="f"
while getopts ":frcsh" opt; do
  case "$opt" in
    f) MODE="f" ;;
    r) MODE="r" ;;
    c) MODE="c" ;;
    s) MODE="s" ;;
    h) usage; exit 0 ;;
    \?) echo "ERROR: invalid option: -$OPTARG" >&2; usage; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$case_dir"

source_openfoam() {
  if [ -n "${WM_PROJECT_DIR:-}" ]; then
    return
  fi

  local bashrc="${OPENFOAM_BASHRC:-/usr/lib/openfoam/openfoam2412/etc/bashrc}"
  if [ ! -f "$bashrc" ]; then
    echo "ERROR: OpenFOAM bashrc not found: $bashrc" >&2
    echo "Set OPENFOAM_BASHRC or source OpenFOAM before running." >&2
    exit 1
  fi

  set +u
  source "$bashrc"
  set -u
}

source_openfoam

default_nprocs() {
  python3 - <<'PY'
import re
from pathlib import Path
p = Path('system/decomposeParDict')
if p.exists():
    m = re.search(r'numberOfSubdomains\s+(\d+)\s*;', p.read_text())
    if m:
        print(m.group(1))
        raise SystemExit
print('8')
PY
}

NPROCS=${1:-$(default_nprocs)}
if [ "$#" -gt 1 ]; then
  echo "ERROR: too many arguments." >&2
  usage
  exit 1
fi
if ! [[ "$NPROCS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: NPROCS must be a positive integer, got '$NPROCS'." >&2
  exit 1
fi

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: required command not found in PATH: $1" >&2
    exit 1
  fi
}

need_cmd python3
need_cmd mpirun
need_cmd myInterFoam
case "$MODE" in
  f) need_cmd checkMesh; need_cmd decomposePar; need_cmd reconstructParMesh; need_cmd reconstructPar ;;
  r|c) need_cmd reconstructParMesh; need_cmd reconstructPar ;;
  s) ;;
esac

if [ "$(id -u)" -eq 0 ]; then
  export OMPI_ALLOW_RUN_AS_ROOT=1
  export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
fi

update_decompose_dict() {
  python3 - "$1" <<'PY'
import re
import sys
from pathlib import Path
p = Path('system/decomposeParDict')
if not p.exists():
    raise SystemExit('ERROR: system/decomposeParDict not found')
s = p.read_text()
s_new, n = re.subn(r'(numberOfSubdomains\s+)\d+(\s*;)', rf'\g<1>{sys.argv[1]}\2', s, count=1)
if n == 0:
    raise SystemExit('ERROR: numberOfSubdomains not found in system/decomposeParDict')
p.write_text(s_new)
PY
}

set_start_from() {
  python3 - "$1" <<'PY'
import re
import sys
from pathlib import Path
target = sys.argv[1]
if target not in {'startTime', 'latestTime'}:
    raise SystemExit(f'ERROR: unsupported startFrom: {target}')
p = Path('system/controlDict')
if not p.exists():
    raise SystemExit('ERROR: system/controlDict not found')
s = p.read_text()
s_new, n = re.subn(r'(^\s*startFrom\s+)\w+(\s*;)', rf'\g<1>{target}\2', s, count=1, flags=re.M)
if n == 0:
    raise SystemExit('ERROR: startFrom not found in system/controlDict')
p.write_text(s_new)
PY
}

clean_root_results() {
  rm -rf 0.* [1-9]* log.* *.log postProcessing dynamicCode case.foam
}

clean_processor_results_keep_mesh() {
  shopt -s nullglob
  local proc
  for proc in processor*; do
    [ -d "$proc" ] || continue
    find "$proc" -mindepth 1 -maxdepth 1 -type d \( -name '0.*' -o -name '[1-9]*' \) -exec rm -rf {} +
    find "$proc" -mindepth 1 -maxdepth 1 -type f \( -name '0.*' -o -name '[1-9]*' \) -delete
  done
  shopt -u nullglob
}

ensure_processors_exist() {
  if ! find . -maxdepth 1 -type d -name 'processor*' | grep -q .; then
    echo "ERROR: no processor* directories found. Run with -f first." >&2
    exit 1
  fi
}

validate_amr_controls() {
  python3 - <<'PY'
import re
from pathlib import Path

p = Path('constant/dynamicMeshDict')
if not p.exists():
    raise SystemExit(0)

s = p.read_text()

def scalar(name):
    m = re.search(rf'\b{name}\s+([^;]+);', s)
    return None if m is None else float(m.group(1))

def switch(name):
    m = re.search(rf'\b{name}\s+(\w+)\s*;', s)
    return None if m is None else m.group(1).lower() in {'true', 'yes', 'on', '1'}

use_grad = switch('useGradIndicator')
lower = scalar('lowerRefineLevel')
unref = scalar('unrefineLevel')

if use_grad and lower is not None and unref is not None and unref >= lower:
    raise SystemExit(
        'ERROR: unrefineLevel %.6g must be lower than lowerRefineLevel %.6g for gradient AMR hysteresis.'
        % (unref, lower)
    )
PY
}

show_failure_tail() {
  local log_file=$1
  echo "---> Last solver log lines from $log_file" >&2
  tail -n 80 "$log_file" >&2 || true
  echo "---> Failure markers" >&2
  grep -n -E 'FOAM FATAL ERROR|MPI_ABORT|Floating point|Signal|Max\(alpha\.water\)|Runtime AMR timing|Selected [0-9]+ cells|Time =' "$log_file" | tail -n 80 >&2 || true
}

reconstruct_all() {
  echo "---> Reconstructing AMR mesh"
  reconstructParMesh -constant > log.reconstructParMesh 2>&1
  echo "---> Reconstructing fields"
  reconstructPar > log.reconstructPar 2>&1
  touch case.foam
}

echo "---> Case directory: $case_dir"
echo "---> Mode: $MODE"
echo "---> MPI ranks: $NPROCS"

validate_amr_controls

case "$MODE" in
  f)
    echo "---> Fresh clean"
    clean_root_results
    find . -maxdepth 1 -type d -name 'processor*' -exec rm -rf {} +
    set_start_from startTime
    update_decompose_dict "$NPROCS"
    echo "---> Running checkMesh"
    checkMesh > log.checkMesh 2>&1
    echo "---> Running decomposePar"
    decomposePar -force > log.decomposePar 2>&1
    ;;
  r)
    ensure_processors_exist
    echo "---> Restart: cleaning results, keeping decomposition"
    clean_root_results
    clean_processor_results_keep_mesh
    set_start_from startTime
    ;;
  c)
    ensure_processors_exist
    echo "---> Continue: keeping results and using latestTime"
    set_start_from latestTime
    ;;
  s)
    ensure_processors_exist
    echo "---> Solve-only: keeping startFrom"
    ;;
esac

LOG_SOLVER="log.myInterFoam.run"
echo "---> Running myInterFoam in parallel"
: > "$LOG_SOLVER"
mpirun --allow-run-as-root -np "$NPROCS" myInterFoam -parallel > "$LOG_SOLVER" 2>&1 &
solver_pid=$!

tail --pid="$solver_pid" -n 0 -F "$LOG_SOLVER" \
  | grep --line-buffered -E '^(Time = |ExecutionTime =|Runtime AMR timing)' || true

set +e
wait "$solver_pid"
solver_status=$?
set -e
if [ "$solver_status" -ne 0 ]; then
  echo "ERROR: myInterFoam failed; see $LOG_SOLVER" >&2
  show_failure_tail "$LOG_SOLVER"
  exit "$solver_status"
fi

reconstruct_all

echo "---> Complete"
