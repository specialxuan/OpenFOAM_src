#!/usr/bin/env bash
set -euo pipefail

CASE_ROOT="/root/Workspace/precice_run"
OUTPUT=""
FORCE=0
DRY_RUN=0
KEEP_DECOMPOSITION=1

usage() {
    cat <<EOF
Usage: $0 [options]

Stage a clean, replayable dual-participant preCICE case without modifying its source.

  --case-root <dir>        workflow root containing pc_case
  --output <file.zip>      output archive (default <case-root>_clean.zip)
  --rebuild-decomposition  omit processor base meshes and repartition on replay
  --force                  overwrite an existing archive
  --dry-run                validate and print the staging plan only
  -h, --help               show this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --case-root) CASE_ROOT="${2:?--case-root needs a value}"; shift 2 ;;
        --output) OUTPUT="${2:?--output needs a value}"; shift 2 ;;
        --rebuild-decomposition) KEEP_DECOMPOSITION=0; shift ;;
        --force) FORCE=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "[PRECICE-PACK ERR] unknown option: $1" >&2; exit 2 ;;
    esac
done

CASE_ROOT="$(realpath -e "$CASE_ROOT")"
CASE_DIR="$CASE_ROOT/pc_case"
NAME="$(basename "$CASE_ROOT")"
OUTPUT="${OUTPUT:-$(dirname "$CASE_ROOT")/${NAME}_clean.zip}"
OUTPUT="$(realpath -m "$OUTPUT")"
if [[ "$OUTPUT" != *.zip ]]; then
    echo "[PRECICE-PACK ERR] --output must end in .zip: $OUTPUT" >&2
    exit 2
fi
case "$OUTPUT" in "$CASE_ROOT"/*) echo "[PRECICE-PACK ERR] output cannot be inside source" >&2; exit 2;; esac
if [ -e "$OUTPUT" ] && [ "$FORCE" -ne 1 ]; then
    echo "[PRECICE-PACK ERR] output exists; use --force: $OUTPUT" >&2
    exit 2
fi
command -v zip >/dev/null 2>&1 || { echo "[PRECICE-PACK ERR] zip not found" >&2; exit 2; }

REQUIRED=(
    "$CASE_ROOT/CaseDescription.md" "$CASE_ROOT/precice_result.mp4"
    "$CASE_DIR/precice-config.xml" "$CASE_DIR/run-coupled.sh"
    "$CASE_DIR/fluid-openfoam/0" "$CASE_DIR/fluid-openfoam/case.foam"
    "$CASE_DIR/fluid-openfoam/constant/polyMesh/points"
    "$CASE_DIR/fluid-openfoam/constant/polyMesh/faces"
    "$CASE_DIR/fluid-openfoam/constant/polyMesh/owner"
    "$CASE_DIR/fluid-openfoam/constant/polyMesh/neighbour"
    "$CASE_DIR/fluid-openfoam/constant/polyMesh/boundary"
    "$CASE_DIR/fluid-openfoam/constant/dynamicMeshDict"
    "$CASE_DIR/fluid-openfoam/system/controlDict"
    "$CASE_DIR/fluid-openfoam/system/preciceDict"
    "$CASE_DIR/fluid-openfoam/system/fvSchemes"
    "$CASE_DIR/fluid-openfoam/system/fvSolution"
    "$CASE_DIR/fluid-openfoam/system/setFieldsDict"
    "$CASE_DIR/fluid-openfoam/system/decomposeParDict"
    "$CASE_DIR/solid-calculix/dam_gate.inp" "$CASE_DIR/solid-calculix/config.yml"
)
for path in "${REQUIRED[@]}"; do
    [ -e "$path" ] || { echo "[PRECICE-PACK ERR] missing: $path" >&2; exit 2; }
done

echo "[PRECICE-PACK] source=$CASE_ROOT"
echo "[PRECICE-PACK] output=$OUTPUT"
echo "[PRECICE-PACK] preserve decomposition=$KEEP_DECOMPOSITION"
echo "[PRECICE-PACK] stage excludes transient results, logs, frames, and postprocess; source remains unchanged"
[ "$DRY_RUN" -eq 1 ] && exit 0

OUTPUT_DIR="$(dirname "$OUTPUT")"
mkdir -p "$OUTPUT_DIR"
TMP="$(mktemp -d "$OUTPUT_DIR/.${NAME}.pack.XXXXXX")"
ARCHIVE_TMP=""
trap 'rm -rf "$TMP"; [ -z "$ARCHIVE_TMP" ] || rm -f "$ARCHIVE_TMP"' EXIT
STAGE="$TMP/$NAME"
STAGE_CASE="$STAGE/pc_case"
mkdir -p "$STAGE_CASE"
cp -a "$CASE_ROOT/CaseDescription.md" "$CASE_ROOT/precice_result.mp4" "$STAGE/"
[ -d "$CASE_ROOT/case_settings" ] && cp -a "$CASE_ROOT/case_settings" "$STAGE/"
cp -a "$CASE_DIR/precice-config.xml" "$CASE_DIR/run-coupled.sh" "$STAGE_CASE/"
cp -a "$CASE_DIR/fluid-openfoam" "$CASE_DIR/solid-calculix" "$STAGE_CASE/"

for path in "$STAGE_CASE/fluid-openfoam"/*; do
    [ -d "$path" ] || continue
    name="$(basename "$path")"
    [[ "$name" =~ ^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$ ]] && [ "$name" != "0" ] && rm -rf "$path"
done
for processor in "$STAGE_CASE/fluid-openfoam"/processor*; do
    [ -d "$processor" ] || continue
    for path in "$processor"/*; do
        [ -d "$path" ] || continue
        name="$(basename "$path")"
        [[ "$name" =~ ^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$ ]] && [ "$name" != "0" ] && rm -rf "$path"
    done
done
find "$STAGE_CASE" -type f \( -name 'log.*' -o -name '*.log' -o -name '*.sta' -o -name '*.cvg' \) -delete
rm -rf "$STAGE_CASE/precice-run" "$STAGE_CASE/postprocess" \
    "$STAGE_CASE/fluid-openfoam/postprocess" "$STAGE/postprocess" "$STAGE/frames"

if [ "$KEEP_DECOMPOSITION" -eq 0 ]; then
    rm -rf "$STAGE_CASE/fluid-openfoam"/processor*
else
    for processor in "$STAGE_CASE/fluid-openfoam"/processor*; do
        [ -d "$processor" ] || continue
        find "$processor" -mindepth 1 -maxdepth 1 ! -name constant -exec rm -rf {} +
        [ -d "$processor/constant/polyMesh" ] || rm -rf "$processor"
    done
fi

cat > "$STAGE/run_case.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE="$ROOT/pc_case"
NPROCS="${NPROCS:-8}"
SOLID_THREADS="${SOLID_THREADS:-8}"
FOAM_BASHRC="${FOAM_BASHRC:-/usr/lib/openfoam/openfoam2412/etc/bashrc}"
[ "${NPROCS##*[!0-9]*}" ] && [ "$NPROCS" -gt 0 ] || { echo "NPROCS must be a positive integer" >&2; exit 2; }
[ "${SOLID_THREADS##*[!0-9]*}" ] && [ "$SOLID_THREADS" -gt 0 ] || { echo "SOLID_THREADS must be a positive integer" >&2; exit 2; }
source "$FOAM_BASHRC"
export OMP_NUM_THREADS="$SOLID_THREADS"
export CCX_NPROC_EQUATION_SOLVER="$SOLID_THREADS"
cd "$CASE/fluid-openfoam"
shopt -s nullglob
checkMesh
setFields
python3 - "$NPROCS" system/decomposeParDict <<'PY'
import re
import sys

path = sys.argv[2]
text = open(path, encoding="utf-8").read()
updated, count = re.subn(r'(^\s*numberOfSubdomains\s+)\d+(\s*;)', r'\g<1>' + sys.argv[1] + r'\2', text, count=1, flags=re.M)
if count != 1:
    raise SystemExit("numberOfSubdomains missing from " + path)
open(path, "w", encoding="utf-8").write(updated)
PY
preserved=(processor*/constant/polyMesh)
if [ "${#preserved[@]}" -eq "$NPROCS" ]; then
    decomposePar -fields
else
    rm -rf processor*
    decomposePar -force
fi
if [ "$(id -u)" -eq 0 ]; then
    MPI=(mpirun --allow-run-as-root -np "$NPROCS")
else
    MPI=(mpirun -np "$NPROCS")
fi
cd "$CASE"
mkdir -p precice-run
(cd solid-calculix && ccx_preCICE -i dam_gate -precice-participant Solid > ../precice-run/log.solid 2>&1) &
solid_pid=$!
(cd fluid-openfoam && "${MPI[@]}" interFoam -parallel > ../precice-run/log.fluid 2>&1) &
fluid_pid=$!
if wait "$fluid_pid"; then fluid_status=0; else fluid_status=$?; fi
if wait "$solid_pid"; then solid_status=0; else solid_status=$?; fi
[ "$fluid_status" -eq 0 ] && [ "$solid_status" -eq 0 ]
EOF
chmod +x "$STAGE/run_case.sh"

cat > "$STAGE/README.md" <<EOF
# $NAME replay package

The package includes the final \`precice_result.mp4\`, the case description, both preCICE participants, \`precice-config.xml\`, \`run-coupled.sh\`, and \`run_case.sh\`. It excludes transient result times, solver logs, restart state, PNG frames, and postprocess outputs.

The default package preserves only each available \`processor*/constant/polyMesh\` base mesh. \`run_case.sh\` uses \`decomposePar -fields\` when the preserved processor count matches \`NPROCS\`; otherwise it removes the stale decomposition and runs \`decomposePar -force\`. A package created with \`--rebuild-decomposition\` has no processor directories and always repartitions.

Run with environment overrides as needed:

\`NPROCS=8 SOLID_THREADS=8 FOAM_BASHRC=/path/to/bashrc ./run_case.sh\`

Dependencies: OpenFOAM v2412 \`interFoam\`, \`ccx_preCICE\`, the preCICE OpenFOAM adapter, MPI, and CalculiX.
EOF

(
    cd "$STAGE"
    find . -type f ! -name MANIFEST.sha256 -print0 | sort -z | xargs -0 sha256sum > MANIFEST.sha256
)
ARCHIVE_TMP="$OUTPUT.tmp.$$.zip"
(cd "$TMP" && zip -qr "$ARCHIVE_TMP" "$NAME")
mv -f "$ARCHIVE_TMP" "$OUTPUT"
ARCHIVE_TMP=""
trap - EXIT
rm -rf "$TMP"
echo "[PRECICE-PACK] archive=$OUTPUT"
echo "[PRECICE-PACK] done; source case was not modified"
