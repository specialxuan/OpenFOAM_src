#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AGARD_ROOT="/root/Workspace/AGARD"
WORK_ROOT="$AGARD_ROOT/Mode-run"
FARFIELD_DECK="$SCRIPT_DIR/mesh/farfield_c3d8.zip"
WING_DECK="$SCRIPT_DIR/mesh/wing_c3d10.zip"
STEPS="1,2"
RESULTS_DIR=""
RESULTS_STEM=""
MODES=10
TIE_TOLERANCE=0.05
WING_TYPE="soft"
FAR_FIELD_MODE="layered"
LAYERS=20
LAYER_EMAX_FACTOR="1e-6"
LAYER_REF_MODULUS="E1"
LAYER_P="0.25"
LAYER_FLOOR="1000"
FAR_FIELD_MODULUS="28.0"
PURE_WING=0
PURE_DAT=""
PURE_FRD=""
MAX_REL=""
MIN_COSINE=""
RUN_CCX=0
CCX_COMMAND="ccx"
DRY_RUN=0

usage() {
    cat <<'EOF'
Usage: run_pipeline.sh [options]

Stages:
  1  assemble the merged mesh deck (far-field C3D8 + wing C3D10)
  2  generate the combined modal deck (optionally --run-ccx)
  3  verify the modal results, then export the FDM CSVs (gate first: no CSV
     is written when a gate fails)

Default stages are 1,2 (no solver launched). Stage 3 needs a solved deck:
run `--steps 1,2,3 --run-ccx`, or point --results-dir at a directory that
already holds <stem>.inp/.dat/.frd.

Pure-wing reference (needed by the stage-3 frequency/identity gates):
  --pure-wing              stage 2 also writes modal/pure_<type>.inp from the
                           wing source deck; with --run-ccx it is also solved,
                           producing modal/pure_<type>.dat/.frd. Stage 3 then
                           defaults --pure-dat/--pure-frd to those files.
  Alternatively pass --pure-dat/--pure-frd for an existing solved pure deck.

Inputs (defaults are the committed meshes next to this script, so no
/root/Workspace checkout is needed; plain, .gz and single-member .zip decks
are all accepted):
  --farfield-deck FILE     far-field C3D8 deck
                           (default <script dir>/mesh/farfield_c3d8.zip)
  --wing-deck FILE         second-order C3D10 wing deck
                           (default <script dir>/mesh/wing_c3d10.zip)
  --work-root DIR          generated-artifact root (default /root/Workspace/AGARD/Mode-run)
  --results-dir DIR        pre-solved deck directory for stage 3
  --results-stem STEM      deck stem inside --results-dir
                           (default combined_<wing-type>_<far-field-mode>)

Controls:
  --steps LIST             comma-separated stages from 1,2,3 (default 1,2)
  --modes N                frequency count (default 10)
  --tie-tolerance X        FSI tie tolerance and FSI-to-SIDE distance gate
                           (default 0.05)
  --wing-type TYPE         wing material: soft (soft_anchored) or
                           hard (hard_e1_0p75) (default soft)
  --far-field-mode MODE    layered (bounded band, default) or uniform
  --layers N               volume-ranked far-field layers (default 20)
  --layer-emax-factor X    E_max = X * reference modulus (default 1e-6)
  --layer-ref-modulus M    E1, E2 or E3 (default E1)
  --layer-p X              E(V)=E_max*(V_min/V)^p exponent (default 0.25)
  --layer-floor X          E_min = E_max / X, the span bound (default 1000)
  --far-field-modulus X    uniform-mode far-field E in Pa (default 28)

Solver / stage 3 pass-through:
  --run-ccx                run CalculiX on the combined deck in stage 2
  --ccx-command CMD        CalculiX executable (default ccx)
  --pure-dat FILE          pure-wing DAT reference for the frequency gate
  --pure-frd FILE          pure-wing FRD reference for the mode-identity gate
  --max-rel X              max |relative frequency difference| (default 0.01)
  --min-cosine X           min mode field cosine (default 0.9)
  --dry-run                print every stage command; write nothing
EOF
}

error() { echo "[AGARD ERR] $*" >&2; }
is_uint() { case "$1" in ''|*[!0-9]*) return 1 ;; *) return 0 ;; esac; }
is_positive() { awk -v value="$1" 'BEGIN { exit !(value + 0 > 0) }'; }

while [ "$#" -gt 0 ]; do
    case "$1" in
        --work-root) WORK_ROOT="${2:?--work-root requires a directory}"; shift 2 ;;
        --farfield-deck) FARFIELD_DECK="${2:?--farfield-deck requires a file}"; shift 2 ;;
        --wing-deck) WING_DECK="${2:?--wing-deck requires a file}"; shift 2 ;;
        --steps) STEPS="${2:?--steps requires a list}"; shift 2 ;;
        --results-dir) RESULTS_DIR="${2:?--results-dir requires a directory}"; shift 2 ;;
        --results-stem) RESULTS_STEM="${2:?--results-stem requires a stem}"; shift 2 ;;
        --modes) MODES="${2:?--modes requires an integer}"; shift 2 ;;
        --tie-tolerance) TIE_TOLERANCE="${2:?--tie-tolerance requires a number}"; shift 2 ;;
        --wing-type) WING_TYPE="${2:?--wing-type requires hard or soft}"; shift 2 ;;
        --far-field-mode) FAR_FIELD_MODE="${2:?--far-field-mode requires layered or uniform}"; shift 2 ;;
        --layers) LAYERS="${2:?--layers requires an integer}"; shift 2 ;;
        --layer-emax-factor) LAYER_EMAX_FACTOR="${2:?--layer-emax-factor requires a number}"; shift 2 ;;
        --layer-ref-modulus) LAYER_REF_MODULUS="${2:?--layer-ref-modulus requires E1, E2 or E3}"; shift 2 ;;
        --layer-p) LAYER_P="${2:?--layer-p requires a number}"; shift 2 ;;
        --layer-floor) LAYER_FLOOR="${2:?--layer-floor requires a number}"; shift 2 ;;
        --far-field-modulus) FAR_FIELD_MODULUS="${2:?--far-field-modulus requires a number}"; shift 2 ;;
        --pure-wing) PURE_WING=1; shift ;;
        --pure-dat) PURE_DAT="${2:?--pure-dat requires a file}"; shift 2 ;;
        --pure-frd) PURE_FRD="${2:?--pure-frd requires a file}"; shift 2 ;;
        --max-rel) MAX_REL="${2:?--max-rel requires a number}"; shift 2 ;;
        --min-cosine) MIN_COSINE="${2:?--min-cosine requires a number}"; shift 2 ;;
        --run-ccx) RUN_CCX=1; shift ;;
        --ccx-command) CCX_COMMAND="${2:?--ccx-command requires an executable}"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) error "unknown option: $1"; usage >&2; exit 2 ;;
    esac
done

NORMALIZED_STEPS="$(printf '%s' "$STEPS" | tr ',' '\n' | sort -nu)"
for step in $NORMALIZED_STEPS; do
    case "$step" in 1|2|3) ;; *) error "invalid stage '$step' (allowed: 1,2,3)"; exit 2 ;; esac
done

is_uint "$MODES" && [ "$MODES" -gt 0 ] || { error "--modes must be a positive integer"; exit 2; }
is_uint "$LAYERS" && [ "$LAYERS" -gt 0 ] || { error "--layers must be a positive integer"; exit 2; }
is_positive "$TIE_TOLERANCE" || { error "--tie-tolerance must be positive"; exit 2; }
is_positive "$LAYER_EMAX_FACTOR" || { error "--layer-emax-factor must be positive"; exit 2; }
is_positive "$LAYER_P" || { error "--layer-p must be positive"; exit 2; }
is_positive "$LAYER_FLOOR" && awk -v value="$LAYER_FLOOR" 'BEGIN { exit !(value > 1) }' || { error "--layer-floor must be greater than 1"; exit 2; }
is_positive "$FAR_FIELD_MODULUS" || { error "--far-field-modulus must be positive"; exit 2; }
case "$WING_TYPE" in hard|soft) ;; *) error "--wing-type must be 'hard' or 'soft'"; exit 2 ;; esac
case "$FAR_FIELD_MODE" in layered|uniform) ;; *) error "--far-field-mode must be 'layered' or 'uniform'"; exit 2 ;; esac
case "$LAYER_REF_MODULUS" in E1|E2|E3) ;; *) error "--layer-ref-modulus must be E1, E2 or E3"; exit 2 ;; esac
if [ -n "$MAX_REL" ]; then is_positive "$MAX_REL" || { error "--max-rel must be positive"; exit 2; }; fi
if [ -n "$MIN_COSINE" ]; then
    is_positive "$MIN_COSINE" && awk -v value="$MIN_COSINE" 'BEGIN { exit !(value <= 1) }' || { error "--min-cosine must be in (0, 1]"; exit 2; }
fi

has_step() { printf '%s\n' "$NORMALIZED_STEPS" | grep -qx "$1"; }
[ "$RUN_CCX" -eq 0 ] || has_step 2 || { error "--run-ccx requires stage 2"; exit 2; }
[ "$PURE_WING" -eq 0 ] || has_step 2 || { error "--pure-wing requires stage 2"; exit 2; }
if has_step 3 && [ -z "$RESULTS_DIR" ] && { ! has_step 2 || [ "$RUN_CCX" -eq 0 ]; }; then
    error "stage 3 without --results-dir requires stages 2,3 and --run-ccx"
    exit 2
fi

STEM="combined_${WING_TYPE}_${FAR_FIELD_MODE}"
[ -n "$RESULTS_STEM" ] || RESULTS_STEM="$STEM"
MESH_DECK="$WORK_ROOT/mesh/merged_base.inp"
MESH_META="$WORK_ROOT/mesh/merged_meta.json"
MODAL_DIR="$WORK_ROOT/modal"
MODAL_DECK="$MODAL_DIR/$STEM.inp"
MODAL_META="$MODAL_DIR/${STEM}_meta.json"
PURE_DECK="$MODAL_DIR/pure_${WING_TYPE}.inp"
if [ -n "$RESULTS_DIR" ]; then
    STAGE3_DECK="$RESULTS_DIR/$RESULTS_STEM.inp"
    STAGE3_DAT="$RESULTS_DIR/$RESULTS_STEM.dat"
    STAGE3_FRD="$RESULTS_DIR/$RESULTS_STEM.frd"
else
    STAGE3_DECK="$MODAL_DECK"
    STAGE3_DAT="$MODAL_DIR/$STEM.dat"
    STAGE3_FRD="$MODAL_DIR/$STEM.frd"
fi
if [ "$PURE_WING" -eq 1 ]; then
    [ -n "$PURE_DAT" ] || PURE_DAT="$MODAL_DIR/pure_${WING_TYPE}.dat"
    [ -n "$PURE_FRD" ] || PURE_FRD="$MODAL_DIR/pure_${WING_TYPE}.frd"
fi

run_step() {
    local step="$1" command
    case "$step" in
        1)
            command="python3 '$SCRIPT_DIR/01_assemble_mesh.py' --farfield-deck '$FARFIELD_DECK' --wing-deck '$WING_DECK' --output '$MESH_DECK' --meta '$MESH_META' --tie-tolerance $TIE_TOLERANCE"
            ;;
        2)
            command="python3 '$SCRIPT_DIR/02_generate_modal_case.py' --input '$MESH_DECK' --output '$MODAL_DECK' --meta '$MODAL_META' --modes $MODES --tie-tolerance $TIE_TOLERANCE --wing-type $WING_TYPE --far-field-mode $FAR_FIELD_MODE"
            if [ "$FAR_FIELD_MODE" = "uniform" ]; then
                command="$command --far-field-modulus $FAR_FIELD_MODULUS"
            else
                command="$command --layers $LAYERS --layer-emax-factor $LAYER_EMAX_FACTOR --layer-ref-modulus $LAYER_REF_MODULUS --layer-p $LAYER_P --layer-floor $LAYER_FLOOR"
            fi
            if [ "$PURE_WING" -eq 1 ]; then
                command="$command --pure-wing --pure-output '$PURE_DECK' --wing-source '$WING_DECK'"
                [ "$RUN_CCX" -eq 1 ] && command="$command --run-pure-ccx"
            fi
            [ "$RUN_CCX" -eq 1 ] && command="$command --run-ccx --ccx-command '$CCX_COMMAND'"
            ;;
        3)
            command="python3 '$SCRIPT_DIR/03_export_modes.py' --input '$STAGE3_DECK' --dat '$STAGE3_DAT' --frd '$STAGE3_FRD' --output-dir '$WORK_ROOT/mode' --modes $MODES"
            if [ -n "$PURE_DAT" ]; then command="$command --pure-dat '$PURE_DAT'"; fi
            if [ -n "$PURE_FRD" ]; then command="$command --pure-frd '$PURE_FRD'"; fi
            if [ -n "$MAX_REL" ]; then command="$command --max-rel $MAX_REL"; fi
            if [ -n "$MIN_COSINE" ]; then command="$command --min-cosine $MIN_COSINE"; fi
            ;;
    esac
    echo "[AGARD-PIPE] stage $step: $command"
    [ "$DRY_RUN" -eq 1 ] && return 0
    mkdir -p "$WORK_ROOT/logs"
    eval "$command" 2>&1 | tee "$WORK_ROOT/logs/stage_${step}.log"
}

echo "[AGARD-PIPE] work-root=$WORK_ROOT steps=$STEPS modes=$MODES tie-tolerance=$TIE_TOLERANCE wing-type=$WING_TYPE far-field=$FAR_FIELD_MODE layers=$LAYERS layer-emax-factor=$LAYER_EMAX_FACTOR layer-ref=$LAYER_REF_MODULUS layer-p=$LAYER_P layer-floor=$LAYER_FLOOR pure-wing=$PURE_WING run-ccx=$RUN_CCX dry-run=$DRY_RUN"
for step in $NORMALIZED_STEPS; do run_step "$step"; done
echo "[AGARD-PIPE] complete"
