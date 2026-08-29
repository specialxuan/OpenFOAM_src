#!/usr/bin/env bash
set -euo pipefail

CASE_ROOT="/root/Workspace/fdm_amr_3"
OUTPUT=""
FORCE=0
DRY_RUN=0
KEEP_DECOMPOSITION=1

usage() {
    cat <<EOF
用法: $0 [选项]

在临时目录中清理并打包可重复运行的 FDM 算例；原始算例不会被修改。

  --case-root <dir>        算例根目录 (默认 $CASE_ROOT)
  --output <file.zip>      输出压缩包 (默认 <case-root>_clean.zip)
  --rebuild-decomposition  不保留 processor*/constant/polyMesh，重跑时重新分区
  --force                  覆盖已存在的压缩包
  --dry-run                只检查并打印计划
  -h, --help               显示帮助
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --case-root) CASE_ROOT="${2:?--case-root 需要参数}"; shift 2 ;;
        --output) OUTPUT="${2:?--output 需要参数}"; shift 2 ;;
        --rebuild-decomposition) KEEP_DECOMPOSITION=0; shift ;;
        --force) FORCE=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "[FDM-PACK ERR] 未知参数: $1" >&2; usage >&2; exit 2 ;;
    esac
done

CASE_ROOT="$(realpath -e "$CASE_ROOT")"
CASE_DIR="$CASE_ROOT/fdm_case"
NAME="$(basename "$CASE_ROOT")"
if [ -z "$OUTPUT" ]; then
    OUTPUT="$(dirname "$CASE_ROOT")/${NAME}_clean.zip"
fi
OUTPUT="$(realpath -m "$OUTPUT")"
if [[ "$OUTPUT" != *.zip ]]; then
    echo "[FDM-PACK ERR] --output 必须使用 .zip 扩展名: $OUTPUT" >&2
    exit 2
fi
if ! command -v zip >/dev/null 2>&1; then
    echo "[FDM-PACK ERR] 找不到 zip 命令" >&2
    exit 2
fi

case "$OUTPUT" in
    "$CASE_ROOT"/*) echo "[FDM-PACK ERR] 输出压缩包不能位于源算例内部" >&2; exit 2 ;;
esac
if [ -e "$OUTPUT" ] && [ "$FORCE" -ne 1 ]; then
    echo "[FDM-PACK ERR] 输出已存在，使用 --force 覆盖: $OUTPUT" >&2
    exit 2
fi

REQUIRED=(
    "$CASE_ROOT/CaseDescription.md"
    "$CASE_ROOT/fdm_result.mp4"
    "$CASE_DIR/0"
    "$CASE_DIR/case.foam"
    "$CASE_DIR/constant/polyMesh/points"
    "$CASE_DIR/constant/polyMesh/faces"
    "$CASE_DIR/constant/polyMesh/owner"
    "$CASE_DIR/constant/polyMesh/neighbour"
    "$CASE_DIR/constant/polyMesh/boundary"
    "$CASE_DIR/constant/dynamicMeshDict"
    "$CASE_DIR/system/controlDict"
    "$CASE_DIR/system/fvSchemes"
    "$CASE_DIR/system/fvSolution"
    "$CASE_DIR/system/setFieldsDict"
    "$CASE_DIR/system/decomposeParDict"
    "$CASE_DIR/mode/FluidNodeCoor.csv"
    "$CASE_DIR/mode/FluidNodeDisp1.csv"
)
for path in "${REQUIRED[@]}"; do
    if [ ! -e "$path" ]; then
        echo "[FDM-PACK ERR] 缺少必要文件: $path" >&2
        exit 2
    fi
done

echo "[FDM-PACK] source=$CASE_ROOT"
echo "[FDM-PACK] output=$OUTPUT"
echo "[FDM-PACK] preserve decomposition=$KEEP_DECOMPOSITION"
echo "[FDM-PACK] remove from staged copy: result times, logs, restart state, processor time results"
if [ "$DRY_RUN" -eq 1 ]; then
    echo "[FDM-PACK] dry-run complete; source unchanged"
    exit 0
fi

OUTPUT_DIR="$(dirname "$OUTPUT")"
mkdir -p "$OUTPUT_DIR"
TMP="$(mktemp -d "$OUTPUT_DIR/.${NAME}.pack.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
STAGE="$TMP/$NAME"
STAGE_CASE="$STAGE/fdm_case"
mkdir -p "$STAGE_CASE/constant" "$STAGE_CASE/system"

cp -a "$CASE_ROOT/CaseDescription.md" "$STAGE/"
cp -a "$CASE_ROOT/fdm_result.mp4" "$STAGE/"
if [ -d "$CASE_ROOT/case_settings" ]; then
    cp -a "$CASE_ROOT/case_settings" "$STAGE/"
fi
cp -a "$CASE_DIR/0" "$STAGE_CASE/"
cp -a "$CASE_DIR/mode" "$STAGE_CASE/"
cp -a "$CASE_DIR/system/." "$STAGE_CASE/system/"
cp -a "$CASE_DIR/case.foam" "$STAGE_CASE/"

for path in "$CASE_DIR/constant"/*; do
    base="$(basename "$path")"
    case "$base" in
        fsiRestart) continue ;;
        polyMesh) cp -a "$path" "$STAGE_CASE/constant/" ;;
        *) cp -a "$path" "$STAGE_CASE/constant/" ;;
    esac
done

if [ "$KEEP_DECOMPOSITION" -eq 1 ]; then
    for proc in "$CASE_DIR"/processor*; do
        [ -d "$proc/constant/polyMesh" ] || continue
        target="$STAGE_CASE/$(basename "$proc")/constant"
        mkdir -p "$target"
        cp -a "$proc/constant/polyMesh" "$target/"
    done
fi

cat > "$STAGE/run_case.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE="$ROOT/fdm_case"
source /usr/lib/openfoam/openfoam2412/etc/bashrc
cd "$CASE"
checkMesh
setFields
if compgen -G 'processor*/constant/polyMesh' >/dev/null; then
    decomposePar -fields
else
    decomposePar -force
fi
if [ "$(id -u)" -eq 0 ]; then
    mpirun --allow-run-as-root -np 8 myInterFoam -parallel
else
    mpirun -np 8 myInterFoam -parallel
fi
EOF
chmod +x "$STAGE/run_case.sh"

cat > "$STAGE/README.md" <<EOF
# $NAME 可重复运行包

包含：最终视频、算例说明、初始网格、初始场、模态文件、OpenFOAM 设置，以及保留的 8 核 processor 基础网格/节点 addressing。

已清理：非零时间结果、求解日志、重建结果日志、\`constant/fsiRestart\`、processor 非零时间结果、PNG 帧和后处理结果。

依赖：OpenFOAM v2412、\`myInterFoam\`、\`libfastDynamicFvMesh.so\`、MPI。

解压后运行：

\`\`\`bash
./run_case.sh
\`\`\`

\`run_case.sh\` 会执行 \`checkMesh\`、\`setFields\`，使用已保留的 processor 基础网格通过 \`decomposePar -fields\` 分配初始场，再运行 8 核 \`myInterFoam -parallel\`。

若使用 \`--rebuild-decomposition\` 打包，则不包含 processor 目录，\`run_case.sh\` 会根据 \`system/decomposeParDict\` 重新分区。
EOF

if [ "$KEEP_DECOMPOSITION" -eq 0 ]; then
    rm -rf "$STAGE_CASE"/processor*
fi

find "$STAGE_CASE" -mindepth 1 -maxdepth 1 -type d -name '[1-9]*' -exec rm -rf {} +
find "$STAGE_CASE" -maxdepth 1 -type f \( -name 'log.*' -o -name '*.log' -o -name 'modal_diagnostics.csv' \) -delete
find "$STAGE_CASE" -mindepth 2 -maxdepth 2 -type d -path '*/processor*/[1-9]*' -exec rm -rf {} +

(
    cd "$STAGE"
    find . -type f ! -name MANIFEST.sha256 -print0 \
        | sort -z \
        | xargs -0 sha256sum > MANIFEST.sha256
)

if find "$STAGE_CASE" -mindepth 1 -maxdepth 1 -type d -name '[1-9]*' | grep -q .; then
    echo "[FDM-PACK ERR] staging still contains nonzero result times" >&2
    exit 1
fi
if find "$STAGE_CASE" -maxdepth 1 -type f -name 'log.*' | grep -q .; then
    echo "[FDM-PACK ERR] staging still contains logs" >&2
    exit 1
fi

ARCHIVE_TMP="$OUTPUT.tmp.$$.zip"
(cd "$TMP" && zip -qr "$ARCHIVE_TMP" "$NAME")
mv -f "$ARCHIVE_TMP" "$OUTPUT"
trap - EXIT
rm -rf "$TMP"
echo "[FDM-PACK] archive=$OUTPUT"
echo "[FDM-PACK] done; source case was not modified"
