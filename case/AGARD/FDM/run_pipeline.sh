#!/usr/bin/env bash
# ============================================================================ #
#  AGARD 445.6 自包含流水线驱动
#
#  阶段:
#    01 解压网格(data/fluid8_mesh.zip)  -> constant/polyMesh
#    02 解压模态(data/modes_4m.zip)     -> mode/ + 节点数门禁
#    03 稳态初场(运行时算)              -> 0/
#    04 组装算例(渲染 myHisa 字典)
#    05 颤振运行(分解 + myHisa)
#    06 后处理(sigma / f / 振幅)
#    07 出图
#
#  默认档位 = 冒烟: 20 步, 8 核, 约 1.5-2 h (含稳态)
# ============================================================================ #
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- environment bootstrap ------------------------------------------------- #
# The solver and the library live in the user tree ($HOME/OpenFOAM/user-v2412),
# which etc/bashrc does not select on its own.
if ! command -v myHisa >/dev/null 2>&1; then
    cat >&2 <<'EOF'
ERROR: myHisa is not on PATH. Load the OpenFOAM environment first, e.g.

  source /usr/lib/openfoam/openfoam2412/etc/bashrc
  . $WM_PROJECT_DIR/bin/tools/change-userdir.sh $HOME/OpenFOAM/user-v2412 -platforms

(the second line selects the user tree that holds myHisa and
libfastDynamicFvMesh.so; etc/bashrc alone points at $HOME/OpenFOAM/$USER-v2412)
EOF
    exit 1
fi
STEPS="1,2,3,4,5,6"
CASE_ROOT=""
MACH="0.96"
PRESSURE="2936"
MODES="2"
AMR="--amr"
END_TIME="0.02"
NPROCS="8"
DRY_RUN=0
SKIP_IF_FRESH=0
FORCE=0

usage() {
    cat <<EOF
Usage: $0 [options]

Stages:
  1  解压网格 (data/fluid8_mesh.zip -> constant/polyMesh)
  2  解压模态 (data/modes_4m.zip -> mode/), 节点数门禁
  3  稳态初场 (运行时算) -> 0/
  4  组装算例 (渲染 myHisa 字典)
  5  颤振运行 (decomposePar + myHisa -parallel)
  6  后处理 (逐周期 p2p, sigma, f)

Options:
  --steps <list>        运行阶段, 逗号分隔 (默认 "$STEPS")
  --case-root <dir>     算例目录 (默认 /root/Workspace/AGARD/FDM-run/M<mach>_q<q>)
  --mach <M>            马赫数 (默认 $MACH)
  --pressure <q>        来流动压 [Pa] (默认 $PRESSURE)
  --modes <n>           叠加模态数, 写成 nActiveModes (默认 $MODES; mode/ 固定 4 阶)
  --no-amr              关闭 AMR
  --end-time <t>        结束时间 (默认 $END_TIME = 20 步; 完整 4 周期用 0.286)
  --nprocs <n>          MPI 进程数 (默认 $NPROCS)
  --force               覆盖已存在的网格/模态
  --dry-run             只打印将执行的命令
  -h, --help            显示本帮助

Examples:
  $0 --dry-run                                   # 预览全流程
  $0 --steps 1,2,4                               # 只生成算例, 不求解
  $0 --steps 5 --end-time 0.286                  # 单点完整 (过夜)
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --steps)        STEPS="${2:?--steps 需要一个参数}"; shift 2 ;;
        --case-root)    CASE_ROOT="${2:?--case-root 需要一个参数}"; shift 2 ;;
        --mach)         MACH="${2:?}"; shift 2 ;;
        --pressure)     PRESSURE="${2:?}"; shift 2 ;;
        --modes)        MODES="${2:?}"; shift 2 ;;
        --no-amr)       AMR="--no-amr"; shift ;;
        --end-time)     END_TIME="${2:?}"; shift 2 ;;
        --nprocs)       NPROCS="${2:?}"; shift 2 ;;
        --force)        FORCE=1; shift ;;
        --dry-run)      DRY_RUN=1; shift ;;
        --skip-if-fresh) SKIP_IF_FRESH=1; shift ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "未知参数: $1" >&2; usage; exit 2 ;;
    esac
done

if [ -z "$CASE_ROOT" ]; then
    CASE_ROOT="/root/Workspace/AGARD/FDM-run/M${MACH}_q${PRESSURE}"
fi

run() {
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [dry-run] $*"
    else
        "$@"
    fi
}

stage() { case ",$STEPS," in *",$1,"*) return 0 ;; *) return 1 ;; esac; }

mkdir -p "$CASE_ROOT"
echo "=============================================================="
echo " case root : $CASE_ROOT"
echo " M = $MACH   q = $PRESSURE Pa   modes = $MODES   AMR = ${AMR#--}"
echo " endTime = $END_TIME   nprocs = $NPROCS   steps = $STEPS"
echo "=============================================================="

if stage 1; then
    echo "--- [1] mesh ---"
    args=(--case-root "$CASE_ROOT")
    [ "$FORCE" -eq 1 ] && args+=(--force)
    run python3 "$SCRIPT_DIR/01_generate_mesh.py" "${args[@]}"
fi

if stage 2; then
    echo "--- [2] modes ---"
    args=(--case-root "$CASE_ROOT" --n-modes 4)
    [ "$FORCE" -eq 1 ] && args+=(--force)
    run python3 "$SCRIPT_DIR/02_prepare_modes.py" "${args[@]}"
fi

if stage 3; then
    echo "--- [3] steady initial field ---"
    run python3 "$SCRIPT_DIR/03_make_steady.py" \
        --case-root "$CASE_ROOT" --mach "$MACH" --pressure "$PRESSURE" \
        --nprocs "$NPROCS"
fi

if stage 4; then
    echo "--- [4] assemble case ---"
    run python3 "$SCRIPT_DIR/04_assemble_case.py" \
        --case-root "$CASE_ROOT" --mach "$MACH" --pressure "$PRESSURE" \
        --modes "$MODES" --end-time "$END_TIME" --nprocs "$NPROCS" "$AMR"
fi

if stage 5; then
    echo "--- [5] flutter run ---"
    run python3 "$SCRIPT_DIR/05_run_flutter.py" \
        --case-root "$CASE_ROOT" --nprocs "$NPROCS" --end-time "$END_TIME"
fi

if stage 6; then
    echo "--- [6] postprocess ---"
    run python3 "$SCRIPT_DIR/06_postprocess.py" --case-root "$CASE_ROOT"
fi

echo "=============================================================="
echo " done: $CASE_ROOT"
