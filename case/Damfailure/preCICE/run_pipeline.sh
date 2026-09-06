#!/usr/bin/env bash
# ============================================================================ #
#  Damfailure preCICE 流水线驱动脚本
#
#  驱动顺序: 01 网格生成 -> 03 组装算例 -> 04 preCICE 双求解器运行 -> 05 视频导出
#            -> 06 位移后处理
#            (preCICE 无独立 02 模态脚本; run_modal_analysis.py 是可选验证)
#
#  支持: 全流程 / 部分流程 (--steps)、命令行参数控制、无参数交互询问、--dry-run
# ============================================================================ #
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --------------------------------------------------------------------------- #
#  默认值
# --------------------------------------------------------------------------- #
CASE_ROOT="/root/Workspace/precice_run"
STEPS="1,3,4,5,6"
RES=115
Z_LAYERS=12
END_TIME=0.5          # 默认完整; 交互/--end-time 可改成短算例 0.02
DELTA_T=0.0005        # = precice time-window-size (必须)
WRITE_INTERVAL=0.01
NPROCS=8
SOLID_THREADS=8
GATE_E=1.0e6          # gate 材料参数, 默认与 FDM 一致 (E=1e6, nu=0, rho=2500)
GATE_NU=0.0
GATE_RHO=2500.0
VIDEO=""              # "" 未定(交互询问) / "yes" / "no"
DRY_RUN=0
SETTINGS_DIR=""
STEP_04_FAILED=0
STEP_04_RC=0

usage() {
    cat <<EOF
用法: $0 [选项]

驱动 Damfailure preCICE 五脚本流水线:
  01 网格生成 -> 03 组装算例 -> 04 preCICE 运行 -> 05 重建/后处理 -> 06 视频导出

选项 (全部可选; 不带任何参数时交互询问):
  --case-root <dir>      算例根目录 (默认 $CASE_ROOT)
  --steps <list>         运行步骤, 逗号分隔 1/3/4/5/6 (默认 $STEPS)
  --res <N>              网格 x 向密度 (传给 01, 默认 $RES)
  --z-layers <N>         z 向层数 (传给 01, 默认 $Z_LAYERS)
  --end-time <T>         仿真时长 (传给 03/04/05, 默认 $END_TIME)
  --delta-t <T>          流体 deltaT, 必须等于 preCICE time-window-size (默认 $DELTA_T)
  --write-interval <T>   输出间隔 (runTime, 传给 03, 默认 $WRITE_INTERVAL)
  --nprocs <N>           流体并行进程数 (传给 03/04, 默认 $NPROCS)
  --solid-threads <N>    固体 OMP 线程数 (传给 04, 默认 $SOLID_THREADS)
  --gate-e <E>           门材料杨氏模量 Pa (非默认时透传给 03, 默认 $GATE_E)
  --gate-nu <V>          门材料泊松比 (非默认时透传给 03, 默认 $GATE_NU)
  --gate-rho <R>         门材料密度 kg/m3 (非默认时透传给 03, 默认 $GATE_RHO)
  --video / --no-video   步骤 6 MP4 行为 (--no-video 仍保留 PNG)
  --settings-dir <dir>   participant-aware settings overlay (传给 03)
  --dry-run              只打印将执行的命令, 不执行
  -h, --help             显示本帮助

示例:
  $0 --dry-run --video                                # 预览全流程 + 视频
  $0 --steps 1,3 --res 30 --z-layers 4               # 只跑网格+组装(小网格)
  $0 --steps 3,4 --case-root /root/Workspace/x       # 只组装+运行
  $0 --steps 6 --case-root /root/Workspace/x         # 只做位移后处理 (需先跑 3/4)
  $0 --end-time 0.02                                  # 短算例全流程
EOF
}

warn()  { echo "[PRECICE WARN] $*" >&2; }
error() { echo "[PRECICE ERR] $*" >&2; }

# --------------------------------------------------------------------------- #
#  参数解析
# --------------------------------------------------------------------------- #
INTERACTIVE=0
[ $# -eq 0 ] && INTERACTIVE=1

while [ $# -gt 0 ]; do
    case "$1" in
        --case-root)       CASE_ROOT="${2:?--case-root 需要一个参数}"; shift 2 ;;
        --steps)           STEPS="${2:?--steps 需要一个参数}"; shift 2 ;;
        --res)             RES="${2:?--res 需要一个参数}"; shift 2 ;;
        --z-layers)        Z_LAYERS="${2:?--z-layers 需要一个参数}"; shift 2 ;;
        --end-time)        END_TIME="${2:?--end-time 需要一个参数}"; shift 2 ;;
        --delta-t)         DELTA_T="${2:?--delta-t 需要一个参数}"; shift 2 ;;
        --write-interval)  WRITE_INTERVAL="${2:?--write-interval 需要一个参数}"; shift 2 ;;
        --nprocs)          NPROCS="${2:?--nprocs 需要一个参数}"; shift 2 ;;
        --solid-threads)   SOLID_THREADS="${2:?--solid-threads 需要一个参数}"; shift 2 ;;
        --gate-e)          GATE_E="${2:?--gate-e 需要一个参数}"; shift 2 ;;
        --gate-nu)         GATE_NU="${2:?--gate-nu 需要一个参数}"; shift 2 ;;
        --gate-rho)        GATE_RHO="${2:?--gate-rho 需要一个参数}"; shift 2 ;;
        --video)           VIDEO=yes; shift ;;
        --no-video)        VIDEO=no; shift ;;
        --settings-dir)    SETTINGS_DIR="${2:?--settings-dir 需要一个参数}"; shift 2 ;;
        --dry-run)         DRY_RUN=1; shift ;;
        -h|--help)         usage; exit 0 ;;
        *)                 error "未知选项: $1"; usage >&2; exit 2 ;;
    esac
done

# --------------------------------------------------------------------------- #
#  无参数: 交互式询问 (回车 = 用默认值)
# --------------------------------------------------------------------------- #
if [ "$INTERACTIVE" -eq 1 ]; then
    echo "[PRECICE-PIPE] 交互模式 (回车使用默认值)"
    read -p "算例根目录 --case-root (默认 $CASE_ROOT): " inp
    if [ -n "${inp:-}" ]; then CASE_ROOT="$inp"; fi

    read -p "网格密度 --res (默认 $RES): " inp
    if [ -n "${inp:-}" ]; then RES="$inp"; fi

    read -p "z 层数 --z-layers (默认 $Z_LAYERS): " inp
    if [ -n "${inp:-}" ]; then Z_LAYERS="$inp"; fi

    read -p "仿真时长 --end-time (默认 $END_TIME): " inp
    if [ -n "${inp:-}" ]; then END_TIME="$inp"; fi

    echo "门材料参数 (回车用默认):"
    read -p "  --gate-e  杨氏模量 Pa (默认 $GATE_E): " inp
    if [ -n "${inp:-}" ]; then GATE_E="$inp"; fi
    read -p "  --gate-nu 泊松比 (默认 $GATE_NU): " inp
    if [ -n "${inp:-}" ]; then GATE_NU="$inp"; fi
    read -p "  --gate-rho 密度 kg/m3 (默认 $GATE_RHO): " inp
    if [ -n "${inp:-}" ]; then GATE_RHO="$inp"; fi

    read -p "并行进程数 --nprocs (默认 $NPROCS): " inp
    if [ -n "${inp:-}" ]; then NPROCS="$inp"; fi

    read -p "导出视频? [y/N]: " inp
    case "${inp:-}" in
        y|Y|yes|Yes|YES) VIDEO=yes ;;
        *)               VIDEO="" ;;   # 未定: 由 steps 是否含 6 决定
    esac
    echo
fi

# 非交互且未显式指定 --video/--no-video 时, VIDEO 保持未定 ("随步骤"):
# 是否合成视频由步骤列表是否含 6 决定 (含 6 自动导出; 不含 6 不导出)

# --------------------------------------------------------------------------- #
#  OpenFOAM 环境 (非强制; 01/04 内部会自行 source bashrc)
# --------------------------------------------------------------------------- #
if [ -f /usr/lib/openfoam/openfoam2412/etc/bashrc ]; then
    # 注意: OpenFOAM bashrc 在 set -e/-u 下 source 可能提前退出, 故临时关闭严格模式
    set +euo pipefail
    # shellcheck disable=SC1091
    source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1
    set -euo pipefail
fi

# Python 子脚本 stdout 改为行缓冲 (否则管道/文件重定向下默认块缓冲 ~8KB,
# run_script 实时过滤会攒批延迟; 04 运行器 flush=True 本身不受影响)
export PYTHONUNBUFFERED=1

# --------------------------------------------------------------------------- #
#  步骤校验 (只接受 1,3,4,5,6)
# --------------------------------------------------------------------------- #
validate_steps() {
    local s
    for s in $(echo "$STEPS" | tr ',' ' '); do
        case "$s" in
            1|3|4|5|6) ;;
            *) error "非法步骤 '$s' (只接受 1,3,4,5,6, 逗号分隔; preCICE 无 02 模态脚本)"; return 1 ;;
        esac
    done
    return 0
}
validate_steps || exit 2

# --------------------------------------------------------------------------- #
#  步骤间依赖检查 (仅警告, 不中断; 下游脚本自身有完整校验)
# --------------------------------------------------------------------------- #
check_deps() {
    local step="$1"
    case "$step" in
        3)
            [ -d "$CASE_ROOT/mesh/fluid-openfoam/constant/polyMesh" ] || \
                warn "步骤 3 需要网格 $CASE_ROOT/mesh/fluid-openfoam/constant/polyMesh (先跑步骤 1)"
            ;;
        4)
            [ -f "$CASE_ROOT/pc_case/precice-config.xml" ] || \
                warn "步骤 4 需要组装算例 $CASE_ROOT/pc_case/precice-config.xml (先跑步骤 3)"
            ;;
        5)
            [ -d "$CASE_ROOT/pc_case/fluid-openfoam" ] || \
                warn "步骤 5 需要组装算例 $CASE_ROOT/pc_case/fluid-openfoam (先跑步骤 3/4)"
            ;;
        6)
            [ -d "$CASE_ROOT/pc_case/fluid-openfoam" ] || \
                warn "步骤 6 需要 preCICE 算例 $CASE_ROOT/pc_case/fluid-openfoam (先跑步骤 3/4)"
            ;;
    esac
}

# --------------------------------------------------------------------------- #
#  门材料参数透传 (仅非默认值追加, 默认与 03 内部默认一致)
# --------------------------------------------------------------------------- #
GATE_ARGS=""
[ "$GATE_E"  != "1.0e6" ]  && GATE_ARGS+=" --gate-e $GATE_E"
[ "$GATE_NU" != "0.0" ]    && GATE_ARGS+=" --gate-nu $GATE_NU"
[ "$GATE_RHO" != "2500.0" ] && GATE_ARGS+=" --gate-rho $GATE_RHO"

# --------------------------------------------------------------------------- #
#  捕获 + 过滤执行器 (步骤 1/3/5/6 及视频合成)
#  子脚本完整 stdout/stderr 落盘到 $1, 控制台只实时显示带 [FDM/[PRECICE 前缀的
#  主要信息行 (log/warn/error); 工具原始输出 (blockMesh/checkMesh/ccx/pvpython/
#  foamToVTK/ffmpeg 等无前缀行) 仅保留在日志, 不透传到控制台.
#  用 awk 而非 grep 过滤: grep 无匹配时 exit 1, 会触发 set -e/pipefail 误判
#  子脚本失败; awk 始终退出 0, 子脚本真实退出码经 PIPESTATUS[0] 取回.
# --------------------------------------------------------------------------- #
run_script() {
    local logf="$1"; shift
    local rc
    mkdir -p "$(dirname "$logf")"
    eval "$*" 2>&1 | tee "$logf" | awk '/^\[(FDM|PRECICE)/ {print; fflush()}'
    rc=${PIPESTATUS[0]}
    return "$rc"
}

# --------------------------------------------------------------------------- #
#  视频行为辅助
# --------------------------------------------------------------------------- #
steps_contains() {
    echo ",$STEPS," | grep -q ",$1,"
}

# --------------------------------------------------------------------------- #
#  执行单个步骤
# --------------------------------------------------------------------------- #
run_step() {
    local step="$1"
    local cmd=""
    case "$step" in
        1) cmd="python3 '$SCRIPT_DIR/01_generate_mesh.py' --res $RES --z-layers $Z_LAYERS --out-dir '$CASE_ROOT/mesh' --view --render" ;;
        3) cmd="python3 '$SCRIPT_DIR/03_assemble_precice_case.py' --mesh-dir '$CASE_ROOT/mesh/fluid-openfoam' --solid-dir '$CASE_ROOT/mesh/solid-calculix' --out-dir '$CASE_ROOT/pc_case' --end-time $END_TIME --delta-t $DELTA_T --write-interval $WRITE_INTERVAL --nprocs $NPROCS"$GATE_ARGS ;;
        4) cmd="python3 '$SCRIPT_DIR/04_run_precice.py' -f --case '$CASE_ROOT/pc_case' --end-time $END_TIME --delta-t $DELTA_T --write-interval $WRITE_INTERVAL --nprocs $NPROCS --solid-threads $SOLID_THREADS" ;;
        5) cmd="python3 '$SCRIPT_DIR/05_postprocess.py' --case '$CASE_ROOT/pc_case' --out-dir '$CASE_ROOT/postprocess'" ;;
        6) cmd="xvfb-run -a pvpython '$SCRIPT_DIR/06_export_precice_video.py' --case '$CASE_ROOT/pc_case' --out '$CASE_ROOT/frames' --mp4-path '$CASE_ROOT/precice_result.mp4'" ;;
    esac

    if [ "$step" -eq 3 ] && [ -n "$SETTINGS_DIR" ]; then
        cmd+=" --settings-dir '$SETTINGS_DIR'"
    fi
    if [ "$step" -eq 6 ] && [ "$VIDEO" = "no" ]; then
        cmd+=" --no-mp4"
    fi
    echo "[PRECICE-PIPE] 步骤 $step: $cmd"
    [ "$DRY_RUN" -eq 1 ] && return 0

    check_deps "$step"
    if [ "$step" -eq 4 ]; then
        # 步骤 4 (运行器): 自带 \r 进度条单行实时刷新, 捕获过滤会破坏刷新, 直接透传
        # (输出已较干净: [PRECICE 04] 前缀行 + 进度条 + solver Time=/ExecutionTime= 行)
        if ! eval "$cmd"; then
            error "步骤 4 执行失败, 返回流水线以继续请求的下游恢复步骤"
            return 1
        fi
    elif [ "$step" -eq 6 ]; then
        # 步骤 6: pvpython 在 xvfb 下渲染完成后可能抛无害 GLXBadContext 导致非零退出码,
        # 但只要帧文件已生成即视为成功。捕获过滤: 只透传 [PRECICE 前缀行, 完整输出落盘。
        if ! run_script "$CASE_ROOT/run_step_6.log" "$cmd"; then
            if ls "$CASE_ROOT"/frames/frame_*.png >/dev/null 2>&1 && { [ "$VIDEO" = "no" ] || [ -s "$CASE_ROOT/precice_result.mp4" ]; }; then
                warn "步骤 6: pvpython 非零但请求的产物完整, 视为 GLX 清理错误"
            else
                error "步骤 6 执行失败: PNG 或 MP4 产物不完整"
                return 1
            fi
        fi
    else
        # 步骤 1/3/5: 捕获过滤, 只透传带 [FDM/[PRECICE 前缀的主要行, 工具原始输出仅落盘
        if ! run_script "$CASE_ROOT/run_step_$step.log" "$cmd"; then
            error "步骤 $step 执行失败, 流水线中止"
            return 1
        fi
    fi
}

# --------------------------------------------------------------------------- #
#  打印配置摘要
# --------------------------------------------------------------------------- #
summary_line() {
    local label="$1"
    local value="$2"
    local padding

    case "$label" in
        算例根目录) padding=1 ;;
        门材料)     padding=5 ;;
        步骤|网格|仿真|并行|视频) padding=7 ;;
        *)          padding=1 ;;
    esac
    printf '[PRECICE-PIPE] %s%*s: %s\n' "$label" "$padding" '' "$value"
}

echo "[PRECICE-PIPE] ===== Damfailure preCICE 流水线 (步骤 $STEPS) ====="
summary_line "算例根目录" "$CASE_ROOT"
summary_line "步骤" "$STEPS"
summary_line "网格" "--res $RES --z-layers $Z_LAYERS"
summary_line "仿真" "endTime=$END_TIME deltaT=$DELTA_T writeInterval=$WRITE_INTERVAL"
summary_line "并行" "$NPROCS 流体进程 / $SOLID_THREADS 固体线程"
summary_line "门材料" "E=$GATE_E nu=$GATE_NU rho=$GATE_RHO${GATE_ARGS:+ (透传: ${GATE_ARGS# }) }"
if [ "$VIDEO" = "yes" ]; then
    VIDEO_STATE="强制是 (--video)"
elif [ "$VIDEO" = "no" ]; then
    VIDEO_STATE="强制否 (--no-video)"
elif steps_contains 6; then
    VIDEO_STATE="步骤6默认合成 MP4 (30 fps)"
else
    VIDEO_STATE="不含步骤6"
fi
summary_line "视频" "$VIDEO_STATE"
[ "$DRY_RUN" -eq 1 ] && echo "[PRECICE-PIPE] DRY-RUN    : 只打印命令, 不执行"
echo

# --------------------------------------------------------------------------- #
#  按步骤顺序执行 (去重 + 数值排序, 保证 1 < 3 < 4 < 5 < 6)
# --------------------------------------------------------------------------- #
for s in $(echo "$STEPS" | tr ',' '\n' | sort -u -n); do
    if [ "$s" -eq 4 ]; then
        set +e
        run_step "$s"
        STEP_04_RC=$?
        set -e
        if [ "$STEP_04_RC" -ne 0 ]; then
            STEP_04_FAILED=1
            error "步骤 4 执行失败 (exit $STEP_04_RC), 继续请求的步骤 5/6"
        fi
    else
        run_step "$s"
    fi
done

echo
echo "[PRECICE-PIPE] 流水线完成"
if [ "$STEP_04_FAILED" -eq 1 ]; then
    error "流水线完成但步骤 4 失败 (exit $STEP_04_RC)"
    exit 1
fi
