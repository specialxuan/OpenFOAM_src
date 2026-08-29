#!/usr/bin/env bash
# ============================================================================ #
#  Damfailure 六脚本流水线驱动脚本
#
#  驱动顺序: 01 网格生成 -> 02 模态分析 -> 03 组装算例 -> 04 FDM 运行
#            -> 05 并行重建/后处理 -> 06 视频导出
#
#  支持: 全流程 / 部分流程 (--steps)、命令行参数控制、无参数交互询问、--dry-run
# ============================================================================ #
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --------------------------------------------------------------------------- #
#  默认值
# --------------------------------------------------------------------------- #
CASE_ROOT="/root/Workspace/fdm_run"
STEPS="1,2,3,4,5,6"
RES=115
Z_LAYERS=12
END_TIME=0.5
WRITE_INTERVAL=50
DELTA_T=""        # 空值: 使用组装设置中的 deltaT
VIDEO=""          # "" 未定(交互询问) / "yes" / "no"
NPROCS=8
AMR="no"          # "yes" / "no" (默认关闭)
SETTINGS_DIR=""   # 可选: 传给 03 的预写设置目录
DRY_RUN=0
STEP_04_FAILED=0
STEP_04_RC=0

usage() {
    cat <<EOF
用法: $0 [选项]

驱动 Damfailure 六脚本流水线:
  01 网格生成 -> 02 模态分析 -> 03 组装算例 -> 04 FDM 运行
  -> 05 并行重建/后处理 -> 06 视频导出

选项 (全部可选; 不带任何参数时交互询问):
  --case-root <dir>      算例根目录 (默认 $CASE_ROOT)
  --steps <list>         运行步骤, 逗号分隔 1-6 (默认 $STEPS), 如 "3,4"; 步骤 6 默认导出 PNG+MP4
  --res <N>              网格 x 向密度 (传给 01, 默认 $RES)
  --z-layers <N>         z 向层数 (传给 01, 默认 $Z_LAYERS)
  --end-time <T>         仿真时长 (传给 03/04, 默认 $END_TIME)
  --write-interval <N>   输出间隔步数 (默认 $WRITE_INTERVAL)
  --delta-t <dt>         时间步长, 传给步骤 3/4 (默认使用设置文件)
  --video / --no-video   步骤 6 是否合成 MP4 (--no-video 仍保留 PNG 帧)
  --nprocs <N>           并行进程数 (传给 03/04, 默认 $NPROCS)
  --amr                  启用动态自适应网格 (dynamic AMR, 传给 03, 默认关闭)
  --settings-dir <dir>   预写算例设置目录, 传给步骤 03 覆盖内嵌默认配置
  --dry-run              只打印将执行的命令, 不执行
  -h, --help             显示本帮助

示例:
  $0 --dry-run --video                                # 预览全流程 + 视频
  $0 --steps 1,2,3 --res 30 --z-layers 4               # 只跑网格+模态+组装(小网格)
  $0 --steps 3,4 --case-root /root/Workspace/x         # 只组装+运行
  $0 --steps 3,4 --amr                                # 组装+运行, 开启动态AMR
  $0 --steps 5 --case-root /root/Workspace/x          # 只做重建+位移后处理 (需先跑 3/4)
EOF
}

warn()  { echo "[FDM WARN] $*" >&2; }
error() { echo "[FDM ERR] $*" >&2; }

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
        --write-interval)  WRITE_INTERVAL="${2:?--write-interval 需要一个参数}"; shift 2 ;;
        --delta-t)          DELTA_T="${2:?--delta-t 需要一个参数}"; shift 2 ;;
        --video)           VIDEO=yes; shift ;;
        --no-video)        VIDEO=no; shift ;;
        --nprocs)          NPROCS="${2:?--nprocs 需要一个参数}"; shift 2 ;;
        --amr)             AMR=yes; shift ;;
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
    echo "[FDM-PIPE] 交互模式 (回车使用默认值)"
    read -p "算例根目录 --case-root (默认 $CASE_ROOT): " inp
    if [ -n "${inp:-}" ]; then CASE_ROOT="$inp"; fi

    read -p "网格密度 --res (默认 $RES): " inp
    if [ -n "${inp:-}" ]; then RES="$inp"; fi

    read -p "z 层数 --z-layers (默认 $Z_LAYERS): " inp
    if [ -n "${inp:-}" ]; then Z_LAYERS="$inp"; fi

    read -p "仿真时长 --end-time (默认 $END_TIME): " inp
    if [ -n "${inp:-}" ]; then END_TIME="$inp"; fi

    read -p "导出视频? [y/N]: " inp
    case "${inp:-}" in
        y|Y|yes|Yes|YES) VIDEO=yes ;;
        *)               VIDEO="" ;;   # 未定: 由 steps 是否含 6 决定
    esac

    read -p "并行进程数 --nprocs (默认 $NPROCS): " inp
    if [ -n "${inp:-}" ]; then NPROCS="$inp"; fi

    read -p "启用AMR (y/N, 默认 no): " inp
    case "${inp:-}" in
        y|Y|yes|Yes|YES) AMR=yes ;;
        *)               AMR=no ;;
    esac

    read -p "预写设置目录 --settings-dir (默认不使用): " inp
    if [ -n "${inp:-}" ]; then SETTINGS_DIR="$inp"; fi
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
#  步骤校验 (只接受 1-6 组合)
# --------------------------------------------------------------------------- #
validate_steps() {
    local s
    for s in $(echo "$STEPS" | tr ',' ' '); do
        case "$s" in
            1|2|3|4|5|6) ;;
            *) error "非法步骤 '$s' (只接受 1-6, 逗号分隔)"; return 1 ;;
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
        2)
            [ -f "$CASE_ROOT/mesh/fluid-openfoam/system/blockMeshDict" ] || \
                warn "步骤 2 需要网格产物 $CASE_ROOT/mesh/fluid-openfoam/system (先跑步骤 1)"
            ;;
        3)
            [ -d "$CASE_ROOT/mesh/fluid-openfoam/constant/polyMesh" ] || \
                warn "步骤 3 需要网格 $CASE_ROOT/mesh/fluid-openfoam/constant/polyMesh (先跑步骤 1)"
            [ -f "$CASE_ROOT/modal_ve/mode/FluidNodeCoor.csv" ] || \
                warn "步骤 3 需要模态文件 $CASE_ROOT/modal_ve/mode/FluidNodeCoor.csv (先跑步骤 2)"
            ;;
        4)
            [ -f "$CASE_ROOT/fdm_case/system/controlDict" ] || \
                warn "步骤 4 需要算例 $CASE_ROOT/fdm_case (先跑步骤 3)"
            ;;
        5)
            [ -d "$CASE_ROOT/fdm_case" ] || \
                warn "步骤 5 需要算例 $CASE_ROOT/fdm_case (先跑步骤 3/4)"
            ;;
        6)
            [ -d "$CASE_ROOT/fdm_case" ] || \
                warn "步骤 6 需要算例 $CASE_ROOT/fdm_case (先跑步骤 5 完成重建)"
            ;;
    esac
}

# --------------------------------------------------------------------------- #
#  捕获 + 过滤执行器 (步骤 1/2/3/5/6 及视频合成)
#  子脚本完整 stdout/stderr 落盘到 $1, 控制台只实时显示带 [FDM 前缀的主要信息行
#  (log/warn/error); 工具原始输出 (blockMesh/checkMesh/ccx/pvpython/foamToVTK/
#  ffmpeg 等无前缀行) 仅保留在日志, 不透传到控制台.
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
        2) cmd="python3 '$SCRIPT_DIR/02_run_modal_virtual_elastic.py' --fluid-dir '$CASE_ROOT/mesh/fluid-openfoam' --out-dir '$CASE_ROOT/modal_ve' --export-mode" ;;
        3) cmd="python3 '$SCRIPT_DIR/03_assemble_fdm_case.py' --mesh-dir '$CASE_ROOT/mesh/fluid-openfoam' --mode-dir '$CASE_ROOT/modal_ve/mode' --out-dir '$CASE_ROOT/fdm_case' --end-time $END_TIME --write-interval $WRITE_INTERVAL --nprocs $NPROCS" ;;
        4) cmd="python3 '$SCRIPT_DIR/04_run_fdm.py' -f --case '$CASE_ROOT/fdm_case' --end-time $END_TIME --write-interval $WRITE_INTERVAL --nprocs $NPROCS" ;;
        5) cmd="python3 '$SCRIPT_DIR/05_postprocess.py' --case '$CASE_ROOT/fdm_case' --out-dir '$CASE_ROOT/postprocess'" ;;
        6)
            cmd="xvfb-run -a pvpython '$SCRIPT_DIR/06_export_fdm_video.py' --case '$CASE_ROOT/fdm_case' --out '$CASE_ROOT/frames' --mp4-path '$CASE_ROOT/fdm_result.mp4'"
            [ "$VIDEO" = "no" ] && cmd+=" --no-mp4"
            ;;
    esac

    # AMR 只影响 03 生成的 constant/dynamicMeshDict (运行期由 fastDynamicFvMesh 执行)
    if [ "$step" -eq 3 ] && [ "$AMR" = "yes" ]; then
        cmd+=" --amr"
    fi
    if [ "$step" -eq 3 ] && [ -n "$SETTINGS_DIR" ]; then
        if [ "$DRY_RUN" -eq 0 ] && [ ! -d "$SETTINGS_DIR" ]; then
            error "--settings-dir 不存在或不是目录: $SETTINGS_DIR"
            return 1
        fi
        cmd+=" --settings-dir '$SETTINGS_DIR'"
    fi
    if [ "$step" -eq 3 ] && [ -n "$DELTA_T" ]; then
        cmd+=" --delta-t '$DELTA_T'"
    fi
    if [ "$step" -eq 4 ] && [ -n "$DELTA_T" ]; then
        cmd+=" --delta-t '$DELTA_T'"
    fi

    echo "[FDM-PIPE] 步骤 $step: $cmd"
    [ "$DRY_RUN" -eq 1 ] && return 0

    check_deps "$step"
    if [ "$step" -eq 4 ]; then
        # 步骤 4 (运行器): 自带 \r 进度条单行实时刷新, 捕获过滤会破坏刷新, 直接透传
        # (输出已较干净: [FDM 04] 前缀行 + 进度条 + solver Time=/ExecutionTime= 行)
        if ! eval "$cmd"; then
            error "步骤 $step 执行失败, 流水线中止"
            return 1
        fi
    elif [ "$step" -eq 6 ]; then
        # ParaView OpenFOAMReader 需要 case.foam 锚点文件; 03 不生成, 沿用仓库 touch 约定
        mkdir -p "$CASE_ROOT/fdm_case"
        [ -f "$CASE_ROOT/fdm_case/case.foam" ] || touch "$CASE_ROOT/fdm_case/case.foam"
        # 步骤 6: 非零退出只有在所需产物完整时才按 GLX teardown 警告处理。
        if ! run_script "$CASE_ROOT/run_step_6.log" "$cmd"; then
            frames_ok=0
            mp4_ok=0
            ls "$CASE_ROOT"/frames/frame_*.png >/dev/null 2>&1 && frames_ok=1
            [ -s "$CASE_ROOT/fdm_result.mp4" ] && mp4_ok=1
            if [ "$frames_ok" -eq 1 ] && { [ "$VIDEO" = "no" ] || [ "$mp4_ok" -eq 1 ]; }; then
                warn "步骤 6: 非零退出但 PNG/MP4 产物完整, 按无害 GLX 清理错误处理"
            else
                error "步骤 6 执行失败: PNG 或 MP4 产物不完整, 详见 run_step_6.log"
                return 1
            fi
        fi
    else
        # 步骤 1/2/3/5: 捕获过滤, 只透传带 [FDM 前缀的主要行, 工具原始输出仅落盘
        if ! run_script "$CASE_ROOT/run_step_$step.log" "$cmd"; then
            error "步骤 $step 执行失败, 流水线中止"
            return 1
        fi
    fi
}

# --------------------------------------------------------------------------- #
#  打印配置摘要
# --------------------------------------------------------------------------- #
echo "[FDM-PIPE] ===== Damfailure 六脚本流水线 (步骤 $STEPS) ====="
echo "[FDM-PIPE] 算例根目录 : $CASE_ROOT"
echo "[FDM-PIPE] 步骤       : $STEPS"
echo "[FDM-PIPE] 网格       : --res $RES --z-layers $Z_LAYERS"
echo "[FDM-PIPE] 仿真       : endTime=$END_TIME writeInterval=$WRITE_INTERVAL"
echo "[FDM-PIPE] 时间步长    : ${DELTA_T:-使用组装设置}"
echo "[FDM-PIPE] 并行       : $NPROCS 进程"
echo "[FDM-PIPE] AMR         : $AMR"
echo "[FDM-PIPE] 设置来源    : ${SETTINGS_DIR:-内嵌默认模板}"
if [ "$VIDEO" = "yes" ]; then
    VIDEO_STATE="步骤6默认合成 MP4 (--video)"
elif [ "$VIDEO" = "no" ]; then
    VIDEO_STATE="仅保留 PNG (--no-video)"
elif steps_contains 6; then
    VIDEO_STATE="步骤6默认合成 MP4 (25 fps)"
else
    VIDEO_STATE="不含步骤6"
fi
echo "[FDM-PIPE] 视频       : $VIDEO_STATE"
[ "$DRY_RUN" -eq 1 ] && echo "[FDM-PIPE] DRY-RUN    : 只打印命令, 不执行"
echo

# --------------------------------------------------------------------------- #
#  按步骤顺序执行 (去重 + 数值排序, 保证 1 < 2 < 3 < 4 < 5 < 6)
# --------------------------------------------------------------------------- #
for s in $(echo "$STEPS" | tr ',' '\n' | sort -u -n); do
    if [ "$s" -eq 4 ]; then
        set +e
        run_step "$s"
        STEP_04_RC=$?
        set -e
        if [ "$STEP_04_RC" -ne 0 ]; then
            STEP_04_FAILED=1
            error "步骤 4 执行失败 (exit $STEP_04_RC), 继续尝试步骤 5/6"
        fi
    else
        run_step "$s"
    fi
done

echo
echo "[FDM-PIPE] 流水线完成"
if [ "$STEP_04_FAILED" -eq 1 ]; then
    error "流水线完成但步骤 4 失败 (exit $STEP_04_RC)"
    exit 1
fi
