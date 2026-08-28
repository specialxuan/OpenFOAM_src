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
VIDEO=""          # "" 未定(交互询问) / "yes" / "no"
NPROCS=8
AMR="no"          # "yes" / "no" (默认关闭)
DRY_RUN=0

usage() {
    cat <<EOF
用法: $0 [选项]

驱动 Damfailure 六脚本流水线:
  01 网格生成 -> 02 模态分析 -> 03 组装算例 -> 04 FDM 运行
  -> 05 并行重建/后处理 -> 06 视频导出

选项 (全部可选; 不带任何参数时交互询问):
  --case-root <dir>      算例根目录 (默认 $CASE_ROOT)
  --steps <list>         运行步骤, 逗号分隔 1-6 (默认 $STEPS), 如 "3,4"; 含 6 则自动导出视频 (帧+MP4)
  --res <N>              网格 x 向密度 (传给 01, 默认 $RES)
  --z-layers <N>         z 向层数 (传给 01, 默认 $Z_LAYERS)
  --end-time <T>         仿真时长 (传给 03/04, 默认 $END_TIME)
  --write-interval <N>   输出间隔步数 (默认 $WRITE_INTERVAL)
  --video / --no-video   覆盖 steps 的视频行为 (默认随 steps 含 6 自动导出; --no-video 强制关闭)
  --nprocs <N>           并行进程数 (传给 03/04, 默认 $NPROCS)
  --amr                  启用动态自适应网格 (dynamic AMR, 传给 03, 默认关闭)
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
        --video)           VIDEO=yes; shift ;;
        --no-video)        VIDEO=no; shift ;;
        --nprocs)          NPROCS="${2:?--nprocs 需要一个参数}"; shift 2 ;;
        --amr)             AMR=yes; shift ;;
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
        *)               VIDEO="" ;;   # 未定: 由 steps 是否含 5 决定
    esac

    read -p "并行进程数 --nprocs (默认 $NPROCS): " inp
    if [ -n "${inp:-}" ]; then NPROCS="$inp"; fi

    read -p "启用AMR (y/N, 默认 no): " inp
    case "${inp:-}" in
        y|Y|yes|Yes|YES) AMR=yes ;;
        *)               AMR=no ;;
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
#  视频合成辅助
# --------------------------------------------------------------------------- #

# 判断步骤列表是否包含某步骤 (如 steps_contains 6)
steps_contains() {
    echo ",$STEPS," | grep -q ",$1,"
}

# ffmpeg 把步骤 6 生成的 PNG 帧合成 MP4 (含帧目录检查 / DRY_RUN 处理)
synthesize_video() {
    ffmpeg_cmd="ffmpeg -y -framerate 10 -i '$CASE_ROOT/frames/frame_%04d.png' -c:v libx264 -pix_fmt yuv420p '$CASE_ROOT/fdm_result.mp4'"
    echo "[FDM-PIPE] 视频合成: $ffmpeg_cmd"
    if [ "$DRY_RUN" -eq 0 ]; then
        if [ ! -d "$CASE_ROOT/frames" ]; then
            error "未找到帧目录 $CASE_ROOT/frames (需要先运行步骤 6)"
            exit 1
        fi
        if ! run_script "$CASE_ROOT/run_video.log" "$ffmpeg_cmd"; then
            error "视频合成失败"
            exit 1
        fi
        echo "[FDM-PIPE] 视频输出: $CASE_ROOT/fdm_result.mp4"
    fi
}

# 合成触发决策 (优先级: --video 强制 > --no-video 强制关闭 > steps 含 6 自动)
maybe_synthesize_video() {
    if [ "$VIDEO" = "yes" ]; then
        synthesize_video
    elif [ "$VIDEO" = "no" ]; then
        :
    elif steps_contains 6; then
        synthesize_video
    fi
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
        6) cmd="xvfb-run -a pvpython '$SCRIPT_DIR/06_export_fdm_video.py' --case '$CASE_ROOT/fdm_case' --out '$CASE_ROOT/frames'" ;;
    esac

    # AMR 只影响 03 生成的 constant/dynamicMeshDict (运行期由 fastDynamicFvMesh 执行)
    if [ "$step" -eq 3 ] && [ "$AMR" = "yes" ]; then
        cmd+=" --amr"
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
        # 步骤 6: pvpython 在 xvfb 下渲染完成后可能抛无害 GLXBadContext 导致非零退出码,
        # 但只要帧文件已生成即视为成功。捕获过滤: 只透传 [FDM 前缀行, 完整输出落盘。
        if ! run_script "$CASE_ROOT/run_step_6.log" "$cmd"; then
            if ls "$CASE_ROOT"/frames/frame_*.png >/dev/null 2>&1; then
                warn "步骤 6: pvpython 退出码非 0 (xvfb GLX 清理错误) 但帧已生成, 视为成功"
            else
                error "步骤 6 执行失败: 未生成任何帧 (pvpython 非零退出)"
                return 1
            fi
        fi
        # 出帧成功后合成视频: steps 含 6 自动; --no-video 强制关闭
        maybe_synthesize_video
    else
        # 步骤 1/2/3/6: 捕获过滤, 只透传带 [FDM 前缀的主要行, 工具原始输出仅落盘
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
echo "[FDM-PIPE] 并行       : $NPROCS 进程"
echo "[FDM-PIPE] AMR         : $AMR"
if [ "$VIDEO" = "yes" ]; then
    VIDEO_STATE="强制是 (--video)"
elif [ "$VIDEO" = "no" ]; then
    VIDEO_STATE="强制否 (--no-video)"
elif steps_contains 6; then
    VIDEO_STATE="随步骤6 (含 6, 自动合成)"
else
    VIDEO_STATE="随步骤6 (不含 6, 不合成)"
fi
echo "[FDM-PIPE] 视频       : $VIDEO_STATE"
[ "$DRY_RUN" -eq 1 ] && echo "[FDM-PIPE] DRY-RUN    : 只打印命令, 不执行"
echo

# --------------------------------------------------------------------------- #
#  按步骤顺序执行 (去重 + 数值排序, 保证 1 < 2 < 3 < 4 < 5 < 6)
# --------------------------------------------------------------------------- #
for s in $(echo "$STEPS" | tr ',' '\n' | sort -u -n); do
    run_step "$s"
done

# --------------------------------------------------------------------------- #
#  视频合成 (步骤 6 生成 PNG 帧后, 用 ffmpeg 合成 MP4)
#  - dry-run 预览: 步骤 6 分支提前返回未触发合成, 这里按决策打印合成命令
#  - 实跑: 步骤 6 分支已合成 (steps 含 6 自动 / --video 强制 / --no-video 关闭);
#    这里仅兜底 --video 且 steps 不含 6 (步骤 6 未运行 -> 无帧, 报错提示)
# --------------------------------------------------------------------------- #
if [ "$DRY_RUN" -eq 1 ]; then
    maybe_synthesize_video
elif [ "$VIDEO" = "yes" ] && ! steps_contains 6; then
    synthesize_video
fi

echo
echo "[FDM-PIPE] 流水线完成"
