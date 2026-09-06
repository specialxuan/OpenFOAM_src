#!/usr/bin/env python3
"""06_export_fdm_video.py - 导出 FDM PNG 帧并合成 MP4

FDM = fastDynamicFvMesh + myInterFoam。读取 case.foam，按 alpha.water 云图着色
（白底 + 蓝色水相），逐时间步保存 PNG，并默认用 ffmpeg 合成 30 fps 的高质量 H.264
帧内编码 MP4。

用法:
  xvfb-run -a pvpython 06_export_fdm_video.py \
       --case <case_dir> --out <out_dir> [--res 1920x1080] [--field alpha.water] \
       [--no-wireframe] [--wire-color 0,0,0] [--framerate 30]
      [--mp4-path <path>] [--no-mp4]

 关键经验 (ParaView 5.11.2):
  - Ubuntu paraview 是 X11 编译版, 必须 xvfb-run 提供 DISPLAY
  - OpenFOAMReader 无 VolumeFields / GetNumberOfTimeSteps, 用 TimestepValues
  - GetDisplayProperties 需要 view 参数
  - 动画场景默认只有 10 帧, 不要用 anim.GoToNext(), 用 GetTimeKeeper().Time 推进
  - 网格线用 ExtractSurface(外表面) + "Surface With Edges", 不叠加独立 Wireframe rep,
    避免 AMR 内部 hex 体对角线在准 2D 投影下形成斜三角线
  - 若 xvfb/OSMesa 下仍出现边线伪影, 用 --no-wireframe 导出纯 Surface 画面
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

from paraview.simple import *


def parse_res(res):
    """'1920x1080' -> (1920, 1080)"""
    try:
        w, h = res.lower().split("x")
        w, h = int(w), int(h)
        assert w > 0 and h > 0
        return w, h
    except Exception:
        sys.stderr.write("[FDM ERR] 无效的 --res '%s', 应为 WxH (如 1920x1080)\n" % res)
        sys.exit(2)


def parse_color(s):
    """'0,0,0' 或 '255,0,0' 或 '0.5,0.5,0.5' -> [r,g,b] (0-1 浮点)"""
    try:
        parts = s.split(",")
        vals = [float(p) for p in parts]
        assert len(vals) == 3
    except Exception:
        sys.stderr.write("[FDM ERR] 无效的 --wire-color '%s', 应为 R,G,B (如 0,0,0 或 255,0,0)\n" % s)
        sys.exit(2)
    # 任一分量 > 1 视为 0-255 范围, 归一化到 0-1
    if any(v > 1.0 for v in vals):
        vals = [v / 255.0 for v in vals]
    return vals


def default_mp4_path(frames_dir):
    return os.path.join(os.path.dirname(os.path.normpath(frames_dir)),
                        "fdm_result.mp4")


def encode_mp4(frames_dir, framerate, mp4_path):
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        sys.stderr.write("[FDM ERR] 找不到 ffmpeg; PNG 帧已保留在 %s\n" % frames_dir)
        return 3
    parent = os.path.dirname(os.path.abspath(mp4_path))
    os.makedirs(parent, exist_ok=True)
    command = [
        ffmpeg, "-y", "-nostdin", "-hide_banner", "-loglevel", "error",
        "-framerate", str(framerate),
        "-i", os.path.join(frames_dir, "frame_%04d.png"),
        "-c:v", "libx264", "-preset", "slow", "-crf", "16",
        "-x264-params", "keyint=1:min-keyint=1:scenecut=0",
        "-pix_fmt", "yuv420p", "-movflags", "+faststart",
        mp4_path,
    ]
    result = subprocess.run(
        command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        text=True, check=False)
    if result.returncode != 0:
        sys.stderr.write("[FDM ERR] ffmpeg 合成失败 (exit %d): %s\n"
                         % (result.returncode, result.stderr.strip()))
        if os.path.isfile(mp4_path):
            os.remove(mp4_path)
        return 3
    print("[FDM 06] MP4 输出: %s (%d fps)" % (mp4_path, framerate))
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="FDM OpenFOAM 结果 -> 逐帧 PNG (alpha.water 云图)"
    )
    ap.add_argument("--case", required=True, help="OpenFOAM case 目录 (含 case.foam)")
    ap.add_argument("--out", default="/root/Workspace/video_frames",
                    help="输出 PNG 目录 (默认 /root/Workspace/video_frames)")
    ap.add_argument("--res", default="1920x1080", help="分辨率 WxH (默认 1920x1080)")
    ap.add_argument("--field", default="alpha.water", help="着色场 (默认 alpha.water)")
    ap.add_argument("--end-time", type=float, default=None,
                    help="渲染到该时间截止 (默认 None = 渲染全部已有时步)")
    ap.add_argument("--wireframe", dest="wireframe", action="store_true", default=True,
                    help="叠加网格线框 (默认开启)")
    ap.add_argument("--no-wireframe", dest="wireframe", action="store_false",
                    help="关闭网格线框")
    ap.add_argument("--wire-color", default="0,0,0",
                    help="线框颜色 R,G,B, 0-255 或 0-1 (默认 0,0,0 黑色)")
    ap.add_argument("--include-t0", dest="include_t0", action="store_true", default=True,
                    help="render the t=0 initial-state frame first (default: on)")
    ap.add_argument("--no-t0", dest="include_t0", action="store_false",
                    help="skip the t=0 initial-state frame (start at first written time step)")
    ap.add_argument("--mp4", dest="mp4", action="store_true", default=True,
                    help="合成 MP4 (默认开启)")
    ap.add_argument("--no-mp4", dest="mp4", action="store_false",
                    help="只输出 PNG 帧, 不合成 MP4")
    ap.add_argument("--framerate", "--fps", type=int, default=30,
                    help="MP4 帧率 (默认 30)")
    ap.add_argument("--mp4-path", default=None,
                    help="MP4 输出路径 (默认 <out 的父目录>/fdm_result.mp4)")
    args = ap.parse_args()

    # pvpython 会用 vtkPythonStdStreamCaptureHelper 替换 sys.stdout, 非 TTY 下其
    # 内部缓冲只在正常 teardown 时冲刷; os._exit(0) 会把它整个丢弃 (flush 无效)。
    # 换回真实 stdout + 行缓冲, 保证成功路径上的 print 日志在 os._exit(0) 前落盘。
    sys.stdout = sys.__stdout__
    sys.stdout.reconfigure(line_buffering=True)

    print("[FDM-PIPE] ===== 步骤 6/6: 视频导出 =====")

    W, H = parse_res(args.res)
    if args.framerate <= 0:
        sys.stderr.write("[FDM ERR] --framerate 必须为正整数\n")
        return 2
    wire_color = parse_color(args.wire_color)
    os.makedirs(args.out, exist_ok=True)
    for stale_frame in glob.glob(os.path.join(args.out, "frame_*.png")):
        os.remove(stale_frame)

    case_foam = os.path.join(args.case, "case.foam")
    if not os.path.exists(case_foam):
        sys.stderr.write("[FDM ERR] 找不到 %s\n" % case_foam)
        return 2

    # ---- 1. 加载 OpenFOAM case ----
    reader = OpenFOAMReader(FileName=case_foam)
    # AMR 粗细过渡会产生五/六边形面；默认分解会生成 pyramid/tetra，
    # Surface With Edges 随后把三角面内部边显示为对角线。
    reader.Decomposepolyhedra = 0
    reader.UpdatePipeline()

    timesteps = list(reader.TimestepValues) if reader.TimestepValues else []

    # OpenFOAMReader 默认 SkipZeroTime=1, 跳过 time=0: TimestepValues 不含 0,
    # 此时 tk.Time=0 会被执行器钳制到首个写盘时步, 帧实为 t=0.01。关闭
    # SkipZeroTime 让 reader 原生加载 0/ 初始场 (默认开启, --no-t0 关闭).
    t0_dir = os.path.join(args.case, "0")
    if (args.include_t0 and os.path.isdir(t0_dir)
            and timesteps and min(timesteps) > 0.0):
        reader.SkipZeroTime = 0
        reader.UpdatePipelineInformation()
        timesteps = list(reader.TimestepValues) if reader.TimestepValues else []
        if 0.0 in timesteps:
            print("[FDM 06] include t=0 initial frame (case has 0/ directory)")
        else:
            sys.stderr.write("[FDM WARN] SkipZeroTime=0 后仍未读到 t=0 时步, 从首个写盘时步开始\n")

    # ---- 2. 按 --end-time 截断 (默认全部) ----
    if args.end_time is not None:
        keep = [t for t in timesteps if t <= args.end_time + 1e-9]
        if not keep:
            sys.stderr.write("[FDM ERR] --end-time %s 小于最小时步 %s, 无帧可渲染\n"
                             % (args.end_time, timesteps[0] if timesteps else "?"))
            return 2
        if len(keep) < len(timesteps):
            print("[FDM 06] --end-time %s: 截断 %d 个时步 -> %d 个"
                  % (args.end_time, len(timesteps), len(keep)))
        timesteps = keep

    n_steps = len(timesteps)
    print("[FDM 06] Number of time steps: %d" % n_steps)
    if n_steps:
        print("[FDM 06] Time range: %s -> %s" % (timesteps[0], timesteps[-1]))
    if n_steps == 0:
        sys.stderr.write("[FDM ERR] 未读到任何时间步\n")
        return 2

    reader.MeshRegions = ["internalMesh"]
    reader.UpdatePipeline()

    # ---- 调试: 确认字段存在, 并判断是 cell 还是 point 数据 ----
    # OpenFOAMReader 5.11 把 volScalarField 作为 CELL 数据加载
    cell_fields = list(reader.CellArrays)
    point_fields = list(reader.PointArrays)
    print("[FDM 06] CellArrays: %s" % ", ".join(cell_fields))
    print("[FDM 06] PointArrays: %s" % ", ".join(point_fields))
    if args.field in point_fields:
        assoc = "POINTS"
    elif args.field in cell_fields:
        assoc = "CELLS"
    else:
        sys.stderr.write("[FDM ERR] 字段 '%s' 不存在, 可用 CellArrays=%s PointArrays=%s\n"
                         % (args.field, cell_fields, point_fields))
        return 2

    # ---- 2. 视图设置 ----
    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    # pvpython 批处理下直接设 Background 不生效 (实测仍是默认 82,87,110 深蓝灰),
    # 必须经调色板加载白底; 该颜色同时是网格外边框和挡板空洞的透出色。
    LoadPalette("WhiteBackground")
    view.Background = [1.0, 1.0, 1.0]          # 白底
    view.OrientationAxesVisibility = 0

    # ---- 3. 显示设置: ExtractSurface(外表面) + Surface With Edges ----
    # 关键: ExtractSurface 只提取外表面边界 patch 的面, 这些面在 x-y/x-z 平面,
    # 边规则正交, 无 AMR 内部 hex 的体对角线斜边 -> 消除视频斜三角线.
    es = ExtractSurface(Input=reader)
    es.UpdatePipeline()
    disp = GetDisplayProperties(es, view)
    disp.Representation = "Surface With Edges" if args.wireframe else "Surface"
    disp.Interpolation = "Flat"                  # 单元色不平滑, VOF 界面更锐利
    disp.Ambient = 1.0                           # 无光照 -> 纯白底 + 纯蓝水相
    disp.Diffuse = 0.0
    disp.EdgeColor = wire_color                  # 仅 Surface With Edges 时可见
    disp.LineWidth = 0.5
    disp.ColorArrayName = [assoc, args.field]
    lut = GetColorTransferFunction(args.field)
    # 白底 + 蓝色水相: alpha=0 -> 白色, alpha=1 -> 蓝色
    lut.RGBPoints = [0.0, 1.0, 1.0, 1.0,   # alpha=0 -> 白
                     1.0, 0.0, 0.2, 1.0]   # alpha=1 -> 蓝
    lut.ColorSpace = "RGB"
    disp.LookupTable = lut
    disp.RescaleTransferFunctionToDataRange(False)
    lut.RescaleTransferFunction(0.0, 1.0)

    # 注: 不再叠加独立 Wireframe representation (旧 3b 段已删除)。独立 Wireframe 会画
    # 所有 z 层 hex 的三角剖分/体斜边, 在准 2D 多 z 层投影下重叠成密集"斜三角线"。
    # "Surface With Edges" 只画 ExtractSurface 外表面边；若软件渲染仍暴露 AMR
    # 表面三角化/悬挂节点边, --no-wireframe 切到纯 Surface 作为严格无斜线路径。

    # ---- 4. 相机: 正对 z 平面 (看到 x-y 平面的溃坝过程) ----
    # 几何: x [0, 0.584], y [0, 0.365], z [0, 0.012] (准二维, z 为薄方向)
    cx, cy, cz = 0.292, 0.1825, 0.006
    view.CameraParallelProjection = 1
    view.CameraPosition = [cx, cy, 1.0]       # 相机在 +z 高位置
    view.CameraFocalPoint = [cx, cy, cz]      # 看向 z=0.006 (中间高度, 显示完整 x-y 水槽)
    view.CameraViewUp = [0, 1, 0]             # y 向上
    # 平行投影缩放: 覆盖整个 x-y 域 (0.584 x 0.365), 留少量边距
    view.CameraParallelScale = 0.21

    # ---- 5. 逐时间步推进 (用 TimeKeeper, 不用动画场景) ----
    tk = GetTimeKeeper()
    for i, t in enumerate(timesteps):
        tk.Time = t
        reader.UpdatePipeline()
        es.UpdatePipeline()
        view.Update()
        fname = os.path.join(args.out, "frame_%04d.png" % i)
        SaveScreenshot(fname, view, ImageResolution=[W, H])
        print("[FDM 06] --- saved %s (t=%.6g)" % (fname, t))

    print("[FDM 06] 帧数: %d  输出: %s" % (n_steps, args.out))
    mp4_status = 0
    if args.mp4:
        mp4_path = args.mp4_path or default_mp4_path(args.out)
        mp4_status = encode_mp4(args.out, args.framerate, mp4_path)
    # pvpython 在 xvfb 下渲染完成后, Python 正常退出时 ParaView/X teardown 会抛
    # 无害 GLXBadContext 导致退出码非 0 (帧已全部生成)。用 os._exit(0) 绕过该
    # X 清理阶段；MP4 必须在此之前完成。
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(mp4_status)


if __name__ == "__main__":
    sys.exit(main())
