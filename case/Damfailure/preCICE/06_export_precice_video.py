#!/usr/bin/env python3
"""参考 FDM 06_export_fdm_video.py, 但读取 preCICE 组装的算例:
  --case <out-dir>  -> <out-dir>/fluid-openfoam/case.foam

按 alpha.water 云图着色 (白底 + 蓝色水相) + 白色底网格线框, 逐时间步渲染 PNG,
脚本默认以 30 fps、帧内压缩的 1080p H.264 在内部调用 ffmpeg 合成 MP4。

用法 (需在 xvfb 下运行 pvpython):
   xvfb-run -a pvpython 06_export_precice_video.py \
       --case /path/to/precice_case --out /path/to/frames [--res 1920x1080] \
       [--end-time 0.02] [--mp4-path /path/to/precice_result.mp4]

 关键经验 (ParaView 5.11.2, 与 FDM 05 相同):
  - OpenFOAMReader 无 GetNumberOfTimeSteps, 用 TimestepValues
  - 动画场景默认只有 10 帧, 用 GetTimeKeeper().Time 推进
   - 网格线用 ExtractSurface + "Surface With Edges", 不叠加独立 Wireframe rep。
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
        sys.stderr.write("[PRECICE ERR] 无效的 --res '%s', 应为 WxH (如 1920x1080)\n" % res)
        sys.exit(2)


def parse_color(s):
    """'0,0,0' 或 '255,0,0' 或 '0.5,0.5,0.5' -> [r,g,b] (0-1 浮点)"""
    try:
        parts = s.split(",")
        vals = [float(p) for p in parts]
        assert len(vals) == 3
    except Exception:
        sys.stderr.write("[PRECICE ERR] 无效的 --wire-color '%s', 应为 R,G,B (如 0,0,0 或 255,0,0)\n" % s)
        sys.exit(2)
    if any(v > 1.0 for v in vals):
        vals = [v / 255.0 for v in vals]
    return vals


def default_mp4_path(frames_dir):
    return os.path.join(os.path.dirname(os.path.normpath(frames_dir)),
                        "precice_result.mp4")


def encode_mp4(frames_dir, framerate, mp4_path):
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        sys.stderr.write("[PRECICE ERR] 找不到 ffmpeg; PNG 帧已保留在 %s\n" % frames_dir)
        return 3
    os.makedirs(os.path.dirname(os.path.abspath(mp4_path)), exist_ok=True)
    command = [
        ffmpeg, "-y", "-nostdin", "-hide_banner", "-loglevel", "error",
        "-framerate", str(framerate),
        "-i", os.path.join(frames_dir, "frame_%04d.png"),
        "-c:v", "libx264", "-preset", "slow", "-crf", "16",
        "-x264-params", "keyint=1:min-keyint=1:scenecut=0",
        "-pix_fmt", "yuv420p", "-movflags", "+faststart",
        mp4_path,
    ]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                            text=True, check=False)
    if result.returncode:
        sys.stderr.write("[PRECICE ERR] ffmpeg 合成失败 (exit %d): %s\n"
                         % (result.returncode, result.stderr.strip()))
        if os.path.isfile(mp4_path):
            os.remove(mp4_path)
        return 3
    print("[PRECICE 06] MP4 输出: %s (%d fps)" % (mp4_path, framerate))
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="preCICE fluid-openfoam 结果 -> 逐帧 PNG (alpha.water 云图)"
    )
    ap.add_argument("--case", required=True,
                    help="preCICE case 根目录 (03 输出, 读取其下 fluid-openfoam/case.foam)")
    ap.add_argument("--out", default="/root/Workspace/precice_frames",
                    help="输出 PNG 目录 (默认 /root/Workspace/precice_frames)")
    ap.add_argument("--res", default="1920x1080", help="分辨率 WxH (默认 1920x1080)")
    ap.add_argument("--field", default="alpha.water", help="着色场 (默认 alpha.water)")
    ap.add_argument("--wireframe", dest="wireframe", action="store_true", default=True,
                    help="叠加网格线框 (默认开启)")
    ap.add_argument("--no-wireframe", dest="wireframe", action="store_false",
                    help="关闭网格线框")
    ap.add_argument("--wire-color", default="0,0,0",
                    help="线框颜色 R,G,B, 0-255 或 0-1 (默认 0,0,0 黑色)")
    ap.add_argument("--end-time", type=float, default=None,
                    help="渲染到该时间截止 (默认 None = 渲染全部已有时步)")
    ap.add_argument("--include-t0", dest="include_t0", action="store_true", default=True,
                    help="先渲染 t=0 初始帧 (默认开启)")
    ap.add_argument("--no-t0", dest="include_t0", action="store_false",
                    help="跳过 t=0 初始帧")
    ap.add_argument("--mp4", dest="mp4", action="store_true", default=True,
                    help="合成 MP4 (默认开启)")
    ap.add_argument("--no-mp4", dest="mp4", action="store_false",
                    help="只输出 PNG 帧, 不合成 MP4")
    ap.add_argument("--framerate", "--fps", type=int, default=30,
                    help="MP4 帧率 (默认 30)")
    ap.add_argument("--mp4-path", default=None,
                    help="MP4 输出路径 (默认 <out 的父目录>/precice_result.mp4)")
    args = ap.parse_args()

    # pvpython 会用 vtkPythonStdStreamCaptureHelper 替换 sys.stdout, 非 TTY 下其
    # 内部缓冲只在正常 teardown 时冲刷; os._exit(0) 会把它整个丢弃 (flush 无效)。
    # 换回真实 stdout + 行缓冲, 保证成功路径上的 print 日志在 os._exit(0) 前落盘。
    sys.stdout = sys.__stdout__
    sys.stdout.reconfigure(line_buffering=True)

    print("[PRECICE-PIPE] ===== 步骤 6: 视频导出 =====")

    W, H = parse_res(args.res)
    if args.framerate <= 0:
        sys.stderr.write("[PRECICE ERR] --framerate 必须为正整数\n")
        return 2
    wire_color = parse_color(args.wire_color)

    case_foam = os.path.join(args.case, "fluid-openfoam", "case.foam")
    if not os.path.exists(case_foam):
        sys.stderr.write("[PRECICE ERR] 找不到 %s\n" % case_foam)
        return 2
    os.makedirs(args.out, exist_ok=True)
    for stale_frame in glob.glob(os.path.join(args.out, "frame_*.png")):
        os.remove(stale_frame)

    # ---- 1. 加载 OpenFOAM case ----
    reader = OpenFOAMReader(FileName=case_foam)
    reader.Decomposepolyhedra = 0
    reader.UpdatePipeline()

    timesteps = list(reader.TimestepValues) if reader.TimestepValues else []

    # OpenFOAMReader 默认 SkipZeroTime=1, 跳过 time=0。关闭后 reader 原生加载 0/
    # 初始场 (默认开启, --no-t0 关闭)。
    t0_dir = os.path.join(args.case, "fluid-openfoam", "0")
    if (args.include_t0 and os.path.isdir(t0_dir)
            and timesteps and min(timesteps) > 0.0):
        reader.SkipZeroTime = 0
        reader.UpdatePipelineInformation()
        timesteps = list(reader.TimestepValues) if reader.TimestepValues else []
        if 0.0 in timesteps:
            print("[PRECICE 05] include t=0 initial frame (case has 0/ directory)")
        else:
            sys.stderr.write("[PRECICE WARN] SkipZeroTime=0 后仍未读到 t=0 时步, 从首个写盘时步开始\n")

    # ---- 2. 按 --end-time 截断 (默认全部) ----
    if args.end_time is not None:
        keep = [t for t in timesteps if t <= args.end_time + 1e-9]
        if not keep:
            sys.stderr.write("[PRECICE ERR] --end-time %s 小于最小时步 %s, 无帧可渲染\n"
                             % (args.end_time, timesteps[0] if timesteps else "?"))
            return 2
        if len(keep) < len(timesteps):
            print("[PRECICE 05] --end-time %s: 截断 %d 个时步 -> %d 个"
                  % (args.end_time, len(timesteps), len(keep)))
        timesteps = keep

    # ---- 2b. 过滤 preCICE 空时间目录 (只有 uniform/, 无着色场文件) ----
    # preCICE adapter 在耦合过程中会残留空的时间目录 (如 0.0095/, 只有 uniform/,
    # 没有 alpha.water 场文件), OpenFOAMReader 的 TimestepValues 会把它们当时步,
    # 导致渲染出全白帧 (蓝像素=0 的假象)。用实际目录名匹配数值, 避免浮点
    # 格式化 (0.01 vs 0.010) 的坑。
    fluid_root = os.path.join(args.case, "fluid-openfoam")
    time_dirs = {}
    for d in glob.glob(os.path.join(fluid_root, "[0-9]*")):
        try:
            time_dirs[float(os.path.basename(d))] = d
        except ValueError:
            pass
    good = []
    skipped = []
    for t in timesteps:
        d = time_dirs.get(t)
        if d is not None and os.path.isfile(os.path.join(d, args.field)):
            good.append(t)
        else:
            skipped.append(t)
    if skipped:
        print("[PRECICE 05] skipped %d empty time-step(s) without '%s': %s"
              % (len(skipped), args.field, ", ".join("%g" % t for t in skipped)))
    timesteps = good
    if not timesteps:
        sys.stderr.write("[PRECICE ERR] 所有时间步都缺少着色场 '%s', 无帧可渲染\n" % args.field)
        return 2

    n_steps = len(timesteps)
    print("[PRECICE 05] Number of time steps: %d" % n_steps)
    if n_steps:
        print("[PRECICE 05] Time range: %s -> %s" % (timesteps[0], timesteps[-1]))
    if n_steps == 0:
        sys.stderr.write("[PRECICE ERR] 未读到任何时间步 (先运行 04_run_precice.py)\n")
        return 2

    reader.MeshRegions = ["internalMesh"]
    reader.UpdatePipeline()

    cell_fields = list(reader.CellArrays)
    point_fields = list(reader.PointArrays)
    print("[PRECICE 05] CellArrays: %s" % ", ".join(cell_fields))
    print("[PRECICE 05] PointArrays: %s" % ", ".join(point_fields))
    if args.field in point_fields:
        assoc = "POINTS"
    elif args.field in cell_fields:
        assoc = "CELLS"
    else:
        sys.stderr.write("[PRECICE ERR] 字段 '%s' 不存在, 可用 CellArrays=%s PointArrays=%s\n"
                         % (args.field, cell_fields, point_fields))
        return 2

    # ---- 3. 视图设置 (白底) ----
    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    # pvpython 批处理下直接设 Background 不生效 (实测仍是默认 82,87,110 深蓝灰),
    # 必须经调色板加载白底; 该颜色同时是网格外边框和挡板空洞的透出色。
    LoadPalette("WhiteBackground")
    view.Background = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 0

    # ---- 4. 显示设置: ExtractSurface + Surface With Edges ----
    es = ExtractSurface(Input=reader)
    es.UpdatePipeline()
    disp = GetDisplayProperties(es, view)
    disp.Representation = "Surface With Edges" if args.wireframe else "Surface"
    disp.Interpolation = "Flat"                  # 单元色不平滑, VOF 界面更锐利
    disp.Ambient = 1.0
    disp.Diffuse = 0.0
    disp.EdgeColor = wire_color                  # 仅 Surface With Edges 时可见
    disp.LineWidth = 0.5
    disp.ColorArrayName = [assoc, args.field]
    lut = GetColorTransferFunction(args.field)
    lut.RGBPoints = [0.0, 1.0, 1.0, 1.0,   # alpha=0 -> 白
                     1.0, 0.0, 0.2, 1.0]   # alpha=1 -> 蓝
    lut.ColorSpace = "RGB"
    disp.LookupTable = lut
    disp.RescaleTransferFunctionToDataRange(False)
    lut.RescaleTransferFunction(0.0, 1.0)


    # ---- 5. 相机: 正对 z 平面 ----
    cx, cy, cz = 0.292, 0.1825, 0.006
    view.CameraParallelProjection = 1
    view.CameraPosition = [cx, cy, 1.0]
    view.CameraFocalPoint = [cx, cy, cz]
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = 0.21

    # ---- 6. 逐时间步推进 ----
    tk = GetTimeKeeper()
    for i, t in enumerate(timesteps):
        tk.Time = t
        reader.UpdatePipeline()
        es.UpdatePipeline()
        view.Update()
        fname = os.path.join(args.out, "frame_%04d.png" % i)
        SaveScreenshot(fname, view, ImageResolution=[W, H])
        print("[PRECICE 05] --- saved %s (t=%.6g)" % (fname, t))

    print("[PRECICE 06] 帧数: %d  输出: %s" % (n_steps, args.out))

    mp4_status = 0
    if args.mp4:
        mp4_status = encode_mp4(args.out, args.framerate,
                                args.mp4_path or default_mp4_path(args.out))
    # pvpython 在 xvfb 下渲染完成后, Python 正常退出时 ParaView/X teardown 会抛
    # 无害 GLXBadContext 导致退出码非 0 (帧已全部生成)。用 os._exit(0) 绕过该
    # X 清理阶段, 干净退出。
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(mp4_status)


if __name__ == "__main__":
    sys.exit(main())
