#!/usr/bin/env python3
"""06_export_fdm_video.py - 导出 FDM PNG 帧并合成 MP4

FDM = fastDynamicFvMesh + myInterFoam。读取 case.foam，按指定场云图着色
（默认仍为 alpha.water：白底 + 蓝色水相），逐时间步保存 PNG，并默认用 ffmpeg
合成 30 fps 的高质量 H.264 帧内编码 MP4。

--field magU 时自动构建 Calculator 标量 mag(U)（U 是矢量场，不能直接着色），
用于本单相满水算例的速度云图；此时标量条（color legend）默认开启。

用法:
  xvfb-run -a pvpython 06_export_fdm_video.py \
       --case <case_dir> --out <out_dir> [--res 1920x1080] [--field alpha.water] \
       [--field magU --vector-field U] [--lut-preset Rainbow_Desaturated] \
       [--color-range 0,0.05] [--color-range-mode auto|current] \
       [--scalar-bar|--no-scalar-bar] [--scalar-bar-title <t>] \
       [--scalar-bar-location 'Lower Right Corner'] [--scalar-bar-format %.4g] \
       [--interpolation auto|flat|gouraud] \
       [--camera-focus 0.298,0.06,0.006] [--camera-scale 0.08] [--zoom-beam] \
       [--no-wireframe] [--wire-color 0,0,0] [--framerate 30] \
       [--mp4-path <path>] [--no-mp4]

 关键经验 (ParaView 5.11.2):
  - Ubuntu paraview 是 X11 编译版, 必须 xvfb-run 提供 DISPLAY
  - OpenFOAMReader 无 VolumeFields / GetNumberOfTimeSteps, 用 TimestepValues
  - GetDisplayProperties 需要 view 参数
  - 动画场景默认只有 10 帧, 不要用 anim.GoToNext(), 用 GetTimeKeeper().Time 推进
  - 时间步进必须显式给 UpdatePipeline(t) 传时间: 客户端读取数组范围
    (PointData[]/GetDataInformation) 在仅设 tk.Time 后会返回陈旧数据, 颜色范围
    预扫描会得到全部相同的范围; reader/filter 链显式 UpdatePipeline(t) 后一致
  - mag(U): ExtractSurface 上 U 同时有 CELL/POINT 副本; Calculator 默认
    AttributeType='Point Data' -> 平滑 Gouraud 速度云图; 关联码 0=POINTS, 1=CELLS
  - 标量条: GetScalarBar(lut, view), WindowLocation 合法值含 'Lower Right Corner'
    (无 'Lower Right'); Position2 在 5.11 已废弃, 不要用
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

MAG_FIELD = "magU"
ALPHA_FIELD = "alpha.water"
FULL_FOCUS = (0.292, 0.18, 0.006)
FULL_SCALE = 0.21
BEAM_FOCUS = (0.298, 0.06, 0.006)   # 悬臂梁 x=[0.292,0.304] y=[0,0.120] 的中心
BEAM_SCALE = 0.08                    # 视野高 ~0.16 m, 宽 ~0.284 m, 梁占满高度
SCALAR_BAR_LOCATIONS = [
    "Any Location", "Lower Left Corner", "Lower Right Corner", "Lower Center",
    "Upper Left Corner", "Upper Right Corner", "Upper Center",
]


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


def parse_vec3(s, opt):
    """'0.298,0.06,0.006' -> [0.298, 0.06, 0.006]"""
    try:
        vals = [float(p) for p in s.split(",")]
        assert len(vals) == 3
        return vals
    except Exception:
        sys.stderr.write("[FDM ERR] 无效的 %s '%s', 应为 X,Y,Z\n" % (opt, s))
        sys.exit(2)


def parse_range(s):
    """'0,0.05' -> (0.0, 0.05)"""
    try:
        parts = [float(p) for p in s.split(",")]
        assert len(parts) == 2 and parts[0] < parts[1]
        return parts[0], parts[1]
    except Exception:
        sys.stderr.write("[FDM ERR] 无效的 --color-range '%s', 应为 'min,max' 且 min<max\n" % s)
        sys.exit(2)


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


def read_array_range(src, name, assoc):
    """读取当前 pipeline 状态下 src 上 name 的分量 0 范围 (min, max)。

    必须在 UpdatePipeline(t) 之后调用; 客户端数组读取不会随 tk.Time 自动刷新。
    """
    arr = src.PointData[name] if assoc == "POINTS" else src.CellData[name]
    return arr.GetRange(0)


def main():
    ap = argparse.ArgumentParser(
        description="FDM OpenFOAM 结果 -> 逐帧 PNG (alpha.water 或 mag(U) 云图)"
    )
    ap.add_argument("--case", required=True, help="OpenFOAM case 目录 (含 case.foam)")
    ap.add_argument("--out", default="/root/Workspace/video_frames",
                    help="输出 PNG 目录 (默认 /root/Workspace/video_frames)")
    ap.add_argument("--res", default="1920x1080", help="分辨率 WxH (默认 1920x1080)")
    ap.add_argument("--field", default=ALPHA_FIELD,
                    help="着色场 (默认 alpha.water); 传 magU 时自动计算 |U| 云图")
    ap.add_argument("--vector-field", default="U",
                    help="--field magU 时求模的矢量场名 (默认 U)")
    ap.add_argument("--lut-preset", default=None,
                    help="ParaView 调色板预设名 (默认: alpha.water 用白->蓝, "
                         "其他场用 Rainbow Desaturated)")
    ap.add_argument("--color-range", default=None,
                    help="固定色标范围 'min,max' (默认自动)")
    ap.add_argument("--color-range-mode", choices=("auto", "current"), default="auto",
                    help="auto: 预扫描全时段固定范围; current: 逐帧自适应 (会闪烁)")
    ap.add_argument("--scalar-bar", dest="scalar_bar", action="store_true", default=None,
                    help="显示标量条/颜色图例 (默认: 非 alpha.water 场开启)")
    ap.add_argument("--no-scalar-bar", dest="scalar_bar", action="store_false",
                    help="关闭标量条")
    ap.add_argument("--scalar-bar-title", default=None,
                    help="标量条标题 (默认: magU -> 'Velocity Magnitude (m/s)')")
    ap.add_argument("--scalar-bar-location", default="Lower Right Corner",
                    choices=SCALAR_BAR_LOCATIONS, help="标量条位置 (默认 Lower Right Corner)")
    ap.add_argument("--scalar-bar-format", default="%.4g", help="标量条标签格式 (默认 %%.4g)")
    ap.add_argument("--interpolation", choices=("auto", "flat", "gouraud"), default="auto",
                    help="表面插值 (默认: alpha.water flat, 其他场点数据 gouraud)")
    ap.add_argument("--camera-focus", default=None,
                    help="相机焦点 X,Y,Z (默认全流域中心 0.292,0.18,0.006)")
    ap.add_argument("--camera-scale", type=float, default=None,
                    help="平行投影尺度 (默认 0.21 全流域; --zoom-beam 时 0.08)")
    ap.add_argument("--zoom-beam", action="store_true",
                    help="梁区放大相机: 焦点 0.298,0.06,0.006, 尺度 0.08")
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

    field_key = args.field.strip()
    mag_mode = field_key.lower() in (MAG_FIELD.lower(), "mag(u)", "|u|")
    color_name = MAG_FIELD if mag_mode else field_key
    if not mag_mode and field_key.upper() == args.vector_field.upper():
        sys.stderr.write("[FDM ERR] '%s' 是矢量场, 请用 --field magU 渲染速度幅值云图\n"
                         % field_key)
        return 2
    fixed_range = parse_range(args.color_range) if args.color_range else None

    # ---- 0. 相机解析 (先算好打印, 避免中途才发现参数错误) ----
    if args.zoom_beam:
        cam_focus = list(BEAM_FOCUS)
        cam_scale = BEAM_SCALE
    else:
        cam_focus = list(FULL_FOCUS)
        cam_scale = FULL_SCALE
    if args.camera_focus is not None:
        cam_focus = parse_vec3(args.camera_focus, "--camera-focus")
    if args.camera_scale is not None:
        cam_scale = args.camera_scale
    if cam_scale <= 0:
        sys.stderr.write("[FDM ERR] --camera-scale 必须为正数\n")
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

    # ---- 3. 显示源: ExtractSurface(外表面) + (可选) Calculator mag(U) ----
    # 关键: ExtractSurface 只提取外表面边界 patch 的面, 这些面在 x-y/x-z 平面,
    # 边规则正交, 无 AMR 内部 hex 的体对角线斜边 -> 消除视频斜三角线.
    es = ExtractSurface(Input=reader)
    es.UpdatePipeline()

    calc = None
    if mag_mode:
        vec = args.vector_field
        es_cells = list(es.CellData.keys())
        es_points = list(es.PointData.keys())
        if vec in es_points:
            vec_assoc = "POINTS"
        elif vec in es_cells:
            vec_assoc = "CELLS"
        else:
            sys.stderr.write("[FDM ERR] 矢量场 '%s' 不存在, ES CellArrays=%s PointArrays=%s\n"
                             % (vec, es_cells, es_points))
            return 2
        calc = Calculator(Input=es)
        # ExtractSurface 上 U 同时存在 CELL/POINT 副本; 默认用点数据获得平滑云图
        calc.AttributeType = "Point Data" if vec_assoc == "POINTS" else "Cell Data"
        calc.ResultArrayName = MAG_FIELD
        calc.Function = "mag(%s)" % vec
        calc.UpdatePipeline()
        assoc = "POINTS" if vec_assoc == "POINTS" else "CELLS"
        print("[FDM 06] magU = mag(%s) on %s data" % (vec, vec_assoc))
    else:
        if field_key in point_fields:
            assoc = "POINTS"
        elif field_key in cell_fields:
            assoc = "CELLS"
        else:
            sys.stderr.write("[FDM ERR] 字段 '%s' 不存在, 可用 CellArrays=%s PointArrays=%s\n"
                             % (field_key, cell_fields, point_fields))
            return 2
        assoc = "POINTS" if assoc == "POINTS" else "CELLS"

    src = calc if calc is not None else es

    # ---- 4. 视图设置 ----
    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    # pvpython 批处理下直接设 Background 不生效 (实测仍是默认 82,87,110 深蓝灰),
    # 必须经调色板加载白底; 该颜色同时是网格外边框和挡板空洞的透出色。
    LoadPalette("WhiteBackground")
    view.Background = [1.0, 1.0, 1.0]          # 白底
    view.OrientationAxesVisibility = 0

    disp = GetDisplayProperties(src, view)
    disp.Representation = "Surface With Edges" if args.wireframe else "Surface"
    # alpha.water 用 Flat (VOF 界面锐利, 保持旧行为); magU 点数据默认 Gouraud 平滑
    if args.interpolation == "auto":
        disp.Interpolation = "Gouraud" if assoc == "POINTS" and color_name != ALPHA_FIELD \
            else "Flat"
    else:
        disp.Interpolation = {"flat": "Flat", "gouraud": "Gouraud"}[args.interpolation]
    disp.Ambient = 1.0                           # 无光照 -> 纯白底 + 纯色场
    disp.Diffuse = 0.0
    disp.EdgeColor = wire_color                  # 仅 Surface With Edges 时可见
    disp.LineWidth = 0.5
    disp.ColorArrayName = [assoc, color_name]
    lut = GetColorTransferFunction(color_name)

    # 调色板: alpha.water 保持白底 + 蓝色水相 (旧行为); 其他场用可读彩虹
    preset = args.lut_preset
    if preset is None and color_name != ALPHA_FIELD:
        preset = "Rainbow Desaturated"
    if preset:
        try:
            lut.ApplyPreset(preset, True)
            print("[FDM 06] LUT preset: %s" % preset)
        except Exception as exc:
            sys.stderr.write("[FDM WARN] LUT 预设 '%s' 不可用 (%s), 使用默认调色板\n"
                             % (preset, exc))
    elif color_name == ALPHA_FIELD:
        # alpha=0 -> 白色, alpha=1 -> 蓝色
        lut.RGBPoints = [0.0, 1.0, 1.0, 1.0,
                         1.0, 0.0, 0.2, 1.0]
        lut.ColorSpace = "RGB"
    disp.LookupTable = lut

    # ---- 5. 颜色范围: 固定 / 全时段预扫描 / 逐帧 ----
    range_lo = range_hi = None
    if fixed_range is not None:
        range_lo, range_hi = fixed_range
        print("[FDM 06] 固定色标范围: %.6g .. %.6g" % (range_lo, range_hi))
    elif color_name == ALPHA_FIELD:
        range_lo, range_hi = 0.0, 1.0
    elif args.color_range_mode == "current":
        print("[FDM 06] 色标范围模式: current (逐帧自适应, 视频会闪烁)")
    else:
        print("[FDM 06] 预扫描色标范围 (%d 帧)..." % n_steps)
        tk = GetTimeKeeper()
        scan_lo, scan_hi = float("inf"), float("-inf")
        for t in timesteps:
            try:
                tk.Time = t
                reader.UpdatePipeline(t)
                es.UpdatePipeline(t)
                if calc is not None:
                    calc.UpdatePipeline(t)
                r = read_array_range(src, color_name, assoc)
            except Exception as exc:
                sys.stderr.write("[FDM WARN] t=%s 范围扫描失败: %s\n" % (t, exc))
                continue
            scan_lo = min(scan_lo, r[0])
            scan_hi = max(scan_hi, r[1])
        if scan_lo == float("inf"):
            sys.stderr.write("[FDM ERR] 色标范围扫描未获得任何数据\n")
            return 2
        # 幅值场下界贴 0 (mag>=0), 避免插值负噪声把色标拖偏
        if mag_mode:
            scan_lo = 0.0
        range_lo, range_hi = scan_lo, scan_hi
        print("[FDM 06] 全时段色标范围: %.6g .. %.6g" % (range_lo, range_hi))
    if range_lo is not None:
        lut.RescaleTransferFunction(range_lo, range_hi)

    # ---- 6. 标量条 (颜色图例) ----
    show_bar = args.scalar_bar
    if show_bar is None:
        show_bar = color_name != ALPHA_FIELD
    if show_bar:
        bar_title = args.scalar_bar_title
        if bar_title is None:
            bar_title = "Velocity Magnitude (m/s)" if mag_mode else color_name
        bar = GetScalarBar(lut, view)
        bar.Visibility = 1
        bar.Title = bar_title
        bar.ComponentTitle = ""
        bar.WindowLocation = args.scalar_bar_location
        bar.LabelFormat = args.scalar_bar_format
        bar.TitleColor = [0, 0, 0]
        bar.LabelColor = [0, 0, 0]
        print("[FDM 06] 标量条: '%s' @ %s" % (bar_title, args.scalar_bar_location))

    # 注: 不再叠加独立 Wireframe representation (旧 3b 段已删除)。独立 Wireframe 会画
    # 所有 z 层 hex 的三角剖分/体斜边, 在准 2D 多 z 层投影下重叠成密集"斜三角线"。
    # "Surface With Edges" 只画 ExtractSurface 外表面边；若软件渲染仍暴露 AMR
    # 表面三角化/悬挂节点边, --no-wireframe 切到纯 Surface 作为严格无斜线路径。

    # ---- 7. 相机: 正对 z 平面 ----
    # 全流域: x [0, 0.584], y [0, 0.365], z [0, 0.012] (准二维, z 为薄方向)
    # 梁区:   x=[0.292,0.304] y=[0,0.120] 的中央悬臂梁 (12 mm 宽)
    view.CameraParallelProjection = 1
    view.CameraPosition = [cam_focus[0], cam_focus[1], 1.0]
    view.CameraFocalPoint = list(cam_focus)
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = cam_scale
    print("[FDM 06] 相机: focus=(%.4g,%.4g,%.4g) scale=%.4g view=%.3gx%.3g m"
          % (cam_focus[0], cam_focus[1], cam_focus[2], cam_scale,
             2.0 * cam_scale * W / float(H), 2.0 * cam_scale))

    # ---- 8. 逐时间步推进 ----
    # 必须显式给每个 source 传时间 t: 否则客户端数据访问可能停留在旧时步,
    # 颜色范围/画面刷新不可靠 (见文件头经验)。
    tk = GetTimeKeeper()
    for i, t in enumerate(timesteps):
        tk.Time = t
        reader.UpdatePipeline(t)
        es.UpdatePipeline(t)
        if calc is not None:
            calc.UpdatePipeline(t)
        if (fixed_range is None and args.color_range_mode == "current"
                and color_name != ALPHA_FIELD):
            r = read_array_range(src, color_name, assoc)
            lut.RescaleTransferFunction(r[0], r[1])
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
