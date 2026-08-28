#!/usr/bin/env python3
"""05_postprocess.py - 重建并后处理 FDM 算例

数据源: fastDynamicFvMesh + myInterFoam 的 FDM 算例 (fdm_case):
  - 基线网格:  <case>/constant/polyMesh/points
  - 写盘时步:  <case>/<t>/polyMesh/points   (t = 0.05, 0.1, ..., 0.5, 共10个,
               动态 mesh 每时步独立保存 points)
位移 = points(t)[node] - baseline[node], 其中 node 是挡板右上角节点,
通过坐标近似匹配定位 (默认 (0.304, 0.08, 0.0), 实测为节点 #40440, 0-indexed)。

输出 CSV: <out-dir>/fdm_gate_topright_disp.csv, 表头 Time,dx,dy,dz,
每行一个写盘时步。默认 --out-dir = <case 的父目录>/postprocess (即算例根目录,
如 --case 传 fdm_case 时 = 其父目录下的 postprocess, 自动创建)。

性能报告 (默认开启, --no-perf 关闭): 从算例日志提取运行性能数据, 输出
<out-dir>/fdm_perf_report.md (markdown, 含 markdown 表格)。数据源:
  - fluid 日志:  <case>/log.interFoam* / <case>/log.myInterFoam* (glob 取最新,
    可用 --fluid-log 显式指定), 提取 Exec 行 / ExecutionTime / ClockTime / 末 Time
                 + FSI timing cumulative (区分 流体求解 vs 动网格+模态 耗时)
  - checkMesh:   <case>/log.checkMesh* (glob 取最新), 提取 cells/points
  - 核数:        Exec 是否含 -parallel + <case>/processor* 目录计数
  - 流水线步骤:  关键产物文件 mtime (mesh/ / modal_ve/ / system/controlDict /
                 <endTime>/polyMesh/points / frames/ / postprocess/*.csv),
                 相邻 mtime 之差估算各步骤耗时 (仅供参考)
所有数值均从日志读取, 读不到打 N/A 并在报告中注明缺失来源, 不中断主流程。

用法:
  python3 05_postprocess.py --case <fdm_case_dir> \
      [--out-dir <dir>] [--node-x 0.304] [--node-y 0.08] [--node-z 0.0] \
      [--res <N>] [--node-tol 0.001] [--perf|--no-perf] [--md]
      [--fluid-log <path>]

注意:
  - FDM 并行运行存在 processor*/ 子目录, 但流水线已 reconstruct 到顶层;
    本脚本优先读顶层 <t>/polyMesh/points, 顶层缺失但 processor0 存在时
    warn 并跳过该时步。
  - AMR 改变点数/编号时, 从该时步 FSI_FLUID patch 的 z=0 前缘点中定位
    顶部带内 x 最大点, 不在全流体网格中做几何回退。
  - 仅用 python3 标准库解析 OpenFOAM 文本 (re + 纯文本), 无第三方依赖。
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import tempfile
import time


def parse_vectors(path):
    """解析 OpenFOAM 文本 points/pointVectorField, 返回 (声明数量, [[x,y,z], ...])

    数据段形如:
        <N>
        (
        (x0 y0 z0)
        ...
        )
    用正则抓所有恰好 3 个浮点的括号三元组, 并按声明数量 N 截断, 从而忽略
    header 噪声与 (pointVectorField 的) boundaryField 附加三元组。
    """
    with open(path) as f:
        text = f.read()
    m = re.search(r'(?m)^\s*(\d+)\s*\n\(\s*$', text)
    n = int(m.group(1)) if m else None
    out = []
    for mm in re.finditer(r'\(([^()]*)\)', text):
        parts = mm.group(1).split()
        if len(parts) == 3:
            try:
                out.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except ValueError:
                pass
    if n is not None:
        out = out[:n]
    return n, out


def locate_node(points, nx, ny, nz, tol=1e-4):
    """按坐标近似匹配定位节点索引; 容差内匹配多个时取最近者并 warn"""
    best, best_d2 = -1, float("inf")
    n_match = 0
    for i, p in enumerate(points):
        if abs(p[0] - nx) < tol and abs(p[1] - ny) < tol and abs(p[2] - nz) < tol:
            n_match += 1
        d2 = (p[0] - nx) ** 2 + (p[1] - ny) ** 2 + (p[2] - nz) ** 2
        if d2 < best_d2:
            best, best_d2 = i, d2
    if best < 0:
        return None
    if n_match != 1:
        sys.stderr.write(
            "[FDM WARN] 坐标 (%g, %g, %g) 容差 ±%g 内匹配到 %d 个点, 取最近者 #%d (%g %g %g)\n"
            % (nx, ny, nz, tol, n_match, best,
               points[best][0], points[best][1], points[best][2]))
    return best


def parse_faces(path):
    faces = []
    with open(path, encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = re.match(r'\s*(\d+)\(([^)]*)\)\s*$', line)
            if match:
                point_ids = [int(value) for value in match.group(2).split()]
                if len(point_ids) == int(match.group(1)):
                    faces.append(point_ids)
    return faces


def read_patch_range(path, patch_name):
    text = read_log(path)
    if text is None:
        return None
    match = re.search(
        r'\b%s\s*\{.*?nFaces\s+(\d+);.*?startFace\s+(\d+);'
        % re.escape(patch_name), text, re.S)
    return (int(match.group(2)), int(match.group(1))) if match else None


def locate_fsi_top_right(poly_mesh_dir, points, z_value, top_band):
    patch_range = read_patch_range(
        os.path.join(poly_mesh_dir, "boundary"), "FSI_FLUID")
    faces_path = os.path.join(poly_mesh_dir, "faces")
    if patch_range is None or not os.path.isfile(faces_path):
        return None
    faces = parse_faces(faces_path)
    start, count = patch_range
    if start + count > len(faces):
        return None
    patch_points = {
        point_id
        for face in faces[start:start + count]
        for point_id in face
        if point_id < len(points)
    }
    front = [
        point_id for point_id in patch_points
        if abs(points[point_id][2] - z_value) <= 1e-8
    ]
    if not front:
        return None
    top_y = max(points[point_id][1] for point_id in front)
    top = [
        point_id for point_id in front
        if points[point_id][1] >= top_y - top_band
    ]
    return max(top, key=lambda point_id: points[point_id][0])


def locate_openfoam_bashrc():
    candidates = [
        "/usr/lib/openfoam/openfoam2412/etc/bashrc",
        "/opt/openfoam2412/etc/bashrc",
    ]
    return next((path for path in candidates if os.path.isfile(path)), None)


def is_amr_case(case):
    path = os.path.join(case, "constant", "dynamicMeshDict")
    text = read_log(path)
    return text is not None and "dynamicRefineFvMeshCoeffs" in text


def numeric_time_names(directory):
    names = set()
    for path in glob.glob(os.path.join(directory, "[0-9]*")):
        if not os.path.isdir(path):
            continue
        try:
            float(os.path.basename(path))
        except ValueError:
            continue
        names.add(os.path.basename(path))
    return names


def needs_reconstruction(case):
    processor_times = numeric_time_names(os.path.join(case, "processor0"))
    return bool(processor_times - numeric_time_names(case))


def run_foam_command(case, bashrc, command, log_name):
    shell = 'source "%s" >/dev/null 2>&1 && %s' % (bashrc, command)
    with open(os.path.join(case, log_name), "w", encoding="utf-8") as stream:
        result = subprocess.run(
            ["bash", "-c", shell], cwd=case, stdout=stream,
            stderr=subprocess.STDOUT, check=False)
    return result.returncode


def reconstruct_parallel_results(case):
    if not needs_reconstruction(case):
        print("[FDM 05] 无待重建并行时步, 跳过 reconstruct")
        return True
    bashrc = locate_openfoam_bashrc()
    if bashrc is None:
        sys.stderr.write("[FDM ERR] 找不到 OpenFOAM bashrc, 无法重建并行结果\n")
        return False
    if is_amr_case(case):
        print("[FDM 05] AMR 并行网格重建: reconstructParMesh")
        if run_foam_command(
                case, bashrc, "reconstructParMesh",
                "log.reconstructParMesh") != 0:
            sys.stderr.write(
                "[FDM ERR] reconstructParMesh 失败, 详见 log.reconstructParMesh\n")
            return False
    print("[FDM 05] 并行场重建: reconstructPar")
    if run_foam_command(
            case, bashrc, "reconstructPar",
            "log.reconstructPar") != 0:
        sys.stderr.write("[FDM ERR] reconstructPar 失败, 详见 log.reconstructPar\n")
        return False
    print("[FDM 05] 并行结果重建完成")
    return True


# ---- 性能报告: 日志读取辅助 (纯标准库, 全部从日志/算例文件读取, 不硬编码) ----

def read_log(filepath):
    """读取文本文件内容; 打不开返回 None (不抛异常)"""
    try:
        with open(filepath, encoding="utf-8", errors="replace") as f:
            return f.read()
    except OSError:
        return None


def find_newest_log(patterns):
    """glob 多个 pattern, 返回 mtime 最新者; 一个都匹配不到返回 None"""
    cands = []
    for pat in patterns:
        cands.extend(glob.glob(pat))
    if not cands:
        return None
    cands.sort(key=os.path.getmtime, reverse=True)
    return cands[0]


def extract_exec_time(text):
    """最后一个 'ExecutionTime = X s  ClockTime = Y s', 返回原始字符串 (保持日志原样);
    未找到返回 (None, None)"""
    ms = re.findall(
        r'ExecutionTime\s*=\s*([\d.]+)\s*s\s+ClockTime\s*=\s*([\d.]+)\s*s', text)
    if ms:
        return ms[-1][0], ms[-1][1]
    return None, None


def extract_count(text, pattern):
    """正则提取整数 (容忍千分位逗号), 失败返回 None"""
    m = re.search(pattern, text)
    if not m:
        return None
    return int(m.group(1).replace(",", ""))


def extract_cells(checkmesh_text):
    """checkMesh 输出 'cells: N' -> 单元数"""
    return extract_count(checkmesh_text, r'cells:\s*([\d,]+)')


def extract_points(checkmesh_text):
    """checkMesh 输出 'points: N' -> 节点数"""
    return extract_count(checkmesh_text, r'points:\s*([\d,]+)')


def extract_last_time(text):
    """日志中最后出现的行首 'Time = X' -> 实际推进到的物理时间 (原始字符串)"""
    ts = re.findall(r'(?m)^Time\s*=\s*([\d.]+)', text)
    return ts[-1] if ts else None


def extract_exec_line(text):
    """'Exec : ...' 行 -> 求解器命令行"""
    m = re.search(r'(?m)^\s*Exec\s*:\s*(.+)$', text)
    return m.group(1).strip() if m else None


def extract_control_key(text, key):
    """controlDict 中 '<key> <value>;' -> float"""
    m = re.search(r'(?m)^\s*%s\s+([\d.eE+-]+)\s*;' % re.escape(key), text)
    return float(m.group(1)) if m else None


def fmt_num(v):
    """None -> 'N/A', 数值 -> 紧凑格式 (%.6g, 去尾零)"""
    if v is None:
        return "N/A"
    s = "%.6g" % v
    return s.rstrip("0").rstrip(".") if "." in s else s


def get_mtime(path):
    """文件/目录 mtime (epoch 秒); 不存在返回 None"""
    try:
        return os.path.getmtime(path)
    except OSError:
        return None


def min_mtime(patterns):
    """一组候选参照物 (glob 模式) 的最小 mtime; 全部缺失返回 (None, None).
    取最小可规避被后续流程改写导致的 mtime 膨胀 (如 controlDict 在运行后被恢复)."""
    best_mt, best_path = None, None
    for pat in patterns:
        for p in glob.glob(pat):
            mt = get_mtime(p)
            if mt is not None and (best_mt is None or mt < best_mt):
                best_mt, best_path = mt, p
    return best_mt, best_path


def fmt_ts(mt):
    """mtime -> 'YYYY-MM-DD HH:MM'; None -> 'N/A'"""
    if mt is None:
        return "N/A"
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(mt))


def fmt_dur(secs):
    """秒 -> 人类可读时长; None -> 'N/A'"""
    if secs is None:
        return "N/A"
    secs = int(round(secs))
    if secs < 60:
        return "%d s" % secs
    if secs < 3600:
        return "%d m %02d s" % (secs // 60, secs % 60)
    if secs < 86400:
        return "%d h %02d m" % (secs // 3600, (secs % 3600) // 60)
    return "%d d %d h" % (secs // 86400, (secs % 86400) // 3600)


def md_cell(v):
    """markdown 表格单元格转义 (竖线转 \\|, None -> 'N/A')"""
    if v is None:
        return "N/A"
    return str(v).replace("|", "\\|")


def extract_fsi_timing(text):
    """fluid 日志中最后一条 'FSI timing' cumulative -> (fluid, fluid_pct,
    mesh, mesh_pct, total); 未找到返回 None"""
    ms = re.findall(
        r'cumulative\{fluid=([\d.]+)\s*\(([\d.]+)%\),\s*mesh=([\d.]+)\s*\(([\d.]+)%\),\s*total=([\d.]+)\}',
        text)
    if ms:
        f, fp, m, mp, tot = ms[-1]
        return float(f), float(fp), float(m), float(mp), float(tot)
    return None


def step_durations(markers):
    """[(no, name, ref_path, mtime), ...] -> [(no, name, ref_path, mtime,
    耗时秒|None)]; 相邻 mtime 之差即步骤耗时估算, 首步/缺失/乱序 -> None"""
    out = []
    prev = None
    for no, name, ref, mt in markers:
        dur = None
        if prev is not None and mt is not None:
            dur = mt - prev
            if dur < 0:
                dur = None  # 时间戳乱序 (如跨天产物被后续改写)
        out.append((no, name, ref, mt, dur))
        if mt is not None:
            prev = mt
    return out


def collect_fdm_step_markers(case, out_dir, end_time):
    """FDM 流水线各步骤参照产物 mtime; 返回 [(no, name, ref_path, mtime), ...]

    步骤参照物 (能拿到就报, 拿不到 mtime=None):
      1 网格 = mesh/fluid-openfoam/constant/polyMesh/points
      2 模态 = modal_ve/mode/FluidNodeCoor.csv
      3 组装 = fdm_case/system/controlDict
      4 运行 = fdm_case/processor0/<endTime>/polyMesh/points (并行) 或顶层时步
      5 后处理 = postprocess/fdm_gate_topright_disp.csv
      6 视频 = frames/frame_*.png
    """
    root = os.path.dirname(os.path.normpath(case))
    mesh_pts = os.path.join(root, "mesh", "fluid-openfoam", "constant",
                            "polyMesh", "points")
    modal_csv = os.path.join(root, "modal_ve", "mode", "FluidNodeCoor.csv")
    # 组装产物: 取 controlDict/decomposeParDict/fvSchemes 的最小 mtime;
    # controlDict 在运行结束后被 04 脚本恢复原值而 mtime 变新, 单独用它会失真.
    ctrl = min_mtime([os.path.join(case, "system", "controlDict"),
                      os.path.join(case, "system", "decomposeParDict"),
                      os.path.join(case, "system", "fvSchemes")])
    if end_time is not None:
        run_pts = os.path.join(case, fmt_num(end_time), "polyMesh", "points")
    else:
        # 未知 endTime 时取最新数值时步的 polyMesh/points
        tdirs = [d for d in glob.glob(os.path.join(case, "[0-9]*"))
                 if os.path.isdir(d)]
        run_pts = (os.path.join(max(tdirs, key=lambda d: float(os.path.basename(d))),
                                "polyMesh", "points")
                   if tdirs else None)
    # 并行运行以 processor0 写盘时刻为准，避免把步骤 05 重建时间误算进步骤 04。
    processor_pts = (os.path.join(case, "processor0", fmt_num(end_time),
                                  "polyMesh", "points")
                     if end_time is not None else None)
    if processor_pts and get_mtime(processor_pts) is not None:
        run_mt, run_ref = get_mtime(processor_pts), processor_pts
    else:
        run_mt, run_ref = (get_mtime(run_pts), run_pts) if run_pts else (None, None)
    if run_mt is None:
        run_mt, run_ref = min_mtime([os.path.join(case, "log.interFoam*"),
                                     os.path.join(case, "log.myInterFoam*")])
    frame_png = os.path.join(root, "frames", "frame_*.png")
    post_csv = os.path.join(out_dir, "fdm_gate_topright_disp.csv")

    markers = [
        (1, "网格生成", mesh_pts, get_mtime(mesh_pts)),
        (2, "模态分析", modal_csv, get_mtime(modal_csv)),
        (3, "算例组装", ctrl[1], ctrl[0]),
        (4, "FSI 运行", run_ref, run_mt),
        (5, "并行重建/后处理", post_csv, get_mtime(post_csv)),
        (6, "视频导出", frame_png, min_mtime([frame_png])[0]),
    ]
    return markers


def write_perf_report(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("[FDM 05] 性能报告: %s" % path)


def build_fdm_perf_report(case, out_dir, csv_path, rows, fluid_log_override):
    """从算例日志提取性能数据并写 <out_dir>/fdm_perf_report.md (markdown)。

    计算口径 (4.2):
      总墙钟      = 日志最后一个 ClockTime
      总 CPU      = 日志最后一个 ExecutionTime
      实际推进物理时间 = 日志最后 'Time = X' (完整跑完时 = endTime)
      单位物理时间耗时 = 总墙钟 / 实际推进物理时间
      单元吞吐率  = 单元数 / 总墙钟 (cells/s)
      求解器内部分阶段 (FDM): 解析最后一条 'FSI timing' cumulative,
        区分 流体求解 与 动网格+模态 (FSI) 的 CPU 耗时 (精确)
      流水线步骤耗时: 关键产物文件 mtime 之差估算 (仅供参考)
    所有数值从日志读取; 读不到填 N/A 并记入备注, 不中断位移 CSV 主流程。
    """
    notes = []

    # ---- 日志定位 (glob 推断, 不依赖绝对路径) ----
    fluid_path = fluid_log_override
    if not fluid_path:
        fluid_path = find_newest_log([
            os.path.join(case, "log.interFoam*"),
            os.path.join(case, "log.myInterFoam*"),
        ])
    checkmesh_path = find_newest_log([os.path.join(case, "log.checkMesh*")])

    fluid_text = read_log(fluid_path) if fluid_path else None
    checkmesh_text = read_log(checkmesh_path) if checkmesh_path else None

    if fluid_text is None:
        notes.append("未找到 fluid 日志 (glob <case>/log.interFoam*/log.myInterFoam*, "
                     "或 --fluid-log)")
    if checkmesh_text is None:
        notes.append("未找到 checkMesh 日志 (<case>/log.checkMesh*)")

    # ---- 网格规模: 优先 checkMesh cells/points ----
    cells = extract_cells(checkmesh_text) if checkmesh_text else None
    points = extract_points(checkmesh_text) if checkmesh_text else None

    # ---- 求解器 / 并行核数: Exec 行 + processor* 目录 ----
    exec_line = extract_exec_line(fluid_text) if fluid_text else None
    proc_dirs = [d for d in glob.glob(os.path.join(case, "processor*"))
                 if os.path.isdir(d)]
    if exec_line and "-parallel" in exec_line:
        nprocs = len(proc_dirs) if proc_dirs else None
        if nprocs is None:
            notes.append("Exec 含 -parallel 但 <case>/processor* 目录缺失")
    else:
        nprocs = 1 if exec_line else None

    # ---- 时间: ExecutionTime/ClockTime 最终值 + 末 Time + controlDict 口径 ----
    exec_time, clock_time = (extract_exec_time(fluid_text)
                             if fluid_text else (None, None))
    last_time = extract_last_time(fluid_text) if fluid_text else None

    control_text = read_log(os.path.join(case, "system", "controlDict"))
    end_time = extract_control_key(control_text, "endTime") if control_text else None
    delta_t = extract_control_key(control_text, "deltaT") if control_text else None

    # 实际推进物理时间 = 日志末 Time (读不到时退化为 endTime, 并备注)
    if last_time is not None:
        advanced = float(last_time)
    elif end_time is not None:
        advanced = end_time
        notes.append("fluid 日志无 'Time =' 记录, 按 endTime %.6g 计" % end_time)
    else:
        advanced = None

    # ---- 求解器内部分阶段: FSI timing cumulative (精确) ----
    fsi = extract_fsi_timing(fluid_text) if fluid_text else None
    if fsi is None and fluid_text is not None:
        notes.append("fluid 日志无 'FSI timing' cumulative 行, 无法区分流体/动网格耗时")

    # ---- 派生指标 ----
    if clock_time is not None and advanced:
        unit_time = float(clock_time) / advanced      # s/1s 物理时间
    else:
        unit_time = None
    if cells is not None and clock_time is not None:
        throughput = cells / float(clock_time)        # cells/s
    else:
        throughput = None

    # ---- 完整度 ----
    if last_time is not None:
        if end_time is not None and float(last_time) >= end_time - 1e-12:
            completeness = "完整跑完 (末 Time = %s, endTime %.6g)" % (last_time, end_time)
        else:
            completeness = ("未跑完: 于 %s 结束 (endTime %s 未到)"
                            % (last_time, fmt_num(end_time)))
    else:
        completeness = "N/A (日志无 'Time =' 记录)"

    # ---- 位移峰值 (从刚写出的位移 CSV 数据读, 可选行) ----
    disp_peak = ""
    if rows:
        t_peak, dx_peak = max(((t, abs(dx)) for t, dx, dy, dz in rows),
                              key=lambda td: td[1])
        disp_peak = ("|dx|max = %.6g (t = %s, dx = %g)"
                     % (dx_peak, fmt_num(t_peak), dx_peak))
    else:
        disp_peak = "N/A (无有效位移时步)"

    # ---- 流水线步骤耗时 (估算, 基于产物 mtime) ----
    markers = collect_fdm_step_markers(case, out_dir, end_time)
    step_rows = step_durations(markers)

    # ---- 组装 markdown 报告 ----
    lines = []
    lines.append("# [FDM] 性能分析报告")
    lines.append("")
    lines.append("> 生成时间: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    lines.append("")
    lines.append("## 一、算例概况")
    lines.append("")
    lines.append("| 项 | 值 |")
    lines.append("|---|---|")
    lines.append("| 算例目录 | %s |" % md_cell(case))
    lines.append("| 网格规模 | %s 单元" % fmt_num(cells)
                 + ("  (%s 节点)" % fmt_num(points) if points is not None else
                    "  (节点数 N/A, log.checkMesh 无 points 行)")
                 + " |")
    lines.append("| 物理时间 | %s s" % fmt_num(advanced)
                 + ("" if end_time is None else "  (endTime %s s)" % fmt_num(end_time))
                 + " |")
    lines.append("| 求解器 | %s |" % md_cell(exec_line if exec_line else "N/A (日志无 Exec 行)"))
    if nprocs is None:
        nproc_cell = "N/A"
    elif nprocs > 1:
        nproc_cell = "%d 核 (-parallel, processor* x %d)" % (nprocs, len(proc_dirs))
    else:
        nproc_cell = "1 核 (串行)"
    lines.append("| 并行核数 | %s |" % nproc_cell)
    lines.append("| 时间步长 | %s s |" % fmt_num(delta_t))
    lines.append("| 完整度 | %s |" % md_cell(completeness))
    lines.append("")
    lines.append("## 二、求解器内部分阶段耗时 (FDM)")
    lines.append("")
    lines.append("| 阶段 | CPU 秒 | 占比 |")
    lines.append("|---|---|---|")
    if fsi:
        f, fp, m, mp, tot = fsi
        lines.append("| 流体求解 (N-S/VOF) | %.2f | %.1f%% |" % (f, fp))
        lines.append("| 动网格+模态 (FSI) | %.2f | %.1f%% |" % (m, mp))
        lines.append("| **合计** | **%.2f** | 100%% |" % tot)
        lines.append("")
        lines.append("> 数据源: %s 的 FSI timing cumulative (最后一条), 精确计时"
                     % md_cell(os.path.basename(fluid_path) if fluid_path else "fluid 日志"))
    else:
        lines.append("| 流体求解 (N-S/VOF) | N/A | N/A |")
        lines.append("| 动网格+模态 (FSI) | N/A | N/A |")
        lines.append("| **合计** | **N/A** | N/A |")
        lines.append("")
        lines.append("> 数据源: fluid 日志无 FSI timing cumulative, 无法区分")
    lines.append("")
    lines.append("## 三、流水线步骤耗时 (估算)")
    lines.append("")
    lines.append("| 步骤 | 内容 | 参照产物 | 完成时刻 | 耗时(估) |")
    lines.append("|---|---|---|---|---|")
    for no, name, ref, mt, dur in step_rows:
        lines.append("| %d | %s | %s | %s | %s |"
                     % (no, name, md_cell(os.path.basename(ref) if ref else "N/A"),
                        fmt_ts(mt), fmt_dur(dur)))
    lines.append("")
    lines.append("> 基于产物文件 mtime 差估算, 仅供参考; 步骤耗时含步骤间的人工/调度间隔,"
                 " 非纯计算耗时。")
    lines.append("")
    lines.append("## 四、总体指标")
    lines.append("")
    lines.append("| 指标 | 值 |")
    lines.append("|---|---|")
    lines.append("| 总 CPU 累计 | %s s |" % (exec_time if exec_time else "N/A"))
    lines.append("| 总墙钟 | %s s |" % (clock_time if clock_time else "N/A"))
    lines.append("| 单位物理时间 | %s s/1s 物理时间 |" % fmt_num(unit_time))
    lines.append("| 单元吞吐率 | %s cells/s |" % fmt_num(throughput))
    lines.append("| 并行效率 | N/A (无串行基线数据, 无法计算加速比; 并行核数 %s) |"
                 % fmt_num(nprocs))
    lines.append("| 位移峰值 | %s |" % md_cell(disp_peak))
    lines.append("")
    lines.append("## 五、备注")
    lines.append("")
    if notes:
        for n in notes:
            lines.append("- %s" % n)
    else:
        lines.append("- 无")
    lines.append("")
    lines.append("- 求解器内部分阶段耗时为精确累计 (FSI timing cumulative);"
                 " 流水线步骤耗时为产物时间戳估算, 仅供参考。")

    # ---- 写报告 (失败仅 warn, 不影响位移 CSV 主流程/返回码) ----
    try:
        write_perf_report(os.path.join(out_dir, "fdm_perf_report.md"), lines)
    except OSError as e:
        sys.stderr.write("[FDM WARN] 写性能报告失败: %s\n" % e)


def write_test_points(path, points):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    lines = [str(len(points)), "("]
    lines.extend("(%s %s %s)" % point for point in points)
    lines.append(")")
    with open(path, "w", encoding="utf-8") as stream:
        stream.write("\n".join(lines) + "\n")


def write_test_fsi_patch(poly_mesh, face_point_ids):
    faces = ["3(%d %d %d)" % ids for ids in face_point_ids]
    with open(os.path.join(poly_mesh, "faces"), "w", encoding="utf-8") as stream:
        stream.write("%d\n(\n%s\n)\n" % (len(faces), "\n".join(faces)))
    with open(os.path.join(poly_mesh, "boundary"), "w", encoding="utf-8") as stream:
        stream.write(
            "1\n(\nFSI_FLUID\n{\n nFaces %d;\n startFace 0;\n}\n)\n"
            % len(faces))


def run_self_test():
    with tempfile.TemporaryDirectory() as tmp:
        case = os.path.join(tmp, "fdm_case")
        out = os.path.join(tmp, "postprocess")
        write_test_points(
            os.path.join(case, "constant", "polyMesh", "points"),
            [(0.0, 0.0, 0.0), (0.304, 0.08, 0.0), (0.5, 0.3, 0.0)])
        write_test_points(
            os.path.join(case, "0.01", "polyMesh", "points"),
            [(0.0, 0.0, 0.0), (0.3041, 0.0801, 0.0),
             (0.302, 0.081, 0.0), (0.306, 0.081, 0.0),
             (0.306, 0.081, 0.012)])
        write_test_fsi_patch(
            os.path.join(case, "0.01", "polyMesh"), [(2, 3, 4)])
        if main(["--case", case, "--out-dir", out, "--no-reconstruct"]) != 0:
            return 1
        csv_path = os.path.join(out, "fdm_gate_topright_disp.csv")
        with open(csv_path, newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        if (len(rows) != 1 or abs(float(rows[0]["dx"]) - 0.002) > 1e-12
                or abs(float(rows[0]["dy"]) - 0.001) > 1e-12):
            return 1
        case2 = os.path.join(tmp, "no_rows", "fdm_case")
        out2 = os.path.join(tmp, "no_rows", "postprocess")
        write_test_points(
            os.path.join(case2, "constant", "polyMesh", "points"),
            [(0.0, 0.0, 0.0), (0.304, 0.08, 0.0)])
        write_test_points(
            os.path.join(case2, "0.01", "polyMesh", "points"),
            [(0.0, 0.0, 0.0), (0.1, 0.1, 0.0), (0.2, 0.2, 0.0)])
        if main(["--case", case2, "--out-dir", out2, "--no-reconstruct"]) != 0:
            return 1
        if os.path.exists(os.path.join(out2, "fdm_gate_topright_disp.csv")):
            return 1
        report = read_log(os.path.join(out2, "fdm_perf_report.md")) or ""
        if "N/A (无有效位移时步)" not in report:
            return 1
    print("[FDM 05] SELF-TEST PASS")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="提取 FDM 挡板右上角位移时间序列并输出 CSV"
    )
    ap.add_argument("--case", required=False,
                    help="FDM 算例目录 (fdm_case 根, 含 constant/polyMesh/points 与 <t>/polyMesh/points)")
    ap.add_argument("--out-dir", default=None,
                    help="输出 CSV 目录 (默认 <case 的父目录>/postprocess, 即算例根目录, 自动创建)")
    ap.add_argument("--node-x", type=float, default=0.304,
                    help="挡板右上角 x (默认 0.304)")
    ap.add_argument("--node-y", type=float, default=0.080,
                    help="挡板右上角 y (默认 0.080)")
    ap.add_argument("--node-z", type=float, default=0.0,
                    help="挡板右上角 z (默认 0.0)")
    ap.add_argument("--res", type=int, default=None,
                    help="预留: 期望网格点数, 与实测不符时 warn (默认不校验)")
    ap.add_argument("--node-tol", type=float, default=5e-4,
                    help="AMR 后 FSI patch 顶部候选带厚度 m (默认 0.0005)")
    ap.add_argument("--perf", dest="perf", action="store_true", default=None,
                    help="生成性能分析报告 (默认开启)")
    ap.add_argument("--no-perf", dest="perf", action="store_false",
                    help="跳过性能分析报告, 只输出位移 CSV")
    ap.add_argument("--md", dest="md", action="store_true", default=False,
                    help="性能报告输出为 markdown (.md, 默认, 无需显式指定)")
    ap.add_argument("--fluid-log", default=None,
                    help="fluid 求解器日志路径 (默认 glob <case>/log.interFoam* 等推断)")
    ap.add_argument("--no-reconstruct", action="store_true",
                    help="跳过 processor* 并行结果重建")
    ap.add_argument("--self-test", action="store_true",
                    help="运行内置 AMR 后处理回归测试并退出")
    args = ap.parse_args(argv)

    if args.self_test:
        return run_self_test()
    if not args.case:
        ap.error("--case is required unless --self-test is used")

    if args.perf is None:
        args.perf = True  # 默认 --perf 开

    print("[FDM-PIPE] ===== 步骤 5: 并行重建与位移后处理 =====")

    case = args.case
    out_dir = args.out_dir if args.out_dir else os.path.join(
        os.path.dirname(os.path.normpath(case)), "postprocess")
    os.makedirs(out_dir, exist_ok=True)

    if not args.no_reconstruct and not reconstruct_parallel_results(case):
        return 2

    # ---- 1. 基线网格 + 节点定位 ----
    base_path = os.path.join(case, "constant", "polyMesh", "points")
    if not os.path.isfile(base_path):
        sys.stderr.write("[FDM ERR] 找不到基线网格 %s\n" % base_path)
        return 2
    n0, base = parse_vectors(base_path)
    print("[FDM 05] 基线网格: %s (%d 点)" % (base_path, len(base)))
    if args.res is not None and len(base) != args.res:
        sys.stderr.write("[FDM WARN] --res %d 与实测网格点数 %d 不符 (仅警告, 不中断)\n"
                         % (args.res, len(base)))

    idx = locate_node(base, args.node_x, args.node_y, args.node_z)
    if idx is None:
        sys.stderr.write("[FDM ERR] 基线网格中找不到挡板右上角 (%g, %g, %g) 附近的节点\n"
                         % (args.node_x, args.node_y, args.node_z))
        return 2
    print("[FDM 05] 挡板右上角节点 #%d: (%g, %g, %g)"
          % (idx, base[idx][0], base[idx][1], base[idx][2]))

    # ---- 2. 收集顶层数值时步目录 ----
    time_dirs = []
    for d in glob.glob(os.path.join(case, "[0-9]*")):
        try:
            t = float(os.path.basename(d))
        except ValueError:
            continue
        if os.path.isdir(d):
            time_dirs.append((t, d))
    time_dirs.sort()
    if not time_dirs:
        sys.stderr.write("[FDM ERR] 算例 %s 下找不到任何数值时步目录 (先跑步骤 3/4)\n" % case)
        return 2

    # ---- 3. 逐时步提取右上角位移 ----
    rows = []
    skipped = []
    for t, d in time_dirs:
        pts_path = os.path.join(d, "polyMesh", "points")
        if not os.path.isfile(pts_path):
            skipped.append(t)
            proc_pts = os.path.join(case, "processor0", os.path.basename(d),
                                    "polyMesh", "points")
            if os.path.isfile(proc_pts):
                sys.stderr.write(
                    "[FDM WARN] 顶层 %s 缺失但 processor0 存在 (疑似未 reconstruct), 跳过 t=%g\n"
                    % (pts_path, t))
            else:
                sys.stderr.write(
                    "[FDM WARN] 顶层 %s 缺失 (t=%g 无网格, 初始态网格在 constant/), 跳过\n"
                    % (pts_path, t))
            continue
        n, pts = parse_vectors(pts_path)
        amr_changed = n0 is not None and n is not None and n != n0
        if amr_changed:
            current_idx = locate_fsi_top_right(
                os.path.join(d, "polyMesh"), pts, args.node_z, args.node_tol)
            if current_idx is None:
                sys.stderr.write(
                    "[FDM WARN] t=%g 无法从 FSI_FLUID patch 定位右上角, 跳过\n" % t)
                skipped.append(t)
                continue
        else:
            current_idx = idx
        current = pts[current_idx]
        dx = current[0] - base[idx][0]
        dy = current[1] - base[idx][1]
        dz = current[2] - base[idx][2]
        rows.append((t, dx, dy, dz))
        print("[FDM 05] t=%-6g node=#%d dx=%.8g dy=%.8g dz=%.8g"
              % (t, current_idx, dx, dy, dz))

    # ---- 4. 写 CSV ----
    csv_path = os.path.join(out_dir, "fdm_gate_topright_disp.csv")
    if rows:
        with open(csv_path, "w") as f:
            f.write("Time,dx,dy,dz\n")
            for t, dx, dy, dz in rows:
                f.write("%.8g,%.8g,%.8g,%.8g\n" % (t, dx, dy, dz))
        print("[FDM 05] 写入 %s (%d 行)" % (csv_path, len(rows)))
    else:
        sys.stderr.write("[FDM WARN] 没有任何有效位移时步, 性能报告将标记为 N/A\n")
    if skipped:
        print("[FDM 05] 跳过 %d 个时步: %s"
              % (len(skipped), ", ".join("%g" % t for t in skipped)))

    # ---- 5. 性能分析报告 (默认开; 缺失日志仅 warn, 不影响返回码) ----
    if args.perf:
        build_fdm_perf_report(case, out_dir, csv_path, rows, args.fluid_log)
    else:
        print("[FDM 05] --no-perf: 跳过性能报告, 仅输出位移 CSV")
    return 0 if rows or args.perf else 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
