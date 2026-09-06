#!/usr/bin/env python3
"""05_postprocess.py - 提取 preCICE 挡板右上角位移时间序列并输出独立 CSV

数据源: preCICE 流体侧算例 (pc_case/fluid-openfoam):
   - 常量网格:  <pc_case>/fluid-openfoam/constant/polyMesh/points
   - 位移场:    <pc_case>/fluid-openfoam/<t>/pointDisplacement
               (OpenFOAM pointVectorField, internalField 为累计位移,
                按点序排列, 无内嵌坐标)
位移 = pointDisplacement(t)[node], 其中 node 是挡板右上角节点,
通过 constant/polyMesh/points 坐标近似匹配定位 (默认 (0.304, 0.08, 0.0),
实测为节点 #40440, 0-indexed; 两算例网格同源同序)。

输出 CSV: <out-dir>/precice_gate_topright_disp.csv, 表头 Time,dx,dy,dz,
每行一个完整时步。默认 --out-dir = 算例根目录/postprocess
(即 <case 的父目录的父目录>/postprocess, 与 FDM 的 fdm_case 同级, 自动创建)。

性能报告 (默认开启, --no-perf 关闭): 从算例日志提取运行性能数据, 输出
<out-dir>/precice_perf_report.md (markdown, 含 markdown 表格)。数据源
(--case 可为 pc_case 或 fluid-openfoam, 日志自动定位):
  - fluid 日志: <case>/precice-run/log.fluid 或 <case 父目录>/precice-run/log.fluid
    (可用 --fluid-log 显式指定), 提取 Exec 行 / ExecutionTime / ClockTime / 末 Time
  - solid 日志: <case>/precice-run/log.solid 或 <case 父目录>/precice-run/log.solid
    (可用 --solid-log 显式指定), 提取 solid 线程数 / convergence 次数 (耦合迭代) /
    "Adapter calling advance" 次数 (窗口) / "Job finished" (完整度判据)
  - solid .sta: <pc_case>/solid-calculix/*.sta, 提取末 INC 步数 / TOT TIME
  - checkMesh:  <case>/log.checkMesh* 或 <case>/fluid-openfoam/log.checkMesh*, 提取 cells
  - 核数:       fluid Exec 是否含 -parallel + processor* 目录计数
  - 流水线步骤:  关键产物文件 mtime (mesh/ / precice-config.xml /
                 precice-run/log.fluid / 视频帧 / postprocess/*.csv),
                 相邻 mtime 之差估算各步骤耗时 (仅供参考)
所有数值均从日志读取, 读不到打 N/A 并在报告中注明缺失来源, 不中断主流程。

用法:
   python3 05_postprocess.py --case <pc_case> \
      [--out-dir <dir>] [--node-x 0.304] [--node-y 0.08] [--node-z 0.0] \
      [--res <N>] [--perf|--no-perf] [--md] [--fluid-log <path>] [--solid-log <path>]

注意:
  - preCICE 耦合过程会残留只有 uniform/ 的中间时步 (如 0.0095, 无
    pointDisplacement), 本脚本以 pointDisplacement 文件存在与否过滤;
    t=0 位移恒为零 (由常量网格承载基线), 从首个写盘时步开始。
   - 并行运行存在 processor*/ 子目录; stage 5 reconstructPar 后优先读顶层
    <t>/pointDisplacement, 顶层缺失但 processor0 存在时 warn 并跳过。
  - 仅用 python3 标准库解析 OpenFOAM 文本 (re + 纯文本), 无第三方依赖。
"""
import argparse
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
    header 噪声与 pointVectorField 的 boundaryField 附加三元组。
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
            "[PRECICE WARN] 坐标 (%g, %g, %g) 容差 ±%g 内匹配到 %d 个点, 取最近者 #%d (%g %g %g)\n"
            % (nx, ny, nz, tol, n_match, best,
               points[best][0], points[best][1], points[best][2]))
    return best


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


def needs_reconstruction(fluid_dir):
    return bool(numeric_time_names(os.path.join(fluid_dir, "processor0"))
                - numeric_time_names(fluid_dir))


def locate_openfoam_bashrc():
    for path in ("/usr/lib/openfoam/openfoam2412/etc/bashrc",
                 "/opt/openfoam2412/etc/bashrc"):
        if os.path.isfile(path):
            return path
    return None


def reconstruct_parallel_results(fluid_dir):
    if not needs_reconstruction(fluid_dir):
        print("[PRECICE 05] 无待重建并行时步, 跳过 reconstructPar")
        return True
    bashrc = locate_openfoam_bashrc()
    if bashrc is None:
        sys.stderr.write("[PRECICE ERR] 找不到 OpenFOAM bashrc, 无法运行 reconstructPar\n")
        return False
    log_path = os.path.join(fluid_dir, "log.reconstructPar")
    command = 'source "%s" >/dev/null 2>&1 && reconstructPar' % bashrc
    print("[PRECICE 05] 并行场重建: reconstructPar")
    with open(log_path, "w", encoding="utf-8") as stream:
        result = subprocess.run(["bash", "-c", command], cwd=fluid_dir,
                                stdout=stream, stderr=subprocess.STDOUT,
                                check=False)
    if result.returncode:
        sys.stderr.write("[PRECICE ERR] reconstructPar 失败, 详见 %s\n" % log_path)
        return False
    print("[PRECICE 05] 并行结果重建完成")
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
    取最小可规避被后续流程改写导致的 mtime 膨胀 (如 precice-config.xml 在
    运行结束后被 04 脚本恢复原值而 mtime 变新)."""
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


def extract_sta_last(text):
    """.sta 文本最后一行数值 -> (INC 步数, TOT TIME); 无则 (None, None)"""
    last = None
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 5:
            try:
                inc = int(parts[1])
                tot = float(parts[4])
            except ValueError:
                continue
            last = (inc, tot)
    return last


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


def collect_precice_step_markers(case, out_dir):
    """preCICE 流水线各步骤参照产物 mtime; 返回 [(no, name, ref_path, mtime), ...]

    步骤参照物 (能拿到就报, 拿不到 mtime=None; preCICE 无模态步骤):
      1 网格 = mesh/fluid-openfoam/constant/polyMesh/points
      3 组装 = pc_case/precice-config.xml
      4 运行 = pc_case/precice-run/log.fluid (或 log.solid)
      5 后处理 = postprocess/precice_gate_topright_disp.csv
      6 视频 = 帧 PNG
    """
    parent = os.path.dirname(os.path.normpath(case))
    root = os.path.dirname(parent)
    mesh_pts = os.path.join(root, "mesh", "fluid-openfoam", "constant",
                            "polyMesh", "points")
    # 组装产物: 取 precice-config.xml/run-coupled.sh/controlDict 的最小 mtime;
    # precice-config.xml 在运行结束后被 04 脚本恢复 max-time 而 mtime 变新,
    # 单独用它会失真.
    cfg = min_mtime([os.path.join(parent, "precice-config.xml"),
                     os.path.join(parent, "run-coupled.sh"),
                     os.path.join(parent, "fluid-openfoam", "system",
                                  "controlDict")])
    run_log = os.path.join(parent, "precice-run", "log.fluid*")
    video_png = os.path.join(root, "precice_frames", "frame_*.png")
    if not glob.glob(video_png):
        video_png = os.path.join(root, "frames", "frame_*.png")
    post_csv = os.path.join(out_dir, "precice_gate_topright_disp.csv")

    markers = [
        (1, "网格生成", mesh_pts, get_mtime(mesh_pts)),
        (3, "算例组装", cfg[1], cfg[0]),
        (4, "FSI 运行", run_log, min_mtime([run_log])[0]),
        (5, "并行重建/后处理", post_csv, get_mtime(post_csv)),
        (6, "视频导出", video_png, min_mtime([video_png])[0]),
    ]
    return markers


def write_perf_report(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("[PRECICE 05] 性能报告: %s" % path)


def build_precice_perf_report(case, out_dir, csv_path, rows,
                              fluid_log_override, solid_log_override):
    """从算例日志提取性能数据并写 <out_dir>/precice_perf_report.md (markdown)。

    计算口径 (4.2):
      总墙钟      = fluid 日志最后一个 ClockTime
      总 CPU      = fluid 日志最后一个 ExecutionTime
      实际推进物理时间 = fluid 日志最后 'Time = X' (崩溃时 < endTime)
      单位物理时间耗时 = 总墙钟 / 实际推进物理时间
      单元吞吐率  = 单元数 / 总墙钟 (cells/s)
      solid 线程  = "Using up to N cpu(s) for the stress calculation"
      耦合迭代数  = solid 日志 "convergence" 出现次数
      solid 推进步数 = solid 日志 "Adapter calling advance" 出现次数
      solid .sta   = 末 INC 步数 / TOT TIME
      求解器内部分阶段 (preCICE): 流体 interFoam ExecutionTime (精确) vs
        固体 ccx 无独立 CPU 计时, 与流体并行启动/结束 (墙钟≈流体, 近似)
      流水线步骤耗时: 关键产物文件 mtime 之差估算 (仅供参考)
    所有数值从日志读取; 读不到填 N/A 并记入备注, 不中断位移 CSV 主流程。
    """
    notes = []
    parent = os.path.dirname(os.path.normpath(case))

    # ---- 日志定位: 兼容 --case 传 pc_case 或 fluid-openfoam 两种布局 ----
    fluid_path = fluid_log_override
    if not fluid_path:
        fluid_path = find_newest_log([
            os.path.join(case, "precice-run", "log.fluid*"),
            os.path.join(parent, "precice-run", "log.fluid*"),
        ])
    solid_path = solid_log_override
    if not solid_path:
        solid_path = find_newest_log([
            os.path.join(case, "precice-run", "log.solid*"),
            os.path.join(parent, "precice-run", "log.solid*"),
        ])
    checkmesh_path = find_newest_log([
        os.path.join(case, "log.checkMesh*"),
        os.path.join(case, "fluid-openfoam", "log.checkMesh*"),
    ])

    fluid_text = read_log(fluid_path) if fluid_path else None
    solid_text = read_log(solid_path) if solid_path else None
    checkmesh_text = read_log(checkmesh_path) if checkmesh_path else None

    if fluid_text is None:
        notes.append("未找到 fluid 日志 (glob <case>/precice-run/log.fluid*, "
                     "或 --fluid-log)")
    if solid_text is None:
        notes.append("未找到 solid 日志 (glob <case>/precice-run/log.solid*, "
                     "或 --solid-log)")
    if checkmesh_text is None:
        notes.append("未找到 checkMesh 日志 (<case>/log.checkMesh*)")

    # ---- 网格规模: 优先 checkMesh cells ----
    cells = extract_cells(checkmesh_text) if checkmesh_text else None

    # ---- fluid 求解器 / 并行核数: Exec 行 + processor* 目录 ----
    exec_line = extract_exec_line(fluid_text) if fluid_text else None
    proc_dirs = [d for d in glob.glob(os.path.join(case, "processor*"))
                 if os.path.isdir(d)]
    proc_dirs += [d for d in glob.glob(os.path.join(case, "fluid-openfoam",
                                                    "processor*"))
                  if os.path.isdir(d)]
    if exec_line and "-parallel" in exec_line:
        fluid_nprocs = len(proc_dirs) if proc_dirs else None
        if fluid_nprocs is None:
            notes.append("Exec 含 -parallel 但无 processor* 目录")
    else:
        fluid_nprocs = 1 if exec_line else None

    # ---- solid: 线程数 / 耦合迭代 / 窗口步数 / .sta 步数 / 完整度 ----
    solid_threads = None
    convergence_count = None
    advance_count = None
    sta_inc, sta_time = None, None
    if solid_text:
        m = re.search(r'Using up to\s+(\d+)\s+cpu\(s\)\s+for the stress calculation',
                      solid_text)
        solid_threads = int(m.group(1)) if m else None
        convergence_count = solid_text.count("convergence")
        advance_count = solid_text.count("Adapter calling advance")
        if solid_threads is None:
            notes.append("solid 日志无 'Using up to N cpu(s) for the stress "
                         "calculation' 行")
    else:
        notes.append("solid 日志缺失, 无法统计线程/迭代/窗口数")

    # .sta 步数 (ccx 瞬态推进步, 非 preCICE 窗口数)
    sta_path = find_newest_log([
        os.path.join(parent, "solid-calculix", "*.sta"),
        os.path.join(case, "solid-calculix", "*.sta"),
    ])
    sta_text = read_log(sta_path) if sta_path else None
    if sta_text:
        sta_inc, sta_time = extract_sta_last(sta_text)
    else:
        notes.append("未找到 solid .sta (glob <pc_case>/solid-calculix/*.sta)")

    # ---- 时间: ExecutionTime/ClockTime 最终值 + 末 Time + controlDict 口径 ----
    exec_time, clock_time = (extract_exec_time(fluid_text)
                             if fluid_text else (None, None))
    last_time = extract_last_time(fluid_text) if fluid_text else None

    control_text = read_log(os.path.join(case, "system", "controlDict"))
    end_time = extract_control_key(control_text, "endTime") if control_text else None
    delta_t = extract_control_key(control_text, "deltaT") if control_text else None

    if last_time is not None:
        advanced = float(last_time)
    elif end_time is not None:
        advanced = end_time
        notes.append("fluid 日志无 'Time =' 记录, 按 endTime %.6g 计" % end_time)
    else:
        advanced = None

    # ---- 派生指标 ----
    if clock_time is not None and advanced:
        unit_time = float(clock_time) / advanced      # s/1s 物理时间
    else:
        unit_time = None
    if cells is not None and clock_time is not None:
        throughput = cells / float(clock_time)        # cells/s
    else:
        throughput = None

    # ---- 完整度: fluid 末 Time vs endTime; solid "Job finished" + fluid 崩溃 ----
    if last_time is not None:
        if end_time is not None and float(last_time) >= end_time - 1e-12:
            completeness = "完整跑完 (末 Time = %s, endTime %.6g)" % (last_time, end_time)
        else:
            completeness = ("未跑完: 于 %s 结束 (endTime %s 未到)"
                            % (last_time, fmt_num(end_time)))
    else:
        completeness = "N/A (fluid 日志无 'Time =' 记录)"

    solid_finished = solid_text is not None and "Job finished" in solid_text
    fluid_crashed = (fluid_text is not None
                     and ("MPI_ABORT" in fluid_text or "FOAM FATAL" in fluid_text))
    if not solid_finished or fluid_crashed:
        crash_note = "preCICE 耦合中断"
        if fluid_crashed and "adjustPhi" in fluid_text:
            crash_note += " (fluid adjustPhi 崩溃)"
        elif fluid_crashed:
            crash_note += " (fluid 日志含 MPI_ABORT/FOAM FATAL)"
        if not solid_finished and solid_text is not None:
            crash_note += "; solid 无 'Job finished'"
        notes.append(crash_note)

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
    markers = collect_precice_step_markers(case, out_dir)
    step_rows = step_durations(markers)

    # ---- 组装 markdown 报告 ----
    lines = []
    lines.append("# [preCICE] 性能分析报告")
    lines.append("")
    lines.append("> 生成时间: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    lines.append("")
    lines.append("## 一、算例概况")
    lines.append("")
    lines.append("| 项 | 值 |")
    lines.append("|---|---|")
    lines.append("| 算例目录 | %s |" % md_cell(case))
    lines.append("| 网格规模 | %s 单元 |" % fmt_num(cells))
    lines.append("| 物理时间 | %s s" % fmt_num(advanced)
                 + ("" if end_time is None else "  (endTime %s s)" % fmt_num(end_time))
                 + " |")
    lines.append("| 求解器 | %s |" % md_cell(exec_line if exec_line else "N/A (日志无 Exec 行)"))
    if fluid_nprocs is None:
        nproc_cell = "fluid N/A 核 + solid %s 线程" % fmt_num(solid_threads)
    else:
        nproc_cell = ("fluid %d 核 (-parallel, processor* x %d) + solid %s 线程"
                      % (fluid_nprocs, len(proc_dirs), fmt_num(solid_threads)))
    lines.append("| 并行核数 | %s |" % nproc_cell)
    lines.append("| 时间步长 | %s s |" % fmt_num(delta_t))
    lines.append("| 完整度 | %s |" % md_cell(completeness))
    lines.append("")
    lines.append("## 二、求解器内部分阶段耗时 (preCICE)")
    lines.append("")
    lines.append("| 阶段 | 数据 | 说明 |")
    lines.append("|---|---|---|")
    if exec_time:
        lines.append("| 流体求解 (interFoam) | %s s (CPU) | log.fluid 最后 ExecutionTime, 精确 |"
                     % exec_time)
    else:
        lines.append("| 流体求解 (interFoam) | N/A | log.fluid 无 ExecutionTime |")
    solid_note = ("无独立 CPU 计时; 与 fluid 并行启动/结束, 墙钟≈流体"
                  + (" (%.0f s)" % float(clock_time) if clock_time else "")
                  + ("; %s 线程" % fmt_num(solid_threads)
                     if solid_threads is not None else "")
                  + ("; .sta 末 %s 步 (TOT TIME %s)"
                     % (fmt_num(sta_inc), fmt_num(sta_time))
                     if sta_inc is not None else ""))
    lines.append("| 固体求解 (ccx) | 无独立 CPU 数据 | %s |" % solid_note)
    lines.append("")
    lines.append("> 数据源: fluid ExecutionTime 精确; 固体无独立 CPU 计时, 与流体并行运行"
                 ", 墙钟≈流体 (近似)。")
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
    lines.append("| 并行效率 | N/A (无串行基线数据, 无法计算加速比; fluid %s 核) |"
                 % fmt_num(fluid_nprocs))
    lines.append("| solid 线程 | %s |" % fmt_num(solid_threads))
    lines.append("| 耦合迭代总数 | %s 次 (solid 日志 'convergence') |"
                 % fmt_num(convergence_count))
    lines.append("| solid 推进步数 | %s 次 (solid 日志 'Adapter calling advance';"
                 " 指 ccx 瞬态推进步数, 非 preCICE 窗口数) |" % fmt_num(advance_count))
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
    lines.append("- 流体/固体分阶段耗时: 流体为精确 ExecutionTime,"
                 " 固体无独立 CPU 计时, 与流体并行 (墙钟≈流体);"
                 " 流水线步骤耗时为产物时间戳估算, 仅供参考。")

    # ---- 写报告 (失败仅 warn, 不影响位移 CSV 主流程/返回码) ----
    try:
        write_perf_report(os.path.join(out_dir, "precice_perf_report.md"), lines)
    except OSError as e:
        sys.stderr.write("[PRECICE WARN] 写性能报告失败: %s\n" % e)


def write_test_vectors(path, vectors):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as stream:
        stream.write("%d\n(\n%s\n)\n" % (
            len(vectors), "\n".join("(%s %s %s)" % vector for vector in vectors)))


def run_self_test():
    with tempfile.TemporaryDirectory() as tmp:
        case_root = os.path.join(tmp, "pc_case")
        fluid_dir = os.path.join(case_root, "fluid-openfoam")
        out_dir = os.path.join(tmp, "postprocess")
        vectors = [(0, 0, 0), (0.304, 0.08, 0)]
        write_test_vectors(os.path.join(fluid_dir, "constant", "polyMesh", "points"), vectors)
        write_test_vectors(os.path.join(fluid_dir, "processor0", "0.01", "pointDisplacement"),
                           [(0, 0, 0), (0.002, 0.001, 0)])
        bashrc = os.path.join(tmp, "bashrc")
        with open(bashrc, "w", encoding="utf-8") as stream:
            stream.write("reconstructPar() { mkdir -p 0.01; cp processor0/0.01/pointDisplacement 0.01/; }\n")
        original_bashrc = locate_openfoam_bashrc
        try:
            globals()["locate_openfoam_bashrc"] = lambda: bashrc
            if main(["--case", case_root, "--out-dir", out_dir, "--no-perf"]) != 0:
                return 1
        finally:
            globals()["locate_openfoam_bashrc"] = original_bashrc
        csv_path = os.path.join(out_dir, "precice_gate_topright_disp.csv")
        text = read_log(csv_path)
        if (text is None or "0.01,0.002,0.001,0" not in text
                or not os.path.isfile(os.path.join(fluid_dir, "log.reconstructPar"))):
            return 1
    print("[PRECICE 05] SELF-TEST PASS")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="提取 preCICE 挡板右上角位移时间序列并输出 CSV"
    )
    ap.add_argument("--case", required=False,
                    help="assembled preCICE case root (contains fluid-openfoam)")
    ap.add_argument("--out-dir", default=None,
                    help="输出 CSV 目录 (默认 算例根目录/postprocess, 即 <case 的父目录的父目录>/postprocess, 自动创建)")
    ap.add_argument("--node-x", type=float, default=0.304,
                    help="挡板右上角 x (默认 0.304)")
    ap.add_argument("--node-y", type=float, default=0.080,
                    help="挡板右上角 y (默认 0.080)")
    ap.add_argument("--node-z", type=float, default=0.0,
                    help="挡板右上角 z (默认 0.0)")
    ap.add_argument("--res", type=int, default=None,
                    help="预留: 期望网格点数, 与实测不符时 warn (默认不校验)")
    ap.add_argument("--perf", dest="perf", action="store_true", default=None,
                    help="生成性能分析报告 (默认开启)")
    ap.add_argument("--no-perf", dest="perf", action="store_false",
                    help="跳过性能分析报告, 只输出位移 CSV")
    ap.add_argument("--md", dest="md", action="store_true", default=False,
                    help="性能报告输出为 markdown (.md, 默认, 无需显式指定)")
    ap.add_argument("--fluid-log", default=None,
                    help="fluid 求解器日志路径 (默认 glob <case>/precice-run/log.fluid* 推断)")
    ap.add_argument("--solid-log", default=None,
                    help="solid 求解器日志路径 (默认 glob <case>/precice-run/log.solid* 推断)")
    ap.add_argument("--no-reconstruct", action="store_true",
                    help="skip reconstructPar when processor results are pending")
    ap.add_argument("--self-test", action="store_true",
                    help="run an isolated fixed-topology extraction test and exit")
    args = ap.parse_args(argv)

    if args.self_test:
        return run_self_test()
    if not args.case:
        ap.error("--case is required unless --self-test is used")

    if args.perf is None:
        args.perf = True  # 默认 --perf 开

    print("[PRECICE-PIPE] ===== 步骤 5: 并行重建与挡板位移后处理 =====")

    case_root = os.path.abspath(args.case)
    case = os.path.join(case_root, "fluid-openfoam")
    if not os.path.isdir(case):
        sys.stderr.write("[PRECICE ERR] --case 缺少 fluid-openfoam: %s\n" % case_root)
        return 2
    if args.out_dir:
        out_dir = args.out_dir
    else:
        out_dir = os.path.join(case_root, "postprocess")
    os.makedirs(out_dir, exist_ok=True)

    if not args.no_reconstruct and not reconstruct_parallel_results(case):
        return 2

    # ---- 1. 常量网格 + 节点定位 (pointDisplacement 无内嵌坐标, 用常量网格定位) ----
    base_path = os.path.join(case, "constant", "polyMesh", "points")
    if not os.path.isfile(base_path):
        sys.stderr.write("[PRECICE ERR] 找不到常量网格 %s (需先用它定位节点)\n" % base_path)
        return 2
    n0, base = parse_vectors(base_path)
    print("[PRECICE 05] 常量网格: %s (%d 点)" % (base_path, len(base)))
    if args.res is not None and len(base) != args.res:
        sys.stderr.write("[PRECICE WARN] --res %d 与实测网格点数 %d 不符 (仅警告, 不中断)\n"
                         % (args.res, len(base)))

    idx = locate_node(base, args.node_x, args.node_y, args.node_z)
    if idx is None:
        sys.stderr.write("[PRECICE ERR] 常量网格中找不到挡板右上角 (%g, %g, %g) 附近的节点\n"
                         % (args.node_x, args.node_y, args.node_z))
        return 2
    print("[PRECICE 05] 挡板右上角节点 #%d: (%g, %g, %g)"
          % (idx, base[idx][0], base[idx][1], base[idx][2]))

    # ---- 2. 收集数值时步目录 (过滤 t=0 与 uniform-only 中间步) ----
    time_dirs = []
    for d in glob.glob(os.path.join(case, "[0-9]*")):
        try:
            t = float(os.path.basename(d))
        except ValueError:
            continue
        if t <= 0.0 or not os.path.isdir(d):
            continue
        time_dirs.append((t, d))
    time_dirs.sort()
    if not time_dirs:
        sys.stderr.write("[PRECICE ERR] 算例 %s 下找不到任何数值时步目录 (先跑步骤 3/4)\n" % case)
        return 2

    # ---- 3. 逐时步读取 pointDisplacement (累计位移, 按点序) ----
    rows = []
    skipped = []
    for t, d in time_dirs:
        pd_path = os.path.join(d, "pointDisplacement")
        if not os.path.isfile(pd_path):
            # uniform-only 中间步 (如 0.0095, 只有 uniform/, 无点位移场) 直接
            # 过滤; 顶层缺失但 processor0 有该场 -> 疑似未 reconstruct, warn 跳过
            skipped.append(t)
            proc_pd = os.path.join(case, "processor0", os.path.basename(d),
                                   "pointDisplacement")
            if os.path.isfile(proc_pd):
                sys.stderr.write(
                    "[PRECICE WARN] 顶层 %s 缺失但 processor0 存在 (疑似未 reconstruct), 跳过 t=%g\n"
                    % (pd_path, t))
            else:
                print("[PRECICE 05] t=%g 无 pointDisplacement (uniform-only 中间步), 过滤"
                      % t)
            continue
        n, disp = parse_vectors(pd_path)
        if n0 is not None and n is not None and n != n0:
            sys.stderr.write(
                "[PRECICE WARN] t=%g pointDisplacement 向量数 %d != 网格点数 %d, 跳过\n"
                % (t, n, n0))
            skipped.append(t)
            continue
        dx, dy, dz = disp[idx]
        rows.append((t, dx, dy, dz))
        print("[PRECICE 05] t=%-6g dx=%.8g dy=%.8g dz=%.8g" % (t, dx, dy, dz))

    # ---- 4. 写 CSV ----
    if not rows:
        sys.stderr.write("[PRECICE ERR] 没有任何有效时步, 未生成数据\n")
        return 2
    csv_path = os.path.join(out_dir, "precice_gate_topright_disp.csv")
    with open(csv_path, "w") as f:
        f.write("Time,dx,dy,dz\n")
        for t, dx, dy, dz in rows:
            f.write("%.8g,%.8g,%.8g,%.8g\n" % (t, dx, dy, dz))
    print("[PRECICE 05] 写入 %s (%d 行)" % (csv_path, len(rows)))
    if skipped:
        print("[PRECICE 05] 跳过 %d 个时步: %s"
              % (len(skipped), ", ".join("%g" % t for t in skipped)))

    # ---- 5. 性能分析报告 (默认开; 缺失日志仅 warn, 不影响返回码) ----
    if args.perf:
        build_precice_perf_report(case, out_dir, csv_path, rows,
                                  args.fluid_log, args.solid_log)
    else:
        print("[PRECICE 05] --no-perf: 跳过性能报告, 仅输出位移 CSV")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
