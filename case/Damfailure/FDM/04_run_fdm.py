#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""04_run_fdm.py

Run the FDM (fastDynamicFvMesh + myInterFoam) dam-break case.

This is a pure-stdlib Python rewrite of the `run_myInterFoam.sh` workflow:

  * three start modes:  fresh (-f) / restart (-r) / continue (-c)
  * serial or MPI-parallel execution (--nprocs N)
  * optional per-run controlDict overrides (--end-time / --write-interval)
  * pre-run validation (tools, mode files, node-count match)
  * built-in expected wall-clock estimate calibrated on measured runs
  * post-run verification ("End" marker, time directory, ExecutionTime)

Usage:
    python3 04_run_fdm.py -f [--nprocs 8] [--end-time 0.5] [--write-interval 50]
    python3 04_run_fdm.py -c --nprocs 8
    python3 04_run_fdm.py --case /path/to/case -f --dry-run

Modes:
    -f, --fresh      clean all results, run checkMesh (+decomposePar), then
                     solve from startTime. Step 05 reconstructs parallel output.
    -r, --restart    keep mesh/decomposition, clean results, solve from
                     startTime (skips checkMesh/decomposePar).
    -c, --continue   keep everything, solve from latestTime.

The solver is started as a blocking child (Popen + live tail of the log), so a
parent-process exit never orphans/kills the run.  Log files are written inside
the case directory:

    log.checkMesh  log.decomposePar  log.interFoam_serial | log.interFoam_parallel_N

Time-estimation model (calibrated on measured dam-failure runs):

    serial   : 51,792 cells, 500 steps (endTime 0.5, deltaT 0.001) -> 9,664 s
    parallel : 181,200 cells, 8 procs, 500 steps                 -> 5,510 s

    serial_cost_per_cell_per_step = 9664 / (51792 * 500)  ~= 3.73e-4 s
    parallel_total = serial_total / (nprocs * 0.77)        (0.77 parallel efficiency)
"""

import argparse
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import time

# --------------------------------------------------------------------------- #
#  Measured benchmark constants (time-estimation calibration)                  #
# --------------------------------------------------------------------------- #

BENCH_CELLS = 51792
BENCH_POINTS = 58175
BENCH_STEPS = 500                    # endTime 0.5 / deltaT 0.001
BENCH_SERIAL_SECONDS = 9664.0        # measured: small mesh, serial, 0.5 s
FULL_RES_CELLS = 181200
FULL_RES_POINTS = 200109
FULL_RES_NPROCS = 8
FULL_RES_SECONDS = 5510.0            # measured: full resolution, 8 procs, 0.5 s

# ~0.187 s/cell over 500 steps -> 3.73e-4 s/cell/step (serial)
SERIAL_COST_PER_CELL_PER_STEP = BENCH_SERIAL_SECONDS / (BENCH_CELLS * BENCH_STEPS)
PARALLEL_EFFICIENCY = 0.77           # calibrated from 8-proc full-resolution run

OPENFOAM_BASHRC_DEFAULT = "/usr/lib/openfoam/openfoam2412/etc/bashrc"

MODE_NAMES = {
    "f": "fresh (clean all, checkMesh + decomposePar, from startTime)",
    "r": "restart (keep mesh/decomposition, clean results, from startTime)",
    "c": "continue (keep all, from latestTime)",
}

# --------------------------------------------------------------------------- #
#  Logging                                                                     #
# --------------------------------------------------------------------------- #


def log(msg=""):
    print(("[FDM 04] " + msg) if msg else "[FDM 04]", flush=True)


def warn(msg):
    print("[FDM WARN] " + msg, flush=True)


def error(msg):
    print("[FDM ERR] " + msg, file=sys.stderr, flush=True)


def raw(msg=""):
    print(msg, flush=True)


# --------------------------------------------------------------------------- #
#  OpenFOAM environment helpers                                                #
# --------------------------------------------------------------------------- #


def locate_bashrc():
    exe = shutil.which("blockMesh")
    if exe:
        proj = os.path.dirname(os.path.dirname(os.path.dirname(exe)))
        rc = os.path.join(proj, "etc", "bashrc")
        if os.path.isfile(rc):
            return rc
    if os.path.isfile(OPENFOAM_BASHRC_DEFAULT):
        return OPENFOAM_BASHRC_DEFAULT
    return None


def resolve_foam_env():
    """Source the OpenFOAM bashrc once and capture the key paths/env vars."""
    rc = locate_bashrc()
    if not rc:
        return None
    keys = ["FOAM_APPBIN", "FOAM_LIBBIN", "FOAM_USER_APPBIN",
            "FOAM_USER_LIBBIN", "WM_PROJECT_DIR", "FOAM_MPI", "WM_MPLIB"]
    cmd = ('source "%s" >/dev/null 2>&1 && '
           'printf "%%s\\n" "$FOAM_APPBIN" "$FOAM_LIBBIN" "$FOAM_USER_APPBIN" '
           '"$FOAM_USER_LIBBIN" "$WM_PROJECT_DIR" "$FOAM_MPI" "$WM_MPLIB"') % rc
    try:
        proc = subprocess.run(["bash", "-c", cmd], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True, timeout=120)
    except Exception:
        return None
    lines = proc.stdout.splitlines()
    if len(lines) < len(keys):
        return None
    env = dict(zip(keys, lines))
    env["bashrc"] = rc
    return env


def foam_shell(cmd, bashrc):
    if bashrc:
        return 'source "%s" >/dev/null 2>&1 && %s' % (bashrc, cmd)
    return cmd


def run_capture(cmd, cwd, bashrc, log_file=None):
    """Run a command through the sourced bashrc, capturing stdout/stderr."""
    shell = foam_shell(cmd, bashrc)
    proc = subprocess.run(["bash", "-c", shell], cwd=cwd,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True)
    output = proc.stdout or ""
    if log_file:
        path = os.path.join(cwd, log_file)
        with open(path, "w", encoding="utf-8", errors="replace") as fh:
            fh.write(output)
    return proc.returncode, output


def run_stream(cmd, cwd, bashrc, log_file, end_time=None, delta_t=None,
               n_steps=None, label="solver"):
    """Run a long-running command, streaming progress to a log + stdout.

    The child is kept as a child of this process (we block and poll), so the
    run is never killed by a parent shell exit.

    Progress mode (end_time and n_steps both provided, e.g. the solver main
    loop): stdout shows a single self-overwriting progress line (``\\r``) with
    current Time / endTime, percentage, a bar graphic, elapsed wall time,
    remaining time and an estimated finish time (ETA).  The child's full stdout
    is still written to the log file, so the complete ``Time =`` /
    ``ExecutionTime =`` records stay available there.

    Non-progress mode (no progress parameters, e.g. check-type commands):
    the old behaviour is kept -- every new line starting with one of
    PROGRESS_PREFIXES is printed to stdout as it appears.
    """
    shell = foam_shell(cmd, bashrc)
    log_path = os.path.join(cwd, log_file)
    progress = (end_time is not None and n_steps is not None)
    with open(log_path, "w", encoding="utf-8", errors="replace") as fh:
        proc = subprocess.Popen(["bash", "-c", shell], cwd=cwd,
                                stdout=fh, stderr=subprocess.STDOUT)
        pos = 0
        t0 = time.time()
        latest = {}
        if progress:
            _render_progress(label, 0.0, end_time, 0.0, t0)
        while True:
            time.sleep(1.0)
            code = proc.poll()
            pos, matched, new_vals = _tail_live(log_path, pos)
            if new_vals:
                latest.update(new_vals)
            if progress:
                wall = time.time() - t0
                cur = latest.get("time", 0.0)
                if code is not None or latest.get("end"):
                    # final read to flush any trailing lines, then freeze the bar
                    time.sleep(0.2)
                    _pos, _matched, nv = _tail_live(log_path, pos)
                    if nv:
                        latest.update(nv)
                    wall = time.time() - t0
                    cur = latest.get("time", 0.0)
                    _render_progress(label, cur, end_time, wall, t0, final=True)
                    break
                _render_progress(label, cur, end_time, wall, t0)
            else:
                for line in matched:
                    print(line, flush=True)
                if code is not None:
                    break
        code = proc.wait()

    if progress:
        cur = latest.get("time", 0.0)
        wall = time.time() - t0
        if code == 0:
            log("=== 求解器完成: 进度到达 %s/%s (步 %d), 耗时 %s"
                % (_fmt_prog(cur), _fmt_prog(end_time), n_steps,
                   _fmt_wall(wall)))
        else:
            error("求解器异常退出 (exit %d), 进度停留 %s/%s; 详见 %s"
                  % (code, _fmt_prog(cur), _fmt_prog(end_time), log_file))
    return code


PROGRESS_PREFIXES = ("Time = ", "ExecutionTime = ", "Courant Number", "End")

BAR_WIDTH = 20


def _tail_live(log_path, pos):
    """Read appended lines since offset pos; return (new_pos, matched, latest)."""
    try:
        size = os.path.getsize(log_path)
    except OSError:
        return pos, [], {}
    if size < pos:
        pos = 0
    if size <= pos:
        return pos, [], {}
    with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
        fh.seek(pos)
        data = fh.read()
    new_pos = pos + len(data)
    matched = []
    latest = {}
    for line in data.splitlines():
        if line.startswith(PROGRESS_PREFIXES):
            matched.append(line)
        if line.startswith("Time = "):
            m = re.search(r'Time\s*=\s*([0-9.eE+-]+)', line)
            if m:
                try:
                    latest["time"] = float(m.group(1))
                except ValueError:
                    pass
        elif line.startswith("ExecutionTime = "):
            m = re.search(r'ExecutionTime\s*=\s*([0-9.eE+-]+)\s*s', line)
            if m:
                try:
                    latest["exec"] = float(m.group(1))
                except ValueError:
                    pass
        elif line.startswith("Courant Number"):
            latest["courant"] = True
        elif line.startswith("End"):
            latest["end"] = True
    return new_pos, matched, latest


def _fmt_prog(value):
    return "%.3f" % value


def _fmt_wall(seconds):
    """m:ss below an hour, h:mm:ss above; '--:--' when unknown."""
    if seconds is None:
        return "--:--"
    s = int(round(seconds))
    if s < 0:
        s = 0
    h = s // 3600
    m = (s % 3600) // 60
    sec = s % 60
    if h > 0:
        return "%d:%02d:%02d" % (h, m, sec)
    return "%d:%02d" % (m, sec)


def _fmt_eta(epoch):
    return time.strftime("%H:%M:%S", time.localtime(epoch))


def _render_progress(label, cur_time, end_time, wall_elapsed, t0, final=False):
    """Single-line bar; final=True ends with a newline, else \\r overwrite."""
    end = end_time or 0.0
    progress = min(1.0, cur_time / end) if end > 0 else 0.0
    filled = int(round(progress * BAR_WIDTH))
    bar = "█" * filled + "░" * (BAR_WIDTH - filled)
    pct = int(round(progress * 100))
    remaining = None
    eta = None
    if progress > 0:
        remaining = wall_elapsed / progress - wall_elapsed
        eta = t0 + wall_elapsed / progress
    line = ("[%s] %s %3d%% | %s/%s s | 已跑 %s | 剩余 %s | 预计完成 ~%s"
            % (label, bar, pct, _fmt_prog(cur_time), _fmt_prog(end),
               _fmt_wall(wall_elapsed), _fmt_wall(remaining),
               _fmt_eta(eta) if eta is not None else "--:--:--"))
    if final:
        print(line, flush=True)
    else:
        print(line, end="\r", flush=True)


# --------------------------------------------------------------------------- #
#  Dictionary parsing / editing                                                #
# --------------------------------------------------------------------------- #


def _dict_match(text, key):
    return re.search(r'^\s*%s\s+([^\s;]+)\s*;' % re.escape(key), text,
                     flags=re.M)


def _dict_value(text, key, default=""):
    m = _dict_match(text, key)
    return m.group(1) if m else default


def _dict_value_float(text, key, default):
    try:
        return float(_dict_value(text, key))
    except (TypeError, ValueError):
        return default


def _dict_value_int(text, key, default):
    try:
        return int(_dict_value(text, key))
    except (TypeError, ValueError):
        return default


def _fmt_num(value):
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value == int(value):
            return str(int(value))
        return format(value, ".10g")
    return str(value)


def read_control_dict(case_dir):
    path = os.path.join(case_dir, "system", "controlDict")
    text = _read_text(path)
    return {
        "application": _dict_value(text, "application", ""),
        "startFrom": _dict_value(text, "startFrom", "startTime"),
        "startTime": _dict_value_float(text, "startTime", 0.0),
        "endTime": _dict_value_float(text, "endTime", 0.5),
        "deltaT": _dict_value_float(text, "deltaT", 0.001),
        "writeControl": _dict_value(text, "writeControl", "timeStep"),
        "writeInterval": _dict_value_int(text, "writeInterval", 1),
    }


def _read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def set_control_dict_value(case_dir, key, value):
    path = os.path.join(case_dir, "system", "controlDict")
    text = _read_text(path)
    new_text, n = re.subn(r'(^\s*%s\s+)[^\s;]+(\s*;)' % re.escape(key),
                          r'\g<1>%s\2' % _fmt_num(value), text, count=1,
                          flags=re.M)
    if n == 0:
        raise RuntimeError("key '%s' not found in system/controlDict" % key)
    with open(path, "w", encoding="utf-8", errors="replace") as fh:
        fh.write(new_text)


def set_start_from(case_dir, target):
    if target not in ("startTime", "latestTime"):
        raise RuntimeError("unsupported startFrom value: %s" % target)
    set_control_dict_value(case_dir, "startFrom", target)


# --------------------------------------------------------------------------- #
#  decomposeParDict handling                                                   #
# --------------------------------------------------------------------------- #

DECOMPOSE_DICT_TEMPLATE = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

numberOfSubdomains %d;

method          scotch;

distributed     no;

roots           ();
"""


def ensure_decompose_dict(case_dir, nprocs):
    """Create a default scotch decomposeParDict, or update numberOfSubdomains."""
    path = os.path.join(case_dir, "system", "decomposeParDict")
    if os.path.isfile(path):
        text = _read_text(path)
        new_text, n = re.subn(r'(numberOfSubdomains\s+)\d+(\s*;)',
                              r'\g<1>%d\2' % nprocs, text, count=1)
        if n == 0:
            raise RuntimeError(
                "numberOfSubdomains not found in system/decomposeParDict")
        with open(path, "w", encoding="utf-8", errors="replace") as fh:
            fh.write(new_text)
        return "update"
    with open(path, "w", encoding="utf-8", errors="replace") as fh:
        fh.write(DECOMPOSE_DICT_TEMPLATE % nprocs)
    return "create"


# --------------------------------------------------------------------------- #
#  Mesh / mode-file counting                                                   #
# --------------------------------------------------------------------------- #

TIME_DIR_RE = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$')


def is_time_dir(name):
    return bool(TIME_DIR_RE.match(name))


def read_poly_mesh_counts(case_dir):
    """Return (n_cells, n_points) from constant/polyMesh headers.

    Primary source is the `note` line OpenFOAM writes into the polyMesh files
    (nPoints/nCells).  Fall back to the point-count line in `points`, then to a
    previously written log.checkMesh.
    """
    poly = os.path.join(case_dir, "constant", "polyMesh")
    n_cells = n_points = None

    owner = os.path.join(poly, "owner")
    if os.path.isfile(owner):
        text = _read_text(owner)
        m = re.search(r'nCells\s*:\s*(\d+)', text)
        if m:
            n_cells = int(m.group(1))
        m = re.search(r'nPoints\s*:\s*(\d+)', text)
        if m:
            n_points = int(m.group(1))

    if n_points is None:
        pts = os.path.join(poly, "points")
        if os.path.isfile(pts):
            text = _read_text(pts)
            m = re.search(r'\n\s*(\d+)\s*\n\s*\(', text)
            if m:
                n_points = int(m.group(1))

    if n_cells is None or n_points is None:
        ccells, cpoints = read_checkmesh_counts(case_dir)
        if n_cells is None:
            n_cells = ccells
        if n_points is None:
            n_points = cpoints

    return n_cells, n_points


def read_checkmesh_counts(case_dir):
    path = os.path.join(case_dir, "log.checkMesh")
    if not os.path.isfile(path):
        return None, None
    text = _read_text(path)
    mc = re.search(r'cells:\s*(\d+)', text)
    mp = re.search(r'points:\s*(\d+)', text)
    return ((int(mc.group(1)) if mc else None),
            (int(mp.group(1)) if mp else None))


def read_mode_node_count(case_dir):
    """Return the node count declared in mode/FluidNodeCoor.csv (2nd field).

    Header format is `dummy, nNode, nMode`.
    """
    path = os.path.join(case_dir, "mode", "FluidNodeCoor.csv")
    if not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
    nums = []
    for part in re.split(r'[,\s]+', first.strip()):
        try:
            nums.append(float(part))
        except ValueError:
            pass
    if len(nums) >= 2:
        return int(nums[1])
    return None


# --------------------------------------------------------------------------- #
#  Time estimation                                                             #
# --------------------------------------------------------------------------- #


def estimate_time(n_cells, nprocs, end_time, delta_t):
    if delta_t is None or delta_t <= 0:
        delta_t = 0.001
    if n_cells is None:
        n_cells = BENCH_CELLS
    n_steps = max(1, int(round(end_time / delta_t)))
    serial_total = n_cells * n_steps * SERIAL_COST_PER_CELL_PER_STEP
    if nprocs <= 1:
        seconds = serial_total
    else:
        seconds = serial_total / (nprocs * PARALLEL_EFFICIENCY)
    return {"n_cells": n_cells, "n_steps": n_steps, "nprocs": nprocs,
            "seconds": seconds, "delta_t": delta_t}


def format_duration(seconds):
    s = int(round(seconds))
    h = s // 3600
    m = (s % 3600) // 60
    sec = s % 60
    if h > 0:
        return "%dh%02dm" % (h, m)
    if m > 0:
        return "%dm%02ds" % (m, sec)
    return "%ds" % sec


def print_estimate(est, n_points, end_time, delta_t):
    log("预期计算时间估算:")
    log("  网格单元: %s" % format(est["n_cells"], ","))
    if n_points is not None:
        log("  网格点数: %s" % format(n_points, ","))
    log("  时间步数: %d (endTime %s, deltaT %s)"
        % (est["n_steps"], _fmt_num(end_time), _fmt_num(delta_t)))
    log("  进程数: %d" % est["nprocs"])
    log("  预计耗时: ~%s s (~%s)"
        % (format(int(round(est["seconds"])), ","), format_duration(est["seconds"])))
    log("  说明: 基于实测标定 (51,792 单元串行 0.5s = 9,664s; "
        "181,200 单元 8 核 = 5,510s)")


# --------------------------------------------------------------------------- #
#  Cleaning                                                                    #
# --------------------------------------------------------------------------- #


def _numeric_dirs(case_dir):
    out = []
    try:
        entries = os.listdir(case_dir)
    except OSError:
        return out
    for name in entries:
        if is_time_dir(name):
            out.append(name)
    return out


def collect_root_cleanup(case_dir):
    """List (without removing) root result dirs / logs to delete.

    Result time directories are removed except the initial '0' (float == 0.0).
    This is a pure read-only pass so --dry-run cannot mutate the case.
    """
    targets = []
    for name in _numeric_dirs(case_dir):
        try:
            if float(name) == 0.0:
                continue
        except ValueError:
            pass
        path = os.path.join(case_dir, name)
        if os.path.isdir(path):
            targets.append(path)
    for pattern in ("log.*", "*.log"):
        for path in glob.glob(os.path.join(case_dir, pattern)):
            targets.append(path)
    modal = os.path.join(case_dir, "modal_diagnostics.csv")
    if os.path.exists(modal):
        targets.append(modal)
    return targets


def collect_processor_cleanup(case_dir, keep_mesh):
    """List processor* paths to delete.

    keep_mesh=False -> whole processor directories.
    keep_mesh=True  -> only per-processor result time dirs (keep constant/ and
                       initial 0/).
    """
    targets = []
    for pdir in glob.glob(os.path.join(case_dir, "processor*")):
        if not os.path.isdir(pdir):
            continue
        if not keep_mesh:
            targets.append(pdir)
            continue
        for entry in os.listdir(pdir):
            if not is_time_dir(entry):
                continue
            try:
                if float(entry) == 0.0:
                    continue
            except ValueError:
                pass
            path = os.path.join(pdir, entry)
            if os.path.isdir(path):
                targets.append(path)
    return targets


def remove_paths(paths):
    """Actually delete the given paths; return their base names."""
    removed = []
    for path in paths:
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)
        else:
            continue
        removed.append(os.path.basename(path))
    return removed


def has_processors(case_dir):
    return any(os.path.isdir(p) for p in glob.glob(os.path.join(case_dir, "processor*")))


# --------------------------------------------------------------------------- #
#  Tool paths / validation                                                     #
# --------------------------------------------------------------------------- #


def solver_path(foam_env):
    return os.path.join(foam_env["FOAM_USER_APPBIN"], "myInterFoam")


def lib_path(foam_env):
    return os.path.join(foam_env["FOAM_USER_LIBBIN"], "libfastDynamicFvMesh.so")


def validate_tools(foam_env, mode, nprocs):
    appbin = foam_env["FOAM_APPBIN"]
    checks = [
        ("myInterFoam", solver_path(foam_env), "x"),
        ("libfastDynamicFvMesh.so", lib_path(foam_env), "f"),
    ]
    if mode == "f":
        checks.append(("checkMesh", os.path.join(appbin, "checkMesh"), "x"))
    if nprocs > 1:
        checks.append(("decomposePar", os.path.join(appbin, "decomposePar"), "x"))
        if not shutil.which("mpirun"):
            error("mpirun not found on PATH (required for parallel runs)")
            return False

    ok = True
    for name, path, kind in checks:
        if not os.path.isfile(path):
            error("%s not found: %s" % (name, path))
            ok = False
        elif kind == "x" and not os.access(path, os.X_OK):
            error("%s not executable: %s" % (name, path))
            ok = False
        else:
            log("OK: %s -> %s" % (name, path))
    return ok


def build_solver_command(foam_env, nprocs, taskset, root):
    solver = solver_path(foam_env)
    if nprocs <= 1:
        cmd = solver
    else:
        mpirun = shutil.which("mpirun") or "mpirun"
        parts = [mpirun]
        if root:
            parts.append("--allow-run-as-root")
        parts.append("-np %d" % nprocs)
        parts.append(solver)
        parts.append("-parallel")
        cmd = " ".join(parts)
    if taskset:
        cmd = "taskset -c %s %s" % (taskset, cmd)
    return cmd


# --------------------------------------------------------------------------- #
#  Post-run verification                                                       #
# --------------------------------------------------------------------------- #


def parse_execution_time(log_path):
    if not os.path.isfile(log_path):
        return None, None
    exec_time = None
    clock_time = None
    with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.search(r'ExecutionTime\s*=\s*([0-9.eE+-]+)\s*s', line)
            if m:
                try:
                    exec_time = float(m.group(1))
                except ValueError:
                    pass
            m = re.search(r'ClockTime\s*=\s*([0-9.eE+-]+)\s*s', line)
            if m:
                try:
                    clock_time = float(m.group(1))
                except ValueError:
                    pass
    return exec_time, clock_time


def find_time_dir(case_dir, value, tol=1e-9):
    for name in _numeric_dirs(case_dir):
        try:
            if abs(float(name) - value) <= tol * max(1.0, abs(value)):
                return name
        except ValueError:
            pass
    return None


def verify_results(case_dir, end_time, log_file, nprocs):
    log_path = os.path.join(case_dir, log_file)
    problems = []

    if os.path.isfile(log_path):
        text = _read_text(log_path)
        if re.search(r'^\s*End\s*$', text, flags=re.M):
            log("OK: solver log '%s' contains 'End' marker" % log_file)
        else:
            problems.append("solver log '%s' missing 'End' marker" % log_file)
    else:
        problems.append("solver log '%s' not found" % log_file)

    end_dir = find_time_dir(case_dir, end_time)
    if end_dir is not None:
        log("OK: end-time directory '%s' exists" % end_dir)
    elif nprocs > 1:
        processor_end = find_time_dir(os.path.join(case_dir, "processor0"), end_time)
        if processor_end is not None:
            log("OK: processor0 end-time directory '%s' exists; "
                "step 05 will reconstruct parallel results" % processor_end)
        else:
            problems.append("no processor0 time directory found for endTime %s"
                            % _fmt_num(end_time))
    else:
        problems.append("no time directory found for endTime %s"
                        % _fmt_num(end_time))

    exec_time, clock_time = parse_execution_time(log_path)
    if exec_time is not None:
        log("ExecutionTime = %s s%s" % (format(exec_time, ".2f"),
            (" (ClockTime = %s s)" % format(clock_time, ".2f")
             if clock_time is not None else "")))
    else:
        problems.append("could not parse ExecutionTime from '%s'" % log_file)

    return problems


# --------------------------------------------------------------------------- #
#  Argument parsing                                                            #
# --------------------------------------------------------------------------- #


def parse_args(argv):
    p = argparse.ArgumentParser(
        prog="04_run_fdm.py",
        description="Run the FDM (fastDynamicFvMesh + myInterFoam) dam-break "
                    "case with fresh/restart/continue modes, optional parallel "
                    "execution, and a built-in expected-time estimate.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python3 04_run_fdm.py -f --nprocs 8\n"
               "  python3 04_run_fdm.py -c --nprocs 8\n"
               "  python3 04_run_fdm.py -f --end-time 0.01 --write-interval 5 --dry-run\n"
    )
    g = p.add_mutually_exclusive_group()
    g.add_argument("-f", "--fresh", dest="mode", action="store_const",
                   const="f", help="fresh start (default)")
    g.add_argument("-r", "--restart", dest="mode", action="store_const",
                   const="r", help="restart: keep mesh/decomposition, clean results")
    g.add_argument("-c", "--continue", dest="mode", action="store_const",
                   const="c", help="continue from latestTime")
    p.set_defaults(mode="f")

    p.add_argument("--case", default=".",
                   help="case directory (default: current directory)")
    p.add_argument("--nprocs", type=int, default=8,
                   help="number of MPI processes (default: 8)")
    p.add_argument("--end-time", type=float, default=None,
                   help="override endTime in controlDict (e.g. 0.5)")
    p.add_argument("--write-interval", type=int, default=None,
                   help="override writeInterval in controlDict")
    p.add_argument("--taskset", default=None,
                   help="optional CPU affinity mask for the solver "
                        "(e.g. '0-7')")
    p.add_argument("--dry-run", action="store_true",
                   help="print plan + time estimate only; do not execute")
    return p.parse_args(argv)


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #


def main(argv):
    args = parse_args(argv)
    case_dir = os.path.abspath(args.case)
    dry = args.dry_run
    mode = args.mode
    nprocs = args.nprocs

    if nprocs < 1:
        error("--nprocs must be a positive integer, got %d" % nprocs)
        return 2

    if not os.path.isfile(os.path.join(case_dir, "system", "controlDict")):
        error("no system/controlDict in case directory: %s" % case_dir)
        return 2

    foam_env = resolve_foam_env()
    if not foam_env:
        error("could not locate OpenFOAM bashrc; cannot run OpenFOAM tools")
        return 2

    root = (os.geteuid() == 0)
    if root:
        os.environ["OMPI_ALLOW_RUN_AS_ROOT"] = "1"
        os.environ["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] = "1"

    # ---- read case configuration -----------------------------------------
    ctrl = read_control_dict(case_dir)
    end_time = args.end_time if args.end_time is not None else ctrl["endTime"]
    delta_t = ctrl["deltaT"]
    write_interval = (args.write_interval if args.write_interval is not None
                      else ctrl["writeInterval"])

    n_cells, n_points = read_poly_mesh_counts(case_dir)
    mode_nodes = read_mode_node_count(case_dir)

    # ---- header -----------------------------------------------------------
    raw("[FDM-PIPE] ===== 步骤 4/6: FDM 运行 =====")
    log("FDM run (fastDynamicFvMesh + myInterFoam)")
    log("--- OpenFOAM bashrc: %s" % foam_env["bashrc"])
    log("--- Case directory: %s" % case_dir)
    log("--- Mode: %s (%s)" % (mode, MODE_NAMES[mode]))
    log("--- Processors: %d (%s)"
        % (nprocs, "parallel" if nprocs > 1 else "serial"))
    if root:
        log("--- Running as root: enabling OpenMPI root override")

    # ---- pre-run validation ----------------------------------------------
    if not validate_tools(foam_env, mode, nprocs):
        return 2

    if mode_nodes is None:
        error("mode/FluidNodeCoor.csv not found in case; required by "
              "fastDynamicFvMesh")
        return 2
    if n_points is not None and mode_nodes != n_points:
        error("mode node count (%d) != polyMesh point count (%d)"
              % (mode_nodes, n_points))
        return 2
    log("OK: mode/FluidNodeCoor.csv node count matches polyMesh points (%d)"
        % mode_nodes)

    # ---- time estimate -----------------------------------------------------
    est = estimate_time(n_cells, nprocs, end_time, delta_t)
    print_estimate(est, n_points, end_time, delta_t)

    # ---- command / plan assembly ------------------------------------------
    solver_cmd = build_solver_command(foam_env, nprocs, args.taskset, root)
    solver_log = ("log.interFoam_parallel_%d" % nprocs
                  if nprocs > 1 else "log.interFoam_serial")
    bashrc = foam_env["bashrc"]

    def show(cmd):
        log("      [dry-run] %s" % cmd)

    # ---- mode-specific preparation ----------------------------------------
    # --dry-run is strictly non-mutating: every file change below is gated on
    # `dry`, and cleanup is split into a read-only collect pass + remove pass.
    if mode == "f":
        log("--- Fresh mode: cleaning all results + processor dirs")
        targets = collect_root_cleanup(case_dir) + \
            collect_processor_cleanup(case_dir, keep_mesh=False)
        if dry:
            log("      [dry-run] would remove: %s"
                % (", ".join(os.path.basename(t) for t in targets)
                   if targets else "(nothing)"))
        else:
            remove_paths(targets)
            log("--- Removed: %s"
                % (", ".join(os.path.basename(t) for t in targets)
                   if targets else "(nothing)"))
        if not dry:
            set_start_from(case_dir, "startTime")
        log("--- Set startFrom = startTime in system/controlDict%s"
            % ("" if not dry else " [dry-run]"))

        if nprocs > 1:
            dict_path = os.path.join(case_dir, "system", "decomposeParDict")
            action = "update" if os.path.isfile(dict_path) else "create"
            if dry:
                log("      [dry-run] decomposeParDict: %s numberOfSubdomains=%d"
                    % (action, nprocs))
            else:
                ensure_decompose_dict(case_dir, nprocs)
                log("--- decomposeParDict: %s (numberOfSubdomains=%d)"
                    % (action, nprocs))

        # checkMesh (fresh only)
        log("--- Running checkMesh -> log.checkMesh")
        check_cmd = os.path.join(foam_env["FOAM_APPBIN"], "checkMesh")
        if dry:
            show(check_cmd)
        else:
            code, out = run_capture(check_cmd, case_dir, bashrc,
                                    log_file="log.checkMesh")
            mesh_ok = code == 0 and "Mesh OK" in out and "Failed" not in out
            if not mesh_ok:
                error("checkMesh failed; see log.checkMesh")
                return 2
            log("--- checkMesh OK")

        # setFields (fresh only; initializes water column from setFieldsDict
        # generated by 03_assemble_fdm_case.py).  Skipped if the dictionary is
        # missing: only a warning, do not crash.
        sf_dict = os.path.join(case_dir, "system", "setFieldsDict")
        if not os.path.isfile(sf_dict):
            warn("system/setFieldsDict not found; "
                 "skipping setFields (water column will NOT be initialized)")
        else:
            log("--- Running setFields -> log.setFields")
            setfields_cmd = os.path.join(foam_env["FOAM_APPBIN"], "setFields")
            if dry:
                show(setfields_cmd)
            else:
                code, out = run_capture(setfields_cmd, case_dir, bashrc,
                                        log_file="log.setFields")
                if code != 0:
                    error("setFields failed; see log.setFields")
                    return 2
                log("--- setFields OK")

        # decomposePar (parallel only)
        if nprocs > 1:
            log("--- Running decomposePar -force -> log.decomposePar")
            decomp_cmd = os.path.join(foam_env["FOAM_APPBIN"],
                                      "decomposePar") + " -force"
            if dry:
                show(decomp_cmd)
            else:
                code, out = run_capture(decomp_cmd, case_dir, bashrc,
                                        log_file="log.decomposePar")
                if code != 0:
                    error("decomposePar failed; see log.decomposePar")
                    return 2
                log("--- decomposePar OK")

    elif mode == "r":
        if nprocs > 1 and not has_processors(case_dir):
            error("restart mode needs processor* dirs; run -f first")
            return 2
        log("--- Restart mode: cleaning results, keeping mesh/decomposition")
        targets = collect_root_cleanup(case_dir) + \
            collect_processor_cleanup(case_dir, keep_mesh=True)
        if dry:
            log("      [dry-run] would remove: %s"
                % (", ".join(os.path.basename(t) for t in targets)
                   if targets else "(nothing)"))
        else:
            remove_paths(targets)
            log("--- Removed: %s"
                % (", ".join(os.path.basename(t) for t in targets)
                   if targets else "(nothing)"))
        if not dry:
            set_start_from(case_dir, "startTime")
        log("--- Set startFrom = startTime in system/controlDict%s"
            % ("" if not dry else " [dry-run]"))

    else:  # "c"
        if nprocs > 1 and not has_processors(case_dir):
            error("continue mode needs processor* dirs; run -f first")
            return 2
        log("--- Continue mode: keeping all results")
        if not dry:
            set_start_from(case_dir, "latestTime")
        log("--- Set startFrom = latestTime in system/controlDict%s"
            % ("" if not dry else " [dry-run]"))

    # apply optional overrides (restored after the run)
    if args.end_time is not None:
        log("--- Override endTime = %s (restored after run)" % _fmt_num(end_time))
        if not dry:
            set_control_dict_value(case_dir, "endTime", end_time)
    if args.write_interval is not None:
        log("--- Override writeControl=timeStep, writeInterval=%d (restored after run)"
            % write_interval)
        if not dry:
            set_control_dict_value(case_dir, "writeControl", "timeStep")
            set_control_dict_value(case_dir, "writeInterval", write_interval)

    # ---- dry run: stop here ----------------------------------------------
    if dry:
        log("--- Dry run: solver would run:")
        show(solver_cmd)
        log("--- Solver log: %s" % solver_log)
        log("Dry run complete (no commands executed).")
        return 0

    # ---- solve -------------------------------------------------------------
    log("--- Starting myInterFoam (%s) -> %s"
        % ("parallel" if nprocs > 1 else "serial", solver_log))
    log("--- 实时进度条 (Time/ExecutionTime/ETA) 如下")
    code = run_stream(solver_cmd, case_dir, bashrc, log_file=solver_log,
                      end_time=end_time, delta_t=delta_t,
                      n_steps=est["n_steps"], label="FDM 04")
    if code != 0:
        error("myInterFoam failed (exit %d); see %s" % (code, solver_log))
        _restore_overrides(case_dir, args, ctrl)
        return code

    # ---- restore optional overrides ---------------------------------------
    _restore_overrides(case_dir, args, ctrl)

    # ---- post-run verification --------------------------------------------
    log("--- Verifying results")
    problems = verify_results(case_dir, end_time, solver_log, nprocs)
    if problems:
        for pr in problems:
            error(pr)
        log("Run finished but verification reported issues (see above).")
        return 3
    exec_t, _clock_t = parse_execution_time(os.path.join(case_dir, solver_log))
    log("endTime: %s  steps: %d  ExecutionTime: %s  (进程: %d)"
        % (_fmt_num(end_time), est["n_steps"],
           ("%ss" % format(exec_t, ".2f")) if exec_t is not None else "?",
           nprocs))
    log("Run complete. Logs: log.checkMesh, log.decomposePar, %s" % solver_log)
    if nprocs > 1:
        log("Parallel results remain under processor*/; run step 05 to reconstruct")
    return 0


def _restore_overrides(case_dir, args, ctrl):
    try:
        if args.end_time is not None:
            set_control_dict_value(case_dir, "endTime", ctrl["endTime"])
        if args.write_interval is not None:
            set_control_dict_value(case_dir, "writeControl",
                                   ctrl.get("writeControl", "timeStep"))
            set_control_dict_value(case_dir, "writeInterval", ctrl["writeInterval"])
    except Exception as exc:
        error("failed to restore controlDict overrides: %s" % exc)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
