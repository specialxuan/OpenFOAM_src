#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""04_run_precice.py

Run the preCICE partitioned dam-break FSI case assembled by
03_assemble_precice_case.py.

The verified baseline zip runs the two solvers as two background subshells
(run-coupled.sh): solid `ccx_preCICE` and fluid `interFoam` are started
simultaneously, both stream into precice-run/, then a single `wait` collects
both exit codes.  This script reproduces that two-child-pair + wait behaviour
in Python (reusing the FDM 04 run_stream progress-bar pattern) so the parent
never orphans the solvers and a single self-overwriting progress bar shows the
fluid Time / endTime / percentage / 已跑 / 剩余 / 预计完成.

  * fluid  : `interFoam`  (original OpenFOAM solver, NOT myInterFoam)
             serial or MPI-parallel (--nprocs N)
  * solid  : `ccx_preCICE -i dam_gate -precice-participant Solid`
             always OMP_NUM_THREADS=4 / CCX_NPROC_EQUATION_SOLVER=4
             (--solid-threads to change)

Modes (mirror FDM 04):
  -f, --fresh      clean all results, run checkMesh (+setFields, +decomposePar
                   for parallel), solve from startTime, then reconstruct.
  -r, --restart    keep mesh/decomposition, clean results, solve from startTime.
  -c, --continue   keep everything, solve from latestTime (fluid only; the
                   preCICE coupling itself always starts from time-window 0).

Post-run verification:
  * fluid log contains the preCICE end marker ("Reached end at") and "End"
  * solid log contains "Job finished"
  * no "MAX ITERATIONS REACHED" / max-iterations error in either log
  * both exit codes are 0

Usage:
    python3 04_run_precice.py --case /path/to/precice_case -f [--nprocs 8] \\
        [--end-time 0.02] [--write-interval 0.01] [--dry-run]

Short run:    --end-time 0.02, deltaT 0.0005 -> 40 coupling time windows.
Full run:     --end-time 0.5 -> 1000 windows (~11 h fluid serial at full res).

Note on the short run: the OpenFOAM preciceAdapter overrides the fluid
controlDict endTime to infinity ("Only preCICE will control the simulation's
endTime"), so --end-time also patches <max-time> in precice-config.xml (the
true coupling stop).  The original <max-time> tag is restored after the run.
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
#  Logging / progress bar                                                      #
# --------------------------------------------------------------------------- #

def log(msg=""):
    print(("[PRECICE 04] " + msg) if msg else "[PRECICE 04]", flush=True)


def warn(msg):
    print("[PRECICE WARN] " + msg, flush=True)


def error(msg):
    print("[PRECICE ERR] " + msg, file=sys.stderr, flush=True)


def raw(msg=""):
    print(msg, flush=True)


PROGRESS_PREFIXES = ("Time = ", "ExecutionTime = ", "Courant Number", "End")
BAR_WIDTH = 20


def _fmt_prog(value):
    return "%.3f" % value


def _fmt_wall(seconds):
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


# --------------------------------------------------------------------------- #
#  OpenFOAM environment helpers                                                #
# --------------------------------------------------------------------------- #

OPENFOAM_BASHRC_DEFAULT = "/usr/lib/openfoam/openfoam2412/etc/bashrc"


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


# --------------------------------------------------------------------------- #
#  Dictionary parsing / editing                                                #
# --------------------------------------------------------------------------- #

def _read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read()


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
    path = os.path.join(case_dir, "fluid-openfoam", "system", "controlDict")
    text = _read_text(path)
    return {
        "application": _dict_value(text, "application", ""),
        "startFrom": _dict_value(text, "startFrom", "startTime"),
        "startTime": _dict_value_float(text, "startTime", 0.0),
        "endTime": _dict_value_float(text, "endTime", 0.02),
        "deltaT": _dict_value_float(text, "deltaT", 0.0005),
        "writeControl": _dict_value(text, "writeControl", "runTime"),
        "writeInterval": _dict_value_int(text, "writeInterval", 1),
    }


def set_control_dict_value(case_dir, key, value):
    path = os.path.join(case_dir, "fluid-openfoam", "system", "controlDict")
    text = _read_text(path)
    new_text, n = re.subn(r'(^\s*%s\s+)[^\s;]+(\s*;)' % re.escape(key),
                          r'\g<1>%s\2' % _fmt_num(value), text, count=1,
                          flags=re.M)
    if n == 0:
        raise RuntimeError("key '%s' not found in fluid-openfoam/system/controlDict"
                           % key)
    with open(path, "w", encoding="utf-8", errors="replace") as fh:
        fh.write(new_text)


def set_start_from(case_dir, target):
    if target not in ("startTime", "latestTime"):
        raise RuntimeError("unsupported startFrom value: %s" % target)
    set_control_dict_value(case_dir, "startFrom", target)


def read_precice_config(config_path):
    """Pull time-window-size / max-time from precice-config.xml."""
    if not os.path.isfile(config_path):
        return {}
    text = _read_text(config_path)
    out = {}
    m = re.search(r'<time-window-size\s+value="([^"]+)"', text)
    if m:
        try:
            out["window"] = float(m.group(1))
        except ValueError:
            pass
    m = re.search(r'<max-time\s+value="([^"]+)"', text)
    if m:
        try:
            out["max_time"] = float(m.group(1))
        except ValueError:
            pass
    m = re.search(r'<max-iterations\s+value="([^"]+)"', text)
    if m:
        try:
            out["max_iter"] = int(m.group(1))
        except ValueError:
            pass
    return out


PRECICE_MAX_TIME_RE = re.compile(r'(<max-time\s+value=")[^"]+(")')


def patch_precice_max_time(config_path, value):
    """Rewrite <max-time value="..."/> in precice-config.xml; return the
    original value string (for restoration) or None if the tag was absent."""
    text = _read_text(config_path)
    m = PRECICE_MAX_TIME_RE.search(text)
    if not m:
        return None
    original = m.group(0)
    new_text, n = PRECICE_MAX_TIME_RE.subn(
        r'\g<1>%s\2' % _fmt_num(value), text, count=1)
    if n != 1:
        return None
    with open(config_path, "w", encoding="utf-8") as fh:
        fh.write(new_text)
    return original


def restore_precice_max_time(config_path, original):
    """Put the original <max-time .../> tag back into precice-config.xml."""
    if original is None:
        return
    text = _read_text(config_path)
    new_text, n = PRECICE_MAX_TIME_RE.subn(lambda _m: original, text, count=1)
    if n != 1:
        warn("could not restore <max-time> in %s (already modified?)"
             % config_path)
        return
    with open(config_path, "w", encoding="utf-8") as fh:
        fh.write(new_text)


# --------------------------------------------------------------------------- #
#  Mesh / time-step helpers                                                    #
# --------------------------------------------------------------------------- #

TIME_DIR_RE = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$')


def is_time_dir(name):
    return bool(TIME_DIR_RE.match(name))


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


def read_poly_mesh_cells(fluid_dir):
    owner = os.path.join(fluid_dir, "constant", "polyMesh", "owner")
    if os.path.isfile(owner):
        m = re.search(r'nCells\s*:\s*(\d+)', _read_text(owner))
        if m:
            return int(m.group(1))
    return None


def collect_root_cleanup(fluid_dir):
    """List (without removing) fluid result time dirs (except initial 0)."""
    targets = []
    for name in _numeric_dirs(fluid_dir):
        try:
            if float(name) == 0.0:
                continue
        except ValueError:
            pass
        path = os.path.join(fluid_dir, name)
        if os.path.isdir(path):
            targets.append(path)
    return targets


def collect_processor_cleanup(fluid_dir, keep_mesh):
    targets = []
    for pdir in glob.glob(os.path.join(fluid_dir, "processor*")):
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


def has_processors(fluid_dir):
    return any(os.path.isdir(p)
               for p in glob.glob(os.path.join(fluid_dir, "processor*")))


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


def ensure_decompose_dict(fluid_dir, nprocs):
    path = os.path.join(fluid_dir, "system", "decomposeParDict")
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
#  Time estimation (calibrated on the verified baseline zip run)               #
# --------------------------------------------------------------------------- #

BENCH_CELLS = 181200
BENCH_STEPS = 1000                  # 0.5 / 0.0005
BENCH_SERIAL_SECONDS = 39185.0      # measured: full-res 181200 cells serial 0.5s
SERIAL_COST_PER_CELL_PER_STEP = BENCH_SERIAL_SECONDS / (BENCH_CELLS * BENCH_STEPS)
PARALLEL_EFFICIENCY = 0.77


def estimate_time(n_cells, nprocs, end_time, delta_t):
    if delta_t is None or delta_t <= 0:
        delta_t = 0.0005
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


# --------------------------------------------------------------------------- #
#  Post-run verification                                                       #
# --------------------------------------------------------------------------- #

def verify_run(case_dir, end_time, delta_t):
    """Check end markers + coupling convergence in precice-run/ logs."""
    run_dir = os.path.join(case_dir, "precice-run")
    fluid_log = os.path.join(run_dir, "log.fluid")
    solid_log = os.path.join(run_dir, "log.solid")
    problems = []

    fluid_text = _read_text(fluid_log) if os.path.isfile(fluid_log) else ""
    solid_text = _read_text(solid_log) if os.path.isfile(solid_log) else ""

    if not os.path.isfile(fluid_log):
        problems.append("precice-run/log.fluid not found")
    elif "Reached end at" in fluid_text:
        m = re.search(r'final time-window:\s*(\d+)', fluid_text)
        nwin = int(m.group(1)) if m else None
        log("OK: fluid reached coupling end (time-windows: %s)"
            % (nwin if nwin is not None else "?"))
        if nwin is not None and end_time > 0 and delta_t > 0:
            want = int(round(end_time / delta_t))
            if nwin < want:
                problems.append("fluid only reached %d coupling windows (expected %d)"
                                % (nwin, want))
    elif "End" in fluid_text:
        log("OK: fluid log contains 'End' marker (preCICE end line not found)")
    else:
        problems.append("fluid log missing both 'Reached end at' and 'End' markers")

    if not os.path.isfile(solid_log):
        problems.append("precice-run/log.solid not found")
    elif "Job finished" in solid_text:
        log("OK: solid log contains 'Job finished'")
    else:
        problems.append("solid log missing 'Job finished' marker")

    both = fluid_text + solid_text
    miter = re.search(r'MAX ITERATIONS REACHED|max-iterations reached',
                      both, flags=re.I)
    if miter:
        problems.append("preCICE reported max coupling iterations reached: %s"
                        % miter.group(0))
    else:
        log("OK: no max-iterations / divergence marker in preCICE logs")

    return problems


def parse_execution_time(log_path):
    if not os.path.isfile(log_path):
        return None, None
    exec_time = clock_time = None
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


# --------------------------------------------------------------------------- #
#  Solver launch (two background children + wait, mirroring run-coupled.sh)     #
# --------------------------------------------------------------------------- #

def build_solver_commands(foam_env, nprocs, taskset, root):
    """Return (fluid_cmd, solid_cmd) strings run inside their subdirectories.

    fluid : interFoam (original, serial or mpirun -np N ... -parallel)
    solid : ccx_preCICE -i dam_gate -precice-participant Solid
    """
    inter_foam = os.path.join(foam_env["FOAM_APPBIN"], "interFoam")
    if nprocs <= 1:
        fluid_cmd = inter_foam
    else:
        mpirun = shutil.which("mpirun") or "mpirun"
        parts = [mpirun]
        if root:
            parts.append("--allow-run-as-root")
        parts.append("-np %d" % nprocs)
        parts.append(inter_foam)
        parts.append("-parallel")
        fluid_cmd = " ".join(parts)
    if taskset:
        fluid_cmd = "taskset -c %s %s" % (taskset, fluid_cmd)
    return fluid_cmd, "ccx_preCICE -i dam_gate -precice-participant Solid"


def run_coupled(case_dir, fluid_cmd, solid_cmd, bashrc, nprocs,
                solid_threads, end_time, label="PRECICE 04"):
    """Launch solid + fluid as two background children, tail log.fluid for the
    progress bar, then wait for both and return (fluid_code, solid_code)."""
    run_dir = os.path.join(case_dir, "precice-run")
    os.makedirs(run_dir, exist_ok=True)
    for name in ("log.fluid", "log.solid"):
        path = os.path.join(run_dir, name)
        if os.path.isfile(path):
            os.remove(path)

    solid_dir = os.path.join(case_dir, "solid-calculix")
    fluid_dir = os.path.join(case_dir, "fluid-openfoam")
    fluid_log = os.path.join(run_dir, "log.fluid")

    # ---- solid background child -------------------------------------------
    solid_shell = foam_shell(
        "export OMP_NUM_THREADS=%d; "
        "export CCX_NPROC_EQUATION_SOLVER=%d; %s" %
        (solid_threads, solid_threads, solid_cmd), bashrc)
    solid_proc = subprocess.Popen(
        ["bash", "-c", solid_shell], cwd=solid_dir,
        stdout=open(os.path.join(run_dir, "log.solid"), "w",
                    encoding="utf-8", errors="replace"),
        stderr=subprocess.STDOUT)

    # ---- fluid background child --------------------------------------------
    fluid_shell = foam_shell(fluid_cmd, bashrc)
    fluid_proc = subprocess.Popen(
        ["bash", "-c", fluid_shell], cwd=fluid_dir,
        stdout=open(fluid_log, "w", encoding="utf-8", errors="replace"),
        stderr=subprocess.STDOUT)

    log("--- started solid (pid %d) + fluid (pid %d)"
        % (solid_proc.pid, fluid_proc.pid))
    log("--- logs: precice-run/log.solid, precice-run/log.fluid")

    # ---- monitor: fluid log drives the progress bar ------------------------
    pos = 0
    t0 = time.time()
    latest = {}
    _render_progress(label, 0.0, end_time, 0.0, t0)
    while True:
        time.sleep(1.0)
        f_code = fluid_proc.poll()
        s_code = solid_proc.poll()
        pos, _matched, new_vals = _tail_live(fluid_log, pos)
        if new_vals:
            latest.update(new_vals)
        wall = time.time() - t0
        cur = latest.get("time", 0.0)
        if f_code is not None or s_code is not None or latest.get("end"):
            time.sleep(0.2)
            _pos, _m2, nv = _tail_live(fluid_log, pos)
            if nv:
                latest.update(nv)
            wall = time.time() - t0
            cur = latest.get("time", 0.0)
            _render_progress(label, cur, end_time, wall, t0, final=True)
            break
        _render_progress(label, cur, end_time, wall, t0)

    fluid_code = fluid_proc.wait()
    solid_code = solid_proc.wait()
    log("=== coupled run finished: fluid exit %d, solid exit %d"
        % (fluid_code, solid_code))
    return fluid_code, solid_code


# --------------------------------------------------------------------------- #
#  Argument parsing + main                                                     #
# --------------------------------------------------------------------------- #

MODE_NAMES = {
    "f": "fresh (clean all, checkMesh + setFields + decomposePar, from startTime)",
    "r": "restart (keep mesh/decomposition, clean results, from startTime)",
    "c": "continue (keep all, from latestTime)",
}


def parse_args(argv):
    p = argparse.ArgumentParser(
        prog="04_run_precice.py",
        description="Run the preCICE partitioned dam-break FSI case "
                    "(interFoam + ccx_preCICE as two background children).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python3 04_run_precice.py --case /path/to/case -f\n"
               "  python3 04_run_precice.py --case /path/to/case -f --nprocs 8\n"
               "  python3 04_run_precice.py --case /path/to/case -f --end-time 0.5 --dry-run\n",
    )
    g = p.add_mutually_exclusive_group()
    g.add_argument("-f", "--fresh", dest="mode", action="store_const",
                   const="f", help="fresh start (default)")
    g.add_argument("-r", "--restart", dest="mode", action="store_const",
                   const="r", help="restart: keep mesh/decomposition, clean results")
    g.add_argument("-c", "--continue", dest="mode", action="store_const",
                   const="c", help="continue from latestTime")
    p.set_defaults(mode="f")

    p.add_argument("--case", required=True,
                   help="assembled preCICE case root (03 output)")
    p.add_argument("--end-time", type=float, default=None,
                   help="override fluid endTime in controlDict (default: short "
                        "0.02; pass 0.5 for the full run)")
    p.add_argument("--write-interval", type=float, default=None,
                   help="override writeInterval (runTime) in controlDict")
    p.add_argument("--nprocs", type=int, default=8,
                   help="fluid MPI processes (default 8)")
    p.add_argument("--solid-threads", type=int, default=8,
                   help="solid OMP_NUM_THREADS / CCX_NPROC_EQUATION_SOLVER "
                        "(default 8)")
    p.add_argument("--taskset", default=None,
                   help="optional CPU affinity mask for the fluid solver "
                        "(e.g. '0-7')")
    p.add_argument("--dry-run", action="store_true",
                   help="print plan + time estimate only; do not execute")
    return p.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    case_dir = os.path.abspath(args.case)
    dry = args.dry_run
    mode = args.mode
    nprocs = args.nprocs
    solid_threads = args.solid_threads

    if nprocs < 1:
        error("--nprocs must be a positive integer, got %d" % nprocs)
        return 2
    if solid_threads < 1:
        error("--solid-threads must be a positive integer, got %d" % solid_threads)
        return 2

    fluid_dir = os.path.join(case_dir, "fluid-openfoam")
    solid_dir = os.path.join(case_dir, "solid-calculix")
    ctrl_path = os.path.join(fluid_dir, "system", "controlDict")
    if not os.path.isfile(ctrl_path):
        error("no fluid-openfoam/system/controlDict in case: %s" % case_dir)
        return 2
    if not os.path.isfile(os.path.join(case_dir, "precice-config.xml")):
        error("no precice-config.xml in case: %s" % case_dir)
        return 2
    if not os.path.isfile(os.path.join(solid_dir, "dam_gate.inp")):
        error("no solid-calculix/dam_gate.inp in case: %s" % case_dir)
        return 2

    foam_env = resolve_foam_env()
    if not foam_env:
        error("could not locate OpenFOAM bashrc; cannot run interFoam")
        return 2

    root = (os.geteuid() == 0)
    if root:
        os.environ["OMPI_ALLOW_RUN_AS_ROOT"] = "1"
        os.environ["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] = "1"

    # ---- read configuration -------------------------------------------------
    ctrl = read_control_dict(case_dir)
    delta_t = ctrl["deltaT"]
    end_time = args.end_time if args.end_time is not None else ctrl["endTime"]
    if end_time <= 0:
        error("endTime must be positive (got %s)" % _fmt_num(end_time))
        return 2
    write_interval = (args.write_interval if args.write_interval is not None
                      else ctrl["writeInterval"])

    precice = read_precice_config(os.path.join(case_dir, "precice-config.xml"))
    if "window" in precice and abs(precice["window"] - delta_t) > 1e-12:
        warn("controlDict deltaT (%s) != precice-config time-window-size (%s); "
             "coupling windows will NOT match fluid steps"
             % (_fmt_num(delta_t), _fmt_num(precice["window"])))
    n_steps = max(1, int(round(end_time / delta_t)))

    n_cells = read_poly_mesh_cells(fluid_dir)
    est = estimate_time(n_cells, nprocs, end_time, delta_t)

    # ---- header -------------------------------------------------------------
    raw("[PRECICE-PIPE] ===== preCICE coupled run =====")
    log("case       : %s" % case_dir)
    log("mode       : %s (%s)" % (mode, MODE_NAMES[mode]))
    log("fluid      : interFoam, %s, deltaT %s, endTime %s, %d windows"
        % ("parallel (%d procs)" % nprocs if nprocs > 1 else "serial",
           _fmt_num(delta_t), _fmt_num(end_time), n_steps))
    log("solid      : ccx_preCICE, OMP/CCX %d threads"
        % solid_threads)
    if n_cells is not None:
        log("mesh       : %d fluid cells" % n_cells)
    log("estimate   : ~%s (endTime %s, %d windows, %d procs)"
        % (format_duration(est["seconds"]), _fmt_num(end_time), n_steps, nprocs))
    if root:
        log("--- Running as root: enabling OpenMPI root override")

    # ---- pre-run validation --------------------------------------------------
    for name, path in (("interFoam", os.path.join(foam_env["FOAM_APPBIN"], "interFoam")),
                       ("ccx_preCICE", shutil.which("ccx_preCICE") or "/usr/bin/ccx_preCICE")):
        if not os.path.isfile(path):
            error("%s not found: %s" % (name, path))
            return 2
    adapter = os.path.join(foam_env["FOAM_USER_LIBBIN"],
                           "libpreciceAdapterFunctionObject.so")
    if not os.path.isfile(adapter):
        error("libpreciceAdapterFunctionObject.so not found: %s" % adapter)
        return 2
    log("OK: interFoam, ccx_preCICE, libpreciceAdapterFunctionObject")

    def show(cmd):
        log("      [dry-run] %s" % cmd)

    # ---- mode-specific preparation -------------------------------------------
    if mode == "f":
        log("--- Fresh mode: cleaning fluid results + processor dirs")
        targets = collect_root_cleanup(fluid_dir) + \
            collect_processor_cleanup(fluid_dir, keep_mesh=False)
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
        log("--- Set startFrom = startTime%s" % ("" if not dry else " [dry-run]"))

        if nprocs > 1:
            dict_path = os.path.join(fluid_dir, "system", "decomposeParDict")
            action = "update" if os.path.isfile(dict_path) else "create"
            if dry:
                log("      [dry-run] decomposeParDict: %s numberOfSubdomains=%d"
                    % (action, nprocs))
            else:
                ensure_decompose_dict(fluid_dir, nprocs)
                log("--- decomposeParDict: %s (numberOfSubdomains=%d)"
                    % (action, nprocs))

        check_cmd = os.path.join(foam_env["FOAM_APPBIN"], "checkMesh")
        log("--- Running checkMesh -> fluid-openfoam/log.checkMesh")
        if dry:
            show(check_cmd)
        else:
            code, out = run_capture(check_cmd, fluid_dir, foam_env["bashrc"],
                                    log_file="log.checkMesh")
            mesh_ok = code == 0 and "Mesh OK" in out and "Failed" not in out
            if not mesh_ok:
                error("checkMesh failed; see fluid-openfoam/log.checkMesh")
                return 2
            log("--- checkMesh OK")

        setfields_cmd = os.path.join(foam_env["FOAM_APPBIN"], "setFields")
        log("--- Running setFields (water column) -> log.setFields")
        if dry:
            show(setfields_cmd)
        else:
            code, out = run_capture(setfields_cmd, fluid_dir, foam_env["bashrc"],
                                    log_file="log.setFields")
            if code != 0:
                error("setFields failed; see fluid-openfoam/log.setFields")
                return 2
            log("--- setFields OK")

        if nprocs > 1:
            decomp_cmd = (os.path.join(foam_env["FOAM_APPBIN"], "decomposePar")
                          + " -force")
            log("--- Running decomposePar -force -> log.decomposePar")
            if dry:
                show(decomp_cmd)
            else:
                code, out = run_capture(decomp_cmd, fluid_dir, foam_env["bashrc"],
                                        log_file="log.decomposePar")
                if code != 0:
                    error("decomposePar failed; see log.decomposePar")
                    return 2
                log("--- decomposePar OK")

    elif mode == "r":
        if nprocs > 1 and not has_processors(fluid_dir):
            error("restart mode needs processor* dirs; run -f first")
            return 2
        log("--- Restart mode: cleaning results, keeping mesh/decomposition")
        targets = collect_root_cleanup(fluid_dir) + \
            collect_processor_cleanup(fluid_dir, keep_mesh=True)
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
        log("--- Set startFrom = startTime%s" % ("" if not dry else " [dry-run]"))

    else:  # "c"
        if nprocs > 1 and not has_processors(fluid_dir):
            error("continue mode needs processor* dirs; run -f first")
            return 2
        log("--- Continue mode: keeping all results")
        warn("continue mode restarts the fluid solver from latestTime but the "
             "preCICE coupling always starts at time-window 0; the fluid 0/ "
             "fields must match the last coupled state for this to be valid")
        if not dry:
            set_start_from(case_dir, "latestTime")
        log("--- Set startFrom = latestTime%s" % ("" if not dry else " [dry-run]"))

    # ---- apply optional overrides (restored after the run) -------------------
    if args.end_time is not None:
        log("--- Override endTime = %s (restored after run)" % _fmt_num(end_time))
        if not dry:
            set_control_dict_value(case_dir, "endTime", end_time)
    if args.write_interval is not None:
        log("--- Override writeInterval = %s (restored after run)"
            % _fmt_num(write_interval))
        if not dry:
            set_control_dict_value(case_dir, "writeInterval", write_interval)

    # ---- preCICE max-time is the real coupling end control ---------------------
    # The preciceAdapter overrides the fluid controlDict endTime to infinity
    # ("Only preCICE will control the simulation's endTime"), so a short run is
    # achieved by patching <max-time> in precice-config.xml, not the endTime.
    # The original tag is restored after the run.
    config_path = os.path.join(case_dir, "precice-config.xml")
    orig_max_time = precice.get("max_time")
    log("--- precice-config.xml: <max-time> %s -> %s (adapter overrides fluid "
        "endTime to infinity; this is the true coupling stop)"
        % (_fmt_num(orig_max_time) if orig_max_time is not None else "?",
           _fmt_num(end_time)))
    orig_tag = None
    if dry:
        log("      [dry-run] patch precice-config.xml <max-time> to %s"
            % _fmt_num(end_time))
    else:
        orig_tag = patch_precice_max_time(config_path, end_time)
        if orig_tag is None:
            error("could not patch <max-time> in precice-config.xml; the "
                  "coupling would run to its configured max-time instead of %s"
                  % _fmt_num(end_time))
            return 2

    # ---- dry run: stop here --------------------------------------------------
    if dry:
        fluid_cmd, solid_cmd = build_solver_commands(foam_env, nprocs,
                                                     args.taskset, root)
        log("--- Dry run: coupled solvers would launch:")
        show("( cd %s && %s > precice-run/log.solid ) &" % (solid_dir, solid_cmd))
        show("( cd %s && %s > precice-run/log.fluid ) &" % (fluid_dir, fluid_cmd))
        if nprocs > 1:
            show(os.path.join(foam_env["FOAM_APPBIN"], "reconstructPar"))
        log("Dry run complete (no commands executed).")
        return 0

    # ---- launch the coupled run ----------------------------------------------
    fluid_cmd, solid_cmd = build_solver_commands(foam_env, nprocs,
                                                 args.taskset, root)
    log("--- Starting coupled run: interFoam (%s) + ccx_preCICE (%d threads)"
        % ("parallel" if nprocs > 1 else "serial", solid_threads))
    log("--- fluid: %s" % fluid_cmd)
    log("--- solid: %s" % solid_cmd)
    fluid_code, solid_code = run_coupled(
        case_dir, fluid_cmd, solid_cmd, foam_env["bashrc"], nprocs,
        solid_threads, end_time)

    try:
        if fluid_code != 0:
            error("interFoam failed (exit %d); see precice-run/log.fluid"
                  % fluid_code)
            return fluid_code
        if solid_code != 0:
            error("ccx_preCICE failed (exit %d); see precice-run/log.solid"
                  % solid_code)
            return solid_code

        # ---- reconstruct (parallel) -------------------------------------------
        if nprocs > 1:
            log("--- Running reconstructPar -> log.reconstructPar")
            recon_cmd = os.path.join(foam_env["FOAM_APPBIN"], "reconstructPar")
            code, out = run_capture(recon_cmd, fluid_dir, foam_env["bashrc"],
                                    log_file="log.reconstructPar")
            if code != 0:
                error("reconstructPar failed; see log.reconstructPar")
                return 2
            log("--- Reconstruction complete")

        # ---- restore optional overrides ----------------------------------------
        _restore_overrides(case_dir, args, ctrl)

        # ---- post-run verification ----------------------------------------------
        log("--- Verifying results")
        problems = verify_run(case_dir, end_time, delta_t)
        if problems:
            for pr in problems:
                error(pr)
            log("Run finished but verification reported issues (see above).")
            return 3

        exec_t, _clock_t = parse_execution_time(
            os.path.join(case_dir, "precice-run", "log.fluid"))
        log("endTime: %s  windows: %d  fluid ExecutionTime: %s"
            % (_fmt_num(end_time), n_steps,
               ("%ss" % format(exec_t, ".2f")) if exec_t is not None else "?"))
        log("Run complete. Logs: precice-run/log.fluid, precice-run/log.solid"
            "%s" % (", log.decomposePar, log.reconstructPar" if nprocs > 1 else ""))
        return 0
    finally:
        if orig_tag is not None:
            restore_precice_max_time(config_path, orig_tag)


def _restore_overrides(case_dir, args, ctrl):
    try:
        if args.end_time is not None:
            set_control_dict_value(case_dir, "endTime", ctrl["endTime"])
        if args.write_interval is not None:
            set_control_dict_value(case_dir, "writeInterval", ctrl["writeInterval"])
    except Exception as exc:
        error("failed to restore controlDict overrides: %s" % exc)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
