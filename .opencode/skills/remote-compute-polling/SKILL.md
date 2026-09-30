---
name: remote-compute-polling
description: "Run a long job (solver/CFD/FEA/batch) on a remote container or server over SSH, launched detached and monitored by periodic polling instead of babysitting. Covers SSH target resolution without hardcoding aliases, read-only pre-flight recon of hardware limits and host load, safe detached launch, 30-minute poll loops with zombie-aware liveness and crash detection, artifact/metric collection, and fair comparison of before/after performance. Use when asked to run/deploy a case on a remote container, '远程计算', 'run it on the remote box over ssh', 'poll every 30 minutes', '轮询', 'make sure it finishes', 'compare efficiency against local run', or any job that outlives a single tool call. Triggers: remote run, remote container, ssh compute, detached job, poll loop, 30-min polling, 远程计算, 轮询, efficiency comparison."
---

# Remote Compute + Polling

Run long jobs on a remote SSH target, detach them so they survive, and poll at a fixed
cadence until they finish. Designed for solver runs (OpenFOAM/CalculiX/preCICE/etc.) but
applies to any long batch job.

## Hard rules

1. **Never hardcode the SSH host or container name.** Resolve it at runtime from
   `~/.ssh/config` (or the name the user says), keep it in a variable, and never let the
   container identity leak into scripts, configs, or saved files.
2. **Never babysit.** Detach the job, then poll on an explicit schedule (default: every
   30 minutes). One poll = one compact status line.
3. **Read-only recon on any host you do not own.** Checking is fine; killing, renaming,
   reprioritizing, or changing another machine's processes is not.
4. **Long jobs are only comparable under equal host conditions.** Record the machine
   identity and its load before/after; a busy host invalidates any comparison.

## Environment conventions (same image, different instances)

Instances are normally launched from the **same image**, so tools, programs, and source
paths are identical across instances; only the instance identity changes.

- Default remote workspace root: **`/root/Workspace/`** — place cases and archives there.
  Override with an env var (`WORKSPACE=<path>`) only when a job needs a different root.
- Paths are stable, so the full layout does not need rediscovering every run — but verify
  once per instance: `hostname` (provenance fingerprint) plus that the expected paths exist.
- Never bake the instance identity (hostname / container id / SSH alias) into scripts,
  configs, or committed files. Pass it as a variable; record it only in run notes.
- If the image may have drifted (newer build than last time), re-run the Step 1 recon.
- Minimal images commonly lack `rsync`, `tmux`, `screen`. Probe with `command -v` before
  designing around any of them (see Steps 2 and 5).

## Step 0 — Resolve the SSH target (no hardcoding)

```bash
grep -E '^Host ' ~/.ssh/config            # list candidate aliases
SSH_HOST="<alias-chosen-by-user>"         # keep in a variable, reuse everywhere
ssh -o BatchMode=yes -o ConnectTimeout=10 "$SSH_HOST" \
  'hostname; uname -a 2>/dev/null; whoami'
```

- If the user names a target ("the remote container", "the big box"), map it to an alias
  from `~/.ssh/config` and confirm if more than one candidate fits.
- Record the returned `hostname`/`uname` — it is your fingerprint for later comparisons.
- The shell on the far side may be Windows (cmd/PowerShell) or Linux. Probe before
  assuming: read its `uname`, `hostname`, and the encoding of its output.

## Step 1 — Read-only pre-flight recon

Learn the machine's speed class and limits **before** interpreting any timing.

```bash
ssh "$SSH_HOST" 'nproc; lscpu | grep -E "^Model name|^CPU\(s\)|^Thread"; free -h;
  df -h /; cat /etc/os-release | head -3'
```

Container-specific limits (a "64-core" container may be capped):

```bash
ssh "$SSH_HOST" 'cat /sys/fs/cgroup/cpu.max; cat /sys/fs/cgroup/cpu.stat | head -4;
  cat /proc/pressure/cpu; grep "^cpu " /proc/stat'
```

- `cpu.max = max 100000` ⇒ no quota; a finite first field ⇒ CPU cap.
- `nr_throttled > 0` ⇒ the container is being throttled.
- `steal` (8th `/proc/stat` field) > 0 ⇒ hypervisor is taking CPU away.
- A quick single-core benchmark gives a portable speed score; run the same loop locally
  for a baseline ratio. Also run it N-way in parallel to see effective parallel throughput.

If the physical host is reachable (and only read-only):

```bash
# Linux host
ssh "$HOST" 'uptime; ps -eo pcpu,pid,comm --sort=-pcpu | head -15'
# Windows host (PowerShell, read-only)
#   total CPU, top processes, and per-process CPU over a sample window
```

Compute effective load: sum per-process CPU deltas over a 10–15 s sample, or use the OS
counter API. Compare against total physical cores. **A host running another big job is the
single most common cause of "impossible" performance regressions.**

## Step 2 — Transfer artifacts safely

Default is `scp` of a single archive — zero dependencies, present wherever OpenSSH is:

```bash
sha256sum local/archive.zip
ssh "$SSH_HOST" 'sha256sum /root/Workspace/archive.zip'   # identical? skip the copy
scp -o BatchMode=yes local/archive.zip "$SSH_HOST":/root/Workspace/
ssh "$SSH_HOST" 'sha256sum /root/Workspace/archive.zip'   # verify after copy
```

Choosing the transfer method:

| Situation | Use | Why |
|---|---|---|
| One archive (the common case) | `scp` | simplest, single integrity check |
| Many small files / whole case tree | `tar czf - -C dir . \| ssh "$SSH_HOST" 'tar xzf - -C /root/Workspace/dest'` | no dependency, one stream, no per-file round trips |
| Repeated incremental sync, resume, `--delete` | `rsync -avP --partial` | the only case that justifies rsync |
| Text-heavy payload (ascii mesh, CSV) | add compression: `scp -C` or the `tar czf` form | large ratio win; skip it for already-zipped archives |

- `rsync` must exist on **both** ends: `command -v rsync` locally and remotely. Minimal
  images usually lack it, and installing per instance is wasted work on ephemeral
  containers — bake it into the image if the workflow depends on it.
- Prefer transferring an archive over the unpacked tree: fewer round trips, one hash to check.
- Always verify by hash on both ends. Re-copying an identical archive wastes time and
  confuses provenance; a matching hash means the transfer step is already done.

## Step 3 — Prepare an isolated run directory

```bash
ssh "$SSH_HOST" 'bash -s' <<'EOF'
set -e
mkdir -p /remote/run_root
cd /remote/run_root
rm -rf case_dir
unzip -q ../archive.zip                  # fresh extract, never patch a dirty tree
cp /remote/run_root/launcher_patched.sh case_dir/launcher.sh   # patched copy lives OUTSIDE
chmod +x case_dir/launcher.sh
EOF
```

- Extract fresh for every run so results are attributable.
- Keep patched launchers **outside** any directory that cleanup will delete.
- Strip stale decomposed/parallel state (`processor*`, caches) before a fresh run.
- If a job writes to a fixed scratch/socket path (e.g. a coupling exchange dir), clear it
  per run, or give each run its own path so parallel attempts cannot collide.

## Step 4 — Verify runtime prerequisites (the silent-load trap)

A job can start, run to completion, and produce **invalid** results because a required
library/plugin silently failed to load. Check for the positive signal, not just absence of
errors.

```bash
grep -m1 "<expected loaded banner>"  <run>/run.log   # e.g. an adapter/plugin load line
grep -m1 "Could not load\|Unknown .* type" <run>/run.log   # the failure signature
```

Typical cause: a library exists in one install tree but the process's search path points
at another (`LD_LIBRARY_PATH`, module path, plugin dir). Verify explicit paths, or set the
search path in the launcher right after the environment is sourced.

## Step 5 — Launch detached (ssh must not hang)

`ssh` does not return while a background child holds the channel. Two rules:

```bash
timeout 10 ssh -o BatchMode=yes "$SSH_HOST" \
  "cd /remote/run_root/case_dir && ENV1=v1 ENV2=v2 \
   setsid nohup bash launcher.sh > run.log 2>&1 < /dev/null & echo launch-issued"
echo "ssh-exit=$?"      # 124 = local timeout fired; the remote job is still alive (by design)
```

1. `setsid nohup … < /dev/null > run.log 2>&1 &` fully detaches.
2. Wrap the ssh call in a local `timeout`; treat exit 124 as success.
3. **Always verify in a separate ssh call**, never in the same one.

Verification checklist (one ssh call, ~2 minutes after launch):

```bash
ssh "$SSH_HOST" 'C=/remote/run_root/case_dir; \
  echo "live=$(ps -C <proc> -o stat= | grep -vc "^Z")"; \
  grep -m1 "<loaded banner>" $C/run.log; \
  grep "<progress pattern>" $C/run.log | tail -1; \
  grep -c "<fatal pattern>" $C/run.log'
```

Confirm: process count matches the requested parallelism, the loaded banner is present,
progress markers advance, and the fatal counter is zero.

### Optional: tmux instead of setsid + nohup

`setsid nohup` is the default for a poll-driven workflow: no dependency, logs land in
files, and Step 6 just greps those files. Reach for `tmux` only when you need a human to
watch or type into the session, side-by-side panes for several participants, scrollback,
or a real TTY (some programs change buffering / progress output without one).

```bash
command -v tmux || echo "tmux not installed (bake it into the image if you want it)"
tmux new-session -d -s run "cd /root/Workspace/case_dir && bash launcher.sh 2>&1 | tee run.log"
tmux capture-pane -p -t run | tail -40      # inspect without attaching
tmux ls; tmux kill-session -t run           # state / stop
```

- Both approaches survive an ssh disconnect; neither survives the container being recreated.
- Keep writing a log file (`tee run.log`) even under tmux: parsing files is more robust than
  scraping panes, so the Step 6 poll loop stays unchanged.
- No tmux and you need a TTY: `script -qec "cmd" run.log` (present in Ubuntu images).
- Install once, in the image, if the team standardizes on tmux — not per ephemeral instance.

## Step 6 — Poll on a fixed cadence (default 30 min)

Run this as a single long-blocking local command; end your turn so it runs unattended.

```bash
while true; do
  sleep 1800                                   # 30 minutes
  out=$(ssh -o BatchMode=yes -o ConnectTimeout=20 "$SSH_HOST" 'bash -s' 2>/dev/null <<'EOF'
LOG=/remote/run_root/case_dir/run.log
# Zombie-aware liveness: pgrep also matches <defunct>, which causes false "still running"
LIVE=$(ps -C <proc> -o stat= 2>/dev/null | grep -v '^Z' | wc -l)
W=$(grep -c '<window/task completed pattern>' "$LOG" 2>/dev/null); W=${W:-0}
F=$(grep -c '<fatal pattern>' "$LOG" 2>/dev/null); F=${F:-0}
if [ "$LIVE" -gt 0 ]; then
  echo "STATUS=RUNNING live=$LIVE done=$W fatal=$F"
  grep '<progress pattern>' "$LOG" | tail -1
  grep '<cumulative-cost pattern>' "$LOG" | tail -1
elif [ "$F" != "0" ]; then
  echo "STATUS=CRASHED fatal=$F"
  grep -m1 -A2 '<fatal pattern>' "$LOG"
else
  echo "STATUS=FINISHED done=$W"
  tail -3 /remote/run_root/case_dir/run.log
  grep '<final-cost pattern>' "$LOG" | tail -1
fi
EOF
)
  echo "[$(date '+%F %T')]"; echo "$out"; echo "---"
  case "$out" in *"STATUS=FINISHED"*|*"STATUS=CRASHED"*) break;; esac
done
```

Rules for a trustworthy poll:

- **Liveness must exclude zombies** (`ps -C name -o stat= | grep -v '^Z'`). `pgrep`/`pgrep -x`
  match `<defunct>` processes, which linger when the parent (e.g. `mpirun`) died — you will
  report "RUNNING" for hours on a dead job.
- Detect terminal states from the **log** as well as the process table (a crash can leave
  orphans); break the loop only on FINISHED/CRASHED.
- Print one compact block per poll: timestamp, status, progress, cumulative cost.
- If progress is frozen across 2–3 consecutive polls, treat it as hung and investigate
  (process state `Z`, blocked syscall, partner died) rather than waiting longer.

## Step 7 — Stopping jobs and cleanup (without shooting yourself)

```bash
pkill -9 -x <procname>              # exact process name — safe
pkill -9 -f "scrip[t].sh"           # bracket trick: the regex no longer matches your own
                                    # ssh command line, which DOES contain the literal name
```

- `pkill -f "script.sh"` sent over ssh **kills your own ssh shell**, because that shell's
  command line contains the pattern. Symptom: the command returns no output at all.
- Prefer PID-based kills for the launcher; leave unrelated users' processes alone.
- Delete only the artifacts the user asked to delete; keep reference runs intact.

## Step 8 — Collect metrics and compare fairly

Extract a normalized cost per unit of work, not just total wall clock:

| Metric | Why |
|---|---|
| wall clock | headline, but host-load sensitive |
| in-solver clock time (e.g. `ClockTime`) | excludes setup/staging |
| in-solver CPU time (e.g. `ExecutionTime`) | **can include busy-wait during coupling — do not read it as pure compute** |
| seconds per task/window | normalizes different amounts of work |
| **seconds per iteration** | best single comparison metric |
| iteration/step count and cap hits | convergence quality — a run that hits an iteration cap did not converge, so its "faster" time is meaningless |
| artifact hashes/curves | result equivalence, not just speed |

Comparability checklist before drawing conclusions:

- [ ] Same machine identity (compare recorded `hostname`). Different container, same host
      can still be comparable — say so explicitly, but prove it from evidence.
- [ ] Host was equally idle for both runs. Check host load/other jobs and their start times.
- [ ] Same input, config, and decomposition/parallel layout, or the difference is explained.
- [ ] Convergence quality matched (no iteration-cap hits, comparable residual history).
- [ ] Any post-processing/IO excluded or matched.

## Pitfall summary (all observed in practice)

1. `ssh` hangs forever when a background child keeps the channel → `timeout` + `setsid nohup` + fd redirection.
2. Libraries silently not loading → the job "succeeds" with wrong physics; grep for the load banner.
3. `pgrep` matches zombies → false "RUNNING" for hours; use `ps -C` + zombie filter.
4. `pkill -f "name"` over ssh kills your own shell → bracket trick or PIDs.
5. Huge unexplained slowdowns → check the **physical host** for other big jobs (read-only) before blaming config.
6. Host resource caps → check cgroup `cpu.max`/throttling, `steal`, and PSI pressure, not just `nproc`.
7. Comparing runs across busy/idle host windows → invalid; re-run or wait for a quiet window.
8. Trusting "more threads = faster" → verify with a short A/B run; oversubscription and
   per-rank thread leakage (a thread-count env var exported globally and inherited by MPI
   ranks) can invert the expected result.
9. Env leakage between participants: exporting a thread count for one participant can
   silently apply to all of them. Set per-participant env explicitly.
10. Designing around a tool the image does not ship (`rsync`, `tmux`, `screen`): probe with
    `command -v` on both ends first; prefer the zero-dependency path (`scp`, `tar | ssh`,
    `setsid nohup`) and add optional tools to the image, not to each ephemeral instance.
11. Hardcoding the instance identity (SSH alias, hostname, container id) into scripts or
    committed files: instances change, paths do not. Keep identity in variables/run notes.
