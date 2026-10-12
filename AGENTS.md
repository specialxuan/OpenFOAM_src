# PROJECT KNOWLEDGE BASE

**Generated:** 2026-10-02 CST
**Updated:** 2026-10-11 CST (HEAD `e4e729d`)
**Branch:** master

## OVERVIEW

OpenFOAM v2412 user-code repository containing a modal-superposition dynamic mesh library, shared FSI controls, custom incompressible/VOF/compressible solver families, and tracked case-automation pipelines (Damfailure, AGARD 445.6, SwingBeam). Build with OpenFOAM `wmake`; there is no repository-wide `Allwmake`, `Allclean`, or CI workflow.

## STRUCTURE

```text
src/
|-- fastDynamicFvMesh/      # dynamic mesh library, AMR/restart/force projection/modal extension, bilingual docs
|-- fsiCoupling/            # header-only FSI/PIMPLE convergence contract
|-- myInterFoam/            # two-phase VOF FSI solver, overset variant, three-phase child
|-- myPimpleFoam/           # incompressible FSI solver, overset and SRF variants
|-- myRhoPimpleFoam/        # compressible FSI solver family and unbuilt reference baseline
|-- myHisa/                 # renamed HiSA v1.12.7 copy with FSI modal-force criterion in its pseudotime loop
|-- case/                   # tracked case generators/pipelines (Python+Bash), NOT wmake code
|   |-- Damfailure/         # dam-break FSI: FDM and preCICE pipelines
|   |-- AGARD/              # AGARD 445.6 wing: Mode/ (modal calibration), FDM/ (case assembly), Mesh/ (deck converters)
|   `-- SwingBeam/          # closed-tank clamped cantilever beam FSI: FDM and preCICE pipelines
|-- DEVELOPMENT.md          # source/generated/local file policy
`-- AGENTS.md               # repository-wide rules
```

`scripts/` was untracked in `effbecd` (v3.3.9) and is gitignored; the live copy now
exists only outside git at `/root/Workspace/_Oldversion_/scripts/` (smoke runner,
mesh/modal converters, staged pipelines).

Nearly every major directory carries its own small `AGENTS.md` (`fastDynamicFvMesh/`,
`fsiCoupling/`, `myPimpleFoam/`, `myInterFoam/` plus `interMixingFoam/`, `myRhoPimpleFoam/`,
`myHisa/`, `case/Damfailure/`, `case/AGARD/` plus `Mode/`, `case/SwingBeam/`). Per-directory
detail (invariants, stage contracts, verified counts) lives there; this file holds only
repo-wide rules.

## WHERE TO LOOK

| Task                       | Location                                                                     | Notes                                                          |
| -------------------------- | ---------------------------------------------------------------------------- | -------------------------------------------------------------- |
| Dynamic mesh runtime       | `fastDynamicFvMesh/fastDynamicFvMesh/`                                       | `update()`, AMR, restart, force assembly, modal extension, diagnostics |
| Mesh library docs          | `fastDynamicFvMesh/README.md`, `fastDynamicFvMesh/Doc/快速动网格程序手册.md` | Must stay synced with library code                             |
| Shared FSI controls        | `fsiCoupling/fsiCoupling.H`                                                  | `fsiCoupling`, `maxFsiIter`, `fsiTolerance`, `fsiDispTolerance`, `fsiDispConsecutiveIters`, `fsiResidual` |
| Incompressible FSI solver  | `myPimpleFoam/`                                                              | Main solver plus overset and SRF variants                      |
| VOF FSI solver             | `myInterFoam/`                                                               | Two-phase main, overset variant, nested three-phase executable |
| Compressible FSI solver    | `myRhoPimpleFoam/`                                                           | Main solver, equation headers, overset variant, baseline copy  |
| HiSA with FSI coupling     | `myHisa/`                                                                    | HiSA copy; FSI keys (incl. displacement OR gate) honored in the pseudotime loop; builds via its own `Allwmake` |
| Damfailure case pipelines  | `case/Damfailure/FDM/`, `case/Damfailure/preCICE/`                           | Numbered Python stages driven by `run_pipeline.sh`             |
| AGARD modal workflow       | `case/AGARD/Mode/`                                                           | Detailed `AGENTS.md`; LFS mesh decks under `Mode/mesh/`        |
| AGARD FDM case assembly    | `case/AGARD/FDM/assemble_fdm_case.py`                                        | Wires converted polyMesh + modal CSVs into a myRhoPimpleFoam + fastDynamicFvMesh case |
| AGARD mesh converters      | `case/AGARD/Mesh/`                                                           | `inp_to_polymesh.py` (CalculiX/Abaqus deck → `constant/polyMesh`), build/smooth tools |
| AGARD validation spec      | `AGARD445.6验证算例规范.md` (repo root, untracked)                            | Conditions/targets/plots for AMR + dynamic-mesh validation; core case Ma 0.96 |
| SwingBeam case pipelines   | `case/SwingBeam/FDM/`, `case/SwingBeam/preCICE/`                             | Closed tank + clamped cantilever beam; six-stage Python stages |
| Smoke validation           | `/root/Workspace/_Oldversion_/scripts/smoke_fastDynamicFvMesh.sh`            | External (untracked); currently blocked, see CONVENTIONS       |
| Mesh/modal/FSI pipeline    | `/root/Workspace/_Oldversion_/scripts/pipeline/`                             | External; staged automation with mesh/mode consistency guard   |
| File policy                | `DEVELOPMENT.md`, `.gitignore`                                               | Source vs generated vs local-only rules                        |

## CODE MAP

No LSP/codegraph tool is wired up for OpenFOAM `.C/.H` in this harness, but each solver
root tracks a clangd `compile_commands.json` (`fastDynamicFvMesh/`, `fsiCoupling/`,
`myInterFoam/`, `myPimpleFoam/`, `myRhoPimpleFoam/`) usable with bear-style tooling.
Grep-verified central nodes:

| Symbol/File                                             | Type     | Location                                                               | Role                                                           |
| ------------------------------------------------------- | -------- | ---------------------------------------------------------------------- | -------------------------------------------------------------- |
| `fastDynamicFvMesh`                                     | class    | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMesh.H`              | `dynamicRefineFvMesh` subclass, modal state, AMR state, caches |
| `fastDynamicFvMesh::update()`                           | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshMain.C`          | Central polymorphic entry; invoked as `mesh.controlledUpdate()` from `pimpleLoopBody.H` (myPimpleFoam, myRhoPimpleFoam) and as `mesh.update()` from `myInterFoam/pimpleLoopBody.H` and `myHisa` `hisaModule.C` |
| `assembleAndRelaxModalForces()`                         | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshUpdate.C`        | Force relaxation and FSI residual publication path             |
| `readControls()`                                        | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshIO.C`            | Runtime dictionary parsing and compatibility checks            |
| `readModeShapes()` / `ensureModeShapesForCurrentMesh()` | methods  | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshModeMapping.C`   | Mode CSV mapping, AMR remap, shared-point sync                 |
| `calcModalForces()`                                     | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshForces.C`        | Pressure/shear/structural force projection                     |
| `solveStructuralDynamics()`                             | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshSolver.C`        | Wilson-Theta modal dynamics                                    |
| `fastDynamicFvMeshModalExtension.C`                     | source   | `fastDynamicFvMesh/fastDynamicFvMesh/`                                 | Modal extension: structural force projection plus displacement/stress/strain recovery snapshots; `modalExtensionEnabled`/`StressEnabled`/`StrainEnabled` all default off |
| `fsiPimpleControl`                                      | class    | `fsiCoupling/fsiCoupling.H`                                            | PIMPLE criteria plus FSI convergence                           |
| `readFsiResidual()`                                     | function | `fsiCoupling/fsiCoupling.H`                                            | Solver-side force-residual read for FSI convergence (with `fsiDisplacementChangeValid()` / `fsiMaxDisplacementChange()` forming the displacement OR gate); used by all four FSI executables (three my*Foam mains + myHisa `pseudotimeControl.C`) |
| `myPimpleFoam.C`                                        | main     | `myPimpleFoam/`                                                        | exe `$(FOAM_USER_APPBIN)`; mesh update inside shared `pimpleLoopBody.H` (included 3× per main) |
| `myInterFoam.C`                                         | main     | `myInterFoam/`                                                         | exe `$(FOAM_USER_APPBIN)`; own `pimpleLoopBody.H` clears stale mesh compression flux on topology change |
| `myRhoPimpleFoam.C`                                     | main     | `myRhoPimpleFoam/`                                                     | exe `$(FOAM_USER_APPBIN)`; branches share `pimpleLoopBody.H`   |
| `interMixingFoam.C`                                     | main     | `myInterFoam/interMixingFoam/`                                         | exe `$(FOAM_APPBIN)`; three-phase VOF with local mixture models |
| `rhoPimpleFoam.C`                                       | main     | `myRhoPimpleFoam/`                                                     | NOT built (no Make entry); upstream baseline reference copy    |

`fsiCoupling.H` is included by the three main FSI solvers (`myPimpleFoam`,
`myInterFoam`, `myRhoPimpleFoam`) and by `myHisa` (`singleRegionSolver.H`,
`hisaModule.C`, `pseudotimeControl.C`); no overset/SRF/three-phase variant carries
FSI wiring.

## CONVENTIONS

- Use `wmake` only. Main targets: `wmake libso fastDynamicFvMesh`, `wmake myPimpleFoam`, `wmake myInterFoam`, `wmake myRhoPimpleFoam`. myHisa builds via its own `myHisa/Allwmake` (src then applications), not a single top-level `wmake`.
- Build nested variants from their own directories with `wmake`; each owns separate `Make/files` and `Make/options` and never inherits the parent's.
- OpenFOAM environment must already be loaded. Main targets install to `$FOAM_USER_LIBBIN`/`$FOAM_USER_APPBIN`; nested variants (`interMixingFoam`, `over*DyMFoam`, `SRFPimpleFoam`) install to `$FOAM_APPBIN` (system tree — intentional, know before you `wmake`).
- Run `ccx` and OpenFOAM tools (`blockMesh`, `checkMesh`, solvers, etc.) only inside a case directory under `/root/Workspace/...` — never in the `src/` repository root. Running them in `src/` drops `.sta/.dat/.cvg` and other artifacts into the repo and pollutes it.
- Scratch files, temporary extraction, tool output, and intermediate files must go under `/root/Workspace/tmp/` — never in the system `/tmp`.
- Baseline smoke archives: the smoke script and `DEVELOPMENT.md` expect `/root/Workspace/case_damfailure.zip` and `case_transient.zip`. **Neither currently exists at `/root/Workspace/`** (only `case_transient.zip` survives under `_Oldversion_/other_zip/`), so the smoke script aborts at its archive check — even with `--dry-run` — until `case_damfailure.zip` is restored there. Never edit baseline archives directly; copy them to temp run dirs before testing.
- Pipeline drivers (`case/*/*/run_pipeline.sh`) support `--dry-run`, `--steps` for partial runs, and `--skip-if-fresh` (mtime gate on generation stages only; solver/post/video stages always run). Default case roots point under `/root/Workspace/...`, not in-tree.
- Video export stages (`06_export_*_video.py`) drive ParaView (`paraview.simple` — needs `pvbatch`/ParaView Python); other pipeline stages are stdlib-only Python.
- New dictionary keys need safe defaults (`readIfPresent`) preserving current behavior. Removing/renaming keys requires a deprecation warning path (canonical example: `refinementFaceMapTolerance` in `fastDynamicFvMeshIO.C`).
- FSI convergence is force-residual based: `fsiResidual < fsiTolerance`; opt-in OR absolute-displacement gate `fsiDispTolerance > 0` (default off). When enabled, the gate must hold for `fsiDispConsecutiveIters` consecutive outer iterations (default `1` = legacy single-step) to reject one-step stagnation false positives. It terminates the partitioned loop directly; in integrated/implicit mode it acts only through the combined PIMPLE/pseudotime AND-veto (so it can stop the loop at the residual criterion—shorter than waiting for force convergence—but cannot terminate the loop on its own). Displacement change is otherwise diagnostic-only. The three `my*Foam` mains wire the gate through `fsiPimpleControl`; myHisa implements the same AND veto inside `pseudotimeControl::criteriaSatisfied()`.
- Changes under `fastDynamicFvMesh/fastDynamicFvMesh/` must update both English and Chinese docs, build the library, and run a smoke path.
- Runtime-control, file-format, or workflow changes must update example case configs and `dynamicMeshDict` guidance.
- `case/AGARD/**/*.zip` is Git LFS (`.gitattributes`); never re-add those archives as plain blobs.

## MACHINES & RUN PLACEMENT

|                       | Local (dev container)                         | Remote `FSIDocker`                       |
| --------------------- | --------------------------------------------- | ---------------------------------------- |
| Form                  | WSL2 container                                | Docker container on WSL2                 |
| CPU                   | i7-14700K (8P+12E); 20 logical cores visible  | 64 cores, no cgroup quota                |
| RAM                   | 31 GB (~27 usable)                            | 251 GB                                   |
| OpenFOAM              | v2412; hisa 1.12.7 + my*Foam built            | v2412; hisa + my*Foam built              |
| Tools                 | full (taskset available)                      | no rsync/tmux/screen; git/unzip/python3  |

- Local parallel runs: **exactly 8 MPI ranks** (= the 8 P-cores). More ranks spill onto E-cores and run slower than 8. WSL2 does not expose P/E topology (`core_type` absent), so do not blind-pin with `taskset`; use `numberOfSubdomains 8;` and `export OMP_NUM_THREADS=1`, and run only one parallel case at a time.
- Send to remote when ANY holds: > ~3M cells, estimated memory > 20 GB, estimated wall time > 12 h, or ranks > 8. The AGARD combined CalculiX modal decks (~37 GB peak RSS) are always remote-class.
- When unsure: run 100–500 steps locally on a coarse mesh / short `endTime`, extrapolate `ExecutionTime`, then decide.
- Remote workflow (run and monitored via the `openfoam-run-monitor` subagent; full procedure: `.opencode/skills/remote-compute-polling/`): resolve the SSH alias from `~/.ssh/config` into a variable — never bake hostnames/container ids into scripts or case files; transfer with `tar czf - -C dir . | ssh "$SSH_HOST" 'tar xzf - -C dest'` (no remote rsync) and verify `sha256sum` on both ends; launch with `setsid nohup ... > run.log 2>&1 < /dev/null &` wrapped in local `timeout` (ssh exit 124 is expected); verify in a separate ssh call; judge liveness with `ps -C <proc> -o stat= | grep -v '^Z'` — the remote PID 1 is `sleep`, so `<defunct>` launchers linger for days without meaning a job is alive; poll every 30 min on a fixed cadence.
- **Every case run must go through the `openfoam-run-monitor` subagent, launched in the background** — solvers, `mpirun`, `decomposePar`/`reconstructPar`, pipeline stages, smoke runs, and long builds: call `subagent(agent="openfoam-run-monitor", background=true)`. Do **not** run these in the foreground, and do **not** launch them as a bare `shell` background job: that shell is reaped when the session ends (observed: a 62-min HiSA run died mid-step at t=0.023 with no solver error). The subagent detaches with `setsid nohup` and keeps a run registry under gitignored `.omo/openfoam-runs/`.
- **Supervise long runs with a `sleep` poll loop**, not by returning early: poll about every 30 min, count only non-zombie processes (`ps -C <proc> -o stat= | grep -v '^Z'`), stop only on a terminal log state (`End`) or a fatal, and emit one compact status line per poll. Remote procedure: `.opencode/skills/remote-compute-polling/`.
- Record hostname, host load, rank count, and solver version for every run; performance comparisons are valid only when both machines were idle.
- `decomposePar` with `method scotch` is **non-reproducible in this environment** (re-running on identical input yields different partitions, hence different rounding). For bit-for-bit A/B comparisons, decompose once and reuse the same `processor*/`, or pin a deterministic decomposition method.

## ANTI-PATTERNS (THIS PROJECT)

- Do not use cmake/make or create parallel build systems.
- Do not edit baseline case archives directly.
- Do not commit `*/lnInclude/`, `*/linux64GccDPInt32Opt/`, `*/.cache/`, `__pycache__/`, `.venv_fdm/`, `.codex`, `.omo/`, `.vscode/mcp.json`, `*.cvg`, or `*.sta`.
- Do not rely on `refinementFaceMapTolerance` as live behavior; it is deprecated, warned, and ignored.
- Do not use geometric fallback for AMR-generated points; topology lineage/mapping rules are mandatory (code-enforced `FatalError` in `fastDynamicFvMeshModeMapping.C`).
- Do not reuse previous time-step mesh compression flux after topology changes.
- Do not claim AMR is unconditionally stable or that refinement can never degrade mesh quality.
- Do not treat stale CalculiX outputs, stale mode files, or old FSI logs as production inputs.
- Do not sync partial `reconstructPar` output into a main case after failure.
- Do not assume `scripts/` exists in-tree; it is untracked. Use `/root/Workspace/_Oldversion_/scripts/` or restore deliberately.

## COMMANDS

```bash
wmake libso fastDynamicFvMesh
wmake myPimpleFoam
wmake myInterFoam
wmake myRhoPimpleFoam
(cd myHisa && ./Allwmake)        # HiSA copy: builds src/ then applications/
case/Damfailure/FDM/run_pipeline.sh --dry-run
case/Damfailure/preCICE/run_pipeline.sh --dry-run
case/AGARD/Mode/run_pipeline.sh --dry-run
case/SwingBeam/FDM/run_pipeline.sh --dry-run
case/SwingBeam/preCICE/run_pipeline.sh --dry-run
/root/Workspace/_Oldversion_/scripts/smoke_fastDynamicFvMesh.sh --dry-run   # currently aborts: case_damfailure.zip missing (see CONVENTIONS)
```

## NOTES

- The external smoke script checks both baseline zips before doing anything (a missing zip aborts even `--dry-run`), builds the library, extracts the zips under `RUN_ROOT`, reads solver applications from `system/controlDict`, and logs under `RUN_ROOT/logs`.
- For AMR/restart runs, prefer reconstructing in an isolated mirror case before syncing reconstructed times back into the main case.
- The external `run_mesh_independence_automation.sh` (`_Oldversion_/scripts/pipeline/`) refuses FSI if `mode/FluidNodeCoor.csv` node count differs from mesh point count unless `--allow-stale-mode` is used for diagnostics.
- No executable unit/regression test suite exists for the C++ code; validation is via builds, smoke cases, pipeline checks, and documented numerical comparisons. `case/AGARD/Mode/tests/`, `case/Damfailure/FDM/test_display_width.py`, and `case/SwingBeam/FDM/test_display_width.py` are Python-only unittest suites.
- Diagram workflow: edit `.drawio` sources first, keep final PNG/JPG only when needed, delete temporary SVG/debug exports.
