# PROJECT KNOWLEDGE BASE

**Generated:** 2026-10-02 CST
**Commit:** 31b9881
**Branch:** master

## OVERVIEW

OpenFOAM v2412 user-code repository containing a modal-superposition dynamic mesh library, shared FSI controls, custom incompressible/VOF/compressible solver families, and tracked case-automation pipelines (Damfailure, AGARD 445.6). Build with OpenFOAM `wmake`; there is no repository-wide `Allwmake`, `Allclean`, or CI workflow.

## STRUCTURE

```text
src/
|-- fastDynamicFvMesh/      # dynamic mesh library, AMR/restart/force projection, bilingual docs
|-- fsiCoupling/            # header-only FSI/PIMPLE convergence contract
|-- myInterFoam/            # two-phase VOF FSI solver, overset variant, three-phase child
|-- myPimpleFoam/           # incompressible FSI solver, overset and SRF variants
|-- myRhoPimpleFoam/        # compressible FSI solver family and unbuilt reference baseline
|-- myHisa/                 # renamed HiSA v1.12.7 copy with FSI modal-force criterion in its pseudotime loop
|-- case/                   # tracked case generators/pipelines (Python+Bash), NOT wmake code
|   |-- Damfailure/         # dam-break FSI: FDM and preCICE pipelines
|   `-- AGARD/              # AGARD 445.6 wing: modal calibration, mesh tools, LFS mesh decks
|-- DEVELOPMENT.md          # source/generated/local file policy
`-- AGENTS.md               # repository-wide rules
```

`scripts/` was untracked in `effbecd` (v3.3.9) and is gitignored; the live copy now
exists only outside git at `/root/Workspace/Oldversion/scripts/` (smoke runner,
mesh/modal converters, staged pipelines).

## WHERE TO LOOK

| Task                       | Location                                                                     | Notes                                                          |
| -------------------------- | ---------------------------------------------------------------------------- | -------------------------------------------------------------- |
| Dynamic mesh runtime       | `fastDynamicFvMesh/fastDynamicFvMesh/`                                       | `update()`, AMR, restart, force assembly, diagnostics          |
| Mesh library docs          | `fastDynamicFvMesh/README.md`, `fastDynamicFvMesh/Doc/快速动网格程序手册.md` | Must stay synced with library code                             |
| Shared FSI controls        | `fsiCoupling/fsiCoupling.H`                                                  | `fsiCoupling`, `maxFsiIter`, `fsiTolerance`, `fsiDispTolerance`, `fsiResidual` |
| Incompressible FSI solver  | `myPimpleFoam/`                                                              | Main solver plus overset and SRF variants                      |
| VOF FSI solver             | `myInterFoam/`                                                               | Two-phase main, overset variant, nested three-phase executable |
| Compressible FSI solver    | `myRhoPimpleFoam/`                                                           | Main solver, equation headers, overset variant, baseline copy  |
| HiSA with FSI coupling     | `myHisa/`                                                                    | HiSA copy; FSI keys honored in pseudotime loop; has AGENTS.md  |
| Damfailure case pipelines  | `case/Damfailure/FDM/`, `case/Damfailure/preCICE/`                           | Numbered Python stages driven by `run_pipeline.sh`             |
| AGARD modal workflow       | `case/AGARD/Mode/`                                                           | Has its own AGENTS.md; LFS mesh decks under `Mode/mesh/`       |
| Smoke validation           | `/root/Workspace/Oldversion/scripts/smoke_fastDynamicFvMesh.sh`              | External (untracked); copies baseline zips, builds lib, runs   |
| Mesh/modal/FSI pipeline    | `/root/Workspace/Oldversion/scripts/pipeline/`                               | External; staged automation with mesh/mode consistency guard   |
| File policy                | `DEVELOPMENT.md`, `.gitignore`                                               | Source vs generated vs local-only rules                        |

## CODE MAP

LSP is not configured for OpenFOAM `.C/.H` in this harness; a codegraph index exists (`.codegraph/`) but no codegraph tool is exposed here. Grep-measured central nodes:

| Symbol/File                                             | Type     | Location                                                               | Refs                        | Role                                                           |
| ------------------------------------------------------- | -------- | ---------------------------------------------------------------------- | --------------------------- | -------------------------------------------------------------- |
| `fastDynamicFvMesh`                                     | class    | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMesh.H`              | 9 library sources           | `dynamicRefineFvMesh` subclass, modal state, AMR state, caches |
| `fastDynamicFvMesh::update()`                           | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshMain.C`          | 7 call sites / 5 exes       | Runtime mesh update orchestration; central polymorphic entry   |
| `assembleAndRelaxModalForces()`                         | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshUpdate.C`        | 1 (Main.C)                  | Force relaxation and FSI residual publication path             |
| `readControls()`                                        | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshIO.C`            | 1 (Main.C)                  | Runtime dictionary parsing and compatibility checks            |
| `readModeShapes()` / `ensureModeShapesForCurrentMesh()` | methods  | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshModeMapping.C`   | 1 each (Main.C / AMR.C)     | Mode CSV mapping, AMR remap, shared-point sync                 |
| `calcModalForces()`                                     | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshForces.C`        | 1 (Update.C)                | Pressure/shear/structural force projection                     |
| `solveStructuralDynamics()`                             | method   | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshSolver.C`        | 1 (Update.C)                | Wilson-Theta modal dynamics                                    |
| `fsiPimpleControl`                                      | class    | `fsiCoupling/fsiCoupling.H`                                            | 3 solver mains              | PIMPLE criteria plus FSI convergence                           |
| `readFsiResidual()`                                     | function | `fsiCoupling/fsiCoupling.H`                                            | 3 solver mains, 6 sites     | Solver-side force residual convergence                         |
| `myPimpleFoam.C`                                        | main     | `myPimpleFoam/`                                                        | exe `$(FOAM_USER_APPBIN)`   | Incompressible FSI solver; mesh via `controlledUpdate()`       |
| `myInterFoam.C`                                         | main     | `myInterFoam/`                                                         | exe `$(FOAM_USER_APPBIN)`   | Two-phase VOF FSI solver                                       |
| `myRhoPimpleFoam.C`                                     | main     | `myRhoPimpleFoam/`                                                     | exe `$(FOAM_USER_APPBIN)`   | Compressible FSI solver; PIMPLE body inlined 3x (no loop hdr)  |
| `interMixingFoam.C`                                     | main     | `myInterFoam/interMixingFoam/`                                         | exe `$(FOAM_APPBIN)`        | Three-phase VOF solver with local mixture models               |
| `rhoPimpleFoam.C`                                       | main     | `myRhoPimpleFoam/`                                                     | NOT built (no Make entry)   | Upstream baseline reference copy; do not treat as a target     |

Only the three main FSI solvers (`myPimpleFoam`, `myInterFoam`, `myRhoPimpleFoam`)
include `fsiCoupling.H`; no overset/SRF/three-phase variant carries FSI wiring.

## CONVENTIONS

- Use `wmake` only. Main targets: `wmake libso fastDynamicFvMesh`, `wmake myPimpleFoam`, `wmake myInterFoam`, `wmake myRhoPimpleFoam`.
- Build nested variants from their own directories with `wmake`; each owns separate `Make/files` and `Make/options` and never inherits the parent's.
- OpenFOAM environment must already be loaded. Main targets install to `$FOAM_USER_LIBBIN`/`$FOAM_USER_APPBIN`; nested variants install to `$FOAM_APPBIN` (system tree — intentional, know before you `wmake`).
- Run `ccx` and OpenFOAM tools (`blockMesh`, `checkMesh`, solvers, etc.) only inside a case directory under `/root/Workspace/...` — never in the `src/` repository root. Running them in `src/` drops `.sta/.dat/.cvg` and other artifacts into the repo and pollutes it.
- Scratch files, temporary extraction, tool output, and intermediate files must go under `/root/Workspace/tmp/` — never in the system `/tmp` (a system directory where everything written is treated as disposable).
- New dictionary keys need safe defaults (`readIfPresent`) preserving current behavior. Removing/renaming keys requires a deprecation warning path (canonical example: `refinementFaceMapTolerance` in `fastDynamicFvMeshIO.C`).
- FSI convergence is force-residual based: `fsiResidual < fsiTolerance`; opt-in OR absolute-displacement gate `fsiDispTolerance > 0` (default off). Displacement change is otherwise diagnostic-only.
- Changes under `fastDynamicFvMesh/fastDynamicFvMesh/` must update both English and Chinese docs, build the library, and run a smoke path.
- Runtime-control, file-format, or workflow changes must update example case configs and `dynamicMeshDict` guidance.
- Baseline archives live at `/root/Workspace/case_damfailure.zip` and `/root/Workspace/case_transient.zip`; copy them to temp run dirs before testing.
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
- Remote workflow (full procedure: `.opencode/skills/remote-compute-polling/`): resolve the SSH alias from `~/.ssh/config` into a variable — never bake hostnames/container ids into scripts or case files; transfer with `tar czf - -C dir . | ssh "$SSH_HOST" 'tar xzf - -C dest'` (no remote rsync) and verify `sha256sum` on both ends; launch with `setsid nohup ... > run.log 2>&1 < /dev/null &` wrapped in local `timeout` (ssh exit 124 is expected); verify in a separate ssh call; judge liveness with `ps -C <proc> -o stat= | grep -v '^Z'` — the remote PID 1 is `sleep`, so `<defunct>` launchers linger for days without meaning a job is alive; poll every 30 min on a fixed cadence.
- Record hostname, host load, rank count, and solver version for every run; performance comparisons are valid only when both machines were idle.

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
- Do not assume `scripts/` exists in-tree; it is untracked. Use `/root/Workspace/Oldversion/scripts/` or restore deliberately.

## COMMANDS

```bash
wmake libso fastDynamicFvMesh
wmake myPimpleFoam
wmake myInterFoam
wmake myRhoPimpleFoam
case/Damfailure/FDM/run_pipeline.sh --dry-run
case/Damfailure/preCICE/run_pipeline.sh --dry-run
case/AGARD/Mode/run_pipeline.sh --dry-run
/root/Workspace/Oldversion/scripts/smoke_fastDynamicFvMesh.sh --dry-run
```

## NOTES

- The external smoke script builds the library, extracts both baseline zips under `RUN_ROOT`, reads solver applications from `system/controlDict`, and logs under `RUN_ROOT/logs`.
- For AMR/restart runs, prefer reconstructing in an isolated mirror case before syncing reconstructed times back into the main case.
- The external `run_mesh_independence_automation.sh` refuses FSI if `mode/FluidNodeCoor.csv` node count differs from mesh point count unless `--allow-stale-mode` is used for diagnostics.
- No executable unit/regression test suite exists for the C++ code; validation is via builds, smoke cases, pipeline checks, and documented numerical comparisons. `case/AGARD/Mode/tests/` and `case/Damfailure/FDM/test_display_width.py` are Python-only unittest suites.
- Diagram workflow: edit `.drawio` sources first, keep final PNG/JPG only when needed, delete temporary SVG/debug exports.
