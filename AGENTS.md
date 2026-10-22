# PROJECT KNOWLEDGE BASE

**Generated:** 2026-10-02 CST
**Updated:** 2026-10-11 CST (HEAD `e4e729d`)
**Branch:** master

## OVERVIEW

OpenFOAM v2412 user-code repository containing a modal-superposition dynamic mesh library, shared FSI controls, custom incompressible/VOF/compressible solver families, and tracked case-automation pipelines (Damfailure, AGARD 445.6, SwingBeam). Build with OpenFOAM `wmake`; there is no repository-wide `Allwmake`, `Allclean`, or CI workflow.

## RESPONSE STYLE

- Default to brief answers: give the conclusion plus the key action/paths only — no extended data dumps, citation lists, or tangents.
- Switch to rich detail (data, references, alternatives, background) only when the user explicitly asks for a detailed answer.

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
|-- README.md               # project overview, build/validation, file policy
`-- AGENTS.md               # repository-wide rules
```

`scripts/` was untracked in `effbecd` (v3.3.9) and is gitignored; the live copy now
exists only outside git at `/root/Workspace/_Oldversion_/scripts/` (mesh/modal
converters, staged pipelines).

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
| SwingBeam case pipelines   | `case/SwingBeam/FDM/`, `case/SwingBeam/preCICE/`                             | Closed tank + clamped cantilever beam; six-stage Python stages |
| Smoke / regression path    | `case/Damfailure/FDM/run_pipeline.sh`                                        | In-repo and self-contained (own mesh + modal stages); `--dry-run` then a short `--steps` run |
| Mesh/modal/FSI pipeline    | `/root/Workspace/_Oldversion_/scripts/pipeline/`                             | External; staged automation with mesh/mode consistency guard   |
| File policy                | `AGENTS.md`, `.gitignore`                                                    | Source vs generated vs local-only rules                        |

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
- On local, run `ccx` and OpenFOAM tools (`blockMesh`, `checkMesh`, solvers, etc.) only inside a case directory under `/root/Workspace/...` — never in the `src/` repository root. On remote hosts, run them only under `$HOME/Special/...`. Running them in `src/` drops `.sta/.dat/.cvg` and other artifacts into the repo and pollutes it.
- Scratch files, temporary extraction, tool output, and intermediate files must go under `/root/Workspace/tmp/` — never in the system `/tmp`.
- Do not install new Python libraries into the system interpreter (no system-wide `pip install`, no `--user`, no distro packages). Install them into the repo-root venv `.venv_fdm/` instead (`.venv_fdm/bin/python -m pip install <package>`) and run Python scripts with `.venv_fdm/bin/python`.
- Smoke/regression validation is **in-repo and self-contained**: `case/Damfailure/FDM/run_pipeline.sh` generates its own mesh and modal stages (no external baseline archive is required) — preview with `--dry-run`, then run a short `--steps` subset.
- Pipeline drivers (`case/*/*/run_pipeline.sh`) support `--dry-run`, `--steps` for partial runs, and `--skip-if-fresh` (mtime gate on generation stages only; solver/post/video stages always run). Default case roots point under `/root/Workspace/...`, not in-tree.
- Video export stages (`06_export_*_video.py`) drive ParaView (`paraview.simple` — needs `pvbatch`/ParaView Python); other pipeline stages are stdlib-only Python.
- New dictionary keys need safe defaults (`readIfPresent`) preserving current behavior. Removing/renaming keys requires a deprecation warning path (canonical example: `refinementFaceMapTolerance` in `fastDynamicFvMeshIO.C`).
- FSI convergence is force-residual based: `fsiResidual < fsiTolerance`; opt-in OR absolute-displacement gate `fsiDispTolerance > 0` (default off). When enabled, the gate must hold for `fsiDispConsecutiveIters` consecutive outer iterations (default `1` = legacy single-step) to reject one-step stagnation false positives. It terminates the partitioned loop directly; in integrated/implicit mode it acts only through the combined PIMPLE/pseudotime AND-veto (so it can stop the loop at the residual criterion—shorter than waiting for force convergence—but cannot terminate the loop on its own). Displacement change is otherwise diagnostic-only. The three `my*Foam` mains wire the gate through `fsiPimpleControl`; myHisa implements the same AND veto inside `pseudotimeControl::criteriaSatisfied()`.
- Changes under `fastDynamicFvMesh/fastDynamicFvMesh/` must update both English and Chinese docs, build the library, and run a smoke path (`case/Damfailure/FDM/run_pipeline.sh --dry-run`, then a short `--steps` FSI run) or document the exact blocker.
- Runtime-control, file-format, or workflow changes must update example case configs and `dynamicMeshDict` guidance.
- `case/AGARD/**/*.zip` is Git LFS (`.gitattributes`); never re-add those archives as plain blobs.

## MACHINES & RUN PLACEMENT

|                       | Local (dev container)                         | `FSIDocker`                                      | `SuperComputer` (Slurm HPC)                          |
| --------------------- | --------------------------------------------- | ------------------------------------------------ | ---------------------------------------------------- |
| Form                  | WSL2 container                                | Docker container on WSL2                         | rented Slurm cluster                                 |
| Access                | local shell                                   | SSH alias `FSIDocker`                            | SSH alias `SuperComputer`                            |
| User / HOME           | `root` / `/root`                              | `root` / `/root`                                 | `xuzili` / `/data/home/xuzili`                       |
| CPU / RAM             | i7-14700K (8P+12E); 20 logical; 31 GB         | 64 cores; 251 GB                                 | `intel36` partition; 36 cores/node; ~93 GB/node      |
| Scheduler             | none                                          | none                                             | Slurm                                                |
| MPI                   | OpenMPI; `mpirun`                              | OpenMPI 4.1.6; `mpirun --allow-run-as-root`      | system OpenMPI 4.0.3; `srun --mpi=pmix_v3`           |
| OpenFOAM              | v2412; hisa 1.12.7 + my*Foam built            | v2412; `myHisa`/`myPimpleFoam`/`myInterFoam`     | v2412 + FDM + `my*Foam` + `myHisa` + preCICE         |
| Environment           | repo-local                                    | `source /usr/lib/openfoam/openfoam2412/etc/bashrc` | `source $HOME/opt/openfoam-precice.env`            |
| Remote write scope    | local repo + `/root/Workspace`                | **only `$HOME/Special`**                         | **only `$HOME/Special`**                             |
| Source (read-only)    | repo                                          | `/root/OpenFOAM/user-v2412/src`                  | `/data/home/xuzili/Special/FDM/src`                  |
| Tools                 | full (taskset available)                      | no rsync/tmux/screen; git/unzip/python3          | Slurm, `srun`, `sbatch`, `sacct`; no rsync           |

- Local parallel runs: **exactly 8 MPI ranks** (= the 8 P-cores). More ranks spill onto E-cores and run slower than 8. WSL2 does not expose P/E topology (`core_type` absent), so do not blind-pin with `taskset`; use `numberOfSubdomains 8;` and `export OMP_NUM_THREADS=1`, and run only one parallel case at a time.
- Send to remote when ANY holds: > ~3M cells, estimated memory > 20 GB, estimated wall time > 12 h, or ranks > 8. The AGARD combined CalculiX modal decks (~37 GB peak RSS) are always remote-class.
- When unsure: run 100–500 steps locally on a coarse mesh / short `endTime`, extrapolate `ExecutionTime`, then decide.
- **Remote write scope is strictly `$HOME/Special`.** On `SuperComputer`, every other `$HOME` directory belongs to group members: do not create, overwrite, or delete there. On `FSIDocker`, treat `/root/Workspace` as legacy/read-only; new cases go under `/root/Special`. Existing source and environment trees outside `Special` are read-only.
- **HPC launcher rule:** on `SuperComputer` use only `srun --mpi=pmix_v3` for MPI. **Never use `mpirun`, and never fall back to it if `srun` fails.** `srun --mpi=openmpi` and `srun --mpi=pmi2` fail on this host; `MpiDefault=none`. If `srun` fails, record the exact error and stop.
- **FSIDocker launcher rule:** there is no Slurm/`srun`; use `mpirun --allow-run-as-root`. Do not copy HPC `srun` commands there.
- **Short jobs (expected <= 10 min): do not poll.** Run synchronously and report the terminal result. HPC: `sbatch --wait` or a direct `srun` in an allocation. FSIDocker: foreground `mpirun` (optionally under `timeout`). If a short job exceeds 10 min, switch to long-job monitoring.
- **Long jobs (> 10 min):** on HPC submit with `sbatch` and poll `squeue`/`sacct` every 30 min (60 min after 5 h); use scheduler state, not PID/zombie logic. On FSIDocker launch detached with `setsid nohup`, then poll every 30/60 min with zombie-aware `ps -C` and log checks. The generic procedure remains in `.opencode/skills/remote-compute-polling/`.
- Transfer with `tar czf - -C dir . | ssh "$SSH_HOST" 'tar xzf - -C dest'` or `scp`; verify `sha256sum` on both ends. No remote `rsync`.
- Reconstruct on a compute node (`srun -n 1 reconstructPar` or a short `sbatch`) for HPC; never on the login node. Pull only root time directories, not `processor*`. On FSIDocker, reconstruct directly in the case under `$HOME/Special`.
- Clean up test cases/temp files only under `$HOME/Special`, and only after the pulled results are verified. On HPC, require `squeue -u xuzili` to be empty first. Never delete source, environment installs, build logs, or legacy `/root/Workspace` without explicit instruction.
- **Every case run must go through the target-specific run subagent:** `openfoam-run-local`, `openfoam-run-fsidocker`, or `openfoam-run-supercomputer`. Short jobs may be run to completion without polling (foreground subagent is acceptable); long jobs must be launched in the background and polled. Do not run cases as a bare `shell` background job.
- Record hostname, host load, rank count, solver version, and launcher for every run; performance comparisons are valid only when both machines were idle.
- `decomposePar` with `method scotch` is **non-reproducible in this environment** (re-running on identical input yields different partitions, hence different rounding). For bit-for-bit A/B comparisons, decompose once and reuse the same `processor*/`, or pin a deterministic decomposition method.

## SUBAGENTS

Project subagents are registered as Markdown files under `.opencode/agents/` (V2-native; the legacy `.opencode/agent/` path is still discovered but no longer used).

| Agent | File | Delegate when | Required deliverable |
| --- | --- | --- | --- |
| `openfoam-code-change` | `.opencode/agents/openfoam-code-change.md` | One scoped patch with a known acceptance test | diff + `wmake` log + smoke evidence |
| `openfoam-analysis` | `.opencode/agents/openfoam-analysis.md` | Write/run an analysis script on specified data | script path + numbers/table + input hash |
| `openfoam-case-assembly` | `.opencode/agents/openfoam-case-assembly.md` | Assemble/verify a runnable case | case dir + dict snippets + `checkMesh` |
| `openfoam-run-local` | `.opencode/agents/openfoam-run-local.md` | Local WSL solver/build/pipeline run | terminal state + timing/residuals/exit code + launcher |
| `openfoam-run-fsidocker` | `.opencode/agents/openfoam-run-fsidocker.md` | FSIDocker solver/build/pipeline run | terminal state + timing/residuals/exit code + launcher |
| `openfoam-run-supercomputer` | `.opencode/agents/openfoam-run-supercomputer.md` | SuperComputer Slurm solver/build/pipeline run | job state + timing/residuals/exit code + launcher |

- All three run agents load the shared rules in `.opencode/skills/openfoam-run-common/SKILL.md`.
- Delegate only mechanical work with explicit pass/fail acceptance criteria. Keep in the main conversation: criteria and trade-off decisions, cross-step causal reasoning, and risky operations (deleting remote originals, changing production data, overwriting the user's uncommitted work).
- Subagents start blank. Every prompt must carry: (1) context — absolute paths, current hashes/IDs, known conclusions, prohibitions (no system Python installs; use `.venv_fdm/`); (2) acceptance criteria — concrete pass/fail, not "make it work"; (3) report format — files changed, commands run, results, what was not verified, no long prose.
- Prefer one subagent per acceptance test: split "patch", "run", and "analyse" rather than asking one agent for all three.
- The agent definitions enforce this too: if the prompt lacks context or acceptance criteria, the subagent must stop and ask instead of guessing.

## FILE POLICY

- Source: `.C/.H` and `Make/` under `fastDynamicFvMesh/fastDynamicFvMesh/`, `myInterFoam/`, `myPimpleFoam/`, `myRhoPimpleFoam/` (and `myHisa/` per its own `AGENTS.md`).
- Tracked docs: `README.md`, `AGENTS.md`, `fastDynamicFvMesh/README.md`, `fastDynamicFvMesh/Doc/快速动网格程序手册.md`; `Plan.md`-style drafts are local and gitignored.
- Never commit build products (`Make/linux64GccDPInt32Opt/`, `.cache/`, `lnInclude/`), local tool config (`.vscode/mcp.json`, `.venv_fdm/`, `__pycache__/`, `.codex`), or generated artifacts (`*.cvg`, `*.sta`, ...).
- Python packages install only into `.venv_fdm/`; pipeline stages are stdlib-only except video export, which needs ParaView's `pvbatch`.

## ANTI-PATTERNS (THIS PROJECT)

- Do not use cmake/make or create parallel build systems.
- Do not edit baseline case archives directly.
- Do not commit `*/lnInclude/`, `*/linux64GccDPInt32Opt/`, `*/.cache/`, `__pycache__/`, `.venv_fdm/`, `.codex`, `.vscode/mcp.json`, `*.cvg`, or `*.sta`.
- Do not rely on `refinementFaceMapTolerance` as live behavior; it is deprecated, warned, and ignored.
- Do not use geometric fallback for AMR-generated points; topology lineage/mapping rules are mandatory (code-enforced `FatalError` in `fastDynamicFvMeshModeMapping.C`).
- Do not reuse previous time-step mesh compression flux after topology changes.
- Do not claim AMR is unconditionally stable or that refinement can never degrade mesh quality.
- Do not treat stale CalculiX outputs, stale mode files, or old FSI logs as production inputs.
- Do not sync partial `reconstructPar` output into a main case after failure.
- Do not assume `scripts/` exists in-tree; it is untracked. Use `/root/Workspace/_Oldversion_/scripts/` or restore deliberately.
- Do not use `mpirun` on `SuperComputer`; do not use `srun` on `FSIDocker`.
- Do not write outside `$HOME/Special` on remote shared hosts; other `$HOME` directories belong to group members.
- Do not poll short jobs (expected <= 10 min); run them synchronously and report the terminal result.
- Do not delete remote test data before the pulled results have been verified, and never delete source/environment/install trees as part of test cleanup.

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
case/Damfailure/FDM/run_pipeline.sh --dry-run   # in-repo smoke/regression path (self-contained)
```

## NOTES

- The smoke/regression path is the in-repo `case/Damfailure/FDM/run_pipeline.sh`: it generates its own mesh and modal stages and then runs the FDM solver, so no external baseline archive is involved. Use `--dry-run` for preflight and `--steps` for a short subset; logs land under the run root.
- For AMR/restart runs, prefer reconstructing in an isolated mirror case before syncing reconstructed times back into the main case.
- The external `run_mesh_independence_automation.sh` (`_Oldversion_/scripts/pipeline/`) refuses FSI if `mode/FluidNodeCoor.csv` node count differs from mesh point count unless `--allow-stale-mode` is used for diagnostics.
- No executable unit/regression test suite exists for the C++ code; validation is via builds, smoke cases, pipeline checks, and documented numerical comparisons. `case/AGARD/Mode/tests/`, `case/Damfailure/FDM/test_display_width.py`, and `case/SwingBeam/FDM/test_display_width.py` are Python-only unittest suites.
