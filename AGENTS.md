# PROJECT KNOWLEDGE BASE

**Generated:** 2026-07-19 CST
**Commit:** 918a4a8
**Branch:** master

## OVERVIEW

OpenFOAM v2412 user-code repository containing a modal-superposition dynamic mesh library, shared FSI controls, custom incompressible/VOF/compressible solver families, and mesh/modal/FSI automation scripts. Build with OpenFOAM `wmake`; there is no repository-wide `Allwmake`, `Allclean`, or CI workflow.

## STRUCTURE

```text
src/
|-- fastDynamicFvMesh/      # dynamic mesh library, AMR/restart/force projection, bilingual docs
|-- fsiCoupling/            # header-only FSI/PIMPLE convergence contract
|-- myInterFoam/            # two-phase VOF FSI solver, overset variant, three-phase child
|-- myPimpleFoam/           # incompressible FSI solver, overset and SRF variants
|-- myRhoPimpleFoam/        # compressible FSI solver family and reference baseline
|-- scripts/                # smoke, mesh conversion, modal export, staged pipelines
|-- DEVELOPMENT.md          # source/generated/local file policy
`-- AGENTS.md               # repository-wide rules
```

## WHERE TO LOOK

| Task                      | Location                                                                     | Notes                                                          |
| ------------------------- | ---------------------------------------------------------------------------- | -------------------------------------------------------------- |
| Dynamic mesh runtime      | `fastDynamicFvMesh/fastDynamicFvMesh/`                                       | `update()`, AMR, restart, force assembly, diagnostics          |
| Mesh library docs         | `fastDynamicFvMesh/README.md`, `fastDynamicFvMesh/Doc/快速动网格程序手册.md` | Must stay synced with library code                             |
| Shared FSI controls       | `fsiCoupling/fsiCoupling.H`                                                  | `fsiCoupling`, `maxFsiIter`, `fsiTolerance`, `fsiResidual`     |
| Incompressible FSI solver | `myPimpleFoam/`                                                              | Main solver plus overset and SRF variants                      |
| VOF FSI solver            | `myInterFoam/`                                                               | Two-phase main, overset variant, nested three-phase executable |
| Compressible FSI solver   | `myRhoPimpleFoam/`                                                           | Main solver, equation headers, overset variant, baseline copy  |
| Smoke validation          | `scripts/smoke_fastDynamicFvMesh.sh`                                         | Copies baseline zips to temp dirs, builds library, runs cases  |
| Mesh/modal/FSI pipeline   | `scripts/pipeline/`                                                          | Staged automation with mesh/mode consistency guard             |
| Mesh conversion tools     | `scripts/mesh/`                                                              | OpenFOAM polyMesh/SU2/STL to CDB or case helpers               |
| Modal tools               | `scripts/modal/`                                                             | CalculiX `.inp/.dat/.frd` to fastDynamicFvMesh CSV bridge      |
| File policy               | `DEVELOPMENT.md`, `.gitignore`                                               | Source vs generated vs local-only rules                        |

## CODE MAP

LSP is not configured for OpenFOAM `.C/.H` in this harness. Codegraph/text search identifies these central nodes.

| Symbol/File                                             | Type     |                                                             Location |           Refs | Role                                                           |
| ------------------------------------------------------- | -------- | -------------------------------------------------------------------: | -------------: | -------------------------------------------------------------- |
| `fastDynamicFvMesh`                                     | class    |            `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMesh.H` |        central | `dynamicRefineFvMesh` subclass, modal state, AMR state, caches |
| `fastDynamicFvMesh::update()`                           | method   |        `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshMain.C` |        central | Runtime mesh update orchestration                              |
| `assembleAndRelaxModalForces()`                         | method   |      `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshUpdate.C` |       internal | Force relaxation and FSI residual publication path             |
| `readControls()`                                        | method   |          `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshIO.C` |       internal | Runtime dictionary parsing and compatibility checks            |
| `readModeShapes()` / `ensureModeShapesForCurrentMesh()` | methods  | `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshModeMapping.C` |       internal | Mode CSV mapping, AMR remap, shared-point sync                 |
| `calcModalForces()`                                     | method   |      `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshForces.C` |       internal | Pressure/shear/structural force projection                     |
| `solveStructuralDynamics()`                             | method   |      `fastDynamicFvMesh/fastDynamicFvMesh/fastDynamicFvMeshSolver.C` |       internal | Wilson-Theta modal dynamics                                    |
| `fsiPimpleControl`                                      | class    |                                          `fsiCoupling/fsiCoupling.H` | 3 solver mains | PIMPLE criteria plus FSI convergence                           |
| `readFsiResidual()`                                     | function |                                          `fsiCoupling/fsiCoupling.H` | 3 solver mains | Solver-side force residual convergence                         |
| `myPimpleFoam.C`                                        | main     |                                                      `myPimpleFoam/` |     executable | Incompressible FSI solver                                      |
| `myInterFoam.C`                                         | main     |                                                       `myInterFoam/` |     executable | Two-phase VOF FSI solver                                       |
| `myRhoPimpleFoam.C`                                     | main     |                                                   `myRhoPimpleFoam/` |     executable | Compressible FSI solver                                        |
| `interMixingFoam.C`                                     | main     |                                       `myInterFoam/interMixingFoam/` |     executable | Three-phase VOF solver with local mixture models               |

## CONVENTIONS

- Use `wmake` only. Main targets: `wmake libso fastDynamicFvMesh`, `wmake myPimpleFoam`, `wmake myInterFoam`, `wmake myRhoPimpleFoam`.
- Build nested variants from their own directories with `wmake`; each owns separate `Make/files` and `Make/options`.
- OpenFOAM environment must already be loaded. Outputs go to `$FOAM_USER_LIBBIN`, `$FOAM_USER_APPBIN`, or variant `FOAM_APPBIN` targets.
- Run `ccx` and OpenFOAM tools (`blockMesh`, `checkMesh`, solvers, etc.) only inside a case directory under `/root/Workspace/...` — never in the `src/` repository root. Running them in `src/` drops `.sta/.dat/.cvg` and other artifacts into the repo and pollutes it.
- New dictionary keys need safe defaults preserving current behavior. Removing/renaming keys requires a deprecation warning path.
- FSI convergence is force-residual based: `fsiResidual < fsiTolerance`. Displacement change is diagnostic-only.
- Changes under `fastDynamicFvMesh/fastDynamicFvMesh/` must update both English and Chinese docs, build the library, and run a smoke path.
- Runtime-control, file-format, or workflow changes must update example case configs and `dynamicMeshDict` guidance.
- Baseline archives live at `/root/Workspace/case_damfailure.zip` and `/root/Workspace/case_transient.zip`; copy them to temp run dirs before testing.

## ANTI-PATTERNS (THIS PROJECT)

- Do not use cmake/make or create parallel build systems.
- Do not edit baseline case archives directly.
- Do not commit `*/lnInclude/`, `*/linux64GccDPInt32Opt/`, `*/.cache/`, `__pycache__/`, `.venv_fdm/`, `.codex`, `.omo/`, `.vscode/mcp.json`, `*.cvg`, or `*.sta`.
- Do not rely on `refinementFaceMapTolerance` as live behavior; it is deprecated, warned, and ignored.
- Do not use geometric fallback for AMR-generated points; topology lineage/mapping rules are mandatory.
- Do not reuse previous time-step mesh compression flux after topology changes.
- Do not claim AMR is unconditionally stable or that refinement can never degrade mesh quality.
- Do not treat stale CalculiX outputs, stale mode files, or old FSI logs as production inputs.
- Do not sync partial `reconstructPar` output into a main case after failure.

## COMMANDS

```bash
wmake libso fastDynamicFvMesh
wmake myPimpleFoam
wmake myInterFoam
wmake myRhoPimpleFoam
scripts/smoke_fastDynamicFvMesh.sh --dry-run
scripts/smoke_fastDynamicFvMesh.sh --end-time 0.01
scripts/pipeline/run_mesh_independence_automation.sh --stage check
```

## NOTES

- `scripts/smoke_fastDynamicFvMesh.sh` builds the library, extracts both baseline zips under `RUN_ROOT`, reads solver applications from `system/controlDict`, and logs under `RUN_ROOT/logs`.
- For AMR/restart runs, prefer reconstructing in an isolated mirror case before syncing reconstructed times back into the main case.
- `scripts/pipeline/run_mesh_independence_automation.sh` refuses FSI if `mode/FluidNodeCoor.csv` node count differs from mesh point count unless `--allow-stale-mode` is used for diagnostics.
- No executable unit/regression test suite exists; validation is via builds, smoke cases, pipeline checks, and documented numerical comparisons.
- Diagram workflow: edit `.drawio` sources first, keep final PNG/JPG only when needed, delete temporary SVG/debug exports.
