# Repository Agent Notes

## Build system (OpenFOAM wmake)

This is an OpenFOAM v2412 user-code repository. Do NOT use cmake/make.
All commands require a loaded OpenFOAM environment (the user's shell has it).

```bash
# Build the shared library (shared by all solvers):
wmake libso fastDynamicFvMesh

# Build solver executables (run from src root or from their directories):
wmake myPimpleFoam
wmake myInterFoam
wmake myRhoPimpleFoam
```

Build configuration lives in `Make/files` and `Make/options` files per component.
Build artifacts go to `$FOAM_USER_LIBBIN` (libraries) and `$FOAM_USER_APPBIN` (executables).

## Source layout

| Directory | Type | Build |
|---|---|---|
| `fastDynamicFvMesh/` | Shared library | `wmake libso fastDynamicFvMesh` |
| `myPimpleFoam/` | Solver executables | `wmake myPimpleFoam` |
| `myInterFoam/` | Solver executables | `wmake myInterFoam` |
| `myRhoPimpleFoam/` | Solver executables | `wmake myRhoPimpleFoam` |
| `fsiCoupling/` | Shared header | (none — `#include`d by all three solver families) |

Each solver directory may contain sub-variants (e.g., `overPimpleDyMFoam/`, `SRFPimpleFoam/`).
`fsiCoupling/fsiCoupling.H` provides the shared `fsiPimpleControl` (extends PIMPLE convergence) and `fsiCouplingControls` struct.

## Baseline case policy

Two baseline case archives exist at `/root/Workspace/`:
- `case_damfailure.zip` — uses `myInterFoam`
- `case_transient.zip` — uses `myPimpleFoam`

For any testing/analysis: copy the archive to a temporary directory, do NOT modify
baselines directly. Use `scripts/smoke_fastDynamicFvMesh.sh` as the standard smoke test runner.

## fastDynamicFvMesh change policy

When modifying ANY file under `fastDynamicFvMesh/fastDynamicFvMesh/`, you MUST:

1. Update `fastDynamicFvMesh/README.md` (English)
2. Update `fastDynamicFvMesh/Doc/快速动网格程序手册.md` (Chinese)
3. Run `wmake libso fastDynamicFvMesh` and confirm it compiles
4. Run at least one smoke test exercising the changed path

If a change modifies runtime controls, file formats, or run workflow, also update
example case configs and `dynamicMeshDict` accordingly.

## Backward compatibility

- New dictionary keys must have safe defaults that preserve existing behavior.
- Do not remove or rename existing dictionary keys without a deprecation warning path.

## Key reference files

- `DEVELOPMENT.md` — file categorization (source vs. generated vs. local), draw.io conventions (tracked in git)
- `Plan.md` — development roadmap and priorities (gitignored, local-only)
- `Plan_1-5_Finished.md` — record of completed priorities 1-5 (gitignored, local-only)
- `AGARD4456_FSI_work_log.md`, `AGARD4456_M0499_FSI_steps.md` — FSI case notes (gitignored, local-only)

## Do NOT commit

- `*/lnInclude/`, `*/linux64GccDPInt32Opt/`, `*/.cache/` (OpenFOAM build artifacts)
- `__pycache__/`, `*.cvg`, `*.sta`
- Local config: `.vscode/mcp.json`, `.codex`, `.venv_fdm/`
