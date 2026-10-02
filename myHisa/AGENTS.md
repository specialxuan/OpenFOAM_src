# myHisa Knowledge

## OVERVIEW

Local copy of the HiSA solver (v1.12.7, density-based coupled implicit compressible solver) with FSI coupling added to its dual-time pseudotime loop, modeled on `myPimpleFoam`. The upstream tree at `/root/OpenFOAM/user-v2412/hisa/` is NEVER modified.

## RENAMES (vs upstream)

- Directory `hisa/` → `myHisa/`; executable `hisa` → `myHisa`; `multiRegionSolver` → `myMultiRegionSolver`; utilities carry a `my` prefix.
- Every library is renamed: `libhisaModule` → `libmyHisaModule`, `libsolverModule` → `libmySolverModule`, `libmodal` → `libmyModal`, etc. Cases loading them must use the `my` names (e.g. `libs ("libmyModal.so")`).
- Upstream-installed `libhisa*`/hisa binary remain intact side by side.

## FSI MODIFICATION (the only functional change vs upstream)

| File | Change |
|---|---|
| `src/modules/hisa/pseudotimeControl/pseudotimeControl.{H,C,I.H}` | Reads `fsiCoupling`/`fsiTolerance` from `constant/dynamicMeshDict` (shared `fsiCoupling.H` contract); `criteriaSatisfied()` ANDs `readFsiResidual(mesh_) < fsiTolerance` for integrated/implicit coupling |
| `applications/singleRegionSolvers/singleRegionSolver.H` | `fsiCoupling partitioned`: outer FSI loop (max `maxFsiIter`) wrapping the pseudotime sub-loop, mesh frozen inside the sub-loop; integrated: per-iteration `FSI/Pseudo Iteration` force-residual log + non-convergence warning; `explicit` (default): byte-identical upstream behavior |
| `src/modules/hisa/hisaModule.{H,C}` | Warn-once validation: integrated expects `moveMeshOuterCorrectors 1`, partitioned expects `0` |
| `src/modules/hisa/Make/options`, `applications/singleRegionSolvers/hisa/Make/options` | `-I../../../../fsiCoupling` (header-only dependency, repo `src/fsiCoupling`) |

## CONVENTIONS

- FSI keys live in `constant/dynamicMeshDict` (`fsiCoupling`, `maxFsiIter`, `fsiTolerance`), same names/defaults as the `my*Foam` family; convergence is force-residual based, displacement change is diagnostic-only.
- Partitioned works because `pseudotimeControl::loop()` resets `corr_` to 0 at sub-loop completion, so the mesh-move guard (`firstIter()`) moves the mesh exactly once per FSI iteration with `moveMeshOuterCorrectors 0`.
- `fsiResidual` is published by `fastDynamicFvMesh::update()`; it lags one pseudotime iteration (same as the `myPimpleFoam` reference).
- Build the whole tree with `./Allwmake` from `src/myHisa/`; original HiSA install is untouched by it.
- `fsiCoupling` keys are still no-ops for the upstream `hisa` binary — only `myHisa` reads them.

## ANTI-PATTERNS

- Do not edit `/root/OpenFOAM/user-v2412/hisa/`; all FSI work happens in this copy.
- Do not build with library names reverted to `libhisa*`; that overwrites the upstream install and silently changes the original `hisa` binary's behavior.
- Do not set `fsiCoupling partitioned` with `moveMeshOuterCorrectors != 0`; the mesh would move inside the frozen sub-loop (warned, not fatal).
- Do not expect `fsiResidual < fsiTolerance` to be reachable when the case uses `couplingRelaxation 1`; tight tolerances may legitimately never converge and the warning path is the correct outcome.

## TESTED

2026-10-02, local 8×P-core (taskset 0-7) parallel runs on the AGARD 445.6 FDM+HiSA case (`/root/Workspace/myhisa_fsi_test/`, fetched from FSIDocker `agard_p2`): explicit = 0 FSI lines (zero regression); integrated = per-iteration residual log, step converged only when `fsiRes < 1e-4` (8.6e-5) else full `nPseudoCorr` + warning; partitioned = 10 FSI iterations/step with ~1 mesh move per FSI iteration (600 pseudo solves vs ~30 mesh updates) + `maxFsiIter` warning.
