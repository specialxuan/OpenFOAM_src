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
| `src/modules/hisa/pseudotimeControl/pseudotimeControl.{H,C,I.H}` | Reads `fsiCoupling`/`fsiTolerance`/`fsiDispTolerance`/`fsiDispConsecutiveIters` from `constant/dynamicMeshDict` (shared `fsiCoupling.H` contract); `criteriaSatisfied()` ANDs `readFsiResidual(mesh_) < fsiTolerance` OR the displacement gate for integrated/implicit coupling; owns the streak/reference-point state (reset at `corr_ == 1`) and exposes `fsiConverged()` |
| `applications/singleRegionSolvers/singleRegionSolver.H` | `fsiCoupling partitioned`: outer FSI loop (max `maxFsiIter`) wrapping the pseudotime sub-loop, mesh frozen inside the sub-loop, convergence `forceOK OR dispOK` with the `fsiDispConsecutiveIters` streak; integrated: per-iteration `FSI/Pseudo Iteration` log + non-convergence warning driven by the pseudotime control's combined flag when the displacement gate is on; `explicit` (default): byte-identical upstream behavior |
| `src/modules/hisa/hisaModule.{H,C}` | Warn-once validation: integrated expects `moveMeshOuterCorrectors 1`, partitioned expects `0` |
| `src/modules/hisa/Make/options`, `applications/singleRegionSolvers/hisa/Make/options` | `-I../../../../fsiCoupling` (header-only dependency, repo `src/fsiCoupling`) |

## CONVENTIONS

- FSI keys live in `constant/dynamicMeshDict` (`fsiCoupling`, `maxFsiIter`, `fsiTolerance`, `fsiDispTolerance`, `fsiDispConsecutiveIters`), same names/defaults as the `my*Foam` family. Convergence is `forceOK OR dispOK`: `forceOK = fsiResidual < fsiTolerance`; `dispOK = (fsiDispTolerance > 0) && maxDiff < fsiDispTolerance` sustained for `fsiDispConsecutiveIters` consecutive iterations (the streak resets on any violation/topology change and at the start of each time step). In integrated/implicit mode this only vetoes the pseudotime residual early-stop (`baseConverged && fsiConverged`); in partitioned mode it terminates the FSI loop directly. Keys default to `0`/`1`, so an unconfigured case is byte-identical to force-only behaviour.
- The displacement gate is a one-iteration lag, exactly like the `fsiResidual` read: at loop call k the control diffs the current mesh points against the reference captured at call k-1 (move k-1's displacement), so the gate and the force residual always see the same iteration.
- The gate state lives in `pseudotimeControl`, which is constructed once per run and reused across time steps (`corr_` resets to 0); the `corr_ <= 1` reset re-captures the reference points and clears the streak each step.
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

2026-10-12, displacement OR-gate verification, local 8 ranks (`OMP_NUM_THREADS=1`) on the AGARD 445.6 soft-wing Fluid8 FDM case (`/root/Workspace/AGARD/Test/agard4456_FDM_AMR_fluid8/case`, run as `myHisa`, AMR off, 505k cells, 4 steps, identical scotch partitions reused for A/B):
- **explicit**: 0 FSI lines, 0 warnings (zero regression); the loop converges at `nPseudoCorrMin` = 20 per step.
- **integrated, gate off**: force residual plateaus at ~5e-3..1.8e-2 and never reaches `fsiTolerance 1e-4`, so the loop runs 20/50/50/50 = 170 pseudotime iterations and warns once per stuck step (4 warnings).
- **integrated + `fsiDispTolerance 1e-9` + `fsiDispConsecutiveIters 2`**: the first mesh move of each step is ~1.5e-7..5.7e-7 while the plateau after 2–3 iterations is 1e-11..4e-10; the gate latches, the loop converges at 20 iterations/step (170 → 80 total), 0 warnings, and the end-of-run forces are bit-identical to the explicit run.
- **partitioned + same gate** (`moveMeshOuterCorrectors 0`): at the first step — where both variants start from the byte-identical decomposition — `fsiDispConsecutiveIters 1` spends 1 outer FSI iteration and `2` spends 2, because the displacement change is already 0 (< 1e-9) so `N=1` accepts the first frozen iteration while `N=2` waits for a second: the streak counter is demonstrated directly. Later steps are not A/B-comparable (the step-1 early exit changes the trajectory; totals 12 vs 17 FSI iterations, 352 s vs 470 s). 0 FSI warnings, 0 `moveMeshOuterCorrectors` warnings.
- Pitfall: running `fsiCoupling partitioned` on a case configured with `moveMeshOuterCorrectors 1` moves the mesh inside the frozen sub-loop (warned once, results skewed) — set it to `0`.
