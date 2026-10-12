# fsiCoupling Knowledge

## OVERVIEW

Header-only solver-side FSI contract shared by the three main solver families.

## WHERE TO LOOK

| Task | Location | Notes |
|---|---|---|
| FSI dictionary reads | `fsiCoupling.H` | Reads `constant/dynamicMeshDict` with `MUST_READ_IF_MODIFIED` |
| PIMPLE stop criteria | `fsiPimpleControl` | ANDs `fsiConverged` into base `pimpleControl::criteriaSatisfied()` — an integrated-mode veto, see CONVENTIONS |
| Coupling mode checks | `isPartitionedFsiCoupling()`, `isIntegratedFsiCoupling()` | `implicit` aliases integrated coupling |
| Force residual read | `readFsiResidual()` | Looks up mesh object `fsiResidual`; fallback value is non-converged |
| Absolute displacement criterion | `fsiCouplingControls.dispTolerance` / `dispConsecutiveIters`, top-level `fsiDispTolerance` / `fsiDispConsecutiveIters` | Opt-in OR gate: sets `fsiConverged` when `maxDiff < dispTolerance` for `dispConsecutiveIters` consecutive outer iterations; defaults `0` / `1` = disabled / legacy single-step |
| Displacement logging | `fsiMaxDisplacementChange()`, `logFsiIteration()` | Per-outer-iteration global max `\|points-points0\|`; used by the OR gate and as diagnostics |

## CONVENTIONS

- No standalone build target; consumers include this header through each solver `Make/options` with `-I../fsiCoupling`.
- Supported `fsiCoupling` values are `explicit`, `partitioned`, `integrated`, and `implicit`.
- Defaults: `fsiCoupling explicit`, `maxFsiIter 10`, `fsiTolerance 1e-4`, `fsiDispTolerance 0` (disabled), `fsiDispConsecutiveIters 1` (legacy single-step gate).
- Convergence flag: `fsiConverged = forceOK || dispOK`, with `forceOK = fsiResidual < fsiTolerance` and `dispOK = (dispTolerance > 0) && fsiDisplacementChangeValid() && maxDiff < dispTolerance` sustained for `dispConsecutiveIters` consecutive outer iterations (`maxDiff` = per-outer-iteration global max `|points-points0|`; the counter resets on any non-satisfying or topology-change iteration).
- How that flag stops the loop is mode-dependent. **Partitioned** runs `while (fsiIter < maxFsiIter && !fsiConverged)`, so `fsiConverged` (force OR displacement) terminates the loop directly. **Integrated/implicit** runs `while (pimple.loop())`, where `fsiPimpleControl::criteriaSatisfied()` is `pimpleControl::criteriaSatisfied() && fsiConverged_` — an AND veto: the flag can only *prevent* PIMPLE's residual-based early stop (forcing iterations until FSI converges), never terminate the loop by itself; with an unmet `residualControl` the loop still runs to `nOuterCorrectors`. In both modes a true `fsiConverged` suppresses the post-loop "FSI force residual did not converge" warning, which is the OR gate's purpose for near-zero modal-force (zero-crossing) spikes.
- `fsiDispConsecutiveIters 1` reproduces the legacy single-step gate bit-for-bit; `2` rejects one-step stagnation/plateau false positives, observable as an extra outer iteration in partitioned mode.
- `fsiDispConsecutiveIters < 1` is clamped to `1` with a warning; values greater than `maxFsiIter` warn that the gate can never trigger.
- Displacement-change checks must tolerate topology changes by comparing point counts before subtracting point fields.
- Keep this file small and header-only unless all solver build files are updated together.

## ANTI-PATTERNS

- Do not change the `fsiResidual` object name without updating `fastDynamicFvMesh` and all solver families.
- Do not make the absolute displacement criterion unconditional; it is opt-in via `fsiDispTolerance > 0` and must stay default-off for bit-identical old behaviour. Keep `fsiDispConsecutiveIters` defaulting to `1` so the legacy single-step gate stays bit-identical. The relative force residual remains the primary gate.
- Do not add solver-family-specific behavior here unless all three main solvers need it.
- Do not remove `implicit` compatibility as an alias for integrated coupling.
