# fsiCoupling Knowledge

## OVERVIEW

Header-only solver-side FSI contract shared by the three main solver families.

## WHERE TO LOOK

| Task | Location | Notes |
|---|---|---|
| FSI dictionary reads | `fsiCoupling.H` | Reads `constant/dynamicMeshDict` with `MUST_READ_IF_MODIFIED` |
| PIMPLE stop criteria | `fsiPimpleControl` | Adds `fsiConverged` to base `pimpleControl::criteriaSatisfied()` |
| Coupling mode checks | `isPartitionedFsiCoupling()`, `isIntegratedFsiCoupling()` | `implicit` aliases integrated coupling |
| Force residual read | `readFsiResidual()` | Looks up mesh object `fsiResidual`; fallback value is non-converged |
| Absolute displacement criterion | `fsiCouplingControls.dispTolerance`, top-level `fsiDispTolerance` | Opt-in OR gate: `maxDiff < dispTolerance` also terminates; default `0` = disabled |
| Displacement logging | `fsiMaxDisplacementChange()`, `logFsiIteration()` | Per-outer-iteration global max `|points-points0|`; used by the OR gate and as diagnostics |

## CONVENTIONS

- No standalone build target; consumers include this header through each solver `Make/options` with `-I../fsiCoupling`.
- Supported `fsiCoupling` values are `explicit`, `partitioned`, `integrated`, and `implicit`.
- Defaults: `fsiCoupling explicit`, `maxFsiIter 10`, `fsiTolerance 1e-4`, `fsiDispTolerance 0` (disabled).
- Integrated/implicit convergence is force-based: `fsiResidual < fsiTolerance`. When top-level `fsiDispTolerance > 0`, the solver ORs an absolute displacement criterion `maxDiff < fsiDispTolerance` (`maxDiff` = per-outer-iteration global max `|points-points0|`), so the loop also stops once the mesh displacement freezes. The OR gate suppresses zero-crossing false warnings of the relative force residual; it is solver-side and default-off.
- Displacement-change checks must tolerate topology changes by comparing point counts before subtracting point fields.
- Keep this file small and header-only unless all solver build files are updated together.

## ANTI-PATTERNS

- Do not change the `fsiResidual` object name without updating `fastDynamicFvMesh` and all solver families.
- Do not make the absolute displacement criterion unconditional; it is opt-in via `fsiDispTolerance > 0` and must stay default-off for bit-identical old behaviour. The relative force residual remains the primary gate.
- Do not add solver-family-specific behavior here unless all three main solvers need it.
- Do not remove `implicit` compatibility as an alias for integrated coupling.
