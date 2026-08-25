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
| Displacement logging | `fsiMaxDisplacementChange()`, `logFsiIteration()` | Diagnostic output only |

## CONVENTIONS

- No standalone build target; consumers include this header through each solver `Make/options` with `-I../fsiCoupling`.
- Supported `fsiCoupling` values are `explicit`, `partitioned`, `integrated`, and `implicit`.
- Defaults: `fsiCoupling explicit`, `maxFsiIter 10`, `fsiTolerance 1e-4`.
- Integrated/implicit convergence must use `fsiResidual < fsiTolerance`.
- Displacement-change checks must tolerate topology changes by comparing point counts before subtracting point fields.
- Keep this file small and header-only unless all solver build files are updated together.

## ANTI-PATTERNS

- Do not change the `fsiResidual` object name without updating `fastDynamicFvMesh` and all solver families.
- Do not make displacement change the convergence criterion again; it is diagnostic-only.
- Do not add solver-family-specific behavior here unless all three main solvers need it.
- Do not remove `implicit` compatibility as an alias for integrated coupling.
