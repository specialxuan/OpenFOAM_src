# fastDynamicFvMesh

`fastDynamicFvMesh` is a modal-superposition dynamic mesh library for OpenFOAM v2412, used for fluid-structure interaction (FSI) style mesh motion.

## What it does

At each mesh `update()`:

1. Computes pressure and viscous-shear traction on configured `fsiPatches`.
2. Projects face loads onto modal coordinates (pressure and shear tracked separately).
3. Applies force relaxation (`couplingRelaxation`) and updates `fsiResidual`.
4. Integrates modal dynamics with Wilson-Theta (`theta`).
5. Reconstructs mesh-point displacement and moves the mesh.
6. On `writeTime()`, writes global modal state files for restart continuity.

## Key runtime behavior

- Force-based coupling in mesh update path:
  - relaxed force is used for structural solve,
  - `fsiResidual` is computed from force update magnitude.
- Restart support:
  - reads persisted modal fields when present (`modeState`, `modeForce`, `appliedModeDisp`, `appliedModeForce`) for non-zero-time continuation/restart only,
  - startup modal mapping now uses current mesh points by default, and only uses `constant/polyMesh/points` when its point count matches the current mesh (prevents `latestTime`+AMR point-count mismatch on restart),
  - with `meshRefinementSupport=true`, startup mapping now synchronizes mapped modal values across shared points in parallel before refinement interpolation, preventing processor-boundary geometry mismatch on continuation runs,
  - writes global modal state files at output times.
- Performance-oriented implementation:
  - precomputed face-mode projection cache,
  - cached FSI patch IDs/references,
  - cached model pointers for turbulence/thermo/transport lookup,
  - per-mode `reduce(sumOp<scalar>())` for parallel force reduction.
- Optional diagnostics and timing:
  - `writeDiagnostics` controls `modal_diagnostics.csv` output (default `false`),
  - `trackTiming` prints fluid-vs-mesh CPU timing at `writeTime()` (default `false`).
- Optional structural forcing:
  - enable by `structuralForceEnabled true`,
  - reads `StructNodeCoor.csv` and `StructNodeDisp1..N.csv` for structure modes,
  - reads `nStructuralForces` files from `structuralForceFilePrefix` (e.g., `StructForce1.csv`),
  - each force file provides target node coordinates and force time-series rows `(time, Fx, Fy, Fz)`,
  - maps listed coordinates to structural nodes using `structuralTargetTolerance`,
  - linearly interpolates force vectors in time and projects them to modal forces.
- Optional refinement-aware FSI coupling:
  - enable by `meshRefinementSupport true`,
  - if `dynamicRefineFvMeshCoeffs` exists in `constant/dynamicMeshDict`, runtime AMR is executed in `fastDynamicFvMesh::update()` before modal-force/mesh-motion steps (same keys as OpenFOAM `dynamicRefineFvMesh`),
  - optional refinement-size floors can be enforced before AMR split selection:
    - `refinementMinCellVolume` (minimum estimated child-cell volume),
    - `refinementMinEdgeLength` (minimum estimated child-edge length),
    - any candidate violating either threshold is removed at refine-candidate stage, so it is excluded from both refine selection and unrefine-protection marking,
    - remaining candidates then follow the normal `dynamicRefineFvMesh` 2:1-consistency expansion,
  - in integrated/implicit FSI loops, topology changes can alter mesh point count within a PIMPLE iteration; solver-side displacement-change checks must guard against `points.size()` changes before computing `gMax(mag(points - points0))`,
- startup mapping supports refined meshes with strict topology interpolation:
  - when a refined startup mesh is detected (`pointLevel>0` exists), geometric CSV mapping is restricted to level-0 points only,
  - AMR-generated points must be recovered by topology rules only (no geometric nearest-neighbour fallback),
  - edge midpoint from 2 endpoint corners (distance weighted),
  - face midpoint from 4 face corners (uniform corner weights),
  - cell midpoint from 8 cell corners (uniform corner weights),
    - startup/restart interpolation runs multiple 2/4/8 passes (bounded) so points that become classifiable after neighboring points are filled can also be recovered,
  - when enabled, split FSI-face forces are aggregated back to reference (pre-refine) faces before modal projection,
  - topology-change updates no longer treat `pointMap` master-label entries as fully mapped values for refinement-added points; these points are deferred to bounded multi-pass 2/4/8 interpolation,
  - when `mapPolyMesh::pointsFromPoints` lineage is available, its source-point averaging is applied before 2/4/8 fallback,
  - parallel topology updates now enforce a globally consistent number of interpolation/synchronization passes (including mixed refine+unrefine steps) to avoid MPI collective call mismatches,
  - modal-force projection is always evaluated on initial reference FSI faces loaded at startup:
    - refine: child-face forces are aggregated to initial reference faces,
    - unrefine: merged-face forces are converted to traction and redistributed to represented initial reference faces using topology-derived area weights,
  - topology mapping uses `mapPolyMesh` connectivity lineage (not geometric nearest-face fallback) for force redistribution.
  - compatibility note: `refinementFaceMapTolerance` is accepted only for backward-compatible dictionaries, validated as `>=0`, emits a deprecation warning, and is otherwise ignored.
  - internal cleanup note: AMR reference-face bookkeeping now keeps only actively used data (reference-face areas, projections, and topology weights); removed unused reference centre/normal cache with no behavior change.

## Build

From repository root:

```bash
wmake libso fastDynamicFvMesh
```

Or inside library directory:

```bash
cd fastDynamicFvMesh
wmake libso
```

## Enable in case

In `system/controlDict`:

```foam
libs ("libfastDynamicFvMesh.so");
```

In `constant/dynamicMeshDict`:

```foam
dynamicFvMesh fastDynamicFvMesh;

fastDynamicFvMeshCoeffs
{
    theta              1.4;
    fsiPatches         ("FLUID_FLUID_STRUCTURE_INTERFACE");

    mappingTolerance   4e-6;
    couplingRelaxation 0.25;    // must be in (0, 1]
    pressureField      p;
    rhoRef             1000;    // optional; required only for kinematic p when rho is unavailable
    pRef               0;
    maxDispChange      1.0;
    structuralForceEnabled   true;
    nStructuralForces        1;
    structuralForceFilePrefix StructForce;
    structuralTargetTolerance 1e-8;
    structNodeCoorFile      StructNodeCoor.csv;
    structNodeDispPrefix    StructNodeDisp;
    meshRefinementSupport   false;
    refinementInterpTolerance 1e-8;
    refinementFaceMapTolerance 0; // deprecated compatibility key; ignored
    refinementMinCellVolume 0;    // optional; 0 disables child-volume floor
    refinementMinEdgeLength 0;    // optional; 0 disables child-edge floor

    // Optional runtime AMR controls (same as dynamicRefineFvMesh)
    // Add dynamicRefineFvMeshCoeffs in the same dynamicMeshDict to enable
    // in-run topology refinement.

    // Optional; size must equal number of modes
    modalMass          (1 1 1 1 1 1 1 1 1 1);
    modalDamp          (0 0 0 0 0 0 0 0 0 0);

    writeDiagnostics   false;
    trackTiming        true;
}
```

## Input files (`mode/`)

Required:

- `FluidNodeCoor.csv`
- `FluidNodeDisp1.csv ... FluidNodeDispN.csv`

Optional:

- `FluidPara.csv` (legacy initial modal velocities; zero if missing)
- `StructNodeCoor.csv`, `StructNodeDisp1.csv ... StructNodeDispN.csv` (required when `structuralForceEnabled true`)
- `StructForce1.csv ... StructForceN.csv` (required when `structuralForceEnabled true`; `N = nStructuralForces`)

## Outputs

- `modal_diagnostics.csv` (only when `writeDiagnostics true`)
  - columns: `Time`, `Force_i`, `PressureForce_i`, `ShearForce_i`, `StructuralForce_i`, `Disp_i`, `Vel_i`, `Acc_i`, `AppliedDisp_i`, `StructuralScale_i`
- restart modal state files written at output times:
  - `modeState`
  - `modeForce`
  - `appliedModeDisp`
  - `appliedModeForce`

In parallel runs, restart modal state files are written by master to:

- `constant/fsiRestart/modeState`
- `constant/fsiRestart/modeForce`
- `constant/fsiRestart/appliedModeDisp`
- `constant/fsiRestart/appliedModeForce`

This keeps restart metadata out of root time directories and avoids AMR/unrefine
`reconstructPar` readUpdate conflicts.

## Timing output

When `trackTiming true`, master prints at write times:

```text
FSI timing [CPU s] at t=...  step{fluid=..., mesh=..., total=...}  cumulative{...}
```

## A/B performance evaluation (2026-03-28)

Compared cases:

- no-efficiency-optimization: `~/Workspace/openfoam_damfailure/`
- with-efficiency-optimization: `~/Workspace/openfoam_damfailure2/`

Same comparison conditions (verified):

- identical `controlDict`, `dynamicMeshDict`, `fvSolution`, `decomposeParDict`;
- same solver (`myInterFoam`), same parallel size (8), same 500 steps.

Result consistency:

- At `t=0.5`, checksums are identical for key outputs:
  - `U`, `p_rgh`, `alpha.water`, `k`, `omega`, `nut`,
  - `modeState`, `modeForce`, `appliedModeDisp`, `appliedModeForce`.

Performance:

- `ExecutionTime`: `4003.59 s -> 3944.01 s` (**-59.58 s**, **+1.49% faster**)
- `FSI timing` cumulative total: `3986.90 s -> 3927.73 s` (**+1.48%**)
- fluid part: `3296.10 s -> 3244.98 s` (**+1.55%**)
- mesh part: `690.80 s -> 682.75 s` (**+1.17%**)

Conclusion:

- For this damfailure benchmark, efficiency improvements provide a small but stable speedup (~1.5%) with no observed result drift.

## Common failure checks

- `Missing required entry 'fsiPatches'`: add `fsiPatches` in `fastDynamicFvMeshCoeffs`.
- `couplingRelaxation` out of range: keep it in `(0, 1]`.
- pressure field not found: set `pressureField` correctly.
- kinematic pressure with no density: provide `rhoRef`, `transportProperties/rho`, or `rho` field.
- zero mapped points: check coordinate frame/units and increase `mappingTolerance`.
- refinement startup interpolation failure (`...points remain unmapped...`):
  - ensure the refined mesh is generated from the same base topology as the mode files,
  - check refinement is hex-like for 2/4/8 corner classification,
  - increase `refinementInterpTolerance` moderately.
- `refinementFaceMapTolerance` warning:
  - this key is deprecated and ignored,
  - topology mapping uses `mapPolyMesh` lineage only (no geometric face-match fallback).

## Commit-by-commit change log

This section maps each repository commit to concrete code/document changes.

### `3524d57` — Alpha v1.0 (2026-03-16)

- Initial repository import.
- Added core components:
  - `fastDynamicFvMesh` library (`.C/.H`, build files, README/docs),
  - `myPimpleFoam` and `myInterFoam` solver variants,
  - legacy reference code under `.oldcode/`.

### `2743529` — Alpha v1.1 update myInterFoam myPimpleFoam (2026-03-16)

- Updated solver-side FSI loop wiring in:
  - `myInterFoam/myInterFoam.C`,
  - `myPimpleFoam/myPimpleFoam.C`.
- `myInterFoam` switched to dictionary default reads for `fsiCoupling`, `maxFsiIter`, `fsiTolerance`.
- Integrated/implicit branch behavior was reshaped around explicit PIMPLE-loop control with forced mesh-motion path (`moveMeshOuterCorrectors` handling).

### `1b67352` — v1.1.1 update documents (2026-03-16)

- Documentation-only commit:
  - updated `fastDynamicFvMesh/README.md`,
  - updated `fastDynamicFvMesh/Doc/快速动网格程序手册.md`,
  - removed `fastDynamicFvMesh/SUMMARY.md`.

### `ed79e23` — v1.2.0 add myRhoPimpleFoam (2026-03-17)

- Added full `myRhoPimpleFoam` solver tree (compressible variant), including:
  - main solver (`myRhoPimpleFoam.C`),
  - equation/header set (`EEqn.H`, `UEqn.H`, `pEqn.H`, etc.),
  - build files and overset variant directory.

### `dbe5395` — v1.2.1 update fastDynamicFvMesh (2026-03-19)

- Major refactor in `fastDynamicFvMesh.C/.H`:
  - split file-reading flow by introducing `readModeFiles(...)`,
  - improved CSV parsing and input validation/reporting,
  - kept `faceDiagnosticsMode` path while removing the separate `writeFaceDiagnostics` switch dependency,
  - cleanup in `readControls`, `patchDensity`, `devRhoReff`, and force assembly path.
- Goal of this commit was robustness/readability of mode-file loading and diagnostics control handling.

### `f080be2` — v1.3.0 fix coupling relaxation bug (2026-03-25)

- Updated `fastDynamicFvMesh` and `myPimpleFoam`:
  - added modal property handling in mesh library (`modalMass` support and consistency checks with mode count),
  - removed face-level diagnostics output path from active implementation,
  - added per-time-step update counter behavior in mesh update flow,
  - integrated/implicit branch in `myPimpleFoam` gained explicit displacement-change convergence checks (`fsiTolerance`) and warning on non-convergence.

### `ca7e7da` — v1.3.1 fix coupling relaxation bugs in myInterFoam and myRhoPimpleFoam (2026-03-25)

- Ported the integrated/implicit convergence-loop fixes from `myPimpleFoam` into:
  - `myInterFoam/myInterFoam.C`,
  - `myRhoPimpleFoam/myRhoPimpleFoam.C`.
- Added per-iteration displacement-change monitoring and non-convergence warnings in these solvers.

### `fe4521c` — v1.3.2 relaxation logic changed to force-based (2026-03-26)

- Core coupling change in `fastDynamicFvMesh`:
  - switched relaxation from displacement-side to force-side,
  - introduced `appliedModeForce_` and `fsiResidual_`,
  - published `fsiResidual` into object registry as `uniformDimensionedScalarField`,
  - used relaxed force for structural solve and removed displacement-relaxation step.
- Solver updates (`myPimpleFoam`, `myInterFoam`, `myRhoPimpleFoam`):
  - integrated/implicit convergence check changed from displacement-difference threshold to `fsiResidual < fsiTolerance`.

### `40fb9b3` — v1.4.0 add restart capability (2026-03-28)

- Added restart persistence mechanics in `fastDynamicFvMesh`:
  - modal fields (`modeState`, `modeForce`, `appliedModeDisp`, `appliedModeForce`) migrated to `IOField` with `READ_IF_PRESENT/AUTO_WRITE`,
  - restart initialization path for reusing saved modal state and avoiding deformation double-application,
  - wrote global modal state files at `writeTime()` into case-root time directory (parallel restart continuity).
- Added runtime controls:
  - `maxDispChange`,
  - `writeDiagnostics` (switchable diagnostics CSV output).
- Added displacement-increment limiter and reduced noisy `Info` output.
- Solver-side integrated/implicit control changed to a custom `fsiPimpleControl` (`criteriaSatisfied() = pimpleCriteria && fsiConverged`).

### `v1.8.0` — restart I/O channel split and force-based convergence unification (2026-04-08)

- `fastDynamicFvMesh` restart persistence no longer writes modal state files into root time directories.
- Restart state is now written/read via `constant/fsiRestart/` to avoid polluting reconstruct targets in AMR+unrefine workflows.
- Legacy fallback remains: if dedicated restart files are absent, startup can still read modal state from current time directory when present.
- `myInterFoam`, `myPimpleFoam`, `myRhoPimpleFoam` FSI convergence decision is unified to force residual:
  - `fsiConverged = (fsiResidual < fsiTolerance)`.
- Displacement change can still be logged diagnostically, but is no longer used as the FSI convergence criterion.

### `0340bff` — v1.4.1 add time tracker (2026-03-28)

- Added optional runtime timing profiler in `fastDynamicFvMesh`:
  - new dictionary switch `trackTiming`,
  - CPU-time accumulation for fluid interval vs mesh update interval,
  - `writeTimingReport()` output at write times:
    - `FSI timing [CPU s] at t=... step{...} cumulative{...}`.

### `ee7826a` — v1.5.0 Efficiency improvement (2026-03-28)

- Performance-focused refactor in `fastDynamicFvMesh.C/.H`:
  - added startup caches for force assembly:
    - cached FSI patch IDs/references (`fsiPatchIDs_`, `fsiPolyPatches_`),
    - cached face-mode projection coefficients (`faceModeProjection_[patch][mode][face]`),
    - reusable pressure scaling storage (`pressureScaleCache_`).
  - added model-pointer cache (`cacheModelPointers()`) for turbulence/thermo/transport lookup in shear evaluation.
- Reduced per-iteration overhead in `calcModalForces()`:
  - removed repeated patch lookup and repeated face shape averaging in hot loop,
  - switched modal force parallel reduction to direct per-mode `reduce(sumOp<scalar>())`.
- Diagnostics I/O path optimization:
  - `writeDiagnostics()` now uses in-memory header state (`diagnosticsHeaderWritten_`) instead of file existence check each call,
  - signature changed from `writeDiagnostics() const` to non-const for header-state update.
