# fastDynamicFvMesh

`fastDynamicFvMesh` is a modal-superposition dynamic mesh library for OpenFOAM v2412, used for fluid-structure interaction (FSI) style mesh motion.

## What it does

At each mesh `update()`:

1. Computes pressure and viscous-shear traction on configured `fsiPatches`.
2. Projects face loads onto modal coordinates (pressure and shear tracked separately).
3. Applies the configured interface coupling acceleration (`couplingAcceleration`, default `constant` force relaxation) and updates `fsiResidual`.
4. Integrates modal dynamics with Wilson-Theta (`theta`).
5. Reconstructs mesh-point displacement and moves the mesh.
6. On `writeTime()`, writes global modal state files for restart continuity.

## Key runtime behavior

- Mode superposition count:
  - `mode/FluidNodeCoor.csv` still declares the number of **available** modes
    `N` (its first data line is `dummy, nNodes, nMode`), and
    `FluidNodeDisp1..N.csv` must all exist.
  - `nActiveModes` (optional, default `N`) superposes only the first `K` modes.
    It must satisfy `1 <= K <= N`; `K > N` aborts before any mode file is read.
  - `modalMass` and `modalDamp` must still contain `N` entries; only the first
    `K` are used. This keeps old dictionaries valid.
  - All runtime arrays, diagnostics, `StructDisp/Stress/Strain.csv`,
    `refinedPointModeShapes`, `fsiReferenceProjection` and the restart store are
    sized to `K`.
  - `K` is part of the restart contract: it is written to
    `constant/fsiRestart/nActiveModes` and a restart with a different `K`
    aborts. Legacy stores without that file are treated as `K = N`.
- Force-based coupling in mesh update path:
  - the raw modal fluid force is turned into the applied force by the configured
    `couplingAcceleration` (`constant` default, or `aitken`/`IQN-ILS`/`IQN-IMVJ`),
  - relaxed force is used for structural solve,
  - `fsiResidual` is computed from force update magnitude.
  - integrated/implicit FSI outer-loop convergence should use
    `fsiResidual < fsiTolerance` (force criterion); displacement change can be
    retained as diagnostics only.
  - the solver-side OR absolute-displacement criterion (`fsiDispTolerance` plus
    `fsiDispConsecutiveIters`, see `fsiCoupling/AGENTS.md`) sets `fsiConverged`
    when the mesh displacement stays frozen for N consecutive outer iterations;
    it is opt-in (default `0` = disabled, `fsiDispConsecutiveIters 1` = legacy
    single-step). In partitioned coupling that flag terminates the loop; in
    integrated/implicit coupling it only vetoes PIMPLE's residual-based early
    stop and cannot end the loop on its own.
- Restart support:
  - reads persisted modal fields when present (`modeState`, `modeForce`, `appliedModeDisp`, `appliedModeForce`) for non-zero-time continuation/restart only. The store is written **globally** (case root) while a relative `IOobject` instance resolves against `time().path()` -- the **processor** directory in a parallel run -- so the read is performed **on the master** from `time().globalPath()` and the (K-sized) fields are then broadcast to every rank; the legacy time-directory fallback uses the same master-read/broadcast pattern. The store is hand-written and does not honour `writePrecision`, so it is now written at 16 significant digits,
  - startup modal mapping matches the **level-0** points against the **undeformed** `constant/polyMesh/points` reference, including on a refined startup mesh where the point counts differ (`hexRef8` keeps the original points first, so base point `pI` corresponds to current level-0 point `pI`). The current mesh points have already been displaced by the modal motion (`<time>/polyMesh/points`) and can exceed `mappingTolerance` even for a short restart interval, which previously left level-0 points unmapped and aborted the run,
  - the mode shapes of AMR-generated points are persisted with the mesh as `<time>/polyMesh/refinedPointModeShapes` (one file per processor, written at full precision) and restored verbatim on restart, so the rebuilt `modeShapes_` reproduce the checkpointed run exactly instead of being re-derived by a nearest-corner search. This is required because `mapPolyMesh::pointsFromPointsMap()` is empty for `hexRef8` refinement: the in-run mapping of refinement-added points is itself a 2/4/8 corner interpolation whose corner selection cannot be reproduced from the mesh alone. Points with no recorded value still fall back to that interpolation,
  - with `meshRefinementSupport=true`, startup mapping now synchronizes mapped modal values across shared points in parallel before refinement interpolation, preventing processor-boundary geometry mismatch on continuation runs,
  - writes global modal state files at output times. At a write time the compact
    store is refreshed on every `update()` call, so the final PIMPLE/FSI
    iteration leaves the end-of-step modal state (the same state as the
    written fields); timing and modal-extension snapshots are still emitted
    once per write time.
- Performance-oriented implementation:
  - precomputed face-mode projection cache,
  - cached FSI patch IDs/references,
  - cached model pointers for turbulence/thermo/transport lookup,
  - per-mode `reduce(sumOp<scalar>())` for parallel force reduction.
- Optional diagnostics and timing:
  - `writeDiagnostics` controls `modal_diagnostics.csv` output (default `false`),
  - `trackTiming` prints fluid-vs-mesh CPU timing at `writeTime()` (default `false`),
  - `trackAmrTiming` prints one runtime-AMR timing line for each
    `dynamicRefineFvMesh::updateTopology()` call (default `false`).
- Optional prescribed modal motion (forced-oscillation / GAF probe):
  - enable by `prescribedMotion true`; the modal ODE is then bypassed and each
    mode follows `q_i(t) = A_i sin(omega t + phi_i)` with
    `omega = 2*pi*prescribedFrequency` (or the mode's own `modeFreq_` when
    `prescribedFrequency <= 0`, the default),
  - `prescribedAmplitude` (list, modal-coordinate units) and `prescribedPhase`
    (list, radians) default to all zeros, so a mode without an entry stays at
    rest; entries beyond the list length are treated as zero,
  - the modal force is still assembled every step, so the generalized
    aerodynamic force of the prescribed motion can be read from
    `modal_diagnostics.csv` (`Force_i` / `PressureForce_i` / `ShearForce_i`),
  - useful for measuring the unsteady GAF at fixed reduced frequency,
    independently of the structural model (frequencies, damping, mode count).
- Optional structural forcing:
  - enable by `structuralForceEnabled true`,
  - reads `StructNodeCoor.csv` and `StructNodeDisp1..N.csv` for structure modes,
  - reads `nStructuralForces` files from `structuralForceFilePrefix` (e.g., `StructForce1.csv`),
  - each force file provides target node coordinates and force time-series rows `(time, Fx, Fy, Fz)`,
  - maps listed coordinates to structural nodes using `structuralTargetTolerance`,
  - linearly interpolates force vectors in time and projects them to modal forces.
- Input parsing and validation:
  - CSV-like input readers accept comma or whitespace separated numeric fields,
  - empty lines and `#` comment lines are ignored for structural force files,
  - parse failures report the input file path and line number where available,
  - dictionary validation now checks `theta >= 1`, positive `mappingTolerance`
    and `maxDispChange`, positive `modalMass` entries, and non-negative
    `modalDamp` entries.
- Transformed coupled patches (e.g. rotational `cyclic` periodicity) require
  eager mesh initialisation:
  - add an empty `solvers` sub-dictionary to `constant/dynamicMeshDict`; the
    `dynamicRefineFvMesh` base then initialises the mesh (and its cyclic
    transformation tensors) before constructing its `hexRef8` engine, which
    otherwise aborts with `Faces do not seem to be correct across coupled
    boundaries`,
  - the refinement engine is still constructed, so such dictionaries must also
    provide `dynamicRefineFvMeshCoeffs { correctFluxes (); dumpLevel false;
    refineInterval 0; }`; `refineInterval 0` disables the runtime AMR path
    (and its mesh-refinement-support warning),
  - dictionaries without `solvers` keep the previous lazy initialisation.
- Optional refinement-aware FSI coupling:
  - enable by `meshRefinementSupport true`,
  - if `dynamicRefineFvMeshCoeffs` exists in `constant/dynamicMeshDict`, runtime AMR is executed in `fastDynamicFvMesh::update()` before modal-force/mesh-motion steps (same keys as OpenFOAM `dynamicRefineFvMesh`),
  - if `dynamicRefineFvMeshCoeffs` exists while `meshRefinementSupport false`,
    a warning is printed because topology changes can invalidate modal-shape and
    FSI face caches,
  - `refinementMappingDiagnostics true` prints per-topology-change mode-shape
    mapping counts for AMR validation (default `false`),
  - optional gradient-based AMR indicator:
    - `useGradIndicator true` switches AMR criteria from raw `field` values to `|grad(gradIndicatorField)|`,
    - `gradIndicatorField` defaults to `alpha.water`,
    - `gradIndicatorMagnitudeField` is optional; when set, `fastDynamicFvMesh`
      reads that registered `volScalarField` as the precomputed gradient
      magnitude and does not call `fvc::grad(gradIndicatorField)` during AMR,
    - when `gradIndicatorMagnitudeField` is not set, the internally computed
      gradient magnitude is cached for the current time index and mesh size, so
      cell refine and point unrefine selection share one gradient calculation,
    - `lowerRefineLevel/upperRefineLevel/unrefineLevel` are then interpreted on gradient magnitude and can be set with explicit hysteresis (`unrefineLevel < lowerRefineLevel`),
    - the local `interMixingFoam` three-phase interface model registers
      `alphaGradMagForAMR` from the `gradAlpha` already computed in
      `mixture.correct()`, which can be used as `gradIndicatorMagnitudeField`,
  - optional refinement-size floors can be enforced before AMR split selection:
    - `refinementMinCellVolume` (minimum estimated child-cell volume),
    - `refinementMinEdgeLength` (minimum estimated child-edge length),
    - any candidate violating either threshold is removed at refine-candidate stage, so it is excluded from both refine selection and unrefine-protection marking,
    - remaining candidates then follow the normal `dynamicRefineFvMesh` 2:1-consistency expansion,
  - optional bottom-layer AMR decision diagnostics:
    - enable by `amrLayerDiagnostics true` in `fastDynamicFvMeshCoeffs`,
    - `amrLayerDiagnosticsZ` sets the target z coordinate (default `0`),
    - `amrLayerDiagnosticsThickness > 0` selects cells within that z-distance;
      when it is `0`, the nearest local/global z layer to the target is used,
    - `amrLayerDiagnosticsNeighborDepth` adds face-neighbor cells around
      bottom-layer/candidate/selected cells,
    - `amrLayerDiagnosticsMaxRows` limits per-rank detailed cell rows without
      changing summary counts,
    - outputs `amr_layer_diagnostics_summary.csv` and
      `amr_layer_diagnostics_cells_proc<N>.csv` with each cell's AMR indicator,
      threshold state, final selection state, and reason such as
      `belowLowerThreshold`, `belowMinCellVolume`, `alreadyAtMaxRefinement`, or
      `candidateButNotSelectedByBaseAMR`,
    - base `dynamicRefineFvMesh` rejections caused by max-cells, protected-cell,
      or 2:1 consistency internals are grouped as
      `candidateButNotSelectedByBaseAMR` unless OpenFOAM base internals are also
      instrumented,
  - in integrated/implicit FSI loops, topology changes can alter mesh point count within a PIMPLE iteration; solver-side displacement-change checks must guard against `points.size()` changes before computing `gMax(mag(points - points0))`,
- startup mapping supports refined meshes with strict topology interpolation:
  - when a refined startup mesh is detected (`pointLevel>0` exists), geometric CSV mapping is restricted to level-0 points only,
  - AMR-generated points must be recovered by topology rules only (no geometric nearest-neighbour fallback),
  - recorded refined-point mode shapes (`<time>/polyMesh/refinedPointModeShapes`) are restored **first** and reproduce the checkpointed values exactly; the 2/4/8 corner rules below are the fallback for points with no recorded value,
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
  - maintenance cleanup note: removed an unused local type alias in `devRhoReff()` to eliminate compiler warning `-Wunused-local-typedefs` with no runtime behavior change.
  - maintainability note: member declarations in `fastDynamicFvMesh.H` are grouped by subsystem (core modal state, timing/diagnostics, structural loading, AMR controls, reference-face mapping) with no behavior change.
  - source-layout note: implementation is split by subsystem into `fastDynamicFvMeshMain.C`, `fastDynamicFvMeshUpdate.C`, `fastDynamicFvMeshIO.C`, `fastDynamicFvMeshModeMapping.C`, `fastDynamicFvMeshForces.C`, `fastDynamicFvMeshSolver.C`, `fastDynamicFvMeshAMR.C`, and `fastDynamicFvMeshDiagnostics.C`; this is a maintainability-only refactor with unchanged runtime keys/behavior.
  - update-orchestration note: `fastDynamicFvMesh::update()` is decomposed into named internal stages for timing, AMR topology update, modal-force relaxation, startup initialization, mesh motion, and write-time output; this keeps the call order and runtime behavior unchanged.

## 界面耦合加速器 (interface coupling accelerators)

`couplingAcceleration` selects how the raw modal fluid force `raw = modeForce_` is
turned into the applied modal force each FSI/PIMPLE outer iteration. All keys are
read with `readIfPresent` and default to the historical constant-relaxation
behaviour, so an existing dictionary behaves exactly as before (verified by the
frozen-decomposition regression; see `case/SwingBeam`).

| Key | Default | Range / tokens | Meaning |
|---|---|---|---|
| `couplingAcceleration` | `constant` | `constant` \| `aitken` \| `IQN-ILS` \| `IQN-IMVJ` | Accelerator selector. `constant` reproduces the historical `F + α(raw-F)` path exactly (including the first-iteration override). |
| `couplingRelaxation` | `1.0` | `(0, 1]` | `constant`: relaxation factor. Others: initial relaxation `α0` on the first iteration of every time window. |
| `iqnMaxHistory` | `10` | `>= 1` | Max current-window least-squares columns kept by IQN-ILS/IMVJ. |
| `iqnFilter` | `QR2` | `none` \| `QR2` \| `exponential` | Column filter. `QR2` is preCICE-faithful (relative orthogonal-complement threshold). `exponential` is an age-cutoff defined by this library (not a preCICE construct). |
| `iqnFilterLimit` | `1e-2` | `> 0`; `< 1` when `exponential` | QR2 relative threshold / exponential decay base. |
| `imvjReuse` | `5` | `>= 0` | Finished time windows whose IMVJ `(W̃, Z)` chunks are retained. `0` degenerates to per-window IQN-ILS. |
| `iqnDiagnostics` | `false` | `true`/`false` | Master-only per-iteration `CouplingAcceleration ...` line (name, active columns, `||r||`). |

- `aitken`: scalar Küttler–Wall secant `α_k = -α_{k-1}(r_{k-1}·(r_k-r_{k-1}))/||r_k-r_{k-1}||²`, clamped to `(0,1]`.
- `IQN-ILS`: current-window interface quasi-Newton with inverse least squares; `F_{k+1} = raw + W c`, `(V)c = -r`.
- `IQN-IMVJ`: adds a cross-window multi-vector inverse-Jacobian approximation with `imvjReuse` retained chunks.
- The first history-live iteration of each window always uses the initial relaxation `couplingRelaxation`; the two startup pseudo-updates (`initializeStartupStateIfNeeded`) never train history.
- Accelerator state is in-memory only (no new restart files); after a restart the accelerator re-seeds from `couplingRelaxation`.
- The published `fsiResidual_` formula is unchanged (`max_m |ΔF|/(|F_new|+1e-12)`) and remains the only published residual.

### Added-mass caveat (measured)

The accelerators are implemented and were verified against their reference formulas,
but on the strongly added-mass-dominated SwingBeam case (closed water tank, fixed-point
gain `|g| ≈ 30–130`) they diverge. Every IQN-ILS/IQN-IMVJ configuration probed —
`α0 ∈ {1e-4, 5e-4, 1e-3}`, `iqnMaxHistory ∈ {5, 10}`, `iqnFilter ∈ {QR2, none}`,
`imvjReuse ∈ {0, 5}` — hit a GAMG SIGFPE between `t=0.001` and `t=0.017 s` at
`endTime 0.05 s` (13/13 cells); Aitken diverges at `t=0.001–0.023 s`. The published
`fsiResidual` stays small mid-window while the internal fixed-point defect `||r||`
grows ≈1.3–2× per step until SIGFPE. preCICE's mature Aitken/IQN-ILS fail on the same
case, so this is a property of the added-mass problem, not of this implementation.

For such cases keep `couplingAcceleration constant` with a small `couplingRelaxation`
(SwingBeam uses `0.004`) and use the solver-side OR displacement criterion
(`fsiDispTolerance`, optionally with `fsiDispConsecutiveIters > 1` to require
consecutive quiet iterations, see `fsiCoupling/AGENTS.md`).

## 模态扩展 (modal extension)

Optional structural post-processing that recovers the structure response from
the same modal data already used to move the fluid mesh. Every switch defaults
to `false`, so an existing case behaves exactly as before: no extra files are
read, no output directory is created, and the default-off path is bit-identical.

- `modalExtensionEnabled` (master switch, default `false`): recovers structural
  displacement. Enabling it also loads `StructNodeCoor.csv` and
  `StructNodeDisp1..N.csv` even when `structuralForceEnabled` is `false`.
- `modalExtensionStressEnabled` (default `false`): additionally recovers
  structural stress. Requires `modalExtensionEnabled true`.
- `modalExtensionStrainEnabled` (default `false`): additionally recovers
  structural strain. Requires `modalExtensionEnabled true`.
- `structNodeStressPrefix` (default `StructNodeStress`): input prefix, e.g.
  `StructNodeStress1.csv`.
- `structNodeStrainPrefix` (default `StructNodeStrain`): input prefix, e.g.
  `StructNodeStrain1.csv`.

### Inputs

When `modalExtensionStressEnabled true`:

- `StructNodeStress1.csv ... StructNodeStressN.csv` (one per mode).

When `modalExtensionStrainEnabled true`:

- `StructNodeStrain1.csv ... StructNodeStrainN.csv` (one per mode).

Each file follows the same contract as `StructNodeDisp`:

- the first numeric line is the frequency, optionally followed by the node and
  mode counts (`frequency, nStructNodes, nMode`) which are validated against
  `StructNodeCoor.csv` when present,
- the body has exactly `nStructNodes` rows in the same ascending node-id order
  as `StructNodeCoor.csv`,
- stress rows carry six components `(SXX,SYY,SZZ,SXY,SYZ,SZX)`, strain rows
  carry six components `(EXX,EYY,EZZ,EXY,EYZ,EZX)`.

Committing fewer rows, or feeding a 3-column displacement file to the 6-column
reader, aborts with a fatal error naming the file and line. Inputs accept comma
or whitespace separators and `#` comment lines like the other CSV readers.

### Outputs

At each `writeTime()`, the master process writes under
`postProcessing/fastDynamicFvMesh/<time>/`:

- `StructDisp.csv` (3 columns: `ux,uy,uz`),
- `StructStress.csv` (6 columns) when stress recovery is on,
- `StructStrain.csv` (6 columns) when strain recovery is on.

Each file uses the FDM CSV contract shared with the CalculiX exporter: a header
`(time, nStructNodes, nMode)` and every data row formatted as scientific
uppercase with 10 decimals and width 20 (`%20.10E`), trailing comma included.
Rows are ordered by the `StructNodeCoor.csv` row order.

### Scope and behavior

- Displacement is superposed from the modal coordinates actually applied to the
  mesh: `u(node) = sum_m appliedModeDisp[m] * structModeShapes[m][node]`.
  Stress and strain use the same superposition over their per-mode files.
- Recovery runs on the master rank only, at write times (the same cadence as
  the restart state files).
- Output is purely derived: there are no restart fields for these quantities.
- Scope is the structural node set only. No stress or strain recovery is
  computed for the virtual-elastic fluid body.

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
    // Optional: superpose only the first K available modes (default: all N).
    // Must be 1 <= nActiveModes <= N, where N is the third value in the
    // FluidNodeCoor.csv header.  modalMass/modalDamp still list N entries.
    // nActiveModes       10;

    theta              1.4;
    fsiPatches         ("FLUID_FLUID_STRUCTURE_INTERFACE");

    mappingTolerance   4e-6;
    couplingRelaxation 0.25;    // must be in (0, 1]
    // Interface coupling accelerator; default constant (historical behaviour).
    // aitken | IQN-ILS | IQN-IMVJ are implemented but diverge on strongly
    // added-mass-dominated cases - see the accelerator section above.
    couplingAcceleration constant;
    iqnMaxHistory      10;
    iqnFilter          QR2;      // none | QR2 | exponential
    iqnFilterLimit     1e-2;
    imvjReuse          5;
    iqnDiagnostics     false;
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

    // Optional modal extension (structural displacement/stress/strain);
    // all switches default false, so these keys may be omitted entirely.
    // modalExtensionEnabled      false;
    // modalExtensionStressEnabled false;
    // modalExtensionStrainEnabled false;
    // structNodeStressPrefix     StructNodeStress;
    // structNodeStrainPrefix     StructNodeStrain;

    meshRefinementSupport   false;
    refinementInterpTolerance 1e-8;
    refinementMappingDiagnostics false;
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
    trackAmrTiming     false;
    amrLayerDiagnostics false;
    amrLayerDiagnosticsZ 0;
    amrLayerDiagnosticsThickness 0;
    amrLayerDiagnosticsNeighborDepth 1;
    amrLayerDiagnosticsMaxRows 0;
    amrLayerDiagnosticsFilePrefix amr_layer_diagnostics;
}

dynamicRefineFvMeshCoeffs
{
    refineInterval     1;
    field              alpha.water;

    // Optional: use |grad(alpha.water)| as AMR indicator
    useGradIndicator   true;
    gradIndicatorField alpha.water;
    // Optional: reuse a registered precomputed magnitude field
    // gradIndicatorMagnitudeField alphaGradMagForAMR;

    // Hysteresis on gradient indicator
    lowerRefineLevel   30;
    upperRefineLevel   500;
    unrefineLevel      8;

    nBufferLayers      3;
    maxRefinement      2;
    maxCells           180000;
}
```

## Input files (`mode/`)

Required:

- `FluidNodeCoor.csv` — its first data line is `dummy, nNodes, nMode`; the third
  value `N` is the number of **available** modes and all
  `FluidNodeDisp1..N.csv` files must exist.
- `FluidNodeDisp1.csv ... FluidNodeDispN.csv` — only the first
  `nActiveModes` are read when the optional key is set.

Optional:

- `FluidPara.csv` (legacy initial modal velocities; zero if missing)
- `StructNodeCoor.csv`, `StructNodeDisp1.csv ... StructNodeDispN.csv` (required when `structuralForceEnabled true` or `modalExtensionEnabled true`)
- `StructForce1.csv ... StructForceN.csv` (required when `structuralForceEnabled true`; `N = nStructuralForces`)
- `StructNodeStress1.csv ... StructNodeStressN.csv` (required when `modalExtensionStressEnabled true`)
- `StructNodeStrain1.csv ... StructNodeStrainN.csv` (required when `modalExtensionStrainEnabled true`)

Structural mode/stress/strain files may provide only the first `K` modes or the
full `N` modes; only the first `K` are read, and the available count must be at
least `K`.

CSV-style inputs accept commas or whitespace as separators. Structural force
files may contain blank lines and `#` comments. Fatal parse errors include the
file path and line number when available.

## Outputs

- `modal_diagnostics.csv` (only when `writeDiagnostics true`)
  - columns: `Time`, `Force_i`, `PressureForce_i`, `ShearForce_i`, `StructuralForce_i`, `Disp_i`, `Vel_i`, `Acc_i`, `AppliedDisp_i`, `StructuralScale_i`
  - the `i` groups cover the active modes (`i = 1..nActiveModes`); with the
    default (no `nActiveModes`) this is the full available set.
  - on a restart/continuation the file is opened in append mode; if it is
    already non-empty the header is not written again.
- modal-extension snapshots written by master at each `writeTime()`, only when
  the matching switch is on (see the "模态扩展 (modal extension)" section):
  - `postProcessing/fastDynamicFvMesh/<time>/StructDisp.csv` (3 columns),
  - `postProcessing/fastDynamicFvMesh/<time>/StructStress.csv` (6 columns),
  - `postProcessing/fastDynamicFvMesh/<time>/StructStrain.csv` (6 columns).
- restart modal state files written at output times:
  - `modeState`
  - `modeForce`
  - `appliedModeDisp`
  - `appliedModeForce`
  - `nActiveModes` (active mode count; prevents changing `K` across restart)

In parallel runs, restart modal state files are written by master to:

- `constant/fsiRestart/modeState`
- `constant/fsiRestart/modeForce`
- `constant/fsiRestart/appliedModeDisp`
- `constant/fsiRestart/appliedModeForce`
- `constant/fsiRestart/nActiveModes`

Legacy stores without `nActiveModes` are interpreted as having been written
with all available modes active (`K = N`). A restart with a different `K` is
rejected because the checkpoint mesh already contains the old modes'
deformation.

This keeps restart metadata out of root time directories and avoids AMR/unrefine
`reconstructPar` readUpdate conflicts. On restart the master reads these files
back from `time().globalPath()/constant/fsiRestart` and broadcasts the modal
state to all ranks -- an `IOobject` read would resolve the relative instance
against `time().path()` (= `processorN`) and always miss the store.

For AMR cases, prefer reconstructing in an isolated mirror case directory first
and only syncing reconstructed time directories back after successful completion.
This prevents partial writes in the main case when `reconstructPar` fails midway
(for example, truncated field files that later trigger reader EOF errors).

## Timing output

When `trackTiming true`, master prints at write times:

```text
FSI timing [CPU s] at t=...  step{fluid=..., mesh=..., total=...}  cumulative{...}
```

When `trackAmrTiming true`, master also prints one line for each runtime AMR
topology-update call:

```text
Runtime AMR timing at t=... [max over ranks]: CPU=... s, Clock=... s, cells ... -> ..., points ... -> ...
```

In parallel, CPU and Clock are reduced with `max` over ranks, while cell and
point counts are global sums.

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

## Per-iteration overhead reduction (v6.3.0, 2026-10-09)

Behavior-preserving optimizations in the per-FSI-iteration mesh-update path;
no dictionary key, output format, or numerical behavior changed:

- `calcModalForces()`: FSI patch face areas are now held as a lightweight
  view instead of a per-call deep copy; the kinematic-pressure density scale
  is computed once per run when no `rho` field is registered (constant
  density), and recomputed only when `rho` exists or the patch face count
  changes (AMR); the reference-face aggregation buffers under
  `meshRefinementSupport` are reused members instead of per-call heap
  allocations; the modal-force parallel reduction now packs each field into a
  single allreduce (`3` collectives instead of `3 x nMode` per iteration).
- `calcStructuralModalForces()`: prescribed structural loading depends only on
  physical time, so it is now evaluated once per time index instead of every
  FSI sub-iteration (removing two `Pstream::broadcast` calls per iteration
  when `structuralForceEnabled` is on).
- `writeDiagnostics()`: `modal_diagnostics.csv` uses a persistent stream
  (opened lazily once, flushed per row) instead of an open/close pair per
  FSI iteration.

Validation (damfailure FDM+AMR case, 8 ranks, 100 steps to `t=0.05`,
`/root/Workspace/Test`): all processor field checksums
(`U`, `p_rgh`, `alpha.water`, `k`, `omega`) and `modal_diagnostics.csv` are
bit-identical before/after; `ExecutionTime` `95.86 s -> 95.31 s` (~-0.6%) on
this small case — the savings scale with FSI iteration count and patch size.

## Common failure checks

- `Entry 'nActiveModes' is K, but ... provides N available modes`: set
  `1 <= nActiveModes <= N`, or omit it to superpose all available modes.
- `modalMass`/`modalDamp` size error: list **all N available modes**, even when
  `nActiveModes < N`; only the first K entries are used.
- `Cannot change nActiveModes across restart`: restart from `t=0`, or use a
  checkpoint written with the same `nActiveModes`.
- `Missing mode shape file ...`: `FluidNodeCoor.csv` declares N available
  modes, so `FluidNodeDisp1..N.csv` must all exist even if only the first K are
  superposed.
- `Missing required entry 'fsiPatches'`: add `fsiPatches` in `fastDynamicFvMeshCoeffs`.
- `couplingRelaxation` out of range: keep it in `(0, 1]`.
- unknown `couplingAcceleration` token: use one of `constant | aitken | IQN-ILS | IQN-IMVJ`.
- `iqnMaxHistory < 1`, `imvjReuse < 0` or `iqnFilterLimit <= 0` are fatal; `iqnFilter` must be `none | QR2 | exponential`, and `exponential` additionally needs `iqnFilterLimit < 1`.
- IQN/Aitken divergence on added-mass-dominated cases is expected (see the accelerator section); use `constant` with a small `couplingRelaxation` there.
- pressure field not found: set `pressureField` correctly.
- kinematic pressure with no density: provide `rhoRef`, `transportProperties/rho`, or `rho` field.
- zero mapped points: check coordinate frame/units and increase `mappingTolerance`.
- `Faces do not seem to be correct across coupled boundaries` with cyclic
  patches: add `solvers { }` and
  `dynamicRefineFvMeshCoeffs { correctFluxes (); dumpLevel false; refineInterval 0; }`
  to `constant/dynamicMeshDict` (eager mesh initialisation, see above).
- refinement startup interpolation failure (`...points remain unmapped...`):
  - ensure the refined mesh is generated from the same base topology as the mode files,
  - check refinement is hex-like for 2/4/8 corner classification,
  - increase `refinementInterpTolerance` moderately.
- `refinementFaceMapTolerance` warning:
  - this key is deprecated and ignored,
  - topology mapping uses `mapPolyMesh` lineage only (no geometric face-match fallback).
- `modalExtensionStressEnabled`/`modalExtensionStrainEnabled` true while
  `modalExtensionEnabled` false: fatal, both require the master switch.
- empty `structNodeStressPrefix` or `structNodeStrainPrefix` while
  `modalExtensionEnabled` true: fatal.
- missing `StructNodeStress{m}.csv` or `StructNodeStrain{m}.csv` while its
  switch is on: fatal, naming the missing file path.
- modal-extension keys present but all switches false: warning that structural
  displacement, stress and strain recovery will not be performed.

## Commit-by-commit change log

This section maps each repository commit to concrete code/document changes.

### v6.7.0 — separate available and active mode counts (2026-11-03)

- `FluidNodeCoor.csv` third value now explicitly means the number of available
  modes `N`; all `FluidNodeDisp1..N.csv` files are still required.
- Added `nActiveModes` in `fastDynamicFvMeshCoeffs` (default `N`, range
  `1 <= K <= N`) to superpose only the first `K` modes.
- Runtime arrays, diagnostics, modal-extension outputs, AMR refined-point
  records, FSI reference projections and the restart store are now K-sized;
  `modalMass`/`modalDamp` remain N-sized by contract.
- Restart stores record `nActiveModes`; legacy stores without it are treated as
  `K = N`. Changing K across restart is rejected.
- `readRefinedPointModeShapes()` and `readFsiReferenceCache()` now run after the
  active count is known, fixing the v6.6 reference-cache restore ordering bug.
- AMR restart cache restore now handles ranks that own zero `fsiPatches` faces
  (`nRef == 0` is valid), keeps the per-rank readers free of collectives, and
  sizes the aggregation buffers so `buildForceProjectionCaches()` does not
  discard and re-snapshot a restored reference set.
- The restart store is refreshed on every `update()` during a write step, so
  the checkpoint contains the end-of-step modal state rather than the
  start-of-step state. `modal_diagnostics.csv` no longer appends a second
  header when restarting into an existing non-empty file.

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

### `b98cf9a` — v1.5.1 update documents (2026-03-28)

- Documentation-focused refresh:
  - updated `fastDynamicFvMesh/README.md`,
  - updated `fastDynamicFvMesh/Doc/快速动网格程序手册.md`.
- `.gitignore` was adjusted alongside the documentation update.

### `e0b42f3` — v2.0.0 add structural force (2026-04-04)

- Added structural-force loading and projection path in `fastDynamicFvMesh`:
  - support for `StructForce*.csv` time-series inputs,
  - coordinate-to-structural-node target mapping with configurable tolerance,
  - modal structural-force accumulation and coupling into structural solve.
- Added/extended runtime controls for structural forcing (`structuralForceEnabled`, `nStructuralForces`, `structuralForceFilePrefix`, `structuralTargetTolerance`, etc.).
- Updated English/Chinese docs for new structural-force workflow and file format.

### `3a5d974` — v3.0.0 add Adaptive Mesh Refinement (2026-04-08)

- Introduced runtime AMR integration in `fastDynamicFvMesh` on top of `dynamicRefineFvMesh`:
  - supports `dynamicRefineFvMeshCoeffs` driven refine/unrefine during `update()`,
  - adds topology-aware mode-shape remapping and cache rebuild for refine/unrefine events.
- Extended build wiring (`fastDynamicFvMesh/Make/options`) and updated docs/config guidance for AMR usage.
- Updated `myInterFoam` integrated-coupling path to stay robust when topology updates change point count.

### `8bf422a` — v3.1.0 refactor restart function & FSI convergence check (2026-04-08)

- `fastDynamicFvMesh` restart persistence no longer writes modal state files into root time directories.
- Restart state is written/read via `constant/fsiRestart/`, with legacy fallback to root-time files for compatibility.
- `myInterFoam`, `myPimpleFoam`, and `myRhoPimpleFoam` unified FSI convergence decision to force residual:
  - `fsiConverged = (fsiResidual < fsiTolerance)`.
- Displacement change is retained as diagnostic output only.

### `0425961` — v3.1.1 fix AMR bugs (2026-04-10)

- Fixed AMR robustness issues in `fastDynamicFvMesh`:
  - hardened refine/unrefine topology mapping for mode-shape continuity,
  - improved parallel AMR consistency in mapping/interpolation flow to avoid unstable behavior across ranks.
- Updated English/Chinese docs to reflect AMR bug fixes and refined behavior notes.

### `01b033a` — v3.1.2 add gradient-based refinement criteria and diagnostics (2026-04-15)

- Added gradient-based AMR indicator controls:
  - `useGradIndicator`,
  - `gradIndicatorField`,
  - threshold interpretation on `|grad(field)|` for refine/unrefine selection.
- Added refinement size-floor controls and related diagnostics:
  - `refinementMinCellVolume`,
  - `refinementMinEdgeLength`,
  - runtime reporting for size-floor candidate filtering.
- Updated docs and compile database (`compile_commands.json`) with the new controls/workflow.

### `9939057` — v3.2.0 split fastDynamicFvMesh (2026-04-15)

- Refactored implementation layout by splitting the original monolithic source into:
  - `fastDynamicFvMeshMain.C`: construction, initialization, and the `update()` entry point,
  - `fastDynamicFvMeshUpdate.C`: internal `update()` stages for timing, topology update, modal-force relaxation, startup initialization, mesh motion, and write-time output,
  - `fastDynamicFvMeshIO.C`: dictionary parsing, restart state I/O, modal CSV input, and structural force input,
  - `fastDynamicFvMeshModeMapping.C`: mode-shape loading, mesh-topology interpolation, shared-point synchronization, and reference-face mapping,
  - `fastDynamicFvMeshForces.C`: fluid/structural modal-force assembly and stress helpers,
  - `fastDynamicFvMeshSolver.C`: Wilson-Theta structural modal dynamics,
  - `fastDynamicFvMeshAMR.C`: runtime AMR indicator, selection, and topology-update hooks,
  - `fastDynamicFvMeshDiagnostics.C`: diagnostics and timing reports.
- Updated build wiring for multi-translation-unit compilation:
  - `fastDynamicFvMesh/Make/files`,
  - `fastDynamicFvMesh/compile_commands.json`,
  - `lnInclude` links for new implementation units.
- Synced documentation with the split:
  - `fastDynamicFvMesh/README.md`,
  - `fastDynamicFvMesh/Doc/快速动网格程序手册.md`.
- Follow-up maintainability refactor: decomposed `fastDynamicFvMesh::update()` into named internal stages in `fastDynamicFvMeshUpdate.C` while preserving the original call order and behavior.

### `09d86c5` - add modal extension (structural displacement/stress/strain recovery) (2026-10-08)

- Consolidated all structural code (force projection, mode/stress/strain
  readers, load-summary logging) into the new translation unit
  `fastDynamicFvMeshModalExtension.C`; the move is logic-, order- and
  message-preserving.
- Added modal-extension dictionary controls, all defaulting to `false`:
  `modalExtensionEnabled` (structural displacement recovery),
  `modalExtensionStressEnabled`, `modalExtensionStrainEnabled`, plus
  `structNodeStressPrefix` (default `StructNodeStress`) and
  `structNodeStrainPrefix` (default `StructNodeStrain`).
- Added the six-component input readers for `StructNodeStress{m}.csv` /
  `StructNodeStrain{m}.csv` and the write-time snapshot writer producing
  `StructDisp.csv`, `StructStress.csv` and `StructStrain.csv` under
  `postProcessing/fastDynamicFvMesh/<time>/`.
- Structural mode files now load when `modalExtensionEnabled` is true, not only
  under `structuralForceEnabled`.
- Added the opt-in `--export-stress` flag to the CalculiX deck generators
  (`*EL FILE` + `S, E` in the `*FREQUENCY` step) and the six-column
  stress/strain CSV export in `case/AGARD/Mode/03_export_modes.py`.
- Updated English/Chinese docs, example `dynamicMeshDict` guidance and the AGENT
  knowledge files for the new switches, files, outputs and failure checks.

### v6.3.0 — reduce per-iteration overhead in mesh update path (2026-10-09)

- `fastDynamicFvMeshForces.C`: face-area view instead of per-call copy;
  constant-density pressure scale computed once (recompute only when a `rho`
  field is registered); reference-face aggregation buffers promoted to reused
  members; modal-force parallel reduction packed to 3 whole-field allreduces.
- `fastDynamicFvMeshModeMapping.C`: pressure-scale cache is invalidated when
  an AMR topology change alters the FSI patch face count (correctness guard);
  the aggregation members are sized alongside the reference-face caches.
- `fastDynamicFvMeshModalExtension.C`: structural modal force is cached per
  time index, removing recompute + two broadcasts per FSI sub-iteration.
- `fastDynamicFvMeshDiagnostics.C`: persistent diagnostics stream replaces
  per-iteration open/close.
- Hygiene: removed the stale `kappa/kappa.C` entry from a myHisa
  `surfaceInterpolationSchemes/Make/files`; marked `myRhoPimpleFoam/rhoPimpleFoam.C`
  as an unbuilt baseline reference; removed a misplaced `-g` from myHisa
  `EXE_INC` flags.
- Validation: damfailure FDM+AMR benchmark (8 ranks, 100 steps) — processor
  field checksums and `modal_diagnostics.csv` bit-identical; `ExecutionTime`
  `95.86 s -> 95.31 s`.
