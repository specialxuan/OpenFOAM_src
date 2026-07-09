# fastDynamicFvMesh Knowledge

## OVERVIEW

Modal-superposition dynamic mesh library for OpenFOAM v2412; highest-blast-radius code in this repo.

## STRUCTURE

| Path | Purpose |
|---|---|
| `fastDynamicFvMesh/fastDynamicFvMesh.H` | Class declaration, state, controls, caches |
| `fastDynamicFvMeshMain.C` | Runtime selection, constructor/init, top-level `update()` |
| `fastDynamicFvMeshUpdate.C` | Update stages, FSI residual publication, timing/output |
| `fastDynamicFvMeshIO.C` | Dictionary parsing, restart I/O, CSV/modal/structural-force reads |
| `fastDynamicFvMeshModeMapping.C` | Mode mapping, AMR topology interpolation, reference-face mapping |
| `fastDynamicFvMeshForces.C` | Pressure/shear/structural modal-force assembly |
| `fastDynamicFvMeshSolver.C` | Wilson-Theta structural dynamics |
| `fastDynamicFvMeshAMR.C` | AMR indicator cache, refine selection, topology hooks |
| `fastDynamicFvMeshDiagnostics.C` | Diagnostics CSV and timing reports |
| `Make/files`, `Make/options` | Library build wiring and OpenFOAM dependencies |
| `README.md`, `Doc/快速动网格程序手册.md` | Required synced docs |

## WHERE TO LOOK

| Task | Location | Notes |
|---|---|---|
| Runtime keys/defaults | `fastDynamicFvMeshIO.C`, `README.md` | Preserve old dictionaries |
| Startup/restart state | `fastDynamicFvMeshIO.C`, `constant/fsiRestart` docs | Fresh `t=0` must not load stale restart state |
| AMR remapping | `fastDynamicFvMeshModeMapping.C`, `fastDynamicFvMeshAMR.C` | Topology lineage only; no geometric fallback for AMR points |
| Force projection | `fastDynamicFvMeshForces.C`, `buildForceProjectionCaches()` | Reference-face aggregation matters under refine/unrefine |
| FSI convergence output | `fastDynamicFvMeshUpdate.C` | Publishes `fsiResidual` read by solvers |
| Timing/diagnostics | `fastDynamicFvMeshDiagnostics.C` | `writeDiagnostics`, `trackTiming`, `trackAmrTiming` |

## CONVENTIONS

- Build from repo root with `wmake libso fastDynamicFvMesh`.
- Keep `Make/files` aligned with every implementation translation unit.
- Add runtime controls with safe defaults; validate dictionary values at read time.
- Keep deprecated keys readable with a warning path when external dictionaries may contain them.
- `fsiResidual` is force-based. Displacement change is diagnostic, not the integrated/implicit FSI convergence criterion.
- Parallel AMR/restart code must keep synchronization counts globally consistent across ranks.

## CHANGE CHECKLIST

- Update `README.md` and `Doc/快速动网格程序手册.md` for any `fastDynamicFvMesh/fastDynamicFvMesh/` code change.
- If controls/file formats/workflow changed, update example `dynamicMeshDict` guidance too.
- Run `wmake libso fastDynamicFvMesh`.
- Run `scripts/smoke_fastDynamicFvMesh.sh` or document the exact blocker.
- For AMR/restart changes, run at least the pipeline `--stage check` and a short FSI path when feasible.

## ANTI-PATTERNS

- Do not reintroduce a monolithic implementation file; the current split is intentional.
- Do not use geometric nearest-neighbour fallback for AMR-generated points.
- Do not load restart state on fresh `t=0` starts.
- Do not treat `refinementFaceMapTolerance` as active behavior.
- Do not write modal restart files into root time directories for parallel runs.
- Do not change runtime file names or dictionary keys without compatibility handling.
