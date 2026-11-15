# fastDynamicFvMesh Knowledge

## OVERVIEW

Modal-superposition dynamic mesh library for OpenFOAM v2412; highest-blast-radius code in this repo. The AMR-free `FdmCore` and the FDM-AMR `fastDynamicFvMesh` facade are separate build units.

## STRUCTURE

The library is split into an AMR-free FDM core and an FDM-AMR facade. The
runtime OpenFOAM type name remains `fastDynamicFvMesh`.

| Path | Purpose |
|---|---|
| `fastDynamicFvMesh/fdmCore.H` | `FdmCore` declaration: modal/FSI state, force projection, structural dynamics, restart; holds a mesh reference, inherits no AMR class |
| `fdmCoreMain.C` | `FdmCore` constructor, type-name helper, `initialiseCoreState()`, `updateFdm()` |
| `fastDynamicFvMesh.H` | `fastDynamicFvMesh` facade: inherits `dynamicRefineFvMesh` (AMR) and `FdmCore`; keeps the historical type name |
| `fastDynamicFvMeshFacade.C` | Runtime selection, facade constructor/`init()`, top-level `update()`, `writeObject()` |
| `fastDynamicFvMeshUpdate.C` | Core update stages, FSI residual publication, timing/output |
| `fastDynamicFvMeshIO.C` | Core dictionary parsing, modal restart I/O, CSV/modal/structural-force reads |
| `fastDynamicFvMeshModeMapping.C` | Core mode mapping and topology interpolation helpers |
| `fastDynamicFvMeshForces.C` | Pressure/shear/structural modal-force assembly |
| `fastDynamicFvMeshModalExtension.C` | Modal extension: structural force projection plus displacement/stress/strain recovery and snapshot output |
| `fastDynamicFvMeshSolver.C` | Wilson-Theta structural dynamics |
| `fastDynamicFvMeshAMR.C` | Facade AMR: indicator cache, refine selection, topology hooks, AMR restart side files (`refinedPointModeShapes`, `fsiReference*`), AMR mode remap |
| `fastDynamicFvMeshDiagnostics.C` | Diagnostics CSV and timing reports |
| `fastDynamicFvMesh/MakeCore/files`, `MakeCore/options` | `libfdmCore` build (no AMR dependency) |
| `Make/files`, `Make/options` | `libfastDynamicFvMesh` facade build; links `-lfdmCore` |
| `Allwmake` | Builds core first, then facade |
| `README.md`, `Doc/快速动网格程序手册.md` | Required synced docs |

## WHERE TO LOOK

| Task | Location | Notes |
|---|---|---|
| Runtime keys/defaults | `fastDynamicFvMeshIO.C`, `README.md` | Preserve old dictionaries |
| Active vs available modes | `fastDynamicFvMeshIO.C` (`readModeFiles`), `fdmCoreMain.C` | `nMode_` = active K; `nAvailableModes_` = `FluidNodeCoor.csv` N; `nActiveModes` key defaults to N; `modalMass/modalDamp` still N |
| Startup/restart state | `fdmCoreMain.C`, `fastDynamicFvMeshIO.C`, `fastDynamicFvMeshAMR.C`, `constant/fsiRestart` docs | Fresh `t=0` must not load stale restart state |
| AMR remapping | `fastDynamicFvMeshAMR.C`, core interpolation in `fastDynamicFvMeshModeMapping.C` | Topology lineage only; no geometric fallback for AMR points |
| Force projection | `fastDynamicFvMeshForces.C`, `buildForceProjectionCaches()` | Reference-face aggregation matters under refine/unrefine |
| Modal-extension controls/outputs | `fastDynamicFvMeshIO.C`, `fastDynamicFvMeshModalExtension.C`, `README.md` | `modalExtensionEnabled`/`StressEnabled`/`StrainEnabled` (all default false); `StructNodeStress`/`Strain` inputs; `postProcessing/fastDynamicFvMesh/<time>/` snapshots; structural nodes only |
| FSI convergence output | `fastDynamicFvMeshUpdate.C` | Publishes `fsiResidual` read by solvers |
| Timing/diagnostics | `fastDynamicFvMeshDiagnostics.C` | `writeDiagnostics`, `trackTiming`, `trackAmrTiming` |

## CONVENTIONS

- Build with `(cd fastDynamicFvMesh && ./Allwmake)`; it builds `libfdmCore.so` first, then `libfastDynamicFvMesh.so`.
- Keep `MakeCore/files` (core) and `Make/files` (facade) aligned with their translation units; core must not include AMR headers or reference `dynamicRefineFvMesh`/`hexRef8` symbols.
- Add runtime controls with safe defaults; validate dictionary values at read time.
- Keep deprecated keys readable with a warning path when external dictionaries may contain them.
- `fsiResidual` is force-based. Displacement change is diagnostic, not the integrated/implicit FSI convergence criterion.
- Parallel AMR/restart code must keep synchronization counts globally consistent across ranks.

## CHANGE CHECKLIST

- Update `README.md` and `Doc/快速动网格程序手册.md` for any `fastDynamicFvMesh/fastDynamicFvMesh/` code change.
- If controls/file formats/workflow changed, update example `dynamicMeshDict` guidance too.
- Run `(cd fastDynamicFvMesh && ./Allwmake)`.
- Check `nm -C $FOAM_USER_LIBBIN/libfdmCore.so | grep -E "dynamicRefineFvMesh|hexRef8|meshCutter"` stays empty.
- Run the in-repo smoke path — `case/Damfailure/FDM/run_pipeline.sh --dry-run`, then a short `--steps` FSI run — or document the exact blocker.
- For AMR/restart changes, run at least the pipeline `--stage check` and a short FSI path when feasible.

## ANTI-PATTERNS

- Do not reintroduce a monolithic implementation file; the current split is intentional.
- Do not use geometric nearest-neighbour fallback for AMR-generated points.
- Do not treat AMR `maxCells` as a hard cap: v2412 base `selectRefineCells` can append all candidates at `(maxCells - N)/7`; check the first refine event's candidate budget.
- Do not load restart state on fresh `t=0` starts.
- Do not treat `refinementFaceMapTolerance` as active behavior.
- Do not write modal restart files into root time directories for parallel runs.
- Do not change runtime file names or dictionary keys without compatibility handling.
