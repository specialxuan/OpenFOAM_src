# Fast Dynamic Mesh for OpenFOAM

## 1. Project Introduction

This library implements a fast dynamic mesh method based on modal superposition, ported from an ANSYS Fluent UDF. It allows for efficient Fluid-Structure Interaction (FSI) simulations by deforming the fluid mesh according to pre-calculated structural modes.

### Key Features

- **Modal Superposition**: Deforms mesh based on mode shapes and generalized coordinates.
- **Two-Way Coupling**: Calculates fluid forces (pressure + viscous shear) to drive the structural dynamics.
- **Parallel Compatible**: Works with MPI parallel decomposition.
- **Diagnostics**: Outputs modal forces, displacements, and optional per-face contributions.

### Quick Usage

1. **Compile**: `wmake libso` inside the source directory.
2. **Load**: Add `libs ("libfastDynamicFvMesh.so");` to `system/controlDict`.
3. **Configure**: Set `dynamicFvMesh fastDynamicFvMesh;` in `constant/dynamicMeshDict` and configure `fastDynamicFvMeshCoeffs`.
4. **Input**: Place mode shape CSV files in the `mode/` directory.

## 2. Bug Fixes and Improvements

The following issues have been fixed and improvements implemented to ensure stability and correctness:

### Library Fixes

- **Dictionary Reading**: Fixed `readControls()` to correctly read `constant/dynamicMeshDict` using `IOdictionary`, replacing the failing `lookupObject` call.
- **Parallel Execution**: Fixed a segmentation fault caused by uninitialized `csvShapes` arrays on slave processors before broadcasting.
- **Path Resolution**: Updated file reading to use absolute paths, ensuring mode files are found correctly in both serial and parallel runs.
- **CSV Parsing**: Rewrote `readModeShapes` to robustly handle trailing commas and whitespace, preventing read failures.
- **Point Mapping**: Relaxed the point mapping tolerance from `1e-12` to `1e-8` to resolve "Mapped 0 points" errors due to floating-point discrepancies.
- **Coupling Relaxation**: Relax mesh displacement every iteration.
- **Reference Pressure**: Remove the effect of refernece pressure.
- **Implicit Update Logic**: Added a check in `update()` to ensure structural dynamics are advanced only once per time step, preventing incorrect acceleration during sub-cycling (e.g., in PIMPLE loops).

### Physics & Scaling

- **Force Scaling**: Fixed negligible forces in incompressible solvers by adding `rhoRef` (reference density) support. Forces are now correctly scaled by density when using kinematic pressure.
- **Pressure Reference**: Added `pRef` parameter to subtract a reference pressure from the field before force calculation.

### Validation Case & Stability

- **Boundary Conditions**: Implemented a robust parabolic inlet using `timeVaryingMappedFixedValue` to avoid `codedFixedValue` compilation issues.
- **Solver Settings**: Added missing `pcorr` (pressure correction) solvers to `system/fvSolution` for valid dynamic mesh motion.
- **Diagnostics**: Implemented detailed runtime diagnostics, including `modal_diagnostics.csv` (forces/displacements per mode) and optional face-level diagnostic files.
