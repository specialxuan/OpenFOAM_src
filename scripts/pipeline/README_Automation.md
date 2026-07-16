# Dam-Failure Mesh-Independence Automation

This folder contains a staged automation workflow for the dam-failure FSI
mesh-independence study.

## Folder layout

- `coarse/`, `medium/`, `fine/`: OpenFOAM cases for the three mesh levels.
- `cdb_inputs/`: fluid CDB-like files converted from OpenFOAM `polyMesh`.
- `calculix_modal_inputs/`: generated CalculiX `.inp` files and metadata.
  Solver outputs are produced here when `--run-ccx` is used.
- `calculix_modal_inputs/archive_stale_pre_fastdynamic_output/`: old CalculiX
  outputs retained for reference. These files were generated before fluid and
  solid nodes were both requested in FRD displacement output, so they must not
  be used as production fastDynamicFvMesh mode sources.
- `analysis/`: lightweight CSV summaries from existing mesh, modal, and FSI
  outputs.
- `quick_tests/`: tiny synthetic fixtures for script-level smoke tests.
- `template_case/`: extracted baseline template used by
  `generate_blockMesh_cases.sh`.

Stale FSI logs generated before matching mode files were available are archived
inside the affected case, for example `coarse/archive_stale_pre_fastdynamic_fsi/`.
They are retained for diagnosis only and are excluded from the current analysis
summary by being kept out of the case root.

## Main script

Use:

```bash
./run_mesh_independence_automation.sh --stage check
```

Available stages:

- `mesh`: regenerate `coarse`, `medium`, and `fine` OpenFOAM cases from
  `/root/Workspace/case_damfailure.zip`.
- `modal`: convert OpenFOAM meshes to CDB, generate CalculiX inputs, optionally
  run CalculiX, export modal summaries, and when `--run-ccx` is used, export
  matching `mode/FluidNode*.csv` and `mode/StructNode*.csv` files.
- `fsi`: run OpenFOAM FSI cases, but only if `mode/FluidNodeCoor.csv` matches
  the current OpenFOAM mesh point count.
- `reconstruct`: reconstruct latest parallel results.
- `analyze`: write `analysis/mesh_independence_summary.csv`.
- `all`: run all stages in order.

## Important guard

The FSI stage intentionally refuses to run when the mode node count differs
from the OpenFOAM mesh point count. This prevents false runs where most mesh
points remain unmapped and stationary.

To regenerate matching mode files:

```bash
./run_mesh_independence_automation.sh --stage modal --levels coarse --run-ccx
```

Then run a short smoke FSI:

```bash
./run_mesh_independence_automation.sh --stage fsi --levels coarse --fsi-end-time 0.001
```

Use `--allow-stale-mode` only for diagnostic experiments, not for production
mesh-independence results.

## Case-local myInterFoam runner

`run_myInterFoam_case.sh` is a reusable runner template for an existing
`myInterFoam` OpenFOAM case. Copy it into a case root and run it there:

```bash
cp scripts/pipeline/run_myInterFoam_case.sh /path/to/case/run_myInterFoam.sh
cd /path/to/case
./run_myInterFoam.sh -f 8
```

Modes are `-f` fresh run, `-r` restart from `startTime`, `-c` continue from
`latestTime`, and `-s` solve with an existing decomposition. All successful
modes reconstruct with `reconstructParMesh -constant` followed by
`reconstructPar`.
