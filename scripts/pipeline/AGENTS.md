# Pipeline Scripts Knowledge

## OVERVIEW

Case-generation, modal-export, FSI-run, reconstruction, and analysis automation for dam-failure/mesh-independence workflows.

## STRUCTURE

| Path | Purpose |
|---|---|
| `README_Automation.md` | Pipeline guide, stages, guards, stale-output policy |
| `run_mesh_independence_automation.sh` | Main coarse/medium/fine staged orchestrator |
| `run_damfailure_pipeline.sh` | Single mesh-level pipeline runner |
| `generate_blockMesh_cases.sh` | Regenerates cases from baseline zip with `blockMesh`, `setFields`, `checkMesh` |

## WHERE TO LOOK

| Task | Location | Notes |
|---|---|---|
| Full staged workflow | `run_mesh_independence_automation.sh` | `check|mesh|modal|fsi|reconstruct|analyze|all` |
| Single-case flow | `run_damfailure_pipeline.sh` | `mesh|modal|export|fsi|all` |
| Mesh-level cases | `generate_blockMesh_cases.sh` | Rebuilds meshes; does not copy old `constant/polyMesh` |
| Modal bridge | `../modal/export_calculix_modes_to_fastdynamic.py` | Creates `mode/FluidNodeCoor.csv` and `FluidNodeDisp*.csv` |
| Mesh converter | `../mesh/openfoam_mesh_to_cdb.py` | Shared CDB input generation |

## CONVENTIONS

- Baseline zip defaults come from `/root/Workspace`; scripts should extract/copy into run dirs, not edit baselines.
- `run_mesh_independence_automation.sh --stage check` is the cheap preflight before heavy FSI runs.
- FSI stage must reject mesh/mode point-count mismatches unless `--allow-stale-mode` is explicitly used for diagnostics.
- `generate_blockMesh_cases.sh` regenerates mesh cases with `blockMesh`, `setFields`, and `checkMesh`.
- Parallel FSI paths use `decomposePar -> mpirun -> reconstructPar`.
- `run_damfailure_pipeline.sh` may temporarily override `system/controlDict` `endTime`; it must restore the original value.

## ANTI-PATTERNS

- Do not copy old meshes forward as a substitute for regenerating with `blockMesh`.
- Do not treat stale CalculiX outputs, stale mode files, or old FSI logs as production inputs.
- Do not sync partial `reconstructPar` output into a main case after failure.
- Do not bypass dependency checks for `python3`, `unzip`, OpenFOAM tools, solver commands, or parallel tools.
