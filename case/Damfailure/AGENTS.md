# case/Damfailure Knowledge

## OVERVIEW

Dam-break FSI case generators: water column collapses onto an elastic gate (CalculiX modal structure + `fastDynamicFvMesh`), with two independent coupling pipelines — modal FDM and preCICE.

## STRUCTURE

| Path | Purpose |
|---|---|
| `CaseDescription.md` | Physics/geometry/materials reference (dimensions, VOF pair, gate `E=1e6 Pa`, 10-mode frequency table, BCs) — authoritative for both pipelines |
| `FDM/run_pipeline.sh` | Six-stage driver: `01` mesh → `02` modal → `03` assemble → `04` FDM run → `05` reconstruct/post → `06` video |
| `FDM/01_generate_mesh.py` … `06_export_fdm_video.py` | Numbered stages; `03_assemble_fdm_case.py` owns the OpenFOAM case settings |
| `FDM/display_width.py`, `FDM/test_display_width.py` | Terminal-width helper and its unittest suite |
| `FDM/package_fdm_case.sh` | Case packaging helper |
| `preCICE/run_pipeline.sh` | Five-stage driver: `01` mesh → `03` assemble → `04` preCICE two-solver run → `05` reconstruct/post → `06` video (no stage 02; `run_modal_analysis.py` is optional validation only) |
| `preCICE/04_run_precice.py` | Coupled fluid+solid preCICE execution |

## CONVENTIONS

- Drivers default to full runs against `/root/Workspace/fdm_run` (FDM) or `/root/Workspace/precice_run` (preCICE); always preview with `--dry-run` first.
- Both drivers accept `--steps`, `--case-root`, `--res`, `--z-layers`, `--end-time`, `--nprocs`, `--settings-dir`, `--video/--no-video`; FDM adds `--amr` (off by default).
- preCICE requires fluid `deltaT == preCICE time-window-size` (default `0.0005`); gate material defaults match the FDM case (`E=1e6, nu=0, rho=2500`).
- Without arguments the drivers ask interactively; all stage commands run under the case root, never in this directory.
- Python stages are stdlib-style scripts; verify edits with `python3 -m unittest FDM/test_display_width.py` (run from `case/Damfailure/`).

## ANTI-PATTERNS

- Do not edit stage numbering casually; `run_pipeline.sh` validates and sorts `--steps` against the `01..06` sequence.
- Do not run the preCICE pipeline with mismatched `deltaT`/time-window; the driver does not repair it for you.
- Do not treat preCICE stage `02` as missing — it intentionally does not exist; only FDM runs the embedded modal stage.
- Do not commit generated case output, logs, or video artifacts; they belong under the case root in `/root/Workspace`.
