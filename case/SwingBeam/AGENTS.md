# case/StructTest Knowledge

## OVERVIEW

Closed, water-filled tank enclosing a centrally-clamped elastic cantilever beam driven by a periodic tip force. Tank `x[0,0.584] y[0,0.365] z[0,0.012] m`; beam `x[0.292,0.304] y[0,0.120] z[0,0.012] m`, `E=1e6 Pa`, `nu=0`, `rho=2500 kg/m3`. Excitation `F(t)=A*sin(2*pi*f*t)` with `A=0.02 N`, `f=2.5 Hz`. Two independent coupling pipelines — modal FDM (`fastDynamicFvMesh`) and direct preCICE (`interFoam` + `ccx_preCICE`). `CaseDescription.md` is authoritative.

## STRUCTURE

| Path | Purpose |
|---|---|
| `CaseDescription.md` | Physics/geometry/materials reference (tank + beam dimensions, `alpha.water=1`, `E=1e6/nu=0/rho=2500`, 10-mode frequency table, periodic tip load, both pipelines) — authoritative for both |
| `FDM/run_pipeline.sh` | Six-stage driver: `01` mesh → `02` modal → `03` assemble → `04` FDM run → `05` reconstruct/post → `06` video |
| `FDM/01_generate_mesh.py` | Fluid `blockMesh` topology plus solid CalculiX hex mesh (beam `6 x 45 x 3`, 1288 nodes / 810 cells) |
| `FDM/02_run_modal_virtual_elastic.py` | CalculiX modal analysis → fluid mode CSVs, and structural `mode/StructNodeCoor.csv` + `mode/StructNodeDisp{m}.csv`; `--excite-force` also writes `mode/StructForce1.csv` |
| `FDM/03_assemble_fdm_case.py` | Owns the OpenFOAM case settings; injects the `structuralForceEnabled` keys (`nStructuralForces 1; structuralForceFilePrefix StructForce; structuralTargetTolerance 1e-3;`) plus `couplingRelaxation` (and `modalMass`/`modalDamp` under `--modal-zeta`) |
| `FDM/04_run_fdm.py` | Runs the `fastDynamicFvMesh` modal-superposition solver |
| `FDM/05_postprocess.py` | Reconstruct + tip-displacement post-processing (`postprocess/fdm_gate_topright_disp.csv`) |
| `FDM/06_export_fdm_video.py` | ParaView frames + MP4 export |
| `FDM/display_width.py`, `FDM/test_display_width.py` | Terminal-width helper and its unittest suite |
| `FDM/package_fdm_case.sh` | FDM case packaging helper |
| `preCICE/run_pipeline.sh` | Five-stage driver: `01` mesh → `03` assemble → `04` preCICE two-solver run → `05` reconstruct/post → `06` video (no stage 02) |
| `preCICE/01_generate_mesh.py` | Fluid + solid meshes and the CalculiX solid deck with the `NtipAux`/`*EQUATION`/`*AMPLITUDE` periodic load |
| `preCICE/03_assemble_precice_case.py` | Assembles `pc_case/` (`precice-config.xml`, `fluid-openfoam`, `solid-calculix`) |
| `preCICE/04_run_precice.py` | Coupled fluid (`interFoam`) + solid (`ccx_preCICE`) preCICE execution |
| `preCICE/05_postprocess.py`, `preCICE/06_export_precice_video.py` | Reconstruct/post + frames/MP4 export |
| `preCICE/run_modal_analysis.py`, `preCICE/package_precice_case.sh` | Optional modal validation and packaging helper |

## CONVENTIONS

- Run roots default to `/root/Workspace/StructTest/fdm` (FDM) and `/root/Workspace/StructTest/precice` (preCICE); always preview with `--dry-run` first.
- Both drivers accept `--case-root`, `--steps`, `--res`, `--z-layers`, `--end-time`, `--write-interval`, `--nprocs`, `--settings-dir`, `--video/--no-video`.
- FDM extras: `--delta-t`, `--no-excite` (stage 02 `--excite-force` is on by default), `--force-amp`, `--force-freq`, `--force-end-time`, `--modal-zeta`, `--coupling-relax`, `--amr`.
- preCICE extras: `--delta-t`, `--gate-e`, `--gate-nu`, `--gate-rho`, `--force-amp`, `--force-freq`, `--force-end-time`, `--solid-threads`.
- FDM `--write-interval` is a step count; preCICE `--write-interval` is runTime, and preCICE `--delta-t` must equal the preCICE time-window-size (default `0.0005`).
- Without arguments the drivers ask interactively; all stage commands run under the case root, never in this directory.

## KEY TUNING

- FDM: the closed, fully water-filled tank makes the beam added-mass-dominated. Stable default is `couplingRelaxation = 0.01`; higher values (`1/0.2/0.05`) blow up in the first time step. Override with `--coupling-relax`.
- preCICE: partitioned coupling requires the 10x stiffer `--gate-e 1.0e7` for stability. At the spec `E=1e6` the coupled run survives ~15 time windows and then aborts with a fluid SIGFPE at about `t=0.008`.
- preCICE additionally pins `p_rgh` TOP to `fixedValue 0`, uses force-only convergence, `correctPhi yes` with PCG/DIC pcorr, and the `NtipAux` + `*EQUATION` MPC for the periodic load.

## ANTI-PATTERNS

- Do not run solvers or OpenFOAM/CalculiX tools in the `src/` repository root; use the run roots under `/root/Workspace/StructTest/`.
- Do not commit generated case output, logs, or video artifacts.
- Do not mix a preCICE `deltaT` with a mismatched time-window; the driver does not repair it.
- Do not treat `E=1e6` as preCICE-stable — use `--gate-e 1.0e7` for the stable smoke.
- Do not edit stage numbering casually; the drivers validate and sort `--steps` against their sequence (FDM `01..06`, preCICE `01,03,04,05,06`).
