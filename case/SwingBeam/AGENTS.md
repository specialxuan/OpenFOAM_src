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
| `FDM/03_assemble_fdm_case.py` | Owns the OpenFOAM case settings; injects the `structuralForceEnabled` keys (`nStructuralForces 1; structuralForceFilePrefix StructForce; structuralTargetTolerance 1e-3;`), top-level `fsiTolerance` (default `1e-3`), `couplingRelaxation` (default `0.7`), `modalMass`/`modalDamp` (default `--modal-zeta 0.03`), `--coupling-acceleration` (default None = library `constant`) and `--fsi-disp-tolerance` (default None = solver `0`); `nOuterCorrectors 100` |
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
- FDM extras: `--delta-t`, `--no-excite` (stage 02 `--excite-force` is on by default), `--force-amp`, `--force-freq`, `--force-end-time`, `--modal-zeta`, `--coupling-relax`, `--coupling-acceleration`, `--fsi-disp-tolerance`, `--amr`.
- preCICE extras: `--delta-t`, `--gate-e`, `--gate-nu`, `--gate-rho`, `--force-amp`, `--force-freq`, `--force-end-time`, `--solid-threads`.
- FDM `--write-interval` is a step count; preCICE `--write-interval` is runTime, and preCICE `--delta-t` must equal the preCICE time-window-size (default `0.0005`).
- Without arguments the drivers ask interactively; all stage commands run under the case root, never in this directory.

## KEY TUNING

- FDM: the closed, fully water-filled tank makes the beam added-mass-dominated. The generated defaults are the requested implicit set: `couplingRelaxation 0.7`, top-level `fsiTolerance 1e-3`, `modalMass (1..1)` + `modalDamp (2*zeta*omega_m)` at `zeta=0.03`, `nOuterCorrectors 100`. Measured on `res=30`/`z=4`, `A=0.1 N`, `dt=1e-3`: `0.7` diverges inside the first step (mesh displacement change ×22 per iteration → GAMG SIGFPE at `t=0.001`); the constant-relaxation stable window is `α ≤ ~6e-3` (`α=0.01`+100 iters crashes at `t=0.029`; `nOuterCorrectors 300` crashes at `t=0.021`; halving `dt` only re-samples the startup window). The library now exposes `couplingAcceleration` (`constant`/`aitken`/`IQN-ILS`/`IQN-IMVJ`) plus IQN keys, but all accelerators diverge on this added-mass case, so `constant` + `--coupling-relax` remains the production lever.
- FDM convergence evidence (best full run, `--coupling-relax 0.004`, 2.0 s): all 2000 steps physically converged — last-iteration mesh displacement change < `1e-6 m` (median `1.1e-9 m`), tip `dx` p2p `24.8 mm` (old `0.01` baseline `25.4 mm`). Force residual median `6.7e-4`, 60.2% of steps `< 1e-3`. The remaining "FSI force residual did not converge" prints are ill-conditioned `|ΔF|/|F|` spikes on near-zero modal forces (mode crossing zero, `|F|~3e-4`, absolute change ~`5e-5`) during startup (`t≤0.011`) and force zero-crossing windows — they do not indicate physical non-convergence.
- FDM production config (Task 5): **`constant` + `--coupling-relax 0.004` + `--fsi-disp-tolerance 1e-6`** (no `--coupling-acceleration`). The OR absolute-displacement criterion suppresses the zero-crossing false warnings (short probe: 50/50 steps, warnings 11→0). Accelerators are implemented and formula-verified but not production-viable here: a bounded probe (13 IQN cells + Aitken; `endTime 0.05 s`, `α0 ∈ {1e-4,5e-4,1e-3}`, `iqnMaxHistory ∈ {5,10}`, `iqnFilter ∈ {QR2,none}`, `imvjReuse ∈ {0,5}`) saw every IQN-ILS/IQN-IMVJ config SIGFPE between `t=0.001–0.017 s` and Aitken between `t=0.001–0.023 s`; the published `fsiResidual` stays low mid-window while the internal `||r||` grows ~1.3–2×/step. preCICE's mature Aitken/IQN-ILS fail on the same case.
- preCICE: partitioned coupling requires the 10x stiffer `--gate-e 1.0e7` for stability. At the spec `E=1e6` the coupled run survives ~15 time windows and then aborts with a fluid SIGFPE at about `t=0.008`.
- preCICE additionally pins `p_rgh` TOP to `fixedValue 0`, uses force-only convergence, `correctPhi yes` with PCG/DIC pcorr, and the `NtipAux` + `*EQUATION` MPC for the periodic load.

## ANTI-PATTERNS

- Do not run solvers or OpenFOAM/CalculiX tools in the `src/` repository root; use the run roots under `/root/Workspace/StructTest/`.
- Do not commit generated case output, logs, or video artifacts.
- Do not mix a preCICE `deltaT` with a mismatched time-window; the driver does not repair it.
- Do not treat `E=1e6` as preCICE-stable — use `--gate-e 1.0e7` for the stable smoke.
- Do not edit stage numbering casually; the drivers validate and sort `--steps` against their sequence (FDM `01..06`, preCICE `01,03,04,05,06`).
