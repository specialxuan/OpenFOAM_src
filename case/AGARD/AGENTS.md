# case/AGARD Knowledge

## OVERVIEW

AGARD 445.6 half-span wing case: CalculiX modal calibration of orthotropic wing materials against Yates 1963 experimental frequencies, plus mesh tooling that feeds the FSI workflows.

## STRUCTURE

| Path | Purpose |
|---|---|
| `CaseDescription.md` | Authoritative calibration doc: C3D10 mesh provenance, material axis mapping (`1,1,0, 0,0,1`), target frequencies, calibrated `hard_e1_0p75` / `soft_anchored` material cards, paste-ready CalculiX blocks |
| `Mode/` | Three-stage modal pipeline (mesh merge → deck generation → gated mode export); has its own detailed AGENTS.md |
| `Mode/mesh/` | Committed input mesh decks (`farfield_c3d8.zip`, `wing_c3d10.zip`, `fluid7/8_c3d8_c3d6.zip`) via Git LFS, with `PROVENANCE.md` |
| `Mesh/` | Fluid-mesh tooling: `build_fluid_mesh.py`, `inp_to_polymesh.py`, `smooth_farfield_polymesh.py` |
| `FDM/assemble_fdm_case.py` | Assembles the AGARD OpenFOAM FDM case |

## CONVENTIONS

- `case/AGARD/**/*.zip` is Git LFS (`.gitattributes`: `filter=lfs diff=lfs merge=lfs -text`); the decks are large and cannot be regenerated on this machine (ICEM/ACIS workflow retired).
- Material constants come only from `CaseDescription.md`; CalculiX cards must use `*ELASTIC, TYPE=ENGINEERING CONSTANTS` with the single `*ORIENTATION` line `1,1,0, 0,0,1` — `TYPE=ORTHO` is forbidden (it interprets values as stiffness coefficients and yields wrong soft materials).
- Generated CalculiX/FRD/DAT/CSV artifacts belong under `/root/Workspace`, never in this directory.
- Detailed stage gates, invariants, and pipeline usage live in `Mode/AGENTS.md` and `Mode/README.md`; read those before touching `Mode/`.

## ANTI-PATTERNS

- Do not re-add mesh deck zips as plain git blobs; LFS routing is mandatory.
- Do not fit materials on the retired first-order C3D4 mesh; it shear-locks and distorts mode ordering — the C3D10 mesh is the calibrated baseline.
- Do not assume calibrated moduli match published paper values; the fit absorbs mesh/boundary/model error and is valid for frequency matching only.
