# AGARD FDM Notes

## Where To Look

- `01_assemble_mesh.py`: far-field C3D8 + wing C3D10 merge, renumbering,
  volume/mid-edge/KD-tree assertions, merged metadata.
- `02_generate_modal_case.py`: element-type-aware FSI surfaces, tie, fixed
  boundaries, calibrated wing material, in-file far-field layering, read-back
  verification.
- `03_export_modes.py`: verification gate first (frequency window, pure-wing
  frequency delta, mode identity), then the legacy FDM CSV contract.
- `run_pipeline.sh`: sorted/validated 3-stage driver, `--farfield-profile`
  (legacy/fluid), `--pure-wing` reference, `--dry-run`.
- `tests/test_workflow.py`: synthetic CLI/format coverage plus an env-gated
  byte-level comparison with the signed-off reference deck archive
  (`/root/Workspace/AGARD/Mode/reference/AGARD4456_reference_decks.zip`).
- `mesh/PROVENANCE.md`: committed input archives, source hashes, far-field
  trimming rule and the rebuild snippet.
- `README.md`: the far-field design rule, gate criteria and measured targets.

## Invariants

- Legacy-profile merge counts: 1,365,551 nodes = 388,465 far-field + 977,086
  wing; elements C3D8 `1..372760`, C3D10 `372761..905133`; wing volume
  `4.4542836953e-3 m^3` (rel `1e-6`); wing NSETs
  `ROOT(193)/TIP(212)/SIDE(510058)/WING(977086)`.  `legacy` is the default
  profile and the signed-off reference archive pins only it; its stage-1
  default is `mesh/farfield_c3d8.zip` and a default run is byte-stable.
- fluid profile (`--farfield-profile fluid`): **no default deck** and no
  pinned far-field counts -- `--farfield-deck` is required and the far-field
  node count, id gaps, per-type element split and `FSI_NSET` size are derived
  from that deck and recorded in `mesh/merged_meta.json` (`far.nodes`,
  `far.node_id_gaps`, `far.node_id_max`, `far.element_type_counts`,
  `far.fsi_nodes`, `far.ignored_sets`).  Structural gates stay: positive
  strictly-ascending node ids rank-renumbered to `1..far_nodes`; all eight
  contract NSETs present and non-empty; `FAR_PART_NSET` exactly the full node
  id set; element connectivity inside and exactly covering that set; only
  `C3D8`/`C3D6` blocks on `FAR_PART_ELSET`, ascending ids per block, any
  source order (the writer groups by type).  Unknown NSET names are ignored
  and recorded, so a new ICEM alias needs no code change; the wing side stays
  pinned exactly as for `legacy`.  Measured members: Fluid7-shaped
  (`mesh/fluid7_c3d8_c3d6.zip`) 280,770 far nodes (4,690 id gaps over
  `1..285460`), 268,384 C3D8 + 748 C3D6, 1,816 `FSI_NSET`, ignored
  `All`/`WING_FAR_HALF_1/WING_FAR_HALF_FARD_1_NSET`, merged 1,257,856 nodes /
  801,505 elements, max FSI-to-`SIDE` 0.017917 m; Fluid8 (522,300 far nodes,
  7,229 id gaps over `1..529529`), 503,504 C3D8 + 1,496 C3D6, 1,817
  `FSI_NSET`, extra ignored `POINTS_NSET`, merged 1,499,386 nodes /
  1,037,373 elements, max FSI-to-`SIDE` 0.016217 m.  The wedges are a thin
  tip/wake sheet inside `FAR_PART_ELSET` taking the bricks' layer modulus (no
  separate material); Fluid8 wedge faces include `S3`/`S4`/`S5` as well as
  `S1`, so no stage assumes `S1`-only.  The combined stem gains an `_fluid`
  suffix.  Measured on the Fluid7-shaped deck, stage 2 yields 1,774 tie faces
  (1,752 + 22, all `S1`), 20 layers of 13,457/13,456 covering all 269,132 far
  elements and 783,086 wing faces, with `E_max = 3518.452897 Pa`.
- Set names crossing stages: far-field `FAR_PART_ELSET`/`FAR_PART_NSET`,
  inner `FSI_NSET`, outer `FRONT/BACK/IN/OUT/TOP/BOTTOM` (`FIXED_*` = outer
  minus FSI), wing `WING` blocks and `ROOT/TIP/SIDE/WING` nodes.  Stage 3
  accepts the aliases `FAR_PART`/`WING_PART` for the old decks.
- C3D10 faces use corner nodes only, table
  `S1=(1,2,3) S2=(1,4,2) S3=(2,4,3) S4=(3,4,1)`; C3D8 faces are the six quads.
- Wing material: `*ELASTIC, TYPE=ENGINEERING CONSTANTS` with the single
  `*ORIENTATION` data line `1,1,0, 0,0,1`; values only from
  `../CaseDescription.md` (`soft_anchored`/`hard_e1_0p75`).  `TYPE=ORTHO`,
  `SOLVER=`, and the isotropic `8.17e9/0.31` wing card are forbidden.
- Far-field: bounded graded band `E(V)=clip(E_max*(V_min/V)^p, E_max/floor,
  E_max)` with `E_max = factor*E_ref` (default `1e-6*E1`), `*DENSITY 0.` and
  Poisson `0.` for every layer; `E_max == factor*E_ref` (rel `1e-12`) and
  `span <= floor` are asserted.
- Frequency deck: `*FREQUENCY` `10,0.334,1.E6` (no `STORAGE=`), `*NODE FILE U`,
  `ROOT,1,3` plus the six `FIXED_*` at DOF 1-3.
- Stage 3 gates: exactly `modes` positive strictly ascending frequencies, the
  1..1e4 Hz sanity window, `max |rel delta| <= 0.01` against the pure wing, and
  mode-identity cosine `>= 0.9` with a diagonal argmax.  A failing gate writes
  no CSV and exits non-zero.
- CSV contract: ascending global node id order, three-field header
  `frequency,count,modes,` with trailing commas, legacy `FluidPara.csv`.
- Inputs: stage 1 (and the stage-2 `--wing-source` reader) opens plain, `.gz`
  or `.zip` decks; a zip must hold exactly one member.  The far-field default
  exists only for `legacy` (`mesh/farfield_c3d8.zip`, the default profile);
  `fluid` requires an explicit `--farfield-deck` (a fluid-family export, e.g.
  the committed `mesh/fluid7_c3d8_c3d6.zip` fixture); the wing default is
  `mesh/wing_c3d10.zip`.  The wing zip is the byte-exact raw source and the
  legacy far-field zip is the trimmed node/C3D8/NSET subset documented in
  `mesh/PROVENANCE.md`.
- Git LFS owns `case/AGARD/Mode/mesh/*.zip`
  (`filter=lfs diff=lfs merge=lfs -text`); never bypass it or re-add the
  archives as plain blobs.
- Never commit derived artifacts: `merged_base.inp`, `combined_*.inp`,
  `pure_*.inp`, or any `*.frd`/`*.dat`/`*.sta`/`*.csv`.  Stages 1-2 rebuild
  them in ~20 s from the zips and they are huge.
- `/root/Workspace/AGARD/Mode/reference/AGARD4456_reference_decks.zip` is the
  read-only signed-off equivalence reference for the env-gated suite; never
  modify it and never re-derive it in place (the suite extracts it under
  `/root/Workspace/tmp/agard-mode-tests/reference`).
- Generated CalculiX/FRD/DAT/CSV/logs and bytecode live under
  `/root/Workspace`, never in this directory.

## Commands

```bash
./run_pipeline.sh --dry-run
./run_pipeline.sh --farfield-profile fluid --farfield-deck /path/to/fluid_icem.inp --dry-run
./run_pipeline.sh --steps 1,2,3 --run-ccx --pure-wing
AGARD_RUN_SLOW_TESTS=1 python3 -m unittest tests/test_workflow.py
PYTHONPYCACHEPREFIX=/root/Workspace/tmp/pycache python3 -m unittest tests/test_workflow.py
bash -n run_pipeline.sh
```

## Anti-Patterns

- Never let the pseudo-elastic far field be both grounded and stiff: a stiff
  grounded medium stiffens the wing (+594.8 % at a uniform 1e5 Pa).
- Never emit an unbounded modulus band: the retired `8.17e6 .. 3.37e-10 Pa`
  20-layer profile decoupled the far field from its clamped boundary and its
  +23.65 % shift was an artifact.
- Never trust exit codes for correctness: the deck assertions, the stage-2
  read-back and the stage-3 gates are the contract.
- Do not treat stale `.dat`/`.frd`/mode files as production inputs, and do not
  run `ccx` as routine validation (each combined run is ~6650-6770 s /
  37.30 GiB).
- Do not modify the read-only reference archive
  (`/root/Workspace/AGARD/Mode/reference/AGARD4456_reference_decks.zip`) or
  write generated artifacts into the repository.
