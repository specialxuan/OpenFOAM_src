# AGARD FDM Modal Workflow

Source-only three-stage pipeline that turns the second-order AGARD 445.6 wing
deck plus the far-field hexahedral deck into a verified CalculiX modal deck and
then into the FDM mode CSVs. The two input meshes are committed here under
`mesh/` (Git LFS); every generated artifact belongs under `/root/Workspace`,
never in this directory.

## Tree

```text
01_assemble_mesh.py        merge far-field C3D8 + wing C3D10 -> mesh-only deck + metadata
02_generate_modal_case.py  element-aware tie surfaces, boundaries, calibrated materials
03_export_modes.py         verification gate FIRST, then FDM CSV export
run_pipeline.sh            sorted/validated 3-stage driver with dry-run
tests/test_workflow.py     synthetic fast suite + env-gated reference-archive equivalence
mesh/*.zip                 committed input meshes (Git LFS)
mesh/PROVENANCE.md         sources, hashes, trimming rule, verify recipe
```

## Inputs

| Archive | Member | Source (read-only) | Content |
| --- | --- | --- | --- |
| `mesh/farfield_c3d8.zip` | `farfield_c3d8.inp` | former AGARD 445.6 assembly export (source revision sha256 `b85e44089d65d3d458b248dac3760ec9e7f824c462be90dabe8279101b059792`); this archive member is the canonical copy | trimmed far-field: 388,465 nodes / 372,760 C3D8 |
| `mesh/wing_c3d10.zip` | `wing_c3d10.inp` | former `/root/Workspace/AGARD/WingMode2/mesh.inp`, tree retired 2026-09-24 — this archive member is now the only copy | 977,086 nodes / 532,373 C3D10, byte-exact raw copy |

The meshes descend from an ICEM/ACIS CAD workflow that no longer exists on
this machine, so they cannot be produced by a simple command; they are
versioned next to the scripts and stored in Git LFS (`.gitattributes`:
`case/AGARD/Mode/mesh/*.zip filter=lfs diff=lfs merge=lfs -text`).
`mesh/PROVENANCE.md` records each source sha256 and size, the exact far-field
trimming rule, the resulting counts, the archive sha256/size, and a
copy-pasteable Python snippet that verifies both archives.

Stages 1-2 accept plain, `.gz` and single-member `.zip` decks transparently;
the pipeline defaults (`run_pipeline.sh`) and the stage-1/stage-2 defaults
point at the two archives above, so `./run_pipeline.sh --steps 1,2` needs no
`/root/Workspace` inputs.

Derived artifacts are regenerable in ~20 s by stages 1-2 and are huge, so
they must NOT be committed: `mesh/merged_base.inp`, `modal/combined_*.inp`,
`modal/pure_*.inp`, and every `*.dat`, `*.frd`, `*.sta` or `*.csv` (solved
results, FDM exports).  They belong under `--work-root` (default
`/root/Workspace/AGARD/Mode-run`).

## Pipeline

| Stage | Input | Output |
| --- | --- | --- |
| 1 | `mesh/farfield_c3d8.zip`, `mesh/wing_c3d10.zip` | `mesh/merged_base.inp` + `mesh/merged_meta.json` |
| 2 | merged deck | `modal/combined_<wing>_<far-field>.inp` + meta (optionally `modal/pure_<wing>.inp`) |
| 3 | solved `<stem>.dat`/`.frd` | `mode/FluidNodeCoor.csv`, `StructNodeCoor.csv`, `FluidNodeDisp*.csv`, `StructNodeDisp*.csv`, `FluidPara.csv` |

Stage 1 renumbers wing nodes by ascending original id into
`388466..1365551` and elements into C3D8 `1..372760`, C3D10 `372761..905133`;
it asserts 1,365,551 nodes, wing volume `4.4542836953e-3 m^3` (rel `1e-6`),
every C3D10 mid-edge node at its corners' midpoint (coordinate-relative `1e-9`)
and all 4045 `FSI_NSET` nodes within the tie tolerance of a wing `SIDE` node
(max measured `0.01831 m` at tolerance `0.05`).

Stage 2 emits `<elem>,S<k>` tie faces with the correct element type: C3D8
quads when all four nodes are in the far-field inner set, C3D10 faces when the
three CORNER nodes are in `SIDE` (face table `S1=(1,2,3) S2=(1,4,2)
S3=(2,4,3) S4=(3,4,1)`; mid-edge nodes belong implicitly).  It writes
`*TIE ... POSITION TOLERANCE=<tol>, ADJUST=NO`, the six `FIXED_*` sets plus
`ROOT` at DOF 1-3, and `*STEP` / `*FREQUENCY 10,0.334,1.E6` (no `STORAGE=`) /
`*NODE FILE U` / `*END STEP`.

## Run

```bash
./run_pipeline.sh --dry-run                  # print all stage commands, write nothing
./run_pipeline.sh                            # stages 1,2 (no solver)
./run_pipeline.sh --steps 1,2,3 --run-ccx --pure-wing
./run_pipeline.sh --wing-type hard
./run_pipeline.sh --far-field-mode uniform --far-field-modulus 28
./run_pipeline.sh --steps 3 --results-dir /root/Workspace/AGARD/Mode-run/modal \
  --results-stem combined_soft_layered      # export an already-solved deck
```

`--pure-wing` builds `modal/pure_<wing-type>.inp` from the wing source deck;
with `--run-ccx` it also solves it, and stage 3 then defaults `--pure-dat` /
`--pure-frd` to those files.  That pure-wing 10-mode reference is what the
stage-3 frequency and identity gates compare against.

## Materials (`../CaseDescription.md` is authoritative)

| Set | Use | E1 (Pa) | E2 | E3 | G12 | G13 | G23 | rho |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `soft_anchored` | `--wing-type soft` (Model No.3) | 3.518452897e9 | 2.673062111e8 | 1.313278384e9 | 3.904665311e8 | 4.201475033e8 | 9.268418440e7 | 413.5 |
| `hard_e1_0p75` | `--wing-type hard` (Model No.1) | 1.013605240e10 | 7.295992132e8 | 1.738397081e9 | 2.341301081e8 | 5.652122962e8 | 8.280724206e8 | 495.5 |

Both use `nu12=0.034`, `nu13=0.033`, `nu23=0.326` and
`*ELASTIC, TYPE=ENGINEERING CONSTANTS` with the single `*ORIENTATION` data line
`1,1,0, 0,0,1` (axis1 span -> E1, axis2 thickness -> E2, axis3 chord -> E3).
`E3/E2 < 1/nu23^2` is asserted (soft 4.913, hard 2.383 against 9.4095).

## Far-field design rule

The pseudo-elastic far field must stay BOTH grounded (connected to its clamped
outer boundary) and soft.  A stiff grounded medium stiffens the enclosed wing:
a uniform `1e5 Pa` far field shifted the 10 wing modes by up to **+594.8 %**.
The old 20-layer band is therefore forbidden: its `8.17e6 .. 3.37e-10 Pa`
profile accidentally DECOUPLED the far field from its `FIXED_*` boundary, so
its `+23.65 %` shift was an artifact, not physics; and an unbounded modulus
span (2.4e16) has no defensible calibration.

The supported rule is bounded and graded:

```text
E(V) = clip(E_max * (V_min / V)**p, E_max / floor, E_max)
E_max = layer-emax-factor * E_ref      (default 1e-6 * E1 = 3518.452897 Pa)
```

* `--far-field-mode layered` (default): 20 volume-ranked layers, one
  `*MATERIAL` + `*ELSET` + `*SOLID SECTION` per layer, `p=0.25`,
  `floor=1000` -> `E_min = 3.518452897 Pa`, span exactly 1000.  Measured soft
  wing target: `+0.2156 %` max frequency delta (mode 1), all 10 modes inside
  the 1 % gate.
* `--far-field-mode uniform`: one massless medium over the whole
  `FAR_PART_ELSET`; the calibrated choice is `E = 28 Pa`, measured
  `+0.9690 %` max delta.

Every layer is massless (`*DENSITY 0.`) with Poisson `0.`.  Stage 2 asserts
`E_max == factor * E_ref` (rel `1e-12`) and `span <= floor`; it never emits the
old `3.37e-10 Pa` outer layers or an unbounded band.

## Stage-3 gate criteria

No CSV is written unless all of these pass:

1. exactly `--modes` frequencies in the `.dat`, all positive, strictly ascending;
2. every frequency inside `--min-frequency .. --max-frequency` (default 1..1e4 Hz);
3. with `--pure-dat`: `max |f_combined/f_pure - 1| <= --max-rel` (default 0.01),
   with the full per-mode table always printed;
4. with `--pure-frd`: each combined mode's wing-field cosine against its best
   pure-wing mode `>= --min-cosine` (default 0.9) and the best match diagonal.

On failure the stage prints `[AGARD ERR]` and exits non-zero with no CSV.  The
verification runs first; the legacy export contract is unchanged (ascending
global node id order, `frequency,count,modes,` three-field header with trailing
commas, plus the legacy `FluidPara.csv`).

## Verified reference and equivalence

`/root/Workspace/AGARD/Mode/reference/AGARD4456_reference_decks.zip`
(read-only, deterministic; documented in `reference/README.md`) is the
signed-off ground truth:
`soft/combined_soft_layered.inp` (sha256 `a4626645…a42a6d`),
`soft/combined_soft_fix.inp` (uniform 28 Pa) and `soft/pure_soft.inp`.
Regenerating stages 1-2 for the soft wing from the committed `mesh/` archives
reproduces the structural content of all three decks **byte-for-byte on the
`*NODE`..EOF data-line projection**
(comment lines, i.e. those starting with `**`, are excluded; the stage-1
projection is a strict byte prefix of `mesh/merged_base.inp` inside that
archive).  Comment
lines are excluded because they carry no solver data and their provenance text
legitimately changes when a producer script is corrected (the material-slot
marker now names `case/AGARD/Mode/02_generate_modal_case.py`); every structural
byte stays pinned.  The pure-wing deck has no stale comment text and keeps a
whole-file sha256 equality.  Real provenance (source paths, sha256, every
measured assertion) lives in `mesh/merged_meta.json`.  The check is automated:

```bash
AGARD_RUN_SLOW_TESTS=1 python3 -m unittest tests/test_workflow.py
```

Reference results: 4,002,822 equations after the tie, 146,920,878 lower-nnz,
0 unmatched tie nodes; each combined run costs ~6650-6770 s wall and 37.30 GiB
peak RSS, dominated by SPOOLES fill-in and essentially modulus-independent
(the old band, uniform and layered all land in the same band).

## Checks

```bash
PYTHONPYCACHEPREFIX=/root/Workspace/tmp/pycache python3 -m unittest tests/test_workflow.py
bash -n run_pipeline.sh
./run_pipeline.sh --dry-run
```
