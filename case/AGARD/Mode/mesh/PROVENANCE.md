# Input mesh provenance

`farfield_c3d8.zip` and `wing_c3d10.zip` are the two input meshes of the AGARD
Mode pipeline.  They are versioned here because the meshes cannot be produced
by a simple command (they descend from an ICEM/ACIS CAD workflow that no
longer exists on this machine), so the pipeline is only reproducible if the
inputs travel with it.  Everything downstream (`merged_base.inp`,
`combined_*.inp`, `pure_*.inp`, solver `.dat`/`.frd`/`.sta`, FDM CSVs) is
derived and must never be committed; stages 1-2 regenerate it in about 20 s.

Both archives are tracked by Git LFS (repository `.gitattributes`:
`case/AGARD/Mode/mesh/*.zip filter=lfs diff=lfs merge=lfs -text`).
Each archive holds exactly one member, with a fixed 1980-01-01 timestamp and
no directory entries; the member bytes are canonical (LF) ASCII.

| Archive | Member | Zip bytes | Zip sha256 |
| --- | --- | --- | --- |
| `farfield_c3d8.zip` | `farfield_c3d8.inp` | 11,198,790 | `cc277631ef48c28820f5155fdbb0d63ab1cb4061658df5ef836564e92aa27b06` |
| `wing_c3d10.zip` | `wing_c3d10.inp` | 34,539,045 | `140779cb3e17fd8f7fb01eade6db4a15cf36e852867124c0e1f873f324dbb7d5` |

## Wing: `wing_c3d10.zip` (raw source copy)

| | |
| --- | --- |
| Source (read-only) | retired 2026-09-24 together with `/root/Workspace/AGARD/WingMode2/`; the committed `wing_c3d10.zip` member below is now the only copy |
| Source bytes / sha256 | 113,871,338 / `071d8d0e0f4e93d8015214349e56755868e9d32fbe61168be41dd6e3b7fce5b0` |
| Member | `wing_c3d10.inp` |
| Member bytes / sha256 | 113,871,338 / `071d8d0e0f4e93d8015214349e56755868e9d32fbe61168be41dd6e3b7fce5b0` |
| Content | 977,086 nodes / 532,373 C3D10 elements / `ROOT`(193) `TIP`(212) `SIDE`(510,058) `WING`(977,086) NSETs / placeholder material and step tail |

The wing member is the **raw source file, byte for byte**.  It is not
reserialized, renumbered, reordered or reformatted in any way; only the zip
container is new.  The member sha256 equals the source sha256.  The source
tree was retired on 2026-09-24, so this committed member is the only copy; the
tree's six-set material calibration record is retained in
`../pure_wing/AGARD4456_wing_mesh_and_calibration_provenance.zip` under
`wing_mode_sets/`.

## Far field: `farfield_c3d8.zip` (trimmed from the assembly)

| | |
| --- | --- |
| Source (read-only) | former AGARD 445.6 assembly export; tree retired, this archive member is now the only copy |
| Source bytes / sha256 | 107,117,815 / `b85e44089d65d3d458b248dac3760ec9e7f824c462be90dabe8279101b059792` |
| Member | `farfield_c3d8.inp` |
| Member bytes / sha256 | 63,974,012 / `bbdc3a1f4a22a56d114819eec4a9e9ce60653051597f2e6e36af4d6aa7fd87d7` |
| Content | 388,465 nodes / 372,760 C3D8 (`FAR_PART_ELSET`) |

Extraction rule (implemented by the snippet below):

* read the source as ASCII lines, stripping the CRLF terminators;
* keep the `*NODE` block restricted to node ids `1..388465`;
* keep the complete `*ELEMENT, TYPE=C3D8, ELSET=FAR_PART_ELSET` block;
* keep exactly the NSETs `FAR_PART_NSET`, `FSI_NSET`, `FRONT_NSET`,
  `BACK_NSET`, `IN_NSET`, `OUT_NSET`, `TOP_NSET`, `BOTTOM_NSET`;
* drop the superseded first-order C3D4 wing block and every wing NSET
  (`ROOT_NSET`, `TIP_NSET`, `SIDE_NSET`, `WING_PART_NSET`);
* emit canonical block headers (`*NODE`,
  `*ELEMENT, TYPE=C3D8, ELSET=FAR_PART_ELSET`, `*NSET, NSET=<name>`) and copy
  every data line verbatim, LF-terminated.

Resulting counts: 388,465 nodes; 372,760 C3D8 elements; NSET sizes
`FAR_PART_NSET` 388,465, `FSI_NSET` 4,045, `FRONT_NSET` 6,461,
`BACK_NSET` 8,330, `IN_NSET` 3,450, `OUT_NSET` 4,850, `TOP_NSET` 2,400,
`BOTTOM_NSET` 2,400.  Stage 1 opens the archive directly and treats this
trimmed inventory as its runtime contract: `Assembly3.inp` is only the
trimming source, so a foreign element type, an unexpected element set or NSET
name, or a node id beyond 388,465 is rejected rather than tolerated.

## Verify the archives

Run this from `case/AGARD/Mode/`.  It verifies both committed members against
their pinned sha256 values and prints both container sha256 values.  Both
source trees (the first-order C3D4 wing export and the AGARD 445.6 assembly
export) were retired 2026-09-24, so `mesh/wing_c3d10.zip` and
`mesh/farfield_c3d8.zip` are the canonical copies of their members and cannot
be rebuilt from source.  Python 3.12 / zlib 1.3 was used here; the member
bytes are zlib-independent, while the container sha256 could in principle
change with a different deflate implementation.

```python
import hashlib
import zipfile
from pathlib import Path

# Both input archives are now canonical: their source trees were retired on
# 2026-09-24, so this recipe verifies the committed members instead of
# rebuilding them.  The far-field member is the trimmed AGARD 445.6 assembly
# export described by the extraction rule above; the wing member is the
# byte-exact raw source copy.
DEST = Path("mesh")
EXPECTED = {
    "farfield_c3d8.zip": (
        "farfield_c3d8.inp",
        "bbdc3a1f4a22a56d114819eec4a9e9ce60653051597f2e6e36af4d6aa7fd87d7",
    ),
    "wing_c3d10.zip": (
        "wing_c3d10.inp",
        "071d8d0e0f4e93d8015214349e56755868e9d32fbe61168be41dd6e3b7fce5b0",
    ),
}
for name, (member, member_sha) in EXPECTED.items():
    with zipfile.ZipFile(DEST / name) as archive:
        assert archive.namelist() == [member]
        assert hashlib.sha256(archive.read(member)).hexdigest() == member_sha
    print(name, hashlib.sha256((DEST / name).read_bytes()).hexdigest())
```

## Equivalence evidence (measured 2026-09-24)

A fresh `./run_pipeline.sh --steps 1,2` with only these archives as inputs
reproduces the signed-off reference deck archive
(`../reference/AGARD4456_reference_decks.zip`):

* stage 1 wrote `merged_base.inp` (153,920,913 B); its `*NODE`..EOF data-line
  projection (153,920,293 B, sha256
  `acdbe222f6a2ba8f1fcf48fd13f34efe41f27a1d3dce1bf09aa6d3e52901aa1b`) is a
  strict byte prefix of the `mesh/merged_base.inp` projection inside the
  reference archive (161,790,564 B, sha256
  `fe61a08b8234566f683fe2886ecc90ce778e5e225ababc8e00d0f0977691a71f`);
* stage 2 wrote `combined_soft_layered.inp` (164,294,649 B) whose projection
  sha256 `74beda49b5671731fa9d7369b2e4aac1d238cccdcbc2d7478b6c846ac365ff76`
  equals the `soft/combined_soft_layered.inp` projection inside the reference
  archive; only comment lines (producer attribution) differ;
* the stage-1 metadata matches the signed-off build report on every
  retained count, set size, volume and FSI-distance value (1,365,551 nodes;
  372,760 C3D8 + 532,373 C3D10; wing volume `4.454283695340668e-3 m^3`;
  max FSI-to-`SIDE` distance `0.018314251661946102 m`).  The raw build report
  was retired with the `Mode3/` tree; the env-gated slow suite pins these same
  values.
