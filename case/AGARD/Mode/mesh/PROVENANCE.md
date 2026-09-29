# Input mesh provenance

`farfield_c3d8.zip` and `wing_c3d10.zip` are the two default input meshes of
the AGARD Mode pipeline; the fluid far-field family (`fluid7_c3d8_c3d6.zip`,
`fluid8_c3d8_c3d6.zip`) holds the optional newer far-field contracts, selected
through the derived `fluid` profile (see below).  They are
versioned here because the meshes cannot be produced by a simple command (they
descend from an ICEM/ACIS CAD workflow that no longer exists on this machine),
so the pipeline is only reproducible if the inputs travel with it.  Everything
downstream (`merged_base.inp`, `combined_*.inp`, `pure_*.inp`, solver
`.dat`/`.frd`/`.sta`, FDM CSVs) is derived and must never be committed; stages
1-2 regenerate it in about 20 s.

All four archives are tracked by Git LFS (repository `.gitattributes`:
`case/AGARD/Mode/mesh/*.zip filter=lfs diff=lfs merge=lfs -text`).
Each archive holds exactly one member, with a fixed 1980-01-01 timestamp and
no directory entries; the member bytes are canonical (LF) ASCII.

| Archive | Member | Zip bytes | Zip sha256 |
| --- | --- | --- | --- |
| `farfield_c3d8.zip` | `farfield_c3d8.inp` | 11,198,790 | `cc277631ef48c28820f5155fdbb0d63ab1cb4061658df5ef836564e92aa27b06` |
| `wing_c3d10.zip` | `wing_c3d10.inp` | 34,539,045 | `140779cb3e17fd8f7fb01eade6db4a15cf36e852867124c0e1f873f324dbb7d5` |
| `fluid7_c3d8_c3d6.zip` | `fluid7_c3d8_c3d6.inp` | 8,046,693 | `f6a6dbb4a2517680eb4bd89c2c4a0bd0aca653b4ead6f3fa50b1df557a03b0fe` |
| `fluid8_c3d8_c3d6.zip` | `fluid8_c3d8_c3d6.inp` | 14,444,106 | `876a5e3186859fe17e0cf6a6315ebc61e11162fe48e4f22e5df67d35f99be2f5` |

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

## Fluid far field family: `fluid7_c3d8_c3d6.zip` and `fluid8_c3d8_c3d6.zip` (raw source copies, CRLF->LF normalized)

The optional newer far-field contracts form one family.  Each instance is a
user-supplied ICEM->ABAQUS fluid export stored raw except for the single
CRLF->LF normalization; the two committed instances are `fluid7_c3d8_c3d6.zip`
(280,770 nodes) and `fluid8_c3d8_c3d6.zip` (522,300 nodes).  Both match the
newer export style: a C3D6 wedge sheet in the same `FAR_PART_ELSET`, a smaller
`FSI_NSET` than the legacy far field, and extra NSET names that the contract
ignores.  The derived `fluid` profile selects an instance with an explicit
deck (`--farfield-profile fluid --farfield-deck mesh/fluid7_c3d8_c3d6.zip` or
`mesh/fluid8_c3d8_c3d6.zip`); the legacy `farfield_c3d8.zip` remains the
default profile and is unchanged.

A zip member must be LF ASCII because stage 1 opens it with
`archive.read(...).decode("ascii")` and applies no newline translation: a `\r`
left in the member would end up in the merged deck.  Storing the normalized
member is therefore what keeps the deck CR-free; CRLF->LF is the single
documented transformation and the raw source is otherwise not renumbered,
reordered or reformatted.

### Instance `fluid7_c3d8_c3d6.zip` (280,770 nodes)

| | |
| --- | --- |
| Source (read-only) | user-supplied ICEM->ABAQUS 4.3.1 fluid export, dated 2026-09-25 20:58 (path deliberately not recorded; the export is identified by the source sha256 below) |
| Source bytes / sha256 | 46,764,219 / `2068e032a217af95339b1788873122e30272ad39cae3bbd50afd8782ea8011e5` |
| Source line endings | pure CRLF: 569,087 `\r\n` pairs, no lone `\r` or `\n` |
| Member | `fluid7_c3d8_c3d6.inp` |
| Member bytes / sha256 | 46,195,132 / `c83e6c8effafd65cbf152d5b6bcbd7d0e55472d1beea5fa92fa92d7958c45a8c` |
| Transformation | CRLF->LF only: the 569,087 `\r\n` pairs become `\n`; the member holds zero `\r` bytes |
| Zip bytes / sha256 | 8,046,693 / `f6a6dbb4a2517680eb4bd89c2c4a0bd0aca653b4ead6f3fa50b1df557a03b0fe` |
| Content | 280,770 nodes (ascending ids `1..285,460` with 4,690 gaps) / 268,384 `C3D8` + 748 `C3D6`, both `ELSET=FAR_PART_ELSET` / NSETs `FAR_PART_NSET` 280,770, `FRONT_NSET` 5,509, `BACK_NSET` 5,661, `IN_NSET` 3,050, `OUT_NSET` 3,650, `TOP_NSET` 2,400, `BOTTOM_NSET` 2,400, `FSI_NSET` 1,816 / ignored aliases `All` and `WING_FAR_HALF_1/WING_FAR_HALF_FARD_1_NSET` (1,312 ids) / placeholder material cards discarded by stage 1 / bbox x -5..6, y 0..7, z -5..5, identical to the legacy far field |

Content inventory (measured 2026-09-25 on the LF member):

* 280,770 nodes, strictly ascending ids `1..285,460`, 4,690 scattered gaps;
* 268,384 `C3D8` then 748 `C3D6`, all in `FAR_PART_ELSET` (source element ids
  23,707..292,882);
* NSET sizes as listed above, including `FSI_NSET` 1,816; the two aliases
  `All` (the `*NODE, NSET=All` set) and
  `WING_FAR_HALF_1/WING_FAR_HALF_FARD_1_NSET` (1,312 ids) are ignored by the
  contract;
* placeholder cards stage 1 discards: `*MATERIAL, NAME=WING` and
  `*ELASTIC, TYPE=ISOTROPIC` (`1, 0.5,`) before `*NODE`, and
  `*SOLID SECTION, ELSET=FAR_PART_ELSET, MATERIAL=WING` after the element
  blocks;
* bounding box x -5..6, y 0..7 (measured minimum `-1.64e-24`, i.e. zero to
  double rounding), z -5..5 -- identical to the legacy far field.

Measured evidence (2026-09-25, from this commit's scripts):

* stage 1 (`--farfield-profile fluid`, `--farfield-deck` this instance):
  merged deck 1,257,856 nodes / 268,384 `C3D8` + 748 `C3D6` + 532,373
  `C3D10` = 801,505 elements; max
  FSI-to-wing-`SIDE` distance `0.017916693293119033 m` at tolerance `0.05`
  with 0 violations;
* stage 2 (soft wing, layered, 20 layers): 1,774 fluid tie faces -- 1,752 on
  `C3D8` bricks (`S1` 192, `S2` 1,500, `S5` 30, `S6` 30) plus the 22 `C3D6`
  wedge faces, all `S1` (the wedge sheet's tie contribution) -- and 783,086
  wing `C3D10` faces; the 20-layer far field covers all 269,132 far elements
  with `E_max = 3518.452897 Pa`;
* equivalence: stage 1 run on the raw CRLF source and on this archive's member
  differs only in the one-line producer comment that stamps the input deck's
  own path (the two invocations named different input files); the `*NODE`..EOF
  data-line
  projection (the repository's equivalence contract) is byte-identical
  (140,178,540 B, sha256
  `f1b8e08c177d7634d3a9e0c37375a54b57aebd2fa920e7be6cb3bedc048194a5`).  When
  the two inputs carry the same basename, the whole merged deck is
  byte-identical (140,179,216 B, sha256
  `363049f29f72df031b0de442cd5def8475f3b3e2032221ed6376637dc716d5dc`).

### Instance `fluid8_c3d8_c3d6.zip` (522,300 nodes)

| | |
| --- | --- |
| Source (read-only) | user-supplied ICEM->ABAQUS 4.3.1 fluid export, dated 2026-09-26 16:45 (path deliberately not recorded; the export is identified by the source sha256 below) |
| Source bytes / sha256 | 87,396,432 / `6bcd7d0b9694415b7945971b6a9366d95df33a0f618f006662cd7bd6f4a6137f` |
| Source line endings | pure CRLF: 1,062,281 `\r\n` pairs, no lone `\r` or `\n` |
| Member | `fluid8_c3d8_c3d6.inp` |
| Member bytes / sha256 | 86,334,151 / `1b05d1fa14bda0ea33812ab251c7b2705c6595d470720f4da53ff92741aeb113` |
| Transformation | CRLF->LF only: the 1,062,281 `\r\n` pairs become `\n`; the member holds zero `\r` bytes |
| Zip bytes / sha256 | 14,444,106 / `876a5e3186859fe17e0cf6a6315ebc61e11162fe48e4f22e5df67d35f99be2f5` |
| Content | 522,300 nodes (ascending ids `1..529,529` with 7,229 gaps) / 503,504 `C3D8` + 1,496 `C3D6`, both `ELSET=FAR_PART_ELSET` / NSETs `FAR_PART_NSET` 522,300, `FRONT_NSET` 8,995, `BACK_NSET` 9,208, `IN_NSET` 4,617, `OUT_NSET` 5,187, `TOP_NSET` 3,420, `BOTTOM_NSET` 3,420, `FSI_NSET` 1,817 / ignored extras `POINTS_NSET` (2 ids), `WING_FAR_HALF_1/WING_FAR_HALF_FARD_1_NSET` (298 ids) and the `*NODE, NSET=All` alias / placeholder material cards discarded by stage 1 / bbox x -5..6, y 0..7, z -5..5, identical to the legacy far field and to `fluid7` |

Content inventory (measured 2026-09-26 on the LF member):

* 522,300 nodes, strictly ascending ids `1..529,529`, 7,229 scattered gaps;
* 503,504 `C3D8` then 1,496 `C3D6`, all in `FAR_PART_ELSET` (source element ids
  37,265..540,768 for the bricks and 35,701..37,196 for the wedges; the two
  blocks leave a 68-id gap between 37,197 and 37,264);
* NSET sizes as listed above, including `FSI_NSET` 1,817; the extra names
  `POINTS_NSET` (2 ids),
  `WING_FAR_HALF_1/WING_FAR_HALF_FARD_1_NSET` (298 ids) and the
  `*NODE, NSET=All` alias are ignored by the contract;
* placeholder cards stage 1 discards: `*MATERIAL, NAME=WING` and
  `*ELASTIC, TYPE=ISOTROPIC` (`1, 0.5,`) before `*NODE, NSET=All`, and
  `*SOLID SECTION, ELSET=FAR_PART_ELSET, MATERIAL=WING` after the element
  blocks;
* bounding box x -5..6, y 0..7 (measured minimum `-2.22e-16`, i.e. zero to
  double rounding), z -5..5 -- identical to the legacy far field and to the
  `fluid7` instance.

Measured evidence (2026-09-26, on the LF member):

* tie feasibility: the 1,817 `FSI_NSET` nodes are at most `0.016217 m` from
  the wing `SIDE` set at tolerance `0.05`, with 0 violations (`fluid7`:
  `0.017916693293119033 m`);
* stage-2 tie surfaces: 1,746 `C3D8` faces plus 40 wedge faces (34 `S1` +
  2 `S3` + 2 `S4` + 2 `S5`) = 1,786 fluid tie faces (`fluid7`: 1,752 + 22
  wedge `S1` = 1,774);
* quality: 0 non-positive Jacobians in the 503,504 bricks; all 1,496 wedges
  consistent under the three-tetrahedron decomposition (the three signed tet
  volumes share one sign); 13,008 bricks below `1e-12 m^3` (stage-2 2x2x2
  Gauss brick volume; thin slivers, same character as the earlier meshes);
* the 1,496 `C3D6` elements form the same thin tip/wake sheet as `fluid7`
  (x `0.8096..1.178`, y `0.762..7`, z `~= +-0.001`; `fluid7` has 748 wedges
  over the same sheet).

### Rebuild/verify recipe (both instances; scratch only, never inside the repository)

```python
import hashlib
import zipfile
from pathlib import Path

# Point SRC at the raw ICEM->ABAQUS fluid export yourself; its path is not
# recorded here by policy.  Identify the source by bytes/sha256: fluid7
# 46,764,219 B / 2068e032... and fluid8 87,396,432 B / 6bcd7d0b....
SRC = Path("<raw ICEM->ABAQUS fluid export>")
STEM = "fluid8_c3d8_c3d6"                # or "fluid7_c3d8_c3d6"
MEMBER = STEM + ".inp"
raw = SRC.read_bytes()
assert raw.count(b"\r") == raw.count(b"\r\n") == raw.count(b"\n")  # pure CRLF
member = raw.replace(b"\r\n", b"\n")
assert b"\r" not in member

info = zipfile.ZipInfo(MEMBER, date_time=(1980, 1, 1, 0, 0, 0))
info.create_system = 3                 # Unix
info.external_attr = 0o100644 << 16    # regular file 0644
info.compress_type = zipfile.ZIP_DEFLATED
with zipfile.ZipFile(STEM + ".zip", "w") as archive:
    archive.writestr(info, member)

print("member", len(member), hashlib.sha256(member).hexdigest())
print("zip   ", Path(STEM + ".zip").stat().st_size,
      hashlib.sha256(Path(STEM + ".zip").read_bytes()).hexdigest())
```

Expected: `fluid7` member 46,195,132 B / `c83e6c8e...`, zip 8,046,693 B /
`f6a6dbb4...`; `fluid8` member 86,334,151 B / `1b05d1fa...`, zip 14,444,106 B /
`876a5e31...`.  With Python 3.12.3 / zlib 1.3 (the environment used here) the
container bytes reproduce exactly -- two builds from the same member are
bit-identical for both instances.  As for the other archives the member is
zlib-independent, while the container sha256 could in principle change with a
different deflate implementation.

## Verify the archives

Run this from `case/AGARD/Mode/`.  It verifies every committed member against
its pinned sha256 value and prints every container sha256 value (the two fluid
instances are additionally reproducible with the rebuild recipe in their
section above).  The legacy source trees (the first-order C3D4 wing export and
the AGARD 445.6 assembly export) were retired 2026-09-24, so
`mesh/wing_c3d10.zip` and `mesh/farfield_c3d8.zip` are the canonical copies of
their members and cannot be rebuilt from source.  Python 3.12 / zlib 1.3 was
used here; the member bytes are zlib-independent, while the container sha256
could in principle change with a different deflate implementation.

```python
import hashlib
import zipfile
from pathlib import Path

# The two legacy input archives are canonical: their source trees were retired
# on 2026-09-24, so this recipe verifies the committed members instead of
# rebuilding them.  The far-field member is the trimmed AGARD 445.6 assembly
# export described by the extraction rule above; the wing member is the
# byte-exact raw source copy.  The two fluid instances are pinned here too
# (their rebuild recipe is in the family section above).
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
    "fluid7_c3d8_c3d6.zip": (
        "fluid7_c3d8_c3d6.inp",
        "c83e6c8effafd65cbf152d5b6bcbd7d0e55472d1beea5fa92fa92d7958c45a8c",
    ),
    "fluid8_c3d8_c3d6.zip": (
        "fluid8_c3d8_c3d6.inp",
        "1b05d1fa14bda0ea33812ab251c7b2705c6595d470720f4da53ff92741aeb113",
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
