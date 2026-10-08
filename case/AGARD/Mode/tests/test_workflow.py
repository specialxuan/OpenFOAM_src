"""Synthetic workflow coverage for the 3-stage AGARD Mode pipeline.

The default suite is fast and writes only under ``/root/Workspace/tmp``; no
CalculiX job is ever launched.  The optional reference-equivalence class is
enabled with ``AGARD_RUN_SLOW_TESTS=1`` and compares the regenerated decks
against the signed-off reference archive at
``/root/Workspace/AGARD/Mode/reference/AGARD4456_reference_decks.zip``
(read-only), which is the strongest available regression guard.  The archive is
extracted lazily into ``/root/Workspace/tmp/agard-mode-tests/reference``; set
``AGARD_REFERENCE_DIR`` to an already-extracted tree to reuse it instead.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile
from dataclasses import replace
from pathlib import Path

import numpy as np

MODE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODE_DIR))
AGARD_ROOT = Path("/root/Workspace/AGARD")
REFERENCE_ARCHIVE = AGARD_ROOT / "Mode" / "reference" / "AGARD4456_reference_decks.zip"
REFERENCE_DIR: Path | None = None  # resolved lazily by the slow suite
SCRATCH = Path("/root/Workspace/tmp/agard-mode-tests")


def load(name: str):
    """Import a stage script whose filename starts with a digit."""
    spec = importlib.util.spec_from_file_location(name, MODE_DIR / f"{name}.py")
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


assemble = load("01_assemble_mesh")
modal = load("02_generate_modal_case")
exporter = load("03_export_modes")


def small_contract() -> assemble.MergeContract:
    """The synthetic two-brick + one-tet contract."""
    return assemble.MergeContract(
        far_nodes=16, far_elems=2, wing_nodes=10, wing_elems=1,
        node_offset=16, elem_offset=2, wing_volume=1.0 / 6.0,
        wing_node_gaps=0, fsi_nodes=4)


CUBE_POINTS = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
               (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
WING_POINTS = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1),
               5: (0.5, 0, 0), 6: (0.5, 0.5, 0), 7: (0, 0.5, 0),
               8: (0, 0, 0.5), 9: (0.5, 0, 0.5), 10: (0, 0.5, 0.5)}
FAR_SET_IDS = {
    "FAR_PART_NSET": list(range(1, 17)),
    "FSI_NSET": [1, 2, 3, 4],
    "FRONT_NSET": [5, 6],
    "BACK_NSET": [7, 8],
    "IN_NSET": [1, 5],
    "OUT_NSET": [2, 6],
    "TOP_NSET": [5, 6, 7, 8],
    "BOTTOM_NSET": [1, 2, 3, 4, 9],
}


def write_farfield(path: Path, scale: float = 0.01) -> None:
    """Far-field deck: two tiny C3D8 bricks and the eight retained NSETs."""
    points = []
    for block, offset in ((1, (0.0, 0.0, 0.0)), (9, (3.0 * scale, 0.0, 0.0))):
        for index, (x, y, z) in enumerate(CUBE_POINTS):
            points.append((block + index,
                           (x * scale + offset[0], y * scale + offset[1],
                            z * scale + offset[2])))
    lines = ["*HEADING\nsynthetic far-field\n*NODE\n"]
    lines.extend(f"{node},{x:.6f},{y:.6f},{z:.6f}\n"
                 for node, (x, y, z) in points)
    lines.append("*ELEMENT, TYPE=C3D8, ELSET=FAR_PART_ELSET\n")
    lines.append("1," + ",".join(str(value) for value in range(1, 9)) + "\n")
    lines.append("2," + ",".join(str(value) for value in range(9, 17)) + "\n")
    for name, ids in FAR_SET_IDS.items():
        lines.append(f"*NSET, NSET={name}\n")
        lines.append(",".join(str(value) for value in ids) + "\n")
    path.write_text("".join(lines), encoding="ascii")


def mixed_contract() -> assemble.MergeContract:
    """The synthetic gapped C3D8+C3D6 far-field contract (fluid shape)."""
    return assemble.MergeContract(
        far_nodes=0, far_elems=0, wing_nodes=10, wing_elems=1,
        node_offset=0, elem_offset=0, wing_volume=1.0 / 6.0,
        wing_node_gaps=0, fsi_nodes=0, profile="fluid",
        far_elem_types=("C3D8", "C3D6"), derive_far=True)


# Source ids 3, 9, 11 are absent: min 1, max 17, count 14, so 3 gaps.
MIXED_FAR_POINTS = {
    1: (0.0, 0.0, 0.0), 2: (0.01, 0.0, 0.0), 4: (0.01, 0.01, 0.0),
    5: (0.0, 0.01, 0.0), 6: (0.0, 0.0, 0.01), 7: (0.01, 0.0, 0.01),
    8: (0.01, 0.01, 0.01), 10: (0.0, 0.01, 0.01),
    12: (3.0, 0.0, 0.0), 13: (3.01, 0.0, 0.0), 14: (3.0, 0.01, 0.0),
    15: (3.0, 0.0, 0.01), 16: (3.01, 0.0, 0.01), 17: (3.0, 0.01, 0.01),
}
MIXED_FAR_SETS = {
    "FAR_PART_NSET": [1, 2, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17],
    "FSI_NSET": [1, 2, 4, 5],
    "FRONT_NSET": [6, 7],
    "BACK_NSET": [8, 10],
    "IN_NSET": [1, 6],
    "OUT_NSET": [2, 7],
    "TOP_NSET": [6, 7, 8, 10],
    "BOTTOM_NSET": [1, 2, 4, 5, 12],
}


def write_mixed_farfield(path: Path,
                         ignored_ids: tuple[int, ...] = (900, 901, 902)) -> None:
    """Gapped C3D8+C3D6 far field plus an ignored NSET alias."""
    lines = ["*HEADING\nsynthetic mixed far-field\n*NODE, NSET=All\n"]
    lines.extend(f"{node},{x:.6f},{y:.6f},{z:.6f}\n"
                 for node, (x, y, z) in MIXED_FAR_POINTS.items())
    lines.append("*ELEMENT, TYPE=C3D8, ELSET=FAR_PART_ELSET\n")
    lines.append("1," + ",".join(str(node) for node in
                                 (1, 2, 4, 5, 6, 7, 8, 10)) + "\n")
    lines.append("*ELEMENT, TYPE=C3D6, ELSET=FAR_PART_ELSET\n")
    lines.append("1," + ",".join(str(node) for node in
                                 (12, 13, 14, 15, 16, 17)) + "\n")
    lines.append("*NSET, NSET=IGNORED_SET\n")
    lines.append(",".join(str(value) for value in ignored_ids) + "\n")
    for name, ids in MIXED_FAR_SETS.items():
        lines.append(f"*NSET, NSET={name}\n")
        lines.append(",".join(str(value) for value in ids) + "\n")
    path.write_text("".join(lines), encoding="ascii")


def write_wing(path: Path) -> None:
    """Second-order wing deck with the placeholder material tail."""
    lines = ["*HEADING\nsynthetic wing\n*NODE\n"]
    lines.extend(f"{node},{x:.6f},{y:.6f},{z:.6f}\n"
                 for node, (x, y, z) in WING_POINTS.items())
    lines.append("*ELEMENT, TYPE=C3D10, ELSET=WING\n")
    lines.append("1," + ",".join(str(value) for value in range(1, 11)) + "\n")
    for name, ids in (("ROOT", [1]), ("TIP", [4]), ("SIDE", [1, 2, 3, 4]),
                      ("WING", list(range(1, 11)))):
        lines.append(f"*NSET, NSET={name}\n")
        lines.append(",".join(str(value) for value in ids) + "\n")
    lines.append("*MATERIAL, NAME=WING\n")
    lines.append("*ELASTIC, TYPE=ENGINEERING CONSTANTS\n")
    lines.append("9.45e9, 6.05e8, 1.01e9, 0.034, 0.033, 0.326, 6.23e8, 8.12e8\n")
    lines.append("2.65e8, 0.\n")
    lines.append("*DENSITY\n495.5,\n")
    lines.append("*ORIENTATION, NAME=WING_AXES\n1,1,0, 0,0,1\n")
    lines.append("*SOLID SECTION, ELSET=WING, MATERIAL=WING, ORIENTATION=WING_AXES\n")
    lines.append(",\n*STEP\n*FREQUENCY\n4,0.334,1.E6\n*BOUNDARY\nROOT,1,3\n*END STEP\n")
    path.write_text("".join(lines), encoding="ascii")


def frd_row(node: int, vector: tuple[float, float, float]) -> str:
    return f" -1{node:10d}{vector[0]:12.5E}{vector[1]:12.5E}{vector[2]:12.5E}\n"


def write_frd(path: Path, mode_vectors: list[dict[int, tuple[float, float, float]]],
              frequencies: list[float]) -> None:
    rows = ["    2C\n"]
    for node, vector in mode_vectors[0].items():
        rows.append(frd_row(node, vector))
    for mode, frequency in enumerate(frequencies, start=1):
        rows.append(f"100CL  101 {frequency} {len(mode_vectors[0])} 2 1MODAL {mode}\n")
        rows.append(f" -4  DISP        {len(mode_vectors[0])}    1\n")
        for node, vector in mode_vectors[mode - 1].items():
            rows.append(frd_row(node, vector))
        rows.append(" -3\n")
    path.write_text("".join(rows), encoding="ascii")


# Verbatim ` -1` records of the T0-captured nodal blocks in
# /root/Workspace/tmp/modal-ext-t0/modal_t0.frd: mode 1 STRESS lines
# 712-717 / TOSTRAIN 961-966, mode 2 STRESS 1701-1706 / TOSTRAIN 1950-1955.
# Copied byte-for-byte because the six `%12.5E` fields pack without
# separators: any column drift must fail the stress tests below.
T0_STRESS_MODE1 = (
    " -1         1-1.20132E+10-5.14855E+09-5.14856E+09-1.41294E+09 6.43629E+02-1.64024E+09\n",
    " -1         2-2.09977E+09-8.99909E+08-8.99908E+08-4.21468E+08 5.10309E+02-7.25513E+08\n",
    " -1         3 7.57655E+09 3.24711E+09 3.24712E+09 7.82593E+08 4.05926E+02-1.14300E+09\n",
    " -1         4-9.27137E+09-3.97347E+09-3.97348E+09-1.64314E+08 1.15574E+02-1.86096E+09\n",
    " -1         5-1.07429E-04-4.61414E-05-4.62234E-05-2.71763E+08-4.37722E-08-1.19995E+09\n",
    " -1         6 9.27137E+09 3.97347E+09 3.97348E+09-1.64314E+08-1.15574E+02-1.86096E+09\n",
)
T0_STRESS_MODE2 = (
    " -1         1 7.57655E+09 3.24712E+09 3.24711E+09 1.14300E+09-4.05926E+02 7.82593E+08\n",
    " -1         2 9.27137E+09 3.97348E+09 3.97347E+09 1.86096E+09 1.15574E+02-1.64314E+08\n",
    " -1         3 1.20132E+10 5.14856E+09 5.14855E+09 1.64024E+09 6.43629E+02-1.41294E+09\n",
    " -1         4-2.09977E+09-8.99908E+08-8.99909E+08 7.25513E+08-5.10309E+02-4.21468E+08\n",
    " -1         5 1.97716E-04 8.47038E-05 8.46451E-05 1.19995E+09 1.76951E-08-2.71763E+08\n",
    " -1         6 2.09977E+09 8.99908E+08 8.99909E+08 7.25513E+08 5.10309E+02-4.21468E+08\n",
)
T0_STRAIN_MODE1 = (
    " -1         1-4.24956E-02-9.75948E-08-1.27271E-07-8.74678E-03 3.98437E-09-1.01539E-02\n",
    " -1         2-7.42775E-03-2.51467E-08-1.64172E-08-2.60909E-03 3.15906E-09-4.49127E-03\n",
    " -1         3 2.68013E-02 4.73824E-08 9.44371E-08 4.84463E-03 2.51288E-09-7.07573E-03\n",
    " -1         4-3.27966E-02-7.24886E-08-1.11033E-07-1.01718E-03 7.15461E-10-1.15202E-02\n",
    " -1         5-3.82053E-16-5.42101E-20-5.01444E-19-1.68234E-03-2.23617E-19-7.42824E-03\n",
    " -1         6 3.27966E-02 7.24886E-08 1.11033E-07-1.01718E-03-7.15461E-10-1.15202E-02\n",
)
T0_STRAIN_MODE2 = (
    " -1         1 2.68013E-02 9.44371E-08 4.73824E-08 7.07573E-03-2.51288E-09 4.84463E-03\n",
    " -1         2 3.27966E-02 1.11033E-07 7.24886E-08 1.15202E-02 7.15461E-10-1.01718E-03\n",
    " -1         3 4.24956E-02 1.27271E-07 9.75948E-08 1.01539E-02 3.98437E-09-8.74678E-03\n",
    " -1         4-7.42775E-03-1.64172E-08-2.51467E-08 4.49127E-03-3.15906E-09-2.60909E-03\n",
    " -1         5 6.99744E-16 3.11708E-19-5.42101E-20 7.42824E-03 1.72795E-19-1.68234E-03\n",
    " -1         6 7.42775E-03 1.64172E-08 2.51467E-08 4.49127E-03 3.15906E-09-2.60909E-03\n",
)
STRESS_DECLARATIONS = (
    " -5  SXX         1    4    1    1\n",
    " -5  SYY         1    4    2    2\n",
    " -5  SZZ         1    4    3    3\n",
    " -5  SXY         1    4    1    2\n",
    " -5  SYZ         1    4    2    3\n",
    " -5  SZX         1    4    3    1\n",
)
STRAIN_DECLARATIONS = (
    " -5  EXX         1    4    1    1\n",
    " -5  EYY         1    4    2    2\n",
    " -5  EZZ         1    4    3    3\n",
    " -5  EXY         1    4    1    2\n",
    " -5  EYZ         1    4    2    3\n",
    " -5  EZX         1    4    3    1\n",
)


def write_frd_with_stress(
        path: Path,
        mode_vectors: list[dict[int, tuple[float, float, float]]],
        stress_rows: tuple[tuple[str, ...], ...],
        strain_rows: tuple[tuple[str, ...], ...],
        frequencies: list[float]) -> None:
    """FRD with the T0 block group per mode: DISP, STRESS, TOSTRAIN, ERROR.

    Every dataset owns a ``100CL`` record and the six-component blocks carry
    the ``-5`` declarations, exactly like the captured ccx output; ERROR holds
    one 25-char component per node and must be skipped by every parser.
    """
    nodes = sorted(mode_vectors[0])
    rows = ["    2C\n"]
    for node in nodes:
        rows.append(frd_row(node, mode_vectors[0][node]))
    for index, frequency in enumerate(frequencies):
        error_rows = tuple(f" -1{node:10d}{999.0:12.5E}\n" for node in nodes)
        blocks = (
            ("DISP", 4, (), tuple(frd_row(node, mode_vectors[index][node])
                                  for node in nodes)),
            ("STRESS", 6, STRESS_DECLARATIONS, stress_rows[index]),
            ("TOSTRAIN", 6, STRAIN_DECLARATIONS, strain_rows[index]),
            ("ERROR", 1, (), error_rows),
        )
        for name, width, declarations, block_rows in blocks:
            header = (f"  100CL  10{index + 1} {frequency} {len(nodes)}"
                      "                     2    1MODAL      1\n")
            rows.append(header)
            rows.append(f" -4  {name:<8} {width}    1\n")
            rows.extend(declarations)
            rows.extend(block_rows)
            rows.append(" -3\n")
    path.write_text("".join(rows), encoding="ascii")


def write_dat(path: Path, frequencies: list[float]) -> None:
    lines = [" MODE NO EIGENVALUE FREQUENCY"]
    lines.extend(f" {index} 1 1 {value}" for index, value in
                 enumerate(frequencies, start=1))
    path.write_text("\n".join(lines) + "\n\n", encoding="ascii")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


NODE_CARD = b"*NODE\n"


def content_from_node(data: bytes) -> bytes:
    """Deck bytes from the first ``*NODE`` card to EOF (header excluded)."""
    return data[data.index(NODE_CARD):]


def header_block(data: bytes) -> bytes:
    """Leading header bytes before the first ``*NODE`` card."""
    return data[:data.index(NODE_CARD)]


def data_line_projection(data: bytes) -> bytes:
    """Deck bytes from the first ``*NODE`` card to EOF, comment lines dropped.

    Comment lines (leading ``**``) carry no solver data, and their provenance
    text legitimately changes when a producer script is corrected (for example
    the material-slot marker).  Comparing only the data lines keeps every
    structural byte pinned -- nodes, elements, sets, surfaces, tie, material
    moduli and step cards -- while allowing accurate comment text.  Do NOT
    "restore" a whole-region equality here: that would re-freeze stale comment
    text, not data.
    """
    region = content_from_node(data)
    return b"".join(line for line in region.splitlines(keepends=True)
                    if not line.strip().startswith(b"**"))


def count_region_lines(data: bytes) -> tuple[int, int]:
    """(comment lines, data lines) in the ``*NODE``..EOF region."""
    lines = content_from_node(data).splitlines(keepends=True)
    comments = sum(1 for line in lines if line.strip().startswith(b"**"))
    return comments, len(lines) - comments


def run_script(script: str, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(MODE_DIR / script), *arguments],
                          text=True, capture_output=True, check=False)


class ScratchCase(unittest.TestCase):
    """Per-test scratch directory outside the repository."""

    def setUp(self) -> None:
        SCRATCH.mkdir(parents=True, exist_ok=True)
        self._temp = tempfile.TemporaryDirectory(dir=SCRATCH)
        self.root = Path(self._temp.name)

    def tearDown(self) -> None:
        self._temp.cleanup()

    def farfield(self) -> Path:
        path = self.root / "far.inp"
        write_farfield(path)
        return path

    def wing(self) -> Path:
        path = self.root / "wing.inp"
        write_wing(path)
        return path

    def merged(self) -> Path:
        far, wing = self.farfield(), self.wing()
        assemble.assemble(far, wing, self.root / "merged.inp",
                          self.root / "merged_meta.json", small_contract())
        return self.root / "merged.inp"

    def combined(self, **overrides) -> Path:
        spec = modal.CaseSpec(**{**{"layers": 2}, **overrides})
        request = modal.BuildRequest(
            merged_deck=self.merged(), output=self.root / "combined.inp",
            meta=self.root / "combined_meta.json", spec=spec)
        modal.generate(request)
        return self.root / "combined.inp"


class AssembleMeshTests(ScratchCase):
    def test_merge_renumbers_when_decks_match_the_contract(self) -> None:
        output, meta = self.root / "merged.inp", self.root / "merged_meta.json"
        report = assemble.assemble(self.farfield(), self.wing(), output, meta,
                                   small_contract())
        self.assertEqual(report["counts"],
                         {"nodes": 26, "c3d8": 2, "c3d10": 1, "elements": 3})
        text = output.read_text(encoding="ascii")
        self.assertIn("*ELEMENT, TYPE=C3D8, ELSET=FAR_PART_ELSET\n", text)
        self.assertIn("*ELEMENT, TYPE=C3D10, ELSET=WING\n", text)
        self.assertIn("3," + ",".join(str(value) for value in range(17, 27)) + "\n",
                      text)
        self.assertIn("*NSET, NSET=FIXED_FRONT\n5,6\n", text)
        self.assertEqual(report["fsi_wing_distance"]["count"], 4)
        self.assertTrue(report["fsi_wing_distance"]["all_within_tolerance"])
        self.assertTrue(json.loads(meta.read_text())["wing"]["volume_m3"] > 0)

    def test_merge_reports_volume_and_mid_edge_errors(self) -> None:
        output = self.root / "merged.inp"
        report = assemble.assemble(self.farfield(), self.wing(), output,
                                   self.root / "meta.json", small_contract())
        self.assertLess(report["wing"]["volume_rel_err"], 1e-6)
        self.assertLessEqual(report["wing"]["mid_edge_rel_coord_max"], 1e-9)
        broken = self.root / "broken_wing.inp"
        text = self.wing().read_text(encoding="ascii").replace(
            "9,0.500000,0.000000,0.500000", "9,0.900000,0.000000,0.500000")
        broken.write_text(text, encoding="ascii")
        with self.assertRaises(assemble.MergeError):
            assemble.assemble(self.farfield(), broken, output,
                              self.root / "meta.json", small_contract())

    def test_merge_rejects_when_fsi_count_differs_from_contract(self) -> None:
        contract = replace(small_contract(), fsi_nodes=3)
        with self.assertRaises(assemble.MergeError):
            assemble.assemble(self.farfield(), self.wing(),
                              self.root / "merged.inp", self.root / "meta.json",
                              contract)

    def test_merge_rejects_when_fsi_nodes_are_far_from_side(self) -> None:
        far = self.root / "far_far.inp"
        write_farfield(far, scale=1.0)
        with self.assertRaises(assemble.MergeError):
            assemble.assemble(far, self.wing(), self.root / "merged.inp",
                              self.root / "meta.json", small_contract())

    def _broken_farfield(self, old: str, new: str) -> Path:
        far = self.root / "far_broken.inp"
        far.write_text(self.farfield().read_text(encoding="ascii").replace(old, new),
                       encoding="ascii")
        return far

    def _reject(self, far: Path) -> assemble.MergeError:
        with self.assertRaises(assemble.MergeError) as caught:
            assemble.assemble(far, self.wing(), self.root / "merged.inp",
                              self.root / "meta.json", small_contract())
        return caught.exception

    def test_merge_rejects_a_foreign_element_type(self) -> None:
        far = self._broken_farfield("TYPE=C3D8", "TYPE=C3D4")
        self.assertIn("element card", str(self._reject(far)))

    def test_merge_rejects_an_unexpected_element_set(self) -> None:
        far = self._broken_farfield("ELSET=FAR_PART_ELSET", "ELSET=OTHER_ELSET")
        self.assertIn("element card", str(self._reject(far)))

    def test_merge_rejects_an_unexpected_nset_name(self) -> None:
        far = self.root / "far_bad_nset.inp"
        far.write_text(self.farfield().read_text(encoding="ascii")
                       + "*NSET, NSET=MYSTERY\n1,\n", encoding="ascii")
        self.assertIn("unexpected far-field NSET", str(self._reject(far)))

    def test_merge_rejects_node_ids_beyond_the_far_field(self) -> None:
        far = self._broken_farfield("\n*ELEMENT, TYPE=C3D8",
                                    "\n17,0.000000,0.000000,0.000000\n"
                                    "*ELEMENT, TYPE=C3D8")
        self.assertIn("node id 17 outside", str(self._reject(far)))

    def test_cli_reports_an_unexpected_nset_as_agard_err(self) -> None:
        far = self.root / "far_bad_nset.inp"
        far.write_text(self.farfield().read_text(encoding="ascii")
                       + "*NSET, NSET=MYSTERY\n1,\n", encoding="ascii")
        result = run_script("01_assemble_mesh.py", "--farfield-deck", str(far),
                            "--wing-deck", str(self.wing()),
                            "--output", str(self.root / "out.inp"),
                            "--meta", str(self.root / "out_meta.json"))
        self.assertEqual(result.returncode, 1)
        self.assertIn("[AGARD ERR]", result.stderr)
        self.assertNotIn("Traceback", result.stderr)


class CompressedInputTests(ScratchCase):
    """Transparent plain / .gz / one-member .zip deck handling."""

    def _merged(self, far: Path, wing: Path, **members: str) -> bytes:
        output = self.root / "merged.inp"
        assemble.assemble(far, wing, output, self.root / "meta.json",
                          small_contract(), **members)
        return output.read_bytes()

    def test_plain_gzip_and_zip_inputs_merge_identically(self) -> None:
        far, wing = self.farfield(), self.wing()
        expected = self._merged(far, wing)
        far_zip = self.root / "far.zip"
        wing_gz = self.root / "wing.inp.gz"
        with zipfile.ZipFile(far_zip, "w") as archive:
            archive.write(far, "far.inp")
        with gzip.open(wing_gz, "wt", encoding="ascii") as handle:
            handle.write(wing.read_text(encoding="ascii"))
        self.assertEqual(data_line_projection(self._merged(far_zip, wing_gz)),
                         data_line_projection(expected))

    def test_named_member_selects_one_entry_of_a_multi_member_zip(self) -> None:
        far, wing = self.farfield(), self.wing()
        multi = self.root / "multi.zip"
        with zipfile.ZipFile(multi, "w") as archive:
            archive.write(far, "far.inp")
            archive.write(wing, "unused.inp")
        self.assertEqual(
            data_line_projection(self._merged(multi, wing,
                                              farfield_member="far.inp")),
            data_line_projection(self._merged(far, wing)))
        with self.assertRaises(assemble.MergeError):
            self._merged(multi, wing)

    def test_corrupt_absent_and_multi_member_zip_fail_clearly(self) -> None:
        far, wing = self.farfield(), self.wing()
        corrupt = self.root / "corrupt.zip"
        corrupt.write_bytes(b"not a zip archive")
        with self.assertRaises(assemble.MergeError):
            self._merged(corrupt, wing)
        with self.assertRaises(assemble.MergeError):
            self._merged(self.root / "missing.zip", wing)
        result = run_script("01_assemble_mesh.py", "--farfield-deck", str(corrupt),
                            "--wing-deck", str(wing),
                            "--output", str(self.root / "out.inp"),
                            "--meta", str(self.root / "out_meta.json"))
        self.assertEqual(result.returncode, 1)
        self.assertIn("[AGARD ERR]", result.stderr)
        self.assertNotIn("Traceback", result.stderr)


class FarfieldProfileTests(ScratchCase):
    """Named-profile far-field parsing: gapped ids, C3D6 block, ignored NSETs."""

    def _mixed(self) -> Path:
        path = self.root / "mixed.inp"
        write_mixed_farfield(path)
        return path

    def _assemble_mixed(self, far: Path) -> tuple[Path, dict]:
        output = self.root / "merged.inp"
        meta = self.root / "merged_meta.json"
        report = assemble.assemble(far, self.wing(), output, meta,
                                   mixed_contract())
        return output, report

    def test_gapped_farfield_rank_renumbers_and_emits_the_c3d6_block(self) -> None:
        output, report = self._assemble_mixed(self._mixed())
        self.assertEqual(report["farfield_profile"], "fluid")
        self.assertEqual(report["counts"], {"nodes": 24, "c3d8": 1, "c3d6": 1,
                                            "c3d10": 1, "elements": 3})
        self.assertEqual(report["far"]["node_id_gaps"], 3)
        self.assertEqual(report["far"]["element_type_counts"],
                         {"C3D8": 1, "C3D6": 1})
        self.assertEqual(report["far"]["element_id_min"], 1)
        self.assertEqual(report["far"]["element_id_max"], 1)
        text = output.read_text(encoding="ascii")
        self.assertIn("*ELEMENT, TYPE=C3D8, ELSET=FAR_PART_ELSET\n"
                      "1,1,2,3,4,5,6,7,8\n", text)
        self.assertIn("*ELEMENT, TYPE=C3D6, ELSET=FAR_PART_ELSET\n"
                      "2,9,10,11,12,13,14\n", text)
        self.assertIn("*ELEMENT, TYPE=C3D10, ELSET=WING\n", text)
        self.assertIn("*NSET, NSET=FAR_PART_NSET\n"
                      "1,2,3,4,5,6,7,8,9,10,11,12,13,14\n", text)
        self.assertIn("*NSET, NSET=FIXED_IN\n5\n", text)
        self.assertIn("*NSET, NSET=BOTTOM_NSET\n1,2,3,4,9\n", text)
        self.assertIn("\n15,0.000000,0.000000,0.000000\n", text)
        self.assertNotIn("IGNORED_SET", text)
        self.assertIn("** far-field profile: fluid\n", text)

    def test_element_blocks_are_accepted_in_any_source_order(self) -> None:
        source = self._mixed().read_text(encoding="ascii")
        c3d8 = ("*ELEMENT, TYPE=C3D8, ELSET=FAR_PART_ELSET\n1,"
                + ",".join(str(node) for node in (1, 2, 4, 5, 6, 7, 8, 10)) + "\n")
        c3d6 = ("*ELEMENT, TYPE=C3D6, ELSET=FAR_PART_ELSET\n1,"
                + ",".join(str(node) for node in (12, 13, 14, 15, 16, 17)) + "\n")
        reordered = self.root / "reordered.inp"
        reordered.write_text(source.replace(c3d8 + c3d6, c3d6 + c3d8),
                             encoding="ascii")
        output, report = self._assemble_mixed(reordered)
        text = output.read_text(encoding="ascii")
        self.assertLess(text.index("*ELEMENT, TYPE=C3D8"),
                        text.index("*ELEMENT, TYPE=C3D6"))
        self.assertEqual(report["far"]["element_type_counts"],
                         {"C3D8": 1, "C3D6": 1})

    def test_derived_profile_records_its_measured_counts_and_ignored_nsets(self) -> None:
        _, report = self._assemble_mixed(self._mixed())
        self.assertEqual(report["far"]["nodes"], 14)
        self.assertEqual(report["far"]["node_id_max"], 17)
        self.assertEqual(report["far"]["elements"], 2)
        self.assertEqual(report["far"]["fsi_nodes"], 4)
        self.assertEqual(report["far"]["ignored_sets"], ["All", "IGNORED_SET"])
        self.assertEqual(report["far"]["node_id_gaps"], 3)

    def test_a_new_nset_alias_needs_no_code_change(self) -> None:
        path = self.root / "new_alias.inp"
        path.write_text(self._mixed().read_text(encoding="ascii")
                        + "*NSET, NSET=BRAND_NEW_ALIAS\n1,2,\n", encoding="ascii")
        output, report = self._assemble_mixed(path)
        self.assertIn("BRAND_NEW_ALIAS", report["far"]["ignored_sets"])
        self.assertNotIn("BRAND_NEW_ALIAS", output.read_text(encoding="ascii"))

    def test_derived_profile_rejects_a_missing_contract_nset(self) -> None:
        path = self.root / "missing_fsi.inp"
        path.write_text(self._mixed().read_text(encoding="ascii").replace(
            "*NSET, NSET=FSI_NSET\n1,2,4,5\n", "*NSET, NSET=FSI_NSET\n"),
            encoding="ascii")
        with self.assertRaises(assemble.MergeError) as caught:
            self._assemble_mixed(path)
        message = str(caught.exception)
        self.assertIn("FSI_NSET", message)
        self.assertIn("empty", message)

    def test_ignored_nset_ids_are_never_validated(self) -> None:
        path = self.root / "mixed_ignored.inp"
        write_mixed_farfield(path, ignored_ids=(900, 901, 902, 903))
        output, report = self._assemble_mixed(path)
        self.assertNotIn("IGNORED_SET", output.read_text(encoding="ascii"))
        self.assertEqual(report["counts"]["nodes"], 24)

    def test_gapped_reference_inside_a_kept_nset_is_rejected(self) -> None:
        source = self._mixed().read_text(encoding="ascii")
        self.assertIn("*NSET, NSET=FSI_NSET\n1,2,4,5\n", source)
        path = self.root / "gapped_nset.inp"
        path.write_text(
            source.replace("*NSET, NSET=FSI_NSET\n1,2,4,5\n",
                           "*NSET, NSET=FSI_NSET\n1,3,4,5\n"),
            encoding="ascii")
        with self.assertRaises(assemble.MergeError) as caught:
            self._assemble_mixed(path)
        message = str(caught.exception)
        self.assertIn("gapped/absent", message)
        self.assertIn("3", message)

    def test_gapped_reference_inside_connectivity_is_rejected(self) -> None:
        source = self._mixed().read_text(encoding="ascii")
        self.assertIn("1,1,2,4,5,6,7,8,10\n", source)
        path = self.root / "gapped_conn.inp"
        path.write_text(
            source.replace("1,1,2,4,5,6,7,8,10\n",
                           "1,1,2,4,5,6,7,3,10\n", 1),
            encoding="ascii")
        with self.assertRaises(assemble.MergeError) as caught:
            self._assemble_mixed(path)
        self.assertIn("gapped/absent", str(caught.exception))

    def test_unsorted_node_ids_are_rejected(self) -> None:
        lines = self._mixed().read_text(encoding="ascii").splitlines(keepends=True)
        start = next(index for index, line in enumerate(lines)
                     if line.startswith("*NODE"))
        lines[start + 1], lines[start + 2] = lines[start + 2], lines[start + 1]
        path = self.root / "unsorted.inp"
        path.write_text("".join(lines), encoding="ascii")
        with self.assertRaises(assemble.MergeError) as caught:
            self._assemble_mixed(path)
        self.assertIn("strictly ascending", str(caught.exception))

    def test_wrong_type_and_wrong_elset_are_rejected(self) -> None:
        source = self._mixed().read_text(encoding="ascii")
        for old, new in (("TYPE=C3D6", "TYPE=C3D4"),
                         ("ELSET=FAR_PART_ELSET", "ELSET=OTHER_ELSET")):
            path = self.root / "broken.inp"
            path.write_text(source.replace(old, new, 1), encoding="ascii")
            with self.assertRaises(assemble.MergeError) as caught:
                self._assemble_mixed(path)
            self.assertIn("element card", str(caught.exception))


class FluidContractTests(unittest.TestCase):
    """The fluid profile derives the far side; legacy stays fully pinned."""

    def test_fluid_contract_derives_the_far_side_and_pins_the_wing(self) -> None:
        contract = assemble.FLUID_CONTRACT
        self.assertEqual(contract.profile, "fluid")
        self.assertTrue(contract.derive_far)
        self.assertIsNone(contract.far_node_id_max)
        self.assertEqual(contract.far_elem_types, ("C3D8", "C3D6"))
        self.assertEqual(contract.wing_nodes, 977086)
        self.assertEqual(contract.wing_elems, 532373)
        self.assertEqual(contract.wing_node_gaps, 12964)
        self.assertEqual(contract.wing_volume, 4.4542836953e-3)

    def test_legacy_contract_keeps_every_pinned_number(self) -> None:
        contract = assemble.LEGACY_CONTRACT
        self.assertFalse(contract.derive_far)
        self.assertEqual(contract.far_nodes, 388465)
        self.assertEqual(contract.far_elems, 372760)
        self.assertEqual(contract.far_node_id_gaps, 0)
        self.assertEqual(contract.far_node_id_max, 388465)
        self.assertEqual(contract.fsi_nodes, 4045)
        self.assertEqual(contract.node_offset, 388465)
        self.assertEqual(contract.elem_offset, 372760)
        self.assertEqual(contract.total_nodes, 1365551)
        self.assertEqual(contract.total_elems, 905133)

    def test_profiles_expose_only_legacy_and_fluid(self) -> None:
        self.assertEqual(sorted(assemble.PROFILES), ["fluid", "legacy"])
        self.assertEqual(assemble.PROFILES["legacy"], assemble.LEGACY_CONTRACT)
        self.assertEqual(assemble.PROFILES["fluid"], assemble.FLUID_CONTRACT)
        self.assertEqual(assemble.AGARD_CONTRACT, assemble.LEGACY_CONTRACT)

    def test_cli_documents_the_profile_choices_and_requires_a_fluid_deck(self) -> None:
        result = run_script("01_assemble_mesh.py", "--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("--farfield-profile", result.stdout)
        self.assertIn("{fluid,legacy}", result.stdout)
        self.assertNotIn("fluid7", result.stdout)
        missing = run_script("01_assemble_mesh.py", "--farfield-profile", "fluid",
                             "--wing-deck", str(MODE_DIR / "mesh" / "wing_c3d10.zip"),
                             "--output", str(SCRATCH / "unused.inp"))
        self.assertEqual(missing.returncode, 2)
        self.assertIn("--farfield-deck", missing.stderr)
        self.assertNotIn("Traceback", missing.stderr)


class SurfaceTests(ScratchCase):
    def test_c3d10_faces_use_corner_table_when_corners_are_in_side(self) -> None:
        mesh = modal.parse_merged(self.merged())
        surfaces = modal.extract_surfaces(mesh)
        self.assertEqual(surfaces.wing_rows.shape[0], 4)
        self.assertEqual(surfaces.wing_flags.tolist(), [1, 2, 3, 4])
        self.assertEqual(surfaces.wing_rows.tolist(), [0, 0, 0, 0])

    def test_c3d8_faces_require_all_four_nodes_in_the_inner_set(self) -> None:
        mesh = modal.parse_merged(self.merged())
        surfaces = modal.extract_surfaces(mesh)
        self.assertEqual(surfaces.fluid_rows.shape[0], 1)
        self.assertEqual(surfaces.fluid_flags.tolist(), [1])

    def test_face_extraction_skips_tets_without_three_corner_side_nodes(self) -> None:
        mesh = modal.parse_merged(self.merged())
        reduced = modal.MergedMesh(
            coords=mesh.coords, far_conn=mesh.far_conn, wing_conn=mesh.wing_conn,
            sets={**mesh.sets, "SIDE": np.array([1], dtype=np.int64)},
            total_nodes=mesh.total_nodes, far_nodes=mesh.far_nodes,
            far_elems=mesh.far_elems)
        with self.assertRaises(modal.DeckError):
            modal.extract_surfaces(reduced)


def wedge_mesh(fsi_ids: list[int], *, brick_fsi: bool = False) -> modal.MergedMesh:
    """Minimal one-brick + one-wedge + one-tet mesh for wedge tie tests.

    Brick deck row 0, wedge deck row 1; the wing tet corners are nodes
    15..18 (all in ``SIDE``) and its mid nodes are 19..24.
    """
    fsi = list(fsi_ids) + ([1, 2, 3, 4] if brick_fsi else [])
    return modal.MergedMesh(
        coords=np.zeros((25, 3)),
        far_conn=np.array([[1, 2, 3, 4, 5, 6, 7, 8]], dtype=np.int64),
        wedge_conn=np.array([[9, 10, 11, 12, 13, 14]], dtype=np.int64),
        wing_conn=np.array([[15, 16, 17, 18, 19, 20, 21, 22, 23, 24]],
                           dtype=np.int64),
        sets={"FSI_NSET": np.array(fsi, dtype=np.int64),
              "SIDE": np.array([15, 16, 17, 18], dtype=np.int64)},
        total_nodes=24, far_nodes=14, far_elems=2)


class WedgeSurfaceTests(unittest.TestCase):
    """C3D6 tie faces: documented CalculiX table and deck-row flags."""

    def test_c3d6_face_table_matches_the_documented_calculix_source(self) -> None:
        self.assertEqual(modal.WEDGE_FACES, (
            ("S1", (0, 1, 2)), ("S2", (3, 4, 5)), ("S3", (0, 1, 4, 3)),
            ("S4", (1, 2, 5, 4)), ("S5", (2, 0, 3, 5))))

    def test_each_wedge_face_flag_is_extracted_on_its_own_nodes(self) -> None:
        for flag, (name, indices) in enumerate(modal.WEDGE_FACES, start=1):
            mesh = wedge_mesh([9 + index for index in indices])
            surfaces = modal.extract_surfaces(mesh)
            # The brick has no FSI face; the wedge is deck row 1 (element id 2).
            self.assertEqual(surfaces.fluid_rows.tolist(), [1], name)
            self.assertEqual(surfaces.fluid_flags.tolist(), [flag], name)

    def test_wedge_face_in_fsi_emits_exactly_one_tie_row_with_its_flag(self) -> None:
        mesh = wedge_mesh([9, 10, 11], brick_fsi=True)
        surfaces = modal.extract_surfaces(mesh)
        self.assertEqual(surfaces.fluid_rows.tolist(), [0, 1])
        self.assertEqual(surfaces.fluid_flags.tolist(), [1, 1])
        wedge_faces = surfaces.fluid_rows >= mesh.brick_elems
        self.assertEqual(surfaces.fluid_rows[wedge_faces].tolist(), [1])
        self.assertEqual(surfaces.fluid_flags[wedge_faces].tolist(), [1])
        buffer = io.StringIO()
        surfaces.write(buffer, 0.05, mesh.far_elems)
        self.assertIn("1,S1\n", buffer.getvalue())
        self.assertIn("2,S1\n", buffer.getvalue())


class WedgeGeometryTests(unittest.TestCase):
    """Exact wedge volumes from the three-tetrahedron decomposition."""

    def test_wedge_volumes_match_an_analytic_prism(self) -> None:
        # Bottom triangle (0,0,0),(1,0,0),(0,1,0) (area 0.5) translated by
        # (0.3,0.2,1.0): an oblique prism of perpendicular height 1 -> V=0.5.
        coords = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
            [0.3, 0.2, 1.0], [1.3, 0.2, 1.0], [0.3, 1.2, 1.0],
        ])
        conn = np.array([[1, 2, 3, 4, 5, 6]], dtype=np.int64)
        volume = modal.wedge_volumes(coords, conn)
        self.assertEqual(volume.shape[0], 1)
        self.assertAlmostEqual(float(volume[0]), 0.5, places=12)

    def test_far_volumes_are_bricks_then_wedges_in_deck_row_order(self) -> None:
        coords = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
            [0.3, 0.2, 1.0], [1.3, 0.2, 1.0], [0.3, 1.2, 1.0],
        ])
        brick = np.array([[1, 2, 3, 4, 5, 6, 7, 8]], dtype=np.int64)
        wedge = np.array([[9, 10, 11, 12, 13, 14]], dtype=np.int64)
        far = modal.far_volumes(coords, brick, wedge)
        self.assertEqual(far.shape[0], 2)
        self.assertAlmostEqual(float(far[0]), 1.0, places=12)
        self.assertAlmostEqual(float(far[1]), 0.5, places=12)


class MixedFarfieldStage2Tests(ScratchCase):
    """Stage-2 parse and tie rows for the derived C3D8+C3D6 fluid shape."""

    def _mixed_merged(self) -> Path:
        far = self.root / "mixed.inp"
        write_mixed_farfield(far)
        output = self.root / "merged.inp"
        assemble.assemble(far, self.wing(), output, self.root / "meta.json",
                          mixed_contract())
        return output

    def test_parse_reads_the_c3d6_block_and_keeps_deck_id_contiguity(self) -> None:
        mesh = modal.parse_merged(self._mixed_merged())
        self.assertEqual(mesh.far_conn.shape, (1, 8))
        self.assertEqual(mesh.wedge_conn.tolist(), [[9, 10, 11, 12, 13, 14]])
        self.assertEqual((mesh.brick_elems, mesh.wedge_elems), (1, 1))
        self.assertEqual((mesh.far_elems, mesh.far_nodes), (2, 14))
        self.assertEqual(mesh.wing_conn.shape, (1, 10))
        surfaces = modal.extract_surfaces(mesh)
        self.assertEqual(surfaces.fluid_rows.tolist(), [0])
        self.assertEqual(surfaces.fluid_flags.tolist(), [1])
        buffer = io.StringIO()
        surfaces.write(buffer, 0.05, mesh.far_elems)
        # Wing C3D10 ids start at the TOTAL far count + 1 (3), not bricks + 1.
        self.assertIn("*SURFACE, NAME=SURF_FSI_WING, TYPE=ELEMENT\n3,S1\n",
                      buffer.getvalue())

    def test_parse_rejects_a_non_contiguous_wedge_element_id(self) -> None:
        merged = self._mixed_merged()
        text = merged.read_text(encoding="ascii")
        self.assertIn("\n2,9,10,11,12,13,14\n", text)
        broken = self.root / "broken.inp"
        broken.write_text(text.replace("\n2,9,10,11,12,13,14\n",
                                       "\n7,9,10,11,12,13,14\n"),
                          encoding="ascii")
        with self.assertRaises(modal.DeckError):
            modal.parse_merged(broken)


class LayerTests(unittest.TestCase):
    def test_layered_config_caps_emax_and_bounds_the_span(self) -> None:
        spec = modal.CaseSpec(layers=5, layer_p=0.25, layer_floor=1000.0)
        volumes = np.geomspace(1e-16, 1.0, 1000)
        config = modal.layered_config(spec, volumes)
        expected = spec.layer_emax_factor * spec.material().E1
        self.assertAlmostEqual(config.E_max, expected, delta=expected * 1e-12)
        self.assertLessEqual(config.contrast, spec.layer_floor * (1.0 + 1e-12))
        self.assertGreaterEqual(config.E_min, config.E_max / spec.layer_floor)
        self.assertEqual(sum(row.count for row in config.table), 1000)
        moduli = [row.E for row in config.table]
        self.assertEqual(moduli, sorted(moduli, reverse=True))

    def test_layered_span_never_collapses_to_the_old_band(self) -> None:
        spec = modal.CaseSpec(layers=20, layer_p=2.5, layer_floor=1000.0)
        volumes = np.geomspace(1e-17, 2.4, 372760)
        config = modal.layered_config(spec, volumes)
        moduli = [row.E for row in config.table]
        self.assertLessEqual(max(moduli) / min(moduli), 1000.0 * (1.0 + 1e-12))
        self.assertGreaterEqual(min(moduli), config.E_max / 1000.0 * (1.0 - 1e-12))
        self.assertLessEqual(max(moduli), spec.layer_emax_factor * spec.material().E1
                             * (1.0 + 1e-12))

    def test_layer_table_floor_marks_only_the_largest_volumes(self) -> None:
        volumes = np.array([1e-16, 1e-8, 1.0, 1e8])
        assignment = np.array([0, 1, 2, 3], dtype=np.int64)
        table = modal.layer_table(volumes, assignment, 0.25, 1000.0, 3518.452897)
        self.assertTrue(table[0].E == 3518.452897)
        self.assertTrue(table[-1].floored)
        self.assertEqual(min(row.E for row in table), 3518.452897 / 1000.0)


class ModalDeckTests(ScratchCase):
    def test_deck_carries_the_verified_card_invariants(self) -> None:
        text = self.combined().read_text(encoding="ascii")
        for invariant in (
                "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
                "3.5184528970e+09",
                "4.1350000000e+02",
                "*ORIENTATION, NAME=WING_AXES\n1,1,0, 0,0,1\n",
                "*SOLID SECTION, ELSET=WING, MATERIAL=soft_anchored, "
                "ORIENTATION=WING_AXES",
                "*FREQUENCY\n10,0.334,1.E6\n",
                "*BOUNDARY\nFIXED_FRONT,1,3\n",
                "ROOT,1,3",
                "*SURFACE, NAME=SURF_FSI_FLUID, TYPE=ELEMENT",
                "*TIE, NAME=FSI_COUPLING, POSITION TOLERANCE=0.05, ADJUST=NO"):
            self.assertIn(invariant, text)
        for absent in ("STORAGE=", "SOLVER=", "TYPE=ORTHO", "*MPC", "8.17e9"):
            self.assertNotIn(absent, text)
        self.assertIn("\n*DENSITY\n0.\n", text)

    def test_hard_wing_uses_the_calibrated_hard_material(self) -> None:
        text = self.combined(wing_type="hard").read_text(encoding="ascii")
        self.assertIn("hard_e1_0p75", text)
        self.assertIn("1.0136052400e+10", text)
        self.assertIn("4.9550000000e+02", text)

    def test_uniform_mode_emits_one_massless_far_field_material(self) -> None:
        deck = self.combined(far_field_mode="uniform", far_field_modulus=28.0)
        text = deck.read_text(encoding="ascii")
        self.assertIn("*MATERIAL, NAME=PFLUID\n*ELASTIC\n2.800000000e+01,0.\n"
                      "*DENSITY\n0.\n", text)
        self.assertNotIn("PFLUID_L", text)
        self.assertEqual(text.count("*MATERIAL"), 2)

    def test_pure_wing_deck_is_built_from_the_wing_source(self) -> None:
        source = self.wing()
        request = modal.BuildRequest(
            merged_deck=self.merged(), output=self.root / "combined.inp",
            meta=self.root / "combined_meta.json",
            spec=modal.CaseSpec(layers=2),
            pure=modal.PureWing(source=source, output=self.root / "pure.inp"))
        report = modal.generate(request)
        text = (self.root / "pure.inp").read_text(encoding="ascii")
        source_text = source.read_text(encoding="ascii")
        copied_prefix = source_text[:source_text.index("*MATERIAL, NAME=WING")]
        self.assertTrue(text.startswith(copied_prefix))
        self.assertIn("*MATERIAL, NAME=soft_anchored", text)
        self.assertIn("*BOUNDARY\nROOT,1,3\n*END STEP\n", text)
        self.assertEqual(report["pure"]["root_nodes"], 1)

    def test_spec_validation_rejects_unbounded_or_unknown_controls(self) -> None:
        for spec in (modal.CaseSpec(modes=0), modal.CaseSpec(layers=0),
                     modal.CaseSpec(layer_floor=1.0),
                     modal.CaseSpec(layer_ref_modulus="E9"),
                     modal.CaseSpec(far_field_mode="shell"),
                     modal.CaseSpec(wing_type="medium"),
                     modal.CaseSpec(far_field_mode="uniform",
                                    far_field_modulus=0.0)):
            with self.assertRaises(modal.DeckError):
                modal.validate_spec(spec)


class ExportGateTests(ScratchCase):
    def _deck(self, aliases: bool = False) -> Path:
        deck = self.root / "case.inp"
        fluid = "FAR_PART" if aliases else "FAR_PART_NSET"
        struct = "WING_PART" if aliases else "WING"
        deck.write_text(
            "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n5,1,1,0\n6,1,0,1\n"
            f"*NSET, NSET={fluid}\n1,2,3\n*NSET, NSET={struct}\n4,5,6\n",
            encoding="ascii")
        return deck

    def _results(self, name: str, combined_freqs: list[float],
                 pure_freqs: list[float]
                 ) -> tuple[Path, Path, Path, Path]:
        combined = [{node: (float(node), float(mode), 0.0) for node in range(1, 7)}
                    for mode in (1, 2)]
        pure = [{node: (float(node + 3), float(mode), 0.0) for node in range(1, 4)}
                for mode in (1, 2)]
        dat = self.root / f"{name}.dat"
        frd = self.root / f"{name}.frd"
        pure_dat = self.root / f"{name}_pure.dat"
        pure_frd = self.root / f"{name}_pure.frd"
        write_dat(dat, combined_freqs)
        write_dat(pure_dat, pure_freqs)
        write_frd(frd, combined, combined_freqs)
        write_frd(pure_frd, pure, pure_freqs)
        return dat, frd, pure_dat, pure_frd

    def _identity_case(self, name: str,
                       combined_fields: list[dict[int, tuple[float, float, float]]],
                       pure_fields: list[dict[int, tuple[float, float, float]]]
                       ) -> tuple[Path, ...]:
        deck = self._deck()
        dat = self.root / f"{name}.dat"
        write_dat(dat, [10.0, 20.0])
        frd = self.root / f"{name}.frd"
        write_frd(frd, combined_fields, [10.0, 20.0])
        pure_frd = self.root / f"{name}_pure.frd"
        write_frd(pure_frd, pure_fields, [10.0, 20.0])
        return deck, dat, frd, dat, pure_frd

    def test_gate_passes_and_exports_the_three_field_contract(self) -> None:
        deck = self._deck()
        dat, frd, pure_dat, pure_frd = self._results("good", [10.0, 20.0],
                                                     [10.0, 20.0])
        output = self.root / "mode"
        result = run_script("03_export_modes.py", "--input", str(deck),
                            "--dat", str(dat), "--frd", str(frd),
                            "--output-dir", str(output), "--modes", "2",
                            "--pure-dat", str(pure_dat),
                            "--pure-frd", str(pure_frd))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("gates passed", result.stdout)
        self.assertEqual((output / "FluidNodeCoor.csv").read_text().splitlines()[0],
                         f"{0.0:20.10E},{3:20.10E},{2:20.10E},")
        self.assertEqual((output / "StructNodeDisp1.csv").read_text().splitlines()[0],
                         f"{10.0:20.10E},{3:20.10E},{2:20.10E},")
        self.assertEqual((output / "FluidPara.csv").read_text(),
                         "19,\n11,\n0.0,\n0.0,\n0.0,\n0.0,\n0.0,\n")

    def test_alias_sets_are_exported_when_merged_names_are_absent(self) -> None:
        deck = self._deck(aliases=True)
        dat, frd, pure_dat, pure_frd = self._results("alias", [10.0, 20.0],
                                                     [10.0, 20.0])
        result = run_script("03_export_modes.py", "--input", str(deck),
                            "--dat", str(dat), "--frd", str(frd),
                            "--output-dir", str(self.root / "mode"),
                            "--modes", "2", "--pure-dat", str(pure_dat),
                            "--pure-frd", str(pure_frd))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.root / "mode" / "StructNodeDisp2.csv").is_file())

    def test_set_lookup_requires_an_exact_nset_name(self) -> None:
        deck = self._deck()
        with self.assertRaises(exporter.ExportError):
            exporter.read_set_ids(deck, "FAR_PART")
        with self.assertRaises(exporter.ExportError):
            exporter.read_set_ids(deck, "WING_PART")

    def test_node_scan_ignores_the_node_file_and_step_cards(self) -> None:
        deck = self._deck()
        deck.write_text(
            deck.read_text(encoding="ascii")
            + "*STEP\n*FREQUENCY\n10,0.334,1.E6\n*NODE FILE\nU\n*END STEP\n",
            encoding="ascii")
        dat, frd, pure_dat, pure_frd = self._results("nodefile", [10.0, 20.0],
                                                     [10.0, 20.0])
        output = self.root / "mode"
        result = run_script("03_export_modes.py", "--input", str(deck),
                            "--dat", str(dat), "--frd", str(frd),
                            "--output-dir", str(output), "--modes", "2",
                            "--pure-dat", str(pure_dat),
                            "--pure-frd", str(pure_frd))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((output / "FluidNodeCoor.csv").is_file())

    def test_gate_aborts_and_writes_nothing_when_modes_are_missing(self) -> None:
        deck = self._deck()
        dat = self.root / "short.dat"
        write_dat(dat, [10.0])
        frd = self.root / "short.frd"
        write_frd(frd, [{node: (1.0, 0.0, 0.0) for node in range(1, 7)}], [10.0])
        output = self.root / "mode"
        result = run_script("03_export_modes.py", "--input", str(deck),
                            "--dat", str(dat), "--frd", str(frd),
                            "--output-dir", str(output), "--modes", "2")
        self.assertEqual(result.returncode, 1)
        self.assertIn("[AGARD ERR]", result.stderr)
        self.assertFalse(output.exists())

    def test_gate_aborts_when_frequency_leaves_the_sanity_window(self) -> None:
        deck = self._deck()
        dat, frd, pure_dat, pure_frd = self._results("window", [10.0, 2.0e4],
                                                     [10.0, 2.0e4])
        output = self.root / "mode"
        result = run_script("03_export_modes.py", "--input", str(deck),
                            "--dat", str(dat), "--frd", str(frd),
                            "--output-dir", str(output), "--modes", "2",
                            "--pure-dat", str(pure_dat),
                            "--pure-frd", str(pure_frd))
        self.assertEqual(result.returncode, 1)
        self.assertFalse(output.exists())

    def test_gate_aborts_when_pure_and_combined_frequencies_diverged(self) -> None:
        deck = self._deck()
        dat, frd, pure_dat, pure_frd = self._results("diverged", [10.0, 20.0],
                                                     [10.0, 25.0])
        output = self.root / "mode"
        result = run_script("03_export_modes.py", "--input", str(deck),
                            "--dat", str(dat), "--frd", str(frd),
                            "--output-dir", str(output), "--modes", "2",
                            "--pure-dat", str(pure_dat),
                            "--pure-frd", str(pure_frd))
        self.assertEqual(result.returncode, 1)
        self.assertIn("max |relative difference|", result.stderr)
        self.assertFalse(output.exists())

    def test_gate_aborts_when_mode_identity_cosine_is_below_threshold(self) -> None:
        wing_one = {4: (1.0, 0.0, 0.0), 5: (1.0, 0.0, 0.0), 6: (1.0, 0.0, 0.0)}
        wing_two = {4: (0.0, 1.0, 0.0), 5: (0.0, 1.0, 0.0), 6: (0.0, 1.0, 0.0)}
        pure_two = {1: (0.0, 0.5, 0.866), 2: (0.0, 0.5, 0.866),
                    3: (0.0, 0.5, 0.866)}
        deck, dat, frd, pure_dat, pure_frd = self._identity_case(
            "cosine", [wing_one, wing_two],
            [{1: (1.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 0.0, 0.0)},
             pure_two])
        output = self.root / "mode"
        result = run_script("03_export_modes.py", "--input", str(deck),
                            "--dat", str(dat), "--frd", str(frd),
                            "--output-dir", str(output), "--modes", "2",
                            "--pure-dat", str(pure_dat),
                            "--pure-frd", str(pure_frd))
        self.assertEqual(result.returncode, 1)
        self.assertIn("field cosine", result.stderr)
        self.assertFalse(output.exists())

    def test_gate_aborts_when_mode_identity_is_off_diagonal(self) -> None:
        wing_one = {4: (1.0, 0.0, 0.0), 5: (1.0, 0.0, 0.0), 6: (1.0, 0.0, 0.0)}
        wing_two = {4: (0.0, 1.0, 0.0), 5: (0.0, 1.0, 0.0), 6: (0.0, 1.0, 0.0)}
        pure_two = {1: (0.0, 0.0, 1.0), 2: (0.0, 0.0, 1.0), 3: (0.0, 0.0, 1.0)}
        deck, dat, frd, pure_dat, pure_frd = self._identity_case(
            "diagonal", [wing_one, wing_two],
            [{1: (1.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 0.0, 0.0)},
             pure_two])
        output = self.root / "mode"
        result = run_script("03_export_modes.py", "--input", str(deck),
                            "--dat", str(dat), "--frd", str(frd),
                            "--output-dir", str(output), "--modes", "2",
                            "--pure-dat", str(pure_dat),
                            "--pure-frd", str(pure_frd))
        self.assertEqual(result.returncode, 1)
        self.assertIn("best-matches", result.stderr)
        self.assertFalse(output.exists())

    def test_gate_window_defaults_match_the_documented_hz_band(self) -> None:
        frequencies = exporter.gate_frequencies(
            self._write_single_dat([9.6, 311.5]), 2, (1.0, 1.0e4))
        self.assertEqual(frequencies, [9.6, 311.5])
        with self.assertRaises(exporter.GateError):
            exporter.gate_frequencies(self._write_single_dat([0.5, 311.5]), 2,
                                      (1.0, 1.0e4))
        with self.assertRaises(exporter.GateError):
            exporter.gate_frequencies(self._write_single_dat([100.0, 50.0]), 2,
                                      (1.0, 1.0e4))

    def _write_single_dat(self, frequencies: list[float]) -> Path:
        path = self.root / f"single_{frequencies[0]}.dat"
        write_dat(path, frequencies)
        return path


class StressExportTests(ScratchCase):
    """T0-pinned STRESS/TOSTRAIN blocks: parsing, export and flag-off no-op."""

    STRESS_NODE4_MODE1 = (-9.27137e9, -3.97347e9, -3.97348e9, -1.64314e8,
                          115.574, -1.86096e9)
    STRESS_NODE4_MODE2 = (-2.09977e9, -8.99908e8, -8.99909e8,
                          7.25513e8, -510.309, -4.21468e8)
    STRAIN_NODE4_MODE1 = (-3.27966e-2, -7.24886e-8, -1.11033e-7,
                          -1.01718e-3, 7.15461e-10, -1.15202e-2)
    LEGACY_FILES = sorted(["FluidNodeCoor.csv", "FluidNodeDisp1.csv",
                           "FluidNodeDisp2.csv", "FluidPara.csv",
                           "StructNodeCoor.csv", "StructNodeDisp1.csv",
                           "StructNodeDisp2.csv"])
    STRESS_FILES = sorted(LEGACY_FILES + ["StructNodeStrain1.csv",
                                          "StructNodeStrain2.csv",
                                          "StructNodeStress1.csv",
                                          "StructNodeStress2.csv"])

    def _deck(self) -> Path:
        deck = self.root / "stress_case.inp"
        deck.write_text(
            "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n5,1,1,0\n6,1,0,1\n"
            "*NSET, NSET=FAR_PART_NSET\n1,2,3\n*NSET, NSET=WING\n4,5,6\n",
            encoding="ascii")
        return deck

    def _displacements(self) -> list[dict[int, tuple[float, float, float]]]:
        return [{node: (float(node), float(mode), 0.0) for node in range(1, 7)}
                for mode in (1, 2)]

    def _run(self, frd: Path, output: Path,
             *extra: str) -> subprocess.CompletedProcess[str]:
        dat = self.root / "stress.dat"
        write_dat(dat, [10.0, 20.0])
        return run_script("03_export_modes.py", "--input", str(self._deck()),
                          "--dat", str(dat), "--frd", str(frd),
                          "--output-dir", str(output), "--modes", "2", *extra)

    def _frd_with_stress(self) -> Path:
        frd = self.root / "stress.frd"
        write_frd_with_stress(frd, self._displacements(),
                              (T0_STRESS_MODE1, T0_STRESS_MODE2),
                              (T0_STRAIN_MODE1, T0_STRAIN_MODE2),
                              [10.0, 20.0])
        return frd

    @staticmethod
    def _row_values(line: str) -> list[float]:
        return [float(value) for value in line.rstrip(",").split(",")]

    def test_export_writes_structural_rows_in_captured_component_order(self) -> None:
        output = self.root / "mode"
        result = self._run(self._frd_with_stress(), output, "--export-stress")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sorted(path.name for path in output.iterdir()),
                         self.STRESS_FILES)
        self.assertFalse((output / "FluidNodeStress1.csv").exists())
        stress = (output / "StructNodeStress1.csv").read_text(
            encoding="ascii").splitlines()
        self.assertEqual(stress[0], f"{10.0:20.10E},{3:20.10E},{2:20.10E},")
        self.assertEqual(len(stress), 1 + 3)
        self.assertEqual(self._row_values(stress[1]),
                         list(self.STRESS_NODE4_MODE1))
        strain = (output / "StructNodeStrain1.csv").read_text(
            encoding="ascii").splitlines()
        self.assertEqual(self._row_values(strain[1]),
                         list(self.STRAIN_NODE4_MODE1))
        second = (output / "StructNodeStress2.csv").read_text(
            encoding="ascii").splitlines()
        self.assertEqual(second[0], f"{20.0:20.10E},{3:20.10E},{2:20.10E},")
        self.assertEqual(self._row_values(second[1]),
                         list(self.STRESS_NODE4_MODE2))

    def test_gate_parser_reads_six_components_and_skips_error_blocks(self) -> None:
        frd = self._frd_with_stress()
        stress = exporter.parse_frd_modes(frd, "STRESS")
        strain = exporter.parse_frd_modes(frd, "TOSTRAIN")
        disp = exporter.parse_frd_modes(frd)
        self.assertEqual([len(stress), len(strain), len(disp)], [2, 2, 2])
        self.assertEqual(stress[0].vectors.shape, (6, 6))
        self.assertEqual(strain[1].vectors.shape, (6, 6))
        self.assertEqual(disp[0].vectors.shape, (6, 3))
        np.testing.assert_array_equal(stress[0].nodes, np.arange(1, 7))
        np.testing.assert_allclose(stress[0].vectors[3],
                                   self.STRESS_NODE4_MODE1, rtol=1e-12)
        np.testing.assert_allclose(strain[0].vectors[3],
                                   self.STRAIN_NODE4_MODE1, rtol=1e-12)
        self.assertFalse(np.any(stress[0].vectors == 999.0))
        self.assertFalse(np.any(strain[1].vectors == 999.0))

    def test_flag_off_is_a_no_op_for_the_legacy_disp_contract(self) -> None:
        frd = self.root / "legacy.frd"
        write_frd(frd, self._displacements(), [10.0, 20.0])
        output = self.root / "mode"
        result = self._run(frd, output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sorted(path.name for path in output.iterdir()),
                         self.LEGACY_FILES)
        rows = "".join(f"{node:20.10E},{1:20.10E},{0:20.10E},\n"
                       for node in (4, 5, 6))
        self.assertEqual((output / "StructNodeDisp1.csv").read_text(
            encoding="ascii"), f"{10.0:20.10E},{3:20.10E},{2:20.10E},\n" + rows)

    def test_a_minus_two_continuation_is_rejected_with_its_line(self) -> None:
        frd = self._frd_with_stress()
        source = frd.read_text(encoding="ascii")
        self.assertEqual(source.count(T0_STRESS_MODE2[0]), 1)
        broken = self.root / "continuation.frd"
        broken.write_text(source.replace(T0_STRESS_MODE2[0],
                                         T0_STRESS_MODE2[0] + " -2extra\n", 1),
                          encoding="ascii")
        with self.assertRaises(exporter.GateError) as caught:
            exporter.parse_frd_modes(broken, "STRESS")
        self.assertRegex(str(caught.exception), r":\d+: STRESS record")
        self.assertIn("-2 continuation", str(caught.exception))
        result = self._run(broken, self.root / "mode", "--export-stress")
        self.assertEqual(result.returncode, 1)
        self.assertIn("[AGARD ERR]", result.stderr)
        self.assertRegex(result.stderr, r"line \d+")

    def test_an_unterminated_block_is_rejected_with_its_line(self) -> None:
        frd = self._frd_with_stress()
        source = frd.read_text(encoding="ascii")
        marker = T0_STRESS_MODE1[-1] + " -3\n"
        self.assertEqual(source.count(marker), 1)
        broken = self.root / "unterminated.frd"
        broken.write_text(source.replace(marker, T0_STRESS_MODE1[-1], 1),
                          encoding="ascii")
        with self.assertRaises(exporter.GateError) as caught:
            exporter.parse_frd_modes(broken, "STRESS")
        self.assertRegex(str(caught.exception), r":\d+: unterminated -4 STRESS")


class PipelineTests(ScratchCase):
    def _pipeline(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(["bash", str(MODE_DIR / "run_pipeline.sh"),
                               *arguments], text=True, capture_output=True,
                              check=False)

    def test_dry_run_orders_stages_forwards_controls_and_writes_nothing(self) -> None:
        work_root = self.root / "work"
        result = self._pipeline("--work-root", str(work_root), "--steps", "3,2,1,2",
                                "--run-ccx",
                                "--modes", "7", "--layers", "12",
                                "--tie-tolerance", "0.2", "--wing-type", "hard",
                                "--far-field-mode", "uniform",
                                "--far-field-modulus", "31",
                                "--pure-wing", "--max-rel", "0.02",
                                "--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertLess(result.stdout.index("stage 1:"), result.stdout.index("stage 2:"))
        self.assertIn("01_assemble_mesh.py", result.stdout)
        self.assertIn("mesh/farfield_c3d8.zip", result.stdout)
        self.assertIn("mesh/wing_c3d10.zip", result.stdout)
        self.assertIn("--modes 7", result.stdout)
        self.assertIn("--wing-type hard", result.stdout)
        self.assertIn("--far-field-mode uniform --far-field-modulus 31", result.stdout)
        self.assertIn("--pure-wing", result.stdout)
        self.assertIn("--max-rel 0.02", result.stdout)
        self.assertIn("wing-type=hard", result.stdout)
        self.assertFalse(work_root.exists())

    def test_dry_run_forwards_layered_and_solver_controls(self) -> None:
        work_root = self.root / "work"
        result = self._pipeline("--work-root", str(work_root), "--steps", "1,2,3",
                                "--run-ccx", "--pure-wing", "--ccx-command", "ccx_mp",
                                "--layers", "18", "--layer-emax-factor", "1e-7",
                                "--layer-ref-modulus", "E2", "--layer-p", "0.3",
                                "--layer-floor", "500", "--min-cosine", "0.95",
                                "--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--layers 18 --layer-emax-factor 1e-7 --layer-ref-modulus E2 "
                      "--layer-p 0.3 --layer-floor 500", result.stdout)
        self.assertIn("--run-pure-ccx", result.stdout)
        self.assertIn("--ccx-command 'ccx_mp'", result.stdout)
        self.assertIn("--min-cosine 0.95", result.stdout)
        self.assertIn("--pure-dat", result.stdout)
        self.assertFalse(work_root.exists())

    def test_stage3_results_dir_route_is_documented_and_forwarded(self) -> None:
        result = self._pipeline("--steps", "3",
                                "--results-dir",
                                "/root/Workspace/AGARD/Mode-run/modal",
                                "--results-stem", "combined_soft_layered",
                                "--modes", "7", "--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--input '/root/Workspace/AGARD/Mode-run/modal/"
                      "combined_soft_layered.inp'", result.stdout)
        self.assertIn("--modes 7", result.stdout)

    def test_dry_run_fluid_profile_needs_the_explicit_deck_and_gains_a_stem(self) -> None:
        work_root = self.root / "work"
        deck = self.root / "fluid8.inp"
        result = self._pipeline("--work-root", str(work_root),
                                "--farfield-profile", "fluid",
                                "--farfield-deck", str(deck), "--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--farfield-profile fluid", result.stdout)
        self.assertIn(f"--farfield-deck '{deck}'", result.stdout)
        self.assertIn("combined_soft_layered_fluid.inp", result.stdout)
        self.assertIn("combined_soft_layered_fluid_meta.json", result.stdout)
        self.assertIn("farfield-profile=fluid", result.stdout)
        self.assertNotIn("mesh/farfield_c3d8.zip", result.stdout)
        self.assertFalse(work_root.exists())

    def test_fluid_profile_without_a_deck_fails_before_running(self) -> None:
        result = self._pipeline("--farfield-profile", "fluid", "--dry-run")
        self.assertEqual(result.returncode, 2)
        self.assertIn("[AGARD ERR]", result.stderr)
        self.assertIn("--farfield-deck", result.stderr)

    def test_dry_run_explicit_deck_wins_over_the_legacy_default(self) -> None:
        work_root = self.root / "work"
        deck = self.root / "custom_far.zip"
        result = self._pipeline("--work-root", str(work_root),
                                "--farfield-profile", "legacy",
                                "--farfield-deck", str(deck), "--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f"--farfield-deck '{deck}'", result.stdout)
        self.assertNotIn("mesh/farfield_c3d8.zip", result.stdout)
        self.assertIn("combined_soft_layered.inp", result.stdout)

    def test_legacy_dry_run_keeps_the_signed_off_command_text(self) -> None:
        work_root = self.root / "work"
        default = self._pipeline("--work-root", str(work_root), "--dry-run")
        explicit = self._pipeline("--work-root", str(work_root),
                                  "--farfield-profile", "legacy", "--dry-run")
        self.assertEqual(default.returncode, 0, default.stderr)
        self.assertEqual(explicit.returncode, 0, explicit.stderr)
        self.assertEqual(default.stdout, explicit.stdout)
        # The legacy stage-1 line is pinned byte-for-byte: the profile flag is
        # only added for a non-legacy profile, and the deck default is the
        # signed-off farfield_c3d8.zip.
        self.assertIn(
            f"[AGARD-PIPE] stage 1: python3 '{MODE_DIR}/01_assemble_mesh.py'"
            f" --farfield-deck '{MODE_DIR}/mesh/farfield_c3d8.zip'"
            f" --wing-deck '{MODE_DIR}/mesh/wing_c3d10.zip'"
            f" --output '{work_root}/mesh/merged_base.inp'"
            f" --meta '{work_root}/mesh/merged_meta.json'"
            f" --tie-tolerance 0.05\n", default.stdout)
        self.assertNotIn("--farfield-profile", default.stdout)
        self.assertIn("combined_soft_layered.inp", default.stdout)
        self.assertNotIn("combined_soft_layered_fluid7", default.stdout)

    def test_help_documents_the_farfield_profiles_and_deck_defaults(self) -> None:
        result = self._pipeline("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("--farfield-profile", result.stdout)
        self.assertIn("fluid", result.stdout)
        self.assertIn("--farfield-deck", result.stdout)
        self.assertIn("mesh/farfield_c3d8.zip", result.stdout)
        self.assertIn("_fluid", result.stdout)
        self.assertNotIn("fluid7", result.stdout)

    def test_pipeline_rejects_invalid_controls_with_exit_code_two(self) -> None:
        for arguments in (("--steps", "4"), ("--modes", "0"), ("--layers", "0"),
                          ("--tie-tolerance", "0"), ("--wing-type", "medium"),
                          ("--far-field-mode", "shell"),
                          ("--farfield-profile", "mystery"),
                          ("--layer-ref-modulus", "E9"),
                          ("--layer-emax-factor", "0"),
                          ("--layer-floor", "1"),
                          ("--min-cosine", "1.5"),
                          ("--max-rel", "0"),
                          ("--steps", "1", "--run-ccx"),
                          ("--steps", "3"),
                          ("--steps", "1,3", "--run-ccx"),
                          ("--steps", "1", "--pure-wing")):
            result = self._pipeline(*arguments, "--dry-run")
            self.assertEqual(result.returncode, 2, arguments)
            self.assertIn("[AGARD ERR]", result.stderr)

    def test_help_documents_the_pure_wing_reference_and_far_field_rule(self) -> None:
        result = self._pipeline("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("--pure-wing", result.stdout)
        self.assertIn("pure_<type>.inp", result.stdout)
        self.assertIn("bounded band", result.stdout)


class RemovedModuleTests(unittest.TestCase):
    def test_no_reference_to_the_deleted_icem_and_layer_modules(self) -> None:
        removed = ("01_convert" + "_icem", "pseudoelastic" + "_layers")
        offenders: list[str] = []
        for path in MODE_DIR.rglob("*"):
            if path.is_dir() or path.name == Path(__file__).name:
                continue
            if path.suffix not in (".py", ".sh", ".md"):
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            if any(name in text for name in removed):
                offenders.append(str(path.relative_to(MODE_DIR)))
        self.assertEqual(offenders, [])


@unittest.skipUnless(os.environ.get("AGARD_RUN_SLOW_TESTS") == "1",
                     "set AGARD_RUN_SLOW_TESTS=1 to check against the "
                     "reference archive")
class ReferenceEquivalenceTests(unittest.TestCase):
    """Byte-level comparison with the signed-off reference decks (read-only)."""

    @classmethod
    def setUpClass(cls) -> None:
        """Resolve the reference lazily so the fast suite extracts nothing."""
        global REFERENCE_DIR
        override = os.environ.get("AGARD_REFERENCE_DIR")
        if override:
            REFERENCE_DIR = Path(override)
            return
        if not REFERENCE_ARCHIVE.is_file():
            return  # setUp() reports the missing reference archive as a skip
        cache = SCRATCH / "reference"
        marker = cache / ".archive.sha256"
        digest = sha256(REFERENCE_ARCHIVE)
        if not (marker.is_file()
                and marker.read_text(encoding="ascii").strip() == digest):
            cache.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(REFERENCE_ARCHIVE) as archive:
                archive.extractall(cache)
            marker.write_text(digest + "\n", encoding="ascii")
        REFERENCE_DIR = cache

    def setUp(self) -> None:
        if REFERENCE_DIR is None or not (
                REFERENCE_DIR / "soft" / "combined_soft_layered.inp").is_file():
            self.skipTest(
                f"reference archive is not available: {REFERENCE_ARCHIVE}")
        SCRATCH.mkdir(parents=True, exist_ok=True)
        self._temp = tempfile.TemporaryDirectory(dir=SCRATCH)
        self.root = Path(self._temp.name)

    def tearDown(self) -> None:
        self._temp.cleanup()

    def _run(self, script: str, *arguments: str) -> None:
        result = run_script(script, *arguments)
        if result.returncode != 0:
            self.fail(f"{script} failed: {result.stderr}")

    def test_soft_layered_deck_matches_the_verified_reference(self) -> None:
        merged = self.root / "merged_base.inp"
        self._run("01_assemble_mesh.py",
                  "--output", str(merged),
                  "--meta", str(self.root / "merged_meta.json"))
        deck = self.root / "combined_soft_layered.inp"
        self._run("02_generate_modal_case.py",
                  "--input", str(merged), "--output", str(deck),
                  "--meta", str(self.root / "combined_meta.json"),
                  "--wing-type", "soft", "--far-field-mode", "layered")
        reference = REFERENCE_DIR / "soft" / "combined_soft_layered.inp"
        reference_bytes = reference.read_bytes()
        deck_bytes = deck.read_bytes()
        # Structural content is the interchange contract.  Compare the
        # *NODE..EOF data-line projection byte-for-byte: comment lines are
        # dropped because the material-slot marker legitimately names its real
        # producer now, and any single differing data byte must still fail.
        self.assertEqual(data_line_projection(deck_bytes),
                         data_line_projection(reference_bytes))
        self.assertIn(b"(rewritten by case/AGARD/Mode/02_generate_modal_case.py)",
                      deck_bytes)
        self.assertNotIn(b"make_cases.py", deck_bytes)
        self.assertNotIn(b"build_merged", deck_bytes)
        header = header_block(deck_bytes)
        self.assertNotIn(b"build_merged", header)
        self.assertNotIn(b"Mode3", header)
        self.assertIn(b"01_assemble_mesh", header)

    def test_stage_one_counts_match_the_verified_build_report(self) -> None:
        merged = self.root / "merged_base.inp"
        meta = self.root / "merged_meta.json"
        self._run("01_assemble_mesh.py", "--output", str(merged),
                  "--meta", str(meta))
        report = json.loads(meta.read_text())
        self.assertEqual(report["counts"], {"nodes": 1365551, "c3d8": 372760,
                                            "c3d10": 532373, "elements": 905133})
        self.assertEqual(report["far"]["fsi_nodes"], 4045)
        self.assertAlmostEqual(report["fsi_wing_distance"]["max_m"],
                               0.018314251661946102, delta=1e-12)
        self.assertLess(report["wing"]["volume_rel_err"], 1e-6)
        reference = REFERENCE_DIR / "mesh" / "merged_base.inp"
        merged_bytes = merged.read_bytes()
        merged_content = data_line_projection(merged_bytes)
        reference_content = data_line_projection(reference.read_bytes())
        # Stage 1 is mesh-only: its *NODE..EOF data-line projection is a strict
        # byte prefix of the signed-off merged deck's projection (the reference
        # appends the material/step tail).  Only the header block and comment
        # text differ.
        self.assertLess(len(merged_content), len(reference_content))
        self.assertEqual(merged_content, reference_content[:len(merged_content)])
        header = header_block(merged_bytes)
        self.assertNotIn(b"build_merged", header)
        self.assertNotIn(b"Mode3", header)
        self.assertIn(b"01_assemble_mesh", header)

    def test_uniform_and_pure_decks_match_the_verified_references(self) -> None:
        merged = self.root / "merged_base.inp"
        self._run("01_assemble_mesh.py", "--output", str(merged),
                  "--meta", str(self.root / "merged_meta.json"))
        uniform = self.root / "combined_soft_uniform.inp"
        pure = self.root / "pure_soft.inp"
        self._run("02_generate_modal_case.py",
                  "--input", str(merged), "--output", str(uniform),
                  "--meta", str(self.root / "uniform_meta.json"),
                  "--wing-type", "soft", "--far-field-mode", "uniform",
                  "--far-field-modulus", "28", "--pure-wing",
                  "--pure-output", str(pure),
                  "--wing-source", str(MODE_DIR / "mesh" / "wing_c3d10.zip"))
        self.assertEqual(data_line_projection(uniform.read_bytes()),
                         data_line_projection(
                             (REFERENCE_DIR / "soft" / "combined_soft_fix.inp").read_bytes()))
        self.assertNotIn(b"make_cases.py", uniform.read_bytes())
        self.assertNotIn(b"build_merged", uniform.read_bytes())
        uniform_header = header_block(uniform.read_bytes())
        self.assertNotIn(b"build_merged", uniform_header)
        self.assertNotIn(b"Mode3", uniform_header)
        self.assertIn(b"01_assemble_mesh", uniform_header)
        self.assertNotIn(b"build_merged", pure.read_bytes())
        # The pure-wing deck has no stale producer text: keep whole-file sha256.
        self.assertEqual(sha256(pure),
                         sha256(REFERENCE_DIR / "soft" / "pure_soft.inp"))


if __name__ == "__main__":
    unittest.main()
