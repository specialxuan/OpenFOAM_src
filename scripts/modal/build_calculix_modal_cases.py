#!/usr/bin/env python3
"""
Build CalculiX modal input decks from ANSYS CDB-like DamFailure.in files.

Workflow:
1. Parse one solid reference CDB (structure mesh and structure node sets).
2. Parse one or more fluid CDB files (fluid mesh and fluid node sets).
3. Combine each fluid mesh with the same solid mesh, then write CalculiX .inp.

This script mirrors the APDL logic in ModeCalculation.apdl:
- pseudo-fluid + solid coupled modal model
- boundary conditions on TOP/BOTTOM/LEFT/RIGHT and FRONT/BACK/SIDE1/SIDE2
- interface compatibility between FSI_FLUID and FSI_SOLID
- pseudo-fluid stiffness scaled by element volume ratio (minV / V)^power
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


FLOAT_RE = re.compile(r"[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[EeDd][+-]?\d+)?")
INT_RE = re.compile(r"[+-]?\d+")
CMBLOCK_RE = re.compile(r"^CMBLOCK,([^,]+),NODE,(\d+)", re.IGNORECASE)


@dataclass(frozen=True)
class Element:
    eid: int
    etype: int
    mat: int
    conn: Tuple[int, int, int, int, int, int, int, int]


@dataclass
class CDBModel:
    nodes: Dict[int, Tuple[float, float, float]]
    elements: List[Element]
    node_sets: Dict[str, List[int]]


def parse_cdb(path: Path) -> CDBModel:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    nodes: Dict[int, Tuple[float, float, float]] = {}
    elements: List[Element] = []
    node_sets: Dict[str, List[int]] = {}

    i = 0
    n = len(lines)
    while i < n:
        raw = lines[i]
        line = raw.strip()
        uline = line.upper()

        if uline.startswith("NBLOCK"):
            i = _parse_nblock(lines, i + 1, nodes)
            continue

        if uline.startswith("EBLOCK"):
            i = _parse_eblock(lines, i + 1, elements)
            continue

        m = CMBLOCK_RE.match(line)
        if m:
            name = m.group(1).strip().upper()
            expected = int(m.group(2))
            i = _parse_cmblock(lines, i + 1, expected, name, node_sets)
            continue

        i += 1

    if not nodes:
        raise ValueError(f"No nodes found in {path}")
    if not elements:
        raise ValueError(f"No elements found in {path}")

    return CDBModel(nodes=nodes, elements=elements, node_sets=node_sets)


def _parse_nblock(lines: Sequence[str], start: int, out_nodes: Dict[int, Tuple[float, float, float]]) -> int:
    i = start
    n = len(lines)

    if i < n and lines[i].strip().startswith("("):
        i += 1

    while i < n:
        line = lines[i]
        s = line.strip()
        us = s.upper()
        if not s:
            i += 1
            continue
        if s.startswith("-1") or us.startswith("N ,") or us.startswith("N,"):
            return i + 1

        head = line[:24]
        ints = INT_RE.findall(head)
        floats = [float(v.replace("D", "E").replace("d", "e")) for v in FLOAT_RE.findall(line[24:])]
        if ints and len(floats) >= 3:
            nid = int(ints[0])
            out_nodes[nid] = (floats[0], floats[1], floats[2])
        i += 1

    return i


def _parse_eblock(lines: Sequence[str], start: int, out_elements: List[Element]) -> int:
    i = start
    n = len(lines)

    if i < n and lines[i].strip().startswith("("):
        i += 1

    while i < n:
        s = lines[i].strip()
        us = s.upper()
        if not s:
            i += 1
            continue
        if s.startswith("-1") or us.startswith("EN ,") or us.startswith("EN,"):
            return i + 1

        vals = [int(v) for v in INT_RE.findall(s)]
        if len(vals) >= 19:
            etype = vals[0]
            mat = vals[1]
            eid = vals[10]
            conn = tuple(vals[11:19])  # C3D8 connectivity in this model
            out_elements.append(Element(eid=eid, etype=etype, mat=mat, conn=conn))  # type: ignore[arg-type]
        i += 1

    return i


def _parse_cmblock(
    lines: Sequence[str],
    start: int,
    expected_count: int,
    name: str,
    out_sets: Dict[str, List[int]],
) -> int:
    i = start
    n = len(lines)

    if i < n and lines[i].strip().startswith("("):
        i += 1

    vals: List[int] = []
    while i < n and len(vals) < expected_count:
        s = lines[i].strip()
        if s:
            vals.extend(int(v) for v in INT_RE.findall(s))
        i += 1

    out_sets[name] = vals[:expected_count]
    return i


def union_sorted(*iterables: Iterable[int]) -> List[int]:
    merged = set()
    for it in iterables:
        merged.update(it)
    return sorted(merged)


def vec_sub(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vec_cross(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vec_dot(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def tet_volume(
    a: Tuple[float, float, float],
    b: Tuple[float, float, float],
    c: Tuple[float, float, float],
    d: Tuple[float, float, float],
) -> float:
    ad = vec_sub(a, d)
    bd = vec_sub(b, d)
    cd = vec_sub(c, d)
    return abs(vec_dot(ad, vec_cross(bd, cd))) / 6.0


def tet_volume_signed(
    a: Tuple[float, float, float],
    b: Tuple[float, float, float],
    c: Tuple[float, float, float],
    d: Tuple[float, float, float],
) -> float:
    ad = vec_sub(a, d)
    bd = vec_sub(b, d)
    cd = vec_sub(c, d)
    return vec_dot(ad, vec_cross(bd, cd)) / 6.0


def hex_volume(coords: Sequence[Tuple[float, float, float]]) -> float:
    if len(coords) != 8:
        raise ValueError("Hex element must have 8 coordinates")

    # 5-tet decomposition for convex HEX8
    tets = [
        (0, 1, 3, 4),
        (1, 2, 3, 6),
        (1, 3, 4, 6),
        (1, 4, 5, 6),
        (3, 4, 6, 7),
    ]
    return sum(tet_volume(coords[i], coords[j], coords[k], coords[l]) for i, j, k, l in tets)


def dist2(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return dx * dx + dy * dy + dz * dz


def remap_ids(ids: Iterable[int], mapping: Dict[int, int]) -> List[int]:
    out: List[int] = []
    for old in ids:
        new = mapping.get(old)
        if new is not None:
            out.append(new)
    return sorted(set(out))


def write_id_chunks(fp, ids: Sequence[int], per_line: int = 16) -> None:
    for i in range(0, len(ids), per_line):
        chunk = ids[i : i + per_line]
        fp.write(", ".join(str(v) for v in chunk) + "\n")


def write_nset(fp, name: str, ids: Sequence[int]) -> None:
    fp.write(f"*NSET, NSET={name}\n")
    if ids:
        write_id_chunks(fp, ids, per_line=16)
    else:
        fp.write("** empty set\n")


def write_elset(fp, name: str, ids: Sequence[int]) -> None:
    fp.write(f"*ELSET, ELSET={name}\n")
    if ids:
        write_id_chunks(fp, ids, per_line=16)
    else:
        fp.write("** empty set\n")


def write_equation(fp, terms: Sequence[Tuple[int, int, float]]) -> None:
    fp.write("*EQUATION\n")
    fp.write(f"{len(terms)}\n")
    flat: List[str] = []
    for nid, dof, coef in terms:
        flat.extend([str(nid), str(dof), f"{coef:.14e}"])
    for i in range(0, len(flat), 9):
        fp.write(", ".join(flat[i : i + 9]) + "\n")


def build_case(
    fluid_model: CDBModel,
    solid_model: CDBModel,
    fluid_cdb_path: Path,
    solid_cdb_path: Path,
    case_name: str,
    out_dir: Path,
    n_modes: int,
    fluid_base_e: float,
    fluid_stiff_power: float,
    fluid_stiff_digits: int,
    ceintf_tol: float,
    interp_k: int,
    solid_e: float = 1.0e6,
    solid_nu: float = 0.0,
    solid_rho: float = 2500.0,
    tet_volumes: Optional[Dict[int, float]] = None,
) -> Path:
    fluid_set_names = [
        "FLUIDDOMAIN",
        "FSI_FLUID",
        "FRONT",
        "BACK",
        "TOP",
        "BOTTOM_FLUID",
        "LEFT",
        "RIGHT",
        "TIE_MASTER",
        "TIE_SLAVE",
    ]
    solid_set_names = [
        "SOLID",
        "FSI_SOLID",
        "BOTTOM_SOLID",
        "SIDE1",
        "SIDE2",
    ]

    fluid_elems_src = [e for e in fluid_model.elements if e.etype == 1]
    solid_elems_src = [e for e in solid_model.elements if e.etype == 2]
    if not fluid_elems_src:
        raise ValueError(f"{case_name}: no fluid elements (etype=1) found")
    if not solid_elems_src:
        raise ValueError(f"{case_name}: no solid elements (etype=2) found")

    fluid_nodes_needed = set()
    for e in fluid_elems_src:
        fluid_nodes_needed.update(e.conn)
    for sname in fluid_set_names:
        fluid_nodes_needed.update(fluid_model.node_sets.get(sname, []))

    solid_nodes_needed = set()
    for e in solid_elems_src:
        solid_nodes_needed.update(e.conn)
    for sname in solid_set_names:
        solid_nodes_needed.update(solid_model.node_sets.get(sname, []))

    fluid_nodes_sorted = sorted(fluid_nodes_needed)
    solid_nodes_sorted = sorted(solid_nodes_needed)

    fluid_map = {old: i + 1 for i, old in enumerate(fluid_nodes_sorted)}
    solid_start = len(fluid_nodes_sorted) + 1
    solid_map = {old: solid_start + i for i, old in enumerate(solid_nodes_sorted)}

    new_nodes: Dict[int, Tuple[float, float, float]] = {}
    for old, new in fluid_map.items():
        if old not in fluid_model.nodes:
            raise ValueError(f"{case_name}: missing fluid node {old}")
        new_nodes[new] = fluid_model.nodes[old]
    for old, new in solid_map.items():
        if old not in solid_model.nodes:
            raise ValueError(f"{case_name}: missing solid node {old}")
        new_nodes[new] = solid_model.nodes[old]

    new_fluid_elems: List[Tuple[int, Tuple[int, ...]]] = []
    new_solid_elems: List[Tuple[int, Tuple[int, ...]]] = []
    next_eid = 1
    for e in fluid_elems_src:
        conn = tuple(fluid_map[n] for n in e.conn)
        new_fluid_elems.append((next_eid, conn))
        next_eid += 1
    for e in solid_elems_src:
        conn = tuple(solid_map[n] for n in e.conn)
        new_solid_elems.append((next_eid, conn))
        next_eid += 1

    # Split fluid elements: hex (5+ unique nodes) vs tet (<=4 unique)
    fluid_hex_elems = []
    fluid_tet_elems = []
    for eid, conn in new_fluid_elems:
        if len(set(conn)) <= 4:
            fluid_tet_elems.append((eid, conn))
        else:
            fluid_hex_elems.append((eid, conn))

    fluid_sets = {name: remap_ids(fluid_model.node_sets.get(name, []), fluid_map) for name in fluid_set_names}
    solid_sets = {name: remap_ids(solid_model.node_sets.get(name, []), solid_map) for name in solid_set_names}

    nset_fluiddomain = fluid_sets.get("FLUIDDOMAIN", [])
    nset_solid = solid_sets.get("SOLID", [])
    nset_fsi_fluid = fluid_sets.get("FSI_FLUID", [])
    nset_fsi_solid = solid_sets.get("FSI_SOLID", [])
    nset_fastdynamic = union_sorted(nset_fluiddomain, nset_solid)

    nset_fix_all = union_sorted(
        fluid_sets.get("TOP", []),
        fluid_sets.get("BOTTOM_FLUID", []),
        solid_sets.get("BOTTOM_SOLID", []),
        fluid_sets.get("LEFT", []),
        fluid_sets.get("RIGHT", []),
    )
    nset_fix_uz = union_sorted(
        fluid_sets.get("FRONT", []),
        fluid_sets.get("BACK", []),
        solid_sets.get("SIDE1", []),
        solid_sets.get("SIDE2", []),
    )

    constrained_fluid_nodes = set(nset_fix_all) | set(nset_fix_uz)
    nset_fsi_fluid_coupled = [nid for nid in nset_fsi_fluid if nid not in constrained_fluid_nodes]

    # Interface equations (CEINTF-like nearest interpolation)
    if not nset_fsi_fluid_coupled or not nset_fsi_solid:
        raise ValueError(f"{case_name}: FSI node sets are empty")

    solid_fsi_coords = [(nid, new_nodes[nid]) for nid in nset_fsi_solid]
    interface_map: Dict[int, List[Tuple[int, float]]] = {}
    tol2 = ceintf_tol * ceintf_tol

    for fnid in nset_fsi_fluid_coupled:
        p = new_nodes[fnid]
        ranked = sorted((dist2(p, q), snid) for snid, q in solid_fsi_coords)
        if not ranked or ranked[0][0] > tol2:
            d = math.sqrt(ranked[0][0]) if ranked else float("inf")
            raise ValueError(
                f"{case_name}: no solid interface node within tolerance for fluid node {fnid}, "
                f"nearest distance={d:.6g}, tol={ceintf_tol:.6g}"
            )

        nearest = ranked[: max(1, min(interp_k, len(ranked)))]
        if nearest[0][0] < 1.0e-24:
            interface_map[fnid] = [(nearest[0][1], 1.0)]
            continue

        invs = []
        for d2, snid in nearest:
            invs.append((snid, 1.0 / max(d2, 1.0e-24)))
        isum = sum(v for _, v in invs)
        interface_map[fnid] = [(snid, val / isum) for snid, val in invs]

    # Pseudo-fluid stiffness scaling by element volume
    fluid_vols: List[Tuple[int, float]] = []
    all_fluid_elems_for_vol = fluid_hex_elems + fluid_tet_elems
    for eid, conn in all_fluid_elems_for_vol:
        if tet_volumes is not None and eid in tet_volumes:
            vol = tet_volumes[eid]
        else:
            coords = [new_nodes[n] for n in conn]
            unique_nodes = len(set(conn))
            if unique_nodes <= 4:
                seen = []
                for n in conn:
                    if n not in seen:
                        seen.append(n)
                a, b, c, d = [new_nodes[n] for n in seen[:4]]
                vol = tet_volume(a, b, c, d)
            else:
                vol = hex_volume(coords)
        if vol <= 0.0:
            raise ValueError(f"{case_name}: non-positive fluid element volume at eid={eid}")
        fluid_vols.append((eid, vol))
    min_vol = min(v for _, v in fluid_vols)

    stiff_groups: Dict[str, Dict[str, object]] = {}
    for eid, vol in fluid_vols:
        ex = fluid_base_e * (min_vol / vol) ** fluid_stiff_power
        key = f"{ex:.{fluid_stiff_digits}e}"
        grp = stiff_groups.get(key)
        if grp is None:
            stiff_groups[key] = {"ex": ex, "eids": [eid]}
        else:
            grp["eids"].append(eid)  # type: ignore[index]

    # Write CalculiX input
    out_dir.mkdir(parents=True, exist_ok=True)
    inp_path = out_dir / f"{case_name}.inp"
    meta_path = out_dir / f"{case_name}.json"

    with inp_path.open("w", encoding="utf-8") as fp:
        fp.write("** ------------------------------------------------------------\n")
        fp.write("** CalculiX modal case generated from ANSYS CDB sources\n")
        fp.write(f"** case_name: {case_name}\n")
        fp.write("** ------------------------------------------------------------\n")

        fp.write("*NODE\n")
        for nid in sorted(new_nodes):
            x, y, z = new_nodes[nid]
            fp.write(f"{nid}, {x:.12e}, {y:.12e}, {z:.12e}\n")

        fp.write("*ELEMENT, TYPE=C3D8, ELSET=E_FLUID_HEX\n")
        for eid, conn in fluid_hex_elems:
            fp.write(f"{eid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}, {conn[4]}, {conn[5]}, {conn[6]}, {conn[7]}\n")

        if fluid_tet_elems:
            fp.write("*ELEMENT, TYPE=C3D4, ELSET=E_FLUID_TET\n")
            for eid, conn in fluid_tet_elems:
                # Write only first 4 unique nodes for C3D4
                seen = []
                for n in conn:
                    if n not in seen:
                        seen.append(n)
                n1, n2, n3, n4 = seen
                # Ensure positive Jacobian: swap n2,n3 if volume is negative
                sv = tet_volume_signed(new_nodes[n1], new_nodes[n2], new_nodes[n3], new_nodes[n4])
                if sv < 0:
                    n2, n3 = n3, n2
                fp.write(f"{eid}, {n1}, {n2}, {n3}, {n4}\n")

        fp.write("*ELEMENT, TYPE=C3D8, ELSET=E_SOLID_ALL\n")
        for eid, conn in new_solid_elems:
            fp.write(f"{eid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}, {conn[4]}, {conn[5]}, {conn[6]}, {conn[7]}\n")

        write_nset(fp, "N_FLUIDDOMAIN", nset_fluiddomain)
        write_nset(fp, "N_SOLID", nset_solid)
        write_nset(fp, "N_FASTDYNAMIC", nset_fastdynamic)
        write_nset(fp, "N_FSI_FLUID", nset_fsi_fluid)
        write_nset(fp, "N_FSI_SOLID", nset_fsi_solid)
        write_nset(fp, "N_FIX_ALL", nset_fix_all)
        write_nset(fp, "N_FIX_UZ", nset_fix_uz)

        for sname, ids in fluid_sets.items():
            write_nset(fp, f"N_{sname}", ids)
        for sname, ids in solid_sets.items():
            write_nset(fp, f"N_{sname}", ids)

        # Fluid stiffness groups -> materials and element sets
        sorted_groups = sorted(stiff_groups.items(), key=lambda kv: float(kv[0]))
        fluid_group_defs: List[Tuple[str, str, float, List[int]]] = []
        for idx, (_, data) in enumerate(sorted_groups, start=1):
            elset_name = f"E_FLUID_G{idx}"
            mat_name = f"MAT_FLUID_G{idx}"
            ex = float(data["ex"])  # type: ignore[arg-type]
            eids = sorted(data["eids"])  # type: ignore[arg-type]
            fluid_group_defs.append((elset_name, mat_name, ex, eids))
            write_elset(fp, elset_name, eids)

        fp.write("*MATERIAL, NAME=MAT_SOLID\n")
        fp.write("*ELASTIC\n")
        fp.write(f"{solid_e:.6e}, {solid_nu:.6e}\n")
        fp.write("*DENSITY\n")
        fp.write(f"{solid_rho:.6e}\n")

        for _, mat_name, ex, _ in fluid_group_defs:
            fp.write(f"*MATERIAL, NAME={mat_name}\n")
            fp.write("*ELASTIC\n")
            fp.write(f"{ex:.12e}, 0.0\n")
            fp.write("*DENSITY\n")
            fp.write("0.0\n")

        fp.write("*SOLID SECTION, ELSET=E_SOLID_ALL, MATERIAL=MAT_SOLID\n")
        fp.write("\n")
        for elset_name, mat_name, _, _ in fluid_group_defs:
            fp.write(f"*SOLID SECTION, ELSET={elset_name}, MATERIAL={mat_name}\n")
            fp.write("\n")
        if fluid_tet_elems:
            fp.write("*SOLID SECTION, ELSET=E_FLUID_TET, MATERIAL=MAT_FLUID_G1\n")
            fp.write("\n")

        fp.write("*BOUNDARY\n")
        fp.write("N_FIX_ALL, 1, 3, 0.0\n")
        fp.write("N_FIX_UZ, 3, 3, 0.0\n")

        for fnid, targets in sorted(interface_map.items()):
            for dof in (1, 2, 3):
                terms = [(fnid, dof, 1.0)]
                for snid, w in targets:
                    terms.append((snid, dof, -w))
                write_equation(fp, terms)

        # TIE constraints: bond TIE_MASTER and TIE_SLAVE surfaces (DISABLED - causes MPC issues)
        # tie_master = fluid_sets.get("TIE_MASTER", [])
        # tie_slave = fluid_sets.get("TIE_SLAVE", [])
        # ... TIE code disabled for now

        fp.write("*STEP\n")
        fp.write("*FREQUENCY\n")
        fp.write(f"{int(n_modes)}\n")
        fp.write("*NODE FILE, NSET=N_FASTDYNAMIC\n")
        fp.write("U\n")
        fp.write("*END STEP\n")

    meta = {
        "case_name": case_name,
        "fluid_nodes": len(fluid_nodes_sorted),
        "solid_nodes": len(solid_nodes_sorted),
        "fastdynamic_output_nodes": len(nset_fastdynamic),
        "fluid_elements": len(fluid_hex_elems) + len(fluid_tet_elems),
        "fluid_hex_elements": len(fluid_hex_elems),
        "fluid_tet_elements": len(fluid_tet_elems),
        "solid_elements": len(new_solid_elems),
        "fluid_stiffness_groups": len(fluid_group_defs),
        "fsi_fluid_nodes": len(nset_fsi_fluid),
        "fsi_fluid_coupled_nodes": len(nset_fsi_fluid_coupled),
        "fsi_solid_nodes": len(nset_fsi_solid),
        "interface_tolerance": ceintf_tol,
        "interface_interp_k": interp_k,
        "fluid_stiff_round_digits": fluid_stiff_digits,
        "n_modes": n_modes,
        "sources": {
            "solid_cdb": str(solid_cdb_path.resolve()),
            "fluid_cdb": str(fluid_cdb_path.resolve()),
        },
    }
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return inp_path


def run_ccx(inp_path: Path, ccx_cmd: str) -> None:
    base = inp_path.with_suffix("")
    subprocess.run([ccx_cmd, "-i", str(base)], check=True, cwd=str(inp_path.parent))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build CalculiX modal .inp files for multiple fluid meshes with one shared solid mesh."
    )
    parser.add_argument("--solid-cdb", required=True, help="Reference DamFailure.in containing the structural mesh")
    parser.add_argument(
        "--fluid-cdb",
        action="append",
        required=True,
        help="Fluid DamFailure.in file. Repeat this argument for multiple fluid meshes.",
    )
    parser.add_argument("--out-dir", default="calculix_modal_cases", help="Output directory for generated .inp files")
    parser.add_argument("--n-modes", type=int, default=10, help="Number of eigenmodes")
    parser.add_argument("--fluid-base-e", type=float, default=0.1, help="Base pseudo-fluid E")
    parser.add_argument(
        "--fluid-stiff-power",
        type=float,
        default=2.5,
        help="Power for pseudo-fluid stiffness scaling: E = baseE*(minV/V)^power",
    )
    parser.add_argument(
        "--fluid-stiff-digits",
        type=int,
        default=12,
        help="Scientific-notation rounding digits used to group pseudo-fluid materials",
    )
    parser.add_argument("--ceintf-tol", type=float, default=0.05, help="Interface coupling tolerance (length)")
    parser.add_argument("--interp-k", type=int, default=4, help="Nearest solid nodes used for each fluid interface node")
    parser.add_argument("--run-ccx", action="store_true", help="Run ccx after generating each .inp")
    parser.add_argument("--ccx-cmd", default="ccx", help="CalculiX executable")
    parser.add_argument("--solid-e", type=float, default=1.0e6, help="Solid Young's modulus [Pa]")
    parser.add_argument("--solid-nu", type=float, default=0.0, help="Solid Poisson ratio")
    parser.add_argument("--solid-rho", type=float, default=2500.0, help="Solid density [kg/m3]")
    parser.add_argument("--tet-volumes-file", help="File with precomputed tetra element volumes (eid vol per line)")
    args = parser.parse_args()

    solid_path = Path(args.solid_cdb).resolve()
    fluid_paths = [Path(p).resolve() for p in args.fluid_cdb]
    out_dir = Path(args.out_dir).resolve()

    solid_model = parse_cdb(solid_path)

    # Load precomputed tetra volumes if provided
    tet_volumes: Dict[int, float] | None = None
    if args.tet_volumes_file:
        tet_volumes = {}
        vpath = Path(args.tet_volumes_file).resolve()
        for line in vpath.read_text().splitlines():
            parts = line.strip().split()
            if len(parts) >= 2:
                tet_volumes[int(parts[0])] = float(parts[1])
        print(f"[INFO] Loaded {len(tet_volumes)} element volumes from {vpath}")

    generated: List[Path] = []
    used_names = set()
    for idx, fpath in enumerate(fluid_paths, start=1):
        fluid_model = parse_cdb(fpath)

        name = fpath.stem
        if name in used_names:
            name = f"{name}_{idx}"
        used_names.add(name)

        inp = build_case(
            fluid_model=fluid_model,
            solid_model=solid_model,
            fluid_cdb_path=fpath,
            solid_cdb_path=solid_path,
            case_name=name,
            out_dir=out_dir,
            n_modes=args.n_modes,
            fluid_base_e=args.fluid_base_e,
            fluid_stiff_power=args.fluid_stiff_power,
            fluid_stiff_digits=args.fluid_stiff_digits,
            ceintf_tol=args.ceintf_tol,
            interp_k=args.interp_k,
            solid_e=args.solid_e,
            solid_nu=args.solid_nu,
            solid_rho=args.solid_rho,
            tet_volumes=tet_volumes,
        )
        generated.append(inp)
        print(f"[OK] generated: {inp}")

    if args.run_ccx:
        ccx_path = shutil.which(args.ccx_cmd)
        if ccx_path is None:
            raise SystemExit(f"ccx executable not found: {args.ccx_cmd}")
        for inp in generated:
            print(f"[RUN] {ccx_path} -i {inp.with_suffix('')}")
            run_ccx(inp, ccx_path)

    print(f"[DONE] total generated cases: {len(generated)}")


if __name__ == "__main__":
    main()
