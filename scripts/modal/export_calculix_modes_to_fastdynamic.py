#!/usr/bin/env python3
"""
Export CalculiX modal displacement results to fastDynamicFvMesh mode CSV files.

Inputs:
  - generated CalculiX .inp from build_calculix_modal_cases.py
  - CalculiX .dat with eigenvalue/frequency table
  - CalculiX .frd with DISP blocks for N_FASTDYNAMIC

Outputs:
  - FluidNodeCoor.csv
  - FluidNodeDisp1.csv ... FluidNodeDispN.csv
  - StructNodeCoor.csv
  - StructNodeDisp1.csv ... StructNodeDispN.csv
  - FluidPara.csv (copied from existing mode dir when available)
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
NSET_RE = re.compile(r"^\*NSET\s*,\s*NSET=([^,\s]+)", re.IGNORECASE)


def parse_inp(path: Path) -> Tuple[Dict[int, Tuple[float, float, float]], Dict[str, List[int]]]:
    nodes: Dict[int, Tuple[float, float, float]] = {}
    nsets: Dict[str, List[int]] = {}
    mode = None
    current_set = None

    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("**"):
            continue

        upper = line.upper()
        if upper.startswith("*NODE"):
            mode = "node"
            current_set = None
            continue

        m = NSET_RE.match(line)
        if m:
            mode = "nset"
            current_set = m.group(1).strip().upper()
            nsets.setdefault(current_set, [])
            continue

        if line.startswith("*"):
            mode = None
            current_set = None
            continue

        if mode == "node":
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 4:
                nid = int(parts[0])
                nodes[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
        elif mode == "nset" and current_set:
            ids = [int(x) for x in re.findall(r"\d+", line)]
            nsets[current_set].extend(ids)

    for name, ids in list(nsets.items()):
        nsets[name] = sorted(set(ids))

    if not nodes:
        raise ValueError(f"No nodes parsed from {path}")
    return nodes, nsets


def parse_dat_frequencies(path: Path) -> List[float]:
    freqs: List[float] = []
    in_table = False
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.strip()
        if "E I G E N V A L U E" in line and "O U T P U T" in line:
            in_table = True
            continue
        if not in_table:
            continue
        if "P A R T I C I P A T I O N" in line:
            break
        parts = line.split()
        if len(parts) == 5 and parts[0].isdigit():
            freqs.append(float(parts[3].replace("D", "E").replace("d", "e")))
    if not freqs:
        raise ValueError(f"No eigenfrequencies parsed from {path}")
    return freqs


def parse_frd_displacements(path: Path) -> Dict[int, Dict[int, Tuple[float, float, float]]]:
    modes: Dict[int, Dict[int, Tuple[float, float, float]]] = {}
    mode_no = None
    in_disp = False
    current: Dict[int, Tuple[float, float, float]] = {}

    for raw in path.read_text(errors="ignore").splitlines():
        if raw.startswith("    1PMODE"):
            vals = re.findall(r"\d+", raw)
            if vals:
                mode_no = int(vals[-1])
            continue

        if raw.startswith(" -4") and "DISP" in raw:
            if mode_no is None:
                raise ValueError(f"DISP block without PMODE in {path}")
            in_disp = True
            current = {}
            continue

        if not in_disp:
            continue

        if raw.startswith(" -3"):
            if mode_no is not None:
                modes[mode_no] = current
            in_disp = False
            current = {}
            continue

        if raw.startswith(" -1"):
            try:
                nid = int(raw[3:13])
            except ValueError:
                parts = raw.split()
                if len(parts) < 2:
                    continue
                nid = int(parts[1])
                payload = " ".join(parts[2:])
            else:
                payload = raw[13:]

            vals = [
                float(x.replace("D", "E").replace("d", "e"))
                for x in FLOAT_RE.findall(payload)
            ]
            if len(vals) >= 3:
                current[nid] = (vals[0], vals[1], vals[2])

    if not modes:
        raise ValueError(f"No displacement blocks parsed from {path}")
    return modes


def write_vec_csv(path: Path, header_a: float, count: int, n_modes: int, rows: Iterable[Tuple[float, float, float]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{header_a:20.10E},{float(count):20.10E},{float(n_modes):20.10E},\n")
        for x, y, z in rows:
            f.write(f"{x:20.10E},{y:20.10E},{z:20.10E},\n")


def write_mode_dir(
    out_dir: Path,
    nodes: Dict[int, Tuple[float, float, float]],
    fluid_ids: List[int],
    solid_ids: List[int],
    freqs: List[float],
    modes: Dict[int, Dict[int, Tuple[float, float, float]]],
    n_modes: int,
    fluid_para_source: Path | None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    missing_modes = [m for m in range(1, n_modes + 1) if m not in modes]
    if missing_modes:
        raise ValueError(f"Missing displacement blocks for modes: {missing_modes}")

    for label, ids in (("fluid", fluid_ids), ("solid", solid_ids)):
        missing_nodes = [
            nid
            for nid in ids
            for mode_no in range(1, n_modes + 1)
            if nid not in modes[mode_no]
        ]
        if missing_nodes:
            unique_missing = sorted(set(missing_nodes))
            sample = ", ".join(str(x) for x in unique_missing[:10])
            raise ValueError(
                f"FRD displacement output is missing {len(unique_missing)} {label} nodes "
                f"(sample: {sample}). Regenerate CalculiX input with N_FASTDYNAMIC output and rerun ccx."
            )

    write_vec_csv(out_dir / "FluidNodeCoor.csv", 0.0, len(fluid_ids), n_modes, (nodes[nid] for nid in fluid_ids))
    write_vec_csv(out_dir / "StructNodeCoor.csv", 0.0, len(solid_ids), n_modes, (nodes[nid] for nid in solid_ids))

    for mode_no in range(1, n_modes + 1):
        freq = freqs[mode_no - 1]
        disp = modes[mode_no]
        write_vec_csv(
            out_dir / f"FluidNodeDisp{mode_no}.csv",
            freq,
            len(fluid_ids),
            n_modes,
            (disp[nid] for nid in fluid_ids),
        )
        write_vec_csv(
            out_dir / f"StructNodeDisp{mode_no}.csv",
            freq,
            len(solid_ids),
            n_modes,
            (disp[nid] for nid in solid_ids),
        )

    fluid_para_target = out_dir / "FluidPara.csv"
    if fluid_para_source and fluid_para_source.exists():
        if fluid_para_source.resolve() != fluid_para_target.resolve():
            shutil.copy2(fluid_para_source, fluid_para_target)
    elif not (out_dir / "FluidPara.csv").exists():
        fluid_para_target.write_text("19,\n11,\n0.0,\n0.0,\n0.0,\n0.0,\n0.0,\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Export CalculiX modal FRD data to fastDynamicFvMesh mode CSV files")
    parser.add_argument("--inp", required=True, help="Generated CalculiX .inp")
    parser.add_argument("--dat", required=True, help="CalculiX .dat")
    parser.add_argument("--frd", required=True, help="CalculiX .frd")
    parser.add_argument("--out-mode-dir", required=True, help="Output mode directory")
    parser.add_argument("--n-modes", type=int, default=10, help="Number of modes to export")
    parser.add_argument("--fluid-para-source", help="Existing FluidPara.csv to copy")
    args = parser.parse_args()

    nodes, nsets = parse_inp(Path(args.inp))
    fluid_ids = nsets.get("N_FLUIDDOMAIN")
    solid_ids = nsets.get("N_SOLID")
    if not fluid_ids:
        raise SystemExit("N_FLUIDDOMAIN not found or empty in INP")
    if not solid_ids:
        raise SystemExit("N_SOLID not found or empty in INP")

    freqs = parse_dat_frequencies(Path(args.dat))
    if len(freqs) < args.n_modes:
        raise SystemExit(f"DAT has only {len(freqs)} modes, requested {args.n_modes}")

    modes = parse_frd_displacements(Path(args.frd))
    write_mode_dir(
        Path(args.out_mode_dir),
        nodes,
        fluid_ids,
        solid_ids,
        freqs,
        modes,
        args.n_modes,
        Path(args.fluid_para_source) if args.fluid_para_source else None,
    )

    print(f"[OK] wrote fastDynamicFvMesh mode files to {Path(args.out_mode_dir).resolve()}")
    print(f"[OK] fluid nodes: {len(fluid_ids)}, solid nodes: {len(solid_ids)}, modes: {args.n_modes}")


if __name__ == "__main__":
    main()
