#!/usr/bin/env python3
"""
Export CalculiX modal eigenvalue output from .dat into CSV/JSON summaries.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List


FIELDS = (
    "mode",
    "eigenvalue",
    "omega_rad_per_time",
    "frequency_cycles_per_time",
    "imag_rad_per_time",
)


def parse_modal_dat(path: Path) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    in_eigen_table = False

    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if "E I G E N V A L U E" in line and "O U T P U T" in line:
            in_eigen_table = True
            continue
        if not in_eigen_table:
            continue
        if "P A R T I C I P A T I O N" in line:
            break

        parts = line.split()
        if len(parts) != 5 or not parts[0].isdigit():
            continue

        mode = int(parts[0])
        values = [float(v.replace("D", "E").replace("d", "e")) for v in parts[1:5]]
        rows.append(
            {
                "mode": mode,
                "eigenvalue": values[0],
                "omega_rad_per_time": values[1],
                "frequency_cycles_per_time": values[2],
                "imag_rad_per_time": values[3],
            }
        )

    if not rows:
        raise ValueError(f"No modal eigenvalue rows found in {path}")
    return rows


def write_single_case_csv(rows: List[Dict[str, float]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(FIELDS))
        writer.writeheader()
        writer.writerows(rows)


def write_multi_case_csv(cases: Dict[str, List[Dict[str, float]]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["case", *FIELDS])
        writer.writeheader()
        for case_name, rows in cases.items():
            for row in rows:
                writer.writerow({"case": case_name, **row})


def main() -> None:
    parser = argparse.ArgumentParser(description="Export modal info from CalculiX .dat file(s)")
    parser.add_argument(
        "--dat",
        action="append",
        required=True,
        help="Path to modal .dat file. Repeat this option for multiple files.",
    )
    parser.add_argument("--csv", help="CSV output path (default: <dat_stem>_modal_summary.csv)")
    parser.add_argument("--json", help="JSON output path (default: <dat_stem>_modal_summary.json)")
    args = parser.parse_args()

    dat_files = [Path(p).resolve() for p in args.dat]
    if len(dat_files) == 1:
        dat_path = dat_files[0]
        case_rows = parse_modal_dat(dat_path)
        csv_out = Path(args.csv).resolve() if args.csv else dat_path.with_name(f"{dat_path.stem}_modal_summary.csv")
        json_out = Path(args.json).resolve() if args.json else dat_path.with_name(f"{dat_path.stem}_modal_summary.json")

        write_single_case_csv(case_rows, csv_out)
        json_out.write_text(json.dumps(case_rows, indent=2), encoding="utf-8")

        print(f"[OK] CSV:  {csv_out}")
        print(f"[OK] JSON: {json_out}")
        print(f"[OK] modes exported: {len(case_rows)}")
        return

    case_data: Dict[str, List[Dict[str, float]]] = {}
    for dat_path in dat_files:
        case_data[dat_path.stem] = parse_modal_dat(dat_path)

    first = dat_files[0]
    csv_out = Path(args.csv).resolve() if args.csv else first.parent / "modal_summary_multi_case.csv"
    json_out = Path(args.json).resolve() if args.json else first.parent / "modal_summary_multi_case.json"

    write_multi_case_csv(case_data, csv_out)
    json_out.write_text(json.dumps(case_data, indent=2), encoding="utf-8")

    print(f"[OK] CSV:  {csv_out}")
    print(f"[OK] JSON: {json_out}")
    print(f"[OK] cases exported: {len(case_data)}")


if __name__ == "__main__":
    main()
