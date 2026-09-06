#!/usr/bin/env python3
"""Convert ICEM-exported ANSYS CDB-like file (NBLOCK/EBLOCK/CMBLOCK) to CalculiX .inp.

Handles:
  NBLOCK,6,SOLID  (3i8,6e16.9)  -> *NODE
  EBLOCK,19,SOLID (19i10)       -> *ELEMENT, TYPE=C3D8 (grouped by TYPE field)
  CMBLOCK,name,NODE,n (8i10)    -> *NSET

Single streaming pass, O(1) memory.
"""
import re
import sys

SRC = "/root/Workspace/Oldversion/AGARD4456/Mesh/Assembly/Assembly.in"
DST = "/root/Workspace/Oldversion/AGARD4456/Mesh/Assembly/Assembly.inp"

INT_RE = re.compile(r"[+-]?\d+")


def sanitize(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_\-]", "_", name.strip())
    if not s or not s[0].isalpha():
        s = "S_" + s
    return s[:80]


def main() -> None:
    n_nodes = 0
    n_elems = 0
    elem_stats = {}   # type -> [count, min_node, max_node]
    nsets = []
    state = None
    fmt_pending = False
    cur_type = None
    cmb_name = None
    cmb_expected = 0
    cmb_got = 0
    cmb_buf = []

    with open(SRC, "r", encoding="ascii", errors="strict") as fin, \
         open(DST, "w", encoding="ascii") as out:
        out.write("** Converted from ICEM ANSYS export Assembly.in\n")
        out.write("** SOLID185 -> C3D8; EBLOCK TYPE field -> ELSET ETYPE<n>\n")

        def flush_cmb():
            nonlocal cmb_buf
            if cmb_name is None:
                return
            out.write(f"*NSET, NSET={cmb_name}\n")
            for i in range(0, len(cmb_buf), 8):
                out.write(",".join(cmb_buf[i:i + 8]) + "\n")
            nsets.append((cmb_name, len(cmb_buf)))
            cmb_buf = []

        for line in fin:
            c0 = line[0:1]
            if state is None:
                if line.startswith("NBLOCK"):
                    state = "nb_fmt"
                elif line.startswith("EBLOCK"):
                    state = "eb_fmt"
                elif line.startswith("CMBLOCK"):
                    parts = line.strip().split(",")
                    if len(parts) >= 4 and parts[2].strip().upper() == "NODE":
                        cmb_name = sanitize(parts[1])
                        cmb_expected = int(parts[3])
                        cmb_got = 0
                        cmb_buf = []
                        state = "cmb_fmt"
                    else:
                        print(f"[skip] non-NODE CMBLOCK: {line.strip()}")
                continue

            if state == "nb_fmt":
                out.write("*NODE\n")
                state = "nb"
                continue
            if state == "eb_fmt":
                state = "eb"
                continue
            if state == "cmb_fmt":
                state = "cmb"
                continue

            if state == "nb":
                if c0.isalpha() or line.startswith("-1") or c0 in "/!":
                    state = None
                    continue
                # (3i8,6e16.9): nid | solid | line | x y z thxy thyz thzx
                nid = line[0:8].strip()
                x = line[24:40].strip()
                y = line[40:56].strip()
                z = line[56:72].strip()
                out.write(f"{nid},{x},{y},{z}\n")
                n_nodes += 1
                continue

            if state == "eb":
                s = line.strip()
                if s == "-1" or s.startswith("EN ") or (c0 not in " 0123456789+-"):
                    state = None
                    cur_type = None
                    continue
                # (19i10): mat type real sec esys death solid shape excl ? eid n1..n8
                typ = int(line[10:20])
                eid = line[100:110].strip()
                if typ != cur_type:
                    out.write(f"*ELEMENT, TYPE=C3D8, ELSET=ETYPE{typ}\n")
                    cur_type = typ
                    elem_stats.setdefault(typ, [0, None, None])
                n1 = int(line[110:120])
                st = elem_stats[typ]
                st[0] += 1
                st[1] = n1 if st[1] is None else min(st[1], n1)
                for i in range(1, 8):
                    v = int(line[110 + i * 10:120 + i * 10])
                    st[1] = min(st[1], v)
                    st[2] = v if st[2] is None else max(st[2], v)
                st[2] = max(st[2], n1)
                out.write(eid + "," + ",".join(
                    line[110 + i * 10:120 + i * 10].strip() for i in range(8)) + "\n")
                n_elems += 1
                continue

            if state == "cmb":
                s = line.strip()
                if s == "-1" or (c0 not in " 0123456789+-"):
                    flush_cmb()
                    cmb_name = None
                    state = None
                    # re-dispatch this line (could be a new CMBLOCK etc.)
                    if line.startswith("CMBLOCK"):
                        parts = line.strip().split(",")
                        if len(parts) >= 4 and parts[2].strip().upper() == "NODE":
                            cmb_name = sanitize(parts[1])
                            cmb_expected = int(parts[3])
                            cmb_got = 0
                            cmb_buf = []
                            state = "cmb_fmt"
                    continue
                vals = INT_RE.findall(line)
                cmb_buf.extend(vals)
                cmb_got += len(vals)
                continue

        flush_cmb()

    print(f"nodes:   {n_nodes}")
    print(f"elements:{n_elems}")
    for typ, (cnt, mn, mx) in sorted(elem_stats.items()):
        print(f"  ETYPE{typ}: {cnt} elems, node range {mn}..{mx}")
    print(f"node sets ({len(nsets)}):")
    for name, cnt in nsets:
        print(f"  {name}: {cnt}")


if __name__ == "__main__":
    main()
