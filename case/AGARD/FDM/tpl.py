#!/usr/bin/env python3
"""Shared access to the zipped case templates.

`templates/myHisa.zip` holds the transient dictionary set, `templates/steady.zip`
the steady-state one. They are shipped as zips so the repository carries two
small files instead of ~20 loose dictionaries.
"""
import os
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
TPL_DIR = os.path.join(HERE, "templates")
MYHISA_ZIP = os.path.join(TPL_DIR, "myHisa.zip")
STEADY_ZIP = os.path.join(TPL_DIR, "steady.zip")


def _open(zip_path, root):
    if not os.path.isfile(zip_path):
        raise SystemExit("template archive missing: %s" % zip_path)
    return zipfile.ZipFile(zip_path), root


def text(zip_path, root, member):
    """Read one template as text, e.g. text(MYHISA_ZIP, 'myHisa', 'system/fvSolution')."""
    z, _ = _open(zip_path, root)
    with z:
        return z.read(root + "/" + member).decode()


def names(zip_path, root, prefix):
    """Members under <root>/<prefix> (files only, no sub-directories)."""
    z, _ = _open(zip_path, root)
    base = root + "/" + prefix
    out = []
    with z:
        for n in z.namelist():
            if not n.startswith(base) or n.endswith("/"):
                continue
            out.append(n[len(base):])
    return sorted(out)


def copy_out(zip_path, root, member, dst):
    """Extract one member to dst (binary safe)."""
    z, _ = _open(zip_path, root)
    with z:
        data = z.read(root + "/" + member)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "wb") as f:
        f.write(data)
