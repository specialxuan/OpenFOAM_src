#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenerate the SVG figures for
`fastDynamicFvMesh/Doc/AMR线性插值下子网格质量与母网格质量的关系.md`.

Dependencies (project venv `.venv_fdm`):
    matplotlib >= 3.11, numpy

Run from the repository root:
    .venv_fdm/bin/python fastDynamicFvMesh/Doc/Figures/amr_mesh_quality/gen_figures.py

SVGs are written next to this script (git-tracked via the .gitignore exception);
QA PNG copies go to /root/Workspace/tmp/amr_mesh_quality_figures/ (not tracked).
"""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
QA_DIR = "/root/Workspace/tmp/amr_mesh_quality_figures"
os.makedirs(QA_DIR, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": ["Noto Sans CJK SC", "Noto Serif CJK SC", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "svg.fonttype": "path",  # portable: text as glyph paths
        "figure.dpi": 110,
    }
)

C_EDGE = "0.25"   # parent wireframe
C_CUT = "0.60"    # cut planes / helper lines
C_CHILD = "C0"    # highlighted child
C_EMID = "C3"     # edge midpoints
C_FMID = "C2"     # face centres
C_CMID = "C1"     # body centre


# --------------------------------------------------------------------------- #
#  orthographic 3D -> 2D projection (2.5D approach; avoids mplot3d z-order)
# --------------------------------------------------------------------------- #
def view_dir(azim_deg, elev_deg):
    a, e = np.radians(azim_deg), np.radians(elev_deg)
    return np.array([np.cos(a) * np.cos(e), np.sin(a) * np.cos(e), np.sin(e)])


def project(pts, azim_deg, elev_deg):
    """Return (screen xy, depth-toward-viewer) for an (...,3) array."""
    pts = np.atleast_2d(np.asarray(pts, dtype=float))
    a, e = np.radians(azim_deg), np.radians(elev_deg)
    sa, ca, se, ce = np.sin(a), np.cos(a), np.sin(e), np.cos(e)
    x2d = -sa * pts[:, 0] + ca * pts[:, 1]
    y2d = -ca * se * pts[:, 0] - sa * se * pts[:, 1] + ce * pts[:, 2]
    depth = ca * ce * pts[:, 0] + sa * ce * pts[:, 1] + se * pts[:, 2]
    return np.column_stack([x2d, y2d]), depth


def draw_edge(ax, p2, a, b, **kw):
    seg = np.array([p2[a], p2[b]])
    ax.plot(seg[:, 0], seg[:, 1], **kw)


def arrow3d(ax, p_from, p_to, color, azim, elev, lw=1.8, ls="-"):
    x1, _ = project(p_from, azim, elev)
    x2, _ = project(p_to, azim, elev)
    ax.add_patch(
        FancyArrowPatch(
            x1[0], x2[0], arrowstyle="-|>", mutation_scale=14,
            color=color, lw=lw, linestyle=ls, zorder=8,
        )
    )


def save(fig, name):
    fig.savefig(os.path.join(HERE, name + ".svg"), bbox_inches="tight")
    fig.savefig(os.path.join(QA_DIR, name + ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# unit-cube topology used by several figures
V_UNIT = np.array(
    [
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ],
    dtype=float,
)
EDGES = [
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
]
FACES = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4],
         [3, 2, 6, 7], [0, 3, 7, 4], [1, 2, 6, 5]]


def hidden_edges(V2faces_normals, vdir):
    vis = [np.dot(n, vdir) > 0 for n in V2faces_normals]

    def hid(edge):
        a, b = edge
        adj = [i for i, f in enumerate(FACES) if a in f and b in f]
        return len(adj) == 2 and not any(vis[i] for i in adj)

    return hid


def face_normals(V):
    return [np.cross(V[f[1]] - V[f[0]], V[f[2]] - V[f[0]]) for f in FACES]


# --------------------------------------------------------------------------- #
#  Fig 1: 2x2x2 subdomain sampling — parent hex, one corner child highlighted
# --------------------------------------------------------------------------- #
def fig_subdomain():
    AZ, EL = 24.0, 18.0
    V = V_UNIT.copy()
    V2, _ = project(V, AZ, EL)
    hid = hidden_edges(face_normals(V), view_dir(AZ, EL))

    fig, ax = plt.subplots(figsize=(6.4, 5.6))

    # parent wireframe
    for e in EDGES:
        draw_edge(ax, V2, *e, color=C_EDGE, lw=1.3,
                  ls=(0, (5, 3)) if hid(e) else "-")

    # mid-plane cut lines (dashed): segments connecting edge midpoints
    em = {}
    for e in EDGES:
        m = (V[e[0]] + V[e[1]]) / 2
        em[e] = m
        em[(e[1], e[0])] = m
    cuts = [
        # xi = 1/2 plane (visible frame on parent faces)
        ((0, 4), (1, 5)), ((1, 5), (2, 6)), ((2, 6), (3, 7)), ((3, 7), (0, 4)),
        # eta = 1/2 plane
        ((0, 1), (3, 2)), ((3, 2), (7, 6)), ((7, 6), (4, 5)), ((4, 5), (0, 1)),
        # zeta = 1/2 plane
        ((0, 3), (1, 2)), ((1, 2), (5, 6)), ((5, 6), (4, 7)), ((4, 7), (0, 3)),
    ]
    for e1, e2 in cuts:
        p1, _ = project(em[e1], AZ, EL)
        p2, _ = project(em[e2], AZ, EL)
        ax.plot([p1[0, 0], p2[0, 0]], [p1[0, 1], p2[0, 1]],
                color=C_CUT, lw=0.9, ls=(0, (3, 2)), zorder=2)

    # child (0,0,0): octant [0, 0.5]^3 in parent coordinates
    s = 0.5
    Vc = s * V_UNIT  # child corners in parent physical coords
    Vc2, _ = project(Vc, AZ, EL)
    for e in EDGES:
        draw_edge(ax, Vc2, *e, color=C_CHILD, lw=2.2, zorder=6)

    # new points of the child
    e_mids = np.array([[s, 0, 0], [0, s, 0], [0, 0, s]])
    f_mids = np.array([[s, s, 0], [s, 0, s], [0, s, s]])
    c_mid = np.array([[s, s, s]])
    for pts, mk, col, lab in [
        (e_mids, "s", C_EMID, "边中点 (2 角点平均)"),
        (f_mids, "^", C_FMID, "面心 (4 角点平均)"),
        (c_mid, "*", C_CMID, "体心 (8 角点平均)"),
    ]:
        p2, _ = project(pts, AZ, EL)
        ax.plot(p2[:, 0], p2[:, 1], mk, ms=9 if mk == "*" else 6,
                color=col, zorder=9, label=lab)

    # label parent corner P000 and the child
    p000, _ = project(V[0], AZ, EL)
    ax.annotate(r"$\mathbf{P}_{000}$", (p000[0, 0], p000[0, 1]),
                xytext=(-16, -14), textcoords="offset points",
                fontsize=11, color=C_EDGE)
    centroid2, _ = project(np.array([[0.25, 0.25, 0.25]]), AZ, EL)
    ax.annotate("角部子单元\n" r"$\mathbf{x}_P(\frac{a+r}{2},\frac{b+s}{2},\frac{c+t}{2})$",
                (centroid2[0, 0], centroid2[0, 1]),
                xytext=(30, -55), textcoords="offset points", fontsize=9.5,
                color=C_CHILD,
                arrowprops=dict(arrowstyle="->", color=C_CHILD, lw=1.1))

    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_title("2×2×2 加密：新点 = 母单元三线性映射的参数子域采样", fontsize=11)
    save(fig, "amr_hexa_subdomain")


# --------------------------------------------------------------------------- #
#  Fig 2: Jacobian / volume scaling — parameter parent -> child local cube
# --------------------------------------------------------------------------- #
def fig_jacobian():
    AZ, EL = 22.0, 16.0
    fig, ax = plt.subplots(figsize=(8.6, 4.6))

    # left: parent parameter cube with (0,0,0) subdomain highlighted
    offL = np.array([0.0, 0.0, 0.0])
    V = V_UNIT + offL
    V2, _ = project(V, AZ, EL)
    hid = hidden_edges(face_normals(V), view_dir(AZ, EL))
    for e in EDGES:
        draw_edge(ax, V2, *e, color=C_EDGE, lw=1.3,
                  ls=(0, (5, 3)) if hid(e) else "-")
    Vc = 0.5 * V_UNIT + offL
    Vc2, _ = project(Vc, AZ, EL)
    for e in EDGES:
        draw_edge(ax, Vc2, *e, color=C_CHILD, lw=2.0, zorder=6)
    ax.annotate("母参数域 $[0,1]^3$\n子域 "
                r"$[\frac{a}{2},\frac{a+1}{2}]\times[\frac{b}{2},\frac{b+1}{2}]\times[\frac{c}{2},\frac{c+1}{2}]$",
                (Vc2[6, 0], Vc2[6, 1]), xytext=(-118, 26),
                textcoords="offset points", fontsize=9.5, color=C_CHILD,
                arrowprops=dict(arrowstyle="->", color=C_CHILD, lw=1.1))

    # mapping arrow to the right
    x_from, _ = project(np.array([[1.15, 0.5, 0.5]]), AZ, EL)
    x_to, _ = project(np.array([[2.35, 0.5, 0.5]]), AZ, EL)
    ax.add_patch(FancyArrowPatch(x_from[0], x_to[0], arrowstyle="-|>",
                                 mutation_scale=20, color="0.15", lw=2.0))
    mid = 0.5 * (x_from[0] + x_to[0])
    ax.text(mid[0], mid[1] + 0.10,
            r"$T_{abc}(r,s,t)=\left(\frac{a+r}{2},\frac{b+s}{2},\frac{c+t}{2}\right)$",
            ha="center", fontsize=10.5)

    # right: child local parameter cube
    offR = np.array([3.1, 0.0, 0.0])
    Vr = V_UNIT + offR
    Vr2, _ = project(Vr, AZ, EL)
    hidr = hidden_edges(face_normals(Vr), view_dir(AZ, EL))
    for e in EDGES:
        draw_edge(ax, Vr2, *e, color=C_CHILD, lw=1.6,
                  ls=(0, (5, 3)) if hidr(e) else "-")
    ax.text(Vr2[2, 0], Vr2[2, 1] + 0.06, "子单元局部参数域 $[0,1]^3$",
            fontsize=9.5, color=C_CHILD, ha="center")

    # scaling facts
    ax.text(
        0.5 * (x_from[0, 0] + x_to[0, 0]), -0.62,
        "链式法则：$\\mathbf{J}_C=\\mathbf{J}_P\\,D T_{abc}=\\frac{1}{2}\\mathbf{J}_P$"
        "\n$\\det\\mathbf{J}_C=\\left(\\frac{1}{2}\\right)^3\\det\\mathbf{J}_P"
        "=\\frac{1}{8}\\det\\mathbf{J}_P$"
        "\n$V_C^{abc}=\\int_{\\Omega_{abc}}\\det\\mathbf{J}_P\\,du\\,dv\\,dw$,"
        "  $\\sum_{abc}V_C^{abc}=V_P$",
        ha="center", fontsize=10.5,
        bbox=dict(boxstyle="round,pad=0.45", fc="0.97", ec="0.75"),
    )

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("子单元是母单元参数子域的重新参数化：Jacobian 缩放 1/2，行列式缩放 1/8",
                 fontsize=11)
    save(fig, "jacobian_scaling")


# --------------------------------------------------------------------------- #
#  Fig 3: non-orthogonality of an internal face (2D projection, affine cells)
# --------------------------------------------------------------------------- #
def fig_nonorth():
    e1 = np.array([1.0, 0.0])
    e2 = np.array([0.5, 1.0])

    def P(xi, eta):
        return xi * e1 + eta * e2

    fig, ax = plt.subplots(figsize=(6.6, 4.6))

    left = [P(0, 0), P(0.5, 0), P(0.5, 1), P(0, 1)]
    right = [P(0.5, 0), P(1, 0), P(1, 1), P(0.5, 1)]
    for cell, col in [(left, C_CHILD), (right, "C7")]:
        poly = np.array(cell + [cell[0]])
        ax.plot(poly[:, 0], poly[:, 1], color=col, lw=1.8)

    # internal face
    f1, f2 = P(0.5, 0), P(0.5, 1)
    ax.plot([f1[0], f2[0]], [f1[1], f2[1]], color="0.15", lw=2.6)
    ax.annotate("内部面 $\\xi=1/2$", (f2[0], f2[1]), xytext=(6, -2),
                textcoords="offset points", fontsize=10)

    # centres and d
    cL = P(0.25, 0.5)
    cR = P(0.75, 0.5)
    ax.plot([cL[0]], [cL[1]], "o", color=C_CHILD, ms=6)
    ax.plot([cR[0]], [cR[1]], "o", color="C7", ms=6)
    ax.annotate("", (cR[0], cR[1]), (cL[0], cL[1]),
                arrowprops=dict(arrowstyle="-|>", color=C_CMID, lw=2.0))
    ax.text(0.5 * (cL[0] + cR[0]), cL[1] + 0.06,
            r"$\mathbf{d}_{int}=\frac{1}{2}\mathbf{e}_\xi$",
            ha="center", fontsize=10.5, color=C_CMID)
    ax.annotate(r"$\mathbf{c}_{own}$", (cL[0], cL[1]), xytext=(-12, -16),
                textcoords="offset points", fontsize=9.5, color=C_CHILD)
    ax.annotate(r"$\mathbf{c}_{nei}$", (cR[0], cR[1]), xytext=(2, -16),
                textcoords="offset points", fontsize=9.5, color="C7")

    # face normal S at face midpoint
    fm = 0.5 * (f1 + f2)
    t = f2 - f1
    n = np.array([t[1], -t[0]])
    n = n / np.linalg.norm(n) * 0.42
    ax.annotate("", (fm[0] + n[0], fm[1] + n[1]), (fm[0], fm[1]),
                arrowprops=dict(arrowstyle="-|>", color=C_FMID, lw=2.0))
    ax.text(fm[0] + n[0] + 0.03, fm[1] + n[1],
            r"$\mathbf{S}_{int}=\frac{1}{4}(\mathbf{e}_\eta\times\mathbf{e}_\zeta)$",
            fontsize=10.5, color=C_FMID)

    # angle arc between d (horizontal) and S
    ang_d = 0.0
    ang_s = np.degrees(np.arctan2(n[1], n[0]))
    arc = Arc(fm, 0.55, 0.55, angle=0, theta1=ang_s, theta2=ang_d,
              color="0.2", lw=1.4)
    ax.add_patch(arc)
    ax.text(fm[0] + 0.33, fm[1] - 0.10, r"$\theta$", fontsize=12)

    ax.text(
        0.52, -0.16,
        "仿射六面体：$\\cos\\theta="
        "\\dfrac{|\\mathbf{e}_\\xi\\cdot(\\mathbf{e}_\\eta\\times\\mathbf{e}_\\zeta)|}"
        "{\\|\\mathbf{e}_\\xi\\|\\,\\|\\mathbf{e}_\\eta\\times\\mathbf{e}_\\zeta\\|}$"
        " 与母面对应方向夹角严格相同",
        transform=ax.transAxes, ha="center", fontsize=9.8,
        bbox=dict(boxstyle="round,pad=0.4", fc="0.97", ec="0.75"),
    )

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("新增内部面的非正交角（仿射单元，2D 投影）", fontsize=11)
    save(fig, "non_orthogonality_internal_face")


# --------------------------------------------------------------------------- #
#  Fig 4: 2D trapezoid worked example
# --------------------------------------------------------------------------- #
def fig_trapezoid():
    def x(xi, eta):
        return np.array([xi, eta * (1 + xi)])

    fig, ax = plt.subplots(figsize=(6.6, 5.4))

    corners = [x(0, 0), x(1, 0), x(1, 1), x(0, 1)]
    poly = np.array(corners + [corners[0]])
    ax.plot(poly[:, 0], poly[:, 1], color=C_EDGE, lw=2.0)
    for (px, py), lab, dx, dy in [
        (x(0, 0), "$P_{00}=(0,0)$", -6, -16),
        (x(1, 0), "$P_{10}=(1,0)$", -6, -16),
        (x(1, 1), "$P_{11}=(1,2)$", 4, 2),
        (x(0, 1), "$P_{01}=(0,1)$", -64, -4),
    ]:
        ax.plot([px], [py], "o", color=C_EDGE, ms=5)
        ax.annotate(lab, (px, py), xytext=(dx, dy),
                    textcoords="offset points", fontsize=10)

    # subdivision lines xi=1/2, eta=1/2
    a1, a2 = x(0.5, 0), x(0.5, 1)
    b1, b2 = x(0, 0.5), x(1, 0.5)
    ax.plot([a1[0], a2[0]], [a1[1], a2[1]], color=C_CUT, lw=1.2, ls=(0, (4, 3)))
    ax.plot([b1[0], b2[0]], [b1[1], b2[1]], color=C_CUT, lw=1.2, ls=(0, (4, 3)))

    # parent: center, top-face midpoint, d_top, n_top
    c = np.array([0.5, 0.75])
    ftop = np.array([0.5, 1.5])
    ax.plot([c[0]], [c[1]], "o", color=C_CHILD, ms=6)
    ax.annotate("$\\mathbf{c}=(0.5,0.75)$", (c[0], c[1]), xytext=(8, -14),
                textcoords="offset points", fontsize=10, color=C_CHILD)
    ax.plot([ftop[0]], [ftop[1]], "s", color=C_FMID, ms=6)
    ax.annotate("$\\mathbf{f}_{top}=(0.5,1.5)$", (ftop[0], ftop[1]),
                xytext=(8, -2), textcoords="offset points",
                fontsize=10, color=C_FMID)
    ax.annotate("", (ftop[0], ftop[1]), (c[0], c[1]),
                arrowprops=dict(arrowstyle="-|>", color=C_CMID, lw=1.8))
    ax.text(0.53, 1.12, "$\\mathbf{d}_{top}$", fontsize=10, color=C_CMID)
    nlen = 0.4
    ax.annotate("", (ftop[0] - nlen / np.sqrt(2), ftop[1] + nlen / np.sqrt(2)),
                (ftop[0], ftop[1]),
                arrowprops=dict(arrowstyle="-|>", color=C_FMID, lw=1.8))
    ax.text(ftop[0] - 0.40, ftop[1] + 0.24, "$\\mathbf{n}_{top}=(-1,1)$",
            fontsize=10, color=C_FMID)
    ax.text(0.30, 1.62, "$\\theta_{top}=45^\\circ$", fontsize=11, color="0.2")

    # internal cut edge xi=1/2, eta in [0,1/2]
    e1, e2 = x(0.5, 0), x(0.5, 0.5)
    ax.plot([e1[0], e2[0]], [e1[1], e2[1]], color="0.15", lw=2.6)
    fmid = 0.5 * (e1 + e2)
    ax.annotate("", (fmid[0] + 0.38, fmid[1]), (fmid[0], fmid[1]),
                arrowprops=dict(arrowstyle="-|>", color=C_FMID, lw=1.8))
    ax.text(fmid[0] + 0.05, fmid[1] + 0.10, "$\\mathbf{n}_{int}=(1,0)$",
            fontsize=10, color=C_FMID)
    cL = x(0.25, 0.25)
    cR = x(0.75, 0.25)
    ax.plot([cL[0]], [cL[1]], "o", color=C_CHILD, ms=6)
    ax.plot([cR[0]], [cR[1]], "o", color="C7", ms=6)
    ax.annotate("$\\mathbf{c}_{left}$", (cL[0], cL[1]), xytext=(-14, -16),
                textcoords="offset points", fontsize=9.5, color=C_CHILD)
    ax.annotate("$\\mathbf{c}_{right}$", (cR[0], cR[1]), xytext=(0, -16),
                textcoords="offset points", fontsize=9.5, color="C7")
    ax.annotate("", (cR[0], cR[1]), (cL[0], cL[1]),
                arrowprops=dict(arrowstyle="-|>", color=C_CMID, lw=1.8))
    ax.text(0.44, 0.40, "$\\mathbf{d}_{int}=(0.5,0.125)$", fontsize=10,
            color=C_CMID)
    ax.text(0.60, 0.16, "$\\theta_{int}\\approx14.0^\\circ$", fontsize=11,
            color="0.2")

    ax.set_xlim(-0.25, 1.45)
    ax.set_ylim(-0.25, 2.2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("2D 验算：$\\mathbf{x}(\\xi,\\eta)=(\\xi,\\,\\eta(1+\\xi))$",
                 fontsize=11)
    save(fig, "trapezoid_2d_check")


# --------------------------------------------------------------------------- #
#  Fig 5 (new): sheared parallelepiped numeric example, theta = 26.57 deg
# --------------------------------------------------------------------------- #
def fig_shear():
    e_xi = np.array([1.0, 0.0, 0.0])
    e_eta = np.array([0.5, 1.0, 0.0])
    e_zeta = np.array([0.0, 0.0, 1.0])

    def P(xi, eta, zeta):
        return xi * e_xi + eta * e_eta + zeta * e_zeta

    V = np.array([P(i, j, k) for i, j, k in
                  [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                   (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]])
    AZ, EL = 24.0, 20.0
    V2, _ = project(V, AZ, EL)
    hid = hidden_edges(face_normals(V), view_dir(AZ, EL))

    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    for e in EDGES:
        draw_edge(ax, V2, *e, color=C_EDGE, lw=1.4,
                  ls=(0, (5, 3)) if hid(e) else "-")

    # internal face xi = 1/2
    fc = np.array([P(0.5, 0, 0), P(0.5, 1, 0), P(0.5, 1, 1), P(0.5, 0, 1)])
    fc2, _ = project(fc, AZ, EL)
    ax.add_patch(Polygon(fc2, closed=True, alpha=0.25, color=C_FMID, zorder=3))
    for i in range(4):
        seg = np.array([fc2[i], fc2[(i + 1) % 4]])
        ax.plot(seg[:, 0], seg[:, 1], color=C_FMID, lw=2.0, zorder=4)
    ax.annotate("内部面 $\\xi=1/2$", (fc2[1, 0], fc2[1, 1]),
                xytext=(4, -12), textcoords="offset points",
                fontsize=10, color=C_FMID)

    # centres + d_int
    cL = P(0.25, 0.5, 0.5)
    cR = P(0.75, 0.5, 0.5)
    cL2, _ = project(cL, AZ, EL)
    cR2, _ = project(cR, AZ, EL)
    ax.plot([cL2[0, 0]], [cL2[0, 1]], "o", color=C_CHILD, ms=6, zorder=9)
    ax.plot([cR2[0, 0]], [cR2[0, 1]], "o", color="C7", ms=6, zorder=9)
    ax.add_patch(FancyArrowPatch(cL2[0], cR2[0], arrowstyle="-|>",
                                 mutation_scale=14, color=C_CMID, lw=2.0,
                                 zorder=9))
    dmid = 0.5 * (cL2[0] + cR2[0])
    ax.text(dmid[0], dmid[1] + 0.05,
            r"$\mathbf{d}_{int}=\frac{1}{2}\mathbf{e}_\xi$",
            fontsize=10.5, color=C_CMID, ha="center")

    # S_int at face centre
    fm = P(0.5, 0.5, 0.5)
    s_dir = np.cross(e_eta, e_zeta)
    s_dir = s_dir / np.linalg.norm(s_dir) * 0.55
    arrow3d(ax, fm, fm + s_dir, C_FMID, AZ, EL, lw=2.2)
    tip, _ = project(fm + s_dir, AZ, EL)
    ax.text(tip[0, 0] + 0.04, tip[0, 1],
            r"$\mathbf{S}_{int}=\frac{1}{4}(\mathbf{e}_\eta\times\mathbf{e}_\zeta)$",
            fontsize=10.5, color=C_FMID)

    ax.text(
        0.5, -0.10,
        "$\\mathbf{e}_\\xi=(1,0,0),\\ \\mathbf{e}_\\eta=(0.5,1,0),"
        "\\ \\mathbf{e}_\\zeta=(0,0,1)$\n"
        "$\\cos\\theta=\\dfrac{0.125}{\\sqrt{0.25^2+0.125^2}\\times0.5}"
        "=0.8944\\ \\Rightarrow\\ \\theta=26.57^\\circ$"
        "（与母面严格相同）",
        transform=ax.transAxes, ha="center", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", fc="0.97", ec="0.75"),
    )

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("数值例：剪切平行六面体的内部面非正交角", fontsize=11)
    save(fig, "affine_shear_example")


# --------------------------------------------------------------------------- #
#  Fig 6 (new): when can child quality be worse than parent?  (three panels)
#  Panel A/B use a verified 3D counterexample (corner detJ>0, interior detJ<0)
#  found by randomized search; the corner coordinates are baked in below.
# --------------------------------------------------------------------------- #
HEX_BAD = np.array(
    [
        [0.1908, 0.5124, 0.2950],
        [1.2917, -0.6501, 0.1057],
        [0.5740, 0.5864, -0.4888],
        [-0.2611, 1.2868, 0.5038],
        [-0.4047, -0.6102, 0.5562],
        [0.9320, -0.3775, 1.0433],
        [1.2576, 1.3079, 1.0134],
        [0.2598, 1.4800, 0.7542],
    ]
)  # corner order: P000,P100,P110,P010,P001,P101,P111,P011

_IJK = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def _tri_map(P, xi, eta, zeta):
    out = np.zeros(3)
    for (i, j, k), p in zip(_IJK, P):
        w = (xi if i else 1 - xi) * (eta if j else 1 - eta) * (zeta if k else 1 - zeta)
        out += w * p
    return out


def _tri_detJ(P, xi, eta, zeta):
    cols = np.zeros((3, 3))
    for (i, j, k), p in zip(_IJK, P):
        wi = xi if i else 1 - xi
        wj = eta if j else 1 - eta
        wk = zeta if k else 1 - zeta
        cols[:, 0] += (1 if i else -1) * wj * wk * p
        cols[:, 1] += (1 if j else -1) * wi * wk * p
        cols[:, 2] += (1 if k else -1) * wi * wj * p
    return np.linalg.det(cols)


def _child_corners(P, a, b, c):
    return np.array([_tri_map(P, (a + i) / 2, (b + j) / 2, (c + k) / 2)
                     for (i, j, k) in _IJK])


def fig_failure_modes():
    P = HEX_BAD
    fig = plt.figure(figsize=(15.5, 5.0))
    axA = fig.add_subplot(1, 3, 1)
    axB = fig.add_subplot(1, 3, 2)
    axC = fig.add_subplot(1, 3, 3)

    # ---- Panel A: 3D warped parent, two inverted children highlighted ------
    AZ, EL = 26.0, 16.0
    V2, _ = project(P, AZ, EL)
    hid = hidden_edges(face_normals(P), view_dir(AZ, EL))
    for e in EDGES:
        draw_edge(axA, V2, *e, color=C_EDGE, lw=1.6,
                  ls=(0, (5, 3)) if hid(e) else "-")

    for (a, b, c) in [(0, 0, 0), (1, 0, 0)]:
        Ck = _child_corners(P, a, b, c)
        Ck2, _ = project(Ck, AZ, EL)
        for e in EDGES:
            draw_edge(axA, Ck2, *e, color="C3", lw=2.2, zorder=7)

    M = _tri_map(P, 0.5, 0.0, 0.0)
    M2, _ = project(M, AZ, EL)
    axA.plot([M2[0, 0]], [M2[0, 1]], "o", color="C3", ms=9, zorder=9)
    axA.annotate("边中点处 $\\det\\mathbf{J}=-0.08$\n(两个子单元共享此角)",
                 (M2[0, 0], M2[0, 1]), xytext=(14, -26),
                 textcoords="offset points", fontsize=9.5, color="C3",
                 arrowprops=dict(arrowstyle="->", color="C3", lw=1.2))
    axA.annotate("$P_{000}$", (V2[0, 0], V2[0, 1]), xytext=(-6, -16),
                 textcoords="offset points", fontsize=10, color=C_EDGE)
    axA.annotate("$P_{111}$", (V2[6, 0], V2[6, 1]), xytext=(2, 8),
                 textcoords="offset points", fontsize=10, color=C_EDGE)
    axA.set_title("机制 1（3D）：角点全合法、内部有负 $\\det\\mathbf{J}$ 口袋\n"
                  "红色为出生即翻转的两个子单元", fontsize=10)
    axA.set_aspect("equal")
    axA.axis("off")

    # ---- Panel B: bar chart, parent corner detJ vs child min detJ ----------
    corner_dets = [_tri_detJ(P, i, j, k) for (i, j, k) in _IJK]
    child_mins = []
    for a in [0, 1]:
        for b in [0, 1]:
            for c in [0, 1]:
                vals = [_tri_detJ(P, a / 2 + u / 2, b / 2 + v / 2, c / 2 + w / 2)
                        for u in np.linspace(0, 1, 9)
                        for v in np.linspace(0, 1, 9)
                        for w in np.linspace(0, 1, 9)]
                child_mins.append(min(vals))

    x = np.arange(8)
    axB.bar(x - 0.21, corner_dets, width=0.4, color="C2", alpha=0.85,
            label="母单元 8 角点 $\\det\\mathbf{J}$（全 $>0$）")
    colors = ["C3" if v < 0 else "0.55" for v in child_mins]
    axB.bar(x + 0.21, child_mins, width=0.4, color=colors, alpha=0.9,
            label="8 个子单元各自 $\\min\\det\\mathbf{J}$")
    axB.axhline(0, color="0.2", lw=1.0)
    axB.set_xticks(x)
    axB.set_xticklabels([str(i + 1) for i in x], fontsize=9)
    axB.set_xlabel("角点 / 子单元编号", fontsize=9.5)
    axB.set_ylabel("$\\det\\mathbf{J}$", fontsize=10)
    axB.legend(fontsize=8.5, frameon=False, loc="upper right")
    axB.set_title("母单元角点全正（体积 $\\approx+1.19$），\n"
                  "但 2 个子单元 $\\min\\det\\mathbf{J}<0$", fontsize=10)
    axB.tick_params(labelsize=9)

    # ---- Panel C: 2:1 transition face on a perfect Cartesian grid ----------
    fine1 = np.array([[0, 0], [0.5, 0], [0.5, 0.5], [0, 0.5]])
    fine2 = fine1 + np.array([0, 0.5])
    coarse = np.array([[0.5, 0], [1.5, 0], [1.5, 1], [0.5, 1]])
    for cell, col, lw in [(fine1, C_CHILD, 1.6), (fine2, C_CHILD, 1.6),
                          (coarse, "C7", 2.2)]:
        poly = np.array(list(cell) + [cell[0]])
        axC.plot(poly[:, 0], poly[:, 1], color=col, lw=lw)

    cF = np.array([0.25, 0.25])
    cC = np.array([1.0, 0.5])
    fC = np.array([0.5, 0.25])
    fP = np.array([0.5, 1.0 / 3.0])  # connector meets face plane at y=1/3

    axC.plot([cF[0]], [cF[1]], "o", color=C_CHILD, ms=6)
    axC.plot([cC[0]], [cC[1]], "o", color="C7", ms=6)
    axC.annotate("细单元中心", (cF[0], cF[1]), xytext=(-16, -16),
                 textcoords="offset points", fontsize=9, color=C_CHILD)
    axC.annotate("粗单元中心", (cC[0], cC[1]), xytext=(-16, -16),
                 textcoords="offset points", fontsize=9, color="C7")
    axC.annotate("", (cC[0], cC[1]), (cF[0], cF[1]),
                 arrowprops=dict(arrowstyle="-|>", color=C_CMID, lw=2.0))
    axC.text(0.62, 0.30, "$\\mathbf{d}$", fontsize=11, color=C_CMID)
    axC.plot([fC[0]], [fC[1]], "s", color=C_FMID, ms=7)
    axC.annotate("面心 $(0.5,\\,0.25)$", (fC[0], fC[1]), xytext=(-92, -14),
                 textcoords="offset points", fontsize=9, color=C_FMID)
    axC.plot([fP[0]], [fP[1]], "x", color="C3", ms=9, mew=2.2)
    axC.annotate("中心连线实际交点\n$(0.5,\\,0.333)$", (fP[0], fP[1]),
                 xytext=(8, 6), textcoords="offset points",
                 fontsize=9, color="C3")
    axC.plot([fC[0], fP[0]], [fC[1], fP[1]], color="C3", lw=2.4)
    axC.annotate("", (fC[0] + 0.28, fC[1]), (fC[0], fC[1]),
                 arrowprops=dict(arrowstyle="-|>", color=C_FMID, lw=1.8))
    axC.text(fC[0] + 0.20, fC[1] - 0.09, "$\\mathbf{S}_f$", fontsize=10,
             color=C_FMID)
    axC.text(
        0.5, -0.16,
        "完美 Cartesian 网格的 2:1 过渡面：\n"
        "$\\theta\\approx18.4^\\circ$，偏斜度 $\\chi\\approx0.105$",
        transform=axC.transAxes, ha="center", fontsize=9.3,
        bbox=dict(boxstyle="round,pad=0.4", fc="0.97", ec="0.75"),
    )
    axC.set_title("机制 2：2:1 过渡面的固有非正交与偏斜", fontsize=10)
    axC.set_xlim(-0.15, 1.7)
    axC.set_ylim(-0.2, 1.15)
    axC.set_aspect("equal")
    axC.axis("off")

    fig.suptitle("什么情况下子单元质量会低于母单元？", fontsize=12)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save(fig, "child_quality_failure_modes")


if __name__ == "__main__":
    fig_subdomain()
    fig_jacobian()
    fig_nonorth()
    fig_trapezoid()
    fig_shear()
    fig_failure_modes()
    print("QA PNGs in", QA_DIR)
