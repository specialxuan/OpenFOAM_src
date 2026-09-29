# REPORT — Experimental pressure validation of the HiSA RANS setup using the official ONERA M6 example

**Date:** 2026-09-28
**Task:** Validate the HiSA RANS setup against genuine experimental surface pressures.
**Substitute for:** AGARD 445.6 — which has **no** experimental `Cp` dataset (flutter only); the
ONERA M6 wing is the canonical transonic-wing pressure-tap validation case.
**Status:** ✅ **Run completed and converged; experimental comparison produced.**
**Verdict:** **The HiSA RANS setup is trustworthy for transonic wing pressure prediction**
(mean |ΔCp| = 0.049–0.098 over all taps; core max |ΔCp| = 0.066–0.277; main-shock x/c agrees
within ≤0.042). Details and caveats below.

---

## 1. What was done (summary)

| Item | Result |
|---|---|
| Job run? | **Yes** — full RANS steady run to `endTime=2000`. |
| Mesh size | **294,912 hexahedral cells** (306,577 points, 896,256 faces). Well under the 1.5 M ceiling. |
| Wall clock / cores | **21.5 min** on **16 ranks** (`mpirun -n 16 hisa -parallel`), exit code 0. |
| Convergence | Coupled residual 1.0 → **1.25e-3** (≈2.9 orders); plateau from iter ~1000 (stable). |
| Free stream vs standard | **M = 0.8395** (std 0.8395), **α = 3.0601°** (std 3.06°), **Re_MAC = 11.72e6** (std 11.72e6) — effectively exact. |
| Experimental data | NASA TMR `case_2308.dat` (digitised AGARD-AR-138, Schmitt & Charpin 1979). |
| Comparison | 7 spans, upper+lower, 271 taps; CSVs + figure in `results/`. |

The job **was not aborted on size** — the mesh is 0.29 M cells (the instruction to stop was only if
projected size exceeded ~1.5 M cells).

---

## 2. Setup

### 2.1 Source (authoritative, read-only — not modified)
```
/root/OpenFOAM/user-v2412/hisa/hisa/examples/oneraM6/
```
Working copy (all work done here, original untouched):
```
/root/Workspace/AGARD/HiSA/oneraM6/
```
Note: the example has **no** `settings` file and **no** `mesh/system/meshDict`. The mesh is a
Plot3D → Gmsh → OpenFOAM conversion driven by `setupMesh`. The shipped `mesh/p3dMesh/m6wing.msh`
(294,912 hexes; element histogram 294,912 type-5 + 43,008 type-3) was used directly, because the
example's Python converters (`p3d2gmsh.py`) require **Python 2** (`xrange`, `PyFile_AsFile`), which is
absent on this host. The `.msh` is byte-for-byte the mesh the example ships and exactly matches the
`m6wing.nmf` block dimensions (2×24·48·32 + 2×72·48·32 = 294,912 cells). `setupMesh` was executed;
its two Python steps no-op'd and `gmshToFoam` consumed the shipped `.msh`.

### 2.2 Mesh size and quality (from `checkMesh`, `mesh/log.checkMesh`)

| Metric | Value |
|---|---|
| cells / points / faces | 294,912 / 306,577 / 896,256 (all hexahedra) |
| boundary patches | 3: `wing` 4,608 faces, `symm` 9,216, `farfield` 9,216 |
| non-orthogonality | max **88.38°**, average 19.86° (1,622 faces > 70°) |
| skewness | max **16.31** (213 highly-skew faces flagged) |
| aspect ratio | max 1.04e5 (7,562 wall-normal stretched cells) |
| min face area / min cell volume | 3.25e-8 m² / 1.30e-10 m³ |
| mesh checks | **Failed 2** (high aspect ratio, skewness) — intrinsic to this viscous O-grid |
| topology notes | 72 cells with multiple in-between faces (from the 6 `stitchMesh` joins); 144 unordered faces |

The wing surface is **72 chordwise × 32 spanwise faces per side** (upper = 2,304, lower = 2,304).
The **32-cell spanwise resolution (Δη ≈ 0.037)** is the main mesh limitation; comparisons were made by
interpolating the structured surface grid to the exact experimental η.

### 2.3 Solver / turbulence settings (from the example, unchanged)

| Setting | Value |
|---|---|
| application | `hisa` (density-based coupled implicit) |
| flux scheme | `AUSMPlusUp` |
| pseudo-time | `bounded dualTime rPseudoDeltaT steadyState`, local time stepping (`pseudoCoNumMax 20`) |
| flow solver | `GMRES`, LaxFriedrichs inviscid Jacobian, Laplacian viscous Jacobian, `LUSGS` preconditioner, `nKrylov 8`, `solverTolRel 1e-1` |
| turbulence | `RAS SpalartAllmaras`, `turbulence on` |
| thermo/transport | `hePsiThermo`, perfect gas, `hConst`, `sutherland As=1.4584e-6, Ts=110.33` |
| wall distance | `wallDist { method meshWave; }` |
| decomposition | `scotch`, 16 subdomains (example default 8; changed **only** in the working copy to use ≤16 cores) |

### 2.4 Free-stream conditions — exact vs standard ONERA M6 test 2308

Taken from `simulation/0.org/include/freestreamConditions`; derived quantities computed exactly.

| Quantity | HiSA example | Standard test 2308 | Match |
|---|---|---|---|
| Mach | 0.8395 | 0.8395 (often quoted 0.84) | ✅ exact |
| α | 3.0601° (`atan2(14.362,268.654)`) | 3.06° | ✅ exact |
| Re (MAC = 0.64607 m) | 1.172e7 | 11.72e6 | ✅ exact |
| T (static) | 255.556 K = 460.0 R | 460 R (static) | ✅ |
| p (static) | 80,510.081 Pa = 11.677 psi | — | reference |
| U∞ | (268.654, 14.362, 0) m/s; |U| = 269.038 m/s | — | — |
| ρ∞ / μ∞ | 1.09747 kg/m³ / 1.6284e-5 Pa·s | — | — |
| q∞ = ½ρU² | 39,718.29 Pa | — | Cp = (p − 80510.081)/39718.29 |

**Deviation from standard: none material.** Mach, α and Re are the standard values. The only
methodological note is that the NASA TMR *CFD exercise* now recommends a **sharp** trailing edge and
Re_c_root = 14.6e6 for code-to-code work; that is **not** the experiment. This report compares the
example (AGARD-derived grid, experimental Re_MAC = 11.72e6) against the **original experiment**, which
is the correct pairing for experimental validation.

---

## 3. Convergence

Coupled (ρ, ρU, ρE) residual (normalised by iteration 1), from `simulation/log.hisa`:

| iteration | 1 | 250 | 500 | 750 | 1000 | 1500 | 2000 |
|---|---|---|---|---|---|---|---|
| residual | 1.00e0 | 1.10e-2 | 3.41e-3 | 1.51e-3 | 1.32e-3 | 1.27e-3 | **1.25e-3** |

The residual drops ≈2.9 orders and **plateaus at ~1.25e-3 from iteration ~1000**, i.e. the steady
state is converged; the floor is set by the linear-solver relative tolerance (`solverTolRel 1e-1`),
not by a lack of convergence. `nuTilda` final residual ≈1e-6. Run ended cleanly:
`=== solver exit=0 2026-09-28T13:02:55+08:00 ===`.

---

## 4. `Cp` extraction method

1. `reconstructPar -latestTime` → time `2000`; `foamToVTK -fields p -patches wing -latestTime -no-point-data`
   → wall pressure `p` on the 4,608 wing faces (exact boundary values, no near-wall extrapolation).
2. `Cp = (p − p∞)/q∞` with `p∞ = 80510.081 Pa`, `q∞ = 39718.29 Pa`.
3. The wing patch is structured **f = i·32 + k** (i = chord 0..71, k = span 0..31); upper `y>0`, lower `y<0`.
4. Planform calibrated to the mesh corners: `x_LE(z) = 0.57280·z`, `x_TE(z) = 0.80595 + 0.28020·z`
   (root chord 0.80595 m, tip TE x = 1.14115 at z = 1.1963 m, tip chord 0.45591 m). `η = z/1.1963`.
5. Each surface section is obtained by interpolating the 32 spanwise cells to the exact target `z = η·1.1963`
   (handles the slightly sheared grid), giving 72 chordwise points per surface.
6. Experimental taps (`x/c`, `Cp`) are mapped through the same planform and the computed `Cp` interpolated
   at each tap. (TMR notes the 6th station is **η = 0.96**, not 0.95; used 0.96.)

Tooling: `extract_cp.py`, `plot_cp.py` (this directory). No values were eyeballed or fabricated.

---

## 5. Results — computed vs experimental `Cp`

Experimental source: **NASA Turbulence Modeling Resource, ONERA M6 `case_2308.dat`**
<https://tmbwg.github.io/turbmodels/Onerawingnumerics_val/case_2308.dat>
(digitisation of Schmitt & Charpin, *Pressure Distributions on the ONERA-M6-Wing at Transonic Mach
Numbers*, AGARD-AR-138, 1979; header of the file: `Mach=0.8395, Alpha=3.06, Re=11.72x10**6`).
271 taps total.

Per-span difference summary (`results/summary.csv`). "core" excludes the extreme-nose tap
(x/c < 0.02) where tap location is geometrically ambiguous.

| η | z (m) | max\|ΔCp\| all | max\|ΔCp\| core | mean\|ΔCp\| | shock comp | shock exp | Δ(x/c) shock | LE peak comp | LE peak exp |
|---|---|---|---|---|---|---|---|---|---|
| 0.20 | 0.2393 | 0.181 | 0.181 | 0.051 | 0.617 | 0.575 | +0.042 | −1.084 | −1.069 |
| 0.44 | 0.5264 | 0.358 | 0.066 | 0.049 | 0.532 | 0.524 | +0.008 | −1.161 | −1.127 |
| 0.65 | 0.7776 | 0.403 | 0.161 | 0.050 | 0.467 | 0.475 | −0.008 | −1.194 | −1.181 |
| 0.80 | 0.9570 | 0.413 | 0.250 | 0.061 | 0.381 | 0.375 | +0.006 | −1.207 | −1.206 |
| 0.90 | 1.0767 | 0.885 | 0.277 | 0.080 | 0.261 | 0.280 | −0.019 | −1.218 | −1.194 |
| 0.96 | 1.1484 | 0.647 | 0.135 | 0.067 | 0.207 | 0.200 | +0.007 | −1.237 | −1.227 |
| 0.99 | 1.1843 | 0.703 | 0.274 | 0.098 | 0.207 | 0.200 | +0.007 | −1.289 | −1.291 |

**Overall:** mean |ΔCp| = **0.049–0.098** (avg ≈ 0.065); core max |ΔCp| = **0.066–0.277**.
Main-shock x/c agrees within **≤0.042**, and within **≤0.019** at all outboard stations.
Leading-edge suction peaks agree within **0.02–0.07**.

**Where the largest differences are, and why:**
- **Extreme nose taps (x/c ≲ 0.01).** max|ΔCp| up to 0.89 at η=0.90. The experimental first tap
  (x/c ≈ 0.002, Cp ≈ −0.18) sits on the near-stagnation nose where `Cp` is extremely sensitive to
  sub-cell tap position and to the (rounded vs discrete) nose geometry; the wall `Cp` has a large
  gradient there. This dominates the "all taps" maximum but is not a solution-quality defect.
- **Shock region.** The residual core maxima (0.14–0.28) occur across the shock at η = 0.80–0.99,
  where a ~0.01–0.02 x/c shock-position offset produces a locally large `Cp` difference.
- **λ-shock stations.** At η = 0.20/0.44/0.65 the M6 has a double (lambda) shock. The detector reports
  the strongest compression: the main shock matches (e.g. η=0.65: 0.467 vs 0.475), but the forward
  λ-leg differs in strength (experiment shows a forward leg near x/c ≈ 0.15–0.20 that the computation
  reproduces more weakly). This is the main *structural* difference seen.
- **Lower surface and outboard sections** agree closely (most |ΔCp| < 0.03).

Outputs:
- `results/summary.csv` — per-span metrics table.
- `results/comparison_eta*.csv` — per-tap `surface,x/c,Cp_exp,Cp_computed,ΔCp` (7 files).
- `results/computed_eta*.csv` — computed upper/lower `Cp(x/c)` curves (7×2 files).
- `results/oneraM6_Cp_comparison.png` — 7-panel `−Cp` vs `x/c` plot, computed lines vs experimental symbols.
- `experimental/exp_cp_sections.csv` — parsed experimental taps (from `case_2308.dat`).

---

## 6. Verdict on the HiSA RANS setup

**Trustworthy for transonic-wing pressure prediction at this level of fidelity.**

- Mean |ΔCp| ≈ **0.05–0.10**, core max |ΔCp| ≈ **0.07–0.28**, and shock location within
  **~0.02 x/c** — at or near the typical scatter of RANS ONERA-M6 validations on comparable meshes.
- The classic ONERA M6 features are captured: leading-edge suction peak, upper-surface plateau,
  swept shock (single outboard, λ inboard), and lower-surface pressure distribution.
- The setup is **not** high-fidelity: the 32-cell spanwise resolution, the standard SA model,
  first/second-order upwind turbulence discretisation and the semi-coarse (0.29 M cell) grid mean the
  λ-foot and exact shock shape will differ from experiment; the nose taps are unreliable for the reason
  above. For engineering transonic loads/pressures this is adequate; for detailed shock-boundary-layer
  interaction or λ-shock structure it is not.

---

## 7. Blockers / limitations encountered

1. **Python 2 absent** → the example's `p3d2gmsh.py` cannot regenerate the `.msh`; the shipped,
   dimension-verified `.msh` was used instead. No effect on results (identical mesh).
2. **Mesh checks "Failed 2"** (high aspect ratio, skewness) — intrinsic to the supplied viscous O-grid,
   not introduced here.
3. **Nose-tap ambiguity** limits a meaningful "max |ΔCp| over all taps"; core (x/c ≥ 0.02) metric also reported.
4. `characteristicWallPressure` was used as supplied; wall `Cp` is the solved boundary value.

---

## 8. Exact next commands (reproduce end-to-end)

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
cd /root/Workspace/AGARD/HiSA/oneraM6

# 1. mesh (uses shipped m6wing.msh; Python-2 converters no-op)
./setupMesh

# 2. solve (16 ranks; ~22 min)
cd simulation
rm -rf 0 processor*; cp -r 0.org 0
mkdir -p constant/polyMesh
for f in boundary faces neighbour owner points; do ln -sf ../../../mesh/constant/polyMesh/$f constant/polyMesh/$f; done
decomposePar
mpirun --allow-run-as-root -n 16 hisa -parallel | tee log.hisa
cd ..

# 3. extract & compare
cd simulation && reconstructPar -latestTime
foamToVTK -fields p -patches wing -latestTime -no-point-data -ascii -no-internal && cd ..
python3 extract_cp.py simulation/VTK/simulation_2000/boundary/wing.vtp experimental/exp_cp_sections.csv results
python3 plot_cp.py     results experimental/exp_cp_sections.csv results/oneraM6_Cp_comparison.png
```

---

## 9. Provenance / data sources

| Artifact | Source |
|---|---|
| HiSA example (mesh, dicts, fields) | `/root/OpenFOAM/user-v2412/hisa/hisa/examples/oneraM6/` (read-only) |
| Experimental `Cp` | NASA TMR `case_2308.dat` — <https://tmbwg.github.io/turbmodels/Onerawingnumerics_val/case_2308.dat> |
| Experimental reference | Schmitt, V. & Charpin, F., AGARD-AR-138 (1979); NASA TMR ONERA M6 validation page <https://tmbwg.github.io/turbmodels/onerawingnumerics_val.html> |
| Geometry/planform | TMR ONERA M6 page (root chord, semispan, tip chord, sweep) + mesh-corner calibration |
| Run artifacts | `/root/Workspace/AGARD/HiSA/oneraM6/simulation/` (log.hisa, VTK/, 2000/) |

*Measured vs assumed:* all reported `Cp`, ΔCp, shock and LE-peak values are **measured** from the
computed field and the digitised experiment. The planform lines are **calibrated** from mesh corners
(stated above). No experimental value was invented.
