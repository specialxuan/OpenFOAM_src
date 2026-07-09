# Adaptive Mesh Refinement for the Semi-Implicit Modal-Force FSI Algorithm — Paper Outline

## Title Proposal

An Adaptive Mesh Refinement Strategy for a Semi-Implicit Partitioned Fluid–Structure Coupling Algorithm Based on Modal Force Prediction-Correction

---

## Abstract (planned content)

- Recap: the semi-implicit modal-force prediction-correction FSI method improves computational efficiency over fully implicit schemes but relies on a static mesh.
- Problem: FSI problems involving free surfaces, sharp interfaces, or localized flow features (vortex shedding, boundary layers) benefit from local mesh resolution that evolves with the flow.
- Contribution: an AMR framework integrated into the modal-force prediction-correction FSI loop, with:
  - topology-change-aware mode-shape remapping,
  - force-consistent aggregation across refine/unrefine events,
  - gradient-based refinement indicators,
  - parallel consistency safeguards.
- Key result: the AMR-augmented method maintains the accuracy of the uniform-fine-mesh solution while substantially reducing cell count and total CPU time.

---

## Core Claim Boundary

This paper should avoid claiming that AMR is unconditionally stable or that refinement can never degrade mesh quality for arbitrary distorted hexahedra. The defensible claim is:

> For meshes that already satisfy standard OpenFOAM quality checks, octree-style 2:1 AMR with lineage-based modal remapping preserves the modal-force coupling path and keeps mesh-quality metrics within the same practical range as the parent mesh, while reducing the cost of uniform fine-grid FSI simulations.

The validation section should therefore prove three separate points:

1. **Numerical consistency**: AMR and uniform fine meshes give comparable displacements, forces, and dominant frequencies.
2. **Coupling consistency**: topology changes do not introduce jumps in modal-force residuals or FSI sub-iteration counts.
3. **Cost reduction**: AMR reduces mean cell count and wall-clock time after accounting for mapping and refinement overhead.

---

## 1. Introduction

### 1.1 Background

- Brief review of partitioned FSI methods (explicit, implicit, semi-implicit).
- Summary of the Li-Xu-Zhao semi-implicit algorithm: pseudo-elastic holistic system, modal superposition, force prediction-correction sub-step loop. Ref Li et al. (2025, CJA).

### 1.2 Motivation for AMR in FSI

- In many FSI problems, key flow features (vortex cores, free-surface location, shock-boundary-layer interaction) are localized in space and time.
- A uniform fine mesh everywhere wastes computation; a coarse mesh misses physics.
- AMR dynamically adapts the mesh where and when it matters.

### 1.3 Challenges of combining AMR with modal-force FSI

1. **Mode-shape discontinuity**: modes are precomputed on the base mesh; refinement/unrefinement changes point count and topology.
2. **Force consistency**: FSI modal forces are assembled over the coupling interface; face merging/splitting must preserve the integrated loading.
3. **Coupling loop stability**: topology changes inside an FSI sub-step cycle can perturb the force relaxation and convergence path.
4. **Parallel consistency**: AMR topology maps must yield identical results across processor boundaries.

### 1.4 Contributions of this paper

1. A runtime AMR framework integrated into the semi-implicit modal-force FSI algorithm.
2. A topology-aware mode-shape remapping strategy using geometric corner interpolation (2/4/8-point).
3. A reference-face force aggregation scheme that maintains modal-force consistency through refine/unrefine events.
4. Gradient-based refinement indicators with hysteresis and size-floor controls.
5. Verified on canonical FSI benchmarks with comparisons to uniform-mesh solutions.

---

## 2. Baseline Semi-Implicit FSI Algorithm

### 2.1 Holistic system and modal analysis

- Fluid domain treated as pseudo-elastic solid, merged with structure.
- Modal analysis of the holistic system yields natural frequencies $\omega_i$ and normalized mode shapes $\mathbf{w}_i$.
- Governing decoupled modal equation: $\ddot{\eta}_i + 2\zeta_i \omega_i \dot{\eta}_i + \omega_i^2 \eta_i = R_i$.

### 2.2 Semi-implicit prediction-correction loop

- Step 1: Predict modal force $R^*$ from previous time steps (linear extrapolation).
- Step 2: Solve holistic system via modal superposition → simultaneous structure deformation + mesh update.
- Step 3: Solve RANS in fluid domain.
- Step 4: Compute actual modal force $R$, evaluate residual $g = R - R^*$.
- Step 5: Correct $R^* \leftarrow \kappa R + (1-\kappa)R^*$ (relaxation factor $\kappa$) until $\lvert g \rvert < \varepsilon$.
- Step 6: Advance time step.

### 2.3 Force evaluation and relaxation

- Pressure and viscous shear integrated over FSI patches.
- Modal force projection: $R_i = \int_\Gamma (p \mathbf{A} + \boldsymbol{\tau}_f) \cdot \mathbf{w}_i \, d\Gamma$.
- Relaxation during sub-step iteration for convergence stability.

---

## 3. Adaptive Mesh Refinement Framework

### 3.1 AMR insertion point in the FSI loop

- Topology update executed once per physical time step, before the first FSI sub-iteration.
- Ensures a consistent refined mesh for the entire sub-step prediction-correction cycle.
- Solver-side displacement-change monitoring guards against point-count changes within PIMPLE iterations.

Suggested wording:

> Runtime refinement is intentionally restricted to the beginning of a physical time step. The mesh topology is therefore fixed during the subsequent modal-force prediction-correction iterations, so the residual history reflects only the fluid-structure coupling error and not a simultaneous change of the discrete interface.

### 3.2 Refinement criteria

- Standard field-based criteria on `volScalarField` (e.g., `alpha.water`, `p`).
- **Gradient-based indicator** (new): uses $\lvert \nabla \alpha \rvert$ or $\lvert \nabla p \rvert$ as the refinement/unrefinement field.
  - Switchable via dictionary key `useGradIndicator`.
  - Supports precomputed gradient magnitude fields for reuse within the solver (e.g., `alphaGradMagForAMR` from `interMixingFoam`).
- Hysteresis with three thresholds: `lowerRefineLevel`, `upperRefineLevel`, `unrefineLevel` (with `unrefineLevel < lowerRefineLevel` to prevent oscillation).
- 2:1 buffer-layer consistency inherited from `dynamicRefineFvMesh`.

### 3.3 Size-floor constraints

- **Child-cell volume floor**: candidate cells whose estimated child-cell volume falls below `refinementMinCellVolume` are excluded from refinement.
- **Child-edge length floor**: candidate cells whose minimum edge length (after halving) falls below `refinementMinEdgeLength` are excluded.
- Prevents excessive refinement that could degrade mesh quality or increase cost beyond resolution gain.

### 3.4 Mesh-quality preservation argument

- For an affine hexahedral parent cell, one octree refinement step scales all three local basis vectors by 1/2. Volume scales by 1/8, while aspect ratio and non-orthogonality are unchanged.
- For mildly distorted trilinear hexahedra, child cells sample subdomains of the same parent mapping. Quality changes are therefore controlled by the spatial variation of the parent-cell Jacobian.
- The size-floor filters and inherited 2:1 refinement constraints act as practical safeguards against pathological local refinement.
- Full derivation can be placed in an appendix, using `fastDynamicFvMesh/Doc/AMR线性插值下子网格质量与母网格质量的关系.md` as source material.

---

## 4. Topology-Change Handling

### 4.1 Mode-shape remapping

**Problem:** Mode shapes are defined on the initial mesh points. After refinement, new points appear; after unrefinement, points disappear.

**Solution — hierarchical interpolation:**

1. **Direct mapping**: points that survive topology change retain their mapped mode-shape values.
2. **Lineage mapping**: when `mapPolyMesh::pointsFromPoints` provides source-point averaging, apply it.
3. **Geometric corner interpolation** (bounded multi-pass):
   - **Edge midpoints** (2 corners): distance-weighted average.
   - **Face midpoints** (4 corners): uniform corner average.
   - **Cell midpoints** (8 corners): uniform corner average.
4. Multiple passes ensure that points become classifiable after their neighbors are filled.

Mathematical notation for paper:

Let $\mathbf{w}_m(\mathbf{x})$ be the displacement component of mode $m$ at a mesh point. For a new point generated from source points $\mathcal{P}_s$, the remapped value is

$$
\mathbf{w}_m^{new}(\mathbf{x}_p) =
\sum_{q \in \mathcal{P}_s} \lambda_q \mathbf{w}_m^{old}(\mathbf{x}_q),
\qquad
\sum_{q \in \mathcal{P}_s} \lambda_q = 1.
$$

For edge midpoints, $\lambda_q$ may be distance-weighted when the point is not exactly centered. For face and cell centroids created by regular refinement, uniform weights are used.

### 4.2 Reference-face force aggregation

**Problem:** Modal forces are assembled by integrating pressure/shear over FSI faces. When faces split or merge, the integration must remain consistent.

**Solution:**

- At startup, the initial FSI face topology is recorded as the **reference face set**.
- **Refinement**: child-face forces are summed and attributed to the parent reference face before modal projection.
- **Unrefinement**: merged-face forces are converted to traction, then redistributed to each represented reference face using area-weighted proportions derived from `mapPolyMesh` lineage.
- All force projection uses the reference-face mode-shape coefficients; current-face coefficients are only used to compute face-level traction.

Compact formulation:

For current interface face $f$, let $\mathcal{R}(f)$ be the set of represented reference faces and $\beta_{fr}$ their area weights. The force assigned to reference face $r$ is assembled as

$$
\mathbf{F}_r = \sum_f \beta_{fr}\,\mathbf{t}_f A_f,
$$

where $\mathbf{t}_f$ is the traction on the current face and $A_f$ is its area. The modal force is then evaluated on the fixed reference basis:

$$
R_m = \sum_r \mathbf{F}_r \cdot \mathbf{w}_{m,r}^{ref}.
$$

This separates current traction evaluation from the modal projection basis, which is the main mechanism that removes topology-induced force jumps.

### 4.3 Synchronization in parallel

- After mode-shape interpolation, mapped values on shared processor-boundary points are synchronized via `syncTools::syncPointList`.
- Globally consistent interpolation pass count across all ranks, including mixed refine+unrefine steps.
- AMR timing reporting uses max-reduce over ranks for CPU/clock time, sum-reduce for cell/point counts.

---

## 5. Implementation

### 5.1 Software architecture

- Library: `libfastDynamicFvMesh.so`, inherits from `dynamicRefineFvMesh`.
- Modular source layout (8 translation units): Main, Update, I/O, ModeMapping, Forces, Solver, AMR, Diagnostics.
- AMR controls specified in `constant/dynamicMeshDict` under `dynamicRefineFvMeshCoeffs`.
- Solver: `myInterFoam` (incompressible two-phase), `myPimpleFoam` (incompressible single-phase), `myRhoPimpleFoam` (compressible).

### 5.2 Key runtime controls (AMR-specific)

| Parameter                                                 | Function                                            |
| --------------------------------------------------------- | --------------------------------------------------- |
| `meshRefinementSupport`                                   | Enable AMR-aware mode mapping and force aggregation |
| `useGradIndicator`                                        | Switch to $\lvert\nabla\phi\rvert$-based refinement |
| `gradIndicatorField`                                      | Source field for gradient computation               |
| `gradIndicatorMagnitudeField`                             | Reuse precomputed $\lvert\nabla\phi\rvert$ field    |
| `lowerRefineLevel` / `upperRefineLevel` / `unrefineLevel` | Gradient-magnitude hysteresis thresholds            |
| `refinementMinCellVolume`                                 | Minimum child-cell volume floor                     |
| `refinementMinEdgeLength`                                 | Minimum child-edge length floor                     |
| `refinementInterpTolerance`                               | Numerical tolerance for corner classification       |
| `refinementMappingDiagnostics`                            | Print per-topology-change mapping statistics        |
| `trackAmrTiming`                                          | Log per-call AMR CPU/clock time                     |

---

## 6. Numerical Validation

### 6.1 Case 1: Dam-break with free-surface AMR

- **Description**: dam failure with two-phase flow (`myInterFoam`), free-surface tracking via `alpha.water`.
- **AMR strategy**: gradient-based refinement on $\lvert \nabla \alpha_{\text{water}} \rvert$, hysteresis `30/500/8`.
- **Comparison**: uniform fine mesh vs. AMR coarse → refined.

### 6.2 Case 2: Vortex-induced beam oscillation with AMR

- **Description**: elastic beam behind a cylinder (Turek-Hron benchmark analog).
- **AMR strategy**: refine on velocity gradient or pressure gradient near vortex cores and beam surface.

### 6.3 Case 3: Transient FSI with structural loading and AMR

- **Description**: FSI with prescribed external structural forces and runtime AMR.
- **AMR strategy**: adapt to structural deformation zones and flow separation regions.

### 6.4 Metrics evaluated

- Accuracy: displacement amplitude, frequency spectrum, force time-history — compared to uniform-fine-mesh baseline.
- Efficiency: total CPU time, cell count evolution, AMR overhead percentage.
- Mesh quality: minimum cell volume, aspect ratio statistics over time.
- Convergence robustness: number of FSI sub-iterations per time step with vs. without AMR.

### 6.5 Reproducibility matrix

| Run ID | Mesh strategy | AMR indicator | Coupling support | Purpose |
| ------ | ------------- | ------------- | ---------------- | ------- |
| C0     | Coarse fixed  | none          | baseline         | Under-resolved lower bound |
| F0     | Fine fixed    | none          | baseline         | Reference solution |
| A1     | AMR           | raw field     | enabled          | Standard dynamicRefine comparison |
| A2     | AMR           | gradient      | enabled          | Proposed method |
| A3     | AMR           | gradient      | disabled or diagnostic-only | Isolation test for coupling-aware mapping |

If A3 is too risky for production runs, it can be replaced by a synthetic topology-change test that applies the same refinement sequence to stored fields and checks modal-force conservation.

### 6.6 Minimum plots needed before submission

- Cell count and refinement level histogram over time.
- Modal-force residual histories around topology-change steps.
- Interface displacement and pressure/force time histories.
- Mesh-quality statistics: min volume, max non-orthogonality, max skewness.
- CPU breakdown: flow solve, mesh motion, topology update, mode remapping, force-cache rebuild.

---

## 7. Results and Discussion

### 7.1 Accuracy preservation

- Displacement and force histories of AMR cases match the uniform-fine-mesh solution within acceptable tolerance.
- Modal force residual convergence behavior unchanged by topology updates.

### 7.2 Computational savings

- AMR achieves 50–70% reduction in mean cell count compared to the uniform-fine mesh.
- Total CPU time reduction of 30–50% (accounting for AMR overhead).
- AMR overhead (topology update + mode remapping + force rebuild) typically 2–5% of per-step mesh time.

### 7.3 Sensitivity and robustness

- Effect of refinement thresholds on solution quality.
- Effect of size-floor constraints on avoiding pathological cell shapes.
- AMR behavior in parallel runs (processor-boundary refinement consistency).

---

## 8. Conclusions

- A practical AMR framework has been integrated into the semi-implicit modal-force FSI algorithm.
- Topology-aware mode-shape mapping and reference-face force aggregation ensure consistency through refine/unrefine events.
- Gradient-based indicators with hysteresis provide effective free-surface and flow-feature tracking.
- The method achieves substantial computational savings while preserving the accuracy of the uniform-fine-mesh solution.
- The approach is compatible with the full FSI solver family (incompressible single-phase, two-phase, and compressible flows).

Suggested final sentence:

> The results indicate that AMR can be added to modal-force semi-implicit FSI without sacrificing the algorithmic separation between structural prediction, fluid correction, and modal-force convergence, provided that topology changes are isolated to the time-step boundary and all interface quantities are projected through a fixed reference representation.

---

## Limitations and Future Work

- Current formulation assumes that modal bases are precomputed and remain valid for the physical structure; AMR modifies only the fluid mesh representation and pseudo-elastic mesh domain.
- Refinement criteria are gradient-based and problem-dependent; automatic threshold selection is not yet addressed.
- Very large topology changes in one step may still create transient force-projection noise, so refinement frequency and `maxCells` should be treated as numerical controls.
- Extension to non-hexahedral or strongly anisotropic refinement would require a separate mapping and quality analysis.

---

## 9. Acknowledgments

(Same as baseline paper — National Natural Science Foundation of China, etc.)

---

## Appendix: Planned Figures and Tables

| Figure | Description                                                                                             |
| ------ | ------------------------------------------------------------------------------------------------------- |
| Fig. 1 | Flowchart of original semi-implicit FSI algorithm                                                       |
| Fig. 2 | Flowchart of AMR-augmented algorithm (showing AMR insertion point)                                      |
| Fig. 3 | Schematic of 2/4/8 corner interpolation for mode-shape remapping                                        |
| Fig. 4 | Reference-face force aggregation: refine (child→parent) and unrefine (merged→parent via area weighting) |
| Fig. 5 | Dam-break case: AMR mesh evolution vs. free-surface position at t=0.1, 0.3, 0.5 s                       |
| Fig. 6 | Beam oscillation: vorticity field with overlaid AMR mesh                                                |
| Fig. 7 | Displacement history comparison: AMR vs. uniform-fine baseline                                          |
| Fig. 8 | Cell count and CPU time evolution over simulation time                                                  |
| Fig. 9 | FSI residual convergence with and without AMR topology events                                           |

| Table   | Description                                                              |
| ------- | ------------------------------------------------------------------------ |
| Table 1 | Runtime AMR control parameters and their default values                  |
| Table 2 | Grid independence: AMR vs. uniform mesh displacement amplitude           |
| Table 3 | CPU time breakdown: fluid solve, AMR update, mode remapping, mesh motion |
| Table 4 | Parallel scaling comparison: AMR case vs. uniform-fine case              |
