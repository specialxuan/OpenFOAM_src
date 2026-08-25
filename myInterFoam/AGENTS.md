# myInterFoam Knowledge

## OVERVIEW

Two-phase VOF FSI solver family with topology-changing mesh support, overset variant, and nested three-phase child solver.

## STRUCTURE

| Path | Purpose |
|---|---|
| `myInterFoam.C` | Custom two-phase VOF FSI solver main |
| `pimpleLoopBody.H` | Header-included alpha, mesh, momentum, pressure sequence |
| `alphaControls.H`, `alphaEqn*.H` | Phase-fraction transport snippets |
| `overInterDyMFoam/` | Separate VOF overset/dynamic executable |
| `interMixingFoam/` | Nested three-phase solver and local mixture models |
| `Make/files`, `Make/options` | Top-level executable target and two-phase deps |

## WHERE TO LOOK

| Task | Location | Notes |
|---|---|---|
| FSI loop behavior | `myInterFoam.C` | Uses `fsiCoupling.H`; topology change must not break point-diff diagnostics |
| VOF solve order | `pimpleLoopBody.H`, `alphaEqn*.H` | Preserve alpha/mesh/flux/pressure ordering |
| Mesh compression flux | `pimpleLoopBody.H` | Clear previous compression flux after topology changes |
| Overset VOF | `overInterDyMFoam/` | Own `Make/`, fringe masks, and overset pressure correction path |
| Three-phase variant | `interMixingFoam/` | Has its own AGENTS.md and local models |
| Build deps | `Make/options` | Includes `../fsiCoupling`; separate from child `interMixingFoam/Make/options` |

## CONVENTIONS

- Use `wmake myInterFoam` from repo root for the top-level two-phase executable.
- Keep FSI convergence force-residual based and displacement change diagnostic-only.
- Preserve topology-change safeguards before subtracting point fields or reusing fluxes.
- Keep `interMixingFoam` self-contained; do not move three-phase-only models into the two-phase parent unless all variants need them.
- Update overset and three-phase `Make/` files separately from the parent solver.

## ANTI-PATTERNS

- Do not reuse previous time-step mesh compression flux after topology changes.
- Do not assume top-level `myInterFoam/Make/options` covers `interMixingFoam`.
- Do not delete local three-phase model sources as duplicate-looking code.
- Do not bypass overset fringe masks or topology guards in `overInterDyMFoam`.
