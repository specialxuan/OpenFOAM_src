# myRhoPimpleFoam Knowledge

## OVERVIEW

Compressible PIMPLE solver family with FSI convergence wiring and dynamic-mesh variants.

## STRUCTURE

| Path | Purpose |
|---|---|
| `myRhoPimpleFoam.C` | Custom compressible FSI solver main |
| `rhoPimpleFoam.C` | Baseline/reference copy kept in-tree |
| `UEqn.H`, `EEqn.H`, `pEqn.H`, `pcEqn.H` | Equation snippets for coupled flow solve |
| `createFields.H`, `createFieldRefs.H` | Field setup and references |
| `correctPhi.H`, `setRDeltaT.H` | Dynamic-mesh/LTS helpers |
| `overRhoPimpleDyMFoam/` | Compressible overset/dynamic variant |
| `Make/files`, `Make/options` | Custom executable target and compressible deps |

## WHERE TO LOOK

| Task | Location | Notes |
|---|---|---|
| FSI loop behavior | `myRhoPimpleFoam.C` | Includes `fsiCoupling.H`; reads force residual |
| Pressure/energy coupling | `pEqn.H`, `pcEqn.H`, `EEqn.H` | Changes can affect convergence and density correction |
| Dynamic mesh correction | `correctPhi.H`, `myRhoPimpleFoam.C` | Topology changes require flux/divergence care |
| Overset variant | `overRhoPimpleDyMFoam/` | Separate executable, separate `Make/` |
| Build deps | `Make/options` | Includes `../fsiCoupling` and thermo/dynamic-mesh libraries |

## CONVENTIONS

- Use `wmake myRhoPimpleFoam` from repo root for the main executable.
- Keep FSI convergence force-residual based: `fsiResidual < fsiTolerance`.
- Keep displacement-change logging diagnostic-only unless changing the shared `fsiCoupling` contract.
- Preserve `rhoPimpleFoam.C` as a reference unless intentionally updating baseline comparison behavior.
- Update variant `Make/files` separately when touching `overRhoPimpleDyMFoam`.

## ANTI-PATTERNS

- Do not make compressible solver changes only in the main file when equation headers own the actual solve step.
- Do not remove `../fsiCoupling` from `Make/options`; the main solver depends on it.
- Do not assume incompressible solver pressure conventions apply directly to this tree.
