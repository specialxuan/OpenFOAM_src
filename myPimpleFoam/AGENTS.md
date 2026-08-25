# myPimpleFoam Knowledge

## OVERVIEW

Incompressible PIMPLE solver family with FSI wiring plus overset and SRF executable variants.

## STRUCTURE

| Path | Purpose |
|---|---|
| `myPimpleFoam.C` | Custom incompressible FSI solver main |
| `pimpleLoopBody.H` | Header-included mesh update, momentum, pressure, turbulence sequence |
| `createFields.H`, `createUfIfPresent.H` | Field setup for moving-mesh solver path |
| `overPimpleDyMFoam/` | Separate incompressible overset/dynamic executable |
| `SRFPimpleFoam/` | Separate rotating-frame solver variant |
| `Make/files`, `Make/options` | Main executable target and FSI include path |

## WHERE TO LOOK

| Task | Location | Notes |
|---|---|---|
| FSI mode behavior | `myPimpleFoam.C` | Partitioned and integrated branches include `fsiCoupling.H` |
| PIMPLE sequence | `pimpleLoopBody.H` | Preserve include order around mesh motion, `UEqn.H`, `pEqn.H` |
| Overset variant | `overPimpleDyMFoam/` | Own `Make/`; not automatically covered by main solver build |
| SRF variant | `SRFPimpleFoam/` | Own field/model assumptions and `Make/` |
| Build deps | `Make/options` | Must keep `../fsiCoupling` include for main FSI solver |

## CONVENTIONS

- Use `wmake myPimpleFoam` from repo root for the main executable.
- Keep FSI convergence force-residual based: `fsiResidual < fsiTolerance`.
- Partitioned FSI freezes mesh motion inside the PIMPLE sub-loop; integrated/implicit FSI moves mesh every PIMPLE loop.
- Keep displacement-change output as diagnostics and reset the reference after topology changes.
- Update variant `Make/files` and `Make/options` independently when touching `overPimpleDyMFoam` or `SRFPimpleFoam`.

## ANTI-PATTERNS

- Do not assume variant `Make/options` inherit from the main solver.
- Do not reorder OpenFOAM include snippets casually; execution order is solver logic.
- Do not remove `../fsiCoupling` from the main solver build path.
- Do not port FSI behavior to overset/SRF variants unless their runtime loops are explicitly updated and validated.
