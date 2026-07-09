# interMixingFoam Knowledge

## OVERVIEW

Nested three-phase VOF solver and local mixture-model stack inside the `myInterFoam` family.

## STRUCTURE

| Path                                                  | Purpose                                                            |
| ----------------------------------------------------- | ------------------------------------------------------------------ |
| `interMixingFoam.C`                                   | Solver main and time-loop orchestration                            |
| `alphaControls.H`, `alphaEqn.H`, `alphaEqnSubCycle.H` | Header-included phase-fraction solve pieces                        |
| `createFields.H`                                      | Field and model construction                                       |
| `incompressibleThreePhaseMixture/`                    | Three-phase mixture implementation                                 |
| `immiscibleIncompressibleThreePhaseMixture/`          | Immiscible mixture specialization                                  |
| `threePhaseInterfaceProperties/`                      | Interface-property model; registers AMR gradient magnitude support |
| `Make/files`, `Make/options`                          | Multi-source executable build wiring                               |

## WHERE TO LOOK

| Task                | Location                         | Notes                                                                |
| ------------------- | -------------------------------- | -------------------------------------------------------------------- |
| Solver control flow | `interMixingFoam.C`              | OpenFOAM header-include style; avoid moving include snippets blindly |
| Alpha transport     | `alphaEqn*.H`                    | Coupled to local mixture/interface classes                           |
| AMR gradient reuse  | `threePhaseInterfaceProperties/` | `alphaGradMagForAMR` can feed `gradIndicatorMagnitudeField`          |
| Build source list   | `Make/files`                     | Includes local model `.C` files before solver main                   |
| Link/include stack  | `Make/options`                   | Distinct from top-level `myInterFoam` options                        |

## CONVENTIONS

- Keep this variant self-contained; do not push three-phase-only model code into top-level `myInterFoam` unless all variants need it.
- Preserve OpenFOAM solver include-snippet ordering around field creation, alpha solve, momentum, pressure, and turbulence correction.
- When adding model sources, update `Make/files` explicitly.
- Treat `alphaGradMagForAMR` as the bridge to `fastDynamicFvMesh` gradient-indicator AMR workflows.

## ANTI-PATTERNS

- Do not assume top-level `myInterFoam/Make/options` covers this nested executable.
- Do not delete local mixture/model sources as duplicate-looking code; they define this solver variant.
- Do not bypass mesh/topology safeguards copied from the main VOF solver family.
