# Repository Agent Notes

## Workspace baseline case policy

1. Workspace root path: `/root/Workspace`.
2. Baseline case archives in Workspace:
   - `case_damfailure.zip`
   - `case_transient.zip`
3. For any new testing/analysis task:
   - default to these two baselines,
   - create a copied new case directory for experiments (do not modify baselines directly),
   - compare analysis results against both baseline cases.
4. For long-running analysis:
   - run monitoring in the background,
   - report status and estimated remaining time every 30 minutes,
   - keep foreground interaction responsive for normal communication.

## fastDynamicFvMesh change policy

When making **any code change** under:

- `fastDynamicFvMesh/fastDynamicFvMesh/`

you must also update documentation in the same change set:

- `fastDynamicFvMesh/README.md` (English)
- `fastDynamicFvMesh/Doc/快速动网格程序手册.md` (Chinese)

## Minimum documentation checklist

For each fastDynamicFvMesh code change:

1. Update runtime options/keys in docs if any dictionary behavior changed.
2. Update I/O description if file formats or outputs changed.
3. Update behavior description if coupling/restart/force/timing logic changed.
4. Keep English and Chinese docs consistent (same feature intent, same defaults).

## Verification reminder

Before finishing fastDynamicFvMesh changes, run:

```bash
wmake libso fastDynamicFvMesh
```

## Case/config synchronization

If a change modifies runtime controls, file formats, or run workflow:

1. Update relevant example case configs (for example `openfoam_transient/constant/dynamicMeshDict`).
2. Update/create matching sample input files when required (for example `mode/StructForce*.csv`).
3. Keep run script help text and docs aligned when script behavior changes
   (for example `run_myPimpleFoam.sh -h`, README, and Chinese manual).

## Backward compatibility

When adding new fastDynamicFvMesh controls:

1. Provide behavior-safe defaults.
2. Preserve existing behavior for cases that do not set new keys.

## Minimum verification evidence

For fastDynamicFvMesh code updates, record evidence from:

1. Library compile success (`wmake libso fastDynamicFvMesh`).
2. At least one short smoke run showing the changed path is exercised.
