#!/usr/bin/env python3
"""Assemble a runnable AGARD 445.6 FDM (fastDynamicFvMesh + myRhoPimpleFoam)
compressible FSI case from a converted fluid polyMesh and a modal CSV set.

The converter ``case/AGARD/Mesh/inp_to_polymesh.py`` turns the CalculiX/Abaqus
far-field deck into ``constant/polyMesh``.  This script only wires that mesh
together with the modal input under ``mode/`` and the compressible FSI settings
(``constant/*``, ``0/*``, ``system/*``) and writes a ``CaseDescription.md``.

The base geometry, mode files and mesh are treated as read-only inputs.  The
only mutation applied to the copied mesh is retyping the ``FRONT`` patch from
``patch`` to ``symmetryPlane`` (the half-model symmetry plane); no node is
moved, merged or removed.

Usage
-----
    python3 assemble_fdm_case.py \\
        --mesh-dir /path/to/mesh \\
        --mode-dir /root/Workspace/AGARD/Mode/Fluid8/mode_soft \\
        --out-dir  /path/to/case \\
        --wing-type soft

Every solver parameter is a documented placeholder (see ``--help`` and the
generated ``CaseDescription.md``); the user confirms them before production.
"""

import argparse
import math
import os
import re
import shutil
import sys

import numpy as np

FOAM_HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2412                                  |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
"""

DIVIDER = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
END_DIVIDER = "// ************************************************************************* //"

# Soft-wing coupled frequencies (Hz), modes 1..10, from the CalculiX solve that
# produced the mode CSVs.  Used only as a cross-check of the mode file headers.
SOFT_FREQ_HZ = (
    9.626082, 39.249460, 49.508190, 95.629780, 125.648100,
    154.359700, 165.313100, 220.347400, 248.963100, 311.560700,
)

MOL_WEIGHT = 28.9
GAMMA = 1.4
R_UNIVERSAL = 8314.46261815324


def _dict_file(object_name, location, body):
    return (
        FOAM_HEADER
        + "FoamFile\n{\n"
        + "    version     2.0;\n"
        + "    format      ascii;\n"
        + "    class       dictionary;\n"
        + '    location    "%s";\n' % location
        + "    object      %s;\n" % object_name
        + "}\n"
        + DIVIDER + "\n\n"
        + body
        + "\n" + END_DIVIDER + "\n"
    )


def _field_file(class_name, object_name, body):
    return (
        FOAM_HEADER
        + "FoamFile\n{\n"
        + "    version     2.0;\n"
        + "    format      ascii;\n"
        + "    class       %s;\n" % class_name
        + '    location    "0";\n'
        + "    object      %s;\n" % object_name
        + "}\n"
        + DIVIDER + "\n\n"
        + body
        + "\n" + END_DIVIDER + "\n"
    )


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def _fmt_num(value):
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value == int(value):
            return str(int(value))
        return format(value, ".10g")
    return str(value)


def _dict_value(text, key, default="N/A"):
    match = re.search(r"(?m)^\s*%s\s+([^;]+);" % re.escape(key), text)
    return match.group(1).strip() if match else default


def read_poly_mesh_point_count(poly_dir):
    points = os.path.join(poly_dir, "points")
    if not os.path.isfile(points):
        return None
    with open(points, "r", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    match = re.search(r"\n\s*(\d+)\s*\n\s*\(", text)
    return int(match.group(1)) if match else None


def read_openfoam_points(path):
    rows = []
    with open(path, "r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped.startswith("("):
                continue
            values = stripped.strip("()").split()
            if len(values) < 3:
                continue
            rows.append((float(values[0]), float(values[1]), float(values[2])))
    return np.asarray(rows, dtype=np.float64)


def read_coor_csv(path):
    with open(path, "r", encoding="ascii", errors="replace") as handle:
        header = handle.readline().split(",")
    declared = int(round(float(header[1])))
    modes = int(round(float(header[2])))
    points = np.loadtxt(path, delimiter=",", skiprows=1, usecols=(0, 1, 2))
    if points.shape != (declared, 3):
        raise SystemExit(
            "[FDM] ERROR: %s holds %s coordinate rows, header declares %d"
            % (path, points.shape, declared)
        )
    return declared, modes, points


def read_mode_frequencies(mode_dir, n_modes):
    freqs = []
    for mode in range(1, n_modes + 1):
        path = os.path.join(mode_dir, "FluidNodeDisp%d.csv" % mode)
        if not os.path.isfile(path):
            raise SystemExit("[FDM] ERROR: missing mode shape file %s" % path)
        with open(path, "r", encoding="ascii", errors="replace") as handle:
            first = handle.readline()
        freqs.append(float(first.split(",")[0]))
    return freqs


def retype_front_symmetry_plane(poly_dir):
    boundary = os.path.join(poly_dir, "boundary")
    with open(boundary, "r", encoding="utf-8") as handle:
        text = handle.read()
    new_text, count = re.subn(
        r"(?ms)(\bFRONT\s*\{.*?type\s+)patch(\s*;)",
        r"\1symmetryPlane\2",
        text,
        count=1,
    )
    if count != 1:
        raise SystemExit(
            "[FDM] ERROR: could not retype FRONT -> symmetryPlane in %s" % boundary
        )
    with open(boundary, "w", encoding="utf-8") as handle:
        handle.write(new_text)


def copy_mode_files(mode_dir, dst_dir, n_modes):
    os.makedirs(dst_dir, exist_ok=True)
    names = ["FluidNodeCoor.csv", "FluidPara.csv"]
    names += ["FluidNodeDisp%d.csv" % mode for mode in range(1, n_modes + 1)]
    for name in names:
        source = os.path.join(mode_dir, name)
        if not os.path.isfile(source):
            if name == "FluidPara.csv":
                print("[FDM] note: optional %s absent; skipped" % name)
                continue
            raise SystemExit("[FDM] ERROR: required mode file missing: %s" % source)
        shutil.copy2(source, os.path.join(dst_dir, name))


def build_modal_damp(freqs, damping_ratio):
    return " ".join(
        _fmt_num(2.0 * damping_ratio * 2.0 * math.pi * freq) for freq in freqs
    )


def build_dynamic_mesh_dict(args, mapping_tolerance, modal_damp):
    common = """    theta              1.4;
    fsiPatches         ("FSI");
    pressureField      p;
    pRef               %s;
    couplingRelaxation %s;
    mappingTolerance   %s;
    modalMass          (1 1 1 1 1 1 1 1 1 1);
    modalDamp          (%s);
    writeDiagnostics   true;
    trackTiming        true;
""" % (
        _fmt_num(args.pressure),
        _fmt_num(args.coupling_relaxation),
        _fmt_num(mapping_tolerance),
        modal_damp,
    )

    if not args.amr:
        body = (
            "dynamicFvMesh fastDynamicFvMesh;\n\n"
            "fsiCoupling integrated;\n"
            "maxFsiIter %d;\n"
            "fsiTolerance %s;\n\n"
            "fastDynamicFvMeshCoeffs\n{\n" % (
                args.max_fsi_iter, _fmt_num(args.fsi_tolerance))
            + common
            + "}\n"
        )
        return body

    amr_keys = (
        "    meshRefinementSupport true;\n"
        "    refinementInterpTolerance 1e-8;\n"
        "    refinementMappingDiagnostics true;\n"
        "    refinementMinCellVolume %s;\n"
        "    refinementMinEdgeLength %s;\n"
        "    trackAmrTiming true;\n" % (
            _fmt_num(args.amr_min_cell_volume), _fmt_num(args.amr_min_edge_length))
    )
    refine = """dynamicRefineFvMeshCoeffs
{
    refineInterval     %d;
    field              p;
    useGradIndicator   true;
    gradIndicatorField p;
    lowerRefineLevel   %s;
    upperRefineLevel   1e10;
    unrefineLevel      %s;
    nBufferLayers      4;
    maxRefinement      %d;
    maxCells           %d;

    correctFluxes
    (
        (phi U)
        (meshPhi none)
    );

    dumpLevel          true;
}
""" % (
        args.amr_interval,
        _fmt_num(args.amr_lower_refine),
        _fmt_num(args.amr_unrefine),
        args.amr_max_refinement,
        args.amr_max_cells,
    )
    body = (
        "dynamicFvMesh fastDynamicFvMesh;\n\n"
        "fsiCoupling integrated;\n"
        "maxFsiIter %d;\n"
        "fsiTolerance %s;\n\n"
        "fastDynamicFvMeshCoeffs\n{\n" % (
            args.max_fsi_iter, _fmt_num(args.fsi_tolerance))
        + common
        + amr_keys
        + "}\n\n"
        + refine
    )
    return body


def write_constant(case_dir, args, mapping_tolerance, modal_damp):
    body = build_dynamic_mesh_dict(args, mapping_tolerance, modal_damp)
    _write(os.path.join(case_dir, "constant", "dynamicMeshDict"),
           _dict_file("dynamicMeshDict", "constant", body))

    thermo = """thermoType
{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}

mixture
{
    specie
    {
        molWeight   %g;
    }
    thermodynamics
    {
        Cp          1005;
        Hf          0;
    }
    transport
    {
        mu          1.82e-05;
        Pr          0.71;
    }
}
""" % MOL_WEIGHT
    _write(os.path.join(case_dir, "constant", "thermophysicalProperties"),
           _dict_file("thermophysicalProperties", "constant", thermo))

    turbulence = """simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
"""
    _write(os.path.join(case_dir, "constant", "turbulenceProperties"),
           _dict_file("turbulenceProperties", "constant", turbulence))

    gravity = """dimensions      [0 1 -2 0 0 0 0];
value           (0 -9.81 0);
"""
    _write(os.path.join(case_dir, "constant", "g"),
           _dict_file("g", "constant", gravity))


FAR = '"(BACK|IN|OUT|TOP|BOTTOM)"'


def write_zero(case_dir, args, freestream):
    u_vec = "(%s 0 0)" % _fmt_num(freestream["U"])
    p = args.pressure
    t = args.temperature
    k0 = freestream["k"]
    omega0 = freestream["omega"]

    p_body = """dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform %s;

boundaryField
{
    FSI
    {
        type            zeroGradient;
    }
    FRONT
    {
        type            symmetryPlane;
    }
    %s
    {
        type            freestreamPressure;
        freestreamValue uniform %s;
    }
}
""" % (_fmt_num(p), FAR, _fmt_num(p))
    _write(os.path.join(case_dir, "0", "p"),
           _field_file("volScalarField", "p", p_body))

    t_body = """dimensions      [0 0 0 1 0 0 0];

internalField   uniform %s;

boundaryField
{
    FSI
    {
        type            zeroGradient;
    }
    FRONT
    {
        type            symmetryPlane;
    }
    %s
    {
        type            inletOutlet;
        inletValue      uniform %s;
        value           uniform %s;
    }
}
""" % (_fmt_num(t), FAR, _fmt_num(t), _fmt_num(t))
    _write(os.path.join(case_dir, "0", "T"),
           _field_file("volScalarField", "T", t_body))

    u_body = """dimensions      [0 1 -1 0 0 0 0];

internalField   uniform %s;

boundaryField
{
    FSI
    {
        type            movingWallVelocity;
        value           uniform (0 0 0);
    }
    FRONT
    {
        type            symmetryPlane;
    }
    %s
    {
        type            freestreamVelocity;
        freestreamValue uniform %s;
    }
}
""" % (u_vec, FAR, u_vec)
    _write(os.path.join(case_dir, "0", "U"),
           _field_file("volVectorField", "U", u_body))

    k_body = """dimensions      [0 2 -2 0 0 0 0];

internalField   uniform %s;

boundaryField
{
    FSI
    {
        type            kqRWallFunction;
        value           uniform %s;
    }
    FRONT
    {
        type            symmetryPlane;
    }
    %s
    {
        type            inletOutlet;
        inletValue      uniform %s;
        value           uniform %s;
    }
}
""" % (_fmt_num(k0), _fmt_num(k0), FAR, _fmt_num(k0), _fmt_num(k0))
    _write(os.path.join(case_dir, "0", "k"),
           _field_file("volScalarField", "k", k_body))

    omega_body = """dimensions      [0 0 -1 0 0 0 0];

internalField   uniform %s;

boundaryField
{
    FSI
    {
        type            omegaWallFunction;
        value           uniform %s;
    }
    FRONT
    {
        type            symmetryPlane;
    }
    %s
    {
        type            inletOutlet;
        inletValue      uniform %s;
        value           uniform %s;
    }
}
""" % (_fmt_num(omega0), _fmt_num(omega0), FAR, _fmt_num(omega0), _fmt_num(omega0))
    _write(os.path.join(case_dir, "0", "omega"),
           _field_file("volScalarField", "omega", omega_body))

    nut_body = """dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    FSI
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    FRONT
    {
        type            symmetryPlane;
    }
    %s
    {
        type            calculated;
        value           uniform 0;
    }
}
""" % FAR
    _write(os.path.join(case_dir, "0", "nut"),
           _field_file("volScalarField", "nut", nut_body))

    alphat_body = """dimensions      [1 -1 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    FSI
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
    FRONT
    {
        type            symmetryPlane;
    }
    %s
    {
        type            calculated;
        value           uniform 0;
    }
}
""" % FAR
    _write(os.path.join(case_dir, "0", "alphat"),
           _field_file("volScalarField", "alphat", alphat_body))


def write_control_dict(case_dir, args):
    body = """application     myRhoPimpleFoam;

libs            ("libfastDynamicFvMesh.so");

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         %s;

deltaT          %s;

writeControl    timeStep;

writeInterval   %d;

purgeWrite      0;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

adjustTimeStep  no;

maxCo           1;

maxDeltaT       1;
""" % (_fmt_num(args.end_time), _fmt_num(args.delta_t), args.write_interval)
    _write(os.path.join(case_dir, "system", "controlDict"),
           _dict_file("controlDict", "system", body))


FV_SCHEMES_BODY = """ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         cellLimited Gauss linear 1;
}

divSchemes
{
    default         none;

    div(phi,U)      Gauss upwind;

    energy          Gauss upwind;
    div(phi,e)      $energy;
    div(phi,K)      $energy;
    div(phi,Ekp)    $energy;

    turbulence      Gauss upwind;
    div(phi,k)      $turbulence;
    div(phi,omega)  $turbulence;

    div(phiv,p)     Gauss upwind;
    div(phid,p)     Gauss upwind;
    div((phi|interpolate(rho)),p) Gauss upwind;
    div(meshPhi,p)  Gauss linear;

    div(((rho*nuEff)*dev2(T(grad(U)))))    Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear limited corrected 0.5;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

wallDist
{
    method          meshWave;
}
"""


def write_fv_solution(case_dir, args):
    body = """solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    "pcorr.*"
    {
        $p;
        tolerance       1e-2;
        relTol          0;
    }

    "(rho|U|k|omega|e)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }

    "(rho|U|k|omega|e)Final"
    {
        $U;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor       yes;
    transonic               yes;
    nOuterCorrectors        %d;
    nCorrectors             2;
    nNonOrthogonalCorrectors 2;
    maxCo                   1;
    maxDeltaT               1;
    rDeltaTSmoothingCoeff   0.1;
    rDeltaTDampingCoeff     0.2;
}

relaxationFactors
{
    equations
    {
        ".*"     1;
    }
}
""" % args.n_outer_correctors
    _write(os.path.join(case_dir, "system", "fvSolution"),
           _dict_file("fvSolution", "system", body))


def write_system(case_dir, args):
    write_control_dict(case_dir, args)
    _write(os.path.join(case_dir, "system", "fvSchemes"),
           _dict_file("fvSchemes", "system", FV_SCHEMES_BODY))
    write_fv_solution(case_dir, args)
    body = """numberOfSubdomains %d;

method          scotch;

distributed     no;

roots           ();
""" % args.nprocs
    _write(os.path.join(case_dir, "system", "decomposeParDict"),
           _dict_file("decomposeParDict", "system", body))


def compute_freestream(args):
    r_gas = R_UNIVERSAL / MOL_WEIGHT
    sound = math.sqrt(GAMMA * r_gas * args.temperature)
    speed = args.mach * sound
    rho = args.pressure / (r_gas * args.temperature)
    k = 1.5 * (speed * 0.01) ** 2
    omega = math.sqrt(k) / (0.09 ** 0.25 * 1.0)
    return {
        "R": r_gas, "a": sound, "U": speed, "rho": rho, "k": k, "omega": omega,
    }


def write_case_description(case_dir, args, freqs, freestream, mapping_tolerance,
                           max_coord_diff, mode_nodes, mesh_points, patch_meta):
    damp_note = (
        "Yates literature value is `0.02` (`C_i = 2*0.02*2*pi*f_i`); the default "
        "here is `%s`." % _fmt_num(args.damping_ratio)
    )
    freq_rows = "\n".join(
        "| %d | `%s` | `%s` |" % (
            mode + 1, _fmt_num(freq), _fmt_num(2.0 * args.damping_ratio * 2.0 * math.pi * freq))
        for mode, freq in enumerate(freqs)
    )
    patch_rows = "\n".join(
        "| `%s` | `%s` | %d |" % (name, ptype, int(n_faces))
        for name, ptype, n_faces in patch_meta
    )
    amr_note = "enabled" if args.amr else "disabled (mesh quality still marginal)"
    text = """# AGARD 445.6 soft-wing FDM compressible FSI case

Assembled by `case/AGARD/FDM/assemble_fdm_case.py` (OpenFOAM v2412,
`myRhoPimpleFoam` + `fastDynamicFvMesh`).

## Sources (read-only inputs)

- Fluid mesh: `Fluid7.inp` (mixed C3D8 + C3D6), converted by
  `case/AGARD/Mesh/inp_to_polymesh.py` into `constant/polyMesh`.
- Modal data: `mode/` copied byte-for-byte from the soft-wing CSV set
  (`FluidNodeCoor.csv`, `FluidNodeDisp1..%d.csv`, `FluidPara.csv`).
  `StructNode*` files are intentionally not copied (structural forcing is off).

## GATES PASSED AT ASSEMBLY

- polyMesh points (written file) = `%d`.
- mode/`FluidNodeCoor.csv` nodes = `%d` (match).
- coordinate cross-check on the WRITTEN `constant/polyMesh/points` vs
  `FluidNodeCoor.csv`: max|diff| = `%s`.
- `mappingTolerance` chosen = `%s` (tight; well below the measured min edge
  `5.946e-6 m`, so no cross-mapping).
- Mesh and mode files were copied without editing geometry or CSV content.
- The only mesh mutation is `FRONT` patch type `patch -> symmetryPlane`.

## Patch map

| patch | type | nFaces |
|---|---:|---:|
%s

`FSI` is the wing surface (`movingWallVelocity` + zero-gradient `p`/`T`);
`FRONT` is the half-model symmetry plane; `BACK/IN/OUT/TOP/BOTTOM` are the
freestream far field.

## PLACEHOLDER parameters (user must confirm)

These are **parameterised placeholders**, not verified production values.

| parameter | value | note |
|---|---:|---|
| `--mach` | `%s` | documented AGARD 445.6 transonic flutter condition |
| `--pressure` | `%s` Pa | freestream static pressure |
| `--temperature` | `%s` K | freestream static temperature |
| `--damping-ratio` | `%s` | %s |
| `--delta-t` | `%s` s | default = `1e-5` |
| `--end-time` | `%s` s | default = `0.05` = half of mode-1 period |
| `--coupling-relaxation` | `%s` | must be in (0, 1] |
| `--fsi-tolerance` | `%s` | force-residual criterion |
| `--max-fsi-iter` | `%d` | explicit/partitioned path |
| `--n-outer-correctors` | `%d` | PIMPLE outer correctors |
| `--write-interval` | `%d` steps | |
| `--nprocs` | `%d` | decomposeParDict |
| `--amr` | `%s` | |

Freestream derived from the placeholders: `a = %s m/s`, `U = %s m/s`,
`rho = %s kg/m3`, `k = %s m2/s2`, `omega = %s 1/s`.
The `Re/ft = 364,600` target corresponds to the documented AGARD 445.6
transonic condition.

## Modal damping

`modalMass = 1.0` for all modes (CalculiX generalized-mass identity).
`modalDamp` is the Wilson-Theta coefficient `C`; with damping ratio
`zeta = %s`, `C_i = 2*zeta*omega_i`, `omega_i = 2*pi*f_i`.

| mode | f_i (Hz) | C_i |
|---:|---:|---:|
%s

## Run

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
checkMesh -case <this case>
myRhoPimpleFoam -case <this case>        # serial smoke
# production: decomposePar -force && mpirun -np %d myRhoPimpleFoam -parallel
```
""" % (
        len(freqs),
        mesh_points,
        mode_nodes,
        "%.3e" % max_coord_diff,
        _fmt_num(mapping_tolerance),
        patch_rows,
        _fmt_num(args.mach),
        _fmt_num(args.pressure),
        _fmt_num(args.temperature),
        _fmt_num(args.damping_ratio),
        damp_note,
        _fmt_num(args.delta_t),
        _fmt_num(args.end_time),
        _fmt_num(args.coupling_relaxation),
        _fmt_num(args.fsi_tolerance),
        args.max_fsi_iter,
        args.n_outer_correctors,
        args.write_interval,
        args.nprocs,
        amr_note,
        _fmt_num(freestream["a"]),
        _fmt_num(freestream["U"]),
        _fmt_num(freestream["rho"]),
        _fmt_num(freestream["k"]),
        _fmt_num(freestream["omega"]),
        _fmt_num(args.damping_ratio),
        freq_rows,
        args.nprocs,
    )
    _write(os.path.join(case_dir, "CaseDescription.md"), text)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Assemble an AGARD 445.6 FDM compressible FSI case from a "
                    "converted polyMesh and a modal CSV set."
    )
    parser.add_argument("--mesh-dir", required=True,
                        help="dir containing constant/polyMesh (converter output)")
    parser.add_argument("--mode-dir", required=True,
                        help="dir containing FluidNodeCoor.csv + FluidNodeDisp*.csv")
    parser.add_argument("--out-dir", required=True, help="output case directory")
    parser.add_argument("--wing-type", choices=("soft", "hard"), default="soft",
                        help="selects the frequency cross-check (default soft)")
    parser.add_argument("--mach", type=float, default=0.96,
                        help="freestream Mach number (default 0.96)")
    parser.add_argument("--pressure", type=float, default=4700.0,
                        help="freestream static pressure Pa (default 4700)")
    parser.add_argument("--temperature", type=float, default=258.01,
                        help="freestream static temperature K (default 258.01)")
    parser.add_argument("--damping-ratio", type=float, default=0.0,
                        help="modal damping ratio zeta; 0.0 default, 0.02 Yates")
    parser.add_argument("--delta-t", type=float, default=1e-5,
                        help="controlDict deltaT seconds (default 1e-5)")
    parser.add_argument("--end-time", type=float, default=0.05,
                        help="controlDict endTime seconds (default 0.05)")
    parser.add_argument("--write-interval", type=int, default=100,
                        help="controlDict writeInterval steps (default 100)")
    parser.add_argument("--nprocs", type=int, default=64,
                        help="decomposeParDict numberOfSubdomains (default 64)")
    parser.add_argument("--coupling-relaxation", type=float, default=1.0,
                        help="fastDynamicFvMesh couplingRelaxation in (0,1] "
                             "(default 1.0)")
    parser.add_argument("--fsi-tolerance", type=float, default=1e-4,
                        help="fsiTolerance force-residual threshold (default 1e-4)")
    parser.add_argument("--max-fsi-iter", type=int, default=10,
                        help="maxFsiIter for partitioned coupling (default 10)")
    parser.add_argument("--n-outer-correctors", type=int, default=2,
                        help="PIMPLE nOuterCorrectors (default 2)")
    parser.add_argument("--mapping-tolerance", type=float, default=None,
                        help="override mappingTolerance; default derives a tight "
                             "value from the measured coordinate cross-check")
    parser.add_argument("--coord-gate-tolerance", type=float, default=1e-9,
                        help="max allowed coordinate diff (default 1e-9)")
    parser.add_argument("--amr", action="store_true",
                        help="enable dynamicRefineFvMesh AMR (default off)")
    parser.add_argument("--amr-interval", type=int, default=1)
    parser.add_argument("--amr-lower-refine", type=float, default=1.0)
    parser.add_argument("--amr-unrefine", type=float, default=0.25)
    parser.add_argument("--amr-max-refinement", type=int, default=1)
    parser.add_argument("--amr-max-cells", type=int, default=600000)
    parser.add_argument("--amr-min-cell-volume", type=float, default=0.0)
    parser.add_argument("--amr-min-edge-length", type=float, default=1e-5)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    print("[FDM] ===== assembling AGARD %s-wing FDM case =====" % args.wing_type)

    mesh_dir = os.path.abspath(args.mesh_dir)
    mode_dir = os.path.abspath(args.mode_dir)
    out_dir = os.path.abspath(args.out_dir)

    src_poly = os.path.join(mesh_dir, "constant", "polyMesh")
    if not os.path.isdir(src_poly):
        raise SystemExit("[FDM] ERROR: --mesh-dir has no constant/polyMesh: %s" % mesh_dir)
    if not os.path.isfile(os.path.join(mode_dir, "FluidNodeCoor.csv")):
        raise SystemExit("[FDM] ERROR: --mode-dir has no FluidNodeCoor.csv: %s" % mode_dir)

    declared_nodes, n_modes, coor = read_coor_csv(os.path.join(mode_dir, "FluidNodeCoor.csv"))
    print("[FDM] mode CSV declares %d nodes, %d modes" % (declared_nodes, n_modes))

    freqs = read_mode_frequencies(mode_dir, n_modes)
    if args.wing_type == "soft":
        expected = np.asarray(SOFT_FREQ_HZ[:n_modes])
        diff = np.max(np.abs(np.asarray(freqs) - expected))
        print("[FDM] soft mode frequencies vs reference: max|diff| = %.3e Hz" % diff)
        if diff > 1e-3:
            raise SystemExit("[FDM] ERROR: mode frequencies do not match the soft set")

    dst_poly = os.path.join(out_dir, "constant", "polyMesh")
    os.makedirs(os.path.dirname(dst_poly), exist_ok=True)
    if os.path.exists(dst_poly):
        shutil.rmtree(dst_poly)
    shutil.copytree(src_poly, dst_poly)
    retype_front_symmetry_plane(dst_poly)
    print("[FDM] copied polyMesh <- %s (FRONT -> symmetryPlane)" % src_poly)

    dst_mode = os.path.join(out_dir, "mode")
    if os.path.exists(dst_mode):
        shutil.rmtree(dst_mode)
    copy_mode_files(mode_dir, dst_mode, n_modes)
    print("[FDM] copied mode/ <- %s (%d mode shapes + coordinates)" % (mode_dir, n_modes))

    mesh_points = read_poly_mesh_point_count(dst_poly)
    if mesh_points is None:
        raise SystemExit("[FDM] ERROR: could not read point count from written polyMesh")
    if mesh_points != declared_nodes:
        raise SystemExit(
            "[FDM] ERROR: mode node count (%d) != polyMesh point count (%d)"
            % (declared_nodes, mesh_points)
        )
    print("[FDM] GATE mode nodes == mesh points == %d" % mesh_points)

    written_points = read_openfoam_points(os.path.join(dst_poly, "points"))
    if written_points.shape != coor.shape:
        raise SystemExit(
            "[FDM] ERROR: written points shape %s != CSV shape %s"
            % (written_points.shape, coor.shape)
        )
    max_coord_diff = float(np.max(np.abs(written_points - coor)))
    if max_coord_diff > args.coord_gate_tolerance:
        raise SystemExit(
            "[FDM] ERROR: coordinate cross-check failed: max|diff| = %.6e > %s"
            % (max_coord_diff, _fmt_num(args.coord_gate_tolerance))
        )
    if args.mapping_tolerance is not None:
        mapping_tolerance = float(args.mapping_tolerance)
    else:
        mapping_tolerance = max(1e-9, 10.0 * max_coord_diff)
    print("[FDM] GATE written points vs FluidNodeCoor.csv max|diff| = %.3e" % max_coord_diff)
    print("[FDM] mappingTolerance = %s" % _fmt_num(mapping_tolerance))

    modal_damp = build_modal_damp(freqs, args.damping_ratio)
    freestream = compute_freestream(args)
    write_constant(out_dir, args, mapping_tolerance, modal_damp)
    write_zero(out_dir, args, freestream)
    write_system(out_dir, args)

    boundary_text = open(os.path.join(dst_poly, "boundary"), encoding="utf-8").read()
    patch_meta = re.findall(
        r"(\w+)\s*\{\s*type\s+(\w+);\s*nFaces\s+(\d+);", boundary_text
    )

    write_case_description(case_dir=out_dir, args=args, freqs=freqs,
                           freestream=freestream,
                           mapping_tolerance=mapping_tolerance,
                           max_coord_diff=max_coord_diff,
                           mode_nodes=declared_nodes, mesh_points=mesh_points,
                           patch_meta=patch_meta)

    print("[FDM] wrote constant/{dynamicMeshDict,thermophysicalProperties,"
          "turbulenceProperties,g}")
    print("[FDM] wrote 0/{p,T,U,k,omega,nut,alphat}")
    print("[FDM] wrote system/{controlDict,fvSchemes,fvSolution,decomposeParDict}")
    print("[FDM] wrote CaseDescription.md")
    print("[FDM] assembled case: %s" % out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
