#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""03_assemble_fdm_case.py

Assemble a complete FDM (fastDynamicFvMesh + myInterFoam) dam-break case from
the outputs of the two preceding pipeline stages:

    01_generate_mesh.py           -> constant/polyMesh (fluid OpenFOAM mesh)
    02_run_modal_virtual_elastic.py -> mode/ (FluidNodeCoor.csv + FluidNodeDisp*.csv)

This script is fully self-contained: every configuration file (constant/*
dynamicMeshDict/transportProperties/turbulenceProperties/g, the 0/ field set,
and the system/ dictionaries) is hard-coded inline.  It does NOT copy any
existing case directory.  It only:

    1. copies constant/polyMesh from --mesh-dir
    2. copies mode/ from --mode-dir
     3. writes every other file from embedded templates, or copies --settings-dir
        when supplied
    4. verifies the mode node count matches the polyMesh point count

Only the Python standard library is used.

Usage:
    python3 03_assemble_fdm_case.py \\
        --mesh-dir /path/to/fluid-openfoam \\
        --mode-dir /path/to/mode \\
        --out-dir  /path/to/fdm_case \\
         [--end-time 0.5] [--write-interval 50] [--delta-t 0.001] [--nprocs 8] \\
        [--amr] [--amr-interval N] [--amr-lower-refine N] [--amr-unrefine N] \\
         [--amr-max-refinement N] [--amr-max-cells N] [--settings-dir DIR]
"""

import argparse
import os
import re
import shutil
import sys

# --------------------------------------------------------------------------- #
#  OpenFOAM file header                                                        #
# --------------------------------------------------------------------------- #

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
CASE_DESCRIPTION_SOURCE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, "CaseDescription.md"))


def _dict_file(object_name, location, body):
    """Full text of a plain dictionary file (constant/* and system/*)."""
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


def _field_file(class_name, object_name, dimensions, internal, boundary):
    """Full text of a 0/ field file (volScalarField / volVectorField)."""
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
        + "dimensions      %s;\n\n" % dimensions
        + "internalField   %s;\n\n" % internal
        + "boundaryField\n{\n"
        + boundary
        + "}\n\n"
        + END_DIVIDER + "\n"
    )


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def _read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def _dict_value(text, key, default="N/A"):
    match = re.search(r"(?m)^\s*%s\s+([^;]+);" % re.escape(key), text)
    return match.group(1).strip() if match else default


def copy_settings(case_dir, settings_dir):
    """Copy prepared settings while preserving the assembled mesh and mode."""
    settings_dir = os.path.abspath(settings_dir)
    if not os.path.isdir(settings_dir):
        raise RuntimeError("settings directory not found: %s" % settings_dir)
    copied = []
    for name in ("constant", "system", "0", "case.foam"):
        source = os.path.join(settings_dir, name)
        if not os.path.exists(source):
            continue
        destination = os.path.join(case_dir, name)
        if os.path.isdir(source):
            for root, _dirs, files in os.walk(source):
                relative = os.path.relpath(root, source)
                excluded = ("polyMesh", "fsiRestart")
                if any(relative == item or relative.startswith(item + os.sep)
                       for item in excluded):
                    continue
                target_root = destination if relative == "." else os.path.join(destination, relative)
                os.makedirs(target_root, exist_ok=True)
                for filename in files:
                    if relative == "." and filename == "CaseDescription.md":
                        continue
                    source_file = os.path.join(root, filename)
                    target_file = os.path.join(target_root, filename)
                    shutil.copy2(source_file, target_file)
                    copied.append(os.path.relpath(target_file, case_dir))
        else:
            shutil.copy2(source, destination)
            copied.append(name)
    if not copied:
        raise RuntimeError("settings directory has no supported settings: %s" % settings_dir)
    return copied


def write_case_description(case_dir, args):
    """Copy the foundational case description and append resolved settings."""
    with open(CASE_DESCRIPTION_SOURCE, "r", encoding="utf-8") as fh:
        foundation = fh.read().replace("\r\n", "\n").rstrip()
    dynamic_mesh = _read_text(os.path.join(case_dir, "constant", "dynamicMeshDict"))
    control = _read_text(os.path.join(case_dir, "system", "controlDict"))
    fv_solution = _read_text(os.path.join(case_dir, "system", "fvSolution"))
    fv_schemes = _read_text(os.path.join(case_dir, "system", "fvSchemes"))
    decompose_path = os.path.join(case_dir, "system", "decomposeParDict")
    decompose = _read_text(decompose_path) if os.path.isfile(decompose_path) else ""
    settings_note = (
        "预写设置目录: `%s`。以下内容为复制后的实际文件。"
        % os.path.abspath(args.settings_dir)
        if args.settings_dir else
        "设置由 `03_assemble_fdm_case.py` 的内嵌模板生成。")
    appendix = """

## 本次组装的求解参数

以下内容由 `03_assemble_fdm_case.py` 根据本次组装实际写入的文件生成。
%s

### 时间与运行

| 参数 | 值 |
|---|---:|
| 求解器 | `myInterFoam` |
| `endTime` | `%s s` |
| `deltaT` | `%s s` |
| `writeControl` | `%s` |
| `writeInterval` | `%s` |
| 并行子域数 | `%s` |
| AMR | `%s` |

### FSI 与 AMR

| 参数 | 值 |
|---|---:|
| `theta` | `%s` |
| `fsiPatches` | `%s` |
| `mappingTolerance` | `%s` |
| `couplingRelaxation` | `%s` |
| `refineInterval` | `%s` |
| `lowerRefineLevel` | `%s` |
| `upperRefineLevel` | `%s` |
| `unrefineLevel` | `%s` |
| `nBufferLayers` | `%s` |
| `maxRefinement` | `%s` |
| `maxCells` | `%s` |
| `refinementMinCellVolume` | `%s` |
| `refinementMinEdgeLength` | `%s` |

### 关键 OpenFOAM 字典

#### `constant/dynamicMeshDict`

```foam
%s
```

#### `system/controlDict`

```foam
%s
```

#### `system/fvSolution`

```foam
%s
```

#### `system/fvSchemes`

```foam
%s
```
""" % (settings_note, _dict_value(control, "endTime"),
         _dict_value(control, "deltaT"), _dict_value(control, "writeControl"),
         _dict_value(control, "writeInterval"),
         _dict_value(decompose, "numberOfSubdomains", "1"),
         "enabled" if "dynamicRefineFvMeshCoeffs" in dynamic_mesh else "disabled",
         _dict_value(dynamic_mesh, "theta"), _dict_value(dynamic_mesh, "fsiPatches"),
         _dict_value(dynamic_mesh, "mappingTolerance"),
         _dict_value(dynamic_mesh, "couplingRelaxation"),
         _dict_value(dynamic_mesh, "refineInterval"),
         _dict_value(dynamic_mesh, "lowerRefineLevel"),
         _dict_value(dynamic_mesh, "upperRefineLevel"),
         _dict_value(dynamic_mesh, "unrefineLevel"),
         _dict_value(dynamic_mesh, "nBufferLayers"),
         _dict_value(dynamic_mesh, "maxRefinement"),
         _dict_value(dynamic_mesh, "maxCells"),
         _dict_value(dynamic_mesh, "refinementMinCellVolume"),
         _dict_value(dynamic_mesh, "refinementMinEdgeLength"),
         dynamic_mesh.rstrip(), control.rstrip(), fv_solution.rstrip(),
         fv_schemes.rstrip())
    case_root = os.path.dirname(os.path.normpath(case_dir))
    _write(os.path.join(case_root, "CaseDescription.md"), foundation + appendix)
    stale_case_description = os.path.join(case_dir, "CaseDescription.md")
    if os.path.isfile(stale_case_description):
        os.remove(stale_case_description)


# --------------------------------------------------------------------------- #
#  constant/ templates                                                         #
# --------------------------------------------------------------------------- #

DYNAMIC_MESH_DICT_BODY = """dynamicFvMesh fastDynamicFvMesh;

fsiCoupling integrated;

fastDynamicFvMeshCoeffs
{
    theta              1.4;
    fsiPatches         ("FSI_FLUID");
    mappingTolerance   4e-6;
    couplingRelaxation 1;
    pressureField      p;
    rhoRef             1;
    pRef               101325;
    writeFaceDiagnostics yes;
    faceDiagnosticsMode 4;
    trackTiming        true;
}
"""

# Blueprint: damfailure_coarse_symm_AMR_fixedDt_refineEveryStep_band
# The four AMR keys are injected into fastDynamicFvMeshCoeffs right before
# trackTiming, matching the blueprint layout (trackTiming then trackAmrTiming).
_AMR_FAST_COEFFS_BODY = """    meshRefinementSupport true;
    refinementInterpTolerance 1e-8;
    refinementMappingDiagnostics true;
    refinementMinCellVolume %s;
    refinementMinEdgeLength %s;
"""

_AMR_REFINE_BODY = """dynamicRefineFvMeshCoeffs
{
    refineInterval     %d;
    field              alpha.water;
    useGradIndicator   true;
    gradIndicatorField alpha.water;
    lowerRefineLevel   %s;
    upperRefineLevel   1e10;
    unrefineLevel      %s;
    nBufferLayers      %d;
    maxRefinement      %d;
    maxCells           %d;

    correctFluxes
    (
        (phi U)
        (rhoPhi none)
        (rho*phi none)
        (alphaPhi0.water none)
        (alphaPhiUn none)
        (nHatf none)
        (meshPhi none)
        (ghf none)
    );

    dumpLevel          true;
}
"""


def build_dynamic_mesh_dict(amr=False, amr_interval=1, amr_lower_refine=12,
                            amr_unrefine=4, amr_nbuffer=4, amr_max_refine=1,
                            amr_max_cells=500000, amr_min_cell_vol=5e-10,
                            amr_min_edge=0.0005):
    """Render constant/dynamicMeshDict body.  Without amr the output is
    byte-identical to the pre-AMR baseline (no meshRefinementSupport, no
    dynamicRefineFvMeshCoeffs block).  With amr, the four refinement-support
    keys are injected into fastDynamicFvMeshCoeffs and the dynamicRefine
    block is appended, mirroring the fixedDt_refineEveryStep_band blueprint."""
    if not amr:
        return DYNAMIC_MESH_DICT_BODY
    anchor = "    trackTiming        true;\n"
    amr_keys = _AMR_FAST_COEFFS_BODY % (_fmt_num(amr_min_cell_vol),
                                        _fmt_num(amr_min_edge))
    fast_body = DYNAMIC_MESH_DICT_BODY.replace(
        anchor,
        amr_keys + anchor + "    trackAmrTiming     true;\n")
    refine_block = _AMR_REFINE_BODY % (
        amr_interval, _fmt_num(amr_lower_refine), _fmt_num(amr_unrefine),
        amr_nbuffer, amr_max_refine, amr_max_cells)
    return fast_body + "\n" + refine_block

TRANSPORT_PROPERTIES_BODY = """phases (water air);

water
{
    transportModel  Newtonian;
    nu              [0 2 -1 0 0 0 0] 1.0048e-06;
    rho             [1 -3 0 0 0 0 0] 998.2;
}

air
{
    transportModel  Newtonian;
    nu              [0 2 -1 0 0 0 0] 1.4607e-05;
    rho             [1 -3 0 0 0 0 0] 1.225;
}

sigma           [1 0 -2 0 0 0 0] 0.07;
"""

TURBULENCE_PROPERTIES_BODY = """simulationType  RAS;

RAS
{
    model           kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
"""

G_BODY = """dimensions      [0 1 -2 0 0 0 0];
value           (0 -9.81 0);
"""


def write_constant(case_dir, amr=False, amr_interval=1, amr_lower_refine=12,
                   amr_unrefine=4, amr_nbuffer=4, amr_max_refine=1,
                   amr_max_cells=500000, amr_min_cell_vol=5e-10,
                   amr_min_edge=0.0005):
    body = build_dynamic_mesh_dict(
        amr=amr, amr_interval=amr_interval, amr_lower_refine=amr_lower_refine,
        amr_unrefine=amr_unrefine, amr_nbuffer=amr_nbuffer,
        amr_max_refine=amr_max_refine, amr_max_cells=amr_max_cells,
        amr_min_cell_vol=amr_min_cell_vol, amr_min_edge=amr_min_edge)
    _write(os.path.join(case_dir, "constant", "dynamicMeshDict"),
           _dict_file("dynamicMeshDict", "constant", body))
    _write(os.path.join(case_dir, "constant", "transportProperties"),
           _dict_file("transportProperties", "constant", TRANSPORT_PROPERTIES_BODY))
    _write(os.path.join(case_dir, "constant", "turbulenceProperties"),
           _dict_file("turbulenceProperties", "constant", TURBULENCE_PROPERTIES_BODY))
    _write(os.path.join(case_dir, "constant", "g"),
           _dict_file("g", "constant", G_BODY))


# --------------------------------------------------------------------------- #
#  0/ field templates                                                          #
# --------------------------------------------------------------------------- #

_SYM_PLANES = """    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
"""


def _write_u(case_dir):
    boundary = """    FSI_FLUID
    {
        type            movingWallVelocity;
        value           uniform (0 0 0);
    }
""" + _SYM_PLANES + """    "(TOP|RIGHT|LEFT|BOTTOM_FLUID)"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
"""
    _write(os.path.join(case_dir, "0", "U"),
           _field_file("volVectorField", "U", "[0 1 -1 0 0 0 0]",
                       "uniform (0 0 0)", boundary))


def _write_p(case_dir):
    boundary = _SYM_PLANES + """    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            calculated;
        value           uniform 0;
    }
"""
    _write(os.path.join(case_dir, "0", "p"),
           _field_file("volScalarField", "p", "[1 -1 -2 0 0 0 0]",
                       "uniform 0", boundary))


def _write_p_rgh(case_dir):
    boundary = _SYM_PLANES + """    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            fixedFluxPressure;
        value           uniform 0;
    }
"""
    _write(os.path.join(case_dir, "0", "p_rgh"),
           _field_file("volScalarField", "p_rgh", "[1 -1 -2 0 0 0 0]",
                       "uniform 0", boundary))


def _write_k(case_dir):
    boundary = _SYM_PLANES + """    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            kqRWallFunction;
        value           uniform 1e-10;
    }
"""
    _write(os.path.join(case_dir, "0", "k"),
           _field_file("volScalarField", "k", "[0 2 -2 0 0 0 0]",
                       "uniform 1e-10", boundary))


def _write_omega(case_dir):
    boundary = _SYM_PLANES + """    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            omegaWallFunction;
        value           uniform 1e-10;
    }
"""
    _write(os.path.join(case_dir, "0", "omega"),
           _field_file("volScalarField", "omega", "[0 0 -1 0 0 0 0]",
                       "uniform 1e-10", boundary))


def _write_nut(case_dir):
    boundary = _SYM_PLANES + """    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
"""
    _write(os.path.join(case_dir, "0", "nut"),
           _field_file("volScalarField", "nut", "[0 2 -1 0 0 0 0]",
                       "uniform 0", boundary))


def _write_alpha_water(case_dir):
    boundary = _SYM_PLANES + """    TOP
    {
        type            zeroGradient;
    }
    LEFT
    {
        type            zeroGradient;
    }
    RIGHT
    {
        type            zeroGradient;
    }
    BOTTOM_FLUID
    {
        type            zeroGradient;
    }
    FSI_FLUID
    {
        type            zeroGradient;
    }
"""
    _write(os.path.join(case_dir, "0", "alpha.water"),
           _field_file("volScalarField", "alpha.water", "[0 0 0 0 0 0 0]",
                       "uniform 0", boundary))


def write_zero_dir(case_dir):
    _write_u(case_dir)
    _write_p(case_dir)
    _write_p_rgh(case_dir)
    _write_k(case_dir)
    _write_omega(case_dir)
    _write_nut(case_dir)
    _write_alpha_water(case_dir)


# --------------------------------------------------------------------------- #
#  system/ templates                                                           #
# --------------------------------------------------------------------------- #

def write_control_dict(case_dir, end_time, delta_t, write_interval):
    body = """application     myInterFoam;

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

maxAlphaCo      1;

maxDeltaT       1;
""" % (_fmt_num(end_time), _fmt_num(delta_t), write_interval)
    _write(os.path.join(case_dir, "system", "controlDict"),
           _dict_file("controlDict", "system", body))


FV_SCHEMES_BODY = """ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         cellLimited leastSquares 1.0;
    grad(U)         cellLimited leastSquares 1.0;
}

divSchemes
{
    default         none;
    div(rho*phi,U)  Gauss linearUpwind grad(U);
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(phi,k)      Gauss linearUpwind default;
    div(phi,omega)  Gauss linearUpwind default;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
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

FV_SOLUTION_BODY = """solvers
{
    "alpha.water.*"
    {
        nAlphaCorr      2;
        nAlphaSubCycles 1;
        cAlpha          1;

        MULESCorr       yes;
        nLimiterIter    3;

        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }

    "p_rgh.*"
    {
        solver          GAMG;
        tolerance       1e-7;
        relTol          0.01;
        smoother        GaussSeidel;
    }

    "p_rghFinal"
    {
        $p_rgh;
        tolerance       1e-7;
        relTol          0;
    }

    "pcorr.*"
    {
        $p_rgh;
        tolerance       1e-2;
        relTol          0;
    }

    "pcorrFinal"
    {
        $pcorr;
        tolerance       1e-2;
        relTol          0;
    }

    "U.*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0;
    }

    "(k|omega|nut).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor   yes;
    nOuterCorrectors    50;
    nCorrectors         1;
    nNonOrthogonalCorrectors 0;
    pRefPoint       (0.392 0.36 0.006);
    pRefValue       101325;

    residualControl
    {
        p_rgh
        {
            tolerance 1e-3;
            relTol 0;
        }
        U
        {
            tolerance 1e-3;
            relTol 0;
        }
        alpha.water
        {
            tolerance 1e-3;
            relTol 0;
        }
        k
        {
            tolerance 1e-3;
            relTol 0;
        }
        omega
        {
            tolerance 1e-3;
            relTol 0;
        }
    }
}

relaxationFactors
{
    fields
    {
        p_rgh           0.3;
        alpha.water     1.0;
    }
    equations
    {
        U               0.7;
        k               0.8;
        omega           0.8;
        ".*"            1.0;
    }
}
"""

SET_FIELDS_DICT_BODY = """defaultFieldValues
(
    volScalarFieldValue alpha.water 0
);

regions
(
    boxToCell
    {
        box (0 0 0) (0.146 0.292 0.012);
        fieldValues
        (
            volScalarFieldValue alpha.water 1
        );
    }
);
"""


def write_system(case_dir, end_time, delta_t, write_interval, nprocs):
    write_control_dict(case_dir, end_time, delta_t, write_interval)
    _write(os.path.join(case_dir, "system", "fvSchemes"),
           _dict_file("fvSchemes", "system", FV_SCHEMES_BODY))
    _write(os.path.join(case_dir, "system", "fvSolution"),
           _dict_file("fvSolution", "system", FV_SOLUTION_BODY))
    _write(os.path.join(case_dir, "system", "setFieldsDict"),
           _dict_file("setFieldsDict", "system", SET_FIELDS_DICT_BODY))
    if nprocs > 1:
        write_decompose_dict(case_dir, nprocs)


def write_decompose_dict(case_dir, nprocs):
    body = """numberOfSubdomains %d;

method          scotch;

distributed     no;

roots           ();
""" % nprocs
    _write(os.path.join(case_dir, "system", "decomposeParDict"),
           _dict_file("decomposeParDict", "system", body))


# --------------------------------------------------------------------------- #
#  Small helpers                                                               #
# --------------------------------------------------------------------------- #

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


# --------------------------------------------------------------------------- #
#  Node / point counting (verification)                                        #
# --------------------------------------------------------------------------- #

def read_mode_node_count(mode_dir):
    """Return the node count declared in mode/FluidNodeCoor.csv (2nd field)."""
    path = os.path.join(mode_dir, "FluidNodeCoor.csv")
    if not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
    nums = []
    for part in first.split(","):
        part = part.strip()
        if not part:
            continue
        try:
            nums.append(float(part))
        except ValueError:
            pass
    if len(nums) >= 2:
        return int(nums[1])
    return None


def read_poly_mesh_point_count(poly_dir):
    """Return nPoints from constant/polyMesh (owner note, then points file)."""
    owner = os.path.join(poly_dir, "owner")
    if os.path.isfile(owner):
        with open(owner, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        m = re.search(r'nPoints\s*:\s*(\d+)', text)
        if m:
            return int(m.group(1))

    points = os.path.join(poly_dir, "points")
    if os.path.isfile(points):
        with open(points, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        m = re.search(r'\n\s*(\d+)\s*\n\s*\(', text)
        if m:
            return int(m.group(1))

    return None


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Assemble a complete FDM dam-break case from "
                    "constant/polyMesh (01_generate_mesh.py) and mode/ "
                    "(02_run_modal_virtual_elastic.py)."
    )
    p.add_argument("--mesh-dir", required=True,
                   help="fluid OpenFOAM dir containing constant/polyMesh")
    p.add_argument("--mode-dir", required=True,
                   help="dir containing FluidNodeCoor.csv + FluidNodeDisp*.csv")
    p.add_argument("--out-dir", required=True,
                   help="output case directory")
    p.add_argument("--settings-dir", default=None,
                   help="prepared settings root containing constant/, system/, 0/; "
                        "copied over embedded defaults without replacing polyMesh")
    p.add_argument("--end-time", type=float, default=0.5,
                   help="controlDict endTime (default 0.5)")
    p.add_argument("--write-interval", type=int, default=50,
                   help="controlDict writeInterval (default 50)")
    p.add_argument("--delta-t", type=float, default=0.001,
                   help="controlDict deltaT (default 0.001)")
    p.add_argument("--nprocs", type=int, default=8,
                   help="number of MPI subdomains for decomposeParDict "
                        "(default 8)")
    p.add_argument("--amr", action="store_true",
                   help="enable dynamic AMR (dynamicRefineFvMeshCoeffs) in "
                        "dynamicMeshDict (default off)")
    p.add_argument("--amr-interval", type=int, default=1,
                   help="AMR refineInterval (steps between mesh refinement) "
                        "(default 1)")
    p.add_argument("--amr-lower-refine", type=float, default=12,
                   help="AMR lowerRefineLevel gradient threshold (default 12)")
    p.add_argument("--amr-unrefine", type=float, default=4,
                   help="AMR unrefineLevel (default 4)")
    p.add_argument("--amr-max-refinement", type=int, default=1,
                   help="AMR maxRefinement levels (default 1)")
    p.add_argument("--amr-max-cells", type=int, default=500000,
                   help="AMR maxCells limit (default 500000)")
    p.add_argument("--amr-min-cell-volume", type=float, default=5e-10,
                   help="AMR minimum child-cell volume (default 5e-10)")
    p.add_argument("--amr-min-edge-length", type=float, default=0.0005,
                   help="AMR minimum child-edge length (default 0.0005 m)")
    return p.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    print("[FDM-PIPE] ===== 步骤 3/5: 算例组装 =====")

    mesh_dir = os.path.abspath(args.mesh_dir)
    mode_dir = os.path.abspath(args.mode_dir)
    out_dir = os.path.abspath(args.out_dir)

    # ---- validate inputs --------------------------------------------------
    src_poly = os.path.join(mesh_dir, "constant", "polyMesh")
    if not os.path.isdir(src_poly):
        print("[FDM ERR] --mesh-dir has no constant/polyMesh: %s" % mesh_dir)
        return 1
    src_mode = mode_dir
    if not os.path.isfile(os.path.join(src_mode, "FluidNodeCoor.csv")):
        print("[FDM ERR] --mode-dir has no FluidNodeCoor.csv: %s" % mode_dir)
        return 1

    # ---- copy mesh + mode --------------------------------------------------
    dst_poly = os.path.join(out_dir, "constant", "polyMesh")
    dst_mode = os.path.join(out_dir, "mode")
    os.makedirs(os.path.dirname(dst_poly), exist_ok=True)
    if os.path.exists(dst_poly):
        shutil.rmtree(dst_poly)
    shutil.copytree(src_poly, dst_poly)
    if os.path.exists(dst_mode):
        shutil.rmtree(dst_mode)
    shutil.copytree(src_mode, dst_mode)
    print("[FDM 03] Copied constant/polyMesh <- %s" % src_poly)
    print("[FDM 03] Copied mode/            <- %s" % src_mode)

    # ---- verify node/point count match ------------------------------------
    mode_nodes = read_mode_node_count(dst_mode)
    mesh_points = read_poly_mesh_point_count(dst_poly)
    if mode_nodes is None:
        print("[FDM ERR] could not parse node count from mode/FluidNodeCoor.csv")
        return 1
    if mesh_points is None:
        print("[FDM ERR] could not determine polyMesh point count")
        return 1
    print("[FDM 03] mode nodes   : %d" % mode_nodes)
    print("[FDM 03] mesh points  : %d" % mesh_points)
    if mode_nodes != mesh_points:
        print("[FDM ERR] mode node count (%d) != polyMesh point count (%d)"
              % (mode_nodes, mesh_points))
        return 1
    print("[FDM 03] OK: mode node count matches polyMesh points")

    # ---- write embedded config --------------------------------------------
    write_constant(out_dir, amr=args.amr, amr_interval=args.amr_interval,
                   amr_lower_refine=args.amr_lower_refine,
                   amr_unrefine=args.amr_unrefine,
                   amr_nbuffer=4,
                   amr_max_refine=args.amr_max_refinement,
                   amr_max_cells=args.amr_max_cells,
                   amr_min_cell_vol=args.amr_min_cell_volume,
                   amr_min_edge=args.amr_min_edge_length)
    if args.amr:
        print("[FDM 03] AMR: enabled (interval=%d, lowerRefine=%s, "
              "unrefine=%s, maxRefinement=%d, maxCells=%d)"
              % (args.amr_interval, _fmt_num(args.amr_lower_refine),
                 _fmt_num(args.amr_unrefine), args.amr_max_refinement,
                 args.amr_max_cells))
    else:
        print("[FDM 03] AMR: disabled")
    write_zero_dir(out_dir)
    write_system(out_dir, args.end_time, args.delta_t, args.write_interval,
                 args.nprocs)
    if args.settings_dir:
        try:
            copied_settings = copy_settings(out_dir, args.settings_dir)
        except (OSError, RuntimeError) as exc:
            print("[FDM ERR] failed to copy --settings-dir: %s" % exc)
            return 1
        print("[FDM 03] Copied prepared settings: %s"
              % ", ".join(copied_settings))
    required_settings = (
        os.path.join(out_dir, "constant", "dynamicMeshDict"),
        os.path.join(out_dir, "system", "controlDict"),
        os.path.join(out_dir, "system", "fvSolution"),
        os.path.join(out_dir, "system", "fvSchemes"),
    )
    missing_settings = [path for path in required_settings if not os.path.isfile(path)]
    if missing_settings:
        print("[FDM ERR] assembled settings missing: %s"
              % ", ".join(missing_settings))
        return 1
    write_case_description(out_dir, args)
    print("[FDM 03] Wrote constant/{dynamicMeshDict,transportProperties,"
          "turbulenceProperties,g}")
    print("[FDM 03] Wrote 0/{U,p,p_rgh,k,omega,nut,alpha.water}")
    print("[FDM 03] Wrote system/{controlDict,fvSchemes,fvSolution,setFieldsDict%s}"
          % (",decomposeParDict" if args.nprocs > 1 else ""))
    print("[FDM 03] Wrote %s with solver parameters"
          % os.path.join(os.path.dirname(os.path.normpath(out_dir)),
                         "CaseDescription.md"))

    # ---- ParaView OpenFOAMReader anchor ------------------------------------
    # An empty case.foam in the case root tells ParaView's OpenFOAMReader to
        # load the case directory (required by 06_export_fdm_video.py).
    case_foam = os.path.join(out_dir, "case.foam")
    if not os.path.isfile(case_foam):
        _write(case_foam, "")
    print("[FDM 03] Wrote case.foam (ParaView anchor)")

    print("")
    print("[FDM 03] 组装完成: %s" % out_dir)
    print("[FDM 03] 节点数匹配: %d=%d" % (mode_nodes, mesh_points))
    print("[FDM 03] endTime=%s deltaT=%s writeInterval=%d nprocs=%d"
          % (_fmt_num(args.end_time), _fmt_num(args.delta_t),
             args.write_interval, args.nprocs))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
