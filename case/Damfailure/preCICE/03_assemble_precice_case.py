#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""03_assemble_precice_case.py

Assemble a complete, runnable preCICE partitioned dam-break FSI case from the
outputs of 01_generate_mesh.py (preCICE variant):

    01_generate_mesh.py  --out-dir <mesh_dir>
        +- fluid-openfoam/constant/polyMesh      (OpenFOAM fluid mesh)
        +- fluid-openfoam/case.foam              (--view)
        `- solid-calculix/dam_gate.inp           (CalculiX C3D8 gate)

This script produces:

    <out-dir>/
    |-- precice-config.xml          (verified template, parallel-implicit + IQN-ILS)
    |-- run-coupled.sh              (two-background-subshell coupled launcher)
    |-- fluid-openfoam/
    |   |-- constant/
    |   |   |-- polyMesh/           (copied from 01)
    |   |   |-- dynamicMeshDict     (dynamicMotionSolverFvMesh + displacementLaplacian)
    |   |   |-- g / transportProperties / turbulenceProperties
    |   |-- system/
    |   |   |-- controlDict         (interFoam, deltaT = precice time-window)
    |   |   |-- preciceDict         (Fluid participant, FSI module, 2 interfaces)
    |   |   |-- fvSchemes / fvSolution / setFieldsDict / decomposeParDict
    |   |-- 0/                      (U alpha.water p_rgh p k omega nut pointDisplacement)
    |   `-- case.foam
    `-- solid-calculix/
        |-- dam_gate.inp            (copied, *MATERIAL overridden to FDM defaults)
        |-- config.yml              (Solid participant adapter config)
        `-- run.sh                  (ccx_preCICE debug launcher, OMP=1)

Every dictionary / field template is taken verbatim from the verified baseline
zip precice_damfailure_fine_noAMR.zip (all values are the zip's exact values).

The solid material follows the FDM baseline (E=1e6, nu=0, rho=2500), NOT the
zip template's 1.4e6/0.4/10000.  Overridable via
--gate-e/--gate-nu/--gate-rho.  Only the *MATERIAL values in dam_gate.inp are
patched; the geometry / *NSET / *BOUNDARY structure is left untouched.

Usage:
    python3 03_assemble_precice_case.py \\
        --mesh-dir  /path/to/fluid-openfoam \\
        --solid-dir /path/to/solid-calculix \\
        --out-dir   /path/to/precice_case \\
        [--end-time 0.02] [--write-interval 0.01] [--delta-t 0.0005] \\
        [--nprocs 1] [--gate-e 1e6] [--gate-nu 0.0] [--gate-rho 2500]
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


def log(msg=""):
    print(("[PRECICE 03] " + msg) if msg else "[PRECICE 03]", flush=True)


def warn(msg):
    print("[PRECICE WARN] " + msg, flush=True)


def error(msg):
    print("[PRECICE ERR] " + msg, file=sys.stderr, flush=True)


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


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


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


def _field_file(class_name, object_name, dimensions, internal, boundary):
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


# --------------------------------------------------------------------------- #
#  Validation helpers                                                          #
# --------------------------------------------------------------------------- #

def _read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def check_fluid_mesh(mesh_dir):
    """Verify the 01 fluid mesh exists and carries the FSI_FLUID patch."""
    poly = os.path.join(mesh_dir, "constant", "polyMesh")
    problems = []
    for name in ("points", "faces", "owner", "neighbour", "boundary"):
        if not os.path.isfile(os.path.join(poly, name)):
            problems.append("fluid polyMesh missing %s (expected in %s)"
                            % (name, poly))
    boundary = os.path.join(poly, "boundary")
    if os.path.isfile(boundary):
        text = _read_text(boundary)
        if not re.search(r'^\s*FSI_FLUID\s*$', text, flags=re.M):
            problems.append("fluid polyMesh boundary has no 'FSI_FLUID' patch")
    return problems


def check_solid_inp(solid_path):
    """Verify the CalculiX deck exists with the GATE material + Ninterface nset."""
    if not os.path.isfile(solid_path):
        return ["solid deck not found: %s" % solid_path]
    text = _read_text(solid_path)
    problems = []
    if not re.search(r'^\*MATERIAL,\s*NAME\s*=\s*GATE\s*$', text, flags=re.M | re.I):
        problems.append("dam_gate.inp missing '*MATERIAL, NAME=GATE'")
    if not re.search(r'^\*ELASTIC\s*$', text, flags=re.M | re.I):
        problems.append("dam_gate.inp missing '*ELASTIC' card")
    if not re.search(r'^\*DENSITY\s*$', text, flags=re.M | re.I):
        problems.append("dam_gate.inp missing '*DENSITY' card")
    if not re.search(r'^\*NSET,\s*NSET\s*=\s*Ninterface\s*$', text, flags=re.M | re.I):
        problems.append("dam_gate.inp missing '*NSET, NSET=Ninterface' "
                        "(preCICE solid adapter patch 'interface' -> 'Ninterface')")
    if not re.search(r'^\*DYNAMIC,\s*ALPHA\s*=\s*0\.0\s*,\s*DIRECT\s*$',
                     text, flags=re.M | re.I):
        problems.append("dam_gate.inp missing '*DYNAMIC, ALPHA=0.0, DIRECT' "
                        "(solid must march at the coupling time-window)")
    return problems


def copy_tree(src, dst):
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


# --------------------------------------------------------------------------- #
#  Material patching for the solid deck                                         #
# --------------------------------------------------------------------------- #

def patch_solid_material(src, dst, e_mod, nu, rho):
    """Copy dam_gate.inp, overriding only the *MATERIAL/*ELASTIC/*DENSITY values.

    The preCICE solid material matches the FDM baseline (E=1e6, nu=0, rho=2500)
    so results are comparable.  Geometry, *NSET and *BOUNDARY cards are copied
    byte-for-byte.
    """
    lines = _read_text(src).splitlines()
    out = []
    pending = None  # "ELASTIC" | "DENSITY"
    patched = {"ELASTIC": False, "DENSITY": False}
    for raw in lines:
        s = raw.strip()
        if pending:
            if not s:
                out.append(raw)
                continue
            if pending == "ELASTIC":
                out.append("%s, %s" % (_fmt_num(e_mod), _fmt_num(nu)))
                patched["ELASTIC"] = True
            else:
                out.append(_fmt_num(rho))
                patched["DENSITY"] = True
            pending = None
            continue
        if s.upper().startswith("*MATERIAL"):
            pending = None
        elif s.upper() == "*ELASTIC":
            pending = "ELASTIC"
        elif s.upper() == "*DENSITY":
            pending = "DENSITY"
        out.append(raw)
    if not (patched["ELASTIC"] and patched["DENSITY"]):
        raise RuntimeError(
            "could not locate *ELASTIC/*DENSITY data lines in %s "
            "(patched=%s); refusing to write an invalid deck" % (src, patched))
    _write(dst, "\n".join(out) + "\n")
    return (e_mod, nu, rho)


# --------------------------------------------------------------------------- #
#  precice-config.xml  (verified zip baseline)                                  #
# --------------------------------------------------------------------------- #

PRECICE_CONFIG_XML = """<?xml version="1.0" encoding="UTF-8" ?>
<precice-configuration>
  <log>
    <sink
      filter="%Severity% &gt; debug and %Rank% = 0"
      format="---[precice] %ColorizedSeverity% %Message%"
      enabled="true" />
  </log>

  <data:vector name="Force" />
  <data:vector name="DisplacementDelta" />

  <mesh name="Fluid-Mesh-Nodes" dimensions="3">
    <use-data name="DisplacementDelta" />
  </mesh>

  <mesh name="Fluid-Mesh-Faces" dimensions="3">
    <use-data name="Force" />
  </mesh>

  <mesh name="Solid-Mesh" dimensions="3">
    <use-data name="Force" />
    <use-data name="DisplacementDelta" />
  </mesh>

  <participant name="Fluid">
    <receive-mesh name="Solid-Mesh" from="Solid" />
    <provide-mesh name="Fluid-Mesh-Nodes" />
    <provide-mesh name="Fluid-Mesh-Faces" />
    <write-data name="Force" mesh="Fluid-Mesh-Faces" />
    <read-data name="DisplacementDelta" mesh="Fluid-Mesh-Nodes" />
    <mapping:rbf-global-direct
      direction="write"
      from="Fluid-Mesh-Faces"
      to="Solid-Mesh"
      constraint="conservative">
      <basis-function:thin-plate-splines />
    </mapping:rbf-global-direct>
    <mapping:rbf-global-direct
      direction="read"
      from="Solid-Mesh"
      to="Fluid-Mesh-Nodes"
      constraint="consistent">
      <basis-function:thin-plate-splines />
    </mapping:rbf-global-direct>
  </participant>

  <participant name="Solid">
    <provide-mesh name="Solid-Mesh" />
    <read-data name="Force" mesh="Solid-Mesh" />
    <write-data name="DisplacementDelta" mesh="Solid-Mesh" />
    <watch-point mesh="Solid-Mesh" name="Gate-Top-Center" coordinate="0.2989999925;0.07999999821;0.0060000001" />
  </participant>

  <m2n:sockets acceptor="Fluid" connector="Solid" exchange-directory=".." />

<coupling-scheme:parallel-implicit>
  <participants first="Fluid" second="Solid" />
  <time-window-size value="0.0005" />
  <max-time value="0.5" />
  <exchange data="Force" mesh="Solid-Mesh" from="Fluid" to="Solid" />
  <exchange data="DisplacementDelta" mesh="Solid-Mesh" from="Solid" to="Fluid" />
  <max-iterations value="50" />
  <min-iterations value="3" />
  <absolute-convergence-measure limit="5e-7" data="DisplacementDelta" mesh="Solid-Mesh" />
  <absolute-convergence-measure limit="1e-3" data="Force" mesh="Solid-Mesh" />
  <acceleration:IQN-ILS>
    <data name="DisplacementDelta" mesh="Solid-Mesh" />
    <data name="Force" mesh="Solid-Mesh" />
    <preconditioner type="residual-sum" />
    <filter type="QR2" limit="1e-1" />
    <initial-relaxation value="0.1" />
    <max-used-iterations value="50" />
    <time-windows-reused value="5" />
  </acceleration:IQN-ILS>
</coupling-scheme:parallel-implicit>
</precice-configuration>
"""


# --------------------------------------------------------------------------- #
#  run-coupled.sh  (verified zip baseline: two background subshells + wait)     #
# --------------------------------------------------------------------------- #

RUN_COUPLED_SH = """#!/bin/bash
# Coupled preCICE run: solid (ccx_preCICE) + fluid (interFoam) started
# simultaneously as background subshells, logs into precice-run/.
# Generated by 03_assemble_precice_case.py from the verified baseline zip.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p precice-run
rm -f precice-run/log.fluid precice-run/log.solid
(
  cd solid-calculix
  export OMP_NUM_THREADS=4
  export CCX_NPROC_EQUATION_SOLVER=4
  ccx_preCICE -i dam_gate -precice-participant Solid > ../precice-run/log.solid 2>&1
) &
solid_pid=$!
(
  cd fluid-openfoam
  interFoam > ../precice-run/log.fluid 2>&1
) &
fluid_pid=$!
wait "$fluid_pid"
fluid_status=$?
wait "$solid_pid"
solid_status=$?
echo "fluid exit status: $fluid_status"
echo "solid exit status: $solid_status"
if [ "$fluid_status" -ne 0 ] || [ "$solid_status" -ne 0 ]; then
  exit 1
fi
"""


# --------------------------------------------------------------------------- #
#  fluid-openfoam/constant/  (zip baseline values)                             #
# --------------------------------------------------------------------------- #

DYNAMIC_MESH_DICT = _dict_file("dynamicMeshDict", "constant", """dynamicFvMesh       dynamicMotionSolverFvMesh;

motionSolverLibs    (fvMotionSolvers);

motionSolver        displacementLaplacian;

diffusivity         quadratic inverseDistance 1(FSI_FLUID);
""")

TRANSPORT_PROPERTIES = _dict_file("transportProperties", "constant", """phases (water air);

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
""")

TURBULENCE_PROPERTIES = _dict_file("turbulenceProperties", "constant", """simulationType  RAS;

RAS
{
    model           kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
""")

G_CONSTANT = _dict_file("g", "constant", """dimensions      [0 1 -2 0 0 0 0];
value           (0 -9.81 0);
""")


# --------------------------------------------------------------------------- #
#  fluid-openfoam/system/  (zip baseline values)                               #
# --------------------------------------------------------------------------- #

def control_dict(end_time, delta_t, write_interval):
    return _dict_file("controlDict", "system", """application     interFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         %s;

deltaT          %s;

writeControl    runTime;

writeInterval   %s;

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

functions
{
    preCICE_Adapter
    {
        type            preciceAdapterFunctionObject;
        libs            ( "libpreciceAdapterFunctionObject.so" );
    }
}
""" % (_fmt_num(end_time), _fmt_num(delta_t), _fmt_num(write_interval)))


PRECICE_DICT = _dict_file("preciceDict", "system", """preciceConfig "../precice-config.xml";

participant Fluid;

modules (FSI);

interfaces
{
    Interface1
    {
        mesh              Fluid-Mesh-Nodes;
        patches           (FSI_FLUID);
        locations         faceNodes;

        readData
        (
            DisplacementDelta
        );

        writeData
        (
        );
    }

    Interface2
    {
        mesh              Fluid-Mesh-Faces;
        patches           (FSI_FLUID);
        locations         faceCenters;

        readData
        (
        );

        writeData
        (
            Force
        );
    }
}

FSI
{
    rho rho [1 -3 0 0 0 0 0] 998.2;
}
""")

FV_SCHEMES = _dict_file("fvSchemes", "system", """ddtSchemes
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
    div(phi,k)      Gauss upwind;
    div(phi,omega)  Gauss upwind;
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
""")

FV_SOLUTION = _dict_file("fvSolution", "system", """solvers
{
    "alpha.water.*"
    {
        nAlphaCorr      2;
        nAlphaSubCycles 4;
        cAlpha          1;
        MULESCorr       yes;
        nLimiterIter    8;
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-08;
        relTol          0;
    }
    "p_rgh.*"
    {
        solver          GAMG;
        tolerance       1e-07;
        relTol          0.01;
        smoother        GaussSeidel;
    }
    "p_rghFinal"
    {
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.001;
        smoother        GaussSeidel;
    }
    "pcorr.*"
    {
        solver          GAMG;
        tolerance       0.01;
        relTol          0;
        smoother        GaussSeidel;
    }
    "pcorrFinal"
    {
        solver          GAMG;
        tolerance       0.01;
        relTol          0;
        smoother        GaussSeidel;
    }
    "(U|cellDisplacement)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0.0001;
        minIter         3;
    }
    "(U|cellDisplacement)Final"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0;
        minIter         3;
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
    momentumPredictor yes;
    nOuterCorrectors 40;
    nCorrectors     1;
    nNonOrthogonalCorrectors 0;
    correctPhi      no;
    pRefPoint       ( 0.392 0.36 0.006 );
    pRefValue       0;
    residualControl
    {
        p_rgh
        {
            tolerance       0.01;
            relTol          0;
        }
        alpha.water
        {
            tolerance       0.001;
            relTol          0;
        }
    }
}

relaxationFactors
{
    fields
    {
        p_rgh           0.3;
        alpha.water     1;
    }
    equations
    {
        U               0.7;
        k               0.8;
        omega           0.8;
        ".*"            1;
    }
}
""")

SET_FIELDS_DICT = _dict_file("setFieldsDict", "system", """defaultFieldValues
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
""")


def decompose_par_dict(nprocs):
    return _dict_file("decomposeParDict", "system", """numberOfSubdomains %d;

method          scotch;

distributed     no;

roots           ();
""" % nprocs)


# --------------------------------------------------------------------------- #
#  fluid-openfoam/0/  (zip baseline fields)                                    #
# --------------------------------------------------------------------------- #

ALPHA_WATER_0 = _field_file(
    "volScalarField", "alpha.water", "[0 0 0 0 0 0 0]",
    "uniform 0",
    """    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    TOP
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
""")

U_0 = _field_file(
    "volVectorField", "U", "[0 1 -1 0 0 0 0]",
    "uniform (0 0 0)",
    """    FSI_FLUID
    {
        type            movingWallVelocity;
        value           uniform (0 0 0);
    }
    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    "(TOP|RIGHT|LEFT|BOTTOM_FLUID)"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
""")

P_RGH_0 = _field_file(
    "volScalarField", "p_rgh", "[1 -1 -2 0 0 0 0]",
    "uniform 0",
    """    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            fixedFluxPressure;
        value           uniform 0;
    }
""")

P_0 = _field_file(
    "volScalarField", "p", "[1 -1 -2 0 0 0 0]",
    "uniform 0",
    """    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            calculated;
        value           uniform 0;
    }
""")

K_0 = _field_file(
    "volScalarField", "k", "[0 2 -2 0 0 0 0]",
    "uniform 1e-10",
    """    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            kqRWallFunction;
        value           uniform 1e-10;
    }
""")

OMEGA_0 = _field_file(
    "volScalarField", "omega", "[0 0 -1 0 0 0 0]",
    "uniform 1e-10",
    """    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            omegaWallFunction;
        value           uniform 1e-10;
    }
""")

NUT_0 = _field_file(
    "volScalarField", "nut", "[0 2 -1 0 0 0 0]",
    "uniform 0",
    """    FRONT
    {
        type            symmetryPlane;
    }
    BACK
    {
        type            symmetryPlane;
    }
    "(TOP|RIGHT|LEFT|BOTTOM_FLUID|FSI_FLUID)"
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
""")

POINT_DISPLACEMENT_0 = _field_file(
    "pointVectorField", "pointDisplacement", "[0 1 0 0 0 0 0]",
    "uniform (0 0 0)",
    """    FRONT
    {
        type            symmetryPlane;
    }

    BACK
    {
        type            symmetryPlane;
    }

    FSI_FLUID
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }

    "(TOP|RIGHT|LEFT|BOTTOM_FLUID)"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
""")


# --------------------------------------------------------------------------- #
#  solid-calculix/  (zip baseline values)                                      #
# --------------------------------------------------------------------------- #

CONFIG_YML = """participants:

    Solid:
        interfaces:
        - nodes-mesh: Solid-Mesh
          patch: interface
          read-data: [Force]
          write-data: [DisplacementDelta]

precice-config-file: ../precice-config.xml
"""

RUN_SH = """#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
ccx_preCICE -i dam_gate -precice-participant Solid
"""


# --------------------------------------------------------------------------- #
#  Argument parsing + main                                                     #
# --------------------------------------------------------------------------- #

def parse_args(argv):
    p = argparse.ArgumentParser(
        prog="03_assemble_precice_case.py",
        description="Assemble a complete preCICE partitioned dam-break FSI case "
                    "from the 01_generate_mesh.py outputs.",
    )
    p.add_argument("--mesh-dir", required=True,
                   help="fluid-openfoam directory produced by 01_generate_mesh.py "
                        "(contains constant/polyMesh)")
    p.add_argument("--solid-dir", required=True,
                   help="solid-calculix directory produced by 01_generate_mesh.py "
                        "(contains dam_gate.inp)")
    p.add_argument("--out-dir", required=True,
                   help="assembled preCICE case root directory")
    p.add_argument("--end-time", type=float, default=0.02,
                   help="fluid endTime in system/controlDict (default 0.02 = short)")
    p.add_argument("--delta-t", type=float, default=0.0005,
                   help="fluid deltaT; must equal the preCICE time-window-size "
                        "(default 0.0005)")
    p.add_argument("--write-interval", type=float, default=0.01,
                   help="fluid writeInterval (runTime), default 0.01")
    p.add_argument("--nprocs", type=int, default=8,
                   help="numberOfSubdomains in decomposeParDict (default 8)")
    p.add_argument("--gate-e", type=float, default=1.0e6,
                   help="gate Young's modulus Pa (FDM baseline 1e6)")
    p.add_argument("--gate-nu", type=float, default=0.0,
                   help="gate Poisson ratio (FDM baseline 0)")
    p.add_argument("--gate-rho", type=float, default=2500.0,
                   help="gate density kg/m3 (FDM baseline 2500)")
    p.add_argument("--precice-max-time", type=float, default=0.5,
                   help="max-time in precice-config.xml (default 0.5)")
    return p.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    mesh_dir = os.path.abspath(args.mesh_dir)
    solid_dir = os.path.abspath(args.solid_dir)
    out_dir = os.path.abspath(args.out_dir)

    if args.delta_t <= 0:
        error("--delta-t must be positive")
        return 2
    if args.nprocs < 1:
        error("--nprocs must be a positive integer")
        return 2
    if args.end_time <= 0:
        error("--end-time must be positive")
        return 2

    log("Assembling preCICE case -> %s" % out_dir)

    # ---- validate inputs --------------------------------------------------
    problems = check_fluid_mesh(mesh_dir) + check_solid_inp(
        os.path.join(solid_dir, "dam_gate.inp"))
    if problems:
        for pr in problems:
            error(pr)
        log("assembly aborted: input validation failed")
        return 1
    log("OK: fluid polyMesh (FSI_FLUID patch) + solid dam_gate.inp (Ninterface)")

    fluid_out = os.path.join(out_dir, "fluid-openfoam")
    solid_out = os.path.join(out_dir, "solid-calculix")

    # ---- copy fluid mesh ---------------------------------------------------
    log("copying constant/polyMesh")
    copy_tree(os.path.join(mesh_dir, "constant", "polyMesh"),
              os.path.join(fluid_out, "constant", "polyMesh"))

    # ---- copy + patch solid deck ------------------------------------------
    log("copying dam_gate.inp (material -> E=%s nu=%s rho=%s)"
        % (_fmt_num(args.gate_e), _fmt_num(args.gate_nu), _fmt_num(args.gate_rho)))
    patch_solid_material(os.path.join(solid_dir, "dam_gate.inp"),
                         os.path.join(solid_out, "dam_gate.inp"),
                         args.gate_e, args.gate_nu, args.gate_rho)

    # ---- root files ---------------------------------------------------------
    _write(os.path.join(out_dir, "precice-config.xml"),
           PRECICE_CONFIG_XML.replace(
               '<max-time value="0.5" />',
               '<max-time value="%s" />' % _fmt_num(args.precice_max_time)))
    _write(os.path.join(out_dir, "run-coupled.sh"), RUN_COUPLED_SH)
    log("wrote precice-config.xml + run-coupled.sh")

    # ---- fluid constant/ ----------------------------------------------------
    _write(os.path.join(fluid_out, "constant", "dynamicMeshDict"),
           DYNAMIC_MESH_DICT)
    _write(os.path.join(fluid_out, "constant", "transportProperties"),
           TRANSPORT_PROPERTIES)
    _write(os.path.join(fluid_out, "constant", "turbulenceProperties"),
           TURBULENCE_PROPERTIES)
    _write(os.path.join(fluid_out, "constant", "g"), G_CONSTANT)
    log("wrote constant/dynamicMeshDict, g, transportProperties, turbulenceProperties")

    # ---- fluid system/ ------------------------------------------------------
    _write(os.path.join(fluid_out, "system", "controlDict"),
           control_dict(args.end_time, args.delta_t, args.write_interval))
    _write(os.path.join(fluid_out, "system", "preciceDict"), PRECICE_DICT)
    _write(os.path.join(fluid_out, "system", "fvSchemes"), FV_SCHEMES)
    _write(os.path.join(fluid_out, "system", "fvSolution"), FV_SOLUTION)
    _write(os.path.join(fluid_out, "system", "setFieldsDict"), SET_FIELDS_DICT)
    _write(os.path.join(fluid_out, "system", "decomposeParDict"),
           decompose_par_dict(args.nprocs))
    log("wrote system/controlDict (deltaT=%s, endTime=%s), preciceDict, "
        "fvSchemes, fvSolution, setFieldsDict, decomposeParDict"
        % (_fmt_num(args.delta_t), _fmt_num(args.end_time)))

    # ---- fluid 0/ ------------------------------------------------------------
    fields = {
        "alpha.water": ALPHA_WATER_0,
        "U": U_0,
        "p_rgh": P_RGH_0,
        "p": P_0,
        "k": K_0,
        "omega": OMEGA_0,
        "nut": NUT_0,
        "pointDisplacement": POINT_DISPLACEMENT_0,
    }
    for name, text in fields.items():
        _write(os.path.join(fluid_out, "0", name), text)
    log("wrote 0/{U, alpha.water, p_rgh, p, k, omega, nut, pointDisplacement}")

    # ---- fluid case.foam (ParaView anchor) -----------------------------------
    _write(os.path.join(fluid_out, "case.foam"),
           "OpenFOAM case marker for ParaView (assembled by 03_assemble_precice_case.py)\n")
    log("wrote fluid-openfoam/case.foam")

    # ---- solid config ---------------------------------------------------------
    _write(os.path.join(solid_out, "config.yml"), CONFIG_YML)
    _write(os.path.join(solid_out, "run.sh"), RUN_SH)
    log("wrote solid-calculix/config.yml + run.sh")

    # ---- summary ---------------------------------------------------------------
    log("=== summary ===")
    log("case root     : %s" % out_dir)
    log("fluid solver  : interFoam (deltaT %s, endTime %s, writeInterval %s)"
        % (_fmt_num(args.delta_t), _fmt_num(args.end_time),
           _fmt_num(args.write_interval)))
    log("solid deck    : dam_gate.inp (E=%s nu=%s rho=%s)"
        % (_fmt_num(args.gate_e), _fmt_num(args.gate_nu), _fmt_num(args.gate_rho)))
    log("precice       : time-window %s, max-time %s, parallel-implicit IQN-ILS"
        % (_fmt_num(args.delta_t), _fmt_num(args.precice_max_time)))
    log("next step     : python3 04_run_precice.py --case %s [--end-time 0.02]"
        % out_dir)
    log("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
