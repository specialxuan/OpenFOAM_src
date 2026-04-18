/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fsiCoupling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote("Transient solver for incompressible, turbulent flow"
                     " of Newtonian fluids on a moving mesh.");

#include "postProcess.H"

#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createDynamicFvMesh.H"
#include "initContinuityErrs.H"
#include "createDyMControls.H"
#include "createFields.H"
#include "createUfIfPresent.H"
#include "CourantNo.H"
#include "setInitialDeltaT.H"

    turbulence->validate();

    if (!LTS)
    {
#include "CourantNo.H"
#include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#include "readDyMControls.H"

        if (LTS)
        {
#include "setRDeltaT.H"
        }
        else
        {
#include "CourantNo.H"
#include "setDeltaT.H"
        }

        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        const fsiCouplingControls fsiControls =
            readFsiCouplingControls(runTime, mesh);

        if (isPartitionedFsiCoupling(fsiControls.coupling))
        {
            label fsiIter = 0;
            bool fsiConverged = false;
            pointField points0 = mesh.points();

            while (fsiIter < fsiControls.maxIter && !fsiConverged)
            {
                fsiIter++;

                // Shadow pimple control for this FSI iteration
                pimpleControl pimple(mesh);

// Update local controls based on new pimple dict
#include "readDyMControls.H"
                moveMeshOuterCorrectors =
                    false; // Force fixed mesh during PIMPLE sub-loop

                while (pimple.loop())
                {
#include "pimpleLoopBody.H"
                }

                const scalar maxDiff = fsiMaxDisplacementChange
                (
                    mesh,
                    points0,
                    fsiIter,
                    "FSI Iteration"
                );
                const scalar fsiRes = readFsiResidual(mesh);
                fsiConverged = (fsiRes < fsiControls.tolerance);
                points0 = mesh.points();

                logFsiIteration("FSI Iteration", fsiIter, fsiRes, maxDiff);
            }
        }
        else if (isIntegratedFsiCoupling(fsiControls.coupling))
        {
            // Use custom PIMPLE loop to enforce FSI convergence
            bool fsiConverged = false;
            label fsiIter = 0;
            fsiPimpleControl pimple(mesh, fsiConverged);

            pointField points0 = mesh.points();

            while (pimple.loop())
            {
                fsiIter++;
#include "readDyMControls.H"
                moveMeshOuterCorrectors = true;

#include "pimpleLoopBody.H"

                const scalar maxDiff = fsiMaxDisplacementChange
                (
                    mesh,
                    points0,
                    fsiIter,
                    "FSI/PIMPLE Iteration"
                );
                const scalar fsiRes = readFsiResidual(mesh);
                points0 = mesh.points();

                logFsiIteration
                (
                    "FSI/PIMPLE Iteration",
                    fsiIter,
                    fsiRes,
                    maxDiff
                );

                fsiConverged = (fsiRes < fsiControls.tolerance);
            }

            if (!fsiConverged)
            {
                 Info << "Warning: FSI force residual did not converge within PIMPLE iterations." << endl;
            }
        }
        else
        {
            pimpleControl pimple(mesh);
            while (pimple.loop())
            {
#include "readDyMControls.H"
#include "pimpleLoopBody.H"
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
