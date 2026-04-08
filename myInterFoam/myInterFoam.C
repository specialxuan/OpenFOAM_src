/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    interFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// Custom PIMPLE control to enforce FSI convergence
class fsiPimpleControl : public pimpleControl
{
    const bool& fsiConverged_;
public:
    fsiPimpleControl(fvMesh& mesh, const bool& fsiConverged, const word& dictName="PIMPLE")
    : pimpleControl(mesh, dictName), fsiConverged_(fsiConverged)
    {}

protected:
    virtual bool criteriaSatisfied()
    {
        return pimpleControl::criteriaSatisfied() && fsiConverged_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

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
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop

        // Read FSI controls
        IOdictionary dynamicMeshDict(
            IOobject(
                "dynamicMeshDict",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
        const word fsiCoupling = dynamicMeshDict.lookupOrDefault<word>("fsiCoupling", "explicit");
        const label maxFsiIter = dynamicMeshDict.lookupOrDefault<label>("maxFsiIter", 10);
        const scalar fsiTolerance = dynamicMeshDict.lookupOrDefault<scalar>("fsiTolerance", 1e-4);

        if (fsiCoupling == "partitioned")
        {
            label fsiIter = 0;
            bool fsiConverged = false;
            pointField points0 = mesh.points();

            while (fsiIter < maxFsiIter && !fsiConverged)
            {
                fsiIter++;

                // Shadow pimple control for this FSI iteration
                pimpleControl pimple(mesh);

                // Update local controls based on new pimple dict
                #include "readDyMControls.H"
                moveMeshOuterCorrectors = false; // Force fixed mesh during PIMPLE sub-loop

                while (pimple.loop())
                {
                    #include "pimpleLoopBody.H"
                }

                // Check convergence
                const pointField& points = mesh.points();
                scalar maxDiff = 0.0;
                if (points.size() == points0.size())
                {
                    maxDiff = max(mag(points - points0));
                }
                else if (Pstream::master())
                {
                    Info << "FSI Iteration " << fsiIter
                         << ": topology change detected on at least one rank"
                         << "; reset displacement-change reference." << endl;
                }
                reduce(maxDiff, maxOp<scalar>());
                scalar fsiRes = 1.0;
                if (mesh.foundObject<uniformDimensionedScalarField>("fsiResidual"))
                {
                    fsiRes = mesh.lookupObject<uniformDimensionedScalarField>("fsiResidual").value();
                }
                fsiConverged = (fsiRes < fsiTolerance);
                points0 = points;

                Info << "FSI Iteration " << fsiIter
                     << ": FSI Force Residual = " << fsiRes
                     << ", max displacement change = " << maxDiff << endl;
            }
        }
        else if (fsiCoupling == "integrated" || fsiCoupling == "implicit") // Option (3)
        {
            // Integrated implicit coupling: Move mesh on every PIMPLE outer corrector
            
            pointField points0 = mesh.points();
            bool fsiConverged = false;
            label fsiIter = 0;
            fsiPimpleControl pimple(mesh, fsiConverged);

            while (pimple.loop())
            {
                fsiIter++;
                #include "readDyMControls.H" // Updates standard flags
                moveMeshOuterCorrectors = true; // Force mesh motion every PIMPLE loop

                #include "pimpleLoopBody.H"

                const pointField& points = mesh.points();
                scalar maxDiff = 0.0;
                if (points.size() == points0.size())
                {
                    maxDiff = max(mag(points - points0));
                }
                else if (Pstream::master())
                {
                    Info << "FSI/PIMPLE Iteration " << fsiIter
                         << ": topology change detected on at least one rank"
                         << "; reset displacement-change reference." << endl;
                }
                reduce(maxDiff, maxOp<scalar>());
                
                scalar fsiRes = 1.0;
                if (mesh.foundObject<uniformDimensionedScalarField>("fsiResidual"))
                {
                    fsiRes = mesh.lookupObject<uniformDimensionedScalarField>("fsiResidual").value();
                }

                Info << "FSI/PIMPLE Iteration " << fsiIter
                     << ": max displacement change = " << maxDiff 
                     << ", FSI Force Residual = " << fsiRes << endl;

                // Check for FSI convergence (using force residual)
                fsiConverged = (fsiRes < fsiTolerance);
                points0 = points;
            }

            if (!fsiConverged)
            {
                 Info << "Warning: FSI force residual did not converge within PIMPLE iterations." << endl;
            }
        }
        else // "explicit" (Option 1) or default
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

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
