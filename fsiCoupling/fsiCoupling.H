/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Shared FSI coupling controls and diagnostics for local solver variants.

\*---------------------------------------------------------------------------*/

#ifndef fsiCoupling_H
#define fsiCoupling_H

#include "fvCFD.H"
#include "pimpleControl.H"
#include "uniformDimensionedFields.H"

class fsiPimpleControl : public pimpleControl
{
    const bool& fsiConverged_;

public:
    fsiPimpleControl
    (
        fvMesh& mesh,
        const bool& fsiConverged,
        const word& dictName = "PIMPLE"
    )
    :
        pimpleControl(mesh, dictName),
        fsiConverged_(fsiConverged)
    {}

protected:
    virtual bool criteriaSatisfied()
    {
        return pimpleControl::criteriaSatisfied() && fsiConverged_;
    }
};


struct fsiCouplingControls
{
    word coupling;
    label maxIter;
    scalar tolerance;
};


inline fsiCouplingControls readFsiCouplingControls
(
    const Time& runTime,
    fvMesh& mesh
)
{
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    fsiCouplingControls controls;
    controls.coupling =
        dynamicMeshDict.lookupOrDefault<word>("fsiCoupling", "explicit");
    controls.maxIter =
        dynamicMeshDict.lookupOrDefault<label>("maxFsiIter", 10);
    controls.tolerance =
        dynamicMeshDict.lookupOrDefault<scalar>("fsiTolerance", 1e-4);

    if
    (
        controls.coupling != "explicit"
     && controls.coupling != "partitioned"
     && controls.coupling != "integrated"
     && controls.coupling != "implicit"
    )
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Unsupported fsiCoupling '" << controls.coupling
            << "'. Valid values are explicit, partitioned, integrated, "
            << "and implicit." << exit(FatalIOError);
    }

    return controls;
}


inline bool isPartitionedFsiCoupling(const word& coupling)
{
    return coupling == "partitioned";
}


inline bool isIntegratedFsiCoupling(const word& coupling)
{
    return coupling == "integrated" || coupling == "implicit";
}


inline scalar readFsiResidual(const fvMesh& mesh)
{
    if (mesh.foundObject<uniformDimensionedScalarField>("fsiResidual"))
    {
        return mesh.lookupObject<uniformDimensionedScalarField>
        (
            "fsiResidual"
        ).value();
    }

    return 1.0;
}


inline scalar fsiMaxDisplacementChange
(
    const fvMesh& mesh,
    const pointField& previousPoints,
    const label iter,
    const word& iterationLabel
)
{
    const pointField& points = mesh.points();
    label topologyChanged = points.size() == previousPoints.size() ? 0 : 1;
    reduce(topologyChanged, maxOp<label>());

    scalar maxDiff = 0.0;
    if (!topologyChanged)
    {
        maxDiff = max(mag(points - previousPoints));
    }

    reduce(maxDiff, maxOp<scalar>());

    if (topologyChanged && Pstream::master())
    {
        Info << iterationLabel << " " << iter
             << ": topology change detected on at least one rank"
             << "; reset displacement-change reference." << endl;
    }

    return maxDiff;
}


inline void logFsiIteration
(
    const word& iterationLabel,
    const label iter,
    const scalar fsiResidual,
    const scalar maxDisplacementChange
)
{
    Info << iterationLabel << " " << iter
         << ": FSI Force Residual = " << fsiResidual
         << ", max displacement change = " << maxDisplacementChange << endl;
}

#endif

// ************************************************************************* //
