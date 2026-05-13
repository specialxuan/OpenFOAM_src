if (pimple.firstIter() || moveMeshOuterCorrectors)
{
    const Switch makeVelocityAmrIndicator =
        runTime.controlDict().lookupOrDefault<Switch>
        (
            "makeVelocityAmrIndicator",
            false
        );

    autoPtr<volScalarField> velocityAmrIndicatorPtr;

    if (makeVelocityAmrIndicator)
    {
        const word velocityAmrIndicatorName =
            runTime.controlDict().lookupOrDefault<word>
            (
                "velocityAmrIndicatorName",
                "amrVortOrGradU"
            );

        const tmp<volTensorField> tGradU(fvc::grad(U));
        const tmp<volVectorField> tVorticity(fvc::curl(U));

        velocityAmrIndicatorPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    velocityAmrIndicatorName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                max(mag(tGradU()), mag(tVorticity()))
            )
        );
    }

    // Do any mesh changes
    mesh.controlledUpdate();

    if (mesh.changing())
    {
        MRF.update();

        if (correctPhi)
        {
            // Calculate absolute flux
            // from the mapped surface velocity
            phi = mesh.Sf() & Uf();

#include "correctPhi.H"

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);
        }

        if (checkMeshCourantNo)
        {
#include "meshCourantNo.H"
        }
    }
}

#include "UEqn.H"

// --- Pressure corrector loop
while (pimple.correct())
{
#include "pEqn.H"
}

if (pimple.turbCorr())
{
    laminarTransport.correct();
    turbulence->correct();
}
