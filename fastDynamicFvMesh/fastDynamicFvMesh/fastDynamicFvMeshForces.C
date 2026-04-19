/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Fast Dynamic Mesh method implementation.

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"
#include "fastDynamicFvMesh.H"
#include "fluidThermo.H"
#include "fvcGrad.H"
#include "mapPolyMesh.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "transportModel.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"
#include "Time.H"
#include "regIOobject.H"
#include "volFields.H"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
void fastDynamicFvMesh::cacheModelPointers()
{
    if (modelPointersCached_)
    {
        return;
    }

    typedef incompressible::turbulenceModel icoTurbModel;
    typedef compressible::turbulenceModel cmpTurbModel;

    if (this->foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        icoTurbPtr_ =
            &this->lookupObject<icoTurbModel>(icoTurbModel::propertiesName);
    }

    if (this->foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        cmpTurbPtr_ =
            &this->lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);
    }

    if (this->foundObject<fluidThermo>(fluidThermo::dictName))
    {
        fluidThermoPtr_ =
            &this->lookupObject<fluidThermo>(fluidThermo::dictName);
    }

    if (this->foundObject<transportModel>("transportProperties"))
    {
        laminarTransportPtr_ =
            &this->lookupObject<transportModel>("transportProperties");
    }

    modelPointersCached_ = true;
}

vector fastDynamicFvMesh::interpolateStructuralForce(
    const scalarField& times, const vectorField& values, const scalar t) const
{
    if (times.size() == 1)
    {
        return values[0];
    }

    if (t <= times.first())
    {
        return values.first();
    }

    if (t >= times.last())
    {
        return values.last();
    }

    for (label i = 0; i < times.size() - 1; ++i)
    {
        const scalar t0 = times[i];
        const scalar t1 = times[i + 1];

        if (t >= t0 && t <= t1)
        {
            const scalar dt = t1 - t0;
            if (dt <= VSMALL)
            {
                return values[i + 1];
            }

            const scalar w = (t - t0) / dt;
            return (1.0 - w) * values[i] + w * values[i + 1];
        }
    }

    return values.last();
}

// patchDensity(patchi, defaultRho)
// Purpose: return a scalarField of density values corresponding to the
// requested patch. Behavior:
//  - If a volScalarField named "rho" exists in the case, the function returns a
//  copy of its boundary field for patchi
//    which may contain per-face or per-cell densities produced by the solver.
//  - Otherwise, a scalarField filled with defaultRho is returned with length
//  equal to the number of faces in the patch.
// Rationale: kinematic pressure fields (p/rho) require a density to convert
// back to dimensional pressure for force computations.
tmp<scalarField> fastDynamicFvMesh::patchDensity(
    const label patchi, const scalar defaultRho) const
{
    if (this->foundObject<volScalarField>("rho"))
    {
        const auto& rhoField = this->lookupObject<volScalarField>("rho");
        return tmp<scalarField>(
            new scalarField(rhoField.boundaryField()[patchi]));
    }

    return tmp<scalarField>(
        new scalarField(this->boundaryMesh()[patchi].size(), defaultRho));
}

// devRhoReff(gradUp, patchi, defaultRho)
// Purpose: compute the deviatoric effective stress tensor for wall-shear
// evaluation on a patch. Implementation notes:
//  - Prefer available turbulence/transport models in order to obtain an
//  effective kinematic viscosity (nuEff or muEff).
//  - When density-dependent viscosity or kinematic-pressure correction is
//  required, obtain per-face density with patchDensity().
//  - The returned tensor field represents -rho * nu_eff * devTwoSymm(gradUp),
//  i.e., the deviatoric viscous stress contribution per face.
//  - This helper centralizes model queries so calcModalForces can obtain shear
//  contributions regardless of solver setup.
tmp<symmTensorField> fastDynamicFvMesh::devRhoReff(const tensorField& gradUp,
    const label patchi, const scalar defaultRho) const
{
    if (icoTurbPtr_)
    {
        const auto& turb = *icoTurbPtr_;

        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>(new symmTensorField(
            -tRho() * turb.nuEff(patchi) * devTwoSymm(gradUp)));
    }
    else if (cmpTurbPtr_)
    {
        const auto& turb = *cmpTurbPtr_;

        return tmp<symmTensorField>(
            new symmTensorField(-turb.muEff(patchi) * devTwoSymm(gradUp)));
    }
    else if (fluidThermoPtr_)
    {
        const auto& thermo = *fluidThermoPtr_;

        return tmp<symmTensorField>(
            new symmTensorField(-thermo.mu(patchi) * devTwoSymm(gradUp)));
    }
    else if (laminarTransportPtr_)
    {
        const auto& laminarT = *laminarTransportPtr_;

        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>(new symmTensorField(
            -tRho() * laminarT.nu(patchi) * devTwoSymm(gradUp)));
    }
    else
    {
        IOdictionary transportProperties(
            IOobject("transportProperties", this->time().constant(), *this,
                IOobject::MUST_READ, IOobject::NO_WRITE));

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);
        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>(
            new symmTensorField(-tRho() * nu.value() * devTwoSymm(gradUp)));
    }
}

// calcModalForces()
// Purpose: compute per-mode generalized forces by projecting fluid pressure and
// wall-shear traction onto mode shapes. Detailed flow:
//  - Zero modal accumulators (modeForce_, modePressureForce_, modeShearForce_).
//  - Fetch the pressure field (supports dimensional pressure or kinematic
//  pressure). If kinematic, determine density (rhoRef or rho field).
//  - Optionally open per-processor face diagnostics CSV (append mode) and write
//  header if first time.
//  - If velocity field 'U' exists, compute its gradient once (fvc::grad(U)) to
//  evaluate viscous stress.
//  - For each configured FSI patch:
//      * Verify patch exists and obtain its face areas and centres.
//      * Compute per-face pressure force = pressure * areaVec (scaled by
//      density for kinematic p).
//      * Compute deviatoric viscous stress and then shear force; remove normal
//      component to leave tangential shear only.
//      * Interpolate mode shape to face by averaging shape vectors at face
//      points.
//      * Project pressure and shear vectors onto the shape vector (dot product)
//      and accumulate into pressure/shear/modal totals.
//      * Optionally write per-face diagnostics lines for the requested mode.
//  - Perform parallel reductions to combine per-processor modal accumulators
//  and (on master) print debug ranges.
void fastDynamicFvMesh::calcModalForces()
{
    cacheModelPointers();

    // Initialize forces
    modeForce_ = 0.0;
    modePressureForce_ = 0.0;
    modeShearForce_ = 0.0;

    if (!this->foundObject<volScalarField>(pressureFieldName_))
    {
        FatalErrorInFunction
            << "Required pressure field '" << pressureFieldName_
            << "' was not found. Configure 'pressureField' in dynamicMeshDict "
            << "or supply the field before fastDynamicFvMesh::update()."
            << exit(FatalError);
    }

    const volScalarField& p =
        this->lookupObject<volScalarField>(pressureFieldName_);
    const dimensionSet kinematicPressureDims(dimPressure / dimDensity);
    const bool dimensionalPressure = (p.dimensions() == dimPressure);
    const bool kinematicPressure = (p.dimensions() == kinematicPressureDims);

    if (!dimensionalPressure && !kinematicPressure)
    {
        FatalErrorInFunction << "Pressure field '" << pressureFieldName_
                             << "' has dimensions " << p.dimensions()
                             << ". Expected either " << dimPressure
                             << " (pressure) or " << kinematicPressureDims
                             << " (kinematic pressure)." << exit(FatalError);
    }

    scalar defaultRho = rhoRef_;

    if (kinematicPressure && defaultRho <= 0)
    {
        IOdictionary transportProperties(
            IOobject("transportProperties", this->time().constant(), *this,
                IOobject::MUST_READ, IOobject::NO_WRITE));

        if (transportProperties.found("rho"))
        {
            transportProperties.lookup("rho") >> defaultRho;
        }
    }

    if (kinematicPressure && defaultRho <= 0 &&
        !this->foundObject<volScalarField>("rho"))
    {
        FatalErrorInFunction
            << "Pressure field '" << pressureFieldName_
            << "' is kinematic, but no density information was found." << nl
            << "Provide 'rhoRef' in dynamicMeshDict, a 'rho' entry in "
            << "constant/transportProperties, or a volScalarField named 'rho'."
            << exit(FatalError);
    }

    tmp<volTensorField> tGradU;
    const volTensorField* gradUPtr = nullptr;

    if (this->foundObject<volVectorField>("U"))
    {
        const volVectorField& U = this->lookupObject<volVectorField>("U");
        tGradU = fvc::grad(U);
        gradUPtr = &tGradU();
    }
    else
    {
        WarningInFunction
            << "Velocity field 'U' not found; wall shear contribution "
            << "to modal forces will be skipped." << endl;
    }

    // Iterate over cached FSI patches
    forAll(fsiPatchIDs_, i)
    {
        const label patchID = fsiPatchIDs_[i];
        const polyPatch& pp = *fsiPolyPatches_[i];
        const fvPatchScalarField& pPatch = p.boundaryField()[patchID];
        const vectorField faceAreas(pp.faceAreas());
        scalarField& pressureScale = pressureScaleCache_[i];
        if (kinematicPressure)
        {
            pressureScale = patchDensity(patchID, defaultRho);
            pressureScaleInitialized_[i] = true;
        }
        else if (!pressureScaleInitialized_[i])
        {
            pressureScale = 1.0;
            pressureScaleInitialized_[i] = true;
        }

        const symmTensorField* devStressPtr = nullptr;
        tmp<symmTensorField> tDevStress;

        if (gradUPtr)
        {
            const tensorField& gradPatch = gradUPtr->boundaryField()[patchID];
            tDevStress = devRhoReff(
                gradPatch, patchID, defaultRho > 0 ? defaultRho : 1.0);
            devStressPtr = &tDevStress();
        }

        if (meshRefinementSupport_ && referenceFsiBuilt_ &&
            referenceFaceModeProjection_.size() == fsiPatches_.size())
        {
            const label nRefFaces = referenceFsiFaceAreas_[i].size();
            vectorField aggregatedPressure(nRefFaces, vector::zero);
            vectorField aggregatedShear(nRefFaces, vector::zero);
            const List<labelList>& faceToRefs = fsiFaceToReferenceFaces_[i];
            const List<scalarField>& faceToRefWeights =
                fsiFaceToReferenceWeights_[i];

            // Loop over current (possibly split) faces and aggregate to
            // reference faces
            forAll(pp, faceI)
            {
                if (faceI >= faceToRefs.size() ||
                    faceI >= faceToRefWeights.size())
                {
                    FatalErrorInFunction
                        << "Invalid refined-face mapping: patch '"
                        << fsiPatches_[i] << "', face " << faceI
                        << " has no reference mapping entry."
                        << exit(FatalError);
                }

                const labelList& refs = faceToRefs[faceI];
                const scalarField& weights = faceToRefWeights[faceI];

                if (refs.size() == 0 || refs.size() != weights.size())
                {
                    FatalErrorInFunction
                        << "Invalid refined/unrefined face mapping on patch '"
                        << fsiPatches_[i] << "', face " << faceI
                        << ": refs=" << refs.size()
                        << ", weights=" << weights.size() << exit(FatalError);
                }

                const vector& areaVec = faceAreas[faceI];
                const scalar areaMag = mag(areaVec);
                const vector faceNormal =
                    areaMag > VSMALL ? areaVec / areaMag : vector::zero;

                const vector pressureForce =
                    (pressureScale[faceI] * pPatch[faceI] - pRef_) * areaVec;
                vector shearForce = vector::zero;

                if (devStressPtr)
                {
                    shearForce = areaVec & (*devStressPtr)[faceI];
                    if (areaMag > VSMALL)
                    {
                        shearForce -= faceNormal * (faceNormal & shearForce);
                    }
                }

                const vector pressureTraction =
                    areaMag > VSMALL ? pressureForce / areaMag : vector::zero;
                const vector shearTraction =
                    areaMag > VSMALL ? shearForce / areaMag : vector::zero;

                forAll(refs, refI)
                {
                    const label refFaceI = refs[refI];
                    if (refFaceI < 0 || refFaceI >= nRefFaces)
                    {
                        FatalErrorInFunction
                            << "Invalid refined/unrefined mapping: patch '"
                            << fsiPatches_[i] << "', face " << faceI
                            << " -> ref " << refFaceI
                            << ", nRefFaces=" << nRefFaces << exit(FatalError);
                    }

                    const scalar w = max(weights[refI], 0.0);
                    const scalar refArea = referenceFsiFaceAreas_[i][refFaceI];
                    aggregatedPressure[refFaceI] +=
                        pressureTraction * (w * refArea);
                    aggregatedShear[refFaceI] += shearTraction * (w * refArea);
                }
            }

            for (label m = 0; m < nMode_; ++m)
            {
                const vectorField& refShape =
                    referenceFaceModeProjection_[i][m];
                forAll(refShape, refFaceI)
                {
                    const scalar pressureModeForce =
                        aggregatedPressure[refFaceI] & refShape[refFaceI];
                    const scalar shearModeForce =
                        aggregatedShear[refFaceI] & refShape[refFaceI];

                    modePressureForce_[m] += pressureModeForce;
                    modeShearForce_[m] += shearModeForce;
                    modeForce_[m] += pressureModeForce + shearModeForce;
                }
            }
        }
        else
        {
            // Loop over faces (legacy path)
            forAll(pp, faceI)
            {
                const vector& areaVec = faceAreas[faceI];
                const vector pressureForce =
                    (pressureScale[faceI] * pPatch[faceI] - pRef_) * areaVec;
                vector shearForce = vector::zero;

                if (devStressPtr)
                {
                    shearForce = areaVec & (*devStressPtr)[faceI];
                    const scalar areaMag = mag(areaVec);

                    if (areaMag > VSMALL)
                    {
                        const vector faceNormal = areaVec / areaMag;
                        shearForce -= faceNormal * (faceNormal & shearForce);
                    }
                }

                for (label m = 0; m < nMode_; ++m)
                {
                    const vector& shapeFace = faceModeProjection_[i][m][faceI];

                    const scalar pressureModeForce = pressureForce & shapeFace;
                    const scalar shearModeForce = shearForce & shapeFace;

                    modePressureForce_[m] += pressureModeForce;
                    modeShearForce_[m] += shearModeForce;
                    modeForce_[m] += pressureModeForce + shearModeForce;
                }
            }
        }
    }

    // Parallel reduction (all-reduce style)
    forAll(modeForce_, i)
    {
        reduce(modeForce_[i], sumOp<scalar>());
        reduce(modePressureForce_[i], sumOp<scalar>());
        reduce(modeShearForce_[i], sumOp<scalar>());
    }
}

void fastDynamicFvMesh::calcStructuralModalForces()
{
    structuralModeForce_ = 0.0;
    structuralForceSignal_ = 0.0;

    if (!structuralForceEnabled_ || structuralForceNodeIDs_.size() == 0)
    {
        return;
    }

    if (Pstream::master())
    {
        const scalar t = this->time().value();

        forAll(structuralForceNodeIDs_, forceI)
        {
            const vector nodalForce =
                interpolateStructuralForce(structuralForceTimes_[forceI],
                    structuralForceValues_[forceI], t);

            structuralForceSignal_ += mag(nodalForce);

            const label nTarget = structuralForceNodeIDs_[forceI].size();
            if (nTarget == 0)
            {
                continue;
            }

            forAll(structuralForceNodeIDs_[forceI], localI)
            {
                const label nodeI = structuralForceNodeIDs_[forceI][localI];

                for (label modeI = 0; modeI < nMode_; ++modeI)
                {
                    structuralModeForce_[modeI] +=
                        nodalForce & structModeShapes_[modeI][nodeI];
                }
            }
        }
    }

    Pstream::broadcast(structuralModeForce_);
    Pstream::broadcast(structuralForceSignal_);
}

} // End namespace Foam
