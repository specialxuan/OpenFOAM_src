/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    fastDynamicFvMesh AMR and topology-mapping implementation.

\*---------------------------------------------------------------------------*/

#include "fastDynamicFvMesh.H"
#include "Pstream.H"
#include "fvcGrad.H"
#include "mapPolyMesh.H"

namespace Foam
{

scalar fastDynamicFvMesh::cellMinEdgeLength(const label celli) const
{
    const cell& cFaces = this->cells()[celli];
    const faceList& fcs = this->faces();
    const pointField& pts = this->points();
    scalar minLen = GREAT;

    forAll(cFaces, cfi)
    {
        const face& f = fcs[cFaces[cfi]];
        const label nfp = f.size();

        for (label i = 0; i < nfp; ++i)
        {
            const point& p0 = pts[f[i]];
            const point& p1 = pts[f[(i + 1) % nfp]];
            minLen = min(minLen, mag(p1 - p0));
        }
    }

    return minLen;
}


bool fastDynamicFvMesh::cellPassesRefinementSizeFloor(const label celli) const
{
    if (refinementMinCellVolume_ <= 0 && refinementMinEdgeLength_ <= 0)
    {
        return true;
    }

    const scalar childVolume = this->V()[celli]/8.0;
    if (refinementMinCellVolume_ > 0 && childVolume < refinementMinCellVolume_)
    {
        return false;
    }

    if (refinementMinEdgeLength_ > 0)
    {
        const scalar minParentEdge = cellMinEdgeLength(celli);
        if (minParentEdge <= VSMALL)
        {
            return false;
        }

        if (0.5*minParentEdge < refinementMinEdgeLength_)
        {
            return false;
        }
    }

    return true;
}


scalarField fastDynamicFvMesh::refinementIndicatorCellField
(
    const scalarField& fallbackCellField
) const
{
    if (!refinementUseGradIndicator_)
    {
        return fallbackCellField;
    }

    if (!this->foundObject<volScalarField>(refinementGradIndicatorField_))
    {
        FatalErrorInFunction
            << "Runtime AMR gradient indicator is enabled, but source field '"
            << refinementGradIndicatorField_
            << "' is not available in the objectRegistry."
            << exit(FatalError);
    }

    const volScalarField& sourceField =
        this->lookupObject<volScalarField>(refinementGradIndicatorField_);
    const tmp<volScalarField> tIndicator = mag(fvc::grad(sourceField));
    return tIndicator().primitiveField();
}


scalarField fastDynamicFvMesh::refinementIndicatorPointField
(
    const scalarField& fallbackPointField
) const
{
    if (!refinementUseGradIndicator_)
    {
        return fallbackPointField;
    }

    if (!this->foundObject<volScalarField>(refinementGradIndicatorField_))
    {
        FatalErrorInFunction
            << "Runtime AMR gradient indicator is enabled, but source field '"
            << refinementGradIndicatorField_
            << "' is not available in the objectRegistry."
            << exit(FatalError);
    }

    const volScalarField& sourceField =
        this->lookupObject<volScalarField>(refinementGradIndicatorField_);
    const tmp<volScalarField> tIndicator = mag(fvc::grad(sourceField));
    return maxCellField(tIndicator());
}


void fastDynamicFvMesh::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    bitSet& candidateCell
) const
{
    const scalarField indicatorCell = refinementIndicatorCellField(vFld);

    dynamicRefineFvMesh::selectRefineCandidates
    (
        lowerRefineLevel,
        upperRefineLevel,
        indicatorCell,
        candidateCell
    );

    if (refinementMinCellVolume_ <= 0 && refinementMinEdgeLength_ <= 0)
    {
        return;
    }

    for (const label celli : candidateCell)
    {
        if (!cellPassesRefinementSizeFloor(celli))
        {
            candidateCell.unset(celli);
        }
    }
}


labelList fastDynamicFvMesh::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const bitSet& candidateCell
) const
{
    if (refinementMinCellVolume_ <= 0 && refinementMinEdgeLength_ <= 0)
    {
        return dynamicRefineFvMesh::selectRefineCells
        (
            maxCells,
            maxRefinement,
            candidateCell
        );
    }

    bitSet sizeFloorFiltered(candidateCell);
    label nRejected = 0;
    const label nCandidates = candidateCell.count();

    for (const label celli : candidateCell)
    {
        if (!cellPassesRefinementSizeFloor(celli))
        {
            sizeFloorFiltered.unset(celli);
            ++nRejected;
        }
    }

    const label totalRejected = returnReduce(nRejected, sumOp<label>());
    const label totalCandidates = returnReduce(nCandidates, sumOp<label>());
    if (Pstream::master() && totalRejected > 0)
    {
        Info<< "Refinement size-floor filter rejected " << totalRejected
            << " candidate cells (refinementMinCellVolume="
            << refinementMinCellVolume_
            << ", refinementMinEdgeLength=" << refinementMinEdgeLength_
            << ")." << endl;
    }
    else if (Pstream::master() && totalCandidates > 0)
    {
        Info<< "Refinement size-floor filter accepted all "
            << totalCandidates << " candidate cells." << endl;
    }

    return dynamicRefineFvMesh::selectRefineCells
    (
        maxCells,
        maxRefinement,
        sizeFloorFiltered
    );
}


labelList fastDynamicFvMesh::selectUnrefinePoints
(
    const scalar unrefineLevel,
    const bitSet& markedCell,
    const scalarField& pFld
) const
{
    const scalarField indicatorPoint = refinementIndicatorPointField(pFld);

    return dynamicRefineFvMesh::selectUnrefinePoints
    (
        unrefineLevel,
        markedCell,
        indicatorPoint
    );
}


void fastDynamicFvMesh::updateMesh(const mapPolyMesh& mpm)
{
    dynamicFvMesh::updateMesh(mpm);

    if (!meshRefinementSupport_)
    {
        return;
    }

    ensureModeShapesForCurrentMesh(mpm);
    buildForceProjectionCaches(&mpm);
}

} // namespace Foam

