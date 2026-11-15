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
#include <fstream>
#include <iomanip>

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

    const scalar childVolume = this->V()[celli] / 8.0;
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

        if (0.5 * minParentEdge < refinementMinEdgeLength_)
        {
            return false;
        }
    }

    return true;
}

void fastDynamicFvMesh::resetAmrLayerDiagnosticsState() const
{
    amrLayerDiagnosticsCandidateCell_.clear();
    amrLayerDiagnosticsIndicatorCell_.setSize(0);
    amrLayerDiagnosticsLower_ = 0.0;
    amrLayerDiagnosticsUpper_ = 0.0;
    amrLayerDiagnosticsLayerCells_ = 0;
    amrLayerDiagnosticsNearbyCells_ = 0;
    amrLayerDiagnosticsCandidates_ = 0;
    amrLayerDiagnosticsSelected_ = 0;
    amrLayerDiagnosticsBelowLower_ = 0;
    amrLayerDiagnosticsAboveUpper_ = 0;
    amrLayerDiagnosticsRejectedSizeFloor_ = 0;
    amrLayerDiagnosticsRejectedBase_ = 0;
    amrLayerDiagnosticsUnrefinePoints_ = 0;
    amrLayerDiagnosticsPendingSummary_ = false;
}

void fastDynamicFvMesh::captureAmrLayerDiagnostics(
    const scalar lowerRefineLevel, const scalar upperRefineLevel,
    const label maxCells, const label maxRefinement,
    const bitSet& candidateCell, const labelList& selectedCells) const
{
    if (!amrLayerDiagnosticsEnabled_)
    {
        return;
    }

    const label nCells = this->nCells();
    if (!nCells)
    {
        resetAmrLayerDiagnosticsState();
        return;
    }

    const scalarField indicatorCell = refinementIndicatorCellField(this->V());
    const vectorField& centres = this->C();
    const labelList& levels = this->meshCutter().cellLevel();
    const bitSet& diagnosticCandidateCell =
        amrLayerDiagnosticsCandidateCell_.size() == nCells
      ? amrLayerDiagnosticsCandidateCell_
      : candidateCell;

    boolList inLayer(nCells, false);
    boolList selected(nCells, false);
    boolList diagnostic(nCells, false);
    labelList neighborDepth(nCells, -1);

    scalar localNearestDist = GREAT;
    forAll(centres, celli)
    {
        localNearestDist = min(
            localNearestDist, mag(centres[celli].z() - amrLayerDiagnosticsZ_));
    }
    const scalar nearestDist = returnReduce(localNearestDist, minOp<scalar>());
    const scalar layerTol = max(SMALL, 1e-9);

    forAll(centres, celli)
    {
        const scalar zDist = mag(centres[celli].z() - amrLayerDiagnosticsZ_);
        if
        (
            amrLayerDiagnosticsThickness_ > 0
          ? zDist <= amrLayerDiagnosticsThickness_
          : mag(zDist - nearestDist) <= layerTol
        )
        {
            inLayer[celli] = true;
            diagnostic[celli] = true;
            neighborDepth[celli] = 0;
        }
    }

    forAll(selectedCells, i)
    {
        const label celli = selectedCells[i];
        if (celli >= 0 && celli < nCells)
        {
            selected[celli] = true;
            diagnostic[celli] = true;
            if (neighborDepth[celli] < 0)
            {
                neighborDepth[celli] = 0;
            }
        }
    }

    for (const label celli : diagnosticCandidateCell)
    {
        if (celli >= 0 && celli < nCells)
        {
            diagnostic[celli] = true;
            if (neighborDepth[celli] < 0)
            {
                neighborDepth[celli] = 0;
            }
        }
    }

    const labelListList& cellNbrs = this->cellCells();
    for (label depth = 0; depth < amrLayerDiagnosticsNeighborDepth_; ++depth)
    {
        bool changed = false;
        forAll(cellNbrs, celli)
        {
            if (neighborDepth[celli] != depth)
            {
                continue;
            }

            forAll(cellNbrs[celli], nbrI)
            {
                const label nbr = cellNbrs[celli][nbrI];
                if (nbr >= 0 && nbr < nCells && neighborDepth[nbr] < 0)
                {
                    neighborDepth[nbr] = depth + 1;
                    diagnostic[nbr] = true;
                    changed = true;
                }
            }
        }

        if (!changed)
        {
            break;
        }
    }

    fileName outPath = this->time().path()
        / (amrLayerDiagnosticsFilePrefix_ + "_cells_proc"
           + Foam::name(Pstream::myProcNo()) + ".csv");
    if (Pstream::parRun())
    {
        outPath = this->time().path().path()
            / (amrLayerDiagnosticsFilePrefix_ + "_cells_proc"
               + Foam::name(Pstream::myProcNo()) + ".csv");
    }

    std::ofstream out(outPath.c_str(), std::ios::app);
    if (!out.good())
    {
        WarningInFunction << "Unable to open " << outPath
                          << " for writing." << endl;
        resetAmrLayerDiagnosticsState();
        return;
    }

    if (!amrLayerDiagnosticsCellHeaderWritten_)
    {
        out << "time,timeIndex,proc,cellI,centerX,centerY,centerZ,cellLevel,"
               "volume,minEdgeLength,indicator,lower,upper,inBottomLayer,"
               "neighborDepth,baseCandidate,passesSizeFloor,selectedForRefine,"
               "reason\n";
        amrLayerDiagnosticsCellHeaderWritten_ = true;
    }

    out << std::setprecision(12);

    const label cachedUnrefinePoints = amrLayerDiagnosticsUnrefinePoints_;
    resetAmrLayerDiagnosticsState();
    amrLayerDiagnosticsUnrefinePoints_ = cachedUnrefinePoints;
    amrLayerDiagnosticsCandidateCell_ = candidateCell;
    amrLayerDiagnosticsIndicatorCell_ = indicatorCell;
    amrLayerDiagnosticsLower_ = lowerRefineLevel;
    amrLayerDiagnosticsUpper_ = upperRefineLevel;

    label rows = 0;
    forAll(diagnostic, celli)
    {
        if (!diagnostic[celli])
        {
            continue;
        }

        const bool isCandidate = diagnosticCandidateCell.test(celli);
        const bool passesSizeFloor = cellPassesRefinementSizeFloor(celli);
        word reason("selected");

        if (selected[celli])
        {
            reason = "selected";
            ++amrLayerDiagnosticsSelected_;
        }
        else if (!isCandidate)
        {
            if (indicatorCell[celli] < lowerRefineLevel)
            {
                reason = "belowLowerThreshold";
                ++amrLayerDiagnosticsBelowLower_;
            }
            else if (indicatorCell[celli] > upperRefineLevel)
            {
                reason = "aboveUpperThreshold";
                ++amrLayerDiagnosticsAboveUpper_;
            }
            else
            {
                reason = "notCandidateByDynamicRefine";
            }
        }
        else if (!passesSizeFloor)
        {
            const scalar childVolume = this->V()[celli] / 8.0;
            if
            (
                refinementMinCellVolume_ > 0
             && childVolume < refinementMinCellVolume_
            )
            {
                reason = "belowMinCellVolume";
            }
            else
            {
                reason = "belowMinEdgeLength";
            }
            ++amrLayerDiagnosticsRejectedSizeFloor_;
        }
        else if (levels[celli] >= maxRefinement)
        {
            reason = "alreadyAtMaxRefinement";
            ++amrLayerDiagnosticsRejectedBase_;
        }
        else
        {
            reason = "candidateButNotSelectedByBaseAMR";
            ++amrLayerDiagnosticsRejectedBase_;
        }

        if (inLayer[celli])
        {
            ++amrLayerDiagnosticsLayerCells_;
        }
        if (neighborDepth[celli] >= 0)
        {
            ++amrLayerDiagnosticsNearbyCells_;
        }
        if (isCandidate)
        {
            ++amrLayerDiagnosticsCandidates_;
        }

        if
        (
            amrLayerDiagnosticsMaxRows_ <= 0
         || rows < amrLayerDiagnosticsMaxRows_
        )
        {
            out << this->time().value() << ',' << this->time().timeIndex()
                << ',' << Pstream::myProcNo() << ',' << celli << ','
                << centres[celli].x() << ',' << centres[celli].y() << ','
                << centres[celli].z() << ',' << levels[celli] << ','
                << this->V()[celli] << ',' << cellMinEdgeLength(celli) << ','
                << indicatorCell[celli] << ',' << lowerRefineLevel << ','
                << upperRefineLevel << ',' << inLayer[celli] << ','
                << neighborDepth[celli] << ',' << isCandidate << ','
                << passesSizeFloor << ',' << selected[celli] << ',' << reason
                << '\n';

            ++rows;
        }
    }

    amrLayerDiagnosticsPendingSummary_ = true;
}

void fastDynamicFvMesh::writeAmrLayerDiagnosticsSummary() const
{
    if (!amrLayerDiagnosticsEnabled_ || !amrLayerDiagnosticsPendingSummary_)
    {
        return;
    }

    const label nLayer = returnReduce(
        amrLayerDiagnosticsLayerCells_, sumOp<label>());
    const label nNearby = returnReduce(
        amrLayerDiagnosticsNearbyCells_, sumOp<label>());
    const label nCandidates = returnReduce(
        amrLayerDiagnosticsCandidates_, sumOp<label>());
    const label nSelected = returnReduce(
        amrLayerDiagnosticsSelected_, sumOp<label>());
    const label nBelowLower = returnReduce(
        amrLayerDiagnosticsBelowLower_, sumOp<label>());
    const label nAboveUpper = returnReduce(
        amrLayerDiagnosticsAboveUpper_, sumOp<label>());
    const label nRejectedSizeFloor = returnReduce(
        amrLayerDiagnosticsRejectedSizeFloor_, sumOp<label>());
    const label nRejectedBase = returnReduce(
        amrLayerDiagnosticsRejectedBase_, sumOp<label>());
    const label nUnrefinePoints = returnReduce(
        amrLayerDiagnosticsUnrefinePoints_, sumOp<label>());

    if (Pstream::master())
    {
        fileName outPath = this->time().path()
            / (amrLayerDiagnosticsFilePrefix_ + "_summary.csv");
        if (Pstream::parRun())
        {
            outPath = this->time().path().path()
                / (amrLayerDiagnosticsFilePrefix_ + "_summary.csv");
        }

        std::ofstream out(outPath.c_str(), std::ios::app);
        if (!out.good())
        {
            WarningInFunction << "Unable to open " << outPath
                              << " for writing." << endl;
            resetAmrLayerDiagnosticsState();
            return;
        }

        if (!amrLayerDiagnosticsSummaryHeaderWritten_)
        {
            out << "time,timeIndex,nLayerCells,nNearbyCells,nCandidate,"
                   "nSelected,nBelowLower,nAboveUpper,nRejectedSizeFloor,"
                   "nRejectedBase,nUnrefinePoints,lower,upper\n";
            amrLayerDiagnosticsSummaryHeaderWritten_ = true;
        }

        out << std::setprecision(12)
            << this->time().value() << ',' << this->time().timeIndex() << ','
            << nLayer << ',' << nNearby << ',' << nCandidates << ','
            << nSelected << ',' << nBelowLower << ',' << nAboveUpper << ','
            << nRejectedSizeFloor << ',' << nRejectedBase << ','
            << nUnrefinePoints << ',' << amrLayerDiagnosticsLower_ << ','
            << amrLayerDiagnosticsUpper_ << '\n';
    }

    resetAmrLayerDiagnosticsState();
}

void fastDynamicFvMesh::clearRefinementGradIndicatorCache() const
{
    refinementGradIndicatorCellCache_.setSize(0);
    refinementGradIndicatorPointCache_.setSize(0);
    refinementGradIndicatorCacheTimeIndex_ = -1;
    refinementGradIndicatorCacheNCells_ = -1;
    refinementGradIndicatorCacheNPoints_ = -1;
    refinementGradIndicatorCacheSourceField_ = word::null;
    refinementGradIndicatorCacheMagnitudeField_ = word::null;
    refinementGradIndicatorCacheValid_ = false;
}

void fastDynamicFvMesh::ensureRefinementGradIndicatorCache() const
{
    if (!refinementUseGradIndicator_)
    {
        clearRefinementGradIndicatorCache();
        return;
    }

    const label timeIndex = this->time().timeIndex();
    const label nCells = this->nCells();
    const label nPoints = this->nPoints();

    if
    (
        refinementGradIndicatorCacheValid_
     && refinementGradIndicatorCacheTimeIndex_ == timeIndex
     && refinementGradIndicatorCacheNCells_ == nCells
     && refinementGradIndicatorCacheNPoints_ == nPoints
     && refinementGradIndicatorCacheSourceField_ == refinementGradIndicatorField_
     && refinementGradIndicatorCacheMagnitudeField_
        == refinementGradIndicatorMagnitudeField_
     && refinementGradIndicatorCellCache_.size() == nCells
     && refinementGradIndicatorPointCache_.size() == nPoints
    )
    {
        return;
    }

    auto cacheIndicator = [&](const volScalarField& indicator)
    {
        if (indicator.size() != nCells)
        {
            FatalErrorInFunction
                << "Runtime AMR gradient indicator field '"
                << indicator.name() << "' has " << indicator.size()
                << " cells, but the mesh has " << nCells << " cells."
                << exit(FatalError);
        }

        refinementGradIndicatorCellCache_ = indicator.primitiveField();
        refinementGradIndicatorPointCache_ = maxCellField(indicator);

        if (refinementGradIndicatorPointCache_.size() != nPoints)
        {
            FatalErrorInFunction
                << "Runtime AMR gradient point indicator built from field '"
                << indicator.name() << "' has "
                << refinementGradIndicatorPointCache_.size()
                << " points, but the mesh has " << nPoints << " points."
                << exit(FatalError);
        }

        refinementGradIndicatorCacheTimeIndex_ = timeIndex;
        refinementGradIndicatorCacheNCells_ = nCells;
        refinementGradIndicatorCacheNPoints_ = nPoints;
        refinementGradIndicatorCacheSourceField_ =
            refinementGradIndicatorField_;
        refinementGradIndicatorCacheMagnitudeField_ =
            refinementGradIndicatorMagnitudeField_;
        refinementGradIndicatorCacheValid_ = true;
    };

    if (!refinementGradIndicatorMagnitudeField_.empty())
    {
        if
        (
            !this->foundObject<volScalarField>
            (
                refinementGradIndicatorMagnitudeField_
            )
        )
        {
            FatalErrorInFunction
                << "Runtime AMR gradient indicator magnitude field '"
                << refinementGradIndicatorMagnitudeField_
                << "' is configured, but it is not available in the "
                << "objectRegistry." << exit(FatalError);
        }

        cacheIndicator
        (
            this->lookupObject<volScalarField>
            (
                refinementGradIndicatorMagnitudeField_
            )
        );
        return;
    }

    if (!this->foundObject<volScalarField>(refinementGradIndicatorField_))
    {
        FatalErrorInFunction
            << "Runtime AMR gradient indicator is enabled, but source field '"
            << refinementGradIndicatorField_
            << "' is not available in the objectRegistry." << exit(FatalError);
    }

    const volScalarField& sourceField =
        this->lookupObject<volScalarField>(refinementGradIndicatorField_);
    const tmp<volScalarField> tIndicator = mag(fvc::grad(sourceField));
    cacheIndicator(tIndicator());
}

scalarField fastDynamicFvMesh::refinementIndicatorCellField(
    const scalarField& fallbackCellField) const
{
    if (!refinementUseGradIndicator_)
    {
        return fallbackCellField;
    }

    ensureRefinementGradIndicatorCache();
    return refinementGradIndicatorCellCache_;
}

scalarField fastDynamicFvMesh::refinementIndicatorPointField(
    const scalarField& fallbackPointField) const
{
    if (!refinementUseGradIndicator_)
    {
        return fallbackPointField;
    }

    ensureRefinementGradIndicatorCache();
    return refinementGradIndicatorPointCache_;
}

void fastDynamicFvMesh::selectRefineCandidates(const scalar lowerRefineLevel,
    const scalar upperRefineLevel, const scalarField& vFld,
    bitSet& candidateCell) const
{
    const scalarField indicatorCell = refinementIndicatorCellField(vFld);

    dynamicRefineFvMesh::selectRefineCandidates(
        lowerRefineLevel, upperRefineLevel, indicatorCell, candidateCell);

    if (amrLayerDiagnosticsEnabled_)
    {
        amrLayerDiagnosticsCandidateCell_ = candidateCell;
        amrLayerDiagnosticsIndicatorCell_ = indicatorCell;
        amrLayerDiagnosticsLower_ = lowerRefineLevel;
        amrLayerDiagnosticsUpper_ = upperRefineLevel;
    }

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

labelList fastDynamicFvMesh::selectRefineCells(const label maxCells,
    const label maxRefinement, const bitSet& candidateCell) const
{
    if (refinementMinCellVolume_ <= 0 && refinementMinEdgeLength_ <= 0)
    {
        const labelList selectedCells = dynamicRefineFvMesh::selectRefineCells(
            maxCells, maxRefinement, candidateCell);
        captureAmrLayerDiagnostics(amrLayerDiagnosticsLower_,
            amrLayerDiagnosticsUpper_, maxCells, maxRefinement, candidateCell,
            selectedCells);
        return selectedCells;
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
        Info << "Refinement size-floor filter rejected " << totalRejected
             << " candidate cells (refinementMinCellVolume="
             << refinementMinCellVolume_
             << ", refinementMinEdgeLength=" << refinementMinEdgeLength_ << ")."
             << endl;
    }
    else if (Pstream::master() && totalCandidates > 0)
    {
        Info << "Refinement size-floor filter accepted all " << totalCandidates
             << " candidate cells." << endl;
    }

    const labelList selectedCells = dynamicRefineFvMesh::selectRefineCells(
        maxCells, maxRefinement, sizeFloorFiltered);
    captureAmrLayerDiagnostics(amrLayerDiagnosticsLower_,
        amrLayerDiagnosticsUpper_, maxCells, maxRefinement, candidateCell,
        selectedCells);
    return selectedCells;
}

labelList fastDynamicFvMesh::selectUnrefinePoints(const scalar unrefineLevel,
    const bitSet& markedCell, const scalarField& pFld) const
{
    const scalarField indicatorPoint = refinementIndicatorPointField(pFld);

    const labelList selectedPoints = dynamicRefineFvMesh::selectUnrefinePoints(
        unrefineLevel, markedCell, indicatorPoint);
    if (amrLayerDiagnosticsEnabled_)
    {
        amrLayerDiagnosticsUnrefinePoints_ = selectedPoints.size();
    }
    return selectedPoints;
}

void fastDynamicFvMesh::updateMesh(const mapPolyMesh& mpm)
{
    dynamicFvMesh::updateMesh(mpm);
    clearRefinementGradIndicatorCache();

    if (!meshRefinementSupport_)
    {
        return;
    }

    ensureModeShapesForCurrentMesh(mpm);
    buildForceProjectionCaches(&mpm);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fastDynamicFvMesh::updateRefinementTopology()
{
    if (runtimeRefinementEnabled_)
    {
        if (!trackAmrTimingEnabled_)
        {
            dynamicRefineFvMesh::updateTopology();
            writeAmrLayerDiagnosticsSummary();
        }
        else
        {
            const label nCellsBefore =
                returnReduce(this->nCells(), sumOp<label>());
            const label nPointsBefore =
                returnReduce(this->nPoints(), sumOp<label>());
            const scalar amrStartCpuTime = this->time().elapsedCpuTime();
            const scalar amrStartClockTime = this->time().elapsedClockTime();

            dynamicRefineFvMesh::updateTopology();
            writeAmrLayerDiagnosticsSummary();

            const scalar localCpu =
                this->time().elapsedCpuTime() - amrStartCpuTime;
            const scalar localClock =
                this->time().elapsedClockTime() - amrStartClockTime;
            const scalar maxCpu = returnReduce(localCpu, maxOp<scalar>());
            const scalar maxClock = returnReduce(localClock, maxOp<scalar>());
            const label nCellsAfter =
                returnReduce(this->nCells(), sumOp<label>());
            const label nPointsAfter =
                returnReduce(this->nPoints(), sumOp<label>());

            if (Pstream::master())
            {
                Info << "Runtime AMR timing at t=" << this->time().timeName()
                     << " [max over ranks]: CPU=" << maxCpu
                     << " s, Clock=" << maxClock << " s, cells "
                     << nCellsBefore << " -> " << nCellsAfter << ", points "
                     << nPointsBefore << " -> " << nPointsAfter << endl;
            }
        }
    }
}

bool fastDynamicFvMesh::refinementTopologyActive() const
{
    return runtimeRefinementEnabled_;
}

const labelList& fastDynamicFvMesh::amrPointLevel() const
{
    return meshCutter_.pointLevel();
}

void fastDynamicFvMesh::readAmrControls(
    const dictionary& fdmDict,
    const dictionary& dynamicMeshDict)
{
    fdmDict.readIfPresent("trackAmrTiming", trackAmrTimingEnabled_);
    fdmDict.readIfPresent(
        "refinementMappingDiagnostics", refinementMappingDiagnostics_);
    scalar refinementFaceMapTolerance = 0.0;
    const bool foundRefinementFaceMapTolerance = fdmDict.readIfPresent(
        "refinementFaceMapTolerance", refinementFaceMapTolerance);
    fdmDict.readIfPresent("refinementMinCellVolume", refinementMinCellVolume_);
    fdmDict.readIfPresent("refinementMinEdgeLength", refinementMinEdgeLength_);
    fdmDict.readIfPresent("amrLayerDiagnostics", amrLayerDiagnosticsEnabled_);
    fdmDict.readIfPresent("amrLayerDiagnosticsZ", amrLayerDiagnosticsZ_);
    fdmDict.readIfPresent(
        "amrLayerDiagnosticsThickness", amrLayerDiagnosticsThickness_);
    fdmDict.readIfPresent(
        "amrLayerDiagnosticsNeighborDepth", amrLayerDiagnosticsNeighborDepth_);
    fdmDict.readIfPresent(
        "amrLayerDiagnosticsMaxRows", amrLayerDiagnosticsMaxRows_);
    fdmDict.readIfPresent(
        "amrLayerDiagnosticsFilePrefix", amrLayerDiagnosticsFilePrefix_);

    if (dynamicMeshDict.found("dynamicRefineFvMeshCoeffs"))
    {
        const dictionary& refineDict =
            dynamicMeshDict.subDict("dynamicRefineFvMeshCoeffs");
        refineDict.readIfPresent(
            "useGradIndicator", refinementUseGradIndicator_);
        refineDict.readIfPresent(
            "gradIndicatorField", refinementGradIndicatorField_);
        refineDict.readIfPresent("gradIndicatorMagnitudeField",
            refinementGradIndicatorMagnitudeField_);
    }

    if (runtimeRefinementEnabled_ && !meshRefinementSupport_ && Pstream::master())
    {
        WarningInFunction
            << "dynamicRefineFvMeshCoeffs is present, but "
            << "'meshRefinementSupport' is false in " << fdmTypeName()
            << "Coeffs. Runtime topology changes can invalidate modal "
            << "shape and FSI face caches; set meshRefinementSupport true "
            << "for AMR-coupled fastDynamicFvMesh cases." << endl;
    }

    if (refinementFaceMapTolerance < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'refinementFaceMapTolerance' must be >= 0 in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.name() << exit(FatalIOError);
    }
    else if (foundRefinementFaceMapTolerance && Pstream::master())
    {
        WarningInFunction
            << "Entry 'refinementFaceMapTolerance' is deprecated and ignored. "
            << "Topology mapping is now fully driven by mapPolyMesh lineage."
            << endl;
    }

    if (refinementMinCellVolume_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'refinementMinCellVolume' must be >= 0 in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.name() << exit(FatalIOError);
    }

    if (refinementMinEdgeLength_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'refinementMinEdgeLength' must be >= 0 in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.name() << exit(FatalIOError);
    }

    if (amrLayerDiagnosticsThickness_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'amrLayerDiagnosticsThickness' must be >= 0 in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.name() << exit(FatalIOError);
    }

    if (amrLayerDiagnosticsNeighborDepth_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'amrLayerDiagnosticsNeighborDepth' must be >= 0 in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.name() << exit(FatalIOError);
    }

    if (amrLayerDiagnosticsMaxRows_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'amrLayerDiagnosticsMaxRows' must be >= 0 in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.name() << exit(FatalIOError);
    }

    if (amrLayerDiagnosticsEnabled_ && amrLayerDiagnosticsFilePrefix_.empty())
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'amrLayerDiagnosticsFilePrefix' cannot be empty when "
            << "'amrLayerDiagnostics' is enabled in sub-dictionary '"
            << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.name()
            << exit(FatalIOError);
    }

    if (refinementUseGradIndicator_ && refinementGradIndicatorField_.empty())
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'gradIndicatorField' cannot be empty when "
            << "'useGradIndicator' is enabled in dynamicRefineFvMeshCoeffs of "
            << dynamicMeshDict.name() << exit(FatalIOError);
    }

    if (Pstream::master() && runtimeRefinementEnabled_ &&
        refinementUseGradIndicator_)
    {
        if (refinementGradIndicatorMagnitudeField_.empty())
        {
            Info << "Runtime AMR indicator uses cached |grad("
                 << refinementGradIndicatorField_
                 << ")| with dynamicRefineFvMesh thresholds." << endl;
        }
        else
        {
            Info << "Runtime AMR indicator uses registered field '"
                 << refinementGradIndicatorMagnitudeField_
                 << "' as |grad(" << refinementGradIndicatorField_
                 << ")| with dynamicRefineFvMesh thresholds." << endl;
        }
    }

    if (Pstream::master() && amrLayerDiagnosticsEnabled_)
    {
        Info << "Runtime AMR bottom-layer diagnostics enabled: z="
             << amrLayerDiagnosticsZ_ << ", thickness="
             << amrLayerDiagnosticsThickness_ << ", neighborDepth="
             << amrLayerDiagnosticsNeighborDepth_ << ", maxRows="
             << amrLayerDiagnosticsMaxRows_ << ", filePrefix='"
             << amrLayerDiagnosticsFilePrefix_ << "'." << endl;
    }
}

void fastDynamicFvMesh::logModeShapeMappingSummary(const label oldNPoints,
    const label newNPoints, const bool anyRefinement,
    const bool shrinkingOrEqualTopology, const label localDirectMapped,
    const label localDeferredAddedPoints, const label localLineageMapped,
    const label localEdgeInterpolated, const label localFaceInterpolated,
    const label localCellInterpolated, const label localUnknown) const
{
    if (!refinementMappingDiagnostics_)
    {
        return;
    }

    const label globalOldNPoints = returnReduce(oldNPoints, sumOp<label>());
    const label globalNewNPoints = returnReduce(newNPoints, sumOp<label>());
    const label globalDirectMapped =
        returnReduce(localDirectMapped, sumOp<label>());
    const label globalDeferredAddedPoints =
        returnReduce(localDeferredAddedPoints, sumOp<label>());
    const label globalLineageMapped =
        returnReduce(localLineageMapped, sumOp<label>());
    const label globalEdgeInterpolated =
        returnReduce(localEdgeInterpolated, sumOp<label>());
    const label globalFaceInterpolated =
        returnReduce(localFaceInterpolated, sumOp<label>());
    const label globalCellInterpolated =
        returnReduce(localCellInterpolated, sumOp<label>());
    const label globalUnknown = returnReduce(localUnknown, sumOp<label>());

    if (Pstream::master())
    {
        Info << "Mode-shape topology mapping diagnostics:" << nl
             << "  topologyChange="
             << (anyRefinement ? "refinement/mixed" : "unrefine/shrink")
             << ", localTopology="
             << (shrinkingOrEqualTopology ? "shrinkingOrEqual" : "growing")
             << nl << "  points old/new=" << globalOldNPoints << "/"
             << globalNewNPoints << nl
             << "  directPointMap=" << globalDirectMapped
             << ", deferredPointMapMasters=" << globalDeferredAddedPoints
             << ", pointsFromPoints=" << globalLineageMapped << nl
             << "  interpolated edge/face/cell=" << globalEdgeInterpolated
             << "/" << globalFaceInterpolated << "/"
             << globalCellInterpolated << nl
             << "  finalUnknown=" << globalUnknown << endl;
    }
}

void fastDynamicFvMesh::ensureModeShapesForCurrentMesh(const mapPolyMesh& mpm)
{
    if (!meshRefinementSupport_ || nMode_ <= 0 || modeShapes_.size() == 0)
    {
        return;
    }

    const label newNPoints = fdmMesh_.points().size();
    const label oldNPoints = modeShapes_[0].size();
    const bool shrinkingOrEqualTopology = (newNPoints <= oldNPoints);
    const bool localPointCountChanged = (newNPoints != oldNPoints);
    const bool anyPointCountChanged =
        returnReduce(label(localPointCountChanged), sumOp<label>()) > 0;

    if (!anyPointCountChanged)
    {
        return;
    }

    const bool anyRefinement =
        returnReduce(label(newNPoints > oldNPoints), sumOp<label>()) > 0;

    const labelList& pointMap = mpm.pointMap();
    const labelList& reversePointMap = mpm.reversePointMap();

    if (pointMap.size() != newNPoints)
    {
        FatalErrorInFunction << "pointMap size (" << pointMap.size()
                             << ") does not match new mesh point count ("
                             << newNPoints << ")." << exit(FatalError);
    }

    List<vectorField> oldModeShapes(modeShapes_);
    modeShapes_.setSize(nMode_);

    forAll(modeShapes_, modeI)
    {
        modeShapes_[modeI].setSize(newNPoints, vector::zero);
    }

    boolList knownMask(newNPoints, false);
    label localDirectMapped = 0;
    label localDeferredAddedPoints = 0;

    for (label pointI = 0; pointI < newNPoints; ++pointI)
    {
        const label oldPoint = pointMap[pointI];
        if (oldPoint < 0 || oldPoint >= oldNPoints)
        {
            continue;
        }

        if (anyRefinement)
        {
            if (oldPoint >= reversePointMap.size())
            {
                FatalErrorInFunction
                    << "reversePointMap size (" << reversePointMap.size()
                    << ") cannot address old point " << oldPoint << "."
                    << exit(FatalError);
            }

            // During refinement, pointMap can assign newly-created points to
            // a single master old point. Those entries must not be treated as
            // fully mapped values; defer them to topology interpolation.
            if (reversePointMap[oldPoint] != pointI)
            {
                ++localDeferredAddedPoints;
                continue;
            }
        }

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            modeShapes_[modeI][pointI] = oldModeShapes[modeI][oldPoint];
        }

        knownMask[pointI] = true;
        ++localDirectMapped;
    }

    const label globalDeferredAddedPoints =
        returnReduce(localDeferredAddedPoints, sumOp<label>());
    if (Pstream::master() && globalDeferredAddedPoints > 0)
    {
        Info << "Deferred " << globalDeferredAddedPoints
             << " refinement-added points from pointMap master labels to "
             << "topology interpolation." << endl;
    }

    // Prefer direct topology lineage for refinement-added points.
    // mapPolyMesh stores source old points for new points in
    // pointsFromPointsMap. Averaging source modal values preserves geometric
    // continuity better than a generic nearest-corner search.
    label localLineageMapped = 0;
    const List<objectMap>& pointsFromPoints = mpm.pointsFromPointsMap();

    forAll(pointsFromPoints, mapI)
    {
        const objectMap& pointMapEntry = pointsFromPoints[mapI];
        const label newPoint = pointMapEntry.index();

        if (newPoint < 0 || newPoint >= newNPoints || knownMask[newPoint])
        {
            continue;
        }

        const labelList& masterPoints = pointMapEntry.masterObjects();
        if (masterPoints.size() == 0)
        {
            continue;
        }

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            modeShapes_[modeI][newPoint] = vector::zero;
        }

        label validMasterCount = 0;
        forAll(masterPoints, masterI)
        {
            const label oldPoint = masterPoints[masterI];
            if (oldPoint < 0 || oldPoint >= oldNPoints)
            {
                continue;
            }

            ++validMasterCount;
            for (label modeI = 0; modeI < nMode_; ++modeI)
            {
                modeShapes_[modeI][newPoint] += oldModeShapes[modeI][oldPoint];
            }
        }

        if (validMasterCount == 0)
        {
            continue;
        }

        const scalar invMasterCount = 1.0 / scalar(validMasterCount);
        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            modeShapes_[modeI][newPoint] *= invMasterCount;
        }

        knownMask[newPoint] = true;
        ++localLineageMapped;
    }

    const label globalLineageMapped =
        returnReduce(localLineageMapped, sumOp<label>());
    if (Pstream::master() && globalLineageMapped > 0)
    {
        Info << "Topology lineage mapped " << globalLineageMapped
             << " refinement-added points from mapPolyMesh::pointsFromPoints."
             << endl;
    }

    syncMappedModeShapes(knownMask);

    label localKnownAfterSync = 0;
    forAll(knownMask, pointI)
    {
        if (knownMask[pointI])
        {
            ++localKnownAfterSync;
        }
    }

    if (localKnownAfterSync == 0 && newNPoints > 0)
    {
        FatalErrorInFunction
            << "No mapped old points found on processor " << Pstream::myProcNo()
            << " while processing topology change; "
            << "cannot interpolate modal shapes on refined mesh."
            << exit(FatalError);
    }

    if (anyRefinement)
    {
        label globalKnownPrev =
            returnReduce(localKnownAfterSync, sumOp<label>());
        label localEdgeInterpolated = 0;
        label localFaceInterpolated = 0;
        label localCellInterpolated = 0;

        for (label passI = 0; passI < 4; ++passI)
        {
            const boolList sourceKnownMask(knownMask);

            if (!shrinkingOrEqualTopology)
            {
                localEdgeInterpolated += interpolateModeShapesFromKnown(
                    knownMask, sourceKnownMask, 2, true, "edge midpoint");
                localFaceInterpolated += interpolateModeShapesFromKnown(
                    knownMask, sourceKnownMask, 4, false, "face midpoint");
                localCellInterpolated += interpolateModeShapesFromKnown(
                    knownMask, sourceKnownMask, 8, false, "cell midpoint");
            }

            syncMappedModeShapes(knownMask);

            label localKnownNow = 0;
            forAll(knownMask, pointI)
            {
                if (knownMask[pointI])
                {
                    ++localKnownNow;
                }
            }

            const label globalKnownNow =
                returnReduce(localKnownNow, sumOp<label>());

            if (globalKnownNow <= globalKnownPrev)
            {
                break;
            }

            globalKnownPrev = globalKnownNow;
        }

        syncMappedModeShapes(knownMask);

        label localUnknown = 0;
        forAll(knownMask, pointI)
        {
            if (!knownMask[pointI])
            {
                ++localUnknown;
            }
        }

        logModeShapeMappingSummary(oldNPoints, newNPoints, anyRefinement,
            shrinkingOrEqualTopology, localDirectMapped,
            localDeferredAddedPoints, localLineageMapped,
            localEdgeInterpolated, localFaceInterpolated,
            localCellInterpolated, localUnknown);

        if (localUnknown > 0)
        {
            FatalErrorInFunction
                << "Refinement interpolation failed for " << localUnknown
                << " local points on processor " << Pstream::myProcNo()
                << " points."
                << " The 2/4/8 topology-corner rules could not classify all "
                << "added points. Adjust refinement controls or "
                << "refinementInterpTolerance." << exit(FatalError);
        }

        return;
    }

    syncMappedModeShapes(knownMask);

    label localUnknown = 0;
    forAll(knownMask, pointI)
    {
        if (!knownMask[pointI])
        {
            ++localUnknown;
        }
    }

    logModeShapeMappingSummary(oldNPoints, newNPoints, anyRefinement,
        shrinkingOrEqualTopology, localDirectMapped,
        localDeferredAddedPoints, localLineageMapped, 0, 0, 0, localUnknown);

    if (localUnknown > 0)
    {
        FatalErrorInFunction
            << (shrinkingOrEqualTopology
                       ? "Topology mapping failed for "
                       : "Refinement interpolation failed for ")
            << localUnknown << " local points on processor "
            << Pstream::myProcNo() << " points."
            << (shrinkingOrEqualTopology
                       ? " Some surviving points were not mapped by pointMap."
                       : " The 2/4/8 topology-corner rules could not classify "
                         "all "
                         "added points. Adjust refinement controls or "
                         "refinementInterpTolerance.")
            << exit(FatalError);
    }
}


} // namespace Foam
