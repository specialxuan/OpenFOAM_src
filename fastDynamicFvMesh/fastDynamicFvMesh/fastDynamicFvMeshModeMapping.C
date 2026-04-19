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
// readModeShapes()
// Purpose: load modal coordinates and modal displacement shapes from the case
// 'mode' directory and map them to mesh points. High-level steps (master
// process):
//  1) Locate mode directory (time()/mode or parent time path when running in
//  parallel). 2) Open FluidNodeCoor.csv and robustly parse the header (accepts
//  comma or space separators). Header contains at least: dummy, nNode, nMode.
//  3) Allocate containers csvPoints and csvShapes and read nCsvNodes
//  coordinates. 4) For each mode m, open FluidNodeDisp(m).csv: read frequency
//  (and optional node/mode counts), then read per-node dx,dy,dz displacements.
//  Allow skip of a non-numeric leading line if needed. 5) Broadcast sizes and
//  data to all processors, then allocate local modeShapes_ sized to localPoints
//  and initialize to zero. 6) Map CSV nodes to local mesh points using a
//  brute-force nearest-neighbor search within mappingTolerance_. 7) Report
//  mapped counts and error out if zero mapped (user must check coordinates or
//  mappingTolerance).
// Notes:
//  - The implementation favors correctness and simplicity (brute-force search)
//  over asymptotic performance; mapping is done once at startup.
//  - Mode frequency values read from files are stored in modeFreq_ and expected
//  to be in Hz (omega computed later as 2*pi*freq).
void fastDynamicFvMesh::readModeShapes()
{
    // 1. Read files on master
    List<point> csvPoints;
    List<List<vector>> csvShapes;
    label nCsvNodes = 0;
    fileName modeDir = this->time().path() / "mode";

    if (Pstream::parRun() && !exists(modeDir))
    {
        modeDir = this->time().path().path() / "mode";
    }

    readModeFiles(modeDir, csvPoints, csvShapes, nCsvNodes);

    // Broadcast sizes
    Pstream::broadcast(nMode_);
    Pstream::broadcast(nCsvNodes);

    // Resize local arrays
    if (!Pstream::master())
        modeFreq_.setSize(nMode_); // Only slaves need resize now

    if (modeMass_.size() == 0)
    {
        modeMass_.setSize(nMode_, 1.0);
    }
    else if (modeMass_.size() != nMode_)
    {
        FatalErrorInFunction << "modalMass size (" << modeMass_.size()
                             << ") does not match number of modes (" << nMode_
                             << ")" << exit(FatalError);
    }

    if (modeDamp_.size() == 0)
    {
        modeDamp_.setSize(nMode_, 0.0);
    }
    else if (modeDamp_.size() != nMode_)
    {
        FatalErrorInFunction << "modalDamp size (" << modeDamp_.size()
                             << ") does not match number of modes (" << nMode_
                             << ")" << exit(FatalError);
    }

    if (nMode_ <= 0)
    {
        FatalErrorInFunction
            << "No fluid modes were loaded from " << modeDir << nl
            << "Ensure mode/FluidNodeCoor.csv and mode/FluidNodeDisp*.csv are "
            << "present and valid." << exit(FatalError);
    }

    modeForce_.setSize(nMode_, 0.0);
    modePressureForce_.setSize(nMode_, 0.0);
    modeShearForce_.setSize(nMode_, 0.0);
    modeForce0_.setSize(nMode_, 0.0);
    modeState_.setSize(nMode_, vector::zero);
    modeState0_.setSize(nMode_, vector::zero);
    initVelocity_.setSize(nMode_, 0.0);
    appliedModeDisp_.setSize(nMode_, 0.0);

    readLegacyParameters(modeDir);

    // Broadcast data
    Pstream::broadcast(modeFreq_);
    Pstream::broadcast(csvPoints);
    Pstream::broadcast(initVelocity_);

    // Resize csvShapes on all procs to ensure loop consistency
    csvShapes.setSize(nMode_);

    forAll(csvShapes, m) Pstream::broadcast(csvShapes[m]);

    // Map to local mesh points.
    // Default to current mesh points. This is required on latestTime restarts
    // after runtime AMR because current point count can differ from
    // constant/polyMesh/points.
    pointField mappingPoints(this->points());
    pointField basePoints;

    // Try to read points from constant/polyMesh (or
    // processor*/constant/polyMesh) Note: in parallel, the polyMesh directory
    // is inside the processor directory. IOobject will handle processor path
    // automatically if we use NO_WRITE.

    // We construct an IOobject to look for "points" in the constant instance.
    IOobject ioPoints("points", this->time().constant(), polyMesh::meshSubDir,
        *this, IOobject::MUST_READ, IOobject::NO_WRITE,
        false // register
    );

    // Custom wrapper to force expecting "vectorField" class name
    // regardless of what Field<vector>::typeName is defined as.
    class vectorFieldIO : public IOField<vector>
    {
    public:
        // Force the type name to match the file header
        virtual const word& type() const
        {
            static const word typeName_val("vectorField");
            return typeName_val;
        }

        vectorFieldIO(const IOobject& io) : IOField<vector>(io) {}
    };

    // Attempt to read with the custom wrapper.
    bool readOk = false;
    try
    {
        // Don't check header via typeHeaderOk as it might check the wrong type
        // Just verify it exists (MUST_READ handles this usually, but we want to
        // catch it)
        if (ioPoints.headerClassName() == "vectorField")
        {
            vectorFieldIO pIO(ioPoints);
            basePoints = pIO;
            readOk = true;
        }
        else if (ioPoints.typeHeaderOk<vectorFieldIO>(true))
        {
            vectorFieldIO pIO(ioPoints);
            basePoints = pIO;
            readOk = true;
        }
    }
    catch (...)
    {
    }

    if (readOk)
    {
        if (basePoints.size() == mappingPoints.size())
        {
            mappingPoints = basePoints;

            if (Pstream::master())
            {
                Info << "  Read " << basePoints.size()
                     << " base points for mapping from "
                     << ioPoints.objectPath() << endl;
            }
        }
        else if (Pstream::master())
        {
            WarningInFunction
                << "constant/polyMesh points count (" << basePoints.size()
                << ") differs from current mesh points count ("
                << mappingPoints.size() << ") on startup/restart. "
                << "Using current mesh points for modal mapping." << endl;
        }
    }
    else if (Pstream::master())
    {
        WarningInFunction
            << "Could not read base points from constant/polyMesh "
            << "(expected class vectorField). Using current mesh points "
            << "for modal mapping." << endl;
    }

    modeShapes_.setSize(nMode_);
    forAll(modeShapes_, m)
    {
        modeShapes_[m].setSize(mappingPoints.size());
        modeShapes_[m] = vector::zero;
    }

    const labelList* pointLevelPtr = nullptr;
    bool topologyOnlyForRefinedPoints = false;
    if (meshRefinementSupport_)
    {
        const labelList& pointLevel = meshCutter_.pointLevel();
        if (pointLevel.size() != mappingPoints.size())
        {
            FatalErrorInFunction
                << "meshCutter pointLevel size (" << pointLevel.size()
                << ") does not match mesh point count (" << mappingPoints.size()
                << ")." << exit(FatalError);
        }

        label localRefinedPointCount = 0;
        forAll(pointLevel, pointI)
        {
            if (pointLevel[pointI] > 0)
            {
                ++localRefinedPointCount;
            }
        }

        const label globalRefinedPointCount =
            returnReduce(localRefinedPointCount, sumOp<label>());
        topologyOnlyForRefinedPoints = (globalRefinedPointCount > 0);
        pointLevelPtr = &pointLevel;

        if (Pstream::master() && topologyOnlyForRefinedPoints)
        {
            Info
                << "Refined startup mesh detected (" << globalRefinedPointCount
                << " points with pointLevel>0). "
                << "Geometric CSV mapping is restricted to level-0 points; "
                << "refined points must be recovered by topology interpolation."
                << endl;
        }
    }

    // Brute force search (simplest correct implementation)
    // Only perform mapping if we have CSV nodes
    label mappedCount = 0;
    boolList mappedMask(mappingPoints.size(), false);
    if (nCsvNodes > 0)
    {
        scalar tolSqr = sqr(mappingTolerance_);

        forAll(mappingPoints, pI)
        {
            if (topologyOnlyForRefinedPoints && (*pointLevelPtr)[pI] > 0)
            {
                continue;
            }

            const point& p = mappingPoints[pI];
            scalar minDistSqr = GREAT;
            label nearestIdx = -1;

            forAll(csvPoints, cI)
            {
                scalar d = magSqr(p - csvPoints[cI]);
                if (d < minDistSqr)
                {
                    minDistSqr = d;
                    nearestIdx = cI;
                }
            }

            if (minDistSqr < tolSqr && nearestIdx != -1)
            {
                mappedCount++;
                mappedMask[pI] = true;
                for (label m = 0; m < nMode_; ++m)
                {
                    modeShapes_[m][pI] = csvShapes[m][nearestIdx];
                }
            }
        }
    }

    if (meshRefinementSupport_)
    {
        syncMappedModeShapes(mappedMask);
    }

    label totalMapped = returnReduce(mappedCount, sumOp<label>());
    const label totalPoints =
        returnReduce(label(mappingPoints.size()), sumOp<label>());

    if (meshRefinementSupport_ && totalMapped > 0 && totalMapped < totalPoints)
    {
        boolList knownMask(mappedMask);

        label prevGlobalKnown = totalMapped;
        for (label passI = 0; passI < 4; ++passI)
        {
            const boolList sourceKnownMask(knownMask);

            interpolateModeShapesFromKnown(
                knownMask, sourceKnownMask, 2, true, "edge midpoint");
            interpolateModeShapesFromKnown(
                knownMask, sourceKnownMask, 4, false, "face midpoint");
            interpolateModeShapesFromKnown(
                knownMask, sourceKnownMask, 8, false, "cell midpoint");

            syncMappedModeShapes(knownMask);

            label localKnownPass = 0;
            forAll(knownMask, pointI)
            {
                if (knownMask[pointI])
                {
                    ++localKnownPass;
                }
            }

            const label totalKnownPass =
                returnReduce(localKnownPass, sumOp<label>());

            if (totalKnownPass <= prevGlobalKnown)
            {
                break;
            }

            prevGlobalKnown = totalKnownPass;
        }

        syncMappedModeShapes(knownMask);

        label localKnown = 0;
        forAll(knownMask, pointI)
        {
            if (knownMask[pointI])
            {
                ++localKnown;
            }
        }

        const label totalKnown = returnReduce(localKnown, sumOp<label>());

        if (Pstream::master())
        {
            Info << "Refinement-aware startup interpolation raised mapped "
                    "points "
                 << "from " << totalMapped << " to " << totalKnown << " out of "
                 << totalPoints << "." << endl;
        }

        if (totalKnown != totalPoints && topologyOnlyForRefinedPoints)
        {
            FatalErrorInFunction
                << "With meshRefinementSupport enabled on a refined startup "
                   "mesh, "
                << (totalPoints - totalKnown)
                << " points remain unmapped after 2/4/8 topology interpolation."
                << nl
                << "No geometric fallback is allowed for AMR-generated points."
                << nl << "Check refinement lineage consistency and "
                << "refinementInterpTolerance." << exit(FatalError);
        }
        else if (totalKnown != totalPoints && Pstream::master())
        {
            WarningInFunction << "With meshRefinementSupport enabled, "
                              << (totalPoints - totalKnown)
                              << " points remain unmapped after 2/4/8 topology "
                                 "interpolation. "
                              << "These points keep zero modal displacement "
                                 "until they can be "
                              << "mapped by later topology updates." << endl;
        }

        totalMapped = totalKnown;
    }

    if (Pstream::master())
    {
        Info << "Mapped " << totalMapped << " points out of " << totalPoints
             << " mesh points." << endl;
    }

    if (totalMapped == 0)
    {
        FatalErrorInFunction
            << "Mapped 0 mesh points from " << modeDir << nl
            << "Check that the mode CSV coordinates correspond to the current "
            << "mesh and that mappingTolerance is large enough."
            << exit(FatalError);
    }

    if (Pstream::master() && totalMapped != totalPoints)
    {
        WarningInFunction
            << "Mapped " << totalMapped << " of " << totalPoints
            << " mesh points. Unmapped points will remain stationary." << endl;
    }

    buildForceProjectionCaches();

    if (structuralForceEnabled_)
    {
        pointField structPoints;
        List<vectorField> structShapes;
        label nStructNodes = 0;
        label nStructModes = 0;

        readStructuralModeFiles(
            modeDir, structPoints, structShapes, nStructNodes, nStructModes);

        Pstream::broadcast(nStructNodes);
        Pstream::broadcast(nStructModes);
        Pstream::broadcast(structPoints);

        if (!Pstream::master())
        {
            structShapes.setSize(nStructModes);
        }

        forAll(structShapes, m)
        {
            Pstream::broadcast(structShapes[m]);
        }

        structNodeCoords_.transfer(structPoints);
        structModeShapes_.transfer(structShapes);
        structuralModeForce_.setSize(nMode_, 0.0);
        readStructuralForceFiles(modeDir);

        if (Pstream::master())
        {
            label totalNodes = 0;
            forAll(structuralForceNodeIDs_, fi)
            {
                totalNodes += structuralForceNodeIDs_[fi].size();
            }

            Info << "Loaded " << nStructuralForces_
                 << " structural force files with mapped node references="
                 << totalNodes << endl;
        }
    }
    else
    {
        structuralModeForce_.setSize(nMode_, 0.0);
        structuralForceNodeIDs_.clear();
        structuralForceTimes_.clear();
        structuralForceValues_.clear();
        structNodeCoords_.clear();
        structModeShapes_.clear();
    }
}

void fastDynamicFvMesh::buildForceProjectionCaches(const mapPolyMesh* mpmPtr)
{
    const polyBoundaryMesh& pbm = this->boundaryMesh();

    const List<List<labelList>> oldFaceToReferenceFaces(
        fsiFaceToReferenceFaces_);
    const List<List<scalarField>> oldFaceToReferenceWeights(
        fsiFaceToReferenceWeights_);
    const bool hasOldReferenceMapping =
        oldFaceToReferenceFaces.size() == fsiPatches_.size() &&
        oldFaceToReferenceWeights.size() == fsiPatches_.size();

    fsiPatchIDs_.setSize(fsiPatches_.size(), -1);
    fsiPolyPatches_.setSize(fsiPatches_.size(), nullptr);
    faceModeProjection_.setSize(fsiPatches_.size());
    pressureScaleCache_.setSize(fsiPatches_.size());
    pressureScaleInitialized_.setSize(fsiPatches_.size(), false);
    fsiFaceToReferenceFaces_.setSize(fsiPatches_.size());
    fsiFaceToReferenceWeights_.setSize(fsiPatches_.size());

    if (!meshRefinementSupport_)
    {
        referenceFsiBuilt_ = false;
        referenceFaceModeProjection_.clear();
        referenceFsiFaceAreas_.clear();
    }
    else if (!referenceFsiBuilt_ ||
        referenceFaceModeProjection_.size() != fsiPatches_.size() ||
        referenceFsiFaceAreas_.size() != fsiPatches_.size())
    {
        referenceFaceModeProjection_.setSize(fsiPatches_.size());
        referenceFsiFaceAreas_.setSize(fsiPatches_.size());
        referenceFsiBuilt_ = false;
    }

    const labelList* oldPatchStartsPtr =
        mpmPtr ? &mpmPtr->oldPatchStarts() : nullptr;
    const labelList* oldPatchSizesPtr =
        mpmPtr ? &mpmPtr->oldPatchSizes() : nullptr;
    const List<objectMap>* facesFromFacesMapPtr =
        mpmPtr ? &mpmPtr->facesFromFacesMap() : nullptr;

    Map<labelList> facesFromFacesByNewFace;
    if (facesFromFacesMapPtr)
    {
        facesFromFacesByNewFace.resize(2 * facesFromFacesMapPtr->size() + 1);
        forAll(*facesFromFacesMapPtr, mapI)
        {
            const objectMap& obj = (*facesFromFacesMapPtr)[mapI];
            facesFromFacesByNewFace.set(obj.index(), obj.masterObjects());
        }
    }

    forAll(fsiPatches_, i)
    {
        const label patchID = pbm.findPatchID(fsiPatches_[i]);
        if (patchID == -1)
        {
            FatalErrorInFunction << "Configured FSI patch '" << fsiPatches_[i]
                                 << "' was not found in boundaryMesh."
                                 << exit(FatalError);
        }

        fsiPatchIDs_[i] = patchID;
        fsiPolyPatches_[i] = &pbm[patchID];

        const polyPatch& pp = *fsiPolyPatches_[i];
        const label nFaces = pp.size();

        pressureScaleCache_[i].setSize(nFaces, 1.0);
        faceModeProjection_[i].setSize(nMode_);

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            vectorField& coeff = faceModeProjection_[i][modeI];
            coeff.setSize(nFaces, vector::zero);

            forAll(pp, faceI)
            {
                const labelList& fPoints = pp[faceI];
                vector shapeFace = vector::zero;

                forAll(fPoints, fp)
                {
                    shapeFace += modeShapes_[modeI][fPoints[fp]];
                }

                if (fPoints.size() > 0)
                {
                    shapeFace /= fPoints.size();
                }

                coeff[faceI] = shapeFace;
            }
        }

        List<labelList>& faceToRefs = fsiFaceToReferenceFaces_[i];
        List<scalarField>& faceToRefWeights = fsiFaceToReferenceWeights_[i];
        faceToRefs.setSize(nFaces);
        faceToRefWeights.setSize(nFaces);

        if (!meshRefinementSupport_)
        {
            forAll(faceToRefs, faceI)
            {
                faceToRefs[faceI].setSize(1, faceI);
                faceToRefWeights[faceI].setSize(1, 1.0);
            }
            continue;
        }

        if (!referenceFsiBuilt_)
        {
            referenceFsiFaceAreas_[i].setSize(nFaces, 0.0);
            referenceFaceModeProjection_[i].setSize(nMode_);

            forAll(pp, faceI)
            {
                const vector area = pp.faceAreas()[faceI];
                const scalar areaMag = mag(area);
                referenceFsiFaceAreas_[i][faceI] = areaMag;
            }

            for (label modeI = 0; modeI < nMode_; ++modeI)
            {
                referenceFaceModeProjection_[i][modeI] =
                    faceModeProjection_[i][modeI];
            }

            forAll(faceToRefs, faceI)
            {
                faceToRefs[faceI].setSize(1, faceI);
                faceToRefWeights[faceI].setSize(1, 1.0);
            }

            continue;
        }

        if (!mpmPtr)
        {
            FatalErrorInFunction << "Reference cache is already initialized, "
                                    "but mapPolyMesh was not provided "
                                 << "while rebuilding topology mappings."
                                 << exit(FatalError);
        }

        if (!hasOldReferenceMapping)
        {
            FatalErrorInFunction
                << "Missing previous topology-to-reference mapping on patch '"
                << fsiPatches_[i] << "' while processing topology update."
                << exit(FatalError);
        }

        const label oldPatchStart = oldPatchStartsPtr->size() > patchID
            ? (*oldPatchStartsPtr)[patchID]
            : -1;
        const label oldPatchSize = oldPatchSizesPtr->size() > patchID
            ? (*oldPatchSizesPtr)[patchID]
            : -1;

        if (oldPatchStart < 0 || oldPatchSize < 0)
        {
            FatalErrorInFunction
                << "Cannot read old patch start/size for patch '"
                << fsiPatches_[i] << "' (patchID=" << patchID << ")."
                << exit(FatalError);
        }

        const List<labelList>& oldPatchRefs = oldFaceToReferenceFaces[i];
        const List<scalarField>& oldPatchWeights = oldFaceToReferenceWeights[i];

        if (oldPatchRefs.size() != oldPatchSize ||
            oldPatchWeights.size() != oldPatchSize)
        {
            FatalErrorInFunction
                << "Old mapping size mismatch on patch '" << fsiPatches_[i]
                << "': " << "oldPatchSize=" << oldPatchSize
                << ", refs=" << oldPatchRefs.size()
                << ", weights=" << oldPatchWeights.size() << exit(FatalError);
        }

        const scalarField& refAreas = referenceFsiFaceAreas_[i];

        label mappedCurrentFaces = 0;
        label oneToManyFaces = 0;

        forAll(pp, faceI)
        {
            labelList& refs = faceToRefs[faceI];
            scalarField& weights = faceToRefWeights[faceI];
            refs.clear();
            weights.clear();

            const label globalFace = pp.start() + faceI;
            DynamicList<label> oldFaces;

            const label oldMasterFace = mpmPtr->faceMap()[globalFace];
            if (oldMasterFace >= 0)
            {
                oldFaces.append(oldMasterFace);
            }

            const auto mergedIter = facesFromFacesByNewFace.cfind(globalFace);
            if (mergedIter.found())
            {
                const labelList& masters = mergedIter();
                forAll(masters, masterI)
                {
                    oldFaces.appendUniq(masters[masterI]);
                }
            }

            if (oldFaces.empty())
            {
                FatalErrorInFunction << "Cannot derive old faces for patch '"
                                     << fsiPatches_[i] << "', face " << faceI
                                     << " (global face " << globalFace << ")."
                                     << exit(FatalError);
            }

            Map<scalar> refAreaAccum(2 * oldFaces.size() + 1);

            forAll(oldFaces, oldI)
            {
                const label oldFace = oldFaces[oldI];
                if (oldFace < oldPatchStart ||
                    oldFace >= oldPatchStart + oldPatchSize)
                {
                    continue;
                }

                const label oldLocalFace = oldFace - oldPatchStart;
                const labelList& prevRefs = oldPatchRefs[oldLocalFace];
                const scalarField& prevWeights = oldPatchWeights[oldLocalFace];

                if (prevRefs.size() != prevWeights.size() || prevRefs.empty())
                {
                    FatalErrorInFunction
                        << "Invalid old mapping entry for patch '"
                        << fsiPatches_[i] << "', old local face "
                        << oldLocalFace << "." << exit(FatalError);
                }

                forAll(prevRefs, refI)
                {
                    const label refFace = prevRefs[refI];
                    if (refFace < 0 || refFace >= refAreas.size())
                    {
                        FatalErrorInFunction << "Out-of-range reference face "
                                             << refFace << " on patch '"
                                             << fsiPatches_[i] << "'."
                                             << exit(FatalError);
                    }

                    const scalar clampedWeight = max(prevWeights[refI], 0.0);
                    const scalar contribArea =
                        max(refAreas[refFace], VSMALL) * clampedWeight;
                    refAreaAccum(refFace, 0.0) += contribArea;
                }
            }

            if (refAreaAccum.empty())
            {
                FatalErrorInFunction
                    << "Topology mapping on patch '" << fsiPatches_[i]
                    << "', face " << faceI
                    << " yielded no reference faces after patch filtering."
                    << exit(FatalError);
            }

            const labelList sortedRefs(refAreaAccum.sortedToc());
            refs.setSize(sortedRefs.size());
            weights.setSize(sortedRefs.size(), 0.0);

            scalar sumArea = 0.0;
            forAll(sortedRefs, refI)
            {
                sumArea += refAreaAccum[sortedRefs[refI]];
            }

            if (sumArea <= VSMALL)
            {
                FatalErrorInFunction
                    << "Non-positive accumulated reference area on patch '"
                    << fsiPatches_[i] << "', face " << faceI << "."
                    << exit(FatalError);
            }

            scalar sumW = 0.0;
            forAll(sortedRefs, refI)
            {
                refs[refI] = sortedRefs[refI];
                weights[refI] = refAreaAccum[sortedRefs[refI]] / sumArea;
                sumW += weights[refI];
            }

            if (sumW <= VSMALL)
            {
                FatalErrorInFunction << "Zero total mapping weight on patch '"
                                     << fsiPatches_[i] << "', face " << faceI
                                     << "." << exit(FatalError);
            }

            if (mag(sumW - 1.0) > 1e-10)
            {
                weights /= sumW;
            }

            ++mappedCurrentFaces;
            if (refs.size() > 1)
            {
                ++oneToManyFaces;
            }
        }

        label globalMappedCurrentFaces = mappedCurrentFaces;
        label globalOneToManyFaces = oneToManyFaces;
        reduce(globalMappedCurrentFaces, sumOp<label>());
        reduce(globalOneToManyFaces, sumOp<label>());

        if (Pstream::master())
        {
            Info << "Topology reference mapping on patch '" << fsiPatches_[i]
                 << "': mapped " << globalMappedCurrentFaces
                 << " current faces, one-to-many faces=" << globalOneToManyFaces
                 << "." << endl;
        }
    }

    if (meshRefinementSupport_ && !referenceFsiBuilt_)
    {
        referenceFsiBuilt_ = true;
        if (Pstream::master())
        {
            Info << "Initialized reference FSI face cache for "
                    "refinement/unrefine-aware "
                 << "force aggregation." << endl;
        }
    }
}

void fastDynamicFvMesh::syncMappedModeShapes(boolList& mappedMask)
{
    if (!meshRefinementSupport_ || !Pstream::parRun() || nMode_ <= 0)
    {
        return;
    }

    if (mappedMask.size() != modeShapes_[0].size())
    {
        FatalErrorInFunction << "mappedMask size (" << mappedMask.size()
                             << ") does not match mode-shape point count ("
                             << modeShapes_[0].size() << ")."
                             << exit(FatalError);
    }

    List<label> ownerProc(mappedMask.size(), -1);
    const label myProc = Pstream::myProcNo();

    forAll(mappedMask, pointI)
    {
        if (mappedMask[pointI])
        {
            ownerProc[pointI] = myProc;
        }
    }

    syncTools::syncPointList(*this, ownerProc, maxEqOp<label>(), label(-1));
    syncTools::syncPointList(*this, mappedMask, orEqOp<bool>(), false);

    for (label modeI = 0; modeI < nMode_; ++modeI)
    {
        vectorField ownedValues(modeShapes_[modeI].size(), vector::zero);

        forAll(ownedValues, pointI)
        {
            if (ownerProc[pointI] == myProc)
            {
                ownedValues[pointI] = modeShapes_[modeI][pointI];
            }
        }

        syncTools::syncPointList(
            *this, ownedValues, plusEqOp<vector>(), vector::zero);
        modeShapes_[modeI].transfer(ownedValues);
    }
}

void fastDynamicFvMesh::syncSharedPointPositions(pointField& pointValues) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (pointValues.size() != this->points().size())
    {
        FatalErrorInFunction << "pointValues size (" << pointValues.size()
                             << ") does not match mesh point count ("
                             << this->points().size() << ")."
                             << exit(FatalError);
    }

    List<label> ownerProc(pointValues.size(), Pstream::myProcNo());
    syncTools::syncPointList(*this, ownerProc, maxEqOp<label>(), label(-1));

    pointField ownedValues(pointValues.size(), point::zero);
    const label myProc = Pstream::myProcNo();

    forAll(ownedValues, pointI)
    {
        if (ownerProc[pointI] == myProc)
        {
            ownedValues[pointI] = pointValues[pointI];
        }
    }

    syncTools::syncPointPositions(
        *this, ownedValues, plusEqOp<point>(), point::zero);
    pointValues.transfer(ownedValues);
}

label fastDynamicFvMesh::interpolateModeShapesFromKnown(boolList& knownMask,
    const boolList& sourceKnownMask, const label cornerCount,
    const bool distanceWeighted, const word& context)
{
    if (cornerCount <= 0)
    {
        return 0;
    }

    const pointField& pts = this->points();
    const labelListList& pCells = this->pointCells();
    const labelListList& cPoints = this->cellPoints();
    label localInterpolated = 0;

    forAll(knownMask, pointI)
    {
        if (knownMask[pointI])
        {
            continue;
        }

        DynamicList<label> candidates;
        const labelList& localCells = pCells[pointI];

        forAll(localCells, cellI)
        {
            const label cI = localCells[cellI];
            const labelList& cellPts = cPoints[cI];

            forAll(cellPts, cpI)
            {
                const label corner = cellPts[cpI];
                if (corner == pointI || !sourceKnownMask[corner])
                {
                    continue;
                }

                bool duplicate = false;
                forAll(candidates, candI)
                {
                    if (candidates[candI] == corner)
                    {
                        duplicate = true;
                        break;
                    }
                }

                if (!duplicate)
                {
                    candidates.append(corner);
                }
            }
        }

        if (candidates.size() < cornerCount)
        {
            continue;
        }

        List<label> order(candidates.size(), -1);
        List<scalar> dist(candidates.size(), GREAT);
        forAll(candidates, candI)
        {
            order[candI] = candI;
            dist[candI] = mag(pts[pointI] - pts[candidates[candI]]);
        }

        std::sort(order.begin(), order.end(),
            [&](const label a, const label b) { return dist[a] < dist[b]; });

        // Select corners by minimizing geometric centroid mismatch at this
        // point (within nearest local candidates), which is more robust than
        // blindly taking the first N nearest corners on distorted AMR meshes.
        const label maxSearchCandidates =
            std::min<label>(order.size(), cornerCount == 2 ? 16 : 12);

        DynamicList<label> corners;
        scalar bestCentroidErr = GREAT;
        scalar bestAvgDist = GREAT;

        if (maxSearchCandidates >= cornerCount)
        {
            DynamicList<label> trial(cornerCount);

            auto updateBest = [&](const DynamicList<label>& labels) {
                point centroid = point::zero;
                scalar avgDist = 0.0;

                forAll(labels, li)
                {
                    const label corner = labels[li];
                    centroid += pts[corner];
                    avgDist += mag(pts[pointI] - pts[corner]);
                }

                centroid /= scalar(labels.size());
                avgDist /= scalar(labels.size());

                const scalar centroidErr = mag(pts[pointI] - centroid);

                if (centroidErr < bestCentroidErr - SMALL ||
                    (mag(centroidErr - bestCentroidErr) <= SMALL &&
                        avgDist < bestAvgDist))
                {
                    bestCentroidErr = centroidErr;
                    bestAvgDist = avgDist;
                    corners.clear();
                    forAll(labels, li)
                    {
                        corners.append(labels[li]);
                    }
                }
            };

            auto chooseCorners = [&](const auto& self, const label start,
                                     const label need) -> void {
                if (need == 0)
                {
                    updateBest(trial);
                    return;
                }

                const label maxI = maxSearchCandidates - need;
                for (label rankI = start; rankI <= maxI; ++rankI)
                {
                    trial.append(candidates[order[rankI]]);
                    self(self, rankI + 1, need - 1);
                    trial.remove();
                }
            };

            chooseCorners(chooseCorners, 0, cornerCount);
        }

        // Safety fallback: retain legacy nearest-corner selection if no
        // valid combination was found in the bounded search.
        if (corners.size() != cornerCount)
        {
            corners.clear();
            for (label pickI = 0; pickI < cornerCount; ++pickI)
            {
                corners.append(candidates[order[pickI]]);
            }
        }

        if (corners.size() != cornerCount)
        {
            continue;
        }

        scalarField weights(cornerCount, 1.0 / scalar(cornerCount));

        if (cornerCount == 2 && distanceWeighted)
        {
            const scalar d0 = mag(pts[pointI] - pts[corners[0]]);
            const scalar d1 = mag(pts[pointI] - pts[corners[1]]);
            const scalar denom = d0 + d1;

            if (denom > VSMALL)
            {
                weights[0] = d1 / denom;
                weights[1] = d0 / denom;
            }
            else
            {
                weights[0] = 0.5;
                weights[1] = 0.5;
            }
        }
        else if (cornerCount == 4 || cornerCount == 8)
        {
            // For hex-like AMR topology (face/cell midpoint creation), equal
            // corner weights preserve straight-line/plane consistency better.
            weights = 1.0 / scalar(cornerCount);
        }
        else
        {
            scalar sumInv = 0.0;
            forAll(corners, cI)
            {
                const scalar d = mag(pts[pointI] - pts[corners[cI]]);
                const scalar inv = 1.0 / max(d, refinementInterpTolerance_);
                weights[cI] = inv;
                sumInv += inv;
            }

            if (sumInv > VSMALL)
            {
                weights /= sumInv;
            }
            else
            {
                weights = 1.0 / scalar(cornerCount);
            }
        }

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            vector value = vector::zero;
            forAll(corners, cI)
            {
                value += weights[cI] * modeShapes_[modeI][corners[cI]];
            }
            modeShapes_[modeI][pointI] = value;
        }

        knownMask[pointI] = true;
        ++localInterpolated;
    }

    if (Pstream::master() && localInterpolated > 0)
    {
        Info << "Refinement interpolation (" << context << ") updated "
             << localInterpolated << " local points using " << cornerCount
             << "-corner topology rule." << endl;
    }

    return localInterpolated;
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

    const label newNPoints = this->points().size();
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


} // End namespace Foam
