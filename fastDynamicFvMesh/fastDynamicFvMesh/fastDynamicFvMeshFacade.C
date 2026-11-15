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

defineTypeNameAndDebug(fastDynamicFvMesh, 0);
addToRunTimeSelectionTable(dynamicFvMesh, fastDynamicFvMesh, IOobject);

namespace
{
    //- Should the mesh geometry be initialised before the refinement engine
    //  is constructed?
    //
    //  The dynamicRefineFvMesh base constructs its hexRef8 refinement engine
    //  in the constructor initialiser list. hexRef8's constructor runs
    //  checkMesh(), which validates the anchors of coupled-patch faces via
    //  syncTools::swapBoundaryFaceList(). For patches that need a
    //  transformation (cyclic rotational/translational), those anchors can
    //  only be compared once the coupled-patch transformation tensors have
    //  been calculated by polyMesh::init(). polyMesh::init() only runs from
    //  the polyMesh constructor when the base is constructed with doInit=true;
    //  without it, transformed cyclic patches abort with
    //  "Faces do not seem to be correct across coupled boundaries".
    //
    //  The eager path makes dynamicMotionSolverListFvMesh construct its
    //  mandatory motion solver (unless a 'solvers' sub-dictionary is present,
    //  which may be empty to construct none). Opt in explicitly through
    //  dynamicMeshDict:
    //    - 'solvers' sub-dictionary present  -> eager initialisation,
    //    - otherwise the legacy lazy path is kept.
    bool eagerMeshInit(const IOobject& io)
    {
        IOdictionary dynamicMeshDict
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                polyMesh::regionName(io.name()),
                io.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        return dynamicMeshDict.found("solvers");
    }
}


fastDynamicFvMesh::fastDynamicFvMesh(const IOobject& io)
    : dynamicRefineFvMesh(io, eagerMeshInit(io)), FdmCore(io, *this),
      runtimeRefinementEnabled_(false),
      refinementMappingDiagnostics_(false), refinementMinCellVolume_(0.0),
      refinementMinEdgeLength_(0.0), refinementUseGradIndicator_(false),
      refinementGradIndicatorField_("alpha.water"),
      refinementGradIndicatorMagnitudeField_(word::null),
      refinementGradIndicatorCellCache_(0),
      refinementGradIndicatorPointCache_(0),
      refinementGradIndicatorCacheTimeIndex_(-1),
      refinementGradIndicatorCacheNCells_(-1),
      refinementGradIndicatorCacheNPoints_(-1),
      refinementGradIndicatorCacheSourceField_(word::null),
      refinementGradIndicatorCacheMagnitudeField_(word::null),
      refinementGradIndicatorCacheValid_(false),
      amrLayerDiagnosticsEnabled_(false), amrLayerDiagnosticsZ_(0.0),
      amrLayerDiagnosticsThickness_(0.0),
      amrLayerDiagnosticsNeighborDepth_(1), amrLayerDiagnosticsMaxRows_(0),
      amrLayerDiagnosticsFilePrefix_("amr_layer_diagnostics"),
      amrLayerDiagnosticsSummaryHeaderWritten_(false),
      amrLayerDiagnosticsCellHeaderWritten_(false),
      amrLayerDiagnosticsCandidateCell_(0),
      amrLayerDiagnosticsIndicatorCell_(0),
      amrLayerDiagnosticsLower_(0.0), amrLayerDiagnosticsUpper_(0.0),
      amrLayerDiagnosticsLayerCells_(0),
      amrLayerDiagnosticsNearbyCells_(0),
      amrLayerDiagnosticsCandidates_(0), amrLayerDiagnosticsSelected_(0),
      amrLayerDiagnosticsBelowLower_(0), amrLayerDiagnosticsAboveUpper_(0),
      amrLayerDiagnosticsRejectedSizeFloor_(0),
      amrLayerDiagnosticsRejectedBase_(0),
      amrLayerDiagnosticsUnrefinePoints_(0),
      amrLayerDiagnosticsPendingSummary_(false), trackAmrTimingEnabled_(false)
{
    init(true);
}

bool fastDynamicFvMesh::init(const bool doInit)
{
    IOdictionary dynamicMeshDict(
        IOobject("dynamicMeshDict", fdmMesh_.time().constant(), *this,
            IOobject::MUST_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER));

    runtimeRefinementEnabled_ =
        dynamicMeshDict.found("dynamicRefineFvMeshCoeffs");

    if (runtimeRefinementEnabled_)
    {
        // refineInterval 0 is dynamicRefineFvMesh's own "no refinement"
        // setting (updateTopology returns immediately); honour it here so
        // such dictionaries do not take the runtime-AMR path.
        runtimeRefinementEnabled_ =
            dynamicMeshDict.subDict("dynamicRefineFvMeshCoeffs")
                .getOrDefault<label>("refineInterval", 1) != 0;
    }

    if (runtimeRefinementEnabled_)
    {
        dynamicRefineFvMesh::init(doInit);
    }
    else
    {
        dynamicFvMesh::init(doInit);
    }

    readControls();
    readAmrControls
    (
        dynamicMeshDict.subDict(fdmTypeName() + "Coeffs"),
        dynamicMeshDict
    );
    initialiseCoreState();

    return true;
}

fastDynamicFvMesh::~fastDynamicFvMesh() {}

bool fastDynamicFvMesh::update()
{
    const label currentTimeIndex = fdmMesh_.time().timeIndex();
    const bool firstIter = isNewTimeStep(currentTimeIndex);

    if (firstIter)
    {
        updateRefinementTopology();
    }

    const bool coreChanged = updateFdm();

    return runtimeRefinementEnabled_ || coreChanged;
}

bool fastDynamicFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok = dynamicRefineFvMesh::writeObject(streamOpt, writeOnProc);

    if (ok && writeOnProc && meshRefinementSupport_)
    {
        ok = writeRefinedPointModeShapes() && ok;
        ok = writeFsiFaceReferenceMap() && ok;
    }

    return ok;
}

} // End namespace Foam
