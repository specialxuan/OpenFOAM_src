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
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fastDynamicFvMesh, 0);
addToRunTimeSelectionTable(dynamicFvMesh, fastDynamicFvMesh, IOobject);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Constructor
// Purpose: initialize the fastDynamicFvMesh instance and load runtime controls
// and mode shapes. Notes:
//  - The initializer list (dynamicFvMesh(io), ...) sets default member values
//  before the body runs.
//  - readControls() parses constant/dynamicMeshDict and validates required
//  entries (fsiPatches, couplingRelaxation, etc.).
//  - readModeShapes() loads modal geometry and displacement shapes from the
//  case 'mode' directory and builds the mapping
//    between CSV node coordinates and local mesh points. This mapping is
//    intentionally brute-force for correctness.
fastDynamicFvMesh::fastDynamicFvMesh(const IOobject& io)
    : dynamicRefineFvMesh(io, false), nMode_(0),
      // Initialize restart containers; dedicated restart I/O is handled
      // explicitly in readRestartStateStore()/writeRestartStateStore().
      modeForce_(
          IOobject("modeForce", io.time().constant() / fileName("fsiRestart"),
              *this, IOobject::NO_READ, IOobject::NO_WRITE)),
      modeState_(
          IOobject("modeState", io.time().constant() / fileName("fsiRestart"),
              *this, IOobject::NO_READ, IOobject::NO_WRITE)),
      initVelocity_(0), appliedModeDisp_(IOobject("appliedModeDisp",
                            io.time().constant() / fileName("fsiRestart"),
                            *this, IOobject::NO_READ, IOobject::NO_WRITE)),
      appliedModeForce_(IOobject("appliedModeForce",
          io.time().constant() / fileName("fsiRestart"), *this,
          IOobject::NO_READ, IOobject::NO_WRITE)),
      fsiResidual_(1.0), theta_(1.4), // Default Wilson-Theta
      mappingTolerance_(4e-6), couplingRelaxation_(1.0),
      pressureFieldName_("p"), rhoRef_(-1.0), pRef_(0.0), maxDispChange_(1.0),
      startupStepCount_(0), lastUpdateTimeIndex_(-1), updateCount_(0),
      writeDiagnosticsEnabled_(false), trackTimingEnabled_(false),
      timingLastUpdateCpuTime_(0.0), timingHasLastUpdate_(false),
      timingFluidCpuAccum_(0.0), timingMeshCpuAccum_(0.0),
      timingFluidCpuAtLastWrite_(0.0), timingMeshCpuAtLastWrite_(0.0),
      lastGlobalWriteTimeIndex_(-1), fsiPatchIDs_(0), fsiPolyPatches_(0),
      faceModeProjection_(0), pressureScaleCache_(0),
      pressureScaleInitialized_(0), diagnosticsHeaderWritten_(false),
      modelPointersCached_(false), icoTurbPtr_(nullptr), cmpTurbPtr_(nullptr),
      fluidThermoPtr_(nullptr), laminarTransportPtr_(nullptr),
      structuralForceEnabled_(false), nStructuralForces_(0),
      structuralForceFilePrefix_("StructForce"),
      structuralTargetTolerance_(1e-8),
      structNodeCoorFile_("StructNodeCoor.csv"),
      structNodeDispPrefix_("StructNodeDisp"), structNodeCoords_(0),
      structModeShapes_(0), structuralForceNodeIDs_(0),
      structuralForceTimes_(0), structuralForceValues_(0),
      structuralModeForce_(0), structuralForceSignal_(0.0),
      meshRefinementSupport_(false), runtimeRefinementEnabled_(false),
      refinementInterpTolerance_(1e-8), refinementMappingDiagnostics_(false),
      refinementMinCellVolume_(0.0), refinementMinEdgeLength_(0.0),
      refinementUseGradIndicator_(false),
      refinementGradIndicatorField_("alpha.water"),
      refinementGradIndicatorMagnitudeField_(word::null),
      refinementGradIndicatorCellCache_(0), refinementGradIndicatorPointCache_(0),
      refinementGradIndicatorCacheTimeIndex_(-1),
      refinementGradIndicatorCacheNCells_(-1),
      refinementGradIndicatorCacheNPoints_(-1),
      refinementGradIndicatorCacheSourceField_(word::null),
      refinementGradIndicatorCacheMagnitudeField_(word::null),
      refinementGradIndicatorCacheValid_(false), referenceFsiBuilt_(false),
      referenceFaceModeProjection_(0), fsiFaceToReferenceFaces_(0),
      fsiFaceToReferenceWeights_(0), referenceFsiFaceAreas_(0)
{
    init(true);
}

bool fastDynamicFvMesh::init(const bool doInit)
{
    IOdictionary dynamicMeshDict(
        IOobject("dynamicMeshDict", this->time().constant(), *this,
            IOobject::MUST_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER));

    runtimeRefinementEnabled_ =
        dynamicMeshDict.found("dynamicRefineFvMeshCoeffs");

    if (runtimeRefinementEnabled_)
    {
        dynamicRefineFvMesh::init(doInit);
    }
    else
    {
        dynamicFvMesh::init(doInit);
    }

    readControls();
    readModeShapes();

    const bool restartStateRead = readRestartStateStore();
    const bool stateRead = restartStateRead && (modeState_.size() > 0);

    if (stateRead && modeState_.size() == nMode_)
    {
        Info << "Restarting fastDynamicFvMesh with " << nMode_
             << " modes from time " << this->time().timeName() << endl;

        startupStepCount_ = 2;

        if (modeForce_.size() != nMode_)
        {
            modeForce_.setSize(nMode_, 0.0);
        }

        if (appliedModeDisp_.size() != nMode_)
        {
            Info << "  appliedModeDisp not read or wrong size; initializing "
                    "from modeState displacement."
                 << endl;
            appliedModeDisp_.setSize(nMode_);
            for (label i = 0; i < nMode_; ++i)
            {
                appliedModeDisp_[i] = modeState_[i].x();
            }
        }

        if (appliedModeForce_.size() != nMode_)
        {
            Info << "  appliedModeForce not read; initializing from modeForce."
                 << endl;
            if (modeForce_.size() == nMode_)
            {
                appliedModeForce_ = modeForce_;
            }
            else
            {
                appliedModeForce_.setSize(nMode_, 0.0);
            }
        }

        modeState0_ = modeState_;
        modeForce0_ = modeForce_;

        if (nMode_ > 0 && mag(appliedModeForce_[0]) < VSMALL &&
            mag(modeForce0_[0]) > VSMALL)
        {
            Info << "  Initializing appliedModeForce from modeForce0 to avoid "
                    "relaxation shock."
                 << endl;
            appliedModeForce_ = modeForce0_;
        }
    }
    else
    {
        if (restartStateRead)
        {
            WarningInFunction << "Read modeState size " << modeState_.size()
                              << " does not match nMode " << nMode_
                              << ". Resetting state." << endl;
        }

        modeForce_.setSize(nMode_, 0.0);
        modeState_.setSize(nMode_, vector::zero);
        appliedModeDisp_.setSize(nMode_, 0.0);
        appliedModeForce_.setSize(nMode_, 0.0);

        startupStepCount_ = 0;
    }

    return true;
}

fastDynamicFvMesh::~fastDynamicFvMesh() {}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// readControls()
// Purpose: parse the dynamicMeshDict sub-dictionary <typeName>Coeffs and
// initialize control parameters. Behavior overview:
//  - Open constant/dynamicMeshDict and obtain the sub-dictionary named
//  "fastDynamicFvMeshCoeffs" (typeName + "Coeffs").
//  - Read optional parameters (theta, mappingTolerance, couplingRelaxation,
//  pressureField, rhoRef).
//  - Read required entry 'fsiPatches' as the list of fluid-solid interface
//  patches.
//  - Validate ranges (e.g., couplingRelaxation must be in (0,1], rhoRef
//  positive if specified).
//  - Handle face-level diagnostics settings (faceDiagnosticsMode
//  conversion from 1-based to 0-based).
// update()
// Purpose: main update hook called by OpenFOAM each timestep to compute modal
// forces, solve structural dynamics, and move the mesh. Sequence:
//  1) Guard against duplicate updates for the same physical time index using
//  lastUpdateTimeIndex_. 2) Compute dt from time().deltaTValue() and save
//  previous modal state/forces for time-integration. 3) Call calcModalForces()
//  to assemble current modal loads (pressure + shear separated for
//  diagnostics). 4) Startup initialization: for the first two update() calls,
//  initialize modeState_ according to legacy conventions (disp=0,
//  vel=initVelocity, acc=modeForce). 5) Solve modal dynamics with
//  solveStructuralDynamics(dt) to update modeState_. 6) Apply
//  couplingRelaxation_ to compute appliedModeDisp_ = old + alpha*(new-old),
//  then compute incremental displacement dDisp = applied-newApplied. 7)
//  Reconstruct point displacements by summing dDisp*shape for each mapped mode
//  shape and call movePoints(newPoints). 8) Write diagnostics and return true
//  to indicate the mesh was updated.
// Notes:
//  - The relaxation step is crucial for stability in explicit partitioned
//  coupling; alpha in (0,1] reduces added-mass instability.
//  - The function writes debug info on the master to help trace
//  divergence/stability issues.
bool fastDynamicFvMesh::update()
{
    const scalar updateStartCpuTime = this->time().elapsedCpuTime();
    accountFluidCpuAtUpdateStart(updateStartCpuTime);

    const label currentTimeIndex = this->time().timeIndex();
    const bool firstIter = beginUpdateStep(currentTimeIndex);

    const scalar dt = this->time().deltaTValue();
    assembleAndRelaxModalForces(firstIter);

    if (initializeStartupStateIfNeeded())
    {
        accountMeshCpuAtUpdateEnd(updateStartCpuTime);
        return updateReturnState();
    }

    solveAndMoveMesh(dt);
    logStructuralLoadSummary();
    accountMeshCpuAtUpdateEnd(updateStartCpuTime);
    writeUpdateOutputs(currentTimeIndex);

    return updateReturnState();
}

} // End namespace Foam
