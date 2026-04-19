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
      refinementInterpTolerance_(1e-8), refinementMinCellVolume_(0.0),
      refinementMinEdgeLength_(0.0), refinementUseGradIndicator_(false),
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

    if (trackTimingEnabled_)
    {
        if (timingHasLastUpdate_)
        {
            const scalar fluidCpu =
                updateStartCpuTime - timingLastUpdateCpuTime_;
            if (fluidCpu > 0)
            {
                timingFluidCpuAccum_ += fluidCpu;
            }
        }
        else
        {
            timingHasLastUpdate_ = true;
        }
    }

    const label currentTimeIndex = this->time().timeIndex();
    const bool firstIter = (currentTimeIndex != lastUpdateTimeIndex_);

    if (runtimeRefinementEnabled_ && firstIter)
    {
        dynamicRefineFvMesh::updateTopology();
    }

    lastUpdateTimeIndex_ = currentTimeIndex;

    if (firstIter)
    {
        updateCount_ = 0;
        // Store previous state only on the first iteration of the time step
        modeForce0_ = modeForce_;
        modeState0_ = modeState_;
    }
    updateCount_++;

    // 1. Calculate time step
    scalar dt = this->time().deltaTValue();

    // 2. Calculate fluid forces and optional structural external loading
    calcModalForces();
    calcStructuralModalForces();

    forAll(modeForce_, modeI)
    {
        modeForce_[modeI] += structuralModeForce_[modeI];
    }

    // Force relaxation and residual calculation
    //   ON the very first iteration of a time step, we are predicting the
    // force. If we just use modeForce0_, we are assuming force doesn't change.
    // But appliedModeForce_ holds the LAST applied force from the PREVIOUS time
    // step's convergence.
    //   HOWEVER, on a restart, appliedModeForce_ might be stale or just read.
    // We should ensure that we don't relax against a zero vector if the force
    // is large.
    //   LET'S rely on the constructor/initialization logic for restart cases.
    // For normal running, appliedModeForce_ carries over.

    fsiResidual_ = 0.0;

    // Temporary storage for next iteration forces
    scalarField nextAppliedForce(nMode_);

    for (label i = 0; i < nMode_; ++i)
    {
        scalar rawFluidForce = modeForce_[i];
        scalar prevRelaxedForce = appliedModeForce_[i];

        // Relax: F_new = F_old + alpha * (F_fluid - F_old)
        // Here F_fluid is the new target, F_old is where we were.
        scalar newRelaxedForce = prevRelaxedForce +
            couplingRelaxation_ * (rawFluidForce - prevRelaxedForce);

        // If this is the first iteration of the time step, and we are
        // RESTARTING (or startup), we might want to trust the computed force
        // more if the previous force is suspiciously zero while the new force
        // is large.
        if (firstIter && mag(prevRelaxedForce) < VSMALL &&
            mag(rawFluidForce) > 1.0)
        {
            // Likely a bad initialization on restart or startup
            newRelaxedForce = rawFluidForce;

            if (Pstream::master() && i == 0)
            {
                Info << "  First iter override: using raw force "
                     << rawFluidForce << " instead of relaxed "
                     << newRelaxedForce << endl;
            }
        }

        // Calculate residual relative to the NEW force magnitude
        // Convergence is when F_fluid matches F_old (correction goes to zero)
        // or when F_new matches F_old (change is zero).
        scalar diff = mag(newRelaxedForce - prevRelaxedForce);
        scalar magForce = mag(newRelaxedForce) + 1e-12; // Avoid div by zero
        scalar res = diff / magForce;

        if (res > fsiResidual_)
            fsiResidual_ = res;

        // Store for next step
        nextAppliedForce[i] = newRelaxedForce;
    }

    // Update state
    appliedModeForce_ = nextAppliedForce;
    modeForce_ = nextAppliedForce; // Use relaxed force for structural dynamics

    // Store FSI residual in object registry for solvers to read
    if (this->foundObject<uniformDimensionedScalarField>("fsiResidual"))
    {
        const_cast<uniformDimensionedScalarField&>(
            this->lookupObject<uniformDimensionedScalarField>("fsiResidual"))
            .value() = fsiResidual_;
    }
    else
    {
        uniformDimensionedScalarField* fsiResPtr =
            new uniformDimensionedScalarField(
                IOobject("fsiResidual", this->time().constant(), *this,
                    IOobject::NO_READ, IOobject::NO_WRITE),
                dimless, fsiResidual_);
        fsiResPtr->store();
    }

    if (startupStepCount_ < 2)
    {
        if (Pstream::master())
        {
            Info << "Initializing Fast Dynamic Mesh state (startup step "
                 << (startupStepCount_ + 1) << " of 2)." << endl;
            if (structuralForceEnabled_)
            {
                label totalTargets = 0;
                forAll(structuralForceNodeIDs_, fi)
                {
                    totalTargets += structuralForceNodeIDs_[fi].size();
                }

                Info << "  Structural load files: " << nStructuralForces_
                     << ", mapped target nodes(total unique-per-file): "
                     << totalTargets << endl;
            }
        }

        for (label i = 0; i < nMode_; ++i)
        {
            scalar F = modeForce_[i];
            scalar v = 0.0;
            if (i < initVelocity_.size())
                v = initVelocity_[i];

            // Legacy Fluent UDF initialisation:
            // D0 = 0, V0 = initVelocity, A0 = F0
            scalar a = F;

            modeState_[i] = vector(0.0, v, a);
            modeState0_[i] = modeState_[i];
            modeForce0_[i] = modeForce_[i];
            appliedModeDisp_[i] = 0.0;
        }

        ++startupStepCount_;
        writeDiagnostics();

        if (trackTimingEnabled_)
        {
            const scalar updateEndCpuTime = this->time().elapsedCpuTime();
            const scalar meshCpu = updateEndCpuTime - updateStartCpuTime;
            if (meshCpu > 0)
            {
                timingMeshCpuAccum_ += meshCpu;
            }
            timingLastUpdateCpuTime_ = updateEndCpuTime;
        }

        return runtimeRefinementEnabled_ || this->moving() ||
            this->topoChanging();
    }

    // 3. Solve dynamics
    if (nMode_ > 0)
    {
        solveStructuralDynamics(dt);
    }

    // 4. Update mesh points
    pointField newPoints = this->points();

    for (label m = 0; m < nMode_; ++m)
    {
        const vectorField& shape = modeShapes_[m];

        // Old displacement applied to the mesh
        const scalar appliedDisp0 = appliedModeDisp_[m];

        // New displacement from structural solver (using relaxed forces)
        const scalar targetDisp = modeState_[m].x();

        // Calculate increment
        scalar dDisp = targetDisp - appliedDisp0;

        // Apply limiter
        if (mag(dDisp) > maxDispChange_)
        {
            if (Pstream::master() && updateCount_ == 1 &&
                m == 0) // Log once per step
            {
                WarningInFunction
                    << "Limiting displacement change for Mode " << m << " from "
                    << dDisp << " to "
                    << (dDisp > 0 ? maxDispChange_ : -maxDispChange_) << endl;
            }
            dDisp = (dDisp > 0) ? maxDispChange_ : -maxDispChange_;
        }

        // Final applied displacement
        const scalar appliedDisp = appliedDisp0 + dDisp;

        // Update the physics state to be consistent with the applied
        // displacement This effectively adds a constraint force implicitly
        modeState_[m].x() = appliedDisp;
        appliedModeDisp_[m] = appliedDisp;

        forAll(newPoints, pI)
        {
            newPoints[pI] += dDisp * shape[pI];
        }
    }

    syncSharedPointPositions(newPoints);
    this->movePoints(newPoints);
    writeDiagnostics();

    if (Pstream::master() && structuralForceEnabled_ && updateCount_ == 1 &&
        nMode_ > 0)
    {
        Info << "Structural modal load summary at t=" << this->time().timeName()
             << ": mode0=" << structuralModeForce_[0]
             << ", signal=" << structuralForceSignal_ << endl;
    }

    if (trackTimingEnabled_)
    {
        const scalar updateEndCpuTime = this->time().elapsedCpuTime();
        const scalar meshCpu = updateEndCpuTime - updateStartCpuTime;
        if (meshCpu > 0)
        {
            timingMeshCpuAccum_ += meshCpu;
        }
        timingLastUpdateCpuTime_ = updateEndCpuTime;
    }

    // Persist restart modal states in a dedicated channel under
    // constant/fsiRestart to avoid mixing restart metadata with reconstructed
    // field time directories.
    if (Pstream::master() && this->time().writeTime() &&
        currentTimeIndex != lastGlobalWriteTimeIndex_)
    {
        lastGlobalWriteTimeIndex_ = currentTimeIndex;

        writeTimingReport();

        writeRestartStateStore();
        Info << "Wrote restart modal state files to " << restartStateInstance()
             << " at t=" << this->time().timeName() << endl;
    }

    return runtimeRefinementEnabled_ || this->moving() || this->topoChanging();
}

} // End namespace Foam
