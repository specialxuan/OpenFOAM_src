/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    FdmCore implementation: modal/FSI state, force assembly, structural
    dynamics, mode mapping and modal restart. No AMR topology code.

\*---------------------------------------------------------------------------*/
#include "Pstream.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"
#include "fdmCore.H"
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

FdmCore::FdmCore
(
    const IOobject& io,
    dynamicMotionSolverListFvMesh& mesh
)
    : fdmMesh_(mesh), nMode_(0),
      nAvailableModes_(0), nActiveModes_(-1),
      // Initialize restart containers; dedicated restart I/O is handled
      // explicitly in readRestartStateStore()/writeRestartStateStore().
      modeForce_(
          IOobject("modeForce", io.time().constant() / fileName("fsiRestart"),
              fdmMesh_, IOobject::NO_READ, IOobject::NO_WRITE)),
      modeState_(
          IOobject("modeState", io.time().constant() / fileName("fsiRestart"),
              fdmMesh_, IOobject::NO_READ, IOobject::NO_WRITE)),
      initVelocity_(0),
      prescribedMotion_(false), prescribedAmplitude_(0), prescribedPhase_(0),
      prescribedFrequency_(-1.0),
      appliedModeDisp_(IOobject("appliedModeDisp",
                            io.time().constant() / fileName("fsiRestart"),
                            fdmMesh_, IOobject::NO_READ, IOobject::NO_WRITE)),
      appliedModeForce_(IOobject("appliedModeForce",
          io.time().constant() / fileName("fsiRestart"), fdmMesh_,
          IOobject::NO_READ, IOobject::NO_WRITE)),
      fsiResidual_(1.0), theta_(1.4), // Default Wilson-Theta
      mappingTolerance_(4e-6), couplingRelaxation_(1.0),
      pressureFieldName_("p"), rhoRef_(-1.0), pRef_(0.0), maxDispChange_(1.0),
      couplingAcceleration_("constant"), iqnMaxHistory_(10),
      iqnFilter_("QR2"), iqnFilterLimit_(1e-2), imvjReuse_(5),
      iqnDiagnostics_(false), aitkenPrevResidual_(0), aitkenAlpha_(1.0),
      aitkenHasPrev_(false), iqnV_(0), iqnW_(0), iqnPrevResidual_(0),
      iqnPrevRaw_(0), iqnHasPrev_(false), accelIterationStarted_(false),
      accelHistoryLive_(false), imvjChunkWtilde_(0), imvjChunkZ_(0),
      imvjChunkCols_(0),
      startupStepCount_(0), lastUpdateTimeIndex_(-1), updateCount_(0),
      writeDiagnosticsEnabled_(false), trackTimingEnabled_(false),
      timingLastUpdateCpuTime_(0.0), timingHasLastUpdate_(false),
      timingFluidCpuAccum_(0.0), timingMeshCpuAccum_(0.0),
      timingFluidCpuAtLastWrite_(0.0), timingMeshCpuAtLastWrite_(0.0),
      lastGlobalWriteTimeIndex_(-1), fsiPatchIDs_(0), fsiPolyPatches_(0),
      faceModeProjection_(0), pressureScaleCache_(0),
      pressureScaleInitialized_(0), diagnosticsHeaderWritten_(false),
      diagStream_(nullptr),
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
      structuralForceTimeIndex_(-1),
      meshRefinementSupport_(false), refinementInterpTolerance_(1e-8),
      referenceFsiBuilt_(false), referenceFaceModeProjection_(0),
      fsiFaceToReferenceFaces_(0), fsiFaceToReferenceWeights_(0),
      referenceFsiFaceAreas_(0), aggregatedPressureRefs_(0),
      aggregatedShearRefs_(0)
{}

const word& FdmCore::fdmTypeName()
{
    static const word name("fastDynamicFvMesh");
    return name;
}

const labelList& FdmCore::amrPointLevel() const
{
    static const labelList empty;
    return empty;
}

void FdmCore::initialiseCoreState()
{
    // readModeShapes() reads the mode header/files first, resolves the active
    // mode count, and only then loads the AMR records and FSI reference cache.
    readModeShapes();

    const bool restartStateRead = readRestartStateStore();
    const bool stateRead = restartStateRead && (modeState_.size() > 0);

    if (stateRead && modeState_.size() == nMode_)
    {
        Info << "Restarting fastDynamicFvMesh with " << nMode_
             << " modes from time " << fdmMesh_.time().timeName() << endl;

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
}

FdmCore::~FdmCore() {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool FdmCore::updateFdm()
{
    const scalar updateStartCpuTime = fdmMesh_.time().elapsedCpuTime();
    accountFluidCpuAtUpdateStart(updateStartCpuTime);

    const label currentTimeIndex = fdmMesh_.time().timeIndex();
    const bool firstIter =
        beginUpdateStep(currentTimeIndex, isNewTimeStep(currentTimeIndex));

    const scalar dt = fdmMesh_.time().deltaTValue();
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
