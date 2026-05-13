/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Fast Dynamic Mesh update-stage implementation.

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "fastDynamicFvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void fastDynamicFvMesh::accountFluidCpuAtUpdateStart(
    const scalar updateStartCpuTime)
{
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
}

bool fastDynamicFvMesh::beginUpdateStep(const label currentTimeIndex)
{
    const bool firstIter = (currentTimeIndex != lastUpdateTimeIndex_);

    if (runtimeRefinementEnabled_ && firstIter)
    {
        if (!trackAmrTimingEnabled_)
        {
            dynamicRefineFvMesh::updateTopology();
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

    lastUpdateTimeIndex_ = currentTimeIndex;

    if (firstIter)
    {
        updateCount_ = 0;
        // Store previous state only on the first iteration of the time step
        modeForce0_ = modeForce_;
        modeState0_ = modeState_;
    }
    updateCount_++;

    return firstIter;
}

void fastDynamicFvMesh::assembleAndRelaxModalForces(const bool firstIter)
{
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

    publishFsiResidual();
}

void fastDynamicFvMesh::publishFsiResidual()
{
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
}

bool fastDynamicFvMesh::initializeStartupStateIfNeeded()
{
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

        return true;
    }

    return false;
}

void fastDynamicFvMesh::solveAndMoveMesh(const scalar dt)
{
    if (nMode_ > 0)
    {
        solveStructuralDynamics(dt);
    }

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
}

void fastDynamicFvMesh::logStructuralLoadSummary()
{
    if (Pstream::master() && structuralForceEnabled_ && updateCount_ == 1 &&
        nMode_ > 0)
    {
        Info << "Structural modal load summary at t=" << this->time().timeName()
             << ": mode0=" << structuralModeForce_[0]
             << ", signal=" << structuralForceSignal_ << endl;
    }
}

void fastDynamicFvMesh::accountMeshCpuAtUpdateEnd(
    const scalar updateStartCpuTime)
{
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
}

void fastDynamicFvMesh::writeUpdateOutputs(const label currentTimeIndex)
{
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
}

bool fastDynamicFvMesh::updateReturnState()
{
    return runtimeRefinementEnabled_ || this->moving() || this->topoChanging();
}

} // End namespace Foam
