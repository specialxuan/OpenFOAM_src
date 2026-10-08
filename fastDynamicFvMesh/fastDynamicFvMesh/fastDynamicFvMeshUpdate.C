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
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace
{
    //- Dot product of two scalar fields (file-local accelerator helper)
    scalar fieldDot(const scalarField& a, const scalarField& b)
    {
        scalar s = 0.0;
        forAll(a, i)
        {
            s += a[i]*b[i];
        }
        return s;
    }

    //- Push one IQN history column to the front (newest first) and cap the
    //  list by dropping the oldest (back) column
    void prependIqnColumn
    (
        List<scalarField>& columns,
        const scalarField& column,
        const label maxColumns
    )
    {
        const label n = columns.size();
        label newSize = n + 1;
        if (newSize > maxColumns)
        {
            newSize = maxColumns;
        }

        List<scalarField> resized(newSize);
        resized[0] = column;
        for (label j = 0; j < newSize - 1 && j < n; ++j)
        {
            resized[j + 1] = columns[j];
        }
        columns.transfer(resized);
    }

    //- Filtered modified Gram-Schmidt QR of the newest-first IQN history.
    //  Column selection follows the frozen none/QR2/exponential rules; Q is
    //  built over the kept columns only, W keeps the column correspondence and
    //  rRows forms the upper-triangular R (rRows[a][b], b >= a).
    //  Returns the number of kept columns (0 => no usable history column).
    label filterIqnColumns
    (
        const List<scalarField>& V,
        const List<scalarField>& W,
        const word& filter,
        const scalar filterLimit,
        List<scalarField>& vCols,
        List<scalarField>& qCols,
        List<scalarField>& rRows,
        List<scalarField>& wCols
    )
    {
        vCols.clear();
        qCols.clear();
        rRows.clear();
        wCols.clear();

        label mUsed = V.size();
        if (filter == "exponential")
        {
            label nCut = label(std::ceil(-std::log(filterLimit)));
            if (nCut < 1)
            {
                nCut = 1;
            }
            if (mUsed > nCut)
            {
                mUsed = nCut;
            }
        }

        // Select columns (newest first), orthogonalising against kept columns
        labelList kept;
        List<scalarField> qSel;
        for (label j = 0; j < mUsed; ++j)
        {
            const scalarField& v = V[j];
            scalarField q(v);
            forAll(qSel, i)
            {
                q -= fieldDot(qSel[i], q)*qSel[i];
            }

            const scalar rho0 = std::sqrt(fieldDot(v, v));
            const scalar rhoOrth = std::sqrt(fieldDot(q, q));

            bool drop = (rhoOrth <= VSMALL);
            if
            (
                !drop
             && filter == "QR2"
             && j > 0
             && rho0 > 0
             && rhoOrth < filterLimit*rho0
            )
            {
                drop = true;
            }

            if (!drop)
            {
                kept.append(j);
                qSel.append(scalarField(q/rhoOrth));
            }
        }

        // Build exact Q/R over the kept columns
        const label mKept = kept.size();
        vCols.setSize(mKept);
        qCols.setSize(mKept);
        wCols.setSize(mKept);
        rRows.setSize(mKept);
        forAll(rRows, a)
        {
            rRows[a] = scalarField(mKept, 0.0);
        }

        for (label a = 0; a < mKept; ++a)
        {
            const scalarField& v = V[kept[a]];
            scalarField q(v);
            for (label b = 0; b < a; ++b)
            {
                const scalar rba = fieldDot(qCols[b], q);
                rRows[b][a] = rba;
                q -= rba*qCols[b];
            }
            rRows[a][a] = std::sqrt(fieldDot(q, q));
            qCols[a] = q/rRows[a][a];
            vCols[a] = v;
            wCols[a] = W[kept[a]];
        }

        return mKept;
    }

    //- Accumulate J_chunks x into out, with J_chunks = sum_q Wtilde_q Z_q.
    //  Wtilde_q is nMode * m_q and Z_q is m_q * nMode, both row-major flat.
    void addImvjChunkMatvec
    (
        scalarField& out,
        const List<scalarField>& chunkWtilde,
        const List<scalarField>& chunkZ,
        const List<label>& chunkCols,
        const scalarField& x
    )
    {
        const label nMode = x.size();

        forAll(chunkWtilde, q)
        {
            const scalarField& wTilde = chunkWtilde[q];
            const scalarField& z = chunkZ[q];
            const label m = chunkCols[q];

            scalarField t(m, 0.0);
            for (label j = 0; j < m; ++j)
            {
                scalar s = 0.0;
                for (label i = 0; i < nMode; ++i)
                {
                    s += z[j*nMode + i]*x[i];
                }
                t[j] = s;
            }

            for (label i = 0; i < nMode; ++i)
            {
                scalar s = 0.0;
                for (label j = 0; j < m; ++j)
                {
                    s += wTilde[i*m + j]*t[j];
                }
                out[i] += s;
            }
        }
    }

    //- Materialise Z = R^-1 Q^T (m x nMode, row-major flat) from the filtered
    //  QR factors of a window's V
    void imvjPseudoInverse
    (
        const List<scalarField>& qCols,
        const List<scalarField>& rRows,
        scalarField& zFlat
    )
    {
        const label m = qCols.size();
        const label nMode = qCols[0].size();

        // R^-1 (rows rInv[a][c])
        List<scalarField> rInv(m, scalarField(m, 0.0));
        for (label c = 0; c < m; ++c)
        {
            for (label a = m - 1; a >= 0; --a)
            {
                scalar rhs = (a == c ? 1.0 : 0.0);
                for (label b = a + 1; b < m; ++b)
                {
                    rhs -= rRows[a][b]*rInv[b][c];
                }
                rInv[a][c] = rhs/rRows[a][a];
            }
        }

        zFlat.setSize(m*nMode);
        for (label j = 0; j < m; ++j)
        {
            for (label i = 0; i < nMode; ++i)
            {
                scalar s = 0.0;
                for (label a = 0; a < m; ++a)
                {
                    s += rInv[j][a]*qCols[a][i];
                }
                zFlat[j*nMode + i] = s;
            }
        }
    }

    //- Drop the oldest (front) retained IMVJ chunk
    void dropOldestImvjChunk
    (
        List<scalarField>& chunkWtilde,
        List<scalarField>& chunkZ,
        List<label>& chunkCols
    )
    {
        const label n = chunkCols.size();

        List<scalarField> newWtilde(n - 1);
        List<scalarField> newZ(n - 1);
        List<label> newCols(n - 1);
        for (label q = 1; q < n; ++q)
        {
            newWtilde[q - 1] = chunkWtilde[q];
            newZ[q - 1] = chunkZ[q];
            newCols[q - 1] = chunkCols[q];
        }

        chunkWtilde.transfer(newWtilde);
        chunkZ.transfer(newZ);
        chunkCols.transfer(newCols);
    }
}

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

    lastUpdateTimeIndex_ = currentTimeIndex;

    if (firstIter)
    {
        updateCount_ = 0;
        // Store previous state only on the first iteration of the time step
        modeForce0_ = modeForce_;
        modeState0_ = modeState_;
        resetCouplingAccelerationState(true);
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

    // Produce the force to apply this iteration. `constant` reproduces the
    // historical per-mode loop bit-for-bit (including the override).
    const scalarField nextAppliedForce =
        applyCouplingAcceleration(modeForce_, firstIter);

    // Residual relative to the NEW force magnitude. Convergence is when
    // F_fluid matches F_old (correction goes to zero) or when F_new matches
    // F_old (change is zero).
    fsiResidual_ = 0.0;
    for (label i = 0; i < nMode_; ++i)
    {
        const scalar prevRelaxedForce = appliedModeForce_[i];
        const scalar newRelaxedForce  = nextAppliedForce[i];
        const scalar diff     = mag(newRelaxedForce - prevRelaxedForce);
        const scalar magForce = mag(newRelaxedForce) + 1e-12; // Avoid div by 0
        const scalar res      = diff / magForce;
        if (res > fsiResidual_)
        {
            fsiResidual_ = res;
        }
    }

    // Update state
    appliedModeForce_ = nextAppliedForce;
    modeForce_ = nextAppliedForce; // Use relaxed force for structural dynamics

    publishFsiResidual();
}

scalarField fastDynamicFvMesh::applyCouplingAcceleration
(
    const scalarField& rawFluidForce,
    const bool firstIter
)
{
    if (couplingAcceleration_ == "constant")
    {
        scalarField next(nMode_);
        for (label i = 0; i < nMode_; ++i)
        {
            const scalar raw = rawFluidForce[i];
            const scalar prevRelaxedForce = appliedModeForce_[i];
            scalar newRelaxedForce =
                prevRelaxedForce
              + couplingRelaxation_*(raw - prevRelaxedForce);

            // If this is the first iteration of the time step, and we are
            // RESTARTING (or startup), trust the computed force more if the
            // previous force is suspiciously zero while the new force is
            // large.
            if (firstIter && mag(prevRelaxedForce) < VSMALL && mag(raw) > 1.0)
            {
                // Likely a bad initialization on restart or startup
                newRelaxedForce = raw;

                if (Pstream::master() && i == 0)
                {
                    Info << "  First iter override: using raw force "
                         << raw << " instead of relaxed "
                         << newRelaxedForce << endl;
                }
            }
            next[i] = newRelaxedForce;
        }
        return next;
    }

    // Non-constant accelerators share the startup gate and the initial
    // relaxation. History is live only after the two startup pseudo-updates
    // (initializeStartupStateIfNeeded skips the structural solve there).
    accelHistoryLive_ = (startupStepCount_ >= 2);

    scalarField next(nMode_);
    for (label i = 0; i < nMode_; ++i)
    {
        next[i] = appliedModeForce_[i]
            + couplingRelaxation_*(rawFluidForce[i] - appliedModeForce_[i]);
    }

    if (!accelHistoryLive_)
    {
        return next;
    }

    // Internal fixed-point residual this iteration
    scalarField r(nMode_);
    for (label i = 0; i < nMode_; ++i)
    {
        r[i] = rawFluidForce[i] - appliedModeForce_[i];
    }

    if (!accelIterationStarted_)
    {
        // First window iteration: initial relaxation, seed history.
        accelIterationStarted_ = true;
        iqnPrevResidual_ = r;
        iqnPrevRaw_ = rawFluidForce;
        iqnHasPrev_ = true;
        aitkenPrevResidual_ = r;
        aitkenAlpha_ = couplingRelaxation_;
        aitkenHasPrev_ = true;

        if (iqnDiagnostics_ && Pstream::master())
        {
            Info << "  CouplingAcceleration " << couplingAcceleration_
                 << ": alpha=" << couplingRelaxation_
                 << ", ||r||=" << std::sqrt(fieldDot(r, r)) << endl;
        }

        return next;
    }

    auto resetIqnWindow = [this]()
    {
        iqnV_.clear();
        iqnW_.clear();
        iqnHasPrev_ = false;
        iqnPrevResidual_.clear();
        iqnPrevRaw_.clear();
        accelIterationStarted_ = false;
    };

    if (couplingAcceleration_ == "aitken")
    {
        scalar alpha = couplingRelaxation_;

        if (aitkenHasPrev_)
        {
            // Kuttler-Wall scalar Aitken:
            //   alpha_k = -alpha_{k-1} (r_{k-1}.(r_k - r_{k-1}))
            //             / ||r_k - r_{k-1}||^2
            const scalarField delta(r - aitkenPrevResidual_);
            const scalar den = fieldDot(delta, delta);
            const scalar dot = fieldDot(aitkenPrevResidual_, delta);

            if (den > 1e-30)
            {
                alpha = -aitkenAlpha_*dot/den;

                if (!std::isfinite(alpha) || alpha <= 0.0)
                {
                    alpha = couplingRelaxation_;
                }

                // Upper clamp only: the added-mass case needs small alpha.
                if (alpha > 1.0)
                {
                    alpha = 1.0;
                }
            }
            else
            {
                // Degenerate residual change: keep the previous relaxation
                alpha = aitkenAlpha_;
            }
        }

        for (label i = 0; i < nMode_; ++i)
        {
            next[i] = appliedModeForce_[i] + alpha*r[i];
        }

        aitkenPrevResidual_ = r;
        aitkenAlpha_ = alpha;

        if (iqnDiagnostics_ && Pstream::master())
        {
            Info << "  CouplingAcceleration aitken: alpha=" << alpha
                 << ", ||r||=" << std::sqrt(fieldDot(r, r)) << endl;
        }

        return next;
    }

    if (couplingAcceleration_ == "IQN-ILS")
    {
        // Append the current secant column at the front (newest first). The
        // first history-enabled iteration seeded r_1/raw_1, so iteration 2 is
        // the first QN step (one column).
        if (iqnHasPrev_)
        {
            const scalarField deltaR(r - iqnPrevResidual_);
            const scalarField deltaRaw(rawFluidForce - iqnPrevRaw_);
            prependIqnColumn(iqnV_, deltaR, iqnMaxHistory_);
            prependIqnColumn(iqnW_, deltaRaw, iqnMaxHistory_);
        }
        iqnPrevResidual_ = r;
        iqnPrevRaw_ = rawFluidForce;
        iqnHasPrev_ = true;

        List<scalarField> vCols;
        List<scalarField> qCols;
        List<scalarField> rRows;
        List<scalarField> wCols;
        const label mKept = filterIqnColumns
        (
            iqnV_, iqnW_, iqnFilter_, iqnFilterLimit_,
            vCols, qCols, rRows, wCols
        );

        if (mKept <= 0)
        {
            // No informative column this iteration: keep initial relaxation
            return next;
        }

        // Least-squares solve V c = -r via the reduced QR (R c = Q^T(-r))
        scalarField c(mKept, 0.0);
        for (label a = mKept - 1; a >= 0; --a)
        {
            scalar rhs = -fieldDot(qCols[a], r);
            for (label b = a + 1; b < mKept; ++b)
            {
                rhs -= rRows[a][b]*c[b];
            }
            c[a] = rhs/rRows[a][a];

            if (!std::isfinite(c[a]))
            {
                resetIqnWindow();
                return next;
            }
        }

        // QN update F_{k+1} = raw_k + W c
        scalarField qnNext(nMode_, 0.0);
        for (label i = 0; i < nMode_; ++i)
        {
            scalar wc = 0.0;
            for (label a = 0; a < mKept; ++a)
            {
                wc += c[a]*wCols[a][i];
            }
            qnNext[i] = rawFluidForce[i] + wc;
        }

        for (label i = 0; i < nMode_; ++i)
        {
            if (!std::isfinite(qnNext[i]))
            {
                resetIqnWindow();
                return next;
            }
        }

        next = qnNext;

        if (iqnDiagnostics_ && Pstream::master())
        {
            Info << "  CouplingAcceleration IQN-ILS: activeColumns=" << mKept
                 << ", history=" << iqnV_.size()
                 << ", ||r||=" << std::sqrt(fieldDot(r, r)) << endl;
        }

        return next;
    }

    if (couplingAcceleration_ == "IQN-IMVJ")
    {
        // Same current-window history append as IQN-ILS
        if (iqnHasPrev_)
        {
            const scalarField deltaR(r - iqnPrevResidual_);
            const scalarField deltaRaw(rawFluidForce - iqnPrevRaw_);
            prependIqnColumn(iqnV_, deltaR, iqnMaxHistory_);
            prependIqnColumn(iqnW_, deltaRaw, iqnMaxHistory_);
        }
        iqnPrevResidual_ = r;
        iqnPrevRaw_ = rawFluidForce;
        iqnHasPrev_ = true;

        List<scalarField> vCols;
        List<scalarField> qCols;
        List<scalarField> rRows;
        List<scalarField> wCols;
        const label mKept = filterIqnColumns
        (
            iqnV_, iqnW_, iqnFilter_, iqnFilterLimit_,
            vCols, qCols, rRows, wCols
        );

        // Newton update F_{k+1} = raw_k + delta,
        //   delta = J_prev(-r_k) + Wtilde_cur Z_cur (-r_k)
        scalarField y(nMode_, 0.0);
        for (label i = 0; i < nMode_; ++i)
        {
            y[i] = -r[i];
        }

        scalarField delta(nMode_, 0.0);
        addImvjChunkMatvec
        (
            delta, imvjChunkWtilde_, imvjChunkZ_, imvjChunkCols_, y
        );

        if (mKept > 0)
        {
            scalarField zFlat;
            imvjPseudoInverse(qCols, rRows, zFlat);

            bool finiteZ = true;
            forAll(zFlat, k)
            {
                if (!std::isfinite(zFlat[k]))
                {
                    finiteZ = false;
                    break;
                }
            }
            if (!finiteZ)
            {
                resetIqnWindow();
                return next;
            }

            // t = Z_cur (-r_k)
            scalarField t(mKept, 0.0);
            for (label j = 0; j < mKept; ++j)
            {
                scalar s = 0.0;
                for (label i = 0; i < nMode_; ++i)
                {
                    s += zFlat[j*nMode_ + i]*y[i];
                }
                t[j] = s;

                if (!std::isfinite(t[j]))
                {
                    resetIqnWindow();
                    return next;
                }
            }

            // delta += (W_cur - J_prev V_cur) t
            for (label j = 0; j < mKept; ++j)
            {
                scalarField jV(nMode_, 0.0);
                addImvjChunkMatvec
                (
                    jV, imvjChunkWtilde_, imvjChunkZ_, imvjChunkCols_,
                    vCols[j]
                );
                for (label i = 0; i < nMode_; ++i)
                {
                    delta[i] += (wCols[j][i] - jV[i])*t[j];
                }
            }
        }
        else if (imvjChunkCols_.size() == 0)
        {
            // Degenerate first window without a usable column: keep the
            // initial relaxation (same fallback as IQN-ILS)
            return next;
        }

        for (label i = 0; i < nMode_; ++i)
        {
            if (!std::isfinite(delta[i]))
            {
                resetIqnWindow();
                return next;
            }
        }

        for (label i = 0; i < nMode_; ++i)
        {
            next[i] = rawFluidForce[i] + delta[i];
        }

        if (iqnDiagnostics_ && Pstream::master())
        {
            Info << "  CouplingAcceleration IQN-IMVJ: activeColumns=" << mKept
                 << ", history=" << iqnV_.size()
                 << ", chunks=" << imvjChunkCols_.size()
                 << ", ||r||=" << std::sqrt(fieldDot(r, r)) << endl;
        }

        return next;
    }

    // All validated couplingAcceleration tokens are dispatched above
    FatalErrorInFunction
        << "Unhandled couplingAcceleration '" << couplingAcceleration_ << "'"
        << exit(FatalError);

    return next;
}

void fastDynamicFvMesh::resetCouplingAccelerationState(const bool newWindow)
{
    // Per-window state is dropped at the first update() of each time step.
    // IQN-IMVJ first commits the finished window's (Wtilde, Z) chunk, computed
    // with the then-current J_prev (the retained chunk list); imvjChunk*_ are
    // retained across windows as a sliding list capped at imvjReuse_.
    if
    (
        newWindow
     && couplingAcceleration_ == "IQN-IMVJ"
     && iqnV_.size() >= 1
    )
    {
        List<scalarField> vCols;
        List<scalarField> qCols;
        List<scalarField> rRows;
        List<scalarField> wCols;
        const label mKept = filterIqnColumns
        (
            iqnV_, iqnW_, iqnFilter_, iqnFilterLimit_,
            vCols, qCols, rRows, wCols
        );

        if (mKept > 0)
        {
            const label nMode = vCols[0].size();

            scalarField zFlat;
            imvjPseudoInverse(qCols, rRows, zFlat);

            // Wtilde = W - J_prev V
            scalarField wTilde(nMode*mKept, 0.0);
            for (label j = 0; j < mKept; ++j)
            {
                scalarField jV(nMode, 0.0);
                addImvjChunkMatvec
                (
                    jV, imvjChunkWtilde_, imvjChunkZ_, imvjChunkCols_,
                    vCols[j]
                );
                for (label i = 0; i < nMode; ++i)
                {
                    wTilde[i*mKept + j] = wCols[j][i] - jV[i];
                }
            }

            bool finite = true;
            forAll(zFlat, k)
            {
                if (!std::isfinite(zFlat[k]))
                {
                    finite = false;
                    break;
                }
            }
            forAll(wTilde, k)
            {
                if (!std::isfinite(wTilde[k]))
                {
                    finite = false;
                    break;
                }
            }

            if (finite)
            {
                imvjChunkWtilde_.append(wTilde);
                imvjChunkZ_.append(zFlat);
                imvjChunkCols_.append(mKept);

                while (imvjChunkCols_.size() > imvjReuse_)
                {
                    dropOldestImvjChunk
                    (
                        imvjChunkWtilde_, imvjChunkZ_, imvjChunkCols_
                    );
                }
            }
        }
    }

    aitkenHasPrev_ = false;
    aitkenPrevResidual_.clear();
    aitkenAlpha_ = 1.0;
    iqnV_.clear();
    iqnW_.clear();
    iqnHasPrev_ = false;
    iqnPrevResidual_.clear();
    iqnPrevRaw_.clear();
    accelIterationStarted_ = false;
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

        writeModalExtensionOutputs();
    }
}

bool fastDynamicFvMesh::updateReturnState()
{
    return runtimeRefinementEnabled_ || this->moving() || this->topoChanging();
}

} // End namespace Foam
