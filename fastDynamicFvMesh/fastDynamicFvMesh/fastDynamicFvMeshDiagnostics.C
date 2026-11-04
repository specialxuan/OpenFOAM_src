/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    fastDynamicFvMesh diagnostics and timing-report implementation.

\*---------------------------------------------------------------------------*/

#include "fastDynamicFvMesh.H"
#include "Pstream.H"
#include <cstdio>
#include <fstream>
#include <iomanip>

namespace Foam
{

// writeDiagnostics()
// Purpose: write per-time-step modal diagnostics to modal_diagnostics.csv on
// the master (or the run directory in parallel runs). Details:
//  - If the diagnostics file does not exist, write a header with columns:
//      Time, Force_i, PressureForce_i, ShearForce_i, Disp_i, Vel_i, Acc_i,
//      AppliedDisp_i for each mode.
//  - Append the current time value and the arrays: modeForce_,
//  modePressureForce_, modeShearForce_, modeState_ components and
//  appliedModeDisp_.
//  - Use high precision formatting to preserve numeric fidelity for
//  post-processing comparisons.
void fastDynamicFvMesh::writeDiagnostics()
{
    if (!writeDiagnosticsEnabled_)
    {
        return;
    }

    if (!Pstream::master())
    {
        return;
    }

    // Time + 10 groups per mode: Force, PressureForce, ShearForce,
    // StructuralForce, Disp, Vel, Acc, AppliedDisp, StructuralScale, RawForce.
    const label expectedColumns = 1 + 10*nMode_;

    if (!diagStream_.valid())
    {
        fileName diagnosticsPath = this->time().path() / "modal_diagnostics.csv";

        if (Pstream::parRun())
        {
            diagnosticsPath =
                this->time().path().path() / "modal_diagnostics.csv";
        }

        // Restart/continuation: inspect the existing header.  A non-empty
        // compatible file keeps its header; an old/incompatible layout is
        // renamed so the new columns cannot be mixed into it.
        bool hasExistingHeader = false;
        label existingColumns = 0;

        {
            std::ifstream existing(diagnosticsPath.c_str());

            if (existing.good())
            {
                std::string header;
                if (std::getline(existing, header) && !header.empty())
                {
                    hasExistingHeader = true;
                    existingColumns = 1;
                    for (const char c : header)
                    {
                        if (c == ',')
                        {
                            ++existingColumns;
                        }
                    }
                }
            }
        }

        if (hasExistingHeader && existingColumns != expectedColumns)
        {
            const fileName backup =
                diagnosticsPath + ".bak" + this->time().timeName();

            if (std::rename(diagnosticsPath.c_str(), backup.c_str()) != 0)
            {
                FatalErrorInFunction
                    << "Cannot rename incompatible diagnostics file "
                    << diagnosticsPath << " to " << backup << exit(FatalError);
            }

            if (Pstream::master())
            {
                WarningInFunction
                    << "Existing " << diagnosticsPath << " has "
                    << existingColumns << " columns, expected "
                    << expectedColumns << "; renamed to " << backup
                    << " and started a new file." << endl;
            }
        }
        else if (hasExistingHeader)
        {
            diagnosticsHeaderWritten_ = true;
        }

        diagStream_.reset(new std::ofstream(diagnosticsPath.c_str(),
            std::ios::app));

        if (!diagStream_->good())
        {
            WarningInFunction << "Unable to open " << diagnosticsPath
                              << " for writing." << endl;
            diagStream_.clear();
            return;
        }
    }

    std::ofstream& diagFile = *diagStream_;

    if (!diagnosticsHeaderWritten_)
    {
        diagFile << "Time";

        for (label i = 0; i < nMode_; ++i)
        {
            diagFile << ",Force_" << (i + 1);
        }

        for (label i = 0; i < nMode_; ++i)
        {
            diagFile << ",PressureForce_" << (i + 1);
        }

        for (label i = 0; i < nMode_; ++i)
        {
            diagFile << ",ShearForce_" << (i + 1);
        }

        for (label i = 0; i < nMode_; ++i)
        {
            diagFile << ",StructuralForce_" << (i + 1);
        }

        for (label i = 0; i < nMode_; ++i)
        {
            diagFile << ",Disp_" << (i + 1) << ",Vel_" << (i + 1) << ",Acc_"
                     << (i + 1) << ",AppliedDisp_" << (i + 1)
                     << ",StructuralScale_" << (i + 1);
        }

        // Appended at the end so existing fixed-column readers keep their
        // original indices.  RawForce is the unrelaxed fluid load
        // (pressure + shear), i.e. the correct signal for GAF post-processing.
        for (label i = 0; i < nMode_; ++i)
        {
            diagFile << ",RawForce_" << (i + 1);
        }

        diagFile << "\n";
        diagnosticsHeaderWritten_ = true;
    }

    diagFile << std::setprecision(12);
    diagFile << this->time().value();

    for (label i = 0; i < nMode_; ++i)
    {
        diagFile << "," << modeForce_[i];
    }

    for (label i = 0; i < nMode_; ++i)
    {
        diagFile << "," << modePressureForce_[i];
    }

    for (label i = 0; i < nMode_; ++i)
    {
        diagFile << "," << modeShearForce_[i];
    }

    for (label i = 0; i < nMode_; ++i)
    {
        diagFile << "," << structuralModeForce_[i];
    }

    const scalar structuralScale =
        structuralForceEnabled_ ? structuralForceSignal_ : 0.0;

    for (label i = 0; i < nMode_; ++i)
    {
        diagFile << "," << modeState_[i].x() << "," << modeState_[i].y() << ","
                 << modeState_[i].z() << "," << appliedModeDisp_[i] << ","
                 << structuralScale;
    }

    for (label i = 0; i < nMode_; ++i)
    {
        diagFile << "," << modePressureForce_[i] + modeShearForce_[i];
    }

    diagFile << "\n";
    diagFile.flush();
}

void fastDynamicFvMesh::writeTimingReport()
{
    if (!trackTimingEnabled_ || !Pstream::master())
    {
        return;
    }

    const scalar totalFluid = timingFluidCpuAccum_;
    const scalar totalMesh = timingMeshCpuAccum_;
    const scalar totalAll = totalFluid + totalMesh;

    const scalar deltaFluid = totalFluid - timingFluidCpuAtLastWrite_;
    const scalar deltaMesh = totalMesh - timingMeshCpuAtLastWrite_;
    const scalar deltaAll = deltaFluid + deltaMesh;

    const scalar totalFluidPct =
        totalAll > VSMALL ? 100.0 * totalFluid / totalAll : 0.0;
    const scalar totalMeshPct =
        totalAll > VSMALL ? 100.0 * totalMesh / totalAll : 0.0;
    const scalar deltaFluidPct =
        deltaAll > VSMALL ? 100.0 * deltaFluid / deltaAll : 0.0;
    const scalar deltaMeshPct =
        deltaAll > VSMALL ? 100.0 * deltaMesh / deltaAll : 0.0;

    Info << "FSI timing [CPU s] at t=" << this->time().timeName()
         << "  step{fluid=" << deltaFluid << " (" << deltaFluidPct << "%)"
         << ", mesh=" << deltaMesh << " (" << deltaMeshPct << "%)"
         << ", total=" << deltaAll << "}" << "  cumulative{fluid=" << totalFluid
         << " (" << totalFluidPct << "%)" << ", mesh=" << totalMesh << " ("
         << totalMeshPct << "%)" << ", total=" << totalAll << "}" << endl;

    timingFluidCpuAtLastWrite_ = totalFluid;
    timingMeshCpuAtLastWrite_ = totalMesh;
}

} // namespace Foam
