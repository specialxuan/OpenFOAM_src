/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Fast Dynamic Mesh modal-extension (structural) implementation.

    Functions moved verbatim from fastDynamicFvMeshForces.C,
    fastDynamicFvMeshIO.C and fastDynamicFvMeshUpdate.C; logic, ordering and
    message text are unchanged.

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "SortableList.H"
#include "fastDynamicFvMesh.H"
#include "Time.H"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace
{

void normalizeCsvLine(string& line)
{
    std::replace(line.begin(), line.end(), ',', ' ');
    std::replace(line.begin(), line.end(), '\t', ' ');
}

bool containsNumericToken(const string& line)
{
    for (const char c : line)
    {
        if ((c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.')
        {
            return true;
        }
    }

    return false;
}

bool isIgnorableCsvLine(const string& line)
{
    for (const char c : line)
    {
        if (c == ' ' || c == '\t' || c == '\r')
        {
            continue;
        }

        return c == '#';
    }

    return true;
}

bool readNextCsvLine(std::ifstream& is, string& line, label& lineNo)
{
    if (std::getline(is, line))
    {
        ++lineNo;
        return true;
    }

    return false;
}

bool readNextCsvDataLine(std::ifstream& is, string& line, label& lineNo)
{
    while (readNextCsvLine(is, line, lineNo))
    {
        if (!isIgnorableCsvLine(line) && containsNumericToken(line))
        {
            return true;
        }
    }

    return false;
}

template<class... Values>
bool parseCsvValues(const string& input, Values&... values)
{
    string line(input);
    normalizeCsvLine(line);
    std::stringstream ss(line);

    bool ok = true;
    using expander = int[];
    (void)expander{0, ((ok = ok && bool(ss >> values)), 0)...};

    return ok;
}

//- Read one 6-component structural mode quantity (stress or strain) for all
//  modes from <prefix><m>.csv into a flat [nodeI*6 + component] layout.
//  Header is (frequency, [nStructNodes, nMode]) like the shape files; each
//  body row must provide six values in the given component order.
void readSixComponentStructuralModeFiles(
    const fileName& modeDir,
    const word& prefix,
    const label nStructNodes,
    const label nMode,
    const word& referenceFileName,
    const char* quantityName,
    const char* componentOrder,
    List<scalarField>& modeValues)
{
    modeValues.setSize(nMode);

    for (label m = 0; m < nMode; ++m)
    {
        const fileName quantityPath =
            modeDir / (prefix + std::to_string(m + 1) + ".csv");
        std::ifstream in(quantityPath);

        if (!in.good())
        {
            FatalErrorInFunction
                << "Cannot open required structural " << quantityName
                << " file " << quantityPath << exit(FatalError);
        }

        string line;
        label lineNo = 0;

        if (!readNextCsvDataLine(in, line, lineNo))
        {
            FatalErrorInFunction
                << "Missing frequency line in structural " << quantityName
                << " file " << quantityPath << exit(FatalError);
        }

        scalar freq = 0.0;
        scalar fileNodeCount = 0.0;
        scalar fileModeCount = 0.0;

        string headerLine(line);
        normalizeCsvLine(headerLine);
        std::stringstream hLine(headerLine);
        if (!(hLine >> freq))
        {
            FatalErrorInFunction
                << "Failed to read frequency from structural " << quantityName
                << " file " << quantityPath << " at line " << lineNo << ": "
                << line << exit(FatalError);
        }

        if ((hLine >> fileNodeCount) && (hLine >> fileModeCount))
        {
            if (label(fileNodeCount) != nStructNodes ||
                label(fileModeCount) != nMode)
            {
                FatalErrorInFunction
                    << "Structural " << quantityName << " file "
                    << quantityPath << " reports " << label(fileNodeCount)
                    << " nodes and " << label(fileModeCount) << " modes, but "
                    << referenceFileName << " reports " << nStructNodes
                    << " nodes and " << nMode << " modes."
                    << exit(FatalError);
            }
        }

        scalarField& nodeValues = modeValues[m];
        nodeValues.setSize(nStructNodes*6, 0.0);

        for (label nodeI = 0; nodeI < nStructNodes; ++nodeI)
        {
            if (!readNextCsvDataLine(in, line, lineNo))
            {
                FatalErrorInFunction
                    << "Unexpected end of file while reading structural "
                    << quantityName << " node " << nodeI << " from "
                    << quantityPath << exit(FatalError);
            }

            scalar s0 = 0.0;
            scalar s1 = 0.0;
            scalar s2 = 0.0;
            scalar s3 = 0.0;
            scalar s4 = 0.0;
            scalar s5 = 0.0;

            if (!parseCsvValues(line, s0, s1, s2, s3, s4, s5))
            {
                FatalErrorInFunction
                    << "Invalid structural " << quantityName
                    << " entry for node " << nodeI << " in " << quantityPath
                    << " at line " << lineNo << ": " << line
                    << ". Expected six values (" << componentOrder << ")."
                    << exit(FatalError);
            }

            nodeValues[nodeI*6 + 0] = s0;
            nodeValues[nodeI*6 + 1] = s1;
            nodeValues[nodeI*6 + 2] = s2;
            nodeValues[nodeI*6 + 3] = s3;
            nodeValues[nodeI*6 + 4] = s4;
            nodeValues[nodeI*6 + 5] = s5;
        }
    }
}

} // End anonymous namespace

vector fastDynamicFvMesh::interpolateStructuralForce(
    const scalarField& times, const vectorField& values, const scalar t) const
{
    if (times.size() == 1)
    {
        return values[0];
    }

    if (t <= times.first())
    {
        return values.first();
    }

    if (t >= times.last())
    {
        return values.last();
    }

    for (label i = 0; i < times.size() - 1; ++i)
    {
        const scalar t0 = times[i];
        const scalar t1 = times[i + 1];

        if (t >= t0 && t <= t1)
        {
            const scalar dt = t1 - t0;
            if (dt <= VSMALL)
            {
                return values[i + 1];
            }

            const scalar w = (t - t0) / dt;
            return (1.0 - w) * values[i] + w * values[i + 1];
        }
    }

    return values.last();
}

void fastDynamicFvMesh::calcStructuralModalForces()
{
    structuralModeForce_ = 0.0;
    structuralForceSignal_ = 0.0;

    if (!structuralForceEnabled_ || structuralForceNodeIDs_.size() == 0)
    {
        return;
    }

    if (Pstream::master())
    {
        const scalar t = this->time().value();

        forAll(structuralForceNodeIDs_, forceI)
        {
            const vector nodalForce =
                interpolateStructuralForce(structuralForceTimes_[forceI],
                    structuralForceValues_[forceI], t);

            structuralForceSignal_ += mag(nodalForce);

            const label nTarget = structuralForceNodeIDs_[forceI].size();
            if (nTarget == 0)
            {
                continue;
            }

            forAll(structuralForceNodeIDs_[forceI], localI)
            {
                const label nodeI = structuralForceNodeIDs_[forceI][localI];

                for (label modeI = 0; modeI < nMode_; ++modeI)
                {
                    structuralModeForce_[modeI] +=
                        nodalForce & structModeShapes_[modeI][nodeI];
                }
            }
        }
    }

    Pstream::broadcast(structuralModeForce_);
    Pstream::broadcast(structuralForceSignal_);
}

void fastDynamicFvMesh::readStructuralModeFiles(const fileName& modeDir,
    pointField& structPoints, List<vectorField>& structShapes,
    label& nStructNodes, label& nStructModes)
{
    if (!Pstream::master())
    {
        return;
    }

    const fileName coorPath = modeDir / structNodeCoorFile_;
    std::ifstream coorFile(coorPath);

    if (!coorFile.good())
    {
        FatalErrorInFunction
            << "Cannot open required structural coordinate file " << coorPath
            << nl << "Set 'structNodeCoorFile' correctly in dynamicMeshDict."
            << exit(FatalError);
    }

    string line;
    label lineNo = 0;
    if (!readNextCsvDataLine(coorFile, line, lineNo))
    {
        FatalErrorInFunction
            << "Missing header line in structural coordinate file " << coorPath
            << exit(FatalError);
    }

    scalar dummy = 0.0;
    scalar nNode = 0.0;
    scalar nMode = 0.0;

    if (!parseCsvValues(line, dummy, nNode, nMode) || nNode <= 0 ||
        nMode <= 0)
    {
        FatalErrorInFunction
            << "Invalid header in structural coordinate file " << coorPath
            << " at line " << lineNo
            << ". Expected: dummy, nNode, nMode with positive counts."
            << exit(FatalError);
    }

    nStructNodes = label(nNode);
    nStructModes = label(nMode);

    if (nStructModes != nMode_)
    {
        FatalErrorInFunction << "Structural mode count (" << nStructModes
                             << ") does not match fluid mode count (" << nMode_
                             << ")." << exit(FatalError);
    }

    structPoints.setSize(nStructNodes);
    for (label i = 0; i < nStructNodes; ++i)
    {
        if (!readNextCsvDataLine(coorFile, line, lineNo))
        {
            FatalErrorInFunction
                << "Unexpected end of file while reading structural node " << i
                << " from " << coorPath << exit(FatalError);
        }

        scalar x = 0.0;
        scalar y = 0.0;
        scalar z = 0.0;

        if (!parseCsvValues(line, x, y, z))
        {
            FatalErrorInFunction
                << "Invalid structural coordinate entry for node " << i
                << " in " << coorPath << " at line " << lineNo << ": "
                << line << exit(FatalError);
        }

        structPoints[i] = point(x, y, z);
    }

    structShapes.setSize(nStructModes);

    for (label m = 0; m < nStructModes; ++m)
    {
        const fileName shapePath =
            modeDir / (structNodeDispPrefix_ + std::to_string(m + 1) + ".csv");
        std::ifstream shapeFile(shapePath);

        if (!shapeFile.good())
        {
            FatalErrorInFunction
                << "Cannot open required structural mode shape file "
                << shapePath << exit(FatalError);
        }

        label shapeLineNo = 0;
        if (!readNextCsvDataLine(shapeFile, line, shapeLineNo))
        {
            FatalErrorInFunction << "Missing frequency line in " << shapePath
                                 << exit(FatalError);
        }

        scalar freq = 0.0;
        scalar fileNodeCount = 0.0;
        scalar fileModeCount = 0.0;

        string headerLine(line);
        normalizeCsvLine(headerLine);
        std::stringstream hShape(headerLine);
        if (!(hShape >> freq))
        {
            FatalErrorInFunction << "Failed to read frequency from "
                                 << shapePath << " at line " << shapeLineNo
                                 << ": " << line << exit(FatalError);
        }

        if ((hShape >> fileNodeCount) && (hShape >> fileModeCount))
        {
            if (label(fileNodeCount) != nStructNodes ||
                label(fileModeCount) != nStructModes)
            {
                FatalErrorInFunction
                    << "Structural mode file " << shapePath << " reports "
                    << label(fileNodeCount) << " nodes and "
                    << label(fileModeCount) << " modes, but "
                    << structNodeCoorFile_ << " reports " << nStructNodes
                    << " nodes and " << nStructModes << " modes."
                    << exit(FatalError);
            }
        }

        structShapes[m].setSize(nStructNodes, vector::zero);

        label dataRow = 0;
        while (dataRow < nStructNodes)
        {
            if (!readNextCsvLine(shapeFile, line, shapeLineNo))
            {
                FatalErrorInFunction
                    << "Unexpected end of file while reading structural node "
                    << dataRow << " from " << shapePath << exit(FatalError);
            }

            if (isIgnorableCsvLine(line))
            {
                continue;
            }

            scalar dx = 0.0;
            scalar dy = 0.0;
            scalar dz = 0.0;

            if (!parseCsvValues(line, dx, dy, dz))
            {
                if (dataRow == 0 && !containsNumericToken(line))
                {
                    continue;
                }

                FatalErrorInFunction
                    << "Invalid structural displacement entry for node "
                    << dataRow << " in " << shapePath << " at line "
                    << shapeLineNo << ": " << line << exit(FatalError);
            }

            structShapes[m][dataRow] = vector(dx, dy, dz);
            ++dataRow;
        }
    }
}

void fastDynamicFvMesh::readStructuralStressStrainFiles(const fileName& modeDir)
{
    structModeStress_.clear();
    structModeStrain_.clear();

    if (!modalExtensionStressEnabled_ && !modalExtensionStrainEnabled_)
    {
        return;
    }

    if (structNodeCoords_.size() == 0)
    {
        FatalErrorInFunction
            << "Structural nodes must be loaded before reading structural "
            << "stress/strain files." << exit(FatalError);
    }

    const label nStructNodes = structNodeCoords_.size();

    if (Pstream::master())
    {
        if (modalExtensionStressEnabled_)
        {
            readSixComponentStructuralModeFiles(
                modeDir, structNodeStressPrefix_, nStructNodes, nMode_,
                structNodeCoorFile_, "stress", "SXX,SYY,SZZ,SXY,SYZ,SZX",
                structModeStress_);
        }

        if (modalExtensionStrainEnabled_)
        {
            readSixComponentStructuralModeFiles(
                modeDir, structNodeStrainPrefix_, nStructNodes, nMode_,
                structNodeCoorFile_, "strain", "EXX,EYY,EZZ,EXY,EYZ,EZX",
                structModeStrain_);
        }
    }

    Pstream::broadcast(structModeStress_);
    Pstream::broadcast(structModeStrain_);
}

void fastDynamicFvMesh::readStructuralForceFiles(const fileName& modeDir)
{
    if (!structuralForceEnabled_)
    {
        structuralForceNodeIDs_.clear();
        structuralForceTimes_.clear();
        structuralForceValues_.clear();
        return;
    }

    if (structNodeCoords_.size() == 0)
    {
        FatalErrorInFunction
            << "Structural nodes must be loaded before reading structural "
            << "force files." << exit(FatalError);
    }

    if (Pstream::master())
    {
        structuralForceNodeIDs_.setSize(nStructuralForces_);
        structuralForceTimes_.setSize(nStructuralForces_);
        structuralForceValues_.setSize(nStructuralForces_);

        for (label forceI = 0; forceI < nStructuralForces_; ++forceI)
        {
            const fileName forcePath = modeDir /
                (structuralForceFilePrefix_ + std::to_string(forceI + 1) +
                    ".csv");
            std::ifstream in(forcePath);

            if (!in.good())
            {
                FatalErrorInFunction << "Cannot open structural force file "
                                     << forcePath << exit(FatalError);
            }

            string line;
            label lineNo = 0;

            if (!readNextCsvDataLine(in, line, lineNo))
            {
                FatalErrorInFunction << "Missing coordinate-count line in "
                                     << forcePath << exit(FatalError);
            }

            label nCoords = -1;
            if (!parseCsvValues(line, nCoords) || nCoords <= 0)
            {
                FatalErrorInFunction << "First numeric line in " << forcePath
                                     << " at line " << lineNo
                                     << " must be positive coordinate count: "
                                     << line
                                     << exit(FatalError);
            }

            DynamicList<label> nodeIDs;
            nodeIDs.reserve(nCoords);

            for (label cI = 0; cI < nCoords; ++cI)
            {
                if (!readNextCsvDataLine(in, line, lineNo))
                {
                    FatalErrorInFunction << "Expected " << nCoords
                                         << " coordinate lines in " << forcePath
                                         << ", but file ended early."
                                         << exit(FatalError);
                }

                scalar x = 0.0;
                scalar y = 0.0;
                scalar z = 0.0;

                if (!parseCsvValues(line, x, y, z))
                {
                    FatalErrorInFunction << "Invalid coordinate line in "
                                         << forcePath << " at line " << lineNo
                                         << ": " << line
                                         << exit(FatalError);
                }

                const point target(x, y, z);
                scalar minDistSqr = GREAT;
                label nearest = -1;

                forAll(structNodeCoords_, nodeI)
                {
                    const scalar d2 = magSqr(structNodeCoords_[nodeI] - target);
                    if (d2 < minDistSqr)
                    {
                        minDistSqr = d2;
                        nearest = nodeI;
                    }
                }

                if (nearest < 0 || minDistSqr > sqr(structuralTargetTolerance_))
                {
                    FatalErrorInFunction
                        << "Cannot map structural force coordinate " << target
                        << " in " << forcePath
                        << " within structuralTargetTolerance="
                        << structuralTargetTolerance_ << "."
                        << exit(FatalError);
                }

                bool duplicate = false;
                forAll(nodeIDs, idI)
                {
                    if (nodeIDs[idI] == nearest)
                    {
                        duplicate = true;
                        break;
                    }
                }

                if (!duplicate)
                {
                    nodeIDs.append(nearest);
                }
            }

            DynamicList<scalar> tList;
            DynamicList<vector> fList;

            while (readNextCsvDataLine(in, line, lineNo))
            {
                scalar tt = 0.0;
                scalar fx = 0.0;
                scalar fy = 0.0;
                scalar fz = 0.0;

                if (!parseCsvValues(line, tt, fx, fy, fz))
                {
                    FatalErrorInFunction << "Invalid time-force row in "
                                         << forcePath << " at line " << lineNo
                                         << ": " << line
                                         << exit(FatalError);
                }

                if (!tList.empty() && tt < tList.last())
                {
                    FatalErrorInFunction << "Force times in " << forcePath
                                         << " must be non-decreasing; line "
                                         << lineNo << " has time " << tt
                                         << " after " << tList.last() << "."
                                         << exit(FatalError);
                }

                tList.append(tt);
                fList.append(vector(fx, fy, fz));
            }

            if (tList.empty())
            {
                FatalErrorInFunction << "No force time-series rows found in "
                                     << forcePath << exit(FatalError);
            }

            structuralForceNodeIDs_[forceI].transfer(nodeIDs);
            structuralForceTimes_[forceI].setSize(tList.size());
            structuralForceValues_[forceI].setSize(fList.size());

            forAll(tList, i)
            {
                structuralForceTimes_[forceI][i] = tList[i];
                structuralForceValues_[forceI][i] = fList[i];
            }

            if (structuralForceNodeIDs_[forceI].size() == 0)
            {
                FatalErrorInFunction
                    << "Structural force file " << forcePath
                    << " did not map to any unique structural nodes."
                    << exit(FatalError);
            }
        }
    }

    Pstream::broadcast(structuralForceNodeIDs_);
    Pstream::broadcast(structuralForceTimes_);
    Pstream::broadcast(structuralForceValues_);
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

} // End namespace Foam

// ************************************************************************* //
