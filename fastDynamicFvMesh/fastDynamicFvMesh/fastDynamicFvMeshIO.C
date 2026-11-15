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
#include "IFstream.H"
#include "IStringStream.H"
#include "regIOobject.H"
#include "volFields.H"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>

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

// Read the list payload of a restart-store file.  Handles both the dedicated
// store (a hand-written FoamFile banner followed by "<n>(...)") and a legacy
// OpenFOAM IOField ("<n>" then "(" ... ")"): skip to the first line whose first
// non-blank character is a digit and parse from there up to the closing divider.
template<class T>
bool readRestartList(const fileName& path, List<T>& result)
{
    if (!isFile(path))
    {
        return false;
    }

    std::ifstream is(path.c_str());

    if (!is.good())
    {
        return false;
    }

    std::string line;
    std::string payload;
    bool started = false;

    while (std::getline(is, line))
    {
        const std::string::size_type i = line.find_first_not_of(" \t\r");

        if (!started)
        {
            if
            (
                i != std::string::npos
             && std::isdigit(static_cast<unsigned char>(line[i]))
            )
            {
                started = true;
            }
        }

        if (started)
        {
            if (i != std::string::npos && line[i] == '/')
            {
                break;      // trailing "// ***" divider
            }

            payload += line;
            payload += '\n';
        }
    }

    if (payload.empty())
    {
        return false;
    }

    IStringStream ps(payload);
    ps >> result;

    return !result.empty();
}


// Broadcast a list that the master read from the global restart store.
template<class T>
void broadcastRestartList(List<T>& f)
{
    label n = f.size();
    Pstream::broadcast(n);

    if (!Pstream::master())
    {
        f.setSize(n);
    }

    Pstream::broadcast(f);
}

} // End anonymous namespace

fileName FdmCore::restartStateInstance() const
{
    return fdmMesh_.time().constant() / fileName("fsiRestart");
}

bool FdmCore::readRestartStateStore()
{
    // Fresh starts at t=0 should not pick up stale restart state files.
    // Only non-zero-time starts (continuations/restarts) are allowed to read
    // dedicated or legacy restart modal state.
    if (fdmMesh_.time().value() <= SMALL)
    {
        return false;
    }

    // The store is written by the master into the GLOBAL case directory
    // (<case>/constant/fsiRestart) -- see writeRestartStateStore().  An IOobject
    // resolves a relative instance against time().path(), which is the
    // PROCESSOR directory in a parallel run, so an IOobject read would look in
    // <case>/processorN/constant/fsiRestart and always miss the store.  Read it
    // on the master from the global path instead and broadcast the (nMode-sized)
    // modal state to every rank.
    const fileName restartDir =
        fdmMesh_.time().globalPath() / restartStateInstance();

    // The mode count is part of the restart contract: a checkpoint written
    // with K_old modes cannot be continued with a different K, because the
    // mesh points already contain K_old modes of deformation.  Legacy stores
    // (before nActiveModes existed) always had K_old = N available modes.
    auto checkRestartActiveModes = [&](const fileName& dir, const char* source)
    {
        scalarField storedActive;
        const bool gotStoredActive =
            readRestartList<scalar>(dir/"nActiveModes", storedActive);

        label storedK = nAvailableModes_;

        if (gotStoredActive)
        {
            if (storedActive.size() != 1)
            {
                FatalErrorInFunction
                    << "Malformed nActiveModes record in " << dir
                    << " (" << source << "): expected one value, found "
                    << storedActive.size() << exit(FatalError);
            }

            storedK = label(storedActive[0]);
        }

        if (storedK < 1 || storedK > nAvailableModes_)
        {
            FatalErrorInFunction
                << "Invalid stored nActiveModes=" << storedK << " in " << dir
                << " (" << source << "); the mode set provides "
                << nAvailableModes_ << " available modes." << exit(FatalError);
        }

        if (storedK != nMode_)
        {
            FatalErrorInFunction
                << "Cannot change nActiveModes across restart: checkpoint was "
                << "written with " << storedK << " active modes, but this run "
                << "requests " << nMode_ << " (available modes: "
                << nAvailableModes_ << ")." << nl
                << "Restart from t=0, or use a checkpoint with the same "
                << "nActiveModes." << exit(FatalError);
        }
    };

    bool loadedAny = false;

    if (Pstream::master())
    {
        // NOTE: every field must be read.  Do NOT chain these with '||' -- the
        // read calls would short-circuit and silently skip the remaining
        // fields (a successful modeForce read would leave modeState unread).
        const bool gotModeForce =
            readRestartList<scalar>(restartDir/"modeForce", modeForce_);
        const bool gotModeState =
            readRestartList<vector>(restartDir/"modeState", modeState_);
        const bool gotAppliedModeDisp =
            readRestartList<scalar>
            (
                restartDir/"appliedModeDisp",
                appliedModeDisp_
            );
        const bool gotAppliedModeForce =
            readRestartList<scalar>
            (
                restartDir/"appliedModeForce",
                appliedModeForce_
            );

        loadedAny =
            gotModeForce || gotModeState
         || gotAppliedModeDisp || gotAppliedModeForce;

        if (loadedAny)
        {
            checkRestartActiveModes(restartDir, "restart store");
        }
    }

    label loadedAnyI = (loadedAny ? 1 : 0);
    Pstream::broadcast(loadedAnyI);
    loadedAny = (loadedAnyI != 0);

    if (loadedAny)
    {
        broadcastRestartList(modeForce_);
        broadcastRestartList(modeState_);
        broadcastRestartList(appliedModeDisp_);
        broadcastRestartList(appliedModeForce_);
        return true;
    }

    // Legacy fallback: modal state written into the case-root time directory by
    // older revisions (v1.4.0).  Same master-read/broadcast treatment.
    const fileName legacyDir =
        fdmMesh_.time().globalPath() / fdmMesh_.time().timeName();

    if (Pstream::master())
    {
        // Same rule as above: read every field, no short-circuiting.
        const bool gotModeForce =
            readRestartList<scalar>(legacyDir/"modeForce", modeForce_);
        const bool gotModeState =
            readRestartList<vector>(legacyDir/"modeState", modeState_);
        const bool gotAppliedModeDisp =
            readRestartList<scalar>
            (
                legacyDir/"appliedModeDisp",
                appliedModeDisp_
            );
        const bool gotAppliedModeForce =
            readRestartList<scalar>
            (
                legacyDir/"appliedModeForce",
                appliedModeForce_
            );

        loadedAny =
            gotModeForce || gotModeState
         || gotAppliedModeDisp || gotAppliedModeForce;

        if (loadedAny)
        {
            checkRestartActiveModes(legacyDir, "legacy time-directory store");
        }
    }

    loadedAnyI = (loadedAny ? 1 : 0);
    Pstream::broadcast(loadedAnyI);
    loadedAny = (loadedAnyI != 0);

    if (loadedAny)
    {
        broadcastRestartList(modeForce_);
        broadcastRestartList(modeState_);
        broadcastRestartList(appliedModeDisp_);
        broadcastRestartList(appliedModeForce_);

        if (Pstream::master())
        {
            Info << "Loaded legacy restart modal state from time directory "
                 << fdmMesh_.time().timeName() << endl;
        }
    }

    return loadedAny;
}

void FdmCore::writeRestartStateStore() const
{
    const fileName restartDir =
        fdmMesh_.time().globalPath() / restartStateInstance();
    mkDir(restartDir);

    auto writeGlobalField = [&](const word& name, const auto& field,
                                const word& className) {
        OFstream os(restartDir / name);

        // Full precision: the restart state is restored verbatim, so the
        // default 6 significant digits silently lose accuracy at every
        // checkpoint.
        os.precision(16);

        os << "/*--------------------------------*- C++ "
              "-*----------------------------------*\\"
           << nl;
        os << "| =========                 |                                   "
              "              |"
           << nl;
        os << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD "
              "Toolbox           |"
           << nl;
        os << "|  \\\\    /   O peration     | Version:  " << OPENFOAM
           << "                                  |" << nl;
        os << "|   \\\\  /    A nd           | Website:  www.openfoam.com      "
              "                |"
           << nl;
        os << "|    \\\\/     M anipulation  |                                 "
              "                |"
           << nl;
        os << "\\*-------------------------------------------------------------"
              "--------------*/"
           << nl;
        os << "FoamFile" << nl;
        os << "{" << nl;
        os << "    version     2.0;" << nl;
        os << "    format      ascii;" << nl;
        os << "    class       " << className << ";" << nl;
        os << "    location    \""
           << restartStateInstance().c_str() << "\";" << nl;
        os << "    object      " << name << ";" << nl;
        os << "}" << nl;
        os << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
              "* * * * * * * //"
           << nl << nl;

        os << field;

        os << nl
           << "// "
              "****************************************************************"
              "********* //"
           << nl;
    };

    writeGlobalField("modeState", modeState_, "vectorField");
    writeGlobalField("modeForce", modeForce_, "scalarField");
    writeGlobalField("appliedModeDisp", appliedModeDisp_, "scalarField");
    writeGlobalField("appliedModeForce", appliedModeForce_, "scalarField");

    // Record the active mode count so a restart cannot silently change K.
    // Stored as a one-entry scalarField to reuse the restart list reader.
    writeGlobalField
    (
        "nActiveModes",
        scalarField(1, scalar(nMode_)),
        "scalarField"
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Destructor
// Purpose: release any resources if required. Current implementation relies on
// RAII and does not perform explicit cleanup.
void FdmCore::readControls()
{
    // Read dictionary
    IOdictionary dynamicMeshDict(
        IOobject("dynamicMeshDict", fdmMesh_.time().constant(), fdmMesh_,
            IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE, false));

    if (!dynamicMeshDict.found(fdmTypeName() + "Coeffs"))
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Missing required sub-dictionary '" << fdmTypeName() << "Coeffs' in "
            << dynamicMeshDict.objectPath() << nl
            << "Add constant/dynamicMeshDict entries for " << fdmTypeName() << '.'
            << exit(FatalIOError);
    }

    const dictionary& fdmDict = dynamicMeshDict.subDict(fdmTypeName() + "Coeffs");

    fdmDict.readIfPresent("theta", theta_);

    if (!fdmDict.found("fsiPatches"))
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Missing required entry 'fsiPatches' in sub-dictionary '"
            << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    fdmDict.lookup("fsiPatches") >> fsiPatches_;

    if (!fsiPatches_.size())
    {
        FatalErrorInFunction
            << "At least one FSI patch name must be provided in "
            << dynamicMeshDict.objectPath() << exit(FatalError);
    }

    fdmDict.readIfPresent("mappingTolerance", mappingTolerance_);
    fdmDict.readIfPresent("couplingRelaxation", couplingRelaxation_);
    fdmDict.readIfPresent("couplingAcceleration", couplingAcceleration_);
    fdmDict.readIfPresent("iqnMaxHistory", iqnMaxHistory_);
    fdmDict.readIfPresent("iqnFilter", iqnFilter_);
    fdmDict.readIfPresent("iqnFilterLimit", iqnFilterLimit_);
    fdmDict.readIfPresent("imvjReuse", imvjReuse_);
    fdmDict.readIfPresent("iqnDiagnostics", iqnDiagnostics_);
    fdmDict.readIfPresent("pressureField", pressureFieldName_);
    const bool foundRhoRef = fdmDict.readIfPresent("rhoRef", rhoRef_);
    fdmDict.readIfPresent("pRef", pRef_);
    fdmDict.readIfPresent("maxDispChange", maxDispChange_);
    fdmDict.readIfPresent("writeDiagnostics", writeDiagnosticsEnabled_);
    fdmDict.readIfPresent("trackTiming", trackTimingEnabled_);
    fdmDict.readIfPresent("structuralForceEnabled", structuralForceEnabled_);
    fdmDict.readIfPresent("nStructuralForces", nStructuralForces_);
    fdmDict.readIfPresent(
        "structuralForceFilePrefix", structuralForceFilePrefix_);
    fdmDict.readIfPresent(
        "structuralTargetTolerance", structuralTargetTolerance_);
    fdmDict.readIfPresent("structNodeCoorFile", structNodeCoorFile_);
    fdmDict.readIfPresent("structNodeDispPrefix", structNodeDispPrefix_);
    fdmDict.readIfPresent("modalExtensionEnabled", modalExtensionEnabled_);
    fdmDict.readIfPresent(
        "modalExtensionStressEnabled", modalExtensionStressEnabled_);
    fdmDict.readIfPresent(
        "modalExtensionStrainEnabled", modalExtensionStrainEnabled_);
    fdmDict.readIfPresent("structNodeStressPrefix", structNodeStressPrefix_);
    fdmDict.readIfPresent("structNodeStrainPrefix", structNodeStrainPrefix_);
    fdmDict.readIfPresent("meshRefinementSupport", meshRefinementSupport_);
    fdmDict.readIfPresent(
        "refinementInterpTolerance", refinementInterpTolerance_);

    if (fdmDict.found("modalMass"))
    {
        fdmDict.lookup("modalMass") >> modeMass_;
    }

    if (fdmDict.found("modalDamp"))
    {
        fdmDict.lookup("modalDamp") >> modeDamp_;
    }

    // Number of modes to superpose (first K of the N available modes).  The
    // upper bound is validated in readModeFiles(), because N is only known
    // after reading mode/FluidNodeCoor.csv.  -1 means "all available modes".
    nActiveModes_ = -1;
    fdmDict.readIfPresent("nActiveModes", nActiveModes_);

    // Prescribed modal motion (forced-oscillation / GAF testing).  All keys
    // are optional and default to "off", so existing cases are unaffected.
    fdmDict.readIfPresent("prescribedMotion", prescribedMotion_);
    fdmDict.readIfPresent("prescribedAmplitude", prescribedAmplitude_);
    fdmDict.readIfPresent("prescribedPhase", prescribedPhase_);
    fdmDict.readIfPresent("prescribedFrequency", prescribedFrequency_);

    if (prescribedMotion_ && Pstream::master())
    {
        Info<< "fastDynamicFvMesh: prescribed modal motion ON"
            << "  amplitude=" << prescribedAmplitude_
            << "  phase=" << prescribedPhase_
            << "  frequency=" << prescribedFrequency_ << " Hz"
            << " (<=0 means each mode's own frequency)" << endl;
    }

    if (couplingRelaxation_ <= 0 || couplingRelaxation_ > 1)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'couplingRelaxation' must be in the range (0, 1] in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if
    (
        couplingAcceleration_ != "constant"
     && couplingAcceleration_ != "aitken"
     && couplingAcceleration_ != "IQN-ILS"
     && couplingAcceleration_ != "IQN-IMVJ"
    )
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'couplingAcceleration' must be one of constant, aitken, "
            << "IQN-ILS or IQN-IMVJ in sub-dictionary '" << fdmTypeName()
            << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (iqnMaxHistory_ < 1)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'iqnMaxHistory' must be >= 1 in sub-dictionary '"
            << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (imvjReuse_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'imvjReuse' must be >= 0 in sub-dictionary '"
            << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (iqnFilterLimit_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'iqnFilterLimit' must be > 0 in sub-dictionary '"
            << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if
    (
        iqnFilter_ != "none"
     && iqnFilter_ != "QR2"
     && iqnFilter_ != "exponential"
    )
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'iqnFilter' must be one of none, QR2 or exponential in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (iqnFilter_ == "exponential" && iqnFilterLimit_ >= 1)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'iqnFilterLimit' must be < 1 when 'iqnFilter' is "
            << "exponential in sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (Pstream::master())
    {
        Info << "Coupling acceleration: " << couplingAcceleration_
             << " (initialRelaxation=" << couplingRelaxation_ << ')' << endl;
    }

    if (theta_ < 1)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'theta' must be >= 1 for Wilson-Theta integration in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (mappingTolerance_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'mappingTolerance' must be positive in sub-dictionary '"
            << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (maxDispChange_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'maxDispChange' must be positive in sub-dictionary '"
            << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    forAll(modeMass_, modeI)
    {
        if (modeMass_[modeI] <= 0)
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'modalMass' must contain positive values; index "
                << modeI << " is " << modeMass_[modeI] << " in "
                << dynamicMeshDict.objectPath() << exit(FatalIOError);
        }
    }

    forAll(modeDamp_, modeI)
    {
        if (modeDamp_[modeI] < 0)
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'modalDamp' must contain non-negative values; index "
                << modeI << " is " << modeDamp_[modeI] << " in "
                << dynamicMeshDict.objectPath() << exit(FatalIOError);
        }
    }

    if (foundRhoRef && rhoRef_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'rhoRef' must be positive in sub-dictionary '" << fdmTypeName()
            << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (structuralTargetTolerance_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'structuralTargetTolerance' must be positive in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (refinementInterpTolerance_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'refinementInterpTolerance' must be positive in "
            << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (structuralForceEnabled_)
    {
        if (nStructuralForces_ <= 0)
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'nStructuralForces' must be > 0 in "
                << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
                << dynamicMeshDict.objectPath() << exit(FatalIOError);
        }

        if (structuralForceFilePrefix_.empty())
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'structuralForceFilePrefix' cannot be empty in "
                << "sub-dictionary '" << fdmTypeName() << "Coeffs' of "
                << dynamicMeshDict.objectPath() << exit(FatalIOError);
        }
    }
    else if
    (
        Pstream::master()
     && (
            fdmDict.found("nStructuralForces")
         || fdmDict.found("structuralForceFilePrefix")
         || fdmDict.found("structuralTargetTolerance")
        )
    )
    {
        WarningInFunction
            << "Structural-force entries are present, but "
            << "'structuralForceEnabled' is false; structural force files "
            << "will not be read." << endl;
    }

    if
    (
        (modalExtensionStressEnabled_ || modalExtensionStrainEnabled_)
     && !modalExtensionEnabled_
    )
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entries 'modalExtensionStressEnabled' and/or "
            << "'modalExtensionStrainEnabled' require 'modalExtensionEnabled' "
            << "to be true in sub-dictionary '" << fdmTypeName() << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (modalExtensionEnabled_)
    {
        if (structNodeStressPrefix_.empty())
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'structNodeStressPrefix' cannot be empty when "
                << "'modalExtensionEnabled' is true in sub-dictionary '"
                << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.objectPath()
                << exit(FatalIOError);
        }

        if (structNodeStrainPrefix_.empty())
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'structNodeStrainPrefix' cannot be empty when "
                << "'modalExtensionEnabled' is true in sub-dictionary '"
                << fdmTypeName() << "Coeffs' of " << dynamicMeshDict.objectPath()
                << exit(FatalIOError);
        }
    }
    else if
    (
        Pstream::master()
     && (
            fdmDict.found("modalExtensionEnabled")
         || fdmDict.found("modalExtensionStressEnabled")
         || fdmDict.found("modalExtensionStrainEnabled")
         || fdmDict.found("structNodeStressPrefix")
         || fdmDict.found("structNodeStrainPrefix")
        )
    )
    {
        WarningInFunction
            << "Modal-extension entries are present, but "
            << "'modalExtensionEnabled' is false; structural displacement, "
            << "stress and strain recovery will not be performed." << endl;
    }
}

// readLegacyParameters(modeDir)
// Purpose: on the master process, attempt to read legacy runtime parameters
// from mode/FluidPara.csv. What it does:
//  - Only the master process attempts to open and parse the file; slaves return
//  early.
//  - The parser is robust to commas and spaces (commas are replaced with spaces
//  before tokenizing).
//  - It first reads two legacy identifier scalars (legacyFsiId, legacyFluidId)
//  then reads an initial velocity per mode if present.
//  - If the file or expected entries are missing, the function warns and leaves
//  initVelocity_ at zero.
void FdmCore::readLegacyParameters(const fileName& modeDir)
{
    if (!Pstream::master())
    {
        return;
    }

    fileName paraPath = modeDir / "FluidPara.csv";
    std::ifstream paraFile(paraPath);

    if (!paraFile.good())
    {
        Info << "Legacy parameter file " << paraPath
             << " not found; using zero initial modal velocities." << endl;
        return;
    }

    label lineNo = 0;

    auto readCsvScalar = [&](scalar& value) -> bool {
        string line;

        while (readNextCsvDataLine(paraFile, line, lineNo))
        {
            if (parseCsvValues(line, value))
            {
                return true;
            }

            WarningInFunction
                << "Skipping malformed scalar row " << lineNo << " in "
                << paraPath << ": " << line << endl;
        }

        return false;
    };

    scalar legacyFsiId = -1;
    scalar legacyFluidId = -1;

    if (!readCsvScalar(legacyFsiId) || !readCsvScalar(legacyFluidId))
    {
        WarningInFunction << "Cannot read legacy ids from " << paraPath
                          << ". Initial modal velocities remain zero." << endl;

        return;
    }

    Info << "Read legacy FluidPara.csv ids (fsi=" << legacyFsiId
         << ", fluid=" << legacyFluidId
         << "). OpenFOAM uses patch names from dynamicMeshDict instead."
         << endl;

    for (label modeI = 0; modeI < nMode_; ++modeI)
    {
        scalar value = 0.0;

        if (!readCsvScalar(value))
        {
            WarningInFunction << "Missing initial velocity for mode " << modeI
                              << " in " << paraPath
                              << ". Remaining initial velocities stay at zero."
                              << endl;

            break;
        }

        initVelocity_[modeI] = value;
    }
}

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
void FdmCore::readModeFiles(const fileName& modeDir,
    List<point>& csvPoints, List<List<vector>>& csvShapes, label& nCsvNodes)
{
    if (!Pstream::master())
    {
        return;
    }

    Info << "Reading mode coordinates..." << endl;

    fileName coorPath = modeDir / "FluidNodeCoor.csv";
    Info << "Trying to open: " << coorPath << endl;

    std::ifstream file(coorPath);
    if (!file.good())
    {
        FatalErrorInFunction
            << "Cannot open required mode coordinate file " << coorPath << nl
            << "Provide mode/FluidNodeCoor.csv in the case root before "
            << "running fastDynamicFvMesh." << exit(FatalError);
    }
    else
    {
        scalar dummy, nNode, nModeHeader;
        string line;
        label lineNo = 0;
        if (!readNextCsvDataLine(file, line, lineNo) ||
            !parseCsvValues(line, dummy, nNode, nModeHeader) || nNode <= 0 ||
            nModeHeader <= 0)
        {
            FatalErrorInFunction
                << "Invalid FluidNodeCoor.csv header in " << coorPath
                << " at line " << lineNo << nl
                << "Expected three comma-separated values with positive "
                << "node and mode counts." << exit(FatalError);
        }

        nCsvNodes = label(nNode);
        nAvailableModes_ = label(nModeHeader);

        // Resolve the requested active mode count.  The header gives the
        // number of available mode files (N); dynamicMeshDict may cap the
        // superposition to the first K modes.  A missing/invalid request is
        // handled here because nMode_ is not known in readControls().
        if (nActiveModes_ < 0)
        {
            nActiveModes_ = nAvailableModes_;
        }

        if (nActiveModes_ < 1 || nActiveModes_ > nAvailableModes_)
        {
            FatalErrorInFunction
                << "Entry 'nActiveModes' is " << nActiveModes_
                << ", but " << coorPath << " provides "
                << nAvailableModes_ << " available modes." << nl
                << "Set 1 <= nActiveModes <= " << nAvailableModes_
                << " (or omit it to superpose all available modes)."
                << exit(FatalError);
        }

        nMode_ = nActiveModes_;

        Info << "  Found " << nCsvNodes << " nodes and " << nAvailableModes_
             << " available modes in CSV; superposing the first " << nMode_
             << " mode(s)." << endl;

        // Keep the FluidNodeCoor header authoritative: all N mode files must
        // exist, even though only the first K are read below.
        for (label m = 0; m < nAvailableModes_; ++m)
        {
            const fileName shapePath =
                modeDir / ("FluidNodeDisp" + std::to_string(m + 1) + ".csv");

            if (!isFile(shapePath))
            {
                FatalErrorInFunction
                    << "Missing mode shape file " << shapePath << nl
                    << "FluidNodeCoor.csv declares " << nAvailableModes_
                    << " available modes, so FluidNodeDisp1.."
                    << nAvailableModes_ << ".csv must all exist."
                    << exit(FatalError);
            }
        }

        modeFreq_.setSize(nMode_); // Resize here on master
        csvPoints.setSize(nCsvNodes);
        csvShapes.setSize(nMode_);
        forAll(csvShapes, m) csvShapes[m].setSize(nCsvNodes);

        for (label i = 0; i < nCsvNodes; ++i)
        {
            scalar x, y, z;
            if (!readNextCsvDataLine(file, line, lineNo))
            {
                FatalErrorInFunction
                    << "Unexpected end of file while reading node " << i
                    << " from " << coorPath << exit(FatalError);
            }
            if (!parseCsvValues(line, x, y, z))
            {
                FatalErrorInFunction << "Invalid coordinate entry for node "
                                     << i << " in " << coorPath << " at line "
                                     << lineNo << ": " << line
                                     << exit(FatalError);
            }

            csvPoints[i] = point(x, y, z);
        }

        // Read mode shapes
        for (label m = 0; m < nMode_; ++m)
        {
            fileName shapePath =
                modeDir / ("FluidNodeDisp" + std::to_string(m + 1) + ".csv");
            std::ifstream mFile(shapePath);

            if (!mFile.good())
            {
                FatalErrorInFunction << "Cannot open required mode shape file "
                                     << shapePath << exit(FatalError);
            }

            string line;
            label shapeLineNo = 0;
            if (!readNextCsvDataLine(mFile, line, shapeLineNo))
            {
                FatalErrorInFunction << "Missing frequency line in "
                                     << shapePath << exit(FatalError);
            }

            scalar freq = 0.0;
            scalar fileNodeCount = 0.0;
            scalar fileModeCount = 0.0;

            string headerLine(line);
            normalizeCsvLine(headerLine);
            std::stringstream ss(headerLine);
            if (!(ss >> freq))
            {
                FatalErrorInFunction << "Failed to read frequency from "
                                     << shapePath << " at line " << shapeLineNo
                                     << ": " << line << exit(FatalError);
            }

            if ((ss >> fileNodeCount) && (ss >> fileModeCount))
            {
                if
                (
                    label(fileNodeCount) != nCsvNodes
                 || label(fileModeCount) < nMode_
                )
                {
                    FatalErrorInFunction
                        << "Mode file " << shapePath << " reports "
                        << label(fileNodeCount) << " nodes and "
                        << label(fileModeCount)
                        << " modes, but FluidNodeCoor.csv reports " << nCsvNodes
                        << " nodes and " << nAvailableModes_
                        << " available modes (nActiveModes=" << nMode_ << ")."
                        << exit(FatalError);
                }
            }

            modeFreq_[m] = freq;
            Info << "  Mode " << m << " Freq: " << freq << endl;

            label dataRow = 0;

            while (dataRow < nCsvNodes)
            {
                if (!readNextCsvLine(mFile, line, shapeLineNo))
                {
                    FatalErrorInFunction
                        << "Unexpected end of file while reading node "
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
                        << "Invalid displacement entry for node " << dataRow
                        << " in " << shapePath << " at line " << shapeLineNo
                        << ": " << line << exit(FatalError);
                }

                csvShapes[m][dataRow] = vector(dx, dy, dz);
                ++dataRow;
            }
        }
    }
}
} // End namespace Foam
