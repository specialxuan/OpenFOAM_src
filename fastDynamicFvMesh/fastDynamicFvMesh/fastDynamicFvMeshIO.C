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

fileName fastDynamicFvMesh::restartStateInstance() const
{
    return this->time().constant() / fileName("fsiRestart");
}

bool fastDynamicFvMesh::readRestartStateStore()
{
    // Fresh starts at t=0 should not pick up stale restart state files.
    // Only non-zero-time starts (continuations/restarts) are allowed to read
    // dedicated or legacy restart modal state.
    if (this->time().value() <= SMALL)
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
        this->time().globalPath() / restartStateInstance();

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
        this->time().globalPath() / this->time().timeName();

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
                 << this->time().timeName() << endl;
        }
    }

    return loadedAny;
}

void fastDynamicFvMesh::writeRestartStateStore() const
{
    const fileName restartDir =
        this->time().globalPath() / restartStateInstance();
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
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Destructor
// Purpose: release any resources if required. Current implementation relies on
// RAII and does not perform explicit cleanup.
void fastDynamicFvMesh::readControls()
{
    // Read dictionary
    IOdictionary dynamicMeshDict(
        IOobject("dynamicMeshDict", this->time().constant(), *this,
            IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE, false));

    if (!dynamicMeshDict.found(typeName + "Coeffs"))
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Missing required sub-dictionary '" << typeName << "Coeffs' in "
            << dynamicMeshDict.objectPath() << nl
            << "Add constant/dynamicMeshDict entries for " << typeName << '.'
            << exit(FatalIOError);
    }

    const dictionary& fdmDict = dynamicMeshDict.subDict(typeName + "Coeffs");

    fdmDict.readIfPresent("theta", theta_);

    if (!fdmDict.found("fsiPatches"))
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Missing required entry 'fsiPatches' in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
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
    fdmDict.readIfPresent("trackAmrTiming", trackAmrTimingEnabled_);
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
    fdmDict.readIfPresent(
        "refinementMappingDiagnostics", refinementMappingDiagnostics_);
    scalar refinementFaceMapTolerance = 0.0;
    const bool foundRefinementFaceMapTolerance = fdmDict.readIfPresent(
        "refinementFaceMapTolerance", refinementFaceMapTolerance);
    fdmDict.readIfPresent("refinementMinCellVolume", refinementMinCellVolume_);
    fdmDict.readIfPresent("refinementMinEdgeLength", refinementMinEdgeLength_);
    fdmDict.readIfPresent("amrLayerDiagnostics", amrLayerDiagnosticsEnabled_);
    fdmDict.readIfPresent("amrLayerDiagnosticsZ", amrLayerDiagnosticsZ_);
    fdmDict.readIfPresent(
        "amrLayerDiagnosticsThickness", amrLayerDiagnosticsThickness_);
    fdmDict.readIfPresent(
        "amrLayerDiagnosticsNeighborDepth", amrLayerDiagnosticsNeighborDepth_);
    fdmDict.readIfPresent(
        "amrLayerDiagnosticsMaxRows", amrLayerDiagnosticsMaxRows_);
    fdmDict.readIfPresent(
        "amrLayerDiagnosticsFilePrefix", amrLayerDiagnosticsFilePrefix_);

    if (dynamicMeshDict.found("dynamicRefineFvMeshCoeffs"))
    {
        const dictionary& refineDict =
            dynamicMeshDict.subDict("dynamicRefineFvMeshCoeffs");
        refineDict.readIfPresent(
            "useGradIndicator", refinementUseGradIndicator_);
        refineDict.readIfPresent(
            "gradIndicatorField", refinementGradIndicatorField_);
        refineDict.readIfPresent("gradIndicatorMagnitudeField",
            refinementGradIndicatorMagnitudeField_);
    }

    if (runtimeRefinementEnabled_ && !meshRefinementSupport_ && Pstream::master())
    {
        WarningInFunction
            << "dynamicRefineFvMeshCoeffs is present, but "
            << "'meshRefinementSupport' is false in " << typeName
            << "Coeffs. Runtime topology changes can invalidate modal "
            << "shape and FSI face caches; set meshRefinementSupport true "
            << "for AMR-coupled fastDynamicFvMesh cases." << endl;
    }

    if (fdmDict.found("modalMass"))
    {
        fdmDict.lookup("modalMass") >> modeMass_;
    }

    if (fdmDict.found("modalDamp"))
    {
        fdmDict.lookup("modalDamp") >> modeDamp_;
    }

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
            << "sub-dictionary '" << typeName << "Coeffs' of "
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
            << "IQN-ILS or IQN-IMVJ in sub-dictionary '" << typeName
            << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (iqnMaxHistory_ < 1)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'iqnMaxHistory' must be >= 1 in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (imvjReuse_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'imvjReuse' must be >= 0 in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (iqnFilterLimit_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'iqnFilterLimit' must be > 0 in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
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
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (iqnFilter_ == "exponential" && iqnFilterLimit_ >= 1)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'iqnFilterLimit' must be < 1 when 'iqnFilter' is "
            << "exponential in sub-dictionary '" << typeName << "Coeffs' of "
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
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (mappingTolerance_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'mappingTolerance' must be positive in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (maxDispChange_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'maxDispChange' must be positive in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
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
            << "Entry 'rhoRef' must be positive in sub-dictionary '" << typeName
            << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (structuralTargetTolerance_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'structuralTargetTolerance' must be positive in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (refinementInterpTolerance_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'refinementInterpTolerance' must be positive in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (refinementFaceMapTolerance < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'refinementFaceMapTolerance' must be >= 0 in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }
    else if (foundRefinementFaceMapTolerance && Pstream::master())
    {
        WarningInFunction
            << "Entry 'refinementFaceMapTolerance' is deprecated and ignored. "
            << "Topology mapping is now fully driven by mapPolyMesh lineage."
            << endl;
    }

    if (refinementMinCellVolume_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'refinementMinCellVolume' must be >= 0 in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (refinementMinEdgeLength_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'refinementMinEdgeLength' must be >= 0 in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (amrLayerDiagnosticsThickness_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'amrLayerDiagnosticsThickness' must be >= 0 in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (amrLayerDiagnosticsNeighborDepth_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'amrLayerDiagnosticsNeighborDepth' must be >= 0 in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (amrLayerDiagnosticsMaxRows_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'amrLayerDiagnosticsMaxRows' must be >= 0 in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (amrLayerDiagnosticsEnabled_ && amrLayerDiagnosticsFilePrefix_.empty())
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'amrLayerDiagnosticsFilePrefix' cannot be empty when "
            << "'amrLayerDiagnostics' is enabled in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    if (refinementUseGradIndicator_ && refinementGradIndicatorField_.empty())
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'gradIndicatorField' cannot be empty when "
            << "'useGradIndicator' is enabled in dynamicRefineFvMeshCoeffs of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (Pstream::master() && runtimeRefinementEnabled_ &&
        refinementUseGradIndicator_)
    {
        if (refinementGradIndicatorMagnitudeField_.empty())
        {
            Info << "Runtime AMR indicator uses cached |grad("
                 << refinementGradIndicatorField_
                 << ")| with dynamicRefineFvMesh thresholds." << endl;
        }
        else
        {
            Info << "Runtime AMR indicator uses registered field '"
                 << refinementGradIndicatorMagnitudeField_
                 << "' as |grad(" << refinementGradIndicatorField_
                 << ")| with dynamicRefineFvMesh thresholds." << endl;
        }
    }

    if (Pstream::master() && amrLayerDiagnosticsEnabled_)
    {
        Info << "Runtime AMR bottom-layer diagnostics enabled: z="
             << amrLayerDiagnosticsZ_ << ", thickness="
             << amrLayerDiagnosticsThickness_ << ", neighborDepth="
             << amrLayerDiagnosticsNeighborDepth_ << ", maxRows="
             << amrLayerDiagnosticsMaxRows_ << ", filePrefix='"
             << amrLayerDiagnosticsFilePrefix_ << "'." << endl;
    }

    if (structuralForceEnabled_)
    {
        if (nStructuralForces_ <= 0)
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'nStructuralForces' must be > 0 in "
                << "sub-dictionary '" << typeName << "Coeffs' of "
                << dynamicMeshDict.objectPath() << exit(FatalIOError);
        }

        if (structuralForceFilePrefix_.empty())
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'structuralForceFilePrefix' cannot be empty in "
                << "sub-dictionary '" << typeName << "Coeffs' of "
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
            << "to be true in sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (modalExtensionEnabled_)
    {
        if (structNodeStressPrefix_.empty())
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'structNodeStressPrefix' cannot be empty when "
                << "'modalExtensionEnabled' is true in sub-dictionary '"
                << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
                << exit(FatalIOError);
        }

        if (structNodeStrainPrefix_.empty())
        {
            FatalIOErrorInFunction(dynamicMeshDict)
                << "Entry 'structNodeStrainPrefix' cannot be empty when "
                << "'modalExtensionEnabled' is true in sub-dictionary '"
                << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
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
void fastDynamicFvMesh::readLegacyParameters(const fileName& modeDir)
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
void fastDynamicFvMesh::readModeFiles(const fileName& modeDir,
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
        scalar dummy, nNode, nMode;
        string line;
        label lineNo = 0;
        if (!readNextCsvDataLine(file, line, lineNo) ||
            !parseCsvValues(line, dummy, nNode, nMode) || nNode <= 0 ||
            nMode <= 0)
        {
            FatalErrorInFunction
                << "Invalid FluidNodeCoor.csv header in " << coorPath
                << " at line " << lineNo << nl
                << "Expected three comma-separated values with positive "
                << "node and mode counts." << exit(FatalError);
        }

        nCsvNodes = label(nNode);
        nMode_ = label(nMode);

        Info << "  Found " << nCsvNodes << " nodes and " << nMode_
             << " modes in CSV." << endl;

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
                if (label(fileNodeCount) != nCsvNodes ||
                    label(fileModeCount) != nMode_)
                {
                    FatalErrorInFunction
                        << "Mode file " << shapePath << " reports "
                        << label(fileNodeCount) << " nodes and "
                        << label(fileModeCount)
                        << " modes, but FluidNodeCoor.csv reports " << nCsvNodes
                        << " nodes and " << nMode_ << " modes."
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


bool fastDynamicFvMesh::writeRefinedPointModeShapes() const
{
    if (!meshRefinementSupport_ || nMode_ <= 0)
    {
        return true;
    }

    const label nPoints = this->points().size();

    if (modeShapes_.size() != nMode_)
    {
        return true;
    }

    forAll(modeShapes_, modeI)
    {
        if (modeShapes_[modeI].size() != nPoints)
        {
            return true;
        }
    }

    const fileName timeDir = this->time().timeName();
    const fileName dir = this->time().path()/timeDir/polyMesh::meshSubDir;

    mkDir(dir);

    const fileName shapesFile = dir/"refinedPointModeShapes";

    const labelList& pointLevel = meshCutter_.pointLevel();

    label nLevel = 0;
    forAll(pointLevel, pI)
    {
        if (pointLevel[pI] > 0)
        {
            ++nLevel;
        }
    }

    std::ofstream os(shapesFile.c_str());

    if (!os.good())
    {
        return false;
    }

    // Full precision: these values are restored verbatim on restart.
    os.precision(17);

    os << "FoamFile" << std::endl;
    os << "{" << std::endl;
    os << "    version     2.0;" << std::endl;
    os << "    format      ascii;" << std::endl;
    os << "    class       pointModeShapes;" << std::endl;
    os << "    location    \"" << (timeDir/"polyMesh").c_str() << "\";"
       << std::endl;
    os << "    object      refinedPointModeShapes;" << std::endl;
    os << "}" << std::endl;
    os << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
          "* * * * * * * //" << std::endl << std::endl;

    // Body:
    //   nBase    <number of level-0 points on this processor>
    //   nTime    <write time>
    //   nMode    <number of modes>
    //   nRefined <number of recorded (pointLevel>0) points>
    //   (
    //     <pointIndex>
    //     <mode0 x y z>
    //     <mode1 x y z>
    //     ...
    //   )
    os << "nBase " << (nPoints - nLevel) << std::endl;
    os << "nTime " << timeDir << std::endl;
    os << "nMode " << nMode_ << std::endl;
    os << "nRefined " << nLevel << std::endl;
    os << "(" << std::endl;

    forAll(pointLevel, pI)
    {
        if (pointLevel[pI] <= 0)
        {
            continue;
        }

        os << pI << std::endl;

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            const vector& s = modeShapes_[modeI][pI];

            os << s.x() << ' ' << s.y() << ' ' << s.z() << std::endl;
        }
    }

    os << ")" << std::endl << std::endl;
    os << "// ***************************************************************"
          "********** //" << std::endl;

    return os.good();
}


label fastDynamicFvMesh::readRefinedPointModeShapes()
{
    recordedRefinedPointIDs_.setSize(0);
    recordedRefinedModeShapes_.setSize(0);
    recordedRefinedNMode_ = 0;

    // NOTE: nMode_ is not known yet -- it is read from the mode CSVs inside
    // readModeShapes(), which runs after this.  The recorded mode count is
    // therefore taken from the file and validated when the values are applied.
    if (!meshRefinementSupport_)
    {
        return 0;
    }

    const fileName shapesFile =
        this->time().path()/this->time().timeName()/polyMesh::meshSubDir
      / "refinedPointModeShapes";

    if (!exists(shapesFile))
    {
        return 0;
    }

    try
    {
        std::ifstream is(shapesFile.c_str());

        if (!is.good())
        {
            return 0;
        }

        // Skip everything up to the end-of-header marker
        std::string line;
        while (std::getline(is, line))
        {
            if (line.find("* * * * * * * * * * * * *") != std::string::npos)
            {
                break;
            }
        }

        std::vector<std::string> tok;

        for (std::string t; is >> t; )
        {
            tok.push_back(t);
        }

        if
        (
            tok.size() < 8
         || tok[0] != "nBase"
         || tok[2] != "nTime"
         || tok[4] != "nMode"
         || tok[6] != "nRefined"
        )
        {
            WarningInFunction
                << "Malformed AMR refined-point mode-shape file " << shapesFile
                << "; ignoring it and falling back to nearest-corner "
                   "interpolation." << endl;

            return 0;
        }

        const label nBaseLocal = std::stol(tok[1]);
        const std::string fileTime = tok[3];
        const label fileNMode = std::stol(tok[5]);
        const label nRefined = std::stol(tok[7]);

        if (fileTime != this->time().timeName())
        {
            WarningInFunction
                << "AMR refined-point mode-shape file " << shapesFile
                << " was written for time " << fileTime
                << " but this run starts from " << this->time().timeName()
                << "; ignoring it." << endl;

            return 0;
        }

        if (fileNMode <= 0)
        {
            WarningInFunction
                << "AMR refined-point mode-shape file " << shapesFile
                << " holds no modes; ignoring it." << endl;

            return 0;
        }

        recordedRefinedNMode_ = fileNMode;

        const label nPoints = this->points().size();
        const labelList& pointLevel = meshCutter_.pointLevel();

        if (pointLevel.size() != nPoints)
        {
            WarningInFunction
                << "pointLevel size (" << pointLevel.size()
                << ") does not match the mesh point count (" << nPoints
                << "); ignoring the recorded refined-point mode shapes." << endl;

            return 0;
        }

        recordedRefinedPointIDs_.setSize(nRefined);
        recordedRefinedModeShapes_.setSize(nRefined);

        const std::size_t perPoint = 1 + std::size_t(3*fileNMode);

        label nLoaded = 0;
        std::size_t k = 8;

        if (k < tok.size() && tok[k] == "(")
        {
            ++k;
        }

        for (label e = 0; e < nRefined; ++e)
        {
            if (k + perPoint > tok.size())
            {
                break;
            }

            const label pI = std::stol(tok[k++]);

            vectorField shapes(fileNMode);

            for (label modeI = 0; modeI < fileNMode; ++modeI)
            {
                const scalar vx = std::stod(tok[k++]);
                const scalar vy = std::stod(tok[k++]);
                const scalar vz = std::stod(tok[k++]);

                shapes[modeI] = vector(vx, vy, vz);
            }

            recordedRefinedPointIDs_[nLoaded] = pI;
            recordedRefinedModeShapes_[nLoaded] = shapes;
            ++nLoaded;
        }

        recordedRefinedPointIDs_.setSize(nLoaded);
        recordedRefinedModeShapes_.setSize(nLoaded);

        label nLevel = 0;
        forAll(pointLevel, pI)
        {
            if (pointLevel[pI] > 0)
            {
                ++nLevel;
            }
        }

        const label globalLoaded = returnReduce(nLoaded, sumOp<label>());
        const label globalLevel = returnReduce(nLevel, sumOp<label>());
        const label globalBase = returnReduce(nBaseLocal, sumOp<label>());
        const label globalPoints = returnReduce(nPoints, sumOp<label>());

        if (Pstream::master())
        {
            Info << "Read AMR refined-point mode shapes from " << shapesFile
                 << nl
                 << "    recorded points : " << globalLoaded
                 << "   refined points : " << globalLevel
                 << "   base points : " << globalBase
                 << "   mesh points : " << globalPoints << endl;
        }

        if (globalLoaded < globalLevel && Pstream::master())
        {
            WarningInFunction
                << (globalLevel - globalLoaded)
                << " refinement-added points were not recorded; they will be "
                   "rebuilt by nearest-corner interpolation." << endl;
        }

        return nLoaded;
    }
    catch (const std::exception& e)
    {
        WarningInFunction
            << "Failed to parse AMR refined-point mode-shape file "
            << shapesFile << " (" << e.what() << "); ignoring it and falling "
               "back to nearest-corner interpolation." << endl;

        recordedRefinedPointIDs_.setSize(0);
        recordedRefinedModeShapes_.setSize(0);

        return 0;
    }
}

bool fastDynamicFvMesh::writeFsiReferenceProjection() const
{
    if (!meshRefinementSupport_ || fsiPatches_.size() == 0 || nMode_ <= 0)
    {
        return true;
    }

    const fileName dir = this->time().path()/this->time().constant();

    mkDir(dir);

    const fileName projFile = dir/"fsiReferenceProjection";

    std::ofstream os(projFile.c_str());

    if (!os.good())
    {
        return false;
    }

    os.precision(17);

    // The reference (pre-refine) faces and their mode projection are invariant
    // for the whole computation, so this file is written once and never
    // rewritten.  Per-processor (each rank owns part of the FSI patch).
    os << "nPatches " << fsiPatches_.size() << std::endl;
    os << "nMode " << nMode_ << std::endl;

    forAll(fsiPatches_, i)
    {
        if
        (
            referenceFsiFaceAreas_.size() != fsiPatches_.size()
         || referenceFaceModeProjection_.size() != fsiPatches_.size()
        )
        {
            continue;
        }

        const scalarField& areas = referenceFsiFaceAreas_[i];
        const label nRef = areas.size();

        os << fsiPatches_[i] << ' ' << nRef << std::endl;

        for (label refI = 0; refI < nRef; ++refI)
        {
            os << areas[refI] << ' ';
        }
        os << std::endl;

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            const vectorField& proj = referenceFaceModeProjection_[i][modeI];

            for (label refI = 0; refI < nRef; ++refI)
            {
                const vector& v = proj[refI];
                os << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
            }
            os << std::endl;
        }
    }

    return os.good();
}


bool fastDynamicFvMesh::writeFsiFaceReferenceMap() const
{
    if (!meshRefinementSupport_ || fsiPatches_.size() == 0)
    {
        return true;
    }

    const fileName dir =
        this->time().path()/this->time().timeName()/polyMesh::meshSubDir;

    mkDir(dir);

    const fileName mapFile = dir/"fsiFaceReferenceMap";

    std::ofstream os(mapFile.c_str());

    if (!os.good())
    {
        return false;
    }

    os.precision(17);

    // Current (possibly split) FSI face -> reference face(s) + area weights.
    // Unlike the projection this depends on the topology, so it is written with
    // every checkpoint.  Per-processor.
    os << "nPatches " << fsiPatches_.size() << std::endl;

    forAll(fsiPatches_, i)
    {
        if
        (
            fsiFaceToReferenceFaces_.size() != fsiPatches_.size()
         || fsiFaceToReferenceWeights_.size() != fsiPatches_.size()
        )
        {
            continue;
        }

        const List<labelList>& refs = fsiFaceToReferenceFaces_[i];
        const List<scalarField>& weights = fsiFaceToReferenceWeights_[i];

        const label nFaces = refs.size();

        os << fsiPatches_[i] << ' ' << nFaces << std::endl;

        for (label faceI = 0; faceI < nFaces; ++faceI)
        {
            const labelList& refsFace = refs[faceI];

            os << refsFace.size();

            forAll(refsFace, refI)
            {
                os << ' ' << refsFace[refI];
            }

            if (faceI < weights.size())
            {
                const scalarField& w = weights[faceI];

                forAll(w, refI)
                {
                    os << ' ' << w[refI];
                }
            }

            os << std::endl;
        }
    }

    return os.good();
}


bool fastDynamicFvMesh::readFsiReferenceCache()
{
    if (!meshRefinementSupport_ || fsiPatches_.size() == 0 || nMode_ <= 0)
    {
        return false;
    }

    try
    {
        // ---- 1) invariant reference projection (written once, constant/) ----
        const fileName projFile =
            this->time().path()/this->time().constant()/"fsiReferenceProjection";

        if (!exists(projFile))
        {
            return false;
        }

        std::ifstream pis(projFile.c_str());

        if (!pis.good())
        {
            return false;
        }

        {
            std::string w;
            label nPatches = -1;
            label nMode = -1;

            pis >> w >> nPatches >> w >> nMode;

            if
            (
                !pis.good()
             || nPatches != label(fsiPatches_.size())
             || nMode != nMode_
            )
            {
                WarningInFunction
                    << "FSI reference projection " << projFile
                    << " is inconsistent with this case (patches "
                    << nPatches << " vs " << fsiPatches_.size()
                    << ", modes " << nMode << " vs " << nMode_
                    << "); ignoring it." << endl;

                return false;
            }
        }

        referenceFsiFaceAreas_.setSize(fsiPatches_.size());
        referenceFaceModeProjection_.setSize(fsiPatches_.size());

        forAll(fsiPatches_, i)
        {
            std::string name;
            label nRef = 0;

            pis >> name >> nRef;

            if (name != std::string(fsiPatches_[i]) || nRef <= 0)
            {
                WarningInFunction
                    << "Malformed FSI reference projection " << projFile
                    << "; ignoring it." << endl;

                referenceFsiFaceAreas_.setSize(0);
                referenceFaceModeProjection_.setSize(0);

                return false;
            }

            scalarField& areas = referenceFsiFaceAreas_[i];
            areas.setSize(nRef);

            for (label refI = 0; refI < nRef; ++refI)
            {
                pis >> areas[refI];
            }

            referenceFaceModeProjection_[i].setSize(nMode_);

            for (label modeI = 0; modeI < nMode_; ++modeI)
            {
                vectorField& proj = referenceFaceModeProjection_[i][modeI];
                proj.setSize(nRef);

                for (label refI = 0; refI < nRef; ++refI)
                {
                    scalar x = 0, y = 0, z = 0;
                    pis >> x >> y >> z;
                    proj[refI] = vector(x, y, z);
                }
            }
        }

        // ---- 2) current-face -> reference-face map for THIS checkpoint ----
        const fileName mapFile =
            this->time().path()/this->time().timeName()/polyMesh::meshSubDir
          / "fsiFaceReferenceMap";

        if (!exists(mapFile))
        {
            referenceFsiFaceAreas_.setSize(0);
            referenceFaceModeProjection_.setSize(0);

            return false;
        }

        std::ifstream mis(mapFile.c_str());

        if (!mis.good())
        {
            referenceFsiFaceAreas_.setSize(0);
            referenceFaceModeProjection_.setSize(0);

            return false;
        }

        {
            std::string w;
            label nPatches = -1;

            mis >> w >> nPatches;

            if (!mis.good() || nPatches != label(fsiPatches_.size()))
            {
                WarningInFunction
                    << "FSI face reference map " << mapFile
                    << " is inconsistent with this case; ignoring it." << endl;

                referenceFsiFaceAreas_.setSize(0);
                referenceFaceModeProjection_.setSize(0);

                return false;
            }
        }

        const polyBoundaryMesh& pbm = this->boundaryMesh();

        fsiFaceToReferenceFaces_.setSize(fsiPatches_.size());
        fsiFaceToReferenceWeights_.setSize(fsiPatches_.size());

        forAll(fsiPatches_, i)
        {
            std::string name;
            label nFaces = 0;

            mis >> name >> nFaces;

            const label patchID = pbm.findPatchID(fsiPatches_[i]);
            const label nFacesNow =
                (patchID >= 0 ? pbm[patchID].size() : -1);

            if (name != std::string(fsiPatches_[i]) || nFaces != nFacesNow)
            {
                WarningInFunction
                    << "FSI face reference map " << mapFile << " does not cover "
                    << "patch '" << fsiPatches_[i] << "' of this mesh ("
                    << nFaces << " vs " << nFacesNow
                    << " faces); ignoring it." << endl;

                referenceFsiFaceAreas_.setSize(0);
                referenceFaceModeProjection_.setSize(0);
                fsiFaceToReferenceFaces_.setSize(0);
                fsiFaceToReferenceWeights_.setSize(0);

                return false;
            }

            List<labelList>& refs = fsiFaceToReferenceFaces_[i];
            List<scalarField>& weights = fsiFaceToReferenceWeights_[i];

            refs.setSize(nFaces);
            weights.setSize(nFaces);

            label nOneToMany = 0;

            for (label faceI = 0; faceI < nFaces; ++faceI)
            {
                label nRefs = 0;
                mis >> nRefs;

                if (nRefs <= 0)
                {
                    WarningInFunction
                        << "Malformed FSI face reference map " << mapFile
                        << "; ignoring it." << endl;

                    referenceFsiFaceAreas_.setSize(0);
                    referenceFaceModeProjection_.setSize(0);
                    fsiFaceToReferenceFaces_.setSize(0);
                    fsiFaceToReferenceWeights_.setSize(0);

                    return false;
                }

                labelList& r = refs[faceI];
                scalarField& w = weights[faceI];

                r.setSize(nRefs);
                w.setSize(nRefs);

                for (label refI = 0; refI < nRefs; ++refI)
                {
                    mis >> r[refI];
                }

                for (label refI = 0; refI < nRefs; ++refI)
                {
                    mis >> w[refI];
                }

                if (nRefs > 1)
                {
                    ++nOneToMany;
                }
            }

            const label globalOneToMany =
                returnReduce(nOneToMany, sumOp<label>());

            if (Pstream::master() && globalOneToMany > 0)
            {
                Info << "Restored FSI reference mapping on patch '"
                     << fsiPatches_[i] << "': " << nFaces << " current faces, "
                     << globalOneToMany << " one-to-many." << endl;
            }
        }

        referenceFsiBuilt_ = true;

        label localRefTotal = 0;
        forAll(referenceFsiFaceAreas_, i)
        {
            localRefTotal += referenceFsiFaceAreas_[i].size();
        }
        const label nRefTotal = returnReduce(localRefTotal, sumOp<label>());

        if (Pstream::master())
        {
            Info << "Restored FSI reference-face cache from " << projFile
                 << " + " << mapFile << nl
                 << "    reference faces : " << nRefTotal
                 << "   (projection is invariant; the current-face mapping is "
                    "from this checkpoint)" << endl;
        }

        return true;
    }
    catch (const std::exception& e)
    {
        WarningInFunction
            << "Failed to parse the FSI reference-face cache ("
            << e.what() << "); ignoring it and re-snapshotting from the "
               "startup mesh." << endl;

        referenceFsiFaceAreas_.setSize(0);
        referenceFaceModeProjection_.setSize(0);
        fsiFaceToReferenceFaces_.setSize(0);
        fsiFaceToReferenceWeights_.setSize(0);

        return false;
    }
}

} // End namespace Foam
