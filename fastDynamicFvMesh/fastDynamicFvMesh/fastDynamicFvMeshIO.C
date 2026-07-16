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

    const fileName restartInstance = restartStateInstance();
    bool loadedAny = false;

    auto readScalarFieldIfPresent = [&](const word& name,
                                        IOField<scalar>& target) -> bool {
        IOobject io(name, restartInstance, *this, IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE, false);

        if (!isFile(io.objectPath()))
        {
            return false;
        }

        IOField<scalar> fld(io);
        target = fld;
        return true;
    };

    auto readVectorFieldIfPresent = [&](const word& name,
                                        IOField<vector>& target) -> bool {
        IOobject io(name, restartInstance, *this, IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE, false);

        if (!isFile(io.objectPath()))
        {
            return false;
        }

        IOField<vector> fld(io);
        target = fld;
        return true;
    };

    loadedAny = readScalarFieldIfPresent("modeForce", modeForce_) || loadedAny;
    loadedAny = readVectorFieldIfPresent("modeState", modeState_) || loadedAny;
    loadedAny = readScalarFieldIfPresent("appliedModeDisp", appliedModeDisp_) ||
        loadedAny;
    loadedAny =
        readScalarFieldIfPresent("appliedModeForce", appliedModeForce_) ||
        loadedAny;

    if (!loadedAny)
    {
        const fileName legacyInstance = this->time().timeName();

        auto readLegacyScalar = [&](const word& name,
                                    IOField<scalar>& target) -> bool {
            IOobject io(name, legacyInstance, *this, IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE, false);

            if (!isFile(io.objectPath()))
            {
                return false;
            }

            IOField<scalar> fld(io);
            target = fld;
            return true;
        };

        auto readLegacyVector = [&](const word& name,
                                    IOField<vector>& target) -> bool {
            IOobject io(name, legacyInstance, *this, IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE, false);

            if (!isFile(io.objectPath()))
            {
                return false;
            }

            IOField<vector> fld(io);
            target = fld;
            return true;
        };

        loadedAny = readLegacyScalar("modeForce", modeForce_) || loadedAny;
        loadedAny = readLegacyVector("modeState", modeState_) || loadedAny;
        loadedAny =
            readLegacyScalar("appliedModeDisp", appliedModeDisp_) || loadedAny;
        loadedAny = readLegacyScalar("appliedModeForce", appliedModeForce_) ||
            loadedAny;

        if (loadedAny && Pstream::master())
        {
            Info << "Loaded legacy restart modal state from time directory "
                 << legacyInstance << endl;
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
        os << "    location    \"" << restartStateInstance() << "\";" << nl;
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

    if (couplingRelaxation_ <= 0 || couplingRelaxation_ > 1)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'couplingRelaxation' must be in the range (0, 1] in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
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

} // End namespace Foam
