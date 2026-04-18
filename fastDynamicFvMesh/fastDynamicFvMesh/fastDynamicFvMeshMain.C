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
    scalar refinementFaceMapTolerance = 0.0;
    const bool foundRefinementFaceMapTolerance = fdmDict.readIfPresent(
        "refinementFaceMapTolerance", refinementFaceMapTolerance);
    fdmDict.readIfPresent("refinementMinCellVolume", refinementMinCellVolume_);
    fdmDict.readIfPresent("refinementMinEdgeLength", refinementMinEdgeLength_);

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
void fastDynamicFvMesh::readModeShapes()
{
    // 1. Read files on master
    List<point> csvPoints;
    List<List<vector>> csvShapes;
    label nCsvNodes = 0;
    fileName modeDir = this->time().path() / "mode";

    if (Pstream::parRun() && !exists(modeDir))
    {
        modeDir = this->time().path().path() / "mode";
    }

    readModeFiles(modeDir, csvPoints, csvShapes, nCsvNodes);

    // Broadcast sizes
    Pstream::broadcast(nMode_);
    Pstream::broadcast(nCsvNodes);

    // Resize local arrays
    if (!Pstream::master())
        modeFreq_.setSize(nMode_); // Only slaves need resize now

    if (modeMass_.size() == 0)
    {
        modeMass_.setSize(nMode_, 1.0);
    }
    else if (modeMass_.size() != nMode_)
    {
        FatalErrorInFunction << "modalMass size (" << modeMass_.size()
                             << ") does not match number of modes (" << nMode_
                             << ")" << exit(FatalError);
    }

    if (modeDamp_.size() == 0)
    {
        modeDamp_.setSize(nMode_, 0.0);
    }
    else if (modeDamp_.size() != nMode_)
    {
        FatalErrorInFunction << "modalDamp size (" << modeDamp_.size()
                             << ") does not match number of modes (" << nMode_
                             << ")" << exit(FatalError);
    }

    if (nMode_ <= 0)
    {
        FatalErrorInFunction
            << "No fluid modes were loaded from " << modeDir << nl
            << "Ensure mode/FluidNodeCoor.csv and mode/FluidNodeDisp*.csv are "
            << "present and valid." << exit(FatalError);
    }

    modeForce_.setSize(nMode_, 0.0);
    modePressureForce_.setSize(nMode_, 0.0);
    modeShearForce_.setSize(nMode_, 0.0);
    modeForce0_.setSize(nMode_, 0.0);
    modeState_.setSize(nMode_, vector::zero);
    modeState0_.setSize(nMode_, vector::zero);
    initVelocity_.setSize(nMode_, 0.0);
    appliedModeDisp_.setSize(nMode_, 0.0);

    readLegacyParameters(modeDir);

    // Broadcast data
    Pstream::broadcast(modeFreq_);
    Pstream::broadcast(csvPoints);
    Pstream::broadcast(initVelocity_);

    // Resize csvShapes on all procs to ensure loop consistency
    csvShapes.setSize(nMode_);

    forAll(csvShapes, m) Pstream::broadcast(csvShapes[m]);

    // Map to local mesh points.
    // Default to current mesh points. This is required on latestTime restarts
    // after runtime AMR because current point count can differ from
    // constant/polyMesh/points.
    pointField mappingPoints(this->points());
    pointField basePoints;

    // Try to read points from constant/polyMesh (or
    // processor*/constant/polyMesh) Note: in parallel, the polyMesh directory
    // is inside the processor directory. IOobject will handle processor path
    // automatically if we use NO_WRITE.

    // We construct an IOobject to look for "points" in the constant instance.
    IOobject ioPoints("points", this->time().constant(), polyMesh::meshSubDir,
        *this, IOobject::MUST_READ, IOobject::NO_WRITE,
        false // register
    );

    // Custom wrapper to force expecting "vectorField" class name
    // regardless of what Field<vector>::typeName is defined as.
    class vectorFieldIO : public IOField<vector>
    {
    public:
        // Force the type name to match the file header
        virtual const word& type() const
        {
            static const word typeName_val("vectorField");
            return typeName_val;
        }

        vectorFieldIO(const IOobject& io) : IOField<vector>(io) {}
    };

    // Attempt to read with the custom wrapper.
    bool readOk = false;
    try
    {
        // Don't check header via typeHeaderOk as it might check the wrong type
        // Just verify it exists (MUST_READ handles this usually, but we want to
        // catch it)
        if (ioPoints.headerClassName() == "vectorField")
        {
            vectorFieldIO pIO(ioPoints);
            basePoints = pIO;
            readOk = true;
        }
        else if (ioPoints.typeHeaderOk<vectorFieldIO>(true))
        {
            vectorFieldIO pIO(ioPoints);
            basePoints = pIO;
            readOk = true;
        }
    }
    catch (...)
    {
    }

    if (readOk)
    {
        if (basePoints.size() == mappingPoints.size())
        {
            mappingPoints = basePoints;

            if (Pstream::master())
            {
                Info << "  Read " << basePoints.size()
                     << " base points for mapping from "
                     << ioPoints.objectPath() << endl;
            }
        }
        else if (Pstream::master())
        {
            WarningInFunction
                << "constant/polyMesh points count (" << basePoints.size()
                << ") differs from current mesh points count ("
                << mappingPoints.size() << ") on startup/restart. "
                << "Using current mesh points for modal mapping." << endl;
        }
    }
    else if (Pstream::master())
    {
        WarningInFunction
            << "Could not read base points from constant/polyMesh "
            << "(expected class vectorField). Using current mesh points "
            << "for modal mapping." << endl;
    }

    modeShapes_.setSize(nMode_);
    forAll(modeShapes_, m)
    {
        modeShapes_[m].setSize(mappingPoints.size());
        modeShapes_[m] = vector::zero;
    }

    const labelList* pointLevelPtr = nullptr;
    bool topologyOnlyForRefinedPoints = false;
    if (meshRefinementSupport_)
    {
        const labelList& pointLevel = meshCutter_.pointLevel();
        if (pointLevel.size() != mappingPoints.size())
        {
            FatalErrorInFunction
                << "meshCutter pointLevel size (" << pointLevel.size()
                << ") does not match mesh point count (" << mappingPoints.size()
                << ")." << exit(FatalError);
        }

        label localRefinedPointCount = 0;
        forAll(pointLevel, pointI)
        {
            if (pointLevel[pointI] > 0)
            {
                ++localRefinedPointCount;
            }
        }

        const label globalRefinedPointCount =
            returnReduce(localRefinedPointCount, sumOp<label>());
        topologyOnlyForRefinedPoints = (globalRefinedPointCount > 0);
        pointLevelPtr = &pointLevel;

        if (Pstream::master() && topologyOnlyForRefinedPoints)
        {
            Info
                << "Refined startup mesh detected (" << globalRefinedPointCount
                << " points with pointLevel>0). "
                << "Geometric CSV mapping is restricted to level-0 points; "
                << "refined points must be recovered by topology interpolation."
                << endl;
        }
    }

    // Brute force search (simplest correct implementation)
    // Only perform mapping if we have CSV nodes
    label mappedCount = 0;
    boolList mappedMask(mappingPoints.size(), false);
    if (nCsvNodes > 0)
    {
        scalar tolSqr = sqr(mappingTolerance_);

        forAll(mappingPoints, pI)
        {
            if (topologyOnlyForRefinedPoints && (*pointLevelPtr)[pI] > 0)
            {
                continue;
            }

            const point& p = mappingPoints[pI];
            scalar minDistSqr = GREAT;
            label nearestIdx = -1;

            forAll(csvPoints, cI)
            {
                scalar d = magSqr(p - csvPoints[cI]);
                if (d < minDistSqr)
                {
                    minDistSqr = d;
                    nearestIdx = cI;
                }
            }

            if (minDistSqr < tolSqr && nearestIdx != -1)
            {
                mappedCount++;
                mappedMask[pI] = true;
                for (label m = 0; m < nMode_; ++m)
                {
                    modeShapes_[m][pI] = csvShapes[m][nearestIdx];
                }
            }
        }
    }

    if (meshRefinementSupport_)
    {
        syncMappedModeShapes(mappedMask);
    }

    label totalMapped = returnReduce(mappedCount, sumOp<label>());
    const label totalPoints =
        returnReduce(label(mappingPoints.size()), sumOp<label>());

    if (meshRefinementSupport_ && totalMapped > 0 && totalMapped < totalPoints)
    {
        boolList knownMask(mappedMask);

        label prevGlobalKnown = totalMapped;
        for (label passI = 0; passI < 4; ++passI)
        {
            const boolList sourceKnownMask(knownMask);

            interpolateModeShapesFromKnown(
                knownMask, sourceKnownMask, 2, true, "edge midpoint");
            interpolateModeShapesFromKnown(
                knownMask, sourceKnownMask, 4, false, "face midpoint");
            interpolateModeShapesFromKnown(
                knownMask, sourceKnownMask, 8, false, "cell midpoint");

            syncMappedModeShapes(knownMask);

            label localKnownPass = 0;
            forAll(knownMask, pointI)
            {
                if (knownMask[pointI])
                {
                    ++localKnownPass;
                }
            }

            const label totalKnownPass =
                returnReduce(localKnownPass, sumOp<label>());

            if (totalKnownPass <= prevGlobalKnown)
            {
                break;
            }

            prevGlobalKnown = totalKnownPass;
        }

        syncMappedModeShapes(knownMask);

        label localKnown = 0;
        forAll(knownMask, pointI)
        {
            if (knownMask[pointI])
            {
                ++localKnown;
            }
        }

        const label totalKnown = returnReduce(localKnown, sumOp<label>());

        if (Pstream::master())
        {
            Info << "Refinement-aware startup interpolation raised mapped "
                    "points "
                 << "from " << totalMapped << " to " << totalKnown << " out of "
                 << totalPoints << "." << endl;
        }

        if (totalKnown != totalPoints && topologyOnlyForRefinedPoints)
        {
            FatalErrorInFunction
                << "With meshRefinementSupport enabled on a refined startup "
                   "mesh, "
                << (totalPoints - totalKnown)
                << " points remain unmapped after 2/4/8 topology interpolation."
                << nl
                << "No geometric fallback is allowed for AMR-generated points."
                << nl << "Check refinement lineage consistency and "
                << "refinementInterpTolerance." << exit(FatalError);
        }
        else if (totalKnown != totalPoints && Pstream::master())
        {
            WarningInFunction << "With meshRefinementSupport enabled, "
                              << (totalPoints - totalKnown)
                              << " points remain unmapped after 2/4/8 topology "
                                 "interpolation. "
                              << "These points keep zero modal displacement "
                                 "until they can be "
                              << "mapped by later topology updates." << endl;
        }

        totalMapped = totalKnown;
    }

    if (Pstream::master())
    {
        Info << "Mapped " << totalMapped << " points out of " << totalPoints
             << " mesh points." << endl;
    }

    if (totalMapped == 0)
    {
        FatalErrorInFunction
            << "Mapped 0 mesh points from " << modeDir << nl
            << "Check that the mode CSV coordinates correspond to the current "
            << "mesh and that mappingTolerance is large enough."
            << exit(FatalError);
    }

    if (Pstream::master() && totalMapped != totalPoints)
    {
        WarningInFunction
            << "Mapped " << totalMapped << " of " << totalPoints
            << " mesh points. Unmapped points will remain stationary." << endl;
    }

    buildForceProjectionCaches();

    if (structuralForceEnabled_)
    {
        pointField structPoints;
        List<vectorField> structShapes;
        label nStructNodes = 0;
        label nStructModes = 0;

        readStructuralModeFiles(
            modeDir, structPoints, structShapes, nStructNodes, nStructModes);

        Pstream::broadcast(nStructNodes);
        Pstream::broadcast(nStructModes);
        Pstream::broadcast(structPoints);

        if (!Pstream::master())
        {
            structShapes.setSize(nStructModes);
        }

        forAll(structShapes, m)
        {
            Pstream::broadcast(structShapes[m]);
        }

        structNodeCoords_.transfer(structPoints);
        structModeShapes_.transfer(structShapes);
        structuralModeForce_.setSize(nMode_, 0.0);
        readStructuralForceFiles(modeDir);

        if (Pstream::master())
        {
            label totalNodes = 0;
            forAll(structuralForceNodeIDs_, fi)
            {
                totalNodes += structuralForceNodeIDs_[fi].size();
            }

            Info << "Loaded " << nStructuralForces_
                 << " structural force files with mapped node references="
                 << totalNodes << endl;
        }
    }
    else
    {
        structuralModeForce_.setSize(nMode_, 0.0);
        structuralForceNodeIDs_.clear();
        structuralForceTimes_.clear();
        structuralForceValues_.clear();
        structNodeCoords_.clear();
        structModeShapes_.clear();
    }
}

void fastDynamicFvMesh::buildForceProjectionCaches(const mapPolyMesh* mpmPtr)
{
    const polyBoundaryMesh& pbm = this->boundaryMesh();

    const List<List<labelList>> oldFaceToReferenceFaces(
        fsiFaceToReferenceFaces_);
    const List<List<scalarField>> oldFaceToReferenceWeights(
        fsiFaceToReferenceWeights_);
    const bool hasOldReferenceMapping =
        oldFaceToReferenceFaces.size() == fsiPatches_.size() &&
        oldFaceToReferenceWeights.size() == fsiPatches_.size();

    fsiPatchIDs_.setSize(fsiPatches_.size(), -1);
    fsiPolyPatches_.setSize(fsiPatches_.size(), nullptr);
    faceModeProjection_.setSize(fsiPatches_.size());
    pressureScaleCache_.setSize(fsiPatches_.size());
    pressureScaleInitialized_.setSize(fsiPatches_.size(), false);
    fsiFaceToReferenceFaces_.setSize(fsiPatches_.size());
    fsiFaceToReferenceWeights_.setSize(fsiPatches_.size());

    if (!meshRefinementSupport_)
    {
        referenceFsiBuilt_ = false;
        referenceFaceModeProjection_.clear();
        referenceFsiFaceAreas_.clear();
    }
    else if (!referenceFsiBuilt_ ||
        referenceFaceModeProjection_.size() != fsiPatches_.size() ||
        referenceFsiFaceAreas_.size() != fsiPatches_.size())
    {
        referenceFaceModeProjection_.setSize(fsiPatches_.size());
        referenceFsiFaceAreas_.setSize(fsiPatches_.size());
        referenceFsiBuilt_ = false;
    }

    const labelList* oldPatchStartsPtr =
        mpmPtr ? &mpmPtr->oldPatchStarts() : nullptr;
    const labelList* oldPatchSizesPtr =
        mpmPtr ? &mpmPtr->oldPatchSizes() : nullptr;
    const List<objectMap>* facesFromFacesMapPtr =
        mpmPtr ? &mpmPtr->facesFromFacesMap() : nullptr;

    Map<labelList> facesFromFacesByNewFace;
    if (facesFromFacesMapPtr)
    {
        facesFromFacesByNewFace.resize(2 * facesFromFacesMapPtr->size() + 1);
        forAll(*facesFromFacesMapPtr, mapI)
        {
            const objectMap& obj = (*facesFromFacesMapPtr)[mapI];
            facesFromFacesByNewFace.set(obj.index(), obj.masterObjects());
        }
    }

    forAll(fsiPatches_, i)
    {
        const label patchID = pbm.findPatchID(fsiPatches_[i]);
        if (patchID == -1)
        {
            FatalErrorInFunction << "Configured FSI patch '" << fsiPatches_[i]
                                 << "' was not found in boundaryMesh."
                                 << exit(FatalError);
        }

        fsiPatchIDs_[i] = patchID;
        fsiPolyPatches_[i] = &pbm[patchID];

        const polyPatch& pp = *fsiPolyPatches_[i];
        const label nFaces = pp.size();

        pressureScaleCache_[i].setSize(nFaces, 1.0);
        faceModeProjection_[i].setSize(nMode_);

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            vectorField& coeff = faceModeProjection_[i][modeI];
            coeff.setSize(nFaces, vector::zero);

            forAll(pp, faceI)
            {
                const labelList& fPoints = pp[faceI];
                vector shapeFace = vector::zero;

                forAll(fPoints, fp)
                {
                    shapeFace += modeShapes_[modeI][fPoints[fp]];
                }

                if (fPoints.size() > 0)
                {
                    shapeFace /= fPoints.size();
                }

                coeff[faceI] = shapeFace;
            }
        }

        List<labelList>& faceToRefs = fsiFaceToReferenceFaces_[i];
        List<scalarField>& faceToRefWeights = fsiFaceToReferenceWeights_[i];
        faceToRefs.setSize(nFaces);
        faceToRefWeights.setSize(nFaces);

        if (!meshRefinementSupport_)
        {
            forAll(faceToRefs, faceI)
            {
                faceToRefs[faceI].setSize(1, faceI);
                faceToRefWeights[faceI].setSize(1, 1.0);
            }
            continue;
        }

        if (!referenceFsiBuilt_)
        {
            referenceFsiFaceAreas_[i].setSize(nFaces, 0.0);
            referenceFaceModeProjection_[i].setSize(nMode_);

            forAll(pp, faceI)
            {
                const vector area = pp.faceAreas()[faceI];
                const scalar areaMag = mag(area);
                referenceFsiFaceAreas_[i][faceI] = areaMag;
            }

            for (label modeI = 0; modeI < nMode_; ++modeI)
            {
                referenceFaceModeProjection_[i][modeI] =
                    faceModeProjection_[i][modeI];
            }

            forAll(faceToRefs, faceI)
            {
                faceToRefs[faceI].setSize(1, faceI);
                faceToRefWeights[faceI].setSize(1, 1.0);
            }

            continue;
        }

        if (!mpmPtr)
        {
            FatalErrorInFunction << "Reference cache is already initialized, "
                                    "but mapPolyMesh was not provided "
                                 << "while rebuilding topology mappings."
                                 << exit(FatalError);
        }

        if (!hasOldReferenceMapping)
        {
            FatalErrorInFunction
                << "Missing previous topology-to-reference mapping on patch '"
                << fsiPatches_[i] << "' while processing topology update."
                << exit(FatalError);
        }

        const label oldPatchStart = oldPatchStartsPtr->size() > patchID
            ? (*oldPatchStartsPtr)[patchID]
            : -1;
        const label oldPatchSize = oldPatchSizesPtr->size() > patchID
            ? (*oldPatchSizesPtr)[patchID]
            : -1;

        if (oldPatchStart < 0 || oldPatchSize < 0)
        {
            FatalErrorInFunction
                << "Cannot read old patch start/size for patch '"
                << fsiPatches_[i] << "' (patchID=" << patchID << ")."
                << exit(FatalError);
        }

        const List<labelList>& oldPatchRefs = oldFaceToReferenceFaces[i];
        const List<scalarField>& oldPatchWeights = oldFaceToReferenceWeights[i];

        if (oldPatchRefs.size() != oldPatchSize ||
            oldPatchWeights.size() != oldPatchSize)
        {
            FatalErrorInFunction
                << "Old mapping size mismatch on patch '" << fsiPatches_[i]
                << "': " << "oldPatchSize=" << oldPatchSize
                << ", refs=" << oldPatchRefs.size()
                << ", weights=" << oldPatchWeights.size() << exit(FatalError);
        }

        const scalarField& refAreas = referenceFsiFaceAreas_[i];

        label mappedCurrentFaces = 0;
        label oneToManyFaces = 0;

        forAll(pp, faceI)
        {
            labelList& refs = faceToRefs[faceI];
            scalarField& weights = faceToRefWeights[faceI];
            refs.clear();
            weights.clear();

            const label globalFace = pp.start() + faceI;
            DynamicList<label> oldFaces;

            const label oldMasterFace = mpmPtr->faceMap()[globalFace];
            if (oldMasterFace >= 0)
            {
                oldFaces.append(oldMasterFace);
            }

            const auto mergedIter = facesFromFacesByNewFace.cfind(globalFace);
            if (mergedIter.found())
            {
                const labelList& masters = mergedIter();
                forAll(masters, masterI)
                {
                    oldFaces.appendUniq(masters[masterI]);
                }
            }

            if (oldFaces.empty())
            {
                FatalErrorInFunction << "Cannot derive old faces for patch '"
                                     << fsiPatches_[i] << "', face " << faceI
                                     << " (global face " << globalFace << ")."
                                     << exit(FatalError);
            }

            Map<scalar> refAreaAccum(2 * oldFaces.size() + 1);

            forAll(oldFaces, oldI)
            {
                const label oldFace = oldFaces[oldI];
                if (oldFace < oldPatchStart ||
                    oldFace >= oldPatchStart + oldPatchSize)
                {
                    continue;
                }

                const label oldLocalFace = oldFace - oldPatchStart;
                const labelList& prevRefs = oldPatchRefs[oldLocalFace];
                const scalarField& prevWeights = oldPatchWeights[oldLocalFace];

                if (prevRefs.size() != prevWeights.size() || prevRefs.empty())
                {
                    FatalErrorInFunction
                        << "Invalid old mapping entry for patch '"
                        << fsiPatches_[i] << "', old local face "
                        << oldLocalFace << "." << exit(FatalError);
                }

                forAll(prevRefs, refI)
                {
                    const label refFace = prevRefs[refI];
                    if (refFace < 0 || refFace >= refAreas.size())
                    {
                        FatalErrorInFunction << "Out-of-range reference face "
                                             << refFace << " on patch '"
                                             << fsiPatches_[i] << "'."
                                             << exit(FatalError);
                    }

                    const scalar clampedWeight = max(prevWeights[refI], 0.0);
                    const scalar contribArea =
                        max(refAreas[refFace], VSMALL) * clampedWeight;
                    refAreaAccum(refFace, 0.0) += contribArea;
                }
            }

            if (refAreaAccum.empty())
            {
                FatalErrorInFunction
                    << "Topology mapping on patch '" << fsiPatches_[i]
                    << "', face " << faceI
                    << " yielded no reference faces after patch filtering."
                    << exit(FatalError);
            }

            const labelList sortedRefs(refAreaAccum.sortedToc());
            refs.setSize(sortedRefs.size());
            weights.setSize(sortedRefs.size(), 0.0);

            scalar sumArea = 0.0;
            forAll(sortedRefs, refI)
            {
                sumArea += refAreaAccum[sortedRefs[refI]];
            }

            if (sumArea <= VSMALL)
            {
                FatalErrorInFunction
                    << "Non-positive accumulated reference area on patch '"
                    << fsiPatches_[i] << "', face " << faceI << "."
                    << exit(FatalError);
            }

            scalar sumW = 0.0;
            forAll(sortedRefs, refI)
            {
                refs[refI] = sortedRefs[refI];
                weights[refI] = refAreaAccum[sortedRefs[refI]] / sumArea;
                sumW += weights[refI];
            }

            if (sumW <= VSMALL)
            {
                FatalErrorInFunction << "Zero total mapping weight on patch '"
                                     << fsiPatches_[i] << "', face " << faceI
                                     << "." << exit(FatalError);
            }

            if (mag(sumW - 1.0) > 1e-10)
            {
                weights /= sumW;
            }

            ++mappedCurrentFaces;
            if (refs.size() > 1)
            {
                ++oneToManyFaces;
            }
        }

        label globalMappedCurrentFaces = mappedCurrentFaces;
        label globalOneToManyFaces = oneToManyFaces;
        reduce(globalMappedCurrentFaces, sumOp<label>());
        reduce(globalOneToManyFaces, sumOp<label>());

        if (Pstream::master())
        {
            Info << "Topology reference mapping on patch '" << fsiPatches_[i]
                 << "': mapped " << globalMappedCurrentFaces
                 << " current faces, one-to-many faces=" << globalOneToManyFaces
                 << "." << endl;
        }
    }

    if (meshRefinementSupport_ && !referenceFsiBuilt_)
    {
        referenceFsiBuilt_ = true;
        if (Pstream::master())
        {
            Info << "Initialized reference FSI face cache for "
                    "refinement/unrefine-aware "
                 << "force aggregation." << endl;
        }
    }
}

void fastDynamicFvMesh::syncMappedModeShapes(boolList& mappedMask)
{
    if (!meshRefinementSupport_ || !Pstream::parRun() || nMode_ <= 0)
    {
        return;
    }

    if (mappedMask.size() != modeShapes_[0].size())
    {
        FatalErrorInFunction << "mappedMask size (" << mappedMask.size()
                             << ") does not match mode-shape point count ("
                             << modeShapes_[0].size() << ")."
                             << exit(FatalError);
    }

    List<label> ownerProc(mappedMask.size(), -1);
    const label myProc = Pstream::myProcNo();

    forAll(mappedMask, pointI)
    {
        if (mappedMask[pointI])
        {
            ownerProc[pointI] = myProc;
        }
    }

    syncTools::syncPointList(*this, ownerProc, maxEqOp<label>(), label(-1));
    syncTools::syncPointList(*this, mappedMask, orEqOp<bool>(), false);

    for (label modeI = 0; modeI < nMode_; ++modeI)
    {
        vectorField ownedValues(modeShapes_[modeI].size(), vector::zero);

        forAll(ownedValues, pointI)
        {
            if (ownerProc[pointI] == myProc)
            {
                ownedValues[pointI] = modeShapes_[modeI][pointI];
            }
        }

        syncTools::syncPointList(
            *this, ownedValues, plusEqOp<vector>(), vector::zero);
        modeShapes_[modeI].transfer(ownedValues);
    }
}

void fastDynamicFvMesh::syncSharedPointPositions(pointField& pointValues) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (pointValues.size() != this->points().size())
    {
        FatalErrorInFunction << "pointValues size (" << pointValues.size()
                             << ") does not match mesh point count ("
                             << this->points().size() << ")."
                             << exit(FatalError);
    }

    List<label> ownerProc(pointValues.size(), Pstream::myProcNo());
    syncTools::syncPointList(*this, ownerProc, maxEqOp<label>(), label(-1));

    pointField ownedValues(pointValues.size(), point::zero);
    const label myProc = Pstream::myProcNo();

    forAll(ownedValues, pointI)
    {
        if (ownerProc[pointI] == myProc)
        {
            ownedValues[pointI] = pointValues[pointI];
        }
    }

    syncTools::syncPointPositions(
        *this, ownedValues, plusEqOp<point>(), point::zero);
    pointValues.transfer(ownedValues);
}

label fastDynamicFvMesh::interpolateModeShapesFromKnown(boolList& knownMask,
    const boolList& sourceKnownMask, const label cornerCount,
    const bool distanceWeighted, const word& context)
{
    if (cornerCount <= 0)
    {
        return 0;
    }

    const pointField& pts = this->points();
    const labelListList& pCells = this->pointCells();
    const labelListList& cPoints = this->cellPoints();
    label localInterpolated = 0;

    forAll(knownMask, pointI)
    {
        if (knownMask[pointI])
        {
            continue;
        }

        DynamicList<label> candidates;
        const labelList& localCells = pCells[pointI];

        forAll(localCells, cellI)
        {
            const label cI = localCells[cellI];
            const labelList& cellPts = cPoints[cI];

            forAll(cellPts, cpI)
            {
                const label corner = cellPts[cpI];
                if (corner == pointI || !sourceKnownMask[corner])
                {
                    continue;
                }

                bool duplicate = false;
                forAll(candidates, candI)
                {
                    if (candidates[candI] == corner)
                    {
                        duplicate = true;
                        break;
                    }
                }

                if (!duplicate)
                {
                    candidates.append(corner);
                }
            }
        }

        if (candidates.size() < cornerCount)
        {
            continue;
        }

        List<label> order(candidates.size(), -1);
        List<scalar> dist(candidates.size(), GREAT);
        forAll(candidates, candI)
        {
            order[candI] = candI;
            dist[candI] = mag(pts[pointI] - pts[candidates[candI]]);
        }

        std::sort(order.begin(), order.end(),
            [&](const label a, const label b) { return dist[a] < dist[b]; });

        // Select corners by minimizing geometric centroid mismatch at this
        // point (within nearest local candidates), which is more robust than
        // blindly taking the first N nearest corners on distorted AMR meshes.
        const label maxSearchCandidates =
            std::min<label>(order.size(), cornerCount == 2 ? 16 : 12);

        DynamicList<label> corners;
        scalar bestCentroidErr = GREAT;
        scalar bestAvgDist = GREAT;

        if (maxSearchCandidates >= cornerCount)
        {
            DynamicList<label> trial(cornerCount);

            auto updateBest = [&](const DynamicList<label>& labels) {
                point centroid = point::zero;
                scalar avgDist = 0.0;

                forAll(labels, li)
                {
                    const label corner = labels[li];
                    centroid += pts[corner];
                    avgDist += mag(pts[pointI] - pts[corner]);
                }

                centroid /= scalar(labels.size());
                avgDist /= scalar(labels.size());

                const scalar centroidErr = mag(pts[pointI] - centroid);

                if (centroidErr < bestCentroidErr - SMALL ||
                    (mag(centroidErr - bestCentroidErr) <= SMALL &&
                        avgDist < bestAvgDist))
                {
                    bestCentroidErr = centroidErr;
                    bestAvgDist = avgDist;
                    corners.clear();
                    forAll(labels, li)
                    {
                        corners.append(labels[li]);
                    }
                }
            };

            auto chooseCorners = [&](const auto& self, const label start,
                                     const label need) -> void {
                if (need == 0)
                {
                    updateBest(trial);
                    return;
                }

                const label maxI = maxSearchCandidates - need;
                for (label rankI = start; rankI <= maxI; ++rankI)
                {
                    trial.append(candidates[order[rankI]]);
                    self(self, rankI + 1, need - 1);
                    trial.remove();
                }
            };

            chooseCorners(chooseCorners, 0, cornerCount);
        }

        // Safety fallback: retain legacy nearest-corner selection if no
        // valid combination was found in the bounded search.
        if (corners.size() != cornerCount)
        {
            corners.clear();
            for (label pickI = 0; pickI < cornerCount; ++pickI)
            {
                corners.append(candidates[order[pickI]]);
            }
        }

        if (corners.size() != cornerCount)
        {
            continue;
        }

        scalarField weights(cornerCount, 1.0 / scalar(cornerCount));

        if (cornerCount == 2 && distanceWeighted)
        {
            const scalar d0 = mag(pts[pointI] - pts[corners[0]]);
            const scalar d1 = mag(pts[pointI] - pts[corners[1]]);
            const scalar denom = d0 + d1;

            if (denom > VSMALL)
            {
                weights[0] = d1 / denom;
                weights[1] = d0 / denom;
            }
            else
            {
                weights[0] = 0.5;
                weights[1] = 0.5;
            }
        }
        else if (cornerCount == 4 || cornerCount == 8)
        {
            // For hex-like AMR topology (face/cell midpoint creation), equal
            // corner weights preserve straight-line/plane consistency better.
            weights = 1.0 / scalar(cornerCount);
        }
        else
        {
            scalar sumInv = 0.0;
            forAll(corners, cI)
            {
                const scalar d = mag(pts[pointI] - pts[corners[cI]]);
                const scalar inv = 1.0 / max(d, refinementInterpTolerance_);
                weights[cI] = inv;
                sumInv += inv;
            }

            if (sumInv > VSMALL)
            {
                weights /= sumInv;
            }
            else
            {
                weights = 1.0 / scalar(cornerCount);
            }
        }

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            vector value = vector::zero;
            forAll(corners, cI)
            {
                value += weights[cI] * modeShapes_[modeI][corners[cI]];
            }
            modeShapes_[modeI][pointI] = value;
        }

        knownMask[pointI] = true;
        ++localInterpolated;
    }

    if (Pstream::master() && localInterpolated > 0)
    {
        Info << "Refinement interpolation (" << context << ") updated "
             << localInterpolated << " local points using " << cornerCount
             << "-corner topology rule." << endl;
    }

    return localInterpolated;
}

void fastDynamicFvMesh::ensureModeShapesForCurrentMesh(const mapPolyMesh& mpm)
{
    if (!meshRefinementSupport_ || nMode_ <= 0 || modeShapes_.size() == 0)
    {
        return;
    }

    const label newNPoints = this->points().size();
    const label oldNPoints = modeShapes_[0].size();
    const bool shrinkingOrEqualTopology = (newNPoints <= oldNPoints);
    const bool localPointCountChanged = (newNPoints != oldNPoints);
    const bool anyPointCountChanged =
        returnReduce(label(localPointCountChanged), sumOp<label>()) > 0;

    if (!anyPointCountChanged)
    {
        return;
    }

    const bool anyRefinement =
        returnReduce(label(newNPoints > oldNPoints), sumOp<label>()) > 0;

    const labelList& pointMap = mpm.pointMap();
    const labelList& reversePointMap = mpm.reversePointMap();

    if (pointMap.size() != newNPoints)
    {
        FatalErrorInFunction << "pointMap size (" << pointMap.size()
                             << ") does not match new mesh point count ("
                             << newNPoints << ")." << exit(FatalError);
    }

    List<vectorField> oldModeShapes(modeShapes_);
    modeShapes_.setSize(nMode_);

    forAll(modeShapes_, modeI)
    {
        modeShapes_[modeI].setSize(newNPoints, vector::zero);
    }

    boolList knownMask(newNPoints, false);
    label localDeferredAddedPoints = 0;

    for (label pointI = 0; pointI < newNPoints; ++pointI)
    {
        const label oldPoint = pointMap[pointI];
        if (oldPoint < 0 || oldPoint >= oldNPoints)
        {
            continue;
        }

        if (anyRefinement)
        {
            if (oldPoint >= reversePointMap.size())
            {
                FatalErrorInFunction
                    << "reversePointMap size (" << reversePointMap.size()
                    << ") cannot address old point " << oldPoint << "."
                    << exit(FatalError);
            }

            // During refinement, pointMap can assign newly-created points to
            // a single master old point. Those entries must not be treated as
            // fully mapped values; defer them to topology interpolation.
            if (reversePointMap[oldPoint] != pointI)
            {
                ++localDeferredAddedPoints;
                continue;
            }
        }

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            modeShapes_[modeI][pointI] = oldModeShapes[modeI][oldPoint];
        }

        knownMask[pointI] = true;
    }

    const label globalDeferredAddedPoints =
        returnReduce(localDeferredAddedPoints, sumOp<label>());
    if (Pstream::master() && globalDeferredAddedPoints > 0)
    {
        Info << "Deferred " << globalDeferredAddedPoints
             << " refinement-added points from pointMap master labels to "
             << "topology interpolation." << endl;
    }

    // Prefer direct topology lineage for refinement-added points.
    // mapPolyMesh stores source old points for new points in
    // pointsFromPointsMap. Averaging source modal values preserves geometric
    // continuity better than a generic nearest-corner search.
    label localLineageMapped = 0;
    const List<objectMap>& pointsFromPoints = mpm.pointsFromPointsMap();

    forAll(pointsFromPoints, mapI)
    {
        const objectMap& pointMapEntry = pointsFromPoints[mapI];
        const label newPoint = pointMapEntry.index();

        if (newPoint < 0 || newPoint >= newNPoints || knownMask[newPoint])
        {
            continue;
        }

        const labelList& masterPoints = pointMapEntry.masterObjects();
        if (masterPoints.size() == 0)
        {
            continue;
        }

        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            modeShapes_[modeI][newPoint] = vector::zero;
        }

        label validMasterCount = 0;
        forAll(masterPoints, masterI)
        {
            const label oldPoint = masterPoints[masterI];
            if (oldPoint < 0 || oldPoint >= oldNPoints)
            {
                continue;
            }

            ++validMasterCount;
            for (label modeI = 0; modeI < nMode_; ++modeI)
            {
                modeShapes_[modeI][newPoint] += oldModeShapes[modeI][oldPoint];
            }
        }

        if (validMasterCount == 0)
        {
            continue;
        }

        const scalar invMasterCount = 1.0 / scalar(validMasterCount);
        for (label modeI = 0; modeI < nMode_; ++modeI)
        {
            modeShapes_[modeI][newPoint] *= invMasterCount;
        }

        knownMask[newPoint] = true;
        ++localLineageMapped;
    }

    const label globalLineageMapped =
        returnReduce(localLineageMapped, sumOp<label>());
    if (Pstream::master() && globalLineageMapped > 0)
    {
        Info << "Topology lineage mapped " << globalLineageMapped
             << " refinement-added points from mapPolyMesh::pointsFromPoints."
             << endl;
    }

    syncMappedModeShapes(knownMask);

    label localKnownAfterSync = 0;
    forAll(knownMask, pointI)
    {
        if (knownMask[pointI])
        {
            ++localKnownAfterSync;
        }
    }

    if (localKnownAfterSync == 0 && newNPoints > 0)
    {
        FatalErrorInFunction
            << "No mapped old points found on processor " << Pstream::myProcNo()
            << " while processing topology change; "
            << "cannot interpolate modal shapes on refined mesh."
            << exit(FatalError);
    }

    if (anyRefinement)
    {
        label globalKnownPrev =
            returnReduce(localKnownAfterSync, sumOp<label>());

        for (label passI = 0; passI < 4; ++passI)
        {
            const boolList sourceKnownMask(knownMask);

            if (!shrinkingOrEqualTopology)
            {
                interpolateModeShapesFromKnown(
                    knownMask, sourceKnownMask, 2, true, "edge midpoint");
                interpolateModeShapesFromKnown(
                    knownMask, sourceKnownMask, 4, false, "face midpoint");
                interpolateModeShapesFromKnown(
                    knownMask, sourceKnownMask, 8, false, "cell midpoint");
            }

            syncMappedModeShapes(knownMask);

            label localKnownNow = 0;
            forAll(knownMask, pointI)
            {
                if (knownMask[pointI])
                {
                    ++localKnownNow;
                }
            }

            const label globalKnownNow =
                returnReduce(localKnownNow, sumOp<label>());

            if (globalKnownNow <= globalKnownPrev)
            {
                break;
            }

            globalKnownPrev = globalKnownNow;
        }
    }

    syncMappedModeShapes(knownMask);

    label localUnknown = 0;
    forAll(knownMask, pointI)
    {
        if (!knownMask[pointI])
        {
            ++localUnknown;
        }
    }

    if (localUnknown > 0)
    {
        FatalErrorInFunction
            << (shrinkingOrEqualTopology
                       ? "Topology mapping failed for "
                       : "Refinement interpolation failed for ")
            << localUnknown << " local points on processor "
            << Pstream::myProcNo() << " points."
            << (shrinkingOrEqualTopology
                       ? " Some surviving points were not mapped by pointMap."
                       : " The 2/4/8 topology-corner rules could not classify "
                         "all "
                         "added points. Adjust refinement controls or "
                         "refinementInterpTolerance.")
            << exit(FatalError);
    }
}

void fastDynamicFvMesh::cacheModelPointers()
{
    if (modelPointersCached_)
    {
        return;
    }

    typedef incompressible::turbulenceModel icoTurbModel;
    typedef compressible::turbulenceModel cmpTurbModel;

    if (this->foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        icoTurbPtr_ =
            &this->lookupObject<icoTurbModel>(icoTurbModel::propertiesName);
    }

    if (this->foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        cmpTurbPtr_ =
            &this->lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);
    }

    if (this->foundObject<fluidThermo>(fluidThermo::dictName))
    {
        fluidThermoPtr_ =
            &this->lookupObject<fluidThermo>(fluidThermo::dictName);
    }

    if (this->foundObject<transportModel>("transportProperties"))
    {
        laminarTransportPtr_ =
            &this->lookupObject<transportModel>("transportProperties");
    }

    modelPointersCached_ = true;
}

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

// patchDensity(patchi, defaultRho)
// Purpose: return a scalarField of density values corresponding to the
// requested patch. Behavior:
//  - If a volScalarField named "rho" exists in the case, the function returns a
//  copy of its boundary field for patchi
//    which may contain per-face or per-cell densities produced by the solver.
//  - Otherwise, a scalarField filled with defaultRho is returned with length
//  equal to the number of faces in the patch.
// Rationale: kinematic pressure fields (p/rho) require a density to convert
// back to dimensional pressure for force computations.
tmp<scalarField> fastDynamicFvMesh::patchDensity(
    const label patchi, const scalar defaultRho) const
{
    if (this->foundObject<volScalarField>("rho"))
    {
        const auto& rhoField = this->lookupObject<volScalarField>("rho");
        return tmp<scalarField>(
            new scalarField(rhoField.boundaryField()[patchi]));
    }

    return tmp<scalarField>(
        new scalarField(this->boundaryMesh()[patchi].size(), defaultRho));
}

// devRhoReff(gradUp, patchi, defaultRho)
// Purpose: compute the deviatoric effective stress tensor for wall-shear
// evaluation on a patch. Implementation notes:
//  - Prefer available turbulence/transport models in order to obtain an
//  effective kinematic viscosity (nuEff or muEff).
//  - When density-dependent viscosity or kinematic-pressure correction is
//  required, obtain per-face density with patchDensity().
//  - The returned tensor field represents -rho * nu_eff * devTwoSymm(gradUp),
//  i.e., the deviatoric viscous stress contribution per face.
//  - This helper centralizes model queries so calcModalForces can obtain shear
//  contributions regardless of solver setup.
tmp<symmTensorField> fastDynamicFvMesh::devRhoReff(const tensorField& gradUp,
    const label patchi, const scalar defaultRho) const
{
    if (icoTurbPtr_)
    {
        const auto& turb = *icoTurbPtr_;

        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>(new symmTensorField(
            -tRho() * turb.nuEff(patchi) * devTwoSymm(gradUp)));
    }
    else if (cmpTurbPtr_)
    {
        const auto& turb = *cmpTurbPtr_;

        return tmp<symmTensorField>(
            new symmTensorField(-turb.muEff(patchi) * devTwoSymm(gradUp)));
    }
    else if (fluidThermoPtr_)
    {
        const auto& thermo = *fluidThermoPtr_;

        return tmp<symmTensorField>(
            new symmTensorField(-thermo.mu(patchi) * devTwoSymm(gradUp)));
    }
    else if (laminarTransportPtr_)
    {
        const auto& laminarT = *laminarTransportPtr_;

        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>(new symmTensorField(
            -tRho() * laminarT.nu(patchi) * devTwoSymm(gradUp)));
    }
    else
    {
        IOdictionary transportProperties(
            IOobject("transportProperties", this->time().constant(), *this,
                IOobject::MUST_READ, IOobject::NO_WRITE));

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);
        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>(
            new symmTensorField(-tRho() * nu.value() * devTwoSymm(gradUp)));
    }
}

// calcModalForces()
// Purpose: compute per-mode generalized forces by projecting fluid pressure and
// wall-shear traction onto mode shapes. Detailed flow:
//  - Zero modal accumulators (modeForce_, modePressureForce_, modeShearForce_).
//  - Fetch the pressure field (supports dimensional pressure or kinematic
//  pressure). If kinematic, determine density (rhoRef or rho field).
//  - Optionally open per-processor face diagnostics CSV (append mode) and write
//  header if first time.
//  - If velocity field 'U' exists, compute its gradient once (fvc::grad(U)) to
//  evaluate viscous stress.
//  - For each configured FSI patch:
//      * Verify patch exists and obtain its face areas and centres.
//      * Compute per-face pressure force = pressure * areaVec (scaled by
//      density for kinematic p).
//      * Compute deviatoric viscous stress and then shear force; remove normal
//      component to leave tangential shear only.
//      * Interpolate mode shape to face by averaging shape vectors at face
//      points.
//      * Project pressure and shear vectors onto the shape vector (dot product)
//      and accumulate into pressure/shear/modal totals.
//      * Optionally write per-face diagnostics lines for the requested mode.
//  - Perform parallel reductions to combine per-processor modal accumulators
//  and (on master) print debug ranges.
void fastDynamicFvMesh::calcModalForces()
{
    cacheModelPointers();

    // Initialize forces
    modeForce_ = 0.0;
    modePressureForce_ = 0.0;
    modeShearForce_ = 0.0;

    if (!this->foundObject<volScalarField>(pressureFieldName_))
    {
        FatalErrorInFunction
            << "Required pressure field '" << pressureFieldName_
            << "' was not found. Configure 'pressureField' in dynamicMeshDict "
            << "or supply the field before fastDynamicFvMesh::update()."
            << exit(FatalError);
    }

    const volScalarField& p =
        this->lookupObject<volScalarField>(pressureFieldName_);
    const dimensionSet kinematicPressureDims(dimPressure / dimDensity);
    const bool dimensionalPressure = (p.dimensions() == dimPressure);
    const bool kinematicPressure = (p.dimensions() == kinematicPressureDims);

    if (!dimensionalPressure && !kinematicPressure)
    {
        FatalErrorInFunction << "Pressure field '" << pressureFieldName_
                             << "' has dimensions " << p.dimensions()
                             << ". Expected either " << dimPressure
                             << " (pressure) or " << kinematicPressureDims
                             << " (kinematic pressure)." << exit(FatalError);
    }

    scalar defaultRho = rhoRef_;

    if (kinematicPressure && defaultRho <= 0)
    {
        IOdictionary transportProperties(
            IOobject("transportProperties", this->time().constant(), *this,
                IOobject::MUST_READ, IOobject::NO_WRITE));

        if (transportProperties.found("rho"))
        {
            transportProperties.lookup("rho") >> defaultRho;
        }
    }

    if (kinematicPressure && defaultRho <= 0 &&
        !this->foundObject<volScalarField>("rho"))
    {
        FatalErrorInFunction
            << "Pressure field '" << pressureFieldName_
            << "' is kinematic, but no density information was found." << nl
            << "Provide 'rhoRef' in dynamicMeshDict, a 'rho' entry in "
            << "constant/transportProperties, or a volScalarField named 'rho'."
            << exit(FatalError);
    }

    tmp<volTensorField> tGradU;
    const volTensorField* gradUPtr = nullptr;

    if (this->foundObject<volVectorField>("U"))
    {
        const volVectorField& U = this->lookupObject<volVectorField>("U");
        tGradU = fvc::grad(U);
        gradUPtr = &tGradU();
    }
    else
    {
        WarningInFunction
            << "Velocity field 'U' not found; wall shear contribution "
            << "to modal forces will be skipped." << endl;
    }

    // Iterate over cached FSI patches
    forAll(fsiPatchIDs_, i)
    {
        const label patchID = fsiPatchIDs_[i];
        const polyPatch& pp = *fsiPolyPatches_[i];
        const fvPatchScalarField& pPatch = p.boundaryField()[patchID];
        const vectorField faceAreas(pp.faceAreas());
        scalarField& pressureScale = pressureScaleCache_[i];
        if (kinematicPressure)
        {
            pressureScale = patchDensity(patchID, defaultRho);
            pressureScaleInitialized_[i] = true;
        }
        else if (!pressureScaleInitialized_[i])
        {
            pressureScale = 1.0;
            pressureScaleInitialized_[i] = true;
        }

        const symmTensorField* devStressPtr = nullptr;
        tmp<symmTensorField> tDevStress;

        if (gradUPtr)
        {
            const tensorField& gradPatch = gradUPtr->boundaryField()[patchID];
            tDevStress = devRhoReff(
                gradPatch, patchID, defaultRho > 0 ? defaultRho : 1.0);
            devStressPtr = &tDevStress();
        }

        if (meshRefinementSupport_ && referenceFsiBuilt_ &&
            referenceFaceModeProjection_.size() == fsiPatches_.size())
        {
            const label nRefFaces = referenceFsiFaceAreas_[i].size();
            vectorField aggregatedPressure(nRefFaces, vector::zero);
            vectorField aggregatedShear(nRefFaces, vector::zero);
            const List<labelList>& faceToRefs = fsiFaceToReferenceFaces_[i];
            const List<scalarField>& faceToRefWeights =
                fsiFaceToReferenceWeights_[i];

            // Loop over current (possibly split) faces and aggregate to
            // reference faces
            forAll(pp, faceI)
            {
                if (faceI >= faceToRefs.size() ||
                    faceI >= faceToRefWeights.size())
                {
                    FatalErrorInFunction
                        << "Invalid refined-face mapping: patch '"
                        << fsiPatches_[i] << "', face " << faceI
                        << " has no reference mapping entry."
                        << exit(FatalError);
                }

                const labelList& refs = faceToRefs[faceI];
                const scalarField& weights = faceToRefWeights[faceI];

                if (refs.size() == 0 || refs.size() != weights.size())
                {
                    FatalErrorInFunction
                        << "Invalid refined/unrefined face mapping on patch '"
                        << fsiPatches_[i] << "', face " << faceI
                        << ": refs=" << refs.size()
                        << ", weights=" << weights.size() << exit(FatalError);
                }

                const vector& areaVec = faceAreas[faceI];
                const scalar areaMag = mag(areaVec);
                const vector faceNormal =
                    areaMag > VSMALL ? areaVec / areaMag : vector::zero;

                const vector pressureForce =
                    (pressureScale[faceI] * pPatch[faceI] - pRef_) * areaVec;
                vector shearForce = vector::zero;

                if (devStressPtr)
                {
                    shearForce = areaVec & (*devStressPtr)[faceI];
                    if (areaMag > VSMALL)
                    {
                        shearForce -= faceNormal * (faceNormal & shearForce);
                    }
                }

                const vector pressureTraction =
                    areaMag > VSMALL ? pressureForce / areaMag : vector::zero;
                const vector shearTraction =
                    areaMag > VSMALL ? shearForce / areaMag : vector::zero;

                forAll(refs, refI)
                {
                    const label refFaceI = refs[refI];
                    if (refFaceI < 0 || refFaceI >= nRefFaces)
                    {
                        FatalErrorInFunction
                            << "Invalid refined/unrefined mapping: patch '"
                            << fsiPatches_[i] << "', face " << faceI
                            << " -> ref " << refFaceI
                            << ", nRefFaces=" << nRefFaces << exit(FatalError);
                    }

                    const scalar w = max(weights[refI], 0.0);
                    const scalar refArea = referenceFsiFaceAreas_[i][refFaceI];
                    aggregatedPressure[refFaceI] +=
                        pressureTraction * (w * refArea);
                    aggregatedShear[refFaceI] += shearTraction * (w * refArea);
                }
            }

            for (label m = 0; m < nMode_; ++m)
            {
                const vectorField& refShape =
                    referenceFaceModeProjection_[i][m];
                forAll(refShape, refFaceI)
                {
                    const scalar pressureModeForce =
                        aggregatedPressure[refFaceI] & refShape[refFaceI];
                    const scalar shearModeForce =
                        aggregatedShear[refFaceI] & refShape[refFaceI];

                    modePressureForce_[m] += pressureModeForce;
                    modeShearForce_[m] += shearModeForce;
                    modeForce_[m] += pressureModeForce + shearModeForce;
                }
            }
        }
        else
        {
            // Loop over faces (legacy path)
            forAll(pp, faceI)
            {
                const vector& areaVec = faceAreas[faceI];
                const vector pressureForce =
                    (pressureScale[faceI] * pPatch[faceI] - pRef_) * areaVec;
                vector shearForce = vector::zero;

                if (devStressPtr)
                {
                    shearForce = areaVec & (*devStressPtr)[faceI];
                    const scalar areaMag = mag(areaVec);

                    if (areaMag > VSMALL)
                    {
                        const vector faceNormal = areaVec / areaMag;
                        shearForce -= faceNormal * (faceNormal & shearForce);
                    }
                }

                for (label m = 0; m < nMode_; ++m)
                {
                    const vector& shapeFace = faceModeProjection_[i][m][faceI];

                    const scalar pressureModeForce = pressureForce & shapeFace;
                    const scalar shearModeForce = shearForce & shapeFace;

                    modePressureForce_[m] += pressureModeForce;
                    modeShearForce_[m] += shearModeForce;
                    modeForce_[m] += pressureModeForce + shearModeForce;
                }
            }
        }
    }

    // Parallel reduction (all-reduce style)
    forAll(modeForce_, i)
    {
        reduce(modeForce_[i], sumOp<scalar>());
        reduce(modePressureForce_[i], sumOp<scalar>());
        reduce(modeShearForce_[i], sumOp<scalar>());
    }
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

// solveStructuralDynamics(dt)
// Purpose: integrate modal equations of motion and update modeState_ for
// displacement, velocity, and acceleration. Method and notes:
//  - A single-degree-of-freedom modal equation is solved per mode: M*x_dd +
//  C*x_d + K*x = F.
//  - Current implementation uses Wilson-Theta single-step integration with
//  theta_ (default 1.4) for stability.
//  - Mass and damping are read from dynamicMeshDict as modalMass and modalDamp.
//  - If not provided, Mass defaults to 1.0 and Damp to 0.0.
//  - Steps per mode:
//      * Extract last-step state (disLast, velLast, accLast) and
//      current/previous modal forces (F_this, F_last).
//      * Compute omega = 2*pi*freq and build an effective stiffness term
//      including Wilson-Theta contributions.
//      * Form the predicted load (includes theta-weighted force and history
//      terms from Mass and Damp).
//      * Solve for displacement at theta (disTheta) and recover acc, vel, dis
//      for the current time step.
//  - Results are stored in modeState_[i] as (disp, vel, acc).
void fastDynamicFvMesh::solveStructuralDynamics(scalar dt)
{
    // const scalar Mass = 1.0;
    // const scalar Damp = 0.0;

    for (label i = 0; i < nMode_; ++i)
    {
        scalar Mass = modeMass_[i];
        scalar Damp = modeDamp_[i];

        scalar dis = 0, vel = 0, acc = 0;

        scalar disLast = modeState0_[i].x();
        scalar velLast = modeState0_[i].y();
        scalar accLast = modeState0_[i].z();

        scalar F_this = modeForce_[i];
        scalar F_last = modeForce0_[i];
        scalar freq = modeFreq_[i];

        // Debug
        /*
        if (Pstream::master() && i==0)
        {
             Info<< "Mode 0: Freq=" << freq << " F_this=" << F_this << "
        F_last=" << F_last << endl;
        }
        */

        scalar omega = 2.0 * constant::mathematical::pi * freq;

        // Ensure non-zero denominator logic
        scalar effectiveK = 6.0 * Mass / sqr(theta_ * dt) +
            3.0 * Damp / (theta_ * dt) + Mass * sqr(omega);

        scalar load = F_last + theta_ * (F_this - F_last);
        load += Mass *
            (6.0 / sqr(theta_ * dt) * disLast + 6.0 / (theta_ * dt) * velLast +
                2.0 * accLast);
        load += Damp *
            (3.0 / (theta_ * dt) * disLast + 2.0 * velLast +
                0.5 * theta_ * dt * accLast);

        scalar disTheta = load / effectiveK;

        scalar tau = theta_ * dt;
        scalar term1 = 6.0 / (sqr(tau) * theta_) * (disTheta - disLast);
        scalar term2 = 6.0 / (sqr(theta_) * dt) * velLast;
        scalar term3 = (1.0 - 3.0 / theta_) * accLast;
        acc = term1 - term2 + term3;

        vel = velLast + 0.5 * dt * (acc + accLast);
        dis = disLast + dt * velLast + sqr(dt) / 6.0 * (acc + 2.0 * accLast);

        /*
        if (Pstream::master() && i==0)
        {
             Info << "Mode 0 Step:"
                  << " dt=" << dt
                  << " theta=" << theta_
                  << " F_this=" << F_this
                  << " load=" << load
                  << " effK=" << effectiveK
                  << " disTheta=" << disTheta
                  << " dis=" << dis
                  << " acc=" << acc << endl;
        }
        */

        modeState_[i] = vector(dis, vel, acc);
    }
}

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
