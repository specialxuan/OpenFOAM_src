/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Fast Dynamic Mesh method based on modal superposition.

\*---------------------------------------------------------------------------*/

#ifndef fastDynamicFvMesh_H // Prevent Multiple Inclusions
#define fastDynamicFvMesh_H

#include "dynamicRefineFvMesh.H"
#include "fieldTypes.H"
#include "uniformDimensionedFields.H"
#include "IOField.H"
#include "fluidThermo.H"
#include "transportModel.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class fastDynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class fastDynamicFvMesh : public dynamicRefineFvMesh
{
    // Private Data

    //- Number of modes
    label nMode_;

    //- Modal frequencies
    scalarField modeFreq_;

    //- Modal mass
    scalarField modeMass_;

    //- Modal damping ratios
    scalarField modeDamp_;

    //- Modal stiffness
    scalarField modeStiff_;

    //- Modal force (current time step)
    //  Now an IOField for restart capability
    IOField<scalar> modeForce_;

    //- Pressure contribution to modal force
    scalarField modePressureForce_;

    //- Viscous wall-shear contribution to modal force
    scalarField modeShearForce_;

    //- Modal force (previous time step)
    scalarField modeForce0_;

    //- Modal displacement (current time step)
    //  Layout: [disp, vel, acc] per mode
    //  Now an IOField for restart capability
    IOField<vector> modeState_;

    //- Modal displacement (previous time step)
    List<vector> modeState0_;

    //- Initial modal velocity
    scalarField initVelocity_;

    //- Relaxed modal displacement actually imposed on the mesh
    //  Now an IOField for restart capability
    IOField<scalar> appliedModeDisp_;

    //- Relaxed modal force from previous inner iteration
    //  Now an IOField for restart capability
    IOField<scalar> appliedModeForce_;

    //- Relative residual of modal force convergence
    scalar fsiResidual_;

    //- Mode shapes (displacement field per mode)
    //  Stores full field for simplicity
    //  modeShapes_[m][pI].x() gives x-displacement for mode m at point pI
    List<vectorField> modeShapes_;

    //- Wilson-Theta parameter
    scalar theta_;

    //- Point-mapping tolerance used for legacy CSV node matching
    scalar mappingTolerance_;

    //- Relaxation applied to structural modal displacement before moving mesh
    scalar couplingRelaxation_;

    //- Pressure field used for modal loading
    word pressureFieldName_;

    //- Reference density for kinematic pressure fields
    scalar rhoRef_;
    scalar pRef_;
    scalar maxDispChange_;

    //- FSI patch names
    wordRes fsiPatches_;

    //- Number of startup calls already treated as Wilson-theta initialisation
    label startupStepCount_;

    //- Last physical time index already processed by update()
    label lastUpdateTimeIndex_;

    //- Count of updates within current time step
    label updateCount_;

    //- Enable/disable modal diagnostics CSV output (disabled by default)
    bool writeDiagnosticsEnabled_;

    //- Enable/disable timing profiler for fluid/mesh CPU usage
    bool trackTimingEnabled_;

    //- CPU time at end of previous update() call
    scalar timingLastUpdateCpuTime_;

    //- Whether timingLastUpdateCpuTime_ is initialized
    bool timingHasLastUpdate_;

    //- Accumulated estimated fluid-solver CPU time
    scalar timingFluidCpuAccum_;

    //- Accumulated mesh-update CPU time
    scalar timingMeshCpuAccum_;

    //- Snapshot of accumulated fluid time at previous writeTime
    scalar timingFluidCpuAtLastWrite_;

    //- Snapshot of accumulated mesh time at previous writeTime
    scalar timingMeshCpuAtLastWrite_;

    //- Last physical time index for which global modal files were written
    label lastGlobalWriteTimeIndex_;

    //- Cached FSI patch ids for faster access during force assembly
    List<label> fsiPatchIDs_;

    //- Cached FSI patch references
    List<const polyPatch*> fsiPolyPatches_;

    //- Cached face-mode projection coefficients [patch][mode][face]
    List<List<vectorField>> faceModeProjection_;

    //- Reused storage of pressure scaling per cached patch
    List<scalarField> pressureScaleCache_;

    //- Whether pressure scale cache is initialized for each cached patch
    List<bool> pressureScaleInitialized_;

    //- Cache writeDiagnostics header status to avoid file existence checks
    bool diagnosticsHeaderWritten_;

    //- Cached model pointers for shear evaluation
    bool modelPointersCached_;
    const incompressible::turbulenceModel* icoTurbPtr_;
    const compressible::turbulenceModel* cmpTurbPtr_;
    const fluidThermo* fluidThermoPtr_;
    const transportModel* laminarTransportPtr_;

    //- Enable/disable prescribed structural nodal loading
    bool structuralForceEnabled_;

    //- Number of structural force files to read
    label nStructuralForces_;

    //- Structural force file prefix (e.g., StructForce -> StructForce1.csv)
    word structuralForceFilePrefix_;

    //- Geometric tolerance for mapping force-file coordinates to structure nodes
    scalar structuralTargetTolerance_;

    //- Structural mode files in mode/ directory
    word structNodeCoorFile_;
    word structNodeDispPrefix_;

    //- Structural nodes and mode shapes loaded from StructNode*.csv
    pointField structNodeCoords_;
    List<vectorField> structModeShapes_;

    //- Mapped structure node ids per force file
    List<labelList> structuralForceNodeIDs_;

    //- Time series (seconds) per force file
    List<scalarField> structuralForceTimes_;

    //- Nodal force vectors per force file (same indexing as structuralForceTimes_)
    List<vectorField> structuralForceValues_;

    //- Structural contribution to modal force (before relaxation)
    scalarField structuralModeForce_;

    //- Current force signal for diagnostics/logging
    scalar structuralForceSignal_;

    //- Enable refinement-aware modal interpolation and FSI force aggregation
    bool meshRefinementSupport_;

    //- Enable runtime topology refinement (dynamicRefineFvMeshCoeffs present)
    bool runtimeRefinementEnabled_;

    //- Geometric tolerance for refinement interpolation classification
    scalar refinementInterpTolerance_;

    //- Optional runtime-AMR minimum child-cell volume floor (<=0 disables)
    scalar refinementMinCellVolume_;

    //- Optional runtime-AMR minimum child-edge length floor (<=0 disables)
    scalar refinementMinEdgeLength_;

    //- Whether reference FSI face cache has been initialized
    bool referenceFsiBuilt_;

    //- Reference face-mode projection coefficients [patch][mode][referenceFace]
    List<List<vectorField>> referenceFaceModeProjection_;

    //- Mapping from current FSI patch faces to represented reference faces
    List<List<labelList>> fsiFaceToReferenceFaces_;

    //- Area weights per current face for represented reference faces
    List<List<scalarField>> fsiFaceToReferenceWeights_;

    //- Reference face areas captured on initial FSI topology
    List<scalarField> referenceFsiFaceAreas_;

    // Private Member Functions

    //- Read modal properties from dictionary and files
    void readControls();

    //- Read mode shapes from external files
    void readModeShapes();

    //- Read legacy runtime parameters from FluidPara.csv
    void readLegacyParameters(const fileName& modeDir);

    //- Read mode files on master
    void readModeFiles(
        const fileName& modeDir,
        List<point>& csvPoints,
        List<List<vector>>& csvShapes,
        label& nCsvNodes);

    //- Read structural node coordinates and mode shapes from StructNode*.csv
    void readStructuralModeFiles(
        const fileName& modeDir,
        pointField& structPoints,
        List<vectorField>& structShapes,
        label& nStructNodes,
        label& nStructModes);

    //- Read structural force files, map coordinates to structure nodes, cache series
    void readStructuralForceFiles(const fileName& modeDir);

    //- Build cached patch IDs, projection and topology-consistent face mapping
    void buildForceProjectionCaches(const mapPolyMesh* mpmPtr = nullptr);

    //- Synchronize mapped modal values on shared points (parallel consistency)
    void syncMappedModeShapes(boolList& mappedMask);

    //- Synchronize moved point positions on shared points (parallel consistency)
    void syncSharedPointPositions(pointField& pointValues) const;

    //- Interpolate mode shapes for currently-unknown points from known corners
    label interpolateModeShapesFromKnown(
        boolList& knownMask,
        const boolList& sourceKnownMask,
        const label cornerCount,
        const bool distanceWeighted,
        const word& context);

    //- Keep mode-shape storage consistent with current mesh topology
    void ensureModeShapesForCurrentMesh(const mapPolyMesh& mpm);

    //- Per-cell minimum edge length (used by runtime refinement size floor)
    scalar cellMinEdgeLength(const label celli) const;

    //- Whether a candidate parent cell satisfies refinement size-floor controls
    bool cellPassesRefinementSizeFloor(const label celli) const;

    //- Cache turbulence/thermo/transport pointers from object registry
    void cacheModelPointers();

    //- Select refinement candidates with optional child-size floor constraints
    virtual void selectRefineCandidates
    (
        const scalar lowerRefineLevel,
        const scalar upperRefineLevel,
        const scalarField& vFld,
        bitSet& candidateCell
    ) const;

    //- Select refinement cells with optional child-size floor constraints
    virtual labelList selectRefineCells
    (
        const label maxCells,
        const label maxRefinement,
        const bitSet& candidateCell
    ) const;

    //- Calculate modal forces from fluid pressure and shear stress
    void calcModalForces();

    //- Calculate structural modal force from prescribed nodal loading
    void calcStructuralModalForces();

    //- Linear interpolation of force time-series
    vector interpolateStructuralForce(
        const scalarField& times,
        const vectorField& values,
        const scalar t) const;

    //- Return patch density used for shear calculations
    tmp<scalarField> patchDensity(const label patchi, scalar defaultRho) const;

    //- Deviatoric effective stress tensor for wall shear evaluation
    tmp<symmTensorField> devRhoReff(
        const tensorField& gradUp, const label patchi, scalar defaultRho) const;

    //- Solve structural dynamics using Wilson-Theta method
    void solveStructuralDynamics(scalar dt);

    //- Write modal force/displacement diagnostics (optional)
    void writeDiagnostics();

    //- Write timing report (optional)
    void writeTimingReport();

    //- Disallow default bitwise copy construct
    fastDynamicFvMesh(const fastDynamicFvMesh&);

    //- Disallow default bitwise assignment
    void operator=(const fastDynamicFvMesh&);

public:
    //- Runtime type information
    TypeName("fastDynamicFvMesh");

    // Constructors

    //- Construct from IOobject
    explicit fastDynamicFvMesh(const IOobject& io);

    //- Destructor
    virtual ~fastDynamicFvMesh();

    // Member Functions

    //- Initialise all non-demand-driven data
    virtual bool init(const bool doInit);

    //- Update the mesh for both mesh motion and topology change
    virtual bool update();

    //- Update custom caches/shapes after topology changes
    virtual void updateMesh(const mapPolyMesh& mpm);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
