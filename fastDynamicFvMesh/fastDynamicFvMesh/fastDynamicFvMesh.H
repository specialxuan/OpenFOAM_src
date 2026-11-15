/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    FDM-AMR runtime facade. Keeps the historical OpenFOAM type name
    'fastDynamicFvMesh' for dictionary and case compatibility while the FDM
    core lives in FdmCore. This class owns only the optional runtime AMR
    topology machinery; the AMR algorithm itself is inherited from
    dynamicRefineFvMesh for now and will move behind an AMR provider in the
    next refactor stage.

\*---------------------------------------------------------------------------*/

#ifndef fastDynamicFvMesh_H
#define fastDynamicFvMesh_H

#include "dynamicRefineFvMesh.H"
#include "fdmCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fastDynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class fastDynamicFvMesh
:
    public dynamicRefineFvMesh,
    public FdmCore
{
    // Protected Data (runtime AMR topology only)

    //- Enable runtime topology refinement (dynamicRefineFvMeshCoeffs present)
    bool runtimeRefinementEnabled_;

    //- Print per-topology-change mode-shape mapping diagnostics
    bool refinementMappingDiagnostics_;

    //- Optional runtime-AMR minimum child-cell volume floor (<=0 disables)
    scalar refinementMinCellVolume_;

    //- Optional runtime-AMR minimum child-edge length floor (<=0 disables)
    scalar refinementMinEdgeLength_;

    //- Use |grad(field)| as runtime-AMR indicator instead of raw field values
    bool refinementUseGradIndicator_;

    //- Source scalar field name used to build |grad(field)| indicator
    word refinementGradIndicatorField_;

    //- Optional registered volScalarField containing precomputed |grad(field)|
    word refinementGradIndicatorMagnitudeField_;

    //- Cached runtime-AMR cell indicator built from |grad(field)|
    mutable scalarField refinementGradIndicatorCellCache_;

    //- Cached runtime-AMR point indicator built from |grad(field)|
    mutable scalarField refinementGradIndicatorPointCache_;

    //- Time index for cached runtime-AMR gradient indicator
    mutable label refinementGradIndicatorCacheTimeIndex_;

    //- Cell count for cached runtime-AMR gradient indicator
    mutable label refinementGradIndicatorCacheNCells_;

    //- Point count for cached runtime-AMR gradient indicator
    mutable label refinementGradIndicatorCacheNPoints_;

    //- Source field name for cached runtime-AMR gradient indicator
    mutable word refinementGradIndicatorCacheSourceField_;

    //- External magnitude field name for cached runtime-AMR gradient indicator
    mutable word refinementGradIndicatorCacheMagnitudeField_;

    //- Whether cached runtime-AMR gradient indicator is valid
    mutable bool refinementGradIndicatorCacheValid_;

    //- Enable bottom-layer runtime-AMR decision diagnostics
    bool amrLayerDiagnosticsEnabled_;

    //- Target z coordinate for bottom-layer AMR diagnostics
    scalar amrLayerDiagnosticsZ_;

    //- Diagnostic layer half-thickness; <=0 selects nearest z layer
    scalar amrLayerDiagnosticsThickness_;

    //- Face-neighbor expansion depth around selected/candidate cells
    label amrLayerDiagnosticsNeighborDepth_;

    //- Maximum per-rank diagnostic cell rows per AMR call; <=0 disables limit
    label amrLayerDiagnosticsMaxRows_;

    //- File prefix for AMR layer diagnostic CSV outputs
    word amrLayerDiagnosticsFilePrefix_;

    //- Header state for AMR layer summary output
    mutable bool amrLayerDiagnosticsSummaryHeaderWritten_;

    //- Header state for AMR layer per-rank cell output
    mutable bool amrLayerDiagnosticsCellHeaderWritten_;

    //- Cached AMR layer candidate mask for current updateTopology call
    mutable bitSet amrLayerDiagnosticsCandidateCell_;

    //- Cached AMR layer indicator values for current updateTopology call
    mutable scalarField amrLayerDiagnosticsIndicatorCell_;

    //- Cached lower refine threshold for current updateTopology call
    mutable scalar amrLayerDiagnosticsLower_;

    //- Cached upper refine threshold for current updateTopology call
    mutable scalar amrLayerDiagnosticsUpper_;

    //- Cached local summary counts for current updateTopology call
    mutable label amrLayerDiagnosticsLayerCells_;
    mutable label amrLayerDiagnosticsNearbyCells_;
    mutable label amrLayerDiagnosticsCandidates_;
    mutable label amrLayerDiagnosticsSelected_;
    mutable label amrLayerDiagnosticsBelowLower_;
    mutable label amrLayerDiagnosticsAboveUpper_;
    mutable label amrLayerDiagnosticsRejectedSizeFloor_;
    mutable label amrLayerDiagnosticsRejectedBase_;

    //- Cached local unrefine point count for current updateTopology call
    mutable label amrLayerDiagnosticsUnrefinePoints_;

    //- True when an AMR layer summary is ready to write
    mutable bool amrLayerDiagnosticsPendingSummary_;

    //- Enable per-runtime-AMR-call timing log output
    bool trackAmrTimingEnabled_;

    // Protected Member Functions (runtime AMR topology only)

    //- Read runtime-AMR dictionary entries
    virtual void readAmrControls(
        const dictionary& fdmDict,
        const dictionary& dynamicMeshDict);

    //- Per-cell minimum edge length (used by runtime refinement size floor)
    scalar cellMinEdgeLength(const label celli) const;

    //- Whether a candidate parent cell satisfies refinement size-floor controls
    bool cellPassesRefinementSizeFloor(const label celli) const;

    //- Clear cached runtime-AMR gradient indicator fields
    void clearRefinementGradIndicatorCache() const;

    //- Ensure cached runtime-AMR gradient indicator fields are current
    void ensureRefinementGradIndicatorCache() const;

    //- Runtime-AMR indicator at cell centers (raw field or |grad(field)|)
    scalarField refinementIndicatorCellField(
        const scalarField& fallbackCellField) const;

    //- Runtime-AMR indicator at points (raw field or |grad(field)|)
    scalarField refinementIndicatorPointField(
        const scalarField& fallbackPointField) const;

    //- Select refinement candidates with optional child-size floor constraints
    virtual void selectRefineCandidates(const scalar lowerRefineLevel,
        const scalar upperRefineLevel, const scalarField& vFld,
        bitSet& candidateCell) const;

    //- Select refinement cells with optional child-size floor constraints
    virtual labelList selectRefineCells(const label maxCells,
        const label maxRefinement, const bitSet& candidateCell) const;

    //- Select unrefine points using the configured runtime-AMR indicator
    virtual labelList selectUnrefinePoints(const scalar unrefineLevel,
        const bitSet& markedCell, const scalarField& pFld) const;

    //- Reset cached AMR layer diagnostic state before/after an AMR call
    void resetAmrLayerDiagnosticsState() const;

    //- Write bottom-layer AMR decision diagnostics for selected/candidate cells
    void captureAmrLayerDiagnostics(const scalar lowerRefineLevel,
        const scalar upperRefineLevel, const label maxCells,
        const label maxRefinement, const bitSet& candidateCell,
        const labelList& selectedCells) const;

    //- Flush cached bottom-layer AMR summary after updateTopology completes
    void writeAmrLayerDiagnosticsSummary() const;

    //- Print AMR/refinement mode-shape mapping diagnostics when enabled
    void logModeShapeMappingSummary(const label oldNPoints,
        const label newNPoints, const bool anyRefinement,
        const bool shrinkingOrEqualTopology, const label localDirectMapped,
        const label localDeferredAddedPoints, const label localLineageMapped,
        const label localEdgeInterpolated, const label localFaceInterpolated,
        const label localCellInterpolated, const label localUnknown) const;

    //- Keep mode-shape storage consistent with current mesh topology
    void ensureModeShapesForCurrentMesh(const mapPolyMesh& mpm);

    //- Hook override: run the runtime topology engine for a new time step
    virtual void updateRefinementTopology();

    //- Hook override: runtime topology engine present and active
    virtual bool refinementTopologyActive() const;

    //- Hook override: AMR point levels from the OpenFOAM refinement engine
    virtual const labelList& amrPointLevel() const;

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

    //- Update topology, then the FDM core
    virtual bool update();

    //- Update custom caches/shapes after topology changes
    virtual void updateMesh(const mapPolyMesh& mpm);

    //- Write the mesh plus AMR checkpoint side files
    virtual bool writeObject
    (
        IOstreamOption streamOpt,
        const bool writeOnProc
    ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
