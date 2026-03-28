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

#include "dynamicFvMesh.H"
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

/*---------------------------------------------------------------------------*\
                       Class fastDynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class fastDynamicFvMesh : public dynamicFvMesh // Inherit from dynamicFvMesh
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

    //- Build cached patch IDs, face-mode projection and pressure scale storage
    void buildForceProjectionCaches();

    //- Cache turbulence/thermo/transport pointers from object registry
    void cacheModelPointers();

    //- Calculate modal forces from fluid pressure and shear stress
    void calcModalForces();

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

    //- Update the mesh for both mesh motion and topology change
    virtual bool update();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
