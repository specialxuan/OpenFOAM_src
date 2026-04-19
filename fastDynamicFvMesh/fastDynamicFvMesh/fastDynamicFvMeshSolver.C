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

} // End namespace Foam
