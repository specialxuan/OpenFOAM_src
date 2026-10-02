/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 1991-2008 OpenCFD Ltd.

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.
Class
    Foam::ausmPlusMFluxScheme

Description
    AUSM+M flux splitting scheme as described by:

    Chen, Shu-sheng & Cai, Fang-jie & Xue, Hai-chao & Wang, Ning & Yan, Chao.
    (2019). An improved AUSM-family scheme with robustness and accuracy for all
    Mach number flows. Applied Mathematical Modelling. 77.
    10.1016/j.apm.2019.09.005.

SourceFiles
    ausmPlusMFluxScheme.C

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa
    Maximilian Maigler
    Maximilian Mudra
        UniBw Munich, Germany

\*---------------------------------------------------------------------------*/

#include "ausmPlusMFluxScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcMeshPhi.H"
#include "fvcReconstruct.H"
#include "fvcSurfaceReconstruct.H"
#include "cellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    defineTypeNameAndDebug(ausmPlusMFluxScheme, 0);
    addToRunTimeSelectionTable(fluxScheme, ausmPlusMFluxScheme, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    ausmPlusMFluxScheme::ausmPlusMFluxScheme(
        const dictionary &dict,
        const fluidThermo &thermo,
        const volScalarField &rho,
        const volVectorField &U,
        const volVectorField &rhoU,
        const volScalarField &rhoE)
        : fluxScheme(typeName, dict),
          mesh_(U.mesh()),
          thermo_(thermo),
          rho_(rho),
          U_(U),
          rhoU_(rhoU),
          rhoE_(rhoE),
          dict_(dict),
          own_(surfaceScalarField(
              IOobject(
                  "own",
                  #if FOUNDATION >= 12
                  mesh_.time().name(),
                  #else
                  mesh_.time().timeName(),
                  #endif
                  mesh_),
              mesh_,
              dimensionedScalar("own", dimless, 1.0))),
          nei_(surfaceScalarField(
              IOobject(
                  "nei",
                  #if FOUNDATION >= 12
                  mesh_.time().name(),
                  #else
                  mesh_.time().timeName(),
                  #endif
                  mesh_),
              mesh_,
              dimensionedScalar("nei", dimless, -1.0)))
    //   ausm_dict_(mesh_.schemesDict().subDict("AUSM"))
    {
        Info << endl;
        Info << "------------------------------------------------------------------------------------------" << endl;
        Info << endl;
        Info << "AUSM+M flux scheme selected!" << endl;

        // sqrMachInf_ = ausm_dict_.lookupOrDefault<scalar>("sqrMachInf", 0.01);
        sqrMachInf_ = dict_.lookupOrDefault<scalar>("sqrMachInf", 0.01);
        // alpha0_ = ausm_dict_.lookupOrDefault<scalar>("alpha0", 3 / 16);
        alpha0_ = dict_.lookupOrDefault<scalar>("alpha0", 3 / 16);

        if (debug)
        {
            Info << "---> DEBUG MODE <---" << endl;
            Info << endl;
            Info << "sqrMachInf_: " << sqrMachInf_ << endl;
            Info << "alpha0_: " << alpha0_ << endl;
        };

        Info << endl;
        Info << "------------------------------------------------------------------------------------------" << endl;
        Info << endl;
    }

    // * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

    ausmPlusMFluxScheme::~ausmPlusMFluxScheme()
    {
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    void Foam::ausmPlusMFluxScheme::calcFlux(
        surfaceScalarField &phi,
        surfaceVectorField &phiUp,
        surfaceScalarField &phiEp,
        surfaceVectorField &Up)
    {
        Info << "AUSM+M: Updating fluxes!" << endl;

        // Cell surface values
        const surfaceScalarField surface_area("surface_area", mesh_.magSf());
        surfaceVectorField face_normal("face_normal", mesh_.Sf() / mesh_.magSf());

        // Velocity of owner and neighbour cell
        tmp<surfaceVectorField> u_L, u_R;
        fvc::surfaceReconstruct(U_, u_L, u_R, RECONSTRUCT_VECTOR_);
        tmp<surfaceScalarField> phi_L = u_L() & mesh_.Sf();
        tmp<surfaceScalarField> phi_R = u_R() & mesh_.Sf();

        // Flux relative to mesh movement
        volVectorField Urel(U_);
        if (mesh_.moving())
        {
            Urel = U_ - fvc::reconstruct(fvc::meshPhi(U_));
            fvc::makeRelative(phi_L.ref(), U_);
            fvc::makeRelative(phi_R.ref(), U_);
        }

        tmp<surfaceVectorField> Urel_L, Urel_R;
        fvc::surfaceReconstruct(Urel, Urel_L, Urel_R, RECONSTRUCT_VECTOR_);
        // Calaculate realtive tangential velocity
        const surfaceScalarField ut_L("ut_L", mag(Urel_L.ref() - (Urel_L.ref() & face_normal) * face_normal));
        const surfaceScalarField ut_R("ut_R", mag(Urel_R.ref() - (Urel_R.ref() & face_normal) * face_normal));
        // Magnitude of normal velocity
        const surfaceScalarField u_normal_L(phi_L.ref() / mesh_.magSf());
        const surfaceScalarField u_normal_R(phi_R.ref() / mesh_.magSf());
        phi_L.clear();
        phi_R.clear();

        // Density
        tmp<surfaceScalarField> rho_L, rho_R;
        fvc::surfaceReconstruct(rho_, rho_L, rho_R, RECONSTRUCT_SCALAR_);
        const surfaceScalarField rho_mean("rho_mean", 0.5 * (rho_L.ref() + rho_R.ref()));

        const surfaceScalarField u_normal_mean("u_normal_mean", 0.5 * (u_normal_L + u_normal_R));

        // Pressure of owner, neighbour cells
        tmp<surfaceScalarField> p_L, p_R;
        fvc::surfaceReconstruct(thermo_.p(), p_L, p_R, RECONSTRUCT_SCALAR_);
        const surfaceScalarField p_mean("p_mean", 0.5 * (p_L.ref() + p_R.ref()));

        // Critical acoustic velocity (Liou 2006)
        // tmp<volScalarField> H = rhoE_ / rho_ + thermo_.p() / rho_;
        // H->rename("H");
        tmp<volScalarField> Hrel = rhoE_ / rho_ + thermo_.p() / rho_;
        Hrel->rename("Hrel");

        if (mesh_.moving())
        {
            Hrel.ref() -= 0.5 * (U_ & U_);
            Hrel.ref() += 0.5 * (Urel & Urel);
        }
        Urel.clear();

        tmp<surfaceScalarField> h_L, h_R;
        fvc::surfaceReconstruct(Hrel, h_L, h_R, RECONSTRUCT_SCALAR_);
        Hrel.clear();

        // tmp<surfaceScalarField> ht_L, ht_R;
        // tmp<volScalarField> Ht = thermo_.Cp() * thermo_.T();
        // fvc::surfaceReconstruct(Ht, ht_L, ht_R, RECONSTRUCT_SCALAR_);
        // Ht.clear();

        tmp<volScalarField> gamma = thermo_.gamma();
        tmp<surfaceScalarField> gamma_L, gamma_R;
        fvc::surfaceReconstruct(gamma.ref(), gamma_L, gamma_R, RECONSTRUCT_SCALAR_);
        surfaceScalarField gamma_mean("gamma_mean", 0.5 * (gamma_L.ref() + gamma_R.ref()));
        gamma.clear();
        gamma_L.clear();
        gamma_R.clear();

        // calculate normal component of enthalpy
        const surfaceScalarField h_normal(
            "h_normal",
            // 0.5 * (h_L.ref() + 0.5 * (magSqr(u_L.ref()) - sqr(ut_L)) + h_R.ref() + 0.5 * (magSqr(u_R.ref()) - sqr(ut_R))));
            0.5 * (h_L.ref() - 0.5 * sqr(ut_L) + h_R.ref() - 0.5 * sqr(ut_R)));
        // 0.5 * (ht_L + 0.5 * (magSqr(u_L) - sqr(ut_L)) + ht_R + 0.5 * (magSqr(u_R) - sqr(ut_R))));

        const surfaceScalarField sqr_c_star(
            "sqr_c_star",
            2.0 * (gamma_mean - 1.0) / (gamma_mean + 1.0) * h_normal);

        const surfaceScalarField c_face(
            "c_face",
            pos0(u_normal_mean) * (sqr_c_star / max(mag(u_normal_L), sqrt(sqr_c_star))) +
                neg(u_normal_mean) * (sqr_c_star / max(mag(u_normal_R), sqrt(sqr_c_star))));

        // Compute Mach numbers of neighbouring cells as in Eq. 18
        const surfaceScalarField Mach_L("Mach_L", u_normal_L / c_face);
        const surfaceScalarField Mach_R("Mach_R", u_normal_R / c_face);
        const surfaceScalarField magnitude_Mach_L("magnitude_Mach_L", mag(Mach_L));
        const surfaceScalarField magnitude_Mach_R("magnitude_Mach_R", mag(Mach_R));

        // Eq. 15: Temporary Mach term for weight functions
        tmp<surfaceScalarField> mach = min(1.0, max(magnitude_Mach_L, magnitude_Mach_R));
        const surfaceScalarField f("f", 0.5 * (1 - cos(PI_ * mach.ref())));
        mach.clear();
        // Eq. 20
        // 0.0 for low mach regime; 1.0 for high mach regime
        const surfaceScalarField f_0("f_0", min(1.0, max(f, sqrMachInf_)));

        // alpha is constant 3/16 (see remark under Eq. 27)
        const scalar alpha = alpha0_;
        // AUSM+up:
        // const surfaceScalarField alpha = alpha0_ * (-4.0 + 5.0 * sqr(f_0));

        // PRESSURE FLUX
        // pressure based sensing function min(p_L / p_R, p_R / p_L) of adjacent faces
        // Eq. 25
        const surfaceScalarField h = calc_h(p_R.ref(), p_L.ref()); // small in shock region

        // Eq. 25
        // 1.0 in shock region; 0.0 without shock
        const surfaceScalarField g("g", 0.5 * (1.0 + cos(PI_ * h)));

        // Eq. 6 psi_L and psi_R
        tmp<surfaceScalarField> p_plus_L =
            pos0(magnitude_Mach_L - 1.0) *
                0.5 * (1.0 + sign(Mach_L)) +
            neg(magnitude_Mach_L - 1.0) *
                (0.25 * sqr(Mach_L + 1.0) * (2.0 - Mach_L) +
                 alpha * Mach_L * sqr(sqr(Mach_L) - 1.0));
        tmp<surfaceScalarField> p_minus_R =
            pos0(magnitude_Mach_R - 1.0) *
                0.5 * (1.0 - sign(Mach_R)) +
            neg(magnitude_Mach_R - 1.0) *
                (0.25 * sqr(Mach_R - 1.0) * (2.0 + Mach_R) -
                 alpha * Mach_R * sqr(sqr(Mach_R) - 1.0));

        // second half of Eq.19
        tmp<surfaceScalarField> dp =
            (p_plus_L.ref() - p_minus_R.ref()) * 0.5 * (p_L.ref() - p_R.ref()) +
            f_0 * ((p_plus_L.ref() + p_minus_R.ref() - 1) * p_mean);

        // Eq. 19
        const surfaceScalarField ps("ps", p_mean + dp.ref());
        dp.clear();

        // Eq. 26 velocity diffusion
        surfaceVectorField pu(
            "pu",
            -g * gamma_mean * p_mean / c_face * p_plus_L.ref() * p_minus_R.ref() * (u_R.ref() - u_L.ref()));
        p_plus_L.clear();
        p_minus_R.clear();
        // Eq. 27
        surfaceVectorField p12("p12", ps * face_normal + pu);
        #if OPENFOAM >= 1712
        p12.setOriented(false);
        #endif

        // MASS FLUX Mach_L_+ + Mach_R_- + Mp
        // Eq. 6
        const surfaceScalarField Mach_plus_L(
            "Mach_plus_L",
            pos0(magnitude_Mach_L - 1.0) *
                    0.5 * (Mach_L + magnitude_Mach_L) +
                neg(magnitude_Mach_L - 1.0) *
                    (0.25 * sqr(Mach_L + 1.0) + 0.125 * sqr(sqr(Mach_L) - 1.0)));
        const surfaceScalarField Mach_minus_R(
            "Mach_minus_R",
            pos0(magnitude_Mach_R - 1.0) *
                    0.5 * (Mach_R - magnitude_Mach_R) +
                neg(magnitude_Mach_R - 1.0) *
                    (-0.25 * sqr(Mach_R - 1.0) - 0.125 * sqr(sqr(Mach_R) - 1.0)));

        // Eq. 14
        tmp<surfaceScalarField> Mp = -0.5 * (1.0 - f) * (p_R.ref() - p_L.ref()) / rho_mean / sqr(c_face) * (1.0 - g);

        // Eq. 10
        const surfaceScalarField Mach_1_2(
            "Mach_1_2",
            Mach_plus_L + Mach_minus_R + Mp.ref());
        Mp.clear();

        // Eq. 5
        surfaceScalarField mdot(
            "mdot",
            c_face * Mach_1_2 * (rho_L * pos0(Mach_1_2) + rho_R * neg(Mach_1_2)));
        #if OPENFOAM >= 1712
        mdot.setOriented(false);
        #endif

        // Face velocity for sigmaDotU (turbulence term)
        Up = surface_area * (u_L.ref() * pos0(Mach_1_2) + u_R.ref() * neg(Mach_1_2));
        #if OPENFOAM >= 1712
        Up.setOriented(true);
        #endif

        phi = surface_area * mdot;
        #if OPENFOAM >= 1712
        phi.setOriented(true);
        #endif

        phiUp =
            surface_area * (0.5 * (mdot + mag(mdot)) * u_L.ref() +
                            0.5 * (mdot - mag(mdot)) * u_R.ref() + p12);
        #if OPENFOAM >= 1712
        phiUp.setOriented(true);
        #endif

        phiEp =
            surface_area * (0.5 * (mdot + mag(mdot)) * h_L.ref() +
                            0.5 * (mdot - mag(mdot)) * h_R.ref());
        #if OPENFOAM >= 1712
        phiEp.setOriented(true);
        #endif

        if (debug)
        {
            Info << max(Up) << endl;
            Info << min(mag(Up)) << endl;
            Info << max(phi) << endl;
            Info << min(mag(phi)) << endl;
            Info << max(phiUp) << endl;
            Info << min(mag(phiUp)) << endl;
            Info << max(phiEp) << endl;
            Info << min(mag(phiEp)) << endl;
            Up.write();
            phiUp.write();
            phiEp.write();
            phi.write();
        }

        // Info << endl;
        Info << "Done updating fluxes!" << endl;
        // Info << endl;
    }

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
    Foam::surfaceScalarField ausmPlusMFluxScheme::calc_h(
        surfaceScalarField p_R,
        surfaceScalarField p_L)
    {
        // Eq. 44: h_k for each surface
        const surfaceScalarField h_k("h_k",
                                     min(p_L / p_R, p_R / p_L));

        // volume field that stores minimum h_k for the cells surfaces
        volScalarField h_min(
            IOobject(
                "h_min",
                #if FOUNDATION >= 12
                mesh_.time().name(),
                #else
                mesh_.time().timeName(),
                #endif
                mesh_),
            mesh_,
            dimensionedScalar("h_min", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0));

        forAll(h_min, celli)
        {
            const labelList &cellFaces = mesh_.cells()[celli]; // get list of faces of current cell
            // (https://www.cfd-online.com/Forums/openfoam-solving/120990-small-great-rootvsmall-what.html)
            // iterate over all faces of celli
            forAll(cellFaces, facei)
            {
                // find minimum h_k of the cell; ONLY FOR INTERNAL BOUNDARIES
                if (cellFaces[facei] < h_k.size())
                {
                    h_min[celli] = min(h_min[celli], h_k[cellFaces[facei]]);
                }
            }
        };

        // Eq. 43: interpolate on cell faces
        const surfaceScalarField h43(
            "h43", min(
                       fvc::surfaceReconstruct(h_min, own_, RECONSTRUCT_SCALAR_),
                       fvc::surfaceReconstruct(h_min, nei_, RECONSTRUCT_SCALAR_)));

        return h43;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
