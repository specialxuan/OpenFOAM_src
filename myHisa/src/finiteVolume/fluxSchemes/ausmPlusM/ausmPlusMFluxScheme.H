/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2016 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2016 Oliver Oxtoby - CSIR, South Africa
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
        UniBw Munich, Germany
    Maximilian Mudra
        UniBw Munich, Germany

\*---------------------------------------------------------------------------*/

#ifndef ausmPlusMFluxScheme_H
#define ausmPlusMFluxScheme_H

#include "fluxScheme.H"
#include "fvcSurfaceReconstruct.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                  Class ausmPlusMFluxScheme Declaration
    \*---------------------------------------------------------------------------*/

    class ausmPlusMFluxScheme
        : public fluxScheme
    {

    protected:
        // Protected data
        const fvMesh &mesh_;
        const fluidThermo &thermo_;
        const volScalarField &rho_;
        const volVectorField &U_;
        const volVectorField &rhoU_;
        const volScalarField &rhoE_;
        const dictionary &dict_;

        //- Directionality for field interpolation
        const surfaceScalarField own_, nei_;

        // const dictionary ausm_dict_;

        // fvSchemes must contain: reconstruct(T) wVanLeer; reconstruct(U) wVanLeerV;
        const string RECONSTRUCT_SCALAR_ = "reconstruct(T)";
        const string RECONSTRUCT_VECTOR_ = "reconstruct(U)";

        const scalar PI_ = constant::mathematical::pi;

    private:
        // Parameters of AUSM+M flux scheme. Can be set in fvSchemes dict.
        scalar alpha0_, sqrMachInf_;

        // pressure based shock sensing function to activate diffusive velocity term
        surfaceScalarField calc_h(
            surfaceScalarField p_R,
            surfaceScalarField p_L);

    public:
        //- Runtime type information
        TypeName("AUSMPlusM");

        // Constructors

        //- Construct from components
        ausmPlusMFluxScheme(
            const dictionary &dict,
            const fluidThermo &thermo,
            const volScalarField &rho,
            const volVectorField &U,
            const volVectorField &rhoU,
            const volScalarField &rhoE);

        // Destructor

        ~ausmPlusMFluxScheme();

        // Member Functions

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux(
            surfaceScalarField &phi,
            surfaceVectorField &phiUp,
            surfaceScalarField &phiEp,
            surfaceVectorField &Up) override;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
