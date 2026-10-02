/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2020-2021 Oliver Oxtoby

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

\*---------------------------------------------------------------------------*/

    // Optional failsafe to avoid -ve rho. Avoid dividing by very small quantity
    #if FOUNDATION >= 10
    const dictionary& dict = mesh.solution().subDict("pseudoTime");
    #else
    const dictionary& dict = mesh.solutionDict().subDict("pseudoTime");
    #endif
    dimensionedScalar rhoMin
    (
        "rhoMin",
        dimDensity,
        dict.lookupOrDefault<scalar>("rhoMin", -GREAT)
    );
    bound(rho, rhoMin);

    // Compute internal fields

    #if FOUNDATION >= 12
    U.internalFieldRef() = rhoU.internalField()/rho.internalField();
    #else
    U.ref() = rhoU.internalField()/rho.internalField();
    #endif

    #if FOUNDATION >= 12
    e.internalFieldRef() = rhoE.internalField()/rho.internalField()
              - 0.5*magSqr(U.internalField()) + eZero_->internalField();
    #else
    e.ref() = rhoE.internalField()/rho.internalField()
              - 0.5*magSqr(U.internalField()) + eZero_->internalField();
    #endif

    // Temperature bound
    dimensionedScalar TMin
    (
        "TMin",
        dimTemperature,
        dict.lookupOrDefault<scalar>("TMin", ROOTSMALL)
    );
    dimensionedScalar TMax
    (
        "TMax",
        dimTemperature,
        dict.lookupOrDefault<scalar>("TMax", GREAT)
    );

    // Bound energy
    volScalarField TBound
    (
        IOobject
        (
            "TBound", 
            #if FOUNDATION >= 12
            runTime.name(),
            #else
            runTime.timeName(),
            #endif
            mesh
        ),
        mesh,
        TMin
    );
    volScalarField eBound(thermo.he(thermo.p(), TBound));

    if (max(neg(e-eBound)).value() > 0.5)
    {
        e = max(e, eBound);
        bounded_ = true;

         Info<< "Bounding " << e.name()
             << " to TMin: " << TMin.value()
             << endl;
    }

    // Only do max bound if it was specified
    if (TMax.value() < GREAT)
    {
        TBound == TMax;
        eBound == thermo.he(thermo.p(), TBound);
        if (max(pos(e-eBound)).value() > 0.5)
        {
            e = min(e, eBound);
            Info<< "Bounding " << e.name()
                << " to TMax: " << TMax.value()
                << endl;
        }
    }

    // Calc T and psi from e
    thermo.correct();

    // Update p from rho
    // Use Newton's method to converge rho = rho(p)
    volScalarField::Internal dRho
    (
        rho.internalField()-scalarVarsPrevIter_[0].internalField()
    );
    #if FOUNDATION >= 12
    p.internalFieldRef() += dRho/psi.internalField();
    #else
    p.ref() += dRho/psi.internalField();
    #endif
    label i = 0;
    dimensionedScalar pSmall("pSmall", p.dimensions(), SMALL);
    while (true)
    {
        // Calculate rho & psi from new p
        thermo.correct();

        //Update p
        volScalarField::Internal deltap
        (
           -(
                thermo.rho()().internalField()-rho.internalField()
            )/psi.internalField()
        );
        #if FOUNDATION >= 12
        p.internalFieldRef() += deltap;
        #else
        p.ref() += deltap;
        #endif

        #if FOUNDATION >= 12
        if (gMax(mag(deltap/stabilise(p.internalField(), pSmall))()) < 1e-5)
        #else
        if (gMax(mag(deltap/stabilise(p.internalField(), pSmall))()) < 1e-5)
        #endif
        {
            break;
        }

        if (i > 100)
        {
            WarningInFunction
                << "Newton's method failed to converge in density update."
                << nl;
            break;
        }
    }

    // Recalc rhoU and rhoE in case rho or e were bounded
    #if FOUNDATION >= 12
    rhoU.internalFieldRef() = rho.internalField()*U.internalField();
    rhoE.internalFieldRef() =
        rho.internalField()*
        (
            e.internalField() - eZero_() + 0.5*magSqr(U.internalField())
        );
    #else
    rhoU.ref() = rho.internalField()*U.internalField();
    rhoE.ref() =
        rho.internalField()*
        (
            e.internalField() - eZero_() + 0.5*magSqr(U.internalField())
        );
    #endif

    // Correct boundary fields
    p.correctBoundaryConditions();
    U.correctBoundaryConditions();
    T.correctBoundaryConditions();

    T.boundaryFieldRef() == max(T.boundaryField(), TMin.value());
    if (TMax.value() < GREAT)
    {
        T.boundaryFieldRef() == min(T.boundaryField(), TMax.value());
    }

    e.boundaryFieldRef() == thermo.he(p, T)->boundaryField();

    thermo.correct();                       // NOTE: Correct psi boundaryField

    tmp<volScalarField> trho = thermo.rho();
    forAll(rho.boundaryField(), patchi)
    {
        // Substitute wall pressure with zero grad to avoid unphysical values
        if
        (
            isA<characteristicWallPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            rho.boundaryFieldRef()[patchi] =
                trho().boundaryField()[patchi].patchInternalField();
        }
        else
        {
            rho.boundaryFieldRef()[patchi] = trho().boundaryField()[patchi];
        }
    }
    rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();
    rhoE.boundaryFieldRef() =
        rho.boundaryField()*
        (
            e.boundaryField()
          - eZero_->boundaryField()
          + 0.5*magSqr(U.boundaryField())
        );
