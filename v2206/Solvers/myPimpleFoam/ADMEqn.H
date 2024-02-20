/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    ADM::ADMEqn

Description
    ADM computations.

\*---------------------------------------------------------------------------*/

#ifndef ADMEqn_H
#define ADMEqn_H

// #include "typeDefs.H"
// #include "ADMno1Parameter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info << "ADM convection-diffusion calculation" << endl;

PtrList<volScalarField>& YPtrs = reaction->Y();

volScalarField Deff = turbulence->nu() / reaction->Sc();

volScalarField rho
(
    IOobject
    (
        "rho",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    reaction->paraPtr()->rho()
);
// Deff.dimensions().reset(Deff.dimensions()*dimDensity);
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(YPtrs, i)
{
    fields.add(YPtrs[i]);
}
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi)")
    )
);


{
    // reaction->correct(Top);
    // Qdot = reaction->Qdot();
    // volScalarField Yt(0.0*Y[0]);

    // volScalarField& Yi = YPtrs[7];
    // Info << Yi.dimensions() << endl;
    // Info << fvm::ddt(rho, Yi)->dimensions() << endl;
    // Info << mvConvection->fvmDiv(phi, Yi)->dimensions() << endl;
    // Info << reaction->R(7)->dimensions() << endl;

    forAll(YPtrs, i)
    {
        volScalarField& Yi = YPtrs[i];

        // fvScalarMatrix YiEqn
        // (
        //     fvm::ddt(rho, Yi)
        //   + mvConvection->fvmDiv(phi, Yi)
        //   - fvm::laplacian(Deff, Yi)
        //  ==
        //     fvOptions(rho, Yi)
        //   + reaction->R(i)
        // );

        fvScalarMatrix YiEqn
        (
            fvm::ddt(Yi)
          + fvm::div(phi, Yi)
          - fvm::laplacian(Deff, Yi)
         ==
            fvOptions(Yi)
          + reaction->R(i)
        );

        YiEqn.relax();

        fvOptions.constrain(YiEqn);

        YiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Yi);
    }

    // Gas Phase ============================================================ //
    // TODO: add rho and make adjustment for compressibility

    // <<< TEST ---------------------------------------------------------------- //
    // volScalarField rho
    // (
    //     IOobject
    //     (
    //         "rho",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar
    //     (
    //         "rho", 
    //         dimensionSet(1,-3,0,0,0,0,0), 
    //         1000
    //     )
    // );
    // volScalarField DMuEff = Deff*rho;
    // <<< TEST ---------------------------------------------------------------- //

    // forAll(GPtrs, i)
    // {
    //     volScalarField& Gi = GPtrs[i];

    //     fvScalarMatrix GiEqn
    //     (
    //         fvm::ddt(Gi)
    //       + fvm::div(phi, Gi)
    //       - fvm::laplacian(Deff, Gi)
    //      ==
    //         fvOptions(Gi)
    //     //   + reaction->RG(i, Gi)
    //     );

    //     GiEqn.relax();

    //     fvOptions.constrain(GiEqn);

    //     GiEqn.solve(mesh.solver("Gi"));

    //     fvOptions.correct(Gi);
    // }

    // Gas Phase ============================================================ //

    // reaction->clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //