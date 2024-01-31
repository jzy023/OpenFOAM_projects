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
// TODO: make it an input from ./constant/ADMpara file
volScalarField Deff = turbulence->nu() / reaction->Sc();

// Create ADM class instance (ADMno1 constructor)
// reaction->R(Yi) becomses the source term
// PtrList<volScalarField>& YPtrs = reaction.Y();
PtrList<volScalarField>& YPtrs = reaction->Y();

{
    reaction->correct(Top);
    // Qdot = reaction->Qdot();
    // volScalarField Yt(0.0*Y[0]);

    forAll(YPtrs, i)
    {
        volScalarField& Yi = YPtrs[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(Yi)
          + fvm::div(phi, Yi)
          - fvm::laplacian(Deff, Yi) 
        //  TODO: add turbulent diffusion
        // - fvm::laplacian(turbulence->muEff(), Yi)
         ==
            fvOptions(Yi)
        //  TODO: add local biochemical reaction as source term
          + reaction->R(i, Yi)
        // + fvOptions(rho, Yi)
        );

        YiEqn.relax();

        fvOptions.constrain(YiEqn);

        YiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Yi);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //