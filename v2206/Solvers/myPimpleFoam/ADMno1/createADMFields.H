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
    ADM::typeDefs

Description
    Data structure definition for ADMno1 parameters .

\*---------------------------------------------------------------------------*/

#ifndef createADMFields_H
#define createADMFields_H

// #include "typeDefs.H"
#include "ADMno1Parameter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

wordList namesSoluable = {"Ssu", "Saa", "Sfa", "Sva", "Sbu", "Spro",
                          "Sac", "Sh2", "Sch4", "SIC", "SIN", "SI", "Sco2"};

wordList namesGaseous = {"Gh2", "Gch4", "Gco2"};

wordList namesParticulate = {"Xc", "Xch", "Xpr", "Xli", "Xsu", "Xaa",
                              "Xfa", "Xc4", "Xpro", "Xac", "Xh2", "XI"};

wordList namesElectrolyte = {"Scat", "San", "ShP", "SohN"};

// making ptr list for all the components
// PtrList<volScalarField>& cmptPtrs(28); // 27 + Sco2 
PtrList<volScalarField> YPtrs(15);
label iNames = 0;

//-  Read soluables

Info<< "Reading ADMno1 initial concentrations for soluables" << endl;

forAll(namesSoluable, i)
{
    YPtrs.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                namesSoluable[i], // IOobject::groupName(namesSoluable[i]),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

iNames += sizeof(namesSoluable);

//-  Read gaseuoses

Info<< "Reading ADMno1 initial concentrations for gaseuoses" << endl;

forAll(namesGaseous, i)
{
    YPtrs.set
    (
        i + iNames,
        new volScalarField
        (
            IOobject
            (
                namesGaseous[i],
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

iNames += sizeof(namesGaseous);

//-  Read particulates

Info<< "Reading ADMno1 initial concentrations for particulates" << endl;

forAll(namesParticulate, i)
{
    YPtrs.set
    (
        i + iNames,
        new volScalarField
        (
            IOobject
            (
                namesParticulate[i],
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

iNames += sizeof(namesParticulate);

//-  Read electrolytes

Info<< "Reading ADMno1 initial concentrations for electrolytes" << endl;

forAll(namesElectrolyte, i)
{
    YPtrs.set
    (
        i + iNames,
        new volScalarField
        (
            IOobject
            (
                namesElectrolyte[i],
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //