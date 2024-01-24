/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
 
\*---------------------------------------------------------------------------*/

#include "ADMno1.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
 
namespace Foam
{
    defineTypeNameAndDebug(ADMno1, 0);
    // defineRunTimeSelectionTable(ADMno1, fvMesh);
}
 
const Foam::word Foam::ADMno1::propertiesName("admno1Properties");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ADMno1::modeCheckErr()
{
    if (!namesOpMode.found(opMode_))
    {
        std::cerr << nl << "--> FOAM FATAL IO ERROR:" << nl
                  << "Unknown operation mode type: " << opMode_
                  << "\n\nValid operation mode types:\n";
        std::cerr << namesOpMode.size() << "(\n";
        forAll(namesOpMode, i)
        {
            std::cerr << namesOpMode[i] << nl;
        }
        std::cerr << ")\n";
        std::exit(1);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ADMno1::ADMno1
// ADMno1::ADMno1
(
    const fvMesh& mesh,
    const IOdictionary& ADMno1Dict
)
:
    IOdictionary(ADMno1Dict),
    opMode_(ADMno1Dict.get<word>("mode")),
    Sc_(ADMno1Dict.lookupOrDefault("Sc", 0.2)),
    R_(ADMno1Dict.lookupOrDefault("R", 0.083145)),
    pH_
    (
        IOobject
        (
            "pH",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "pHDefault", 
            dimensionSet(0,0,0,0,0,0,0), 
            ADMno1Dict.lookupOrDefault("pH", 7.26)
        )
    ),
    Scat_(ADMno1Dict.lookupOrDefault("Scat", 0.00)),
    San_(ADMno1Dict.lookupOrDefault("San", 0.0052))
{

    Info<< "\nSelecting ADM no1 operation mode " << opMode_ << endl;
    
    modeCheckErr();

    para_.setOpMode(namesOpMode.find(opMode_));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Main substances concentration initialization

    Info<< "Reading ADM no1 initial concentrations for soluables" << endl;

    label iNames = 0;
    label nSpecies = namesSoluable.size() + namesParticulate.size();

    YPtrs_.resize(nSpecies);

    forAll(namesSoluable, i)
    {
        YPtrs_.set
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

    iNames += namesSoluable.size();

    //-  Read particulates

    Info<< "Reading ADMno1 initial concentrations for particulates" << endl;

    forAll(namesParticulate, i)
    {
        YPtrs_.set
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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //-  Gaseuoses initialization

    Info<< "Reading ADM no1 initial concentrations for gaseuoses" << endl;

    GPtrs_.resize(namesGaseous.size());

    forAll(namesGaseous, i)
    {
        GPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    namesGaseous[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    namesGaseous[i] + "Default", 
                    dimensionSet(0,0,0,0,0,0,0),
                    para_.getGini(i)
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //-  Electrolytes (calculated) initialization

    Info<< "Initializing concentrations for electrolyte" << endl;

    EPtrs_.resize(namesElectrolyte.size());

    forAll(namesElectrolyte, i)
    {
        EPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    namesElectrolyte[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    namesElectrolyte[i] + "Default", 
                    dimensionSet(0,0,0,0,0,0,0), 
                    para_.getEini(i)
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //-  Medians initialization

    Info<< "Initializing concentrations for medians" << endl;

    MPtrs_.resize(namesMedians.size());

    forAll(namesMedians, i)
    {
        MPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    namesMedians[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    namesMedians[i] + "Default", 
                    dimensionSet(0,0,0,0,0,0,0), 
                    para_.getMini(i)
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Set up temperal concentration fields

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Inhibition coeffs initialization

    IPtrs_.resize(8);

    for (int i = 0; i < 8; i++)
    {
        IPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Inh",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE// TODO: choose if writing
                ),
                mesh,
                dimensionedScalar
                (
                    "Inh", 
                    dimensionSet(0,0,0,0,0,0,0), 
                    1e-20
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Reaction rate initialization
    
    RRPtrs_.resize(19);

    for (int i = 0; i < 19; i++)
    {
        RRPtrs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "RRs",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "RRs", 
                    dimensionSet(0,0,0,0,0,0,0), 
                    1e-20
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading field operational temperature" << endl;

    volScalarField Top_
    (
        IOobject
        (
            "Top",
            mesh.time().timeName(),  // or runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "TopDefault", 
            dimensionSet(0,0,0,1,0,0,0), 
            para_.getTbase()
        )
    );
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
 
Foam::autoPtr<Foam::ADMno1> Foam::ADMno1::New
(
    const fvMesh& mesh
)
{
    IOdictionary ADMno1Dict
    (
        IOobject
        (
            propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    // TODO: do it properly!!! with virtual destructors and constructor hash tables 
    // new keywaord is not gonna last!
    ADMno1* reactionPtr = new ADMno1(mesh, ADMno1Dict);
    return autoPtr<ADMno1>(reactionPtr);
}


// ************************************************************************* //