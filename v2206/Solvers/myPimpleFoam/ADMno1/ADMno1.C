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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
    Sc_(ADMno1Dict.lookupOrDefault("Sc", 0.2))
{

    Info<< "\nSelecting ADM no1 operation mode " << opMode_ << endl;
    
    modeCheckErr();

    admParameters_.setOpMode(namesOpMode.find(opMode_));
    

    //- Recreate "createADMFields.H"

    Info<< "Reading ADM no1 initial concentrations for soluables" << endl;

    label iNames = 0;
    label nSpecies = namesSoluable.size() + namesGaseous.size();
                // TODO: add other species
                //    + namesParticulate.size() + namesMedians.size();

    YPtrs_.resize(nSpecies);

    forAll(namesSoluable, i)
    {
        // Info<< namesSoluable[i] << endl;
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

    //-  Read gaseuoses

    Info<< "Reading ADM no1 initial concentrations for gaseuoses" << endl;

    forAll(namesGaseous, i)
    {
        YPtrs_.set
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

    // iNames += sizeof(namesGaseous);

    // //-  Read particulates

    // Info<< "Reading ADMno1 initial concentrations for particulates" << endl;

    // forAll(namesParticulate, i)
    // {
    //     YPtrs_.set
    //     (
    //         i + iNames,
    //         new volScalarField
    //         (
    //             IOobject
    //             (
    //                 namesParticulate[i],
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::MUST_READ,
    //                 IOobject::AUTO_WRITE
    //             ),
    //             mesh
    //         )
    //     );
    // }

    // iNames += sizeof(namesParticulate);

    // //-  Read electrolytes

    // Info<< "Reading ADMno1 initial concentrations for electrolytes" << endl;

    // forAll(namesElectrolyte, i)
    // {
    //     YPtrs_.set
    //     (
    //         i + iNames,
    //         new volScalarField
    //         (
    //             IOobject
    //             (
    //                 namesElectrolyte[i],
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::MUST_READ,
    //                 IOobject::AUTO_WRITE
    //             ),
    //             mesh
    //         )
    //     );
    // }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Set up temperal concentration fields

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Info<< "Reading ADMno1 initial operating temperaturesbased on runMode" << endl;

    // TODO: if statement from runMode
    // scalar TopDefault_ = admParameters_.getTbase();

    dimensionedScalar TopDefault("TopDefault", dimensionSet(0,0,0,1,0,0,0), admParameters_.getTbase()); 
    Info<< "Reading field operational temperature" << endl;
    volScalarField Top
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
        TopDefault
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
            IOobject::MUST_READ_IF_MODIFIED,
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