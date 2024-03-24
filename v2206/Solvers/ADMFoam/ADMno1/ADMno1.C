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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ADMno1::ADMno1
(
    const fvMesh& mesh,
    const IOdictionary& ADMno1Dict
)
:
    IOdictionary(ADMno1Dict),
    para_(ADMno1Dict.get<word>("mode")),
    Sc_(ADMno1Dict.lookupOrDefault("Sc", 0.2)),
    R_(ADMno1Dict.lookupOrDefault("R", 0.083145)),
    KP_(ADMno1Dict.lookupOrDefault("Kpip", 5e4)),
    Vfrac_(ADMno1Dict.lookupOrDefault("Vfrac", 0.088235)),
    Pext_
    (
        "Pext", 
        dimPressure,
        ADMno1Dict.lookupOrDefault("Pext", 1.013)
    ),
    fac_
    (
        IOobject
        (
            "fac",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "facDefault", 
            dimless, 
            ADMno1Dict.lookupOrDefault("fac", 1)
        )
    ),
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
            dimless, 
            ADMno1Dict.lookupOrDefault("pH", 7.26)
        )
    ),
    ShP_
    (
        IOobject
        (
            "ShP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
           "ShP_Default", 
            dimMass,
            para_.Pini()
        )
    ),
    Scat_
    (
        "Scat",
        dimMass/dimVolume, //TODO
        ADMno1Dict.lookupOrDefault("Scat", 0.00)
    ),
    San_
    (
        "San",
        dimMass/dimVolume, //TODO
        ADMno1Dict.lookupOrDefault("San", 0.0052)
    )
{

    Info<< "\nSelecting ADM no1 operation mode " << ADMno1Dict.get<word>("mode") << endl;

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

    //- Initializing derivatives

    dYPtrs_.resize(nSpecies);

    for(label i = 0; i < nSpecies; i++)
    {
        dYPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "d" + YPtrs_[i].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    YPtrs_[0].dimensions()/dimTime, 
                    Zero
                )
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
                    YPtrs_[0].dimensions(),
                    para_.Gini(i)
                )
            )
        );
    }
    
    //- Initializing derivatives

    dGPtrs_.resize(namesGaseous.size());

    for(label i = 0; i < namesGaseous.size(); i++)
    {
        dGPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "d" + GPtrs_[i].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    GPtrs_[0].dimensions()/dimTime, 
                    Zero
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
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    namesMedians[i] + "Default", 
                    YPtrs_[0].dimensions(),
                    para_.Mini(i)
                )
            )
        );
    }

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
                    "Inh" + Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE// TODO: choose if writing
                ),
                mesh,
                dimensionedScalar
                (
                    dimless, 
                    Zero
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Kinetic rate initialization
    
    KRPtrs_.resize(19);

    for (int i = 0; i < 19; i++)
    {
        KRPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRs" + Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    dYPtrs_[0].dimensions(), 
                    Zero
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Gas trasfer rate initialization
    
    GRPtrs_.resize(3);

    for (int i = 0; i < 3; i++)
    {
        GRPtrs_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "GRs" + Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    dGPtrs_[0].dimensions(), 
                    Zero
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // reset dimensions 
    para_.setParaDim(YPtrs_[0].dimensions());
    ShP_.dimensions().reset(YPtrs_[0].dimensions());
    Scat_.dimensions().reset(YPtrs_[0].dimensions());
    San_.dimensions().reset(YPtrs_[0].dimensions());

    // 
    nIaa_ = 3.0 / (para_.pHL().ULaa - para_.pHL().LLaa);  // aa
    nIac_ = 3.0 / (para_.pHL().ULac - para_.pHL().LLac);  // ac
    nIh2_ = 3.0 / (para_.pHL().ULh2 - para_.pHL().LLh2);  // h2

    // >>> TEST
    // dYPtrs_[7].field() = 0.01*2.5055e-7;

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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ADMno1::clear()
{
    forAll(dYPtrs_, i)
    {
        dYPtrs_[i] *= 0.0;
    }

    forAll(dGPtrs_, i)
    {
        dGPtrs_[i] *= 0.0;
    }
}

void Foam::ADMno1::correct(volScalarField& T)
{
    //- calculate thermal factor
    // thermalFac(T);

    //- calculate gas phase transfer rates
    gasPhaseRate(T);

    //- calculate gas exit rates
    gasSourceRate(T);

    //- calculate raction rates
    kineticRate();

    //- Acid-base calculations
    calcShp();

    //- calculate dY with STOI
    sulSourceRate();

}


tmp<fvScalarMatrix> Foam::ADMno1::R
(
    label i
) const
{
    DimensionedField<scalar, volMesh> dY = dYPtrs_[i];

        tmp<fvScalarMatrix> tSu
        (
            new fvScalarMatrix
            (
                YPtrs_[i],
                dY.dimensions()*dimVolume
                // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
            )
        );

    fvScalarMatrix& Su = tSu.ref();
    
    // https://www.openfoam.com/documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
    Su += dY; 

    return tSu;
}; 


tmp<fvScalarMatrix> Foam::ADMno1::RG
(
    label i
) const
{
    DimensionedField<scalar, volMesh> dG = dGPtrs_[i];

        tmp<fvScalarMatrix> tSu
        (
            new fvScalarMatrix
            (
                GPtrs_[i],
                dG.dimensions()*dimVolume
                // dimMass/dimTime // <- for compressible flow and uses fvm::ddt(rho, Yi)
            )
        );

    fvScalarMatrix& Su = tSu.ref();
    
    // https://www.openfoam.com/documentation/guides/latest/api/fvMatrix_8C_source.html#l01708
    Su += dG; 

    return tSu;
}; 

// ************************************************************************* //