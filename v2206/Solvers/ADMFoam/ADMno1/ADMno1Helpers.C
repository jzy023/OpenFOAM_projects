/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Inhibition calculations

volScalarField::Internal Foam::ADMno1::calcInhibition
(
    volScalarField Y, 
    dimensionedScalar denum
)
{
    return 1.0/ (1.0 + (Y.internalField() / denum));
}

volScalarField::Internal Foam::ADMno1::calcInhibitionHP
(
    volScalarField::Internal Shp,
    const dimensionedScalar UL,
    const dimensionedScalar LL,
    const dimensionedScalar n
)
{
    dimensionedScalar Kph = pow(10, -0.5 * (UL + LL));
    Kph.dimensions().reset(Shp.dimensions());

    return pow(Kph, n) / (pow(Shp, n) + pow(Kph, n));
}


//- Kinetic rate calculations

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k,
    volScalarField X
)
{
    return k * X;
}

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k, 
    volScalarField S,
    const dimensionedScalar K,
    volScalarField X,
    volScalarField::Internal I
)
{
    return k * (S / (K + S.internalField())) * X.internalField() * I;
}

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k, 
    volScalarField S1,
    const dimensionedScalar K,
    volScalarField X,
    volScalarField S2,
    volScalarField::Internal I
)
{
    return k * (S1.internalField() / (K + S1.internalField())) * X.internalField() * 
           (1.0 / (1.0 + (S2.internalField() / S1.internalField()))) * I;
}


//- Acid-base calculations

volScalarField::Internal Foam::ADMno1::fSion
(
    const dimensionedScalar Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return Kax * Sx / (Kax + Shp);
}

volScalarField::Internal Foam::ADMno1::fSion
(
    const volScalarField Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return Kax * Sx / (Kax + Shp);
}

volScalarField::Internal Foam::ADMno1::dfSion
(
    const dimensionedScalar Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return - Kax * Sx / ((Kax + Shp) * (Kax + Shp));
}

volScalarField::Internal Foam::ADMno1::dfSion
(
    const volScalarField Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return - Kax.internalField() * Sx / ((Kax + Shp) * (Kax + Shp));
}


//- Components source term calculations

volScalarField::Internal Foam::ADMno1::concPerComponent
(
    label j,
    const admPara para,
    PtrList<volScalarField::Internal> KRPtrs
)
{
    volScalarField::Internal dY
    (
        IOobject
        (
            "dY",
            KRPtrs[0].mesh().time().timeName(),
            KRPtrs[0].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        KRPtrs[0].mesh(),
        dimensionedScalar
        (
           "dY_Default", 
            KRPtrs[0].dimensions(),
            Zero
        )
    );

    // TODO: room for optimization (dont loop through all 19 elements since a lot of them are 0s)
    for (int i = 0; i < 19; i++) {
        dY += para.STOI[i][j] * KRPtrs[i] * para.DTOS(); //check if it works
    }
    return dY;
}

// ************************************************************************* //