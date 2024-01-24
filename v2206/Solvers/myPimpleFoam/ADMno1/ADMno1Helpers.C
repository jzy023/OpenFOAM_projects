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

//- Inhibition calculations=

volScalarField Foam::ADMno1::calcInhibition
(
    volScalarField Y, 
    scalar denom
)
{
    return 1 / (1 + (Y / denom));
}

volScalarField Foam::ADMno1::calcInhibitionHP
(
    volScalarField Shp,
    scalar UL,
    scalar LL,
    label n
)
{
    scalar Kph = pow(10, -0.5 * (UL + LL));
    return pow(Kph, n) / (pow(Shp, n) + pow(Kph, n));
}

//- Reaction rate calculations

volScalarField Foam::ADMno1::calcRho
(
    scalar k, 
    volScalarField X
)
{
    return k * X;
}

volScalarField Foam::ADMno1::calcRho
(
    scalar k, 
    volScalarField S,
    scalar K,
    volScalarField X,
    volScalarField I
)
{
    return k * (S / (K + S)) * X * I;
}

volScalarField Foam::ADMno1::calcRho
(
    scalar k, 
    volScalarField S1,
    scalar K,
    volScalarField X,
    volScalarField S2,
    volScalarField I
)
{
    return k * (S1 / (K + S1)) * X * (1 / (1 + (S2 / S1))) * I;
}




// ************************************************************************* //