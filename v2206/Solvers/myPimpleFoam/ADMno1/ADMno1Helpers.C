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
    volScalarField* YPtr, 
    scalar denom
)
{
    return 1 / (1 + (*YPtr / denom));
}

volScalarField Foam::ADMno1::calcInhibitionHP
(
    volScalarField* ShpPtr,
    scalar UL,
    scalar LL,
    label n
)
{
    scalar Kph = pow(10, -0.5 * (UL + LL));
    return pow(Kph, n) / (pow(*ShpPtr, n) + pow(Kph, n));
}

//- Reaction rate calculations

volScalarField Foam::ADMno1::calcRho
(
    scalar k, 
    volScalarField* XPtr
)
{
    return k * *XPtr;
}

volScalarField Foam::ADMno1::calcRho
(
    scalar k, 
    volScalarField* SPtr,
    scalar K,
    volScalarField* XPtr,
    scalar I
)
{
    return k * (*SPtr / (K + *SPtr)) * *XPtr * I;
}

volScalarField Foam::ADMno1::calcRho
(
    scalar k, 
    volScalarField* SPtr1,
    scalar K,
    volScalarField* XPtr,
    volScalarField* SPtr2,
    scalar I
)
{
    return k * (*SPtr1 / (K + *SPtr1)) * *XPtr * (1 / (1 + (*SPtr2 / *SPtr1))) * I;
}




// ************************************************************************* //