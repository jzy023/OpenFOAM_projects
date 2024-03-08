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

volScalarField Foam::ADMno1::calcInhibition
(
    volScalarField Y, 
    dimensionedScalar denum
)
{
    return 1 / (1 + (Y / denum));
}

volScalarField Foam::ADMno1::calcInhibitionHP
(
    volScalarField Shp,
    dimensionedScalar UL,
    dimensionedScalar LL,
    dimensionedScalar n
)
{
    dimensionedScalar Kph = pow(10, -0.5 * (UL + LL));
    Kph.dimensions().reset(Shp.dimensions());

    return pow(Kph, n) / (pow(Shp, n) + pow(Kph, n));
}

//- Kinetic rate calculations

volScalarField Foam::ADMno1::calcRho
(
    dimensionedScalar k,
    volScalarField X
)
{
    return k * X;
}

volScalarField Foam::ADMno1::calcRho
(
    dimensionedScalar k, 
    volScalarField S,
    dimensionedScalar K,
    volScalarField X,
    volScalarField I
)
{
    return k * (S / (K + S)) * X * I;
}

volScalarField Foam::ADMno1::calcRho
(
    dimensionedScalar k, 
    volScalarField S1,
    dimensionedScalar K,
    volScalarField X,
    volScalarField S2,
    volScalarField I
)
{
    return k * (S1 / (K + S1)) * X * (1 / (1 + (S2 / S1))) * I;
}

//- Components source term calculations

// volScalarField Foam::ADMno1::concPerComponent
// (
//     const admPara *paraPtr, 
//     int *jPtr
// )
// {
//     data_type concComponent = 0;
//     for (int i = 0; i < 19; i++) {
//         data_type temp = (paraPtr->STOI[i][*jPtr]) * (rho[i]); //check if it works
//         concComponent += temp;
//     }
//     *jPtr += 1;
//     return concComponent;
// }

//- Acid-base calculations

volScalarField Foam::ADMno1::fSion
(
    dimensionedScalar Kax,
    volScalarField Sx,
    volScalarField Shp
)
{
    return Kax * Sx / (Kax + Shp);
}

volScalarField Foam::ADMno1::fShp()
{
    ETempPtrs_[0] = fSion // SvaN
    (
        para_.Ka().va,
        YPtrs_[3],
        ETempPtrs_[6]
    );

    ETempPtrs_[1] = fSion // SbuN
    (
        para_.Ka().bu,
        YPtrs_[4],
        ETempPtrs_[6]
    ); 

    ETempPtrs_[2] = fSion // SproN
    (
        para_.Ka().pro,
        YPtrs_[5],
        ETempPtrs_[6]
    ); 

    ETempPtrs_[3] = fSion // SacN
    (
        para_.Ka().ac,
        YPtrs_[6],
        ETempPtrs_[6]
    );

    ETempPtrs_[4] = fSion // Shco3N
    (
        para_.Ka().co2,
        YPtrs_[9], // SIC
        ETempPtrs_[6]
    ); 

    MPtrs_[1] = fSion // Snh3
    (
        para_.Ka().IN,
        YPtrs_[10], // SIN
        ETempPtrs_[6]
    );

    // calc SohN
    // TODO: the original ADMno1 is quite inconsistent with the dimensions
    // TODO: maybe reverse the dimensionsScalar to scalar in para_?
    volScalarField SohN = para_.Ka().W / ETempPtrs_[6]; 
    SohN.dimensions().reset(ETempPtrs_[6].dimensions()); 
    ETempPtrs_[5] = SohN;

    //     Scat_ - San_ + ShP           - SohN          + (SIN - Snh3)                 
    return Scat_ - San_ + ETempPtrs_[6] - ETempPtrs_[5] + (YPtrs_[10] - MPtrs_[1]) - 
           ETempPtrs_[4] - ETempPtrs_[3]/64 - ETempPtrs_[2]/112 - ETempPtrs_[1]/160 - ETempPtrs_[0]/208;
    //     Shco3N        - SacN/64          - SproN/112         - SbuN/160          - SvaN/208

    // return Scat_ - San_ - ETempPtrs_[4] - ETempPtrs_[3]/64 - ETempPtrs_[2]/112 - ETempPtrs_[1]/160 - ETempPtrs_[0]/208;
           
}

// ************************************************************************* //