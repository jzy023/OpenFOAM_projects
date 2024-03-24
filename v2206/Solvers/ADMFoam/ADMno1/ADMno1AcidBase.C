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

volScalarField Foam::ADMno1::fShp
(
    volScalarField& ShpTemp
)
{
    ETempPtrs_[0] = fSion // SvaN
    (
        para_.Ka().va,
        YPtrs_[3],
        ShpTemp
    );

    ETempPtrs_[1] = fSion // SbuN
    (
        para_.Ka().bu,
        YPtrs_[4],
        ShpTemp
    ); 

    ETempPtrs_[2] = fSion // SproN
    (
        para_.Ka().pro,
        YPtrs_[5],
        ShpTemp
    ); 

    ETempPtrs_[3] = fSion // SacN
    (
        para_.Ka().ac,
        YPtrs_[6],
        ShpTemp
    );

    ETempPtrs_[4] = fSion // Shco3N
    (
        para_.Ka().co2,
        YPtrs_[9], // SIC
        ShpTemp
    ); 

    MPtrs_[1] = fSion // Snh3
    (
        para_.Ka().IN,
        YPtrs_[10], // SIN
        ShpTemp
    );

    // calc SohN
    // TODO: the original ADMno1 is quite inconsistent with the dimensions
    // TODO: maybe reverse the dimensionsScalar to scalar in para_?
    volScalarField SohN = para_.Ka().W / ShpTemp; 
    SohN.dimensions().reset(ShpTemp.dimensions()); 
    ETempPtrs_[5] = SohN;

    //     Scat_ - San_ + ShP           - SohN          + (SIN - Snh3)                 
    return Scat_ - San_ + ShpTemp - ETempPtrs_[5] + (YPtrs_[10] - MPtrs_[1]) - 
    //     Shco3N        - SacN/64            - SproN/112           - SbuN/160            - SvaN/208
           ETempPtrs_[4] - ETempPtrs_[3]/64.0 - ETempPtrs_[2]/112.0 - ETempPtrs_[1]/160.0 - ETempPtrs_[0]/208.0;
           
}

volScalarField Foam::ADMno1::dfShp
(
    volScalarField& ShpTemp
)
{
    volScalarField dSvaN = dfSion
    (
        para_.Ka().va,
        YPtrs_[3],
        ShpTemp
    );

    volScalarField dSbuN = dfSion
    (
        para_.Ka().bu,
        YPtrs_[4],
        ShpTemp
    );

    volScalarField dSproN = dfSion
    (
        para_.Ka().pro,
        YPtrs_[5],
        ShpTemp
    );

    volScalarField dSacN = dfSion
    (
        para_.Ka().ac,
        YPtrs_[6],
        ShpTemp
    );

    volScalarField dShco3N = dfSion
    (
        para_.Ka().co2,
        YPtrs_[9], // SIC
        ShpTemp
    );

    volScalarField dSnh3 = dfSion // Snh3
    (
        para_.Ka().IN,
        YPtrs_[10], // SIN
        ShpTemp
    );

    // calc SohN
    // TODO: the original ADMno1 is quite inconsistent with the dimensions
    // TODO: maybe reverse the dimensionsScalar to scalar in para_?
    volScalarField dSohN = - para_.Ka().W  / (ShpTemp * ShpTemp);
    dSohN.dimensions().reset(dSvaN.dimensions());

    dimensionedScalar uniField
    (
        // YPtrs_[0].dimensions(),
        dimless,
        One
    );

    return uniField - dSnh3 - dShco3N - dSacN/64.0 - dSproN/112.0 - dSbuN/160.0 - dSvaN/208.0 - dSohN;

}

void Foam::ADMno1::calcShp()
{

    //TODO: IO dictionary for these parameters
    scalar tol = 1e-8;
    label nIter = 1e3;
    label i = 0;

    // initial value of x, E and dEdx
    volScalarField x = ETempPtrs_[6];
    volScalarField x0 = ETempPtrs_[6];
    volScalarField E = ETempPtrs_[6];
    volScalarField dE = ETempPtrs_[6];
    
    while
    (
        max(mag(E.field())) > tol &&
        i < nIter
    )
    {
        E.field() = fShp(x).field();
        dE.field() = dfShp(x).field();
        // store old value x
        x0 = x;
        // update x
        x.field() = x0.field() - E.field()/dE.field();
        // false check
        // if (x.field())
        // {
        //     /* code */
        // }
        i++;
    };

    Info << "Newton-Raphson:\tSolving for Sh+, No Interations " << i << endl;
    ETempPtrs_[6] = x;
}


// ************************************************************************* //