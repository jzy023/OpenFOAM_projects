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

volScalarField::Internal Foam::ADMno1::fShp
(
    volScalarField::Internal& ShpTemp
)
{
    volScalarField::Internal SvaN = fSion
    (
        para_.Ka().va,
        YPtrs_[3].internalField(),
        ShpTemp
    );

    volScalarField::Internal SbuN = fSion
    (
        para_.Ka().bu,
        YPtrs_[4].internalField(),
        ShpTemp
    ); 

    volScalarField::Internal SproN = fSion
    (
        para_.Ka().pro,
        YPtrs_[5].internalField(),
        ShpTemp
    ); 

    volScalarField::Internal SacN = fSion
    (
        para_.Ka().ac,
        YPtrs_[6].internalField(),
        ShpTemp
    );

    volScalarField::Internal Shco3N = fSion
    (
        para_.Ka().co2,
        YPtrs_[9].internalField(), // SIC
        ShpTemp
    ); 

    // Sco2
    MPtrs_[0].ref() = YPtrs_[9].internalField() - Shco3N; 

    // Snh3
    MPtrs_[1].ref() = fSion
    (
        para_.Ka().IN,
        YPtrs_[10].internalField(), // SIN
        ShpTemp
    );

    // calc SohN
    // TODO: the original ADMno1 is quite inconsistent with the dimensions
    // TODO: maybe reverse the dimensionsScalar to scalar in para_?
    volScalarField::Internal SohN = para_.Ka().W / ShpTemp; 
    SohN.dimensions().reset(ShpTemp.dimensions()); 
    
        // Scat_ - San_ + ShP     - SohN + (SIN                        - Snh3)                 
    return Scat_ - San_ + ShpTemp - SohN + (YPtrs_[10].internalField() - MPtrs_[1].internalField()) - 
           Shco3N - SacN/64.0 - SproN/112.0 - SbuN/160.0 - SvaN/208.0;
           
}

volScalarField::Internal Foam::ADMno1::dfShp
(
    volScalarField::Internal& ShpTemp
)
{
    volScalarField::Internal dSvaN = dfSion
    (
        para_.Ka().va,
        YPtrs_[3].internalField(),
        ShpTemp
    );

    volScalarField::Internal dSbuN = dfSion
    (
        para_.Ka().bu,
        YPtrs_[4].internalField(),
        ShpTemp
    );

    volScalarField::Internal dSproN = dfSion
    (
        para_.Ka().pro,
        YPtrs_[5].internalField(),
        ShpTemp
    );

    volScalarField::Internal dSacN = dfSion
    (
        para_.Ka().ac,
        YPtrs_[6].internalField(),
        ShpTemp
    );

    volScalarField::Internal dShco3N = dfSion
    (
        para_.Ka().co2,
        YPtrs_[9].internalField(), // SIC
        ShpTemp
    );

    volScalarField::Internal dSnh3 = dfSion // Snh3
    (
        para_.Ka().IN,
        YPtrs_[10].internalField(), // SIN
        ShpTemp
    );

    // calc SohN
    // TODO: the original ADMno1 is quite inconsistent with the dimensions
    // TODO: maybe reverse the dimensionsScalar to scalar in para_?
    volScalarField::Internal dSohN = - para_.Ka().W / (ShpTemp * ShpTemp);
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
    scalar tol = 1e-16;
    label nIter = 1e3;
    label i = 0;

    // initial value of x, E and dEdx
    volScalarField::Internal x = ShP_;   // x = Shp
    volScalarField::Internal E = ShP_;   // E = dShp/dt
    volScalarField::Internal dE = ShP_;  // dE = (dShp/dt)/dShp
    
    do
    {
        E.field() = fShp(x).field();
        dE.field() = dfShp(x).field();
        x.field() = x.field() - E.field()/dE.field();
        // false check
        if( min(x.field()) < 0 )
        {
            std::cerr << nl << "--> FOAM FATAL IO ERROR:" << nl
                      << "Proton (H+) concentration below Zero\n";
            std::exit(1);
        }
        i++;
    }
    while
    (
        max(mag(E.field())) > tol &&
        i < nIter
    );

    Info << "Newton-Raphson:\tSolving for Sh+, No Interations " << i // << endl;
         << ", min Shp: " << min(x.field()) << ", max Shp: " << max(x.field()) << endl;

    // ShP
    ShP_ = x;

    // Sco2
    MPtrs_[0].ref() = YPtrs_[9].internalField() - fSion
    (
        para_.Ka().co2,
        YPtrs_[9].internalField(), // SIC
        ShP_
    ); 

    // Snh3
    MPtrs_[1].ref() = fSion
    (
        para_.Ka().IN,
        YPtrs_[10].internalField(), // SIN
        ShP_
    );
}


// ************************************************************************* //