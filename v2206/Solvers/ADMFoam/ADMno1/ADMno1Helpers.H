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

// #include "fvCFD.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Inhibition calculations

volScalarField calcInhibition
(
    volScalarField Y, 
    dimensionedScalar denom
);

volScalarField calcInhibitionHP
(
    volScalarField Shp,
    dimensionedScalar UL,
    dimensionedScalar LL,
    dimensionedScalar n
);

//- Kinetic rate calculations

volScalarField calcRho
(
    dimensionedScalar k, 
    volScalarField X
);

volScalarField calcRho
(
    dimensionedScalar k, 
    volScalarField S,
    dimensionedScalar K,
    volScalarField X,
    volScalarField I
);

volScalarField calcRho
(
    dimensionedScalar k, 
    volScalarField S1,
    dimensionedScalar K,
    volScalarField X,
    volScalarField S2,
    volScalarField I
);

//- Components source term calculations



//- Acid-base calculations

volScalarField fSion
(
    dimensionedScalar Kax,
    volScalarField Sx,
    volScalarField Shp
);

volScalarField dfSion
(
    dimensionedScalar Kax,
    volScalarField Sx,
    volScalarField Shp
);


// ************************************************************************* //