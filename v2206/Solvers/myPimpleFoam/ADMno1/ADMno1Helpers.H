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

scalar calcInhibition
(
    scalar num,
    scalar denom
);

volScalarField calcInhibitionH
(
    volScalarField* ShpPtr,
    scalar UL,
    scalar LL,
    label n
);

volScalarField calcInhibitionPH
(
    volScalarField* phPtr,
    scalar UL,
    scalar LL,
    label n
);

//- Reaction rate calculations

volScalarField calcRho
(
    scalar k, 
    volScalarField* XPtr
);

volScalarField calcRho
(
    scalar k, 
    volScalarField* SPtr,
    scalar K,
    volScalarField* XPtr,
    scalar I
);

volScalarField calcRho
(
    scalar k, 
    volScalarField* SPtr1,
    scalar K,
    volScalarField* XPtr,
    volScalarField* SPtr2,
    scalar I
);

// ************************************************************************* //