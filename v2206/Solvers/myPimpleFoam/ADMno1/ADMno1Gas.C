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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ADMno1::gasR(volScalarField& Top)
{
    fac_ = (1.0 / para_.getTbase() - 1.0 / Top) / (100.0 * R_);

    GRPtrs_[0] = para_.kLa * (YPtrs_[7] - R_ * Top * GPtrs_[0] * para_.KH.h2 * exp(-4180.0 * fac_));

    GRPtrs_[1] = para_.kLa * (YPtrs_[8] - R_ * Top * GPtrs_[1] * para_.KH.ch4 * exp(-14240.0 * fac_));

    GRPtrs_[2] = para_.kLa * (YPtrs_[9] - R_ * Top * GPtrs_[2] * para_.KH.co2 * exp(-19410.0 * fac_));

}

// ************************************************************************* //