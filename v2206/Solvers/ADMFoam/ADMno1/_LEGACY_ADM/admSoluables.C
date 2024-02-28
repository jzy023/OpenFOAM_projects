/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

Class
    ADM::admSoluables

Description
    C file for declairing AMDno1 class

\*---------------------------------------------------------------------------*/

#include "scalar.H"
#include "admSoluables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::admSoluables<scalar>::vsType::typeName = "admS";

template<>
const char* const Foam::Vector<float>::vsType::typeName = "floatAdmS";

template<>
const char* const Foam::Vector<double>::vsType::typeName = "doubleAdmS";

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#undef  defineTraits
#define defineTraits(Type, Prefix)                                            \
                                                                              \
    template<>                                                                \
    const char* const Foam::admSoluables<Type>::vsType::componentNames[] =    \
    {                                                                         \
        "su_", "aa_", "fa_", "va_", "bu_", "pro_",                            \
        "ac_", "h2_", "ch4_", "IC_", "IN_",  "I_"                             \
    };                                                                        \

// defineTraits(float, floatScalar);
// defineTraits(double, doubleScalar);

#undef defineTraits

// ************************************************************************* //