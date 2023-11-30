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
    H file for declairing AMDno1 class

\*---------------------------------------------------------------------------*/

#ifndef admSoluables_H
#define admSoluables_H

#include "typeDefs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class admSoluables
{
private:
    IOdictionary    dict_; // might not nedded
    scalar          soluables_[12];
public:
    admSoluables
    (
        const IOdictionary  &io, 
        const double        soluables[12]
    )
    :
        dict_(io)
    {
        for (int i = 0; i < 12; ++i) {
            soluables_[i] = soluables[i];
        }
    };
};

namespace Foam
{
template<>
class pTraits<admSoluables>
{
    admSoluables p_;
public:
        // Appropriate public typedefs such as cmptType, etc. 
        // Look into the Scalar.H and see how the partial specialization is 
        // done for the Scalar, and do the same for your type. 
    //- Equivalent type of labels used for valid component indexing
    typedef admSoluables cmptType;

    typedef label labelType;
 
    // Member Constants
 
    //- Dimensionality of space
    static constexpr direction dim = 3;

    //- Rank of Scalar is 0
    static constexpr direction rank = 0;

    //- Number of components in Scalar is 1
    static constexpr direction nComponents = 1;
};
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //