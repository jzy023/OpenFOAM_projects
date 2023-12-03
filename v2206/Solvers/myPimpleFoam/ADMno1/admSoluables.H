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
#include "contiguous.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class List;
 
/*---------------------------------------------------------------------------*\
                           Class Vector Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class admSoluables
:
    public VectorSpace<admSoluables<Cmpt>, Cmpt, 12>
{
public:

    // Typedefs
 
        //- Equivalent type of labels used for valid component indexing
        typedef admSoluables<label> labelAdmSoluables;
 
 
    // Member Constants
 
        //- Rank of Vector is 1
        static constexpr direction rank = 1;
 
 
    //- Component labeling enumeration
    enum components { su_, aa_, fa_, va_, bu_, pro_, ac_, h2_, ch4_, IC_, IN_, I_ };

    // Generated Methods
 
        //- Default construct
        admSoluables() = default;
 
        //- Copy construct
        admSoluables(const admSoluables&) = default;
 
        //- Copy assignment
        admSoluables& operator=(const admSoluables&) = default;

    // Constructors
 
        //- Construct initialized to zero
        inline admSoluables(const Foam::zero);
 
        //- Copy construct from VectorSpace of the same rank
        template<class Cmpt2>
        inline admSoluables(const VectorSpace<admSoluables<Cmpt2>, Cmpt2, 12>& sols);
 
        //- Construct from three components
        inline admSoluables(const Cmpt& Ssu, const Cmpt& Saa, const Cmpt& Sfa,
                            const Cmpt& Sva, const Cmpt& Sbu, const Cmpt& Spro,
                            const Cmpt& Sac, const Cmpt& Sh2, const Cmpt& Sch4,
                            const Cmpt& SIC, const Cmpt& SIN, const Cmpt& SI
                            );
 
        //- Construct from Istream
        inline explicit admSoluables(Istream& is);

    // Member Functions
    
        //- Access to the vector su component
        inline const Cmpt& su() const;
 
        //- Access to the vector aa component
        inline const Cmpt& aa() const;
 
        //- Access to the vector fa component
        inline const Cmpt& fa() const;

        //- Access to the vector va component
        inline const Cmpt& va() const;
 
        //- Access to the vector bu component
        inline const Cmpt& bu() const;

        //- Access to the vector pro component
        inline const Cmpt& pro() const;

        //- Access to the vector ac component
        inline const Cmpt& ac() const;

        //- Access to the vector h2 component
        inline const Cmpt& h2() const;

        //- Access to the vector ch4 component
        inline const Cmpt& ch4() const;

        //- Access to the vector IC component
        inline const Cmpt& IC() const;

        //- Access to the vector IN component
        inline const Cmpt& IN() const;

        //- Access to the vector I component
        inline const Cmpt& I() const;
 
    //- Access to the vector su component
        inline Cmpt& su();
 
        //- Access to the vector aa component
        inline Cmpt& aa();
 
        //- Access to the vector fa component
        inline Cmpt& fa();

        //- Access to the vector va component
        inline Cmpt& va();
 
        //- Access to the vector bu component
        inline Cmpt& bu();

        //- Access to the vector pro component
        inline Cmpt& pro();

        //- Access to the vector ac component
        inline Cmpt& ac();

        //- Access to the vector h2 component
        inline Cmpt& h2();

        //- Access to the vector ch4 component
        inline Cmpt& ch4();

        //- Access to the vector IC component
        inline Cmpt& IC();

        //- Access to the vector IN component
        inline Cmpt& IN();

        //- Access to the vector I component
        inline Cmpt& I();
 
        //- Inplace normalise the vector by its magnitude
        //  For small magnitudes (less than ROOTVSMALL) set to zero.
        //  Will not be particularly useful for a vector of labels
        // inline Vector<Cmpt>& normalise(const scalar tol = ROOTVSMALL);
 
        // //- Inplace removal of components that are collinear to the given
        // //- unit vector.
        // inline Vector<Cmpt>& removeCollinear(const Vector<Cmpt>& unitVec);
 
 
        // //- Return *this (used for point which is a typedef to Vector<scalar>.
        // inline const Vector<Cmpt>& centre
        // (
        //     const UList<Vector<Cmpt>>&  /* (unused) */
        // ) const;


};


// Typedefs
typedef admSoluables<scalar>    admS;
typedef admSoluables<float>     floatAdmS;
typedef admSoluables<double>    doubleAdmS;
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
}

#include "admSoluablesI.H"

#endif

// ************************************************************************* //