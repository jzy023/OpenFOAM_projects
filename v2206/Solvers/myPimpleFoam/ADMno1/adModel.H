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
    ADM::adModel

Description
    Anaerobic Digestion Model No.1 class.


\*---------------------------------------------------------------------------*/

#ifndef adModel_H
#define adModel_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "ODESystem.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "adModelParameter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// namespace ADM
// {

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class adModel
{

private:

    //- Reference to the field of specie concentrations
    PtrList<volScalarField>& YPtrs_;

    //- List of reaction rate per specie [kg/m3/s]
    PtrList<volScalarField::Internal> RR_;

    //- Temporary concentration field
    mutable scalarField c_;
 
    //- Temporary rate-of-change of concentration field
    mutable scalarField dcdt_;

public:
 
    //- Run time type information
    TypeName("adModel");

    // Constructors


    // Destructor


    // Member Functions
    //- Return const access to the chemical source terms for specie, i
    inline const volScalarField::Internal& RR ( const label i ) const;

    //- Return non const access to chemical source terms [kg/m3/s]
    virtual volScalarField::Internal& RR ( const label i );

    // combustionModels/laminar.C -> template<class ReactionThermo>
                                  // Foam::tmp<Foam::fvScalarMatrix>
                                  // Foam::combustionModels::laminar<ReactionThermo>::R(volScalarField& Y) const

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //