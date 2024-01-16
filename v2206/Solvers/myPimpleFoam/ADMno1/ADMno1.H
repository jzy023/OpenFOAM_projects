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
    ADM::ADMno1
    >>> inspired from reactingFoam/CombustionModel<psiReactionThermo>


Description
    Anaerobic Digestion ADModel No.1 class.


\*---------------------------------------------------------------------------*/

#ifndef ADMno1_H
#define ADMno1_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "fvCFD.H"
#include "autoPtr.H"
#include "ODESystem.H"
#include "volFields.H"
#include "simpleMatrix.H"

#include "typeDefs.H"
#include "ADMno1Parameter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class ADMno1
// inherite from other classes
:
    public IOdictionary
{

private:

    //- Thermodynamical condition
    label runMode_;

    //- ADM no1 parameters
    admPara admParameters_;

    //- Reference to the field of specie concentrations
    PtrList<volScalarField> YPtrs_;

    //- List of reaction rate per specie [kg/m3/s]
    PtrList<volScalarField::Internal> RR_;

    //- Temporary concentration fields
    // PtrList<volScalarField>& YPtrsTemps_;
    mutable scalarField c_;
 
    //- Temporary rate-of-change of concentration field
    mutable scalarField dcdt_;

    //- Species names for soluable
    //- !!! Sh2 is calculated iteratively 
    const wordList namesSoluable = {"Ssu", "Saa", "Sfa", "Sva", "Sbu", "Spro",
                                    "Sac", "Sh2", "Sch4", "SIC", "SIN", "SI", "Scat", "San"};

    //- Species names for gaseouse
    const wordList namesGaseous = {"Gh2", "Gch4", "Gco2"};

    //- Species names for particulate
    const wordList namesParticulate = {"Xc", "Xch", "Xpr", "Xli", "Xsu", "Xaa",
                                       "Xfa", "Xc4", "Xpro", "Xac", "Xh2", "XI"};

    //- Species names for acid-base medians
    //- !!! All calculated iteratively
    const wordList namesAcidBase = {"SvaN", "SbuN", "SproN", "SacN", "Shco3N", "SohN"
                                        "Snh3", "ShP"};

    //- calculated by simple subtractions
    const wordList namesMedians = {"Sco2", "Snh4"};


    // const word dictName("");

public:
 
    //- Run time type information
    //> runTimeSelection mechanism: 
    //> https://openfoamwiki.net/index.php/OpenFOAM_guide/runTimeSelection_mechanism
    TypeName("ADMno1");


    // Declare run-time constructor selection table

    // declareRunTimeSelectionTable
    // (
    //     autoPtr,
    //     ADMno1,
    //     fvMesh,
    //     (const fvMesh& mesh, label runMode),
    //     (mesh, runMode)
    // );
 

    //- Default name of the turbulence properties dictionary

    static const word propertiesName;


    // Constructors

    ADMno1
    ( 
        const fvMesh& mesh,
        label runMode 
    );

    ADMno1
    ( 
        const fvMesh& mesh,
        const IOdictionary& ADMno1Dict
    );


    // Selectors

    template<class ADModel>
    static autoPtr<ADModel> New
    (
        const fvMesh&
    );


    static autoPtr<ADMno1> New
    (
        const fvMesh&
    );


    // Destructor

    ~ADMno1() = default;

    // Member Functions
    //- Return YPtrs_
    // inline const PtrList<volScalarField>& Y() const;
    inline PtrList<volScalarField>& Y();


    //- Return const access to the chemical source terms for specie, i
    inline const volScalarField::Internal& RR ( const label i ) const;


    //- Return non const access to chemical source terms [kg/m3/s]
    // virtual volScalarField::Internal& RR ( const label i );


    // combustionModels/laminar.C -> template<class ReactionThermo>
                                  // Foam::tmp<Foam::fvScalarMatrix>
                                  // Foam::combustionModels::laminar<ReactionThermo>::R(volScalarField& Y) const

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#include "ADMno1I.H"

#ifdef NoRepository
    #include "ADMno1Templates.C"
#endif
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //