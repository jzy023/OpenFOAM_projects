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

    //- Operation mode
    word opMode_;

    //- Schmidt number
    scalar Sc_;

    //- ADM no1 parameters
    admPara para_;

    //- Ihibition factors
    scalarList nIh_ {3};

    //- Gas constant
    scalar R_;

    //- Thermal condition factor
    volScalarField fac_;				

    //- pH value
    volScalarField pH_;

    //- List of pointers to the field of specie concentrations
    PtrList<volScalarField> YPtrs_;

    //- List of pointers to the field of temperal specie concentrations
    PtrList<scalarField> YTempPtrs_;

    //- List of pointers to the field of specie concentrations derivatives
    PtrList<volScalarField::Internal> dYPtrs_;

    //- List of pointers to the field of gas concentrations
    PtrList<volScalarField> GPtrs_;

    //- List of pointers to the field of electrolytes concentrations
    PtrList<volScalarField> EPtrs_;

    //- List of pointers to the field of medians concentrations
    PtrList<volScalarField> MPtrs_;

    //- List of pointers to inhibition 
    PtrList<volScalarField> IPtrs_;

    //- List of pointers to kinetic rate per specie [kg/m3/s] and [mol/m3/s]
    PtrList<volScalarField> KRPtrs_;

    //- List of pointers to gas trasfer rate
    PtrList<volScalarField> GRPtrs_;

    //- Temporary concentration fields
    // PtrList<volScalarField>& YPtrsTemps_;
    mutable scalarField c_;
 
    //- Temporary rate-of-change of concentration field
    mutable scalarField dcdt_;

    //- Operation modes (thermal conditions for biomass)
    const wordList namesOpMode = {"Meso", "MesoSolid", "Thermo"};  

    //- Species names for soluable
    //- !!! Sh2 is calculated iteratively 
    const wordList namesSoluable = {"Ssu", "Saa", "Sfa", "Sva", "Sbu", "Spro",
                                    "Sac", "Sh2", "Sch4", "SIC", "SIN", "SI"};

    //- Species names for particulate
    const wordList namesParticulate = {"Xc", "Xch", "Xpr", "Xli", "Xsu", "Xaa",
                                       "Xfa", "Xc4", "Xpro", "Xac", "Xh2", "XI"};

    //- Species names for gaseouse
    const wordList namesGaseous = {"Gh2", "Gch4", "Gco2"};

    //- Species names for acid-base medians
    //- !!! All calculated iteratively
    const wordList namesElectrolyte = {"SvaN", "SbuN", "SproN", "SacN", "Shco3N", "SohN",
                                       "ShP"};

    //- calculated by simple subtractions
    const wordList namesMedians = {"Sco2", "Snh3", "Snh4"};

    //- inert substances
    const wordList namesInerts = {"Scat", "San"};

    //- Cation concentration (inert)
    scalar Scat_;

    //- Anion concentration (inert)
    scalar San_;

    //- Pirvate Helper Functions
    #include "ADMno1Helpers.H"

    //- Calculate kinetic rate
    void kineticRate(volScalarField& Top);

    //-
    void calcSh2();


    


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
    //     (const fvMesh& mesh, word runMode),
    //     (mesh, runMode)
    // );


    //- Default name of the turbulence properties dictionary

    static const word propertiesName;


    // Constructors

    ADMno1
    ( 
        const fvMesh& mesh,
        word runMode 
    );

    ADMno1
    ( 
        const fvMesh& mesh,
        const IOdictionary& ADMno1Dict
    );


    // Selectors

    // template<class ADModel>
    // static autoPtr<ADModel> New (const fvMesh&);

    static autoPtr<ADMno1> New (const fvMesh&);


    // Destructor

    ~ADMno1() = default;


    // Inline Member Functions
    //- Return Schmidt number
    inline scalar Sc();

    //- Return a pointer of admno1 parameters
    inline admPara* paraPtr();
    
    //- Return YPtrs_
    inline PtrList<volScalarField>& Y();

    //- Return dYPtrs_
    inline PtrList<volScalarField::Internal>& dY();


    // Inline Member Functions
    //- Print Error messages
    void modeCheckErr();

    //- Correct components consumption rate
    void correct(volScalarField& Top);

    //- Return components consumption rate matrix, i.e. source term.
    tmp<fvScalarMatrix> R(label i, volScalarField& Yi) const;
    // void R(label i, volScalarField& Yi) const;

   
    //- Return const access to the chemical source terms for specie, i
    // inline const volScalarField::Internal& RR ( const label i ) const;


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
// #ifdef NoRepository
//     #include "ADMno1Templates.C"
// #endif
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //