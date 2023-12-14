/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::propertiesReader
 
Description
    ADM properties reader for OpenFOAM format
 
SourceFiles
    propertiesReader.C
 
\*---------------------------------------------------------------------------*/
 
#ifndef propertiesReader_H
#define propertiesReader_H

#include "typeInfo.H"
#include "specieElement.H"
#include "Reaction.H"
#include "ReactionList.H"
#include "runTimeSelectionTables.H"
#include "fileName.H"
#include "labelList.H"
#include "speciesTable.H"
#include "atomicWeights.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
// namespace Foam
// {
 
/*---------------------------------------------------------------------------*\
                       Class foamChemistry Declaration
\*---------------------------------------------------------------------------*/

class propertiesReader
{

    //- runMode of ADMno1
    label mode_;

    //- Table of species
    speciesTable& speciesTable_;

    //- Table of species composition
    speciesCompositionTable speciesComposition_;

public:
 
    //- Runtime type information
    TypeName("chemistryReader");

// Constructors
 
    //- Construct from foamChemistry and thermodynamics file names
    propertiesReader
    (
        const fileName& reactionsFileName,
        speciesTable& species,
        const fileName& thermoFileName
    );

    //- Construct by getting foamChemistry and thermodynamics file names
    //- from dictionary
    propertiesReader
    (
        const dictionary& thermoDict,
        speciesTable& species
    );

// Member Functions

    //- ADMno1 run mode
    const label mode() const {
        return mode_;
    }

    //- Table of species
    const speciesTable& species() const
    {
        return speciesTable_;
    }

    //- Table of species composition
    const speciesCompositionTable& specieComposition() const
    {
        return speciesComposition_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
// } // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* //