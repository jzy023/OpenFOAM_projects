/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
 
#include "adModel.H" 

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class model>
Foam::autoPtr<model> Foam::adModel::New
(
    const fvMesh& mesh,
    label runMode
)
{
    // TODO: do it properly!!! with virtual destructors and constructor hash tables 
    // new keywaord is not gonna last!

    // IOdictionary thermoDict
    // (
    //     IOobject
    //     (
    //         phasePropertyName(dictName, phaseName),
    //         mesh.time().constant(),
    //         mesh,
    //         IOobject::MUST_READ_IF_MODIFIED,
    //         IOobject::NO_WRITE,
    //         false
    //     )
    // );
 
    // auto* ctorPtr = getThermoOrDie<Thermo>
    // (
    //     thermoDict,
    //     *(Thermo::fvMeshConstructorTablePtr_)
    // );

    auto* reactionPtr = new model(mesh, runMode);
    return autoPtr<model>(reactionPtr);
}

// ************************************************************************* //