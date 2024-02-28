/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
 
#ifndef makeAdModelTypes_H
#define makeAdModelTypes_H
 
#include "addToRunTimeSelectionTable.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
// namespace Foam
// {
 
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
// #define typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie)          \
//                                                                                \
//     typedef                                                                    \
//         Transport                                                              \
//         <                                                                      \
//             species::thermo                                                    \
//             <                                                                  \
//                 Thermo                                                         \
//                 <                                                              \
//                     EqnOfState                                                 \
//                     <                                                          \
//                         Specie                                                 \
//                     >                                                          \
//                 >,                                                             \
//                 Type                                                           \
//             >                                                                  \
//         >                                                                      \
//         Transport##Type##Thermo##EqnOfState##Specie
 

// #define addThermoPhysicsThermo(BaseThermo,CThermoMixtureThermoPhys)            \
//                                                                                \
//     addToRunTimeSelectionTable                                                 \
//     (                                                                          \
//         BaseThermo,                                                            \
//         CThermoMixtureThermoPhys,                                              \
//         fvMesh                                                                 \
//     );                                                                         \
//                                                                                \
//     addToRunTimeSelectionTable                                                 \
//     (                                                                          \
//         BaseThermo,                                                            \
//         CThermoMixtureThermoPhys,                                              \
//         fvMeshDictPhase                                                        \
//     );
 
 
// #define makeThermoPhysicsThermo(BaseThermo,CThermo,Mixture,ThermoPhys)         \
//                                                                                \
//     defineThermoPhysicsThermo(BaseThermo, CThermo, Mixture, ThermoPhys);       \
//                                                                                \
//     addThermoPhysicsThermo(BaseThermo, CThermo##Mixture##ThermoPhys)
 
 
// #define makeThermoPhysicsThermos(BaseThermo,CThermo,Mixture,ThermoPhys)        \
//                                                                                \
//     defineThermoPhysicsThermo(BaseThermo, CThermo, Mixture, ThermoPhys);       \
//                                                                                \
//                                                                                \
//     addThermoPhysicsThermo(basicThermo, CThermo##Mixture##ThermoPhys);         \
//     addThermoPhysicsThermo(fluidThermo, CThermo##Mixture##ThermoPhys);         \
//     addThermoPhysicsThermo(BaseThermo, CThermo##Mixture##ThermoPhys)
 
 
// #define makeThermo(BaseThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
//                                                                                \
//     typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);             \
//                                                                                \
//     makeThermoPhysicsThermo                                                    \
//     (                                                                          \
//         BaseThermo,                                                            \
//         CThermo,                                                               \
//         Mixture,                                                               \
//         Transport##Type##Thermo##EqnOfState##Specie                            \
//     )
 
 
// #define makeThermos(BaseThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
//                                                                                \
//     typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);             \
//                                                                                \
//     makeThermoPhysicsThermos                                                   \
//     (                                                                          \
//         BaseThermo,                                                            \
//         CThermo,                                                               \
//         Mixture,                                                               \
//         Transport##Type##Thermo##EqnOfState##Specie                            \
//     )
 
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
// } // End namespace Foam
 
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* //