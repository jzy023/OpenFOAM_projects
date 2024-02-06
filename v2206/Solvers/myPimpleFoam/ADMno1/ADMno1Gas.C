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

void Foam::ADMno1::GasPhaseRate(volScalarField& Top)
{
    //- Thermal condition factor

    volScalarField TopDummy(Top);

    TopDummy.dimensions().reset(dimless);

    fac_ = (1.0 / para_.Tbase().value() - 1.0 / TopDummy) / (100.0 * R_);

    GRPtrs_[0] = para_.kLa() * (YPtrs_[7] - R_ * TopDummy * GPtrs_[0] * para_.KH().h2 * exp(-4180.0 * fac_));

    GRPtrs_[1] = para_.kLa() * (YPtrs_[8] - R_ * TopDummy * GPtrs_[1] * para_.KH().ch4 * exp(-14240.0 * fac_));

    GRPtrs_[2] = para_.kLa() * (YPtrs_[9] - R_ * TopDummy * GPtrs_[2] * para_.KH().co2 * exp(-19410.0 * fac_));

}


void Foam::ADMno1::GasExitRate(volScalarField& Top)
{

    // field of cell volume for mesh 
    scalarField volMeshField = GPtrs_[0].mesh().V().field();            

    // particle scaled gas volume
    scalarField volGas = volMeshField * fracG_; 

    // particle scaled liquid volume
    scalarField volLiq = volMeshField * fracL_;

    // particle scaled pipe resistance
    volScalarField kp = 0.0*fac_;

    kp.dimensions().reset(dimVolume/dimTime/dimPressure);

    kp.field() = KP_ * (volMeshField / (volGas + volLiq));  

    //- gas pressure
    volScalarField TopDummy(Top);

    TopDummy.dimensions().reset(dimless); 

    volScalarField Ph2o = GPtrs_[0];

    Ph2o.field() = para_.KH().h2 * exp(5290.0 * fac_ * 100 * R_);

    volScalarField Pgas = (GPtrs_[0] / 16.0 + GPtrs_[1] / 64.0 + GPtrs_[2]) * R_ * TopDummy + Ph2o;

    Pgas.dimensions().reset(dimPressure);

    // TODO: 1.013 -> atmospheric pressure!
    dimensionedScalar Patm( dimPressure, 1.013 );

    volScalarField qGasLocal = kp * (Pgas - Patm); 

	qGasLocal.min(0.0);

    dGPtrs_[0].field() = GRPtrs_[0].field() * volLiq / volGas - GPtrs_[0].field() * qGasLocal.field() / volGas;

    dGPtrs_[1].field() = GRPtrs_[1].field() * volLiq / volGas - GPtrs_[1].field() * qGasLocal.field() / volGas;

    dGPtrs_[2].field() = GRPtrs_[2].field() * volLiq / volGas - GPtrs_[2].field() * qGasLocal.field() / volGas;
    
};


// ************************************************************************* //