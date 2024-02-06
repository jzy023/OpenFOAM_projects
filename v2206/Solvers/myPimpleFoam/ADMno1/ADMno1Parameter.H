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
    ADM::ADMno1Parameter

Description
    H file for declairing AMDno1 parameter class

\*---------------------------------------------------------------------------*/

#ifndef ADMno1Parameter_H
#define ADMno1Parameter_H

// #include "typeDefs.H"
#include "ADMno1Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class admPara
{
private:

    //- Operation modes (thermal conditions for biomass)
    const wordList namesOpMode = {"Meso", "MesoSolid", "Thermo"}; 

    //- Error message
    int errMessage;

    //- Day-to-second conversion
    const scalar DTOS_ = 1.0/86400.0;

    //- Base operational temperature
    dimensionedScalar Tbase_;   

    //- Operational temperature
    dimensionedScalar Top_;   

    //- Decay rate [1/day]
    decayRate kDec_; 

    //- Yield of biomass (Y)
	yieldBiomass yB_;

    //- Yield of product (f)
	yieldProduct yP_;

    //- Carbon content
	carbonContent CC_;	

    //- Inhibition factor for h2 and nh3
	inhibitionParaI KI_;

    //- Inhibition factor
	inhibitionParaS KS_;

	//- Acid-base kinematic para [M/day]
    dimensionedScalar kAB_;

    //- Acid-base equilibrium para [M]
	admAcidBaseEq Ka_;

    //- Henry's Law coeffecients
	gasHenryCoe KH_;

    // Gas transfer rate coefficients [d^-1]
    dimensionedScalar kLa_;  

    // pH Limits
    limPH pHL_;

    //- Nitrogen content
	nitrogenContent NC_;

    // List<scalar> Yini;  // main

    List<scalar> Gini_;   // gas

    List<scalar> Eini_;   // electrolytes

    List<scalar> Mini_;   // medians

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    int printErrMessage(word runMode);

    dimensionedScalar defineTbase(word runMode);

    dimensionedScalar defineTop(word runMode);

	decayRate defineRC(word runMode);

    yieldBiomass defineYields(word runMode);

    inhibitionParaI defineKI(word runMode);

    inhibitionParaS defineKS(word runMode);

	void defineSTOI();

	// void defineAcidBaseDAE();

    void listResizing();

    void defineInitialState(word runMode);

public:

    //- Biochemical rate coefficients (nu_ij)
    List<List<dimensionedScalar>> STOI;

    //- Some acied base stuff
	List<List<dimensionedScalar>> abDAE;

    //- Constructor

    admPara() : 
        errMessage
        (
            printErrMessage("Meso")
        ),
        Tbase_
        (
            defineTbase("Meso")
        ),
        Top_
        (
            defineTop("Meso")
        ),
        kDec_
        (
            defineRC("Meso")
        ),
        yB_
        (
            defineYields("Meso")
        ),
        yP_
        (   
            0.10, 0.25, 0.20,
            0.20, 0.25, 0.95,
            0.19, 0.13, 0.27,
            0.41, 0.06, 0.23,
            0.26, 0.05, 0.40
        ),
        CC_
        (
            0.2786, 0.03, 0.0313,
            0.03, 0.022, 0.03,
            0.0313, 0.03, 0.0217,
            0.025, 0.0268, 0.0313,
            0.0313, 0.024, 0.0156
        ),
        KI_
        (
            defineKI("Meso")
        ),
        KS_
        (
            defineKS("Meso")
        ),
        kAB_
        (
            dimMoles/dimTime, 1e8
        ),
        Ka_
        (
            1.380e-5, 1.514e-5, 1.318e-5,
            1.738e-5, 4.467e-7, 5.623e-10, 1e-14
        ),
        KH_
        (
            7.384654293536963e-04,  // h2
            0.001161902733673,      // ch4
            0.027146692900075,      // co2
            0.031300000000000       // h2o
        ),
    	kLa_
        (
            dimless/dimTime, 200
        ),
        pHL_
        (
            5.5, 4.0, 7.0, 6.0, 6.0, 5.0
        ),
        NC_
        (
            0.007, 0.005714
        )
    {
        defineInitialState("Meso");
    	defineSTOI();
    	// defineAcidBaseDAE();
    };

    admPara(word runMode);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    scalar DTOS() const { return DTOS_; }

    const dimensionedScalar Tbase() const { return Tbase_; }

    const dimensionedScalar Top() const { return Top_; }

	const decayRate kDec() const { return kDec_; }

	const yieldBiomass yB() const { return yB_; }

	const yieldProduct yP() const { return yP_; }

	const carbonContent CC() const { return CC_; }

    const inhibitionParaI KI() const { return KI_; }

    const inhibitionParaS KS() const { return KS_; }

    const dimensionedScalar kAB() const { return kAB_; }

    const admAcidBaseEq Ka() const { return Ka_; }

    const gasHenryCoe KH() const { return KH_; }

    const dimensionedScalar kLa() const { return kLa_; }

    const limPH pHL() const { return pHL_; }

    const nitrogenContent NC() const { return NC_; }

    scalar Gini(label idx) const { return Gini_[idx]; }

    scalar Eini(label idx) const { return Eini_[idx]; }

    scalar Mini(label idx) const { return Mini_[idx]; }
    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //