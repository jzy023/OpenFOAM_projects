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

#include "typeDefs.H"
#include "ADMno1Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class admPara
{
private:

    scalar Tbase_;      // base operational temperature

    // yieldBiomassTest yB;
	yieldBiomass yB_;	// yield of biomass (Y)
	yieldProduct yP_;	// yield of product (f)
	carbonContent CC_;	// carbon content

	scalar N_aa_;  // Nitrogen content
	scalar N_bac_; // Nitrogen content Biomass

	void defineRCMeso();
	void defineRCMesoSolid();
	void defineRCThermo();

	void defineYieldsMeso();
	void defineYieldsThermo();

    void defineKI_Meso();

	void defineKIMeso();
	void defineKIMesoSolid();
	void defineKIThermo();

	void defineKSMeso();
	void defineKSMesoSolid();
	void defineKSThermo();

	void defineSTOI();
	void defineAcidBaseDAE();

    List<scalar> Yini;  // main
    List<scalar> Gini;   // gas
    List<scalar> Eini;   // electrolytes
    List<scalar> Mini;   // medians

    void listResizing();

    void defineInitialStateMeso();      //TODO: need work
    void defineInitialStateMesoSolid();
    void defineInitialStateThermo();


public:
	std::vector<std::vector<scalar>> STOI;
	std::vector<std::vector<scalar>> abDAE;
	reactionRateConstant RC; // k_dis; k_dec; k_hyd; k_m;

	inhibitionParaI K_I;	// 
	inhibitionParaS K_S;	// 
	admAcidBaseKinetic kAB;	// 
	admAcidBaseEq Ka;		// 
	gasHenryCoe KH;			// Henry's Law
    scalar kLa;       	    // gas transfer rate coefficients [d^-1]

    // pH Limits
	scalar pH_UL_aa;
	scalar pH_LL_aa;

	scalar pH_UL_ac;
	scalar pH_LL_ac;

	scalar pH_UL_h2;
	scalar pH_LL_h2;

    gasHenryCoeTest test;
    gasHenryCoeTest defineTest()
    {
        return gasHenryCoeTest(0.1, 0.1, 0.1, 0.1);
    };

    admPara() : 
        Tbase_(35),
        // Yields of product
        yP_({0.10, 0.25, 0.20,
             0.20, 0.25, 0.95,
             0.19, 0.13, 0.27,
             0.41, 0.06, 0.23,
             0.26, 0.05, 0.40}),
        // Carbon Content
        CC_({0.2786, 0.03, 0.0313,
             0.03, 0.022, 0.03,
             0.0313, 0.03, 0.0217,
             0.025, 0.0268, 0.0313,
             0.0313, 0.024, 0.0156}),
        // Nitrogen Content
        N_aa_(0.007), N_bac_(0.005714),
        // Acid Base Kinetics
        kAB({1e10, 1e10, 1e10,
             1e10, 1e10, 1e10}),
        // Acide base Equilibrium Para
        Ka({1.380e-5, 1.514e-5, 1.318e-5,
            1.738e-5, 4.467e-7, 5.623e-10,
            1e-14}),
        // Henry's Law Coefficients
        KH({7.384654293536963e-04, // h2
            0.001161902733673,     // ch4
            0.027146692900075,     // co2
            0.031300000000000}),   // h2o

        // Gas Transfer Coefficients
        // kL({200, 200, 200}),
    	kLa(200),
        // pH bounds
        pH_UL_aa(5.5), pH_LL_aa(4),
        pH_UL_ac(7), pH_LL_ac(6),
        pH_UL_h2(6), pH_LL_h2(5),

        test(defineTest())
    {
    	defineRCMeso();
    	defineYieldsMeso();
    	defineKIMeso();
    	defineKSMeso();
    	defineSTOI();
    	defineAcidBaseDAE();
    };

	// member functions

    void setOpMode(label idxMode);

	const reactionRateConstant getRC() const { return RC; }

	const yieldBiomass getYieldBio() const { return yB_; }

	const yieldProduct getYieldPro() const { return yP_; }

	const carbonContent getCarbonCon() const { return CC_; }

    scalar getTbase() const { return Tbase_; }

    scalar getYini(label idx) const { return Yini[idx]; }
    scalar getGini(label idx) const { return Gini[idx]; }
    scalar getEini(label idx) const { return Eini[idx]; }
    scalar getMini(label idx) const { return Mini[idx]; }
    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //