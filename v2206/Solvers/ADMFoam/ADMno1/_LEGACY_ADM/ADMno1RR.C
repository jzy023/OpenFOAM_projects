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

void Foam::ADMno1::KineticRate(volScalarField& Top)
{

    //- Thermal condition factor
    // TODO: fix dimensions!

    volScalarField TopDummy(Top);

    TopDummy.dimensions().reset(dimless);

    fac_ = (1.0 / para_.Tbase()() - 1.0 / TopDummy) / (100.0 * R_);

    //- Inhibiitons

    IPtrs_[0] = calcInhibitionHP // aa
    (
        EPtrs_.last(), // ShP
        para_.pH_UL_aa, 
        para_.pH_LL_aa,
        nIh_[0]
    );

    IPtrs_[1] = calcInhibitionHP // ac
    (
        EPtrs_.last(), // ShP
        para_.pH_UL_ac, 
        para_.pH_LL_ac,
        nIh_[1]
    );

    IPtrs_[2] = calcInhibitionHP // h2
    (
        EPtrs_.last(), // ShP
        para_.pH_UL_h2, 
        para_.pH_LL_h2,
        nIh_[2]
    );

    // >>> TODO: which one is correct ???
    // IPtrs_[3] = calcInhibition // IN
    // (
    //     YPtrs_[10], // SIN
    //     para_.K_S.IN
    // );

    // TODO: make an overload
    dimensionedScalar KS
    (
        dimMass/dimVolume,
        para_.K_S.IN
    );
    IPtrs_[3] = 1 / (1 + (KS / YPtrs_[10]));
    

	IPtrs_[4] = calcInhibition // h2fa
    (
        YPtrs_[7], // Sh2
        para_.K_I.h2fa
    );

	IPtrs_[5] = calcInhibition // h2c4
    (
        YPtrs_[7], // Sh2
        para_.K_I.h2c4
    );

	IPtrs_[6] = calcInhibition // h2pro
    (
        YPtrs_[7], // Sh2
        para_.K_I.h2pro
    );

	IPtrs_[7] = calcInhibition // nh3
    (
        MPtrs_[1], // Snh3
        para_.K_I.nh3
    );

    //- Kinetic rates

    // KRPtrs_[0] = calcRho
    // (
    //     para_.RC.dis,
    //     YPtrs_[12] // Xc
    // );

    // KRPtrs_[1] = calcRho
    // (
    //     para_.RC.hyd_ch,
    //     YPtrs_[13] // Xch
    // );

    // KRPtrs_[2] = calcRho
    // (
    //     para_.RC.hyd_pr,
    //     YPtrs_[14] // Xpr
    // );

    // KRPtrs_[3] = calcRho
    // (
    //     para_.RC.hyd_li,
    //     YPtrs_[15] // Xli
    // );

    // KRPtrs_[4] = calcRho
    // (
    //     para_.RC.m_su,
    //     YPtrs_[0], // Ssu
    //     para_.K_S.su,
    //     YPtrs_[16], // Xsu
    //     IPtrs_[0] * IPtrs_[3] // Iphaa*IIN
    // );

    // KRPtrs_[5] = calcRho
    // (
    //     para_.RC.m_aa,
    //     YPtrs_[1], // Saa
    //     para_.K_S.aa,
    //     YPtrs_[17], // Xaa
    //     IPtrs_[0] * IPtrs_[3] // Iphaa*IIN
    // );

    // KRPtrs_[6] = calcRho
    // (
    //     para_.RC.m_fa,
    //     YPtrs_[2], // Sfa
    //     para_.K_S.fa,
    //     YPtrs_[18], // Xfa
    //     IPtrs_[0] * IPtrs_[3] * IPtrs_[4] //Iphaa*IIN*Ih2fa
    // );

    // KRPtrs_[7] = calcRho
    // (
    //     para_.RC.m_c4,
    //     YPtrs_[3], // Sva
    //     para_.K_S.c4,
    //     YPtrs_[19], // Xc4
    //     YPtrs_[4],  // Sbu
    //     IPtrs_[0] * IPtrs_[3] * IPtrs_[5] //Iphaa*IIN*Ih2c4
    // );

    // KRPtrs_[8] = calcRho
    // (
    //     para_.RC.m_c4,
    //     YPtrs_[4], // Sbu
    //     para_.K_S.c4,
    //     YPtrs_[19], // Xc4
    //     YPtrs_[3],  // Sva
    //     IPtrs_[0] * IPtrs_[3] * IPtrs_[5] //Iphaa*IIN*Ih2c4
    // );

    // KRPtrs_[9] = calcRho
    // (
    //     para_.RC.m_pro,
    //     YPtrs_[5], // Spro
    //     para_.K_S.pro,
    //     YPtrs_[20], // Xpro
    //     IPtrs_[0] * IPtrs_[3] * IPtrs_[6]  //Iphaa*IIN*Ih2pro
    // );

    // KRPtrs_[10] = calcRho
    // (
    //     para_.RC.m_ac,
    //     YPtrs_[6], // Sac
    //     para_.K_S.ac,
    //     YPtrs_[21], // Xac
    //     IPtrs_[1] * IPtrs_[3] * IPtrs_[7] // Iphac*IIN*Inh3
    // );

	// // >>> in Rosen et al. implementation, no intermediate used for S_h2
	// KRPtrs_[11] = calcRho
    // (
    //     para_.RC.m_h2,
    //     YPtrs_[7], // Sh2
    //     para_.K_S.h2,
    //     YPtrs_[22], // Xh2
    //     IPtrs_[2] * IPtrs_[3] // Iphh2*IIN
    // );

	// KRPtrs_[12] = calcRho
    // (
    //     para_.RC.dec_xsu,
    //     YPtrs_[16] // Xsu
    // );

    // KRPtrs_[13] = calcRho
    // (
    //     para_.RC.dec_xaa,
    //     YPtrs_[17] // Xaa
    // );

    // KRPtrs_[14] = calcRho
    // (
    //     para_.RC.dec_xfa,
    //     YPtrs_[18] // Xfa
    // );

    // KRPtrs_[15] = calcRho
    // (
    //     para_.RC.dec_xc4,
    //     YPtrs_[19] // Xc4
    // );

    // KRPtrs_[16] = calcRho
    // (
    //     para_.RC.dec_xpro,
    //     YPtrs_[20] // Xpro
    // );

    // KRPtrs_[17] = calcRho
    // (
    //     para_.RC.dec_xac,
    //     YPtrs_[21] // Xac
    // );

    // KRPtrs_[18] = calcRho
    // (
    //     para_.RC.dec_xh2,
    //     YPtrs_[22] // Xh2
    // );

    // GRPtrs_[0] = para_.kLa * (YPtrs_[7] - R_ * TopDummy * GPtrs_[0] * para_.KH.h2 * exp(-4180.0 * fac_));

    // GRPtrs_[1] = para_.kLa * (YPtrs_[8] - R_ * TopDummy * GPtrs_[1] * para_.KH.ch4 * exp(-14240.0 * fac_));

    // GRPtrs_[2] = para_.kLa * (YPtrs_[9] - R_ * TopDummy * GPtrs_[2] * para_.KH.co2 * exp(-19410.0 * fac_));

}

// ************************************************************************* //