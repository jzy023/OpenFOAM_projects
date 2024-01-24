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

void Foam::ADMno1::RR()
{
    //- Inhibiitons

    scalar n_aa = 3.0 / (para_.pH_UL_aa - para_.pH_LL_aa);
    scalar n_ac = 3.0 / (para_.pH_UL_ac - para_.pH_LL_ac);
    scalar n_h2 = 3.0 / (para_.pH_UL_h2 - para_.pH_LL_h2);

    IPtrs_[0] = calcInhibitionHP // aa
    (
        EPtrs_.last(), // ShP
        para_.pH_UL_aa, 
        para_.pH_LL_aa,
        n_aa
    );

    IPtrs_[1] = calcInhibitionHP // ac
    (
        EPtrs_.last(), // ShP
        para_.pH_UL_ac, 
        para_.pH_LL_ac,
        n_ac
    );

    IPtrs_[2] = calcInhibitionHP // h2
    (
        EPtrs_.last(), // ShP
        para_.pH_UL_h2, 
        para_.pH_LL_h2,
        n_h2
    );

    // >>> TODO: which one is correct ???
    // IPtrs_[3] = calcInhibition // IN
    // (
    //     YPtrs_[10], // SIN
    //     para_.K_S.IN
    // );

    IPtrs_[3] = 1 / (1 + (para_.K_S.IN / YPtrs_[10]));
    

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

    //- Raction rates

    RRPtrs_[0] = calcRho
    (
        para_.RC.dis,
        YPtrs_[12] // Xc
    );

    RRPtrs_[1] = calcRho
    (
        para_.RC.hyd_ch,
        YPtrs_[13] // Xch
    );

    RRPtrs_[2] = calcRho
    (
        para_.RC.hyd_pr,
        YPtrs_[14] // Xpr
    );

    RRPtrs_[3] = calcRho
    (
        para_.RC.hyd_li,
        YPtrs_[15] // Xli
    );

    RRPtrs_[4] = calcRho
    (
        para_.RC.m_su,
        YPtrs_[0], // Ssu
        para_.K_S.su,
        YPtrs_[16], // Xsu
        IPtrs_[0] * IPtrs_[3] // Iphaa*IIN
    );

    RRPtrs_[5] = calcRho
    (
        para_.RC.m_aa,
        YPtrs_[1], // Saa
        para_.K_S.aa,
        YPtrs_[17], // Xaa
        IPtrs_[0] * IPtrs_[3] // Iphaa*IIN
    );

    RRPtrs_[6] = calcRho
    (
        para_.RC.m_fa,
        YPtrs_[2], // Sfa
        para_.K_S.fa,
        YPtrs_[18], // Xfa
        IPtrs_[0] * IPtrs_[3] * IPtrs_[4] //Iphaa*IIN*Ih2fa
    );

    RRPtrs_[7] = calcRho
    (
        para_.RC.m_c4,
        YPtrs_[3], // Sva
        para_.K_S.c4,
        YPtrs_[19], // Xc4
        YPtrs_[4],  // Sbu
        IPtrs_[0] * IPtrs_[3] * IPtrs_[5] //Iphaa*IIN*Ih2c4
    );

    RRPtrs_[8] = calcRho
    (
        para_.RC.m_c4,
        YPtrs_[4], // Sbu
        para_.K_S.c4,
        YPtrs_[19], // Xc4
        YPtrs_[3],  // Sva
        IPtrs_[0] * IPtrs_[3] * IPtrs_[5] //Iphaa*IIN*Ih2c4
    );

    RRPtrs_[9] = calcRho
    (
        para_.RC.m_pro,
        YPtrs_[5], // Spro
        para_.K_S.pro,
        YPtrs_[20], // Xpro
        IPtrs_[0] * IPtrs_[3] * IPtrs_[6]  //Iphaa*IIN*Ih2pro
    );

    RRPtrs_[10] = calcRho
    (
        para_.RC.m_ac,
        YPtrs_[6], // Sac
        para_.K_S.ac,
        YPtrs_[21], // Xac
        IPtrs_[1] * IPtrs_[3] * IPtrs_[7] // Iphac*IIN*Inh3
    );

	// >>> in Rosen et al. implementation, no intermediate used for S_h2
	RRPtrs_[11] = calcRho
    (
        para_.RC.m_h2,
        YPtrs_[7], // Sh2
        para_.K_S.h2,
        YPtrs_[22], // Xh2
        IPtrs_[2] * IPtrs_[3] // Iphh2*IIN
    );

	RRPtrs_[12] = calcRho
    (
        para_.RC.dec_xsu,
        YPtrs_[16] // Xsu
    );

    RRPtrs_[13] = calcRho
    (
        para_.RC.dec_xaa,
        YPtrs_[17] // Xaa
    );

    RRPtrs_[14] = calcRho
    (
        para_.RC.dec_xfa,
        YPtrs_[18] // Xfa
    );

    RRPtrs_[15] = calcRho
    (
        para_.RC.dec_xc4,
        YPtrs_[19] // Xc4
    );

    RRPtrs_[16] = calcRho
    (
        para_.RC.dec_xpro,
        YPtrs_[20] // Xpro
    );

    RRPtrs_[17] = calcRho
    (
        para_.RC.dec_xac,
        YPtrs_[21] // Xac
    );

    RRPtrs_[18] = calcRho
    (
        para_.RC.dec_xh2,
        YPtrs_[22] // Xh2
    );

}

// ************************************************************************* //