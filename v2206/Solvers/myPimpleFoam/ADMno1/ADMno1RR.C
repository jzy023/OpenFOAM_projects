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

    // IPtrs_[0] = calcInhibitionH // aa
    // (
    //     &(EPtrs_.last()), // ShP
    //     para_.pH_UL_aa, 
    //     para_.pH_LL_aa, n_aa
    // );

    // IPtrs_[1] = calcInhibitionH // ac
    // (
    //     &(EPtrs_.last()), // ShP
    //     para_.pH_UL_ac, 
    //     para_.pH_LL_ac, n_ac
    // );

    // IPtrs_[2] = calcInhibitionH // h2
    // (
    //     &(EPtrs_.last()), // ShP
    //     para_.pH_UL_h2, 
    //     para_.pH_LL_h2, n_h2
    // );

    // // >>> ???
    // IPtrs_[3] = calcInhibition // IN
    // (
    //     &(YPtrs_[10]), // SIN
    //     para_.K_S.IN
    // );

	// IPtrs_[4] = calcInhibition // h2fa
    // (
    //     &(YPtrs_[7]), // Sh2
    //     para_.K_I.h2fa
    // );

	// IPtrs_[5] = calcInhibition // h2c4
    // (
    //     &(YPtrs_[7]), // Sh2
    //     para_.K_I.h2c4
    // );

	// IPtrs_[6] = calcInhibition // h2pro
    // (
    //     &(YPtrs_[7]), // Sh2
    //     para_.K_I.h2pro
    // );

	// IPtrs_[7] = calcInhibition // nh3
    // (
    //     &(MPtrs_[1]), // Snh3
    //     para_.K_I.nh3
    // );



    // //- Raction rates
    // RRPtrs_[0] = calcRho(para_.RC.dis, Conc_.X_c);
    // RRPtrs_[1] = calcRho(para_.RC.hyd_ch, Conc_.X_ch);
    // RRPtrs_[2] = calcRho(para_.RC.hyd_pr, Conc_.X_pr);
    // RRPtrs_[3] = calcRho(para_.RC.hyd_li, Conc_.X_li);

    // RRPtrs_[4] = calcRho(para_.RC.m_su, Conc_.S_su, para_.K_S.su, Conc_.X_su, I.phaa * I.IN);
    // RRPtrs_[5] = calcRho(para_.RC.m_aa, Conc_.S_aa, para_.K_S.aa, Conc_.X_aa, I.phaa * I.IN);
    // RRPtrs_[6] = calcRho(para_.RC.m_fa, Conc_.S_fa, para_.K_S.fa, Conc_.X_fa, I.phaa * I.IN * I.h2fa);
    // RRPtrs_[7] = calcRho(para_.RC.m_c4, Conc_.S_va, para_.K_S.c4, Conc_.X_c4, Conc_.S_bu, I.phaa * I.IN * I.h2c4);
    // RRPtrs_[8] = calcRho(para_.RC.m_c4, Conc_.S_bu, para_.K_S.c4, Conc_.X_c4, Conc_.S_va, I.phaa * I.IN * I.h2c4);
    // RRPtrs_[9] = calcRho(para_.RC.m_pro, Conc_.S_pro, para_.K_S.pro, Conc_.X_pro, I.phaa * I.IN * I.h2pro);
    // RRPtrs_[10] = calcRho(para_.RC.m_ac, Conc_.S_ac, para_.K_S.ac, Conc_.X_ac, I.phac * I.IN * I.nh3);

	// // >>> in Rosen et al. implementation, no intermediate used for S_h2
	// RRPtrs_[11] = calcRho(para_.RC.m_h2, Conc_.S_h2, para_.K_S.h2, Conc_.X_h2, I.phh2 * I.IN);

	// RRPtrs_[12] = calcRho(para_.RC.dec_xsu, Conc_.X_su);
    // RRPtrs_[13] = calcRho(para_.RC.dec_xaa, Conc_.X_aa);
    // RRPtrs_[14] = calcRho(para_.RC.dec_xfa, Conc_.X_fa);
    // RRPtrs_[15] = calcRho(para_.RC.dec_xc4, Conc_.X_c4);
    // RRPtrs_[16] = calcRho(para_.RC.dec_xpro, Conc_.X_pro);
    // RRPtrs_[17] = calcRho(para_.RC.dec_xac, Conc_.X_ac);
    // RRPtrs_[18] = calcRho(para_.RC.dec_xh2, Conc_.X_h2);

}

// ************************************************************************* //