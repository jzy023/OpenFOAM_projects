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

void Foam::ADMno1::kineticRate()
{

    //- Inhibiitons

    IPtrs_[0] = calcInhibitionHP // aa
    (
        ShP_,
        para_.pHL().ULaa, 
        para_.pHL().LLaa,
        nIaa_
    );

    IPtrs_[1] = calcInhibitionHP // ac
    (
        ShP_,
        para_.pHL().ULac, 
        para_.pHL().LLac,
        nIac_
    );

    IPtrs_[2] = calcInhibitionHP // h2
    (
        ShP_,
        para_.pHL().ULh2, 
        para_.pHL().LLh2,
        nIh2_
    );

    // >>> TODO: which one is correct ???
    // IPtrs_[3] = calcInhibition // IN
    // (
    //     YPtrs_[10], // SIN
    //     para_.KS().IN
    // );

    IPtrs_[3] = 1.0 / (1.0 + (para_.KS().IN / YPtrs_[10]));
    

	IPtrs_[4] = calcInhibition // h2fa
    (
        YPtrs_[7], // Sh2
        para_.KI().h2fa
    );

	IPtrs_[5] = calcInhibition // h2c4
    (
        YPtrs_[7], // Sh2
        para_.KI().h2c4
    );

	IPtrs_[6] = calcInhibition // h2pro
    (
        YPtrs_[7], // Sh2
        para_.KI().h2pro
    );

	IPtrs_[7] = calcInhibition // nh3
    (
        MPtrs_[1], // Snh3
        para_.KI().nh3
    );

    //- Kinetic rates

    KRPtrs_[0] = calcRho
    (
        para_.kDec().dis,
        YPtrs_[12] // Xc
    );

    KRPtrs_[1] = calcRho
    (
        para_.kDec().hyd_ch,
        YPtrs_[13] // Xch
    );

    KRPtrs_[2] = calcRho
    (
        para_.kDec().hyd_pr,
        YPtrs_[14] // Xpr
    );

    KRPtrs_[3] = calcRho
    (
        para_.kDec().hyd_li,
        YPtrs_[15] // Xli
    );

    KRPtrs_[4] = calcRho
    (
        para_.kDec().m_su,
        YPtrs_[0], // Ssu
        para_.KS().su,
        YPtrs_[16], // Xsu
        IPtrs_[0] * IPtrs_[3] // Iphaa*IIN
    );

    KRPtrs_[5] = calcRho
    (
        para_.kDec().m_aa,
        YPtrs_[1], // Saa
        para_.KS().aa,
        YPtrs_[17], // Xaa
        IPtrs_[0] * IPtrs_[3] // Iphaa*IIN
    );

    KRPtrs_[6] = calcRho
    (
        para_.kDec().m_fa,
        YPtrs_[2], // Sfa
        para_.KS().fa,
        YPtrs_[18], // Xfa
        IPtrs_[0] * IPtrs_[3] * IPtrs_[4] //Iphaa*IIN*Ih2fa
    );

    KRPtrs_[7] = calcRho
    (
        para_.kDec().m_c4,
        YPtrs_[3], // Sva
        para_.KS().c4,
        YPtrs_[19], // Xc4
        YPtrs_[4],  // Sbu
        IPtrs_[0] * IPtrs_[3] * IPtrs_[5] //Iphaa*IIN*Ih2c4
    );

    KRPtrs_[8] = calcRho
    (
        para_.kDec().m_c4,
        YPtrs_[4], // Sbu
        para_.KS().c4,
        YPtrs_[19], // Xc4
        YPtrs_[3],  // Sva
        IPtrs_[0] * IPtrs_[3] * IPtrs_[5] //Iphaa*IIN*Ih2c4
    );

    KRPtrs_[9] = calcRho
    (
        para_.kDec().m_pro,
        YPtrs_[5], // Spro
        para_.KS().pro,
        YPtrs_[20], // Xpro
        IPtrs_[0] * IPtrs_[3] * IPtrs_[6]  //Iphaa*IIN*Ih2pro
    );

    KRPtrs_[10] = calcRho
    (
        para_.kDec().m_ac,
        YPtrs_[6], // Sac
        para_.KS().ac,
        YPtrs_[21], // Xac
        IPtrs_[1] * IPtrs_[3] * IPtrs_[7] // Iphac*IIN*Inh3
    );

	// >>> in Rosen et al. implementation, no intermediate used for S_h2
	KRPtrs_[11] = calcRho
    (
        para_.kDec().m_h2,
        YPtrs_[7], // Sh2
        para_.KS().h2,
        YPtrs_[22], // Xh2
        IPtrs_[2] * IPtrs_[3] // Iphh2*IIN
    );

	KRPtrs_[12] = calcRho
    (
        para_.kDec().dec_xsu,
        YPtrs_[16] // Xsu
    );

    KRPtrs_[13] = calcRho
    (
        para_.kDec().dec_xaa,
        YPtrs_[17] // Xaa
    );

    KRPtrs_[14] = calcRho
    (
        para_.kDec().dec_xfa,
        YPtrs_[18] // Xfa
    );

    KRPtrs_[15] = calcRho
    (
        para_.kDec().dec_xc4,
        YPtrs_[19] // Xc4
    );

    KRPtrs_[16] = calcRho
    (
        para_.kDec().dec_xpro,
        YPtrs_[20] // Xpro
    );

    KRPtrs_[17] = calcRho
    (
        para_.kDec().dec_xac,
        YPtrs_[21] // Xac
    );

    KRPtrs_[18] = calcRho
    (
        para_.kDec().dec_xh2,
        YPtrs_[22] // Xh2
    );
}


void Foam::ADMno1::sulSourceRate()
{
    for(label j = 0; j < 7; j++)
    {
        for (int i = 0; i < 19; i++)
        {
            dYPtrs_[j] += para_.STOI[i][j] * KRPtrs_[i] * para_.DTOS(); //check if it works
        }
    }

    for(label j = 8; j < YPtrs_.size(); j++)
    {
        for (int i = 0; i < 19; i++)
        {
            dYPtrs_[j] += para_.STOI[i][j] * KRPtrs_[i] * para_.DTOS();
        }
    }
    
    //- calculate dSh2 iteratively
    // TODO: implement it! with Rosen et al.
    // RSh2(); 

    //- calculate with STOI and gas transer
    // dYPtrs_[7] -= GRPtrs_[0]; // Sh2 - Gh2

    dYPtrs_[8] -= GRPtrs_[1]; // Sch4 - Gch4

    // Sco2!
    dYPtrs_[9] -= GRPtrs_[2]; // SIC - Gco2

}

// ************************************************************************* //