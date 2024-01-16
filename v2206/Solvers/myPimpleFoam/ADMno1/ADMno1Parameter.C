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
    C file for defining AMDno1 parameter class functions

\*---------------------------------------------------------------------------*/

#include "ADMno1Parameter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

admPara::admPara(int admMode) 
: 
    admMode_(admMode),
    // Yields of product
    yP_({0.10, 0.20, 0.20,
         0.20, 0.30, 0.95,
         0.19, 0.13, 0.27,
         0.41, 0.06, 0.23,
         0.26, 0.05, 0.40}),
    // si_xc;	xi_sc;	ch_xc;
    // pr_xc;	li_xc;	fa_li;
    // h2_su;	bu_su;	pro_su;
    // ac_su;	h2_aa;	va_aa;
    // bu_aa;	pro_aa;	ac_aa;

    // Carbon Content
    CC_({0.02786, 0.03, 0.0313,
         0.03, 0.022, 0.03,
         0.0313, 0.03, 0.0217,
         0.025, 0.0268, 0.0313,
         0.0313, 0.024, 0.0156}),
    // xc;		si;		ch;
    // pr;		li;		xi;
    // su;		aa;		fa;
    // bu;		pro;	ac;
    // bac;		va;		ch4;

    // Nitrogen Content
    N_aa_(0.007), N_bac_(0.005714),

    // Acid Base Kinetics
    kAB({1e10, 1e10, 1e10,
         1e10, 1e10, 1e10}),

    // Acide base Equilibrium Para
    Ka({1.380e-5, 1.514e-5, 1.318e-5,
        1.738e-5, 4.467e-7, 5.623e-10,
        1e-14}),
    // va;	bu;	  pro;
    // ac;	co2;  IN;
    // Kw

    // Henry's Law Coefficients
    KH({7.384654293536963e-04, // h2
        0.001161902733673,     // ch4
        0.027146692900075,     // co2
        0.031300000000000}),   // h2o

    // Gas Transfer Coefficients
    // TODO: simplify
    // kL({200, 200, 200}),
    kLa(200),

    pH_UL_aa(5.5), pH_LL_aa(4),
    pH_UL_ac(7), pH_LL_ac(6),
    pH_UL_h2(6), pH_LL_h2(5)
{
	if (admMode == 0)
	{
		defineRCMeso();
		defineYieldsMeso();
		defineKIMeso();
		defineKSMeso();
	}
	else if (admMode == 1)
	{
		defineRCMesoSolid();
		defineYieldsMeso();
		defineKIMeso();
		defineKSMesoSolid();
	}
	else if (admMode == 2)
	{
		defineRCThermo();
		defineYieldsThermo();
		defineKIThermo();
		defineKSThermo();
	}
	else
	{
		throw "admMode undefined";
	}

	/* ========================= Stoichematic for ==========================
	>>> calculate Stoichematic table based on yields 
	======================================================================== */
	defineSTOI();
	defineAcidBaseDAE();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void admPara::defineRCMeso()
{
	RC.dis = 0.4;
	RC.hyd_ch = 0.25;
	RC.hyd_pr = 0.2;
	RC.hyd_li = 0.1;

	RC.m_su = 30.0;
	RC.m_aa = 50.0;
	RC.m_fa = 6.0;
	RC.m_c4 = 20.0;
	RC.m_pro = 13.0;
	RC.m_ac = 8.0;
	RC.m_h2 = 35.0;

	RC.dec_xsu = 0.02;
	RC.dec_xaa = 0.02;
	RC.dec_xfa = 0.02;
	RC.dec_xc4 = 0.02;
	RC.dec_xpro = 0.02;
	RC.dec_xac = 0.02;
	RC.dec_xh2 = 0.02;
}

void admPara::defineRCMesoSolid()
{
	RC.dis = 0.5;
	RC.hyd_ch = 10.0;
	RC.hyd_pr = 10.0;
	RC.hyd_li = 10.0;

	RC.m_su = 30.0;
	RC.m_aa = 50.0;
	RC.m_fa = 6.0;
	RC.m_c4 = 20.0;
	RC.m_pro = 13.0;
	RC.m_ac = 8.0;
	RC.m_h2 = 35.0;

	RC.dec_xsu = 0.02;
	RC.dec_xaa = 0.02;
	RC.dec_xfa = 0.02;
	RC.dec_xc4 = 0.02;
	RC.dec_xpro = 0.02;
	RC.dec_xac = 0.02;
	RC.dec_xh2 = 0.02;
}

void admPara::defineRCThermo()
{
	RC.dis = 1.0;
	RC.hyd_ch = 10.0;
	RC.hyd_pr = 10.0;
	RC.hyd_li = 10.0;

	RC.m_su = 70.0;
	RC.m_aa = 70.0;
	RC.m_fa = 10.0;
	RC.m_c4 = 30.0;
	RC.m_pro = 20.0;
	RC.m_ac = 16.0;
	RC.m_h2 = 35.0;

	RC.dec_xsu = 0.04;
	RC.dec_xaa = 0.04;
	RC.dec_xfa = 0.04;
	RC.dec_xc4 = 0.04;
	RC.dec_xpro = 0.04;
	RC.dec_xac = 0.04;
	RC.dec_xh2 = 0.04;
}

void admPara::defineYieldsMeso()
{
	yB_.su = 0.10;
	yB_.aa = 0.08;
	yB_.fa = 0.06;
	yB_.c4 = 0.06;
	yB_.pro = 0.04;
	yB_.ac = 0.05;
	yB_.h2 = 0.06;
}

void admPara::defineYieldsThermo()
{
	yB_.su = 0.10;
	yB_.aa = 0.08;
	yB_.fa = 0.06;
	yB_.c4 = 0.06;
	yB_.pro = 0.05;
	yB_.ac = 0.05;
	yB_.h2 = 0.06;
}

void admPara::defineKIMeso()
{
	K_I.h2fa = 5.0e-6;
	K_I.h2c4 = 1.0e-5;
	K_I.h2pro = 3.5e-6;
	K_I.nh3 = 0.0018;
}

void admPara::defineKIThermo()
{
	// K_I.h2fa = N/A !!! from ADM no1 book
	K_I.h2fa = 5.0e-6;  // from Rosen ADM-BSM DIGESTERPAR[55]
	K_I.h2c4 = 3.0e-5;
	K_I.h2pro = 1e-5;
	K_I.nh3 = 0.011;
}

void admPara::defineKSMeso()
{
	K_S.IN = 1e-4;
	K_S.nh3 = 1e-4;
	K_S.su = 0.5;
	K_S.aa = 0.3;
	K_S.fa = 0.4;
	K_S.c4 = 0.3;
	K_S.pro = 0.3;
	K_S.ac = 0.15;
	K_S.h2 = 2.5e-5;
}

void admPara::defineKSMesoSolid()
{
	K_S.IN = 1e-4;
	K_S.nh3 = 1e-4;
	K_S.su = 0.5;
	K_S.aa = 0.3;
	K_S.fa = 0.4;
	K_S.c4 = 0.2;
	K_S.pro = 0.1;
	K_S.ac = 0.15;
	K_S.h2 = 7e-6;
}

void admPara::defineKSThermo()
{
	K_S.IN = 1e-4;
	K_S.nh3 = 1e-4;
	K_S.su = 1.0;
	K_S.aa = 0.3;
	K_S.fa = 0.4;
	K_S.c4 = 0.4;
	K_S.pro = 0.3;
	K_S.ac = 0.3;
	K_S.h2 = 5e-5;
}

void admPara::defineSTOI()
{
	STOI.resize(19);
	for (int i = 0; i < 19; i++)
	{
		STOI[i].resize(24);
	}
	// rows for processes and columns for substrates (visco[j,i] from the table)
	STOI[0][11] = yP_.si_xc;
	STOI[1][0] = 1;
	STOI[2][1] = 1;
	STOI[3][0] = 1 - yP_.fa_li;
	STOI[3][2] = yP_.fa_li;

	STOI[4][0] = -1;
	STOI[4][4] = (1 - yB_.su) * yP_.bu_su;
	STOI[4][5] = (1 - yB_.su) * yP_.pro_su;
	STOI[4][6] = (1 - yB_.su) * yP_.ac_su;
	STOI[4][7] = (1 - yB_.su) * yP_.h2_su;
	STOI[4][10] = -yB_.su * N_bac_;

	STOI[5][1] = -1;
	STOI[5][3] = (1 - yB_.aa) * yP_.va_aa;
	STOI[5][4] = (1 - yB_.aa) * yP_.bu_aa;
	STOI[5][5] = (1 - yB_.aa) * yP_.pro_aa;
	STOI[5][6] = (1 - yB_.aa) * yP_.ac_aa;
	STOI[5][7] = (1 - yB_.aa) * yP_.h2_aa;
	STOI[5][10] = N_aa_ - yB_.aa * N_bac_;

	STOI[6][2] = -1;
	STOI[6][6] = (1 - yB_.fa) * 0.7;
	STOI[6][7] = (1 - yB_.fa) * 0.3;
	STOI[6][10] = -yB_.fa * N_bac_;

	STOI[7][3] = -1;
	STOI[7][5] = (1 - yB_.c4) * 0.54;
	STOI[7][6] = (1 - yB_.c4) * 0.31;
	STOI[7][7] = (1 - yB_.c4) * 0.15;
	STOI[7][10] = -yB_.c4 * N_bac_;

	STOI[8][4] = -1;
	STOI[8][6] = (1 - yB_.c4) * 0.8;
	STOI[8][7] = (1 - yB_.c4) * 0.2;
	STOI[8][10] = -yB_.c4 * N_bac_;

	STOI[9][5] = -1;
	STOI[9][6] = (1 - yB_.pro) * 0.57;
	STOI[9][7] = (1 - yB_.pro) * 0.43;
	STOI[9][10] = -yB_.pro * N_bac_;

	STOI[10][6] = -1;
	STOI[10][8] = (1 - yB_.ac);
	STOI[10][10] = -yB_.ac * N_bac_;

	STOI[11][7] = -1;
	STOI[11][8] = (1 - yB_.h2);
	STOI[11][10] = -yB_.h2 * N_bac_;

	STOI[0][12] = -1;
	STOI[0][13] = yP_.ch_xc;
	STOI[0][14] = yP_.pr_xc;
	STOI[0][15] = yP_.li_xc;
	STOI[0][23] = yP_.xi_xc;

	STOI[1][13] = -1;
	STOI[2][14] = -1;
	STOI[3][15] = -1;

	STOI[4][16] = yB_.su;
	STOI[5][17] = yB_.aa;
	STOI[6][18] = yB_.fa;
	STOI[7][19] = yB_.c4;

	STOI[8][19] = yB_.c4;
	STOI[9][20] = yB_.pro;
	STOI[10][21] = yB_.ac;
	STOI[11][22] = yB_.h2;

	for (int i = 12; i < 19; i++)
	{
		STOI[i][12] = 1;
		STOI[i][i + 4] = -1;
	}

	//for loop for Carbon content Using Rosen
	STOI[0][9] = -(STOI[0][11] * CC_.si + STOI[0][12] * CC_.xc + STOI[0][13] * CC_.ch +
				   STOI[0][14] * CC_.pr + STOI[0][15] * CC_.li + STOI[0][23] * CC_.xi);

	STOI[1][9] = STOI[1][0] * CC_.su +  STOI[1][13] * CC_.ch;

	STOI[2][9] = STOI[2][1] * CC_.aa + STOI[2][14] * CC_.pr;

	STOI[3][9] = STOI[3][0] * CC_.su + STOI[3][2] * CC_.fa + STOI[3][15] * CC_.li;    

	STOI[4][9] = -(STOI[4][0] * CC_.su + STOI[4][4] * CC_.bu + STOI[4][5] * CC_.pro +
				   STOI[4][6] * CC_.ac + STOI[4][16] * CC_.bac);

	STOI[5][9] = -(STOI[5][1] * CC_.aa + STOI[5][3] * CC_.va + STOI[5][4] * CC_.bu +
				   STOI[5][5] * CC_.pro + STOI[5][6] * CC_.ac + STOI[5][17] * CC_.bac);

	STOI[6][9] = STOI[6][2] * CC_.fa + STOI[6][6] * CC_.ac + STOI[6][18] * CC_.bac;

	STOI[7][9] = STOI[7][3] * CC_.va + STOI[7][5] * CC_.pro + STOI[7][6] * CC_.ac + STOI[7][19] * CC_.bac;

	STOI[8][9] = STOI[8][4] * CC_.bu + STOI[8][6] * CC_.ac + STOI[6][19] * CC_.bac;

	STOI[9][9] = -(STOI[9][5] * CC_.pro + STOI[9][6] * CC_.ac + STOI[9][20] * CC_.bac);

	STOI[10][9] = -(STOI[10][6] * CC_.ac + STOI[10][8] * CC_.ch4 + STOI[10][21] * CC_.bac);

	STOI[11][9] = -(STOI[11][8] * CC_.ch4 + STOI[11][22] * CC_.bac);

	STOI[12][9] = STOI[12][12] * CC_.xc +  STOI[12][16] * CC_.bac;
	//Extra CC for decay processes
	// STOI[13][9] = STOI[13][12] * CC_.xc + STOI[13][17]* CC_.bac;
	// STOI[14][9] = STOI[14][12] * CC_.xc + STOI[14][18]* CC_.bac;
	// STOI[15][9] = STOI[15][12] * CC_.xc + STOI[15][19]* CC_.bac;
	// STOI[16][9] = STOI[16][12] * CC_.xc + STOI[16][20]* CC_.bac;
	// STOI[17][9] = STOI[17][12] * CC_.xc + STOI[17][21]* CC_.bac;
	// STOI[18][9] = STOI[18][12] * CC_.xc + STOI[18][22]* CC_.bac;
	
}

void admPara::defineAcidBaseDAE()
{
	abDAE.resize(9);
	for (int i = 0; i < 9; i++)
	{
		abDAE[i].resize(9);
	}

	for (int i = 0; i < 8; i++)
	{
		abDAE[i][i] = 1;
	}

	abDAE[7][6] = 1;

	abDAE[8][0] = -1;		//OH-
	abDAE[8][1] = static_cast<scalar> (-1 / 208); //vaN
	abDAE[8][2] = static_cast<scalar> (-1 / 160); //buN
	abDAE[8][3] = static_cast<scalar> (-1 / 112); //proN
	abDAE[8][4] = static_cast<scalar> (-1 / 64);	//acN
	abDAE[8][5] = -1;		//hco3N
	abDAE[8][6] = 0;		//nh3
	abDAE[8][7] = 1;		//nh+
	abDAE[8][8] = 1;
}


// ************************************************************************* //