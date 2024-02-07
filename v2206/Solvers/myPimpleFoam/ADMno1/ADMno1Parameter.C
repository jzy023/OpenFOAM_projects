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

admPara::admPara
(
    word runMode
)
:
    errMessage
    (
        printErrMessage(runMode)
    ),
    Tbase_
    (
        defineTbase(runMode)
    ),
    Top_
    (
        defineTop(runMode)
    ),
    kDec_
    (
        defineRC(runMode)
    ),
    yB_
    (
        defineYields(runMode)
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
        defineKI(runMode)
    ),
    KS_
    (
        defineKS(runMode)
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
    defineInitialState(runMode);
	defineSTOI();
	// defineAcidBaseDAE();
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int admPara::printErrMessage
(
    word runMode
)
{
    if (!namesOpMode.found(runMode)) 
    {
        std::cerr << nl << "--> FOAM FATAL IO ERROR:" << nl
                  << "Unknown operation mode type: " << runMode
                  << "\n\nValid operation mode types:\n";
        std::cerr << namesOpMode.size() << "(\n";
        forAll(namesOpMode, i)
        {
            std::cerr << namesOpMode[i] << nl;
        }
        std::cerr << ")\n";
        std::exit(1);
    }

    return 0;
}

dimensionedScalar admPara::defineTbase
(
    word runMode
)
{
    if(runMode == "Thermo")
    { // TODO: need check
        return dimensionedScalar(dimTemperature, 55);
    }
    else
    {
        return dimensionedScalar(dimTemperature, 25);
    }
}

dimensionedScalar admPara::defineTop
(
    word runMode
)
{
    if(runMode == "Thermo")
    {
        return dimensionedScalar(dimTemperature, 55);
    }
    else
    {
        return dimensionedScalar(dimTemperature, 35);
    }
}

decayRate admPara::defineRC
(
    word runMode
)
{
    if(runMode == "Meso")
    {
        return decayRate
        (
            0.4,    // dis
            0.25,   // hyd_ch
            0.2,    // hyd_pr
            0.1,    // hyd_li
            30.0,   // m_su
            50.0,   // m_aa
            6.0,    // m_fa
            20.0,   // m_c4
            13.0,   // m_pro
            8.0,    // m_ac
            35.0,   // m_h2
            0.02,   // dec_xsu
            0.02,   // dec_xaa
            0.02,   // dec_xfa
            0.02,   // dec_xc4
            0.02,   // dec_xpro
            0.02,   // dec_xac
            0.02    // dec_xh2
        );
    }
    else if(runMode == "MesoSolid")
    {
        return decayRate
        (
            0.5,    // dis
            10.0,   // hyd_ch
            10.0,   // hyd_pr
            10.0,   // hyd_li
            30.0,   // m_su
            50.0,   // m_aa
            6.0,    // m_fa
            20.0,   // m_c4
            13.0,   // m_pro
            8.0,    // m_ac
            35.0,   // m_h2
            0.02,   // dec_xsu
            0.02,   // dec_xaa
            0.02,   // dec_xfa
            0.02,   // dec_xc4
            0.02,   // dec_xpro
            0.02,   // dec_xac
            0.02    // dec_xh2
        );
    }
    else
    {
        return decayRate
        (
            1.0,    // dis
            10.0,   // hyd_ch
            10.0,   // hyd_pr
            10.0,   // hyd_li
            70.0,   // m_su
            70.0,   // m_aa
            10.0,   // m_fa
            30.0,   // m_c4
            20.0,   // m_pro
            16.0,   // m_ac
            35.0,   // m_h2
            0.04,   // dec_xsu
            0.04,   // dec_xaa
            0.04,   // dec_xfa
            0.04,   // dec_xc4
            0.04,   // dec_xpro
            0.04,   // dec_xac
            0.04    // dec_xh2
        );
    }
}

yieldBiomass admPara::defineYields
(
    word runMode
)
{
    if(runMode == "Meso")
    {
        return yieldBiomass
        (
            0.10,
            0.08,
            0.06,
            0.06,
            0.04,
            0.05,
            0.06
        );
    }
    else
    {
        return yieldBiomass
        (
            0.10,
            0.08,
            0.06,
            0.06,
            0.05,
            0.05,
            0.06
        );
    }
}

inhibitionParaI admPara::defineKI
(
    word runMode
)
{
    if(runMode == "Meso")
    {
        return inhibitionParaI
        (
            5.0e-6, // h2fa
	        1.0e-5, // h2c4
	        3.5e-6, // h2pro
	        0.0018  // nh3
        );
    }
    else if(runMode == "MesoSolid")
    {   // TODO: check value!
        return inhibitionParaI
        (
            5.0e-6, // h2fa
	        1.0e-5, // h2c4
	        3.5e-6, // h2pro
	        0.0018  // nh3
        );
    }
    else
    {
        return inhibitionParaI
        (   // from Rosen ADM-BSM DIGESTERPAR[55]
            5.0e-6, // h2fa
	        3.0e-5, // h2c4
	        1e-5,   // h2pro
	        0.011   // nh3
        );
    }
}

inhibitionParaS admPara::defineKS
(
    word runMode
)
{
    if(runMode == "Meso")
    {
        return inhibitionParaS
        (
            1e-4,   // IN
            1e-4,   // nh3
            0.5,    // su
            0.3,    // aa
            0.4,    // fa
            0.3,    // c4
            0.3,    // pro
            0.15,   // ac
            2.5e-5  // h2
        );
    }
    else if(runMode == "MesoSolid")
    {   // TODO: check value!
        return inhibitionParaS
        (
            1e-4,   // IN
            1e-4,   // nh3
            0.5,    // su
            0.3,    // aa
            0.4,    // fa
            0.2,    // c4
            0.1,    // pro
            0.15,   // ac
            7e-6    // h2
        );
    }
    else
    {
        return inhibitionParaS
        (
            1e-4,   // IN
            1e-4,   // nh3
            1.0,    // su
            0.3,    // aa
            0.4,    // fa
            0.4,    // c4
            0.3,    // pro
            0.3,    // ac
            5e-5    // h2
        );
    }
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
	STOI[4][10] = -yB_.su * NC_.bac;

	STOI[5][1] = -1;
	STOI[5][3] = (1 - yB_.aa) * yP_.va_aa;
	STOI[5][4] = (1 - yB_.aa) * yP_.bu_aa;
	STOI[5][5] = (1 - yB_.aa) * yP_.pro_aa;
	STOI[5][6] = (1 - yB_.aa) * yP_.ac_aa;
	STOI[5][7] = (1 - yB_.aa) * yP_.h2_aa;
	STOI[5][10] = NC_.aa - yB_.aa * NC_.bac;

	STOI[6][2] = -1;
	STOI[6][6] = (1 - yB_.fa) * 0.7;
	STOI[6][7] = (1 - yB_.fa) * 0.3;
	STOI[6][10] = -yB_.fa * NC_.bac;

	STOI[7][3] = -1;
	STOI[7][5] = (1 - yB_.c4) * 0.54;
	STOI[7][6] = (1 - yB_.c4) * 0.31;
	STOI[7][7] = (1 - yB_.c4) * 0.15;
	STOI[7][10] = -yB_.c4 * NC_.bac;

	STOI[8][4] = -1;
	STOI[8][6] = (1 - yB_.c4) * 0.8;
	STOI[8][7] = (1 - yB_.c4) * 0.2;
	STOI[8][10] = -yB_.c4 * NC_.bac;

	STOI[9][5] = -1;
	STOI[9][6] = (1 - yB_.pro) * 0.57;
	STOI[9][7] = (1 - yB_.pro) * 0.43;
	STOI[9][10] = -yB_.pro * NC_.bac;

	STOI[10][6] = -1;
	STOI[10][8] = (1 - yB_.ac);
	STOI[10][10] = -yB_.ac * NC_.bac;

	STOI[11][7] = -1;
	STOI[11][8] = (1 - yB_.h2);
	STOI[11][10] = -yB_.h2 * NC_.bac;

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

// void admPara::defineAcidBaseDAE()
// {
// 	abDAE.resize(9);
// 	for (int i = 0; i < 9; i++)
// 	{
// 		abDAE[i].resize(9);
// 	}

// 	for (int i = 0; i < 8; i++)
// 	{
// 		abDAE[i][i] = 1;
// 	}

// 	abDAE[7][6] = 1;

// 	abDAE[8][0] = -1;		//OH-
// 	abDAE[8][1] = static_cast<scalar> (-1 / 208); //vaN
// 	abDAE[8][2] = static_cast<scalar> (-1 / 160); //buN
// 	abDAE[8][3] = static_cast<scalar> (-1 / 112); //proN
// 	abDAE[8][4] = static_cast<scalar> (-1 / 64);	//acN
// 	abDAE[8][5] = -1;		//hco3N
// 	abDAE[8][6] = 0;		//nh3
// 	abDAE[8][7] = 1;		//nh+
// 	abDAE[8][8] = 1;
// };

void admPara::listResizing()
{
    Gini_.resize(3);
    Eini_.resize(7);
    Mini_.resize(3); 
}

void admPara::defineInitialState(word runMode)
{
    listResizing();

    if(runMode == "Meso")
    {// TODO: update values
        Gini_[0] = 1.1032e-5;    // G_h2
        Gini_[1] = 1.6535;       // G_ch4
        Gini_[2] = 0.0135;       // G_co2

        Eini_[0] = 0.012284;     // S_vaN
        Eini_[1] = 0.013953;     // S_buN
        Eini_[2] = 0.017511;     // S_proN
        Eini_[3] = 0.089035;     // S_acN
        Eini_[4] = 0.08568;      // S_hco3N
        Eini_[5] = 0;            // S_ohN
        Eini_[6] = 5.4562e-8;    // S_hP

        Mini_[0] = 0.00942;      // S_co2
        Mini_[1] = 0.001884;     // S_nh3
        Mini_[2] = 0.092584;     // S_nh4
    }
    else if(runMode == "MesoSolid")
    {
        Gini_[0] = 1.1032e-5;    // G_h2
        Gini_[1] = 1.6535;       // G_ch4
        Gini_[2] = 0.0135;       // G_co2

        Eini_[0] = 0.012284;     // S_vaN
        Eini_[1] = 0.013953;     // S_buN
        Eini_[2] = 0.017511;     // S_proN
        Eini_[3] = 0.089035;     // S_acN
        Eini_[4] = 0.08568;      // S_hco3N
        Eini_[5] = 0;            // S_ohN
        Eini_[6] = 5.4562e-8;    // S_hP

        Mini_[0] = 0.00942;      // S_co2
        Mini_[1] = 0.001884;     // S_nh3
        Mini_[2] = 0.092584;     // S_nh4
    }
    else
    {// TODO: update values
        Gini_[0] = 1.1032e-5;    // G_h2
        Gini_[1] = 1.6535;       // G_ch4
        Gini_[2] = 0.0135;       // G_co2

        Eini_[0] = 0.012284;     // S_vaN
        Eini_[1] = 0.013953;     // S_buN
        Eini_[2] = 0.017511;     // S_proN
        Eini_[3] = 0.089035;     // S_acN
        Eini_[4] = 0.08568;      // S_hco3N
        Eini_[5] = 0;            // S_ohN
        Eini_[6] = 5.4562e-8;    // S_hP

        Mini_[0] = 0.00942;      // S_co2
        Mini_[1] = 0.001884;     // S_nh3
        Mini_[2] = 0.092584;     // S_nh4
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //