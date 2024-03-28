/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

//- Inhibition calculations

volScalarField::Internal Foam::ADMno1::calcInhibition
(
    volScalarField Y, 
    dimensionedScalar denum
)
{
    // return 1.0/ (1.0 + (Y.internalField() / denum));
    return 1.0/ (1.0 + (Y / denum));
}

volScalarField::Internal Foam::ADMno1::dCalcInhibition
(
    volScalarField Y, 
    dimensionedScalar denom
)
{
    return (-1 / (denom * (1 + (Y/denom)) * (1 + (Y/denom))));
}

volScalarField::Internal Foam::ADMno1::calcInhibitionHP
(
    volScalarField::Internal Shp,
    const dimensionedScalar UL,
    const dimensionedScalar LL,
    const dimensionedScalar n
)
{
    dimensionedScalar Kph = pow(10, -0.5 * (UL + LL));
    Kph.dimensions().reset(Shp.dimensions());

    return pow(Kph, n) / (pow(Shp, n) + pow(Kph, n));
}


//- Kinetic rate calculations

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k,
    volScalarField X
)
{
    return k * X.internalField();
}

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k, 
    volScalarField S,
    const dimensionedScalar K,
    volScalarField X,
    volScalarField::Internal I
)
{
    return k * (S.internalField() / (K + S.internalField())) * X.internalField() * I;
}

volScalarField::Internal Foam::ADMno1::calcRho
(
    const dimensionedScalar k, 
    volScalarField S1,
    const dimensionedScalar K,
    volScalarField X,
    volScalarField S2,
    volScalarField::Internal I
)
{
    return k * (S1.internalField() / (K + S1.internalField())) * X.internalField() * 
           (1.0 / (1.0 + (S2.internalField() / S1.internalField()))) * I;
}


//- Components source term calculations

volScalarField::Internal Foam::ADMno1::concPerComponent
(
    label j,
    const admPara para,
    PtrList<volScalarField::Internal> KRPtrs
)
{
    volScalarField::Internal dY
    (
        IOobject
        (
            "dY",
            KRPtrs[0].mesh().time().timeName(),
            KRPtrs[0].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        KRPtrs[0].mesh(),
        dimensionedScalar
        (
           "dY_Default", 
            KRPtrs[0].dimensions(),
            Zero
        )
    );

    // TODO: room for optimization (dont loop through all 19 elements since a lot of them are 0s)
    for (int i = 0; i < 19; i++) {
        dY += para.STOI[i][j] * KRPtrs[i] * para.DTOS(); //check if it works
    }
    return dY;
}


//- Sh2 calculations

volScalarField::Internal Foam::ADMno1::fSh2
(
    const surfaceScalarField &flux,
    volScalarField Sh2Temp
)
{
	volScalarField::Internal I_h2fa = calcInhibition // h2_fa
    (
        Sh2Temp,
        para_.KI().h2fa
    );

	volScalarField::Internal I_h2c4 = calcInhibition // h2_c4
    (
        Sh2Temp,
        para_.KI().h2c4
    );

	volScalarField::Internal I_h2pro = calcInhibition // h2_pro
    (
        Sh2Temp,
        para_.KI().h2pro
    );

    PtrList<volScalarField::Internal> KRPtrs_temp = KRPtrs_;

    KRPtrs_temp[6] = KRPtrs_temp[6]/IPtrs_[4]*I_h2fa;
    KRPtrs_temp[7] = KRPtrs_temp[7]/IPtrs_[5]*I_h2c4;
    KRPtrs_temp[8] = KRPtrs_temp[8]/IPtrs_[5]*I_h2c4;
    KRPtrs_temp[9] = KRPtrs_temp[9]/IPtrs_[6]*I_h2pro;

    KRPtrs_temp[11] = calcRho
    (
        para_.kDec().m_h2,
        Sh2Temp,
        para_.KS().h2,
        YPtrs_[22], // Xh2
        IPtrs_[2] * IPtrs_[3] // Iphh2*IIN
    );

    return concPerComponent(7, para_, KRPtrs_temp); // + convection - fGasRhoH2(paraPtr, Sh2);

    // return I_h2pro;
    
    // rho_temp[6] = calcRho(paraPtr->RC.m_fa, ConcTemp_.S_fa, paraPtr->K_S.fa, ConcTemp_.X_fa, I.phaa * I.IN * I_h2fa_temp);
    // rho_temp[7] = calcRho(paraPtr->RC.m_c4, ConcTemp_.S_va, paraPtr->K_S.c4, ConcTemp_.X_c4, ConcTemp_.S_bu, I.phaa * I.IN * I_h2c4_temp);
    // rho_temp[8] = calcRho(paraPtr->RC.m_c4, ConcTemp_.S_bu, paraPtr->K_S.c4, ConcTemp_.X_c4, ConcTemp_.S_va, I.phaa * I.IN * I_h2c4_temp);
    // rho_temp[9] = calcRho(paraPtr->RC.m_pro, ConcTemp_.S_pro, paraPtr->K_S.pro, ConcTemp_.X_pro, I.phaa * I.IN * I_h2pro_temp);
    // rho_temp[11] = calcRho(paraPtr->RC.m_h2, Sh2, paraPtr->K_S.h2, ConcTemp_.X_h2, I.phh2 * I.IN);

    // int j = 7; // Sh2 as the Conc_[7] to be used in STOI matrix
    // return concPerComponentNR(paraPtr, &j, rho_temp) + qLocal*(CInflow_.S_h2 - Sh2) - fGasRhoH2(paraPtr, Sh2);
}

volScalarField::Internal Foam::ADMno1::dfSh2
(
    const surfaceScalarField &flux,
    volScalarField Sh2Temp
)
{
    volScalarField::Internal I_h2fa = dCalcInhibition // h2_fa
    (
        Sh2Temp,
        para_.KI().h2fa
    );

	volScalarField::Internal I_h2c4 = dCalcInhibition // h2_c4
    (
        Sh2Temp,
        para_.KI().h2c4
    );

	volScalarField::Internal I_h2pro = dCalcInhibition // h2_pro
    (
        Sh2Temp,
        para_.KI().h2pro
    );

    PtrList<volScalarField::Internal> KRPtrs_temp = KRPtrs_;

    KRPtrs_temp[6] = KRPtrs_temp[6]/IPtrs_[4]*I_h2fa;
    KRPtrs_temp[7] = KRPtrs_temp[7]/IPtrs_[5]*I_h2c4;
    KRPtrs_temp[8] = KRPtrs_temp[8]/IPtrs_[5]*I_h2c4;
    KRPtrs_temp[9] = KRPtrs_temp[9]/IPtrs_[6]*I_h2pro;

    // KRPtrs_temp[11] = calcRho
    // (
    //     para_.kDec().m_h2,
    //     Sh2Temp,
    //     para_.KS().h2,
    //     YPtrs_[22], // Xh2
    //     IPtrs_[2] * IPtrs_[3] // Iphh2*IIN
    // );

    return I_h2pro; // + convection

    // rho_temp[6] = calcRho(paraPtr->RC.m_fa, ConcTemp_.S_fa, paraPtr->K_S.fa, ConcTemp_.X_fa, I.phaa * I.IN * I_h2fa_temp);
    // rho_temp[7] = calcRho(paraPtr->RC.m_c4, ConcTemp_.S_va, paraPtr->K_S.c4, ConcTemp_.X_c4, ConcTemp_.S_bu, I.phaa * I.IN * I_h2c4_temp);
    // rho_temp[8] = calcRho(paraPtr->RC.m_c4, ConcTemp_.S_bu, paraPtr->K_S.c4, ConcTemp_.X_c4, ConcTemp_.S_va, I.phaa * I.IN * I_h2c4_temp);
    // rho_temp[9] = calcRho(paraPtr->RC.m_pro, ConcTemp_.S_pro, paraPtr->K_S.pro, ConcTemp_.X_pro, I.phaa * I.IN * I_h2pro_temp);
    // data_type temp = paraPtr->K_S.h2 + Sh2;
    // rho_temp[11] = (paraPtr->RC.m_h2 * ConcTemp_.X_h2 * I.phh2 * I.IN * paraPtr->K_S.h2) / (temp * temp);
    
    // int j = 7; // Sh2 as the Conc_[7] to be used in STOI matrix
    // return concPerComponentNR(paraPtr, &j, rho_temp) - qLocal - dfGasRhoH2(paraPtr);
}


//- Acid-base calculations

volScalarField::Internal Foam::ADMno1::fSion
(
    const dimensionedScalar Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return Kax * Sx / (Kax + Shp);
}

volScalarField::Internal Foam::ADMno1::fSion
(
    const volScalarField Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return Kax * Sx / (Kax + Shp);
}

volScalarField::Internal Foam::ADMno1::dfSion
(
    const dimensionedScalar Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return - Kax * Sx / ((Kax + Shp) * (Kax + Shp));
}

volScalarField::Internal Foam::ADMno1::dfSion
(
    const volScalarField Kax,
    const volScalarField::Internal& Sx,
    volScalarField::Internal Shp
)
{
    return - Kax.internalField() * Sx / ((Kax + Shp) * (Kax + Shp));
}


// ************************************************************************* //