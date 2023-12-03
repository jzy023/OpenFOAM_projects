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
    ADM::admSoluables

Description
    H file for defining inline functions of AMDno1 class

\*---------------------------------------------------------------------------*/
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
template<class Cmpt>
inline Foam::admSoluables<Cmpt>::admSoluables(const Foam::zero)
:
    admSoluables::vsType(Zero)
{}
 
 
template<class Cmpt>
template<class Cmpt2>
inline Foam::admSoluables<Cmpt>::admSoluables
(
    const VectorSpace<admSoluables<Cmpt2>, Cmpt2, 12>& sol
)
:
    admSoluables::vsType(sol)
{}
 
 
template<class Cmpt>
inline Foam::admSoluables<Cmpt>::admSoluables
(
    const Cmpt& Ssu, const Cmpt& Saa, const Cmpt& Sfa,
    const Cmpt& Sva, const Cmpt& Sbu, const Cmpt& Spro,
    const Cmpt& Sac, const Cmpt& Sh2, const Cmpt& Sch4,
    const Cmpt& SIC, const Cmpt& SIN, const Cmpt& SI
)
{
    this->v_[su_]    = Ssu;
    this->v_[aa_]    = Saa;
    this->v_[fa_]    = Sfa;
    this->v_[va_]    = Sva;
    this->v_[bu_]    = Sbu;
    this->v_[pro_]   = Spro;
    this->v_[ac_]    = Sac;
    this->v_[h2_]    = Sh2;
    this->v_[ch4_]   = Sch4;
    this->v_[IC_]    = SIC;
    this->v_[IN_]    = SIN;
    this->v_[I_]     = SI;
}


template<class Cmpt>
inline Foam::admSoluables<Cmpt>::admSoluables(Istream& is)
:
    admSoluables::vsType(is)
{}

 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::su() const
{
    return this->v_[su_];
}
 
template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::aa() const
{
    return this->v_[aa_];
}
 
template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::fa() const
{
    return this->v_[fa_];
}

template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::va() const
{
    return this->v_[va_];
}
 
template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::bu() const
{
    return this->v_[bu_];
}
 
template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::pro() const
{
    return this->v_[pro_];
}

template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::ac() const
{
    return this->v_[ac_];
}
 
template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::h2() const
{
    return this->v_[h2_];
}
 
template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::ch4() const
{
    return this->v_[ch4_];
}

template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::IC() const
{
    return this->v_[IC_];
}
 
template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::IN() const
{
    return this->v_[IN_];
}
 
template<class Cmpt>
inline const Cmpt& Foam::admSoluables<Cmpt>::I() const
{
    return this->v_[I_];
}

template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::su()
{
    return this->v_[su_];
}
 
template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::aa()
{
    return this->v_[aa_];
}
 
template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::fa()
{
    return this->v_[fa_];
}

template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::va()
{
    return this->v_[va_];
}
 
template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::bu()
{
    return this->v_[bu_];
}
 
template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::pro()
{
    return this->v_[pro_];
}

template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::ac()
{
    return this->v_[ac_];
}
 
template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::h2()
{
    return this->v_[h2_];
}
 
template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::ch4()
{
    return this->v_[ch4_];
}

template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::IC()
{
    return this->v_[IC_];
}
 
template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::IN()
{
    return this->v_[IN_];
}
 
template<class Cmpt>
inline Cmpt& Foam::admSoluables<Cmpt>::I()
{
    return this->v_[I_];
}

// ************************************************************************* //