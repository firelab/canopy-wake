/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "myStochasticDispersionRAS.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myStochasticDispersionRAS<CloudType>::myStochasticDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner)
{}


template<class CloudType>
Foam::myStochasticDispersionRAS<CloudType>::myStochasticDispersionRAS
(
    const myStochasticDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myStochasticDispersionRAS<CloudType>::~myStochasticDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::myStochasticDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    Random& rnd = this->owner().rndGen();

    const scalar C_0 = 4;

    const scalar k = this->kPtr_->primitiveField()[celli];
    const scalar epsilon =
        this->epsilonPtr_->primitiveField()[celli] + rootVSmall;
    
    const scalar sigma2 = sqrt(2*k/3.0);    // the variance, so stddev^2
    
    const vector randNum( rnd.scalarNormal(), rnd.scalarNormal(), rnd.scalarNormal() );
    
    //std::cout << "dt = \"" << dt << "\"" << std::endl;
    //std::cout << "tTurb = \"" << tTurb << "\"" << std::endl;
    //std::cout << "k = \"" << k << "\", epps = \"" << epsilon << "\"" << std::endl;
    //std::cout << "Uc = (" << Uc[0] << "," << Uc[1] << "," << Uc[2] << ")" << std::endl;
    //std::cout << "U = (" << U[0] << "," << U[1] << "," << U[2] << ")" << std::endl;
    //std::cout << "UTurb = (" << UTurb[0] << "," << UTurb[1] << "," << UTurb[2] << ")" << std::endl;
    //std::cout << "randNum = (" << randNum[0] << "," << randNum[1] << "," << randNum[2] << ")" << std::endl;
    
    // this is the equation for homogeneous isotropic turbulence as an explicit scheme from https://baileylab.ucdavis.edu/publications/Bailey_2017_BLM.pdf
    // term1 is the explicit scheme term, term2 is the drift term, term3 is what I call the wiggle term
    // notice that sigma2 is the variance, or stddev^2 (sigma is stddev), the weiner process has mean 0 variance dt, so it is sqrt(dt) not dt in the eqn
    UTurb = UTurb - C_0*epsilon/(2*sigma2)*UTurb*dt + sqrt(C_0*epsilon*dt)*randNum;
    
    tTurb += dt;
    if ( UTurb[0] == 0 && UTurb[1] == 0 && UTurb[2] )
    {
        tTurb = 0;
    }
    
    //std::cout << "new tTurb = \"" << tTurb << "\"" << std::endl;
    //std::cout << "new UTurb = (" << UTurb[0] << "," << UTurb[1] << "," << UTurb[2] << ")\n" << std::endl;

    return Uc + UTurb;
}


// ************************************************************************* //
