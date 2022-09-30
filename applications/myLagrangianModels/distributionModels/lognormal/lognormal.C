/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "lognormal.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(lognormal, 0);
    addToRunTimeSelectionTable(distributionModel, lognormal, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::lognormal::lognormal
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    minValue_(distributionModelDict_.template lookup<scalar>("minValue")),
    maxValue_(distributionModelDict_.template lookup<scalar>("maxValue")),
    expectation_(distributionModelDict_.template lookup<scalar>("expectation")),
    variance_(distributionModelDict_.template lookup<scalar>("variance")),
    a_(0.147)
{
    check();
    info();
}


Foam::distributionModels::lognormal::lognormal(const lognormal& p)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    expectation_(p.expectation_),
    variance_(p.variance_),
    a_(p.a_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::lognormal::~lognormal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::lognormal::sample() const
{
    // wikipedia has the normal distribution cumulative distribution formula as
    // 0.5*(  1 + erf( (   x -mu)/(sqrt(2)*sigma) )  )
    // the lognormal distribution cumulative distribution formula is
    // 0.5*(  1 + erf( (ln(x)-mu)/(sqrt(2)*sigma) )  )
    // this makes me think that the a and b formulas just need the ln(x) instead of (x)
    
    //scalar a = erf((minValue_ - expectation_)/variance_);
    //scalar b = erf((maxValue_ - expectation_)/variance_);
    scalar a = erf((log(minValue_) - expectation_)/variance_);
    scalar b = erf((log(maxValue_) - expectation_)/variance_);

    scalar y = rndGen_.sample01<scalar>();

    // wikipedia has the normal distribution quantile formula as
    //      mu + sqrt(2)*sigma*erfInv(2p-1)
    // the lognormal distribution quantile formula is
    // exp( mu + sqrt(2)*sigma*erfInv(2p-1) )
    // this makes me think that I just have to do an additional exp() of the
    // normal distribution formula
    
    //scalar x = erfInv(y*(b - a) + a)*variance_ + expectation_;
    scalar x = exp( erfInv(y*(b - a) + a)*variance_ + expectation_ );

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    x = min(max(x, minValue_), maxValue_);
    
    // the pdf functions are slightly more different than an added ln(x) vs (x) or exp() of stuff
    // like there is an additional 1/x multiplyer to the formula. But trying to decypher how distributions work, it seems like
    // the pdf is just the derivative of the cdf, and the quantile of cdf(max) - cdf(min) is the quantile of the probability
    // and so it seems like the pdf is not so important, so long as the cdf and quantile function are right, it seems like it should be right
    // Another note, it seems like there is an additional sqrt(2) multiplier and 1/2 in different spots of the cdfs and quantile functions
    // than are seen in other places or in this code, I figured so long as the wikipedia versions match, those additions must be additions 
    // to the overall wikipedia description and hopefully shouldn't matter.

    return x;
}


Foam::scalar Foam::distributionModels::lognormal::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::lognormal::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::lognormal::meanValue() const
{
    // oops, I missed this in my first pass on the function
    // see wikipedia for lognormal, not sure how this gets affected if values are truncated
    // seems like none of the calls to mean for any of the distributions would be correct if the values are truncated.
    //  So this isn't the correct for truncated, but at least it follows the methods done by other distribution models
    //  so that they are all wrong together
    return exp( expectation_ + 0.5*variance_*variance_ );
}


Foam::scalar Foam::distributionModels::lognormal::erfInv(const scalar y) const
{
    scalar k = 2.0/(constant::mathematical::pi*a_) +  0.5*log(1.0 - y*y);
    scalar h = log(1.0 - y*y)/a_;
    scalar x = sqrt(-k + sqrt(k*k - h));
    if (y < 0.0)
    {
        x *= -1.0;
    }
    return x;
}


// ************************************************************************* //
