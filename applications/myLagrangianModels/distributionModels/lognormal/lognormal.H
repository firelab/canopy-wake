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

Class
    Foam::distributionModels::lognormal

Description
    A lognormal distribution model

    \verbatim
        model = strength * exp(-0.5*((x - expectation)/variance)^2 )
    \endverbatim

    strength only has meaning if there's more than one distribution model

SourceFiles
    lognormal.C

\*---------------------------------------------------------------------------*/

#ifndef lognormal_H
#define lognormal_H

#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                           Class lognormal Declaration
\*---------------------------------------------------------------------------*/

class lognormal
:
    public distributionModel
{
    // Private Data


        //- Distribution minimum
        scalar minValue_;

        //- Distribution maximum
        scalar maxValue_;


        // Model coefficients

            scalar expectation_;
            scalar variance_;

            scalar a_;


public:

    //- Runtime type information
    TypeName("lognormal");


    // Constructors

        //- Construct from components
        lognormal(const dictionary& dict, Random& rndGen);

        //- Construct copy
        lognormal(const lognormal& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new lognormal(*this));
        }


    //- Destructor
    virtual ~lognormal();


    // Member Functions

        //- Sample the distributionModel
        virtual scalar sample() const;

        //- Return the minimum value
        virtual scalar minValue() const;

        //- Return the maximum value
        virtual scalar maxValue() const;

        //- Return the mean value
        virtual scalar meanValue() const;

        virtual scalar erfInv(const scalar y) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
