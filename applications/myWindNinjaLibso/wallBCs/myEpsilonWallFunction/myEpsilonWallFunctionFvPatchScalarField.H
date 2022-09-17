/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::myEpsilonWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulence dissipation wall function
    condition for high Reynolds number, turbulent flow cases.

    The condition can be applied to wall boundaries for which it
    - calculates \c epsilon and \c G
    - inserts near wall epsilon values directly into the epsilon equation
        to act as a constraint

    where

    \vartable
        epsilon | turbulence dissipation field
        G       | turbulence generation field
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            myEpsilonWallFunction;
    }
    \endverbatim

See also
    Foam::fixedInternalValueFvPatchField

SourceFiles
    myEpsilonWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef myEpsilonWallFunctionFvPatchScalarField_H
#define myEpsilonWallFunctionFvPatchScalarField_H

#include "fixedInternalValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class momentumTransportModel;

/*---------------------------------------------------------------------------*\
           Class myEpsilonWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class myEpsilonWallFunctionFvPatchScalarField
:
    public fixedInternalValueFvPatchField<scalar>
{
protected:

    // Protected data

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();


public:

    //- Runtime type information
    TypeName("myEpsilonWallFunction");


    // Constructors

        //- Construct from patch and internal field
        myEpsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myEpsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  myEpsilonWallFunctionFvPatchScalarField
        //  onto a new patch
        myEpsilonWallFunctionFvPatchScalarField
        (
            const myEpsilonWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        myEpsilonWallFunctionFvPatchScalarField
        (
            const myEpsilonWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new myEpsilonWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        myEpsilonWallFunctionFvPatchScalarField
        (
            const myEpsilonWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new myEpsilonWallFunctionFvPatchScalarField(*this, iF)
            );
        }

    //- Destructor
    virtual ~myEpsilonWallFunctionFvPatchScalarField()
    {}


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate(const Pstream::commsTypes);


        // I-O

            //- Write
            virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
