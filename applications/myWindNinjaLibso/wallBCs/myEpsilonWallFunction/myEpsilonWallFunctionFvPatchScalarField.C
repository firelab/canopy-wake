/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

\*---------------------------------------------------------------------------*/

#include "myEpsilonWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

#include "wallFvPatch.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::myEpsilonWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("myEpsilonWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from patch and internal field
Foam::myEpsilonWallFunctionFvPatchScalarField::
myEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}

//- Construct from patch, internal field and dictionary
Foam::myEpsilonWallFunctionFvPatchScalarField::
myEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{
    checkType();
}

//- Construct by mapping given
//  myEpsilonWallFunctionFvPatchScalarField
//  onto a new patch
Foam::myEpsilonWallFunctionFvPatchScalarField::
myEpsilonWallFunctionFvPatchScalarField
(
    const myEpsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{
    checkType();
}

//- Copy constructor
Foam::myEpsilonWallFunctionFvPatchScalarField::
myEpsilonWallFunctionFvPatchScalarField
(
    const myEpsilonWallFunctionFvPatchScalarField& ptf
)
:
    fixedInternalValueFvPatchField<scalar>(ptf),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{
    checkType();
}

//- Construct and return a clone
// implemented in .H file

//- Copy constructor setting internal field reference
Foam::myEpsilonWallFunctionFvPatchScalarField::
myEpsilonWallFunctionFvPatchScalarField
(
    const myEpsilonWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, iF),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{
    checkType();
}

//- Construct and return a clone setting internal field reference
// implemented in .H file


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myEpsilonWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    
    const label patchI = patch().index();

    const momentumTransportModel& turbulence =
        db().lookupObject<momentumTransportModel>("momentumTransport");
    //const momentumTransportModel& turbulence =
    //    db().lookupObject<momentumTransportModel>
    //    (
    //        IOobject::groupName
    //        (
    //            momentumTransportModel::typeName,
    //            internalField().group()
    //        )
    //    );
    
    const scalarField& y = turbulence.y()[patchI];

    const scalar Cmu25 = pow025(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);
    
    //volScalarField& G =
    //    const_cast<volScalarField&>
    //    (
    //        db().lookupObject<volScalarField>
    //        (
    //            turbulence.type() + ".G"
    //        )
    //    );
    DimensionedField<scalar, volMesh>& G =
        const_cast<DimensionedField<scalar, volMesh>&>
        (
            db().lookupObject<DimensionedField<scalar, volMesh>>(turbulence.GName())
        );

    //DimensionedField<scalar, volMesh>& epsilon =
    //    const_cast<DimensionedField<scalar, volMesh>&>
    //    (
    //        //dimensionedInternalField()
    //        internalField()
    //    );
    DimensionedField<scalar, volMesh>& epsilon = 
        const_cast<DimensionedField<scalar, volMesh>&>(internalField());

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tnu = turbulence.nu();
    const scalarField& nuw = tnu().boundaryField()[patchI];
    //const tmp<scalarField> tnuw = turbulence.nu(patchI);
    //const scalarField& nuw = tnuw();

    const tmp<volScalarField> tnut = turbulence.nut();
    const volScalarField& nut = tnut();
    const scalarField& nutw = nut.boundaryField()[patchI];
    //const nutWallFunctionFvPatchScalarField& nutw =
    //    nutWallFunctionFvPatchScalarField::nutw(turbulence, patchi);

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set epsilon and G
    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        epsilon[faceCellI] = Cmu75*pow(k[faceCellI], 1.5)/(kappa_*y[faceI]);

        G[faceCellI] =
            (nutw[faceI] + nuw[faceI])
           *magGradUw[faceI]
           *Cmu25*sqrt(k[faceCellI])
           /(kappa_*y[faceI]);
    }
    
    
    fixedInternalValueFvPatchField<scalar>::updateCoeffs();
}

void Foam::myEpsilonWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedInternalValueFvPatchField<scalar>::evaluate(commsType);
}

void Foam::myEpsilonWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchField<scalar>::write(os);
    writeEntry(os, "Cmu", Cmu_);
    writeEntry(os, "kappa", kappa_);
    writeEntry(os, "E", E_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        myEpsilonWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
