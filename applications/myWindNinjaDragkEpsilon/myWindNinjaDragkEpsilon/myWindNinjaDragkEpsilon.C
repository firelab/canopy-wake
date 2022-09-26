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

#include "myWindNinjaDragkEpsilon.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::myWindNinjaDragkEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    
	// my stuff added to the initializer list
	// coefficients
	betaP
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaP",
            this->coeffDict_,
            1.0
        )
    ),
	betaD
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaD",
            this->coeffDict_,
            5.03
        )
    ),
	Cepps4
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cepps4",
            this->coeffDict_,
            0.78
        )
    ),
	Cepps5
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cepps5",
            this->coeffDict_,
            0.78
        )
    ),
	// fields
	Cd
    (
        IOobject
        (
            "Cd",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
	LAD
    (
        IOobject
        (
            "LAD",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
	// end stuff added to the initializer list
    
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::sigma() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigma",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - this->nut_*twoSymm(fvc::grad(this->U_)), // notice the missing dev() wrapper on symmTensor
            k_.boundaryField().types()  // notice the patch types aren't overwritten for those without an equivalent for symmTensor
        )
    );
    
    /* OpenFOAM 8 version of the above code */
    /*tmp<volScalarField> tk(k());

    // Get list of patchField type names from k
    wordList patchFieldTypes(tk().boundaryField().types());

    // For k patchField types which do not have an equivalent for symmTensor
    // set to calculated
    forAll(patchFieldTypes, i)
    {
        if
        (
           !fvPatchField<symmTensor>::patchConstructorTablePtr_
                ->found(patchFieldTypes[i])
        )
        {
            patchFieldTypes[i] = calculatedFvPatchField<symmTensor>::typeName;
        }
    }

    return volSymmTensorField::New
    (
        IOobject::groupName("R", this->alphaRhoPhi_.group()),
        ((2.0/3.0)*I)*tk() - (this->nut_)*dev(twoSymm(fvc::grad(this->U_))),
        patchFieldTypes
    );*/
}

template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::devTau() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devTau",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           (-(this->alpha_*this->rho_*this->nuEff()))*dev(twoSymm(fvc::grad(this->U_))) // notice this is the same
        )
    );
    
    /* OpenFOAM 8 version of the above code */
    /*return volSymmTensorField::New
    (
        IOobject::groupName("devTau", this->alphaRhoPhi_.group()),
        (-(this->alpha_*this->rho_*this->nuEff()))
       *dev(twoSymm(fvc::grad(this->U_)))
    );*/
}

template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::divDevTau
(
    volVectorField& U
) const
{
    // before, I did an info/cout to see which of these functions is called, and whether
    // they were replacing the inherited functions properly. Looks like they are 
    // inheriting properly, but this was the only function actually used and 
    // therefore called in my specific simulation test case
    
    return
    (
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)    // notice this is the same
      - fvc::div(this->alpha_*this->rho_*this->nuEff()*dev(T(fvc::grad(U))))    // notice this uses dev instead of dev2, requiring a different entry in fvSchemes
    );
    
    /* OpenFOAM 8 version of the above code */
    /*return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    );*/
}

template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(this->alpha_*rho*this->nuEff(), U)    // notice that this is the same
      - fvc::div(this->alpha_*rho*this->nuEff()*dev(T(fvc::grad(U))))    // notice that this is using dev instead of dev2
    );
    
    /* OpenFOAM 8 version of the above code */
    /*return
    (
      - fvc::div((this->alpha_*rho*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*rho*this->nuEff(), U)
    );*/
}


template<class BasicMomentumTransportModel>
bool myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        // added scalars to read
		betaP.readIfPresent(this->coeffDict());
		betaD.readIfPresent(this->coeffDict());
		Cepps4.readIfPresent(this->coeffDict());
		Cepps5.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void myWindNinjaDragkEpsilon<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();
    
    volScalarField::Internal G(this->GName(), nut*2*magSqr(symm(fvc::grad(U))) );
    // OpenFOAM 8 version of the above code
    //tmp<volTensorField> tgradU = fvc::grad(U);
    //volScalarField::Internal G
    //(
    //    this->GName(),
    //    nut()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    //);
    //tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + Cd()*LAD()*Cepps4*betaP*pow3(mag(U()))*epsilon_()/k_() - fvm::Sp(   Cd()*LAD()*Cepps5*betaD*mag(U()) , epsilon_   )		// line added. Had a lot of trouble configuring this. In the end I rearranged stuff a bunch and split it into two terms, one that is always positive (explicit) and one that is always negative (implicit). Tested it and whether you do the explicit terms with no wrapper, or with the fvm::Su( , epsilon ) wrapper, it behaves the same and the results did not change. I decided to use the simpler version.
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + Cd()*LAD()*betaP*pow3(mag(U())) - fvm::Sp(   Cd()*LAD()*betaD*mag(U()) , k_   )		// line added. Had a lot of trouble configuring this. In the end I rearranged stuff a bunch and split it into two terms, one that is always positive (explicit) and one that is always negative (implicit). Tested it and whether you do the explicit terms with no wrapper, or with the fvm::Su( , k ) wrapper, it behaves the same and the results did not change. I decided to use the simpler version.
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
