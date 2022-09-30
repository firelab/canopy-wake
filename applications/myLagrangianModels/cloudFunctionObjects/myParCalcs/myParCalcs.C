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

#include "myParCalcs.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::myParCalcs<CloudType>::calc_totalVals()
{
    // this should be called by each constructor, just need to get any of these values once,
    // at the beginning before any simulation stuff
    
    // need to access each of the individual injectors using this->owner().injectors() following methods shown in InjectionModelList.C
    // for now, just set nParTotal_ to 1000 and massTotal_ to 1000 kg.
    //nParTotal_ = 1000;
    //massTotal_ = 1000;
    
    
    nParTotal_ = 0.0;   // start nParTotal_ at zero and add the total nPar to be released from each injector over all simulation timesteps
    massTotal_ = 0.0;   // start massTotal_ at zero and add the total mass to be released from each injector over all simulation timesteps
    InjectionModelList<CloudType>& injectorList = this->owner().injectors();    // get the reference to the list of injectors
    forAll(injectorList, i)     // loop over each injector in the injector list
    {
        // set the reference to the current injector
        InjectionModel<CloudType>& currentInjectionModel = injectorList.operator[](i);
        
        // get the total number of particles/parcels to inject over the simulation for the current injector
        // so get the start and stop times of the injection, and calculate it from parcelsToInject
        scalar currentInjectionTimeStart = currentInjectionModel.timeStart();
        scalar currentInjectionTimeEnd = currentInjectionModel.timeEnd();
        // for some injections, timeStart() matches timeEnd() as SOI and SOI rather than SOI and SOI+duration
        // and some parcelsToInject do (timeEnd - timeStart)*parcelsPerSecond which would be zero if this were the case
        // but it SHOULD be all right, all the injections that do (timeEnd - timeStart)*parcelsPerSecond seem to have their
        // timeStart() and timeEnd() functions returning SOI and SOI+duration so that won't happen for them,
        // and the ones that DO have this timeEnd - timeStart as zero seem to use a different calculation that avoids this problem
        // just be aware that it could be a problem, I only studied three different injection models, this sample may not hold over all injection models
        scalar nTotal = currentInjectionModel.parcelsToInject(currentInjectionTimeStart,currentInjectionTimeEnd);
        
        // increment the total number of particles/parcels to be released using the number of particles/parcels to be released over the simulation for the current injector
        nParTotal_ += nTotal;
        
        // get the total mass set by the user for the current injector
        scalar mt = currentInjectionModel.massTotal();
        
        // increment the total mass using the massTotal given by the current injector
        // note that this value may be an arbitrarily set value, even a value of zero, but it is the best we've got.
        //  The number of particles is much more reliable, but without knowing the particle diameters, which aren't consistently set
        //  till after or during injection, nPars can't be converted into a total mass. VolumeTotal is set and used more frequently than
        //  massTotal, but the value of it is not controlled from the constant properties file: some injections set it using particle diameters,
        //  but most set it using an integration of the flow rate profile. So the volume total is not so straightforward to set than mass total
        //  because mass total is set as a single set number instead of being calculated from a profile and sometimes set from diameters.
        //  The important thing is that the user knows what value they set for massTotal for a given injection.
        massTotal_ += mt;
    }
    
    // if nParTotal_ or massTotal_ are 0, need to warn, and set them to 1, to avoid a divide by 0 error later
    if ( nParTotal_ == 0 )
    {
        std::cout << "    !!! Warning !!! myParCalcs nParTotal_ found to be 0, setting it to 1 to avoid divide by 0 problems in calculations"
                  << std::endl;
        nParTotal_ = 1;
    }
    if ( massTotal_ == 0 )
    {
        std::cout << "    !!! Warning !!! myParCalcs massTotal_ found to be 0, setting it to 1 to avoid divide by 0 problems in calculations"
                  << std::endl;
        massTotal_ = 1;
    }
}


// * * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::myParCalcs<CloudType>::write()
{
    // make each of the pointers to the calculated data write
    // to the output files designated during constructor time
    // this function is a virtual function called at each write time
    // so for most of these variables this means writing to the
    // time directories at a given write time
    
    // write for the calculated val data pointers
    dataPtr_numDensity->write();
    dataPtr_conc->write();
    dataPtr_volFrac->write();
    dataPtr_massFrac->write();
    
    // write for the cloud data pointers
    dataPtr_cloud_theta->write();
    dataPtr_cloud_alpha->write();
    dataPtr_cloud_rhoEff->write();
}


template<class CloudType>
void Foam::myParCalcs<CloudType>::initializeDataPointers()
{
    // this is called by each type of constructor, and it expects each constructor
    // to have already initialized each data pointer as a null pointer
    // now the data for each data pointer is initialized and the data pointers are
    // finally set to point to the designated data for all the calculations
    // as the data is initialized, the data storage is given a name and where to write to
    
    
    const fvMesh& mesh = this->owner().mesh();

    // initialize and set the calculated val data pointers
    
    dataPtr_numDensity.reset
    (
        new volScalarField
        (
            IOobject
            (
                "numDensity",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    dataPtr_conc.reset
    (
        new volScalarField
        (
            IOobject
            (
                "conc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    dataPtr_volFrac.reset
    (
        new volScalarField
        (
            IOobject
            (
                "volFrac",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    dataPtr_massFrac.reset
    (
        new volScalarField
        (
            IOobject
            (
                "massFrac",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    
    // initialize and set the cloud data pointers
    
    dataPtr_cloud_theta.reset
    (
        new volScalarField
        (
            IOobject
            (
                this->owner().name() + "Theta",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    dataPtr_cloud_alpha.reset
    (
        new volScalarField
        (
            IOobject
            (
                this->owner().name() + "Alpha",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    dataPtr_cloud_rhoEff.reset
    (
        new volScalarField
        (
            IOobject
            (
                this->owner().name() + "RhoEff",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myParCalcs<CloudType>::myParCalcs
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    dataPtr_numDensity(nullptr),
    dataPtr_conc(nullptr),
    dataPtr_volFrac(nullptr),
    dataPtr_massFrac(nullptr),
    dataPtr_cloud_theta(nullptr),
    dataPtr_cloud_alpha(nullptr),
    dataPtr_cloud_rhoEff(nullptr),
    nParTotal_(1),  // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    massTotal_(1)   // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
{
    // define the pointers to the calculation data right from the start
    // then later just need to reset the variables each iteration of calculations
    initializeDataPointers();
    
    // calculate the total nPar and mass information
    calc_totalVals();
}


template<class CloudType>
Foam::myParCalcs<CloudType>::myParCalcs
(
    const myParCalcs<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    dataPtr_numDensity(nullptr),
    dataPtr_conc(nullptr),
    dataPtr_volFrac(nullptr),
    dataPtr_massFrac(nullptr),
    dataPtr_cloud_theta(nullptr),
    dataPtr_cloud_alpha(nullptr),
    dataPtr_cloud_rhoEff(nullptr),
    nParTotal_(1),  // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    massTotal_(1)   // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
{
    // define the pointers to the calculation data right from the start
    // then later just need to reset the variables each iteration of calculations
    initializeDataPointers();
    
    // calculate the total nPar and mass information
    calc_totalVals();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myParCalcs<CloudType>::~myParCalcs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::myParCalcs<CloudType>::preEvolve()
{
    // reset the values of the variables to zero for the next iteration of calculations
    
    // initialize and set the calculated val data pointers
    dataPtr_numDensity->primitiveFieldRef() = 0.0;
    dataPtr_conc->primitiveFieldRef() = 0.0;
    dataPtr_volFrac->primitiveFieldRef() = 0.0;
    dataPtr_massFrac->primitiveFieldRef() = 0.0;
    
    // initialize and set the cloud data pointers
    dataPtr_cloud_theta->primitiveFieldRef() = 0.0;
    dataPtr_cloud_alpha->primitiveFieldRef() = 0.0;
    dataPtr_cloud_rhoEff->primitiveFieldRef() = 0.0;
}

template<class CloudType>
void Foam::myParCalcs<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point& position0,
    bool& keepParticle
)
{
    // do the summation calculations for each calculated variable
    
    // do the summation calculations for number density
    volScalarField& numDensity = dataPtr_numDensity();
    numDensity[p.cell()] += p.nParticle();
    
    // do the summation calculations for concentration
    volScalarField& conc = dataPtr_conc();
    conc[p.cell()] += p.nParticle()*p.mass();
    
    // do the summation calculations for volme fraction
    volScalarField& volFrac = dataPtr_volFrac();
    volFrac[p.cell()] += p.nParticle()*p.volume();
    
    // do the summation calculations for mass fraction
    volScalarField& massFrac = dataPtr_massFrac();
    massFrac[p.cell()] += p.nParticle()*p.mass();
}

template<class CloudType>
void Foam::myParCalcs<CloudType>::postEvolve()
{
    // do the final calculation for each calculated variable
    
    const fvMesh& mesh = this->owner().mesh();
    
    // do the final calculation for number density
    volScalarField& numDensity = dataPtr_numDensity();
    numDensity.primitiveFieldRef() /= mesh.V();
    
    // do the final calculation for concentration
    volScalarField& conc = dataPtr_conc();
    conc.primitiveFieldRef() /= mesh.V();
    
    // do the final calculation for volume fraction
    volScalarField& volFrac = dataPtr_volFrac();
    volFrac.primitiveFieldRef() /= mesh.V();
    
    // do the final calculation for mass fraction
    volScalarField& massFrac = dataPtr_massFrac();
    massFrac.primitiveFieldRef() /= massTotal_;
    
    
    // get the values for the cloud data pointers, by calling the functions from the cloud
    
    // get the values of the cloud theta, the cloud particle volume fraction field
    volScalarField& theta = dataPtr_cloud_theta();
    theta.primitiveFieldRef() += this->owner().theta();
    
    // get the values of the cloud alpha, the cloud particle mass fraction field
    volScalarField& alpha = dataPtr_cloud_alpha();
    alpha.primitiveFieldRef() += this->owner().alpha();
    
    // get the values of the cloud rhoEff, the cloud particle effective density field (looks more like a concentration to me)
    volScalarField& rhoEff = dataPtr_cloud_rhoEff();
    rhoEff.primitiveFieldRef() += this->owner().rhoEff();
    

    CloudFunctionObject<CloudType>::postEvolve();
}

// ************************************************************************* //
