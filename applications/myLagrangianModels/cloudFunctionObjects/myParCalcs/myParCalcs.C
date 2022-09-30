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
    
    // turns out that while the injector information exists at constructor time, the injector information isn't filled until AFTER constructor
    // time, so this actually needs called after constructor time but before the calcs. Unfortunately that isn't possible, the closest
    // thing is to call it every time even though these values only need calculated once.
    
    nParTotal_ = 0.0;   // start nParTotal_ at zero and add the total nPar to be released over all simulation timesteps from each injector
    massTotal_ = 0.0;   // start massTotal_ at zero and add the total mass to be released over all simulation timesteps from each injector
    InjectionModelList<CloudType>& injectorList = this->owner().injectors();    // get the reference to the list of injectors
    forAll(injectorList, i)     // loop over each injector in the injector list
    {
        // set the reference to the current injector
        InjectionModel<CloudType>& currentInjectionModel = injectorList.operator[](i);
        
        // get the total number of particles/parcels to inject over the simulation for the current injector
        // so get the start and stop times of the injection, and calculate it from parcelsToInject
        scalar currentInjectionTimeStart = currentInjectionModel.timeStart();
        scalar currentInjectionTimeEnd = currentInjectionModel.timeEnd();
        //std::cout << "timeStart = \"" << currentInjectionTimeStart << "\", timeEnd \"" << currentInjectionTimeEnd << "\"" << std::endl;
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
    // I thought about setting the value to 1000 instead, so you don't get weird mass fractions greater than 1, but since that acts as a secondary
    // warning that something is going on, I think the value of 1 is better.
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
    dataPtr_numPar->write();
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
    
    dataPtr_numPar.reset
    (
        new volScalarField
        (
            IOobject
            (
                "numPar",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
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


template<class CloudType>
void Foam::myParCalcs<CloudType>::resetVars()
{
    // reset the values of the variables to zero for the next iteration of calculations
    
    // initialize and set the calculated val data pointers
    dataPtr_numPar->primitiveFieldRef() = 0.0;
    dataPtr_numDensity->primitiveFieldRef() = 0.0;
    dataPtr_conc->primitiveFieldRef() = 0.0;
    dataPtr_volFrac->primitiveFieldRef() = 0.0;
    dataPtr_massFrac->primitiveFieldRef() = 0.0;
    
    // initialize and set the cloud data pointers
    dataPtr_cloud_theta->primitiveFieldRef() = 0.0;
    dataPtr_cloud_alpha->primitiveFieldRef() = 0.0;
    dataPtr_cloud_rhoEff->primitiveFieldRef() = 0.0;
    
    // reset the pastIter variables
    pastIter_origId = 0.0;
    pastIter_numPar = 0.0;
    pastIter_numDensity = 0.0;
    pastIter_conc = 0.0;
    pastIter_volFrac = 0.0;
    pastIter_massFrac = 0.0;
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
    dataPtr_numPar(nullptr),
    dataPtr_numDensity(nullptr),
    dataPtr_conc(nullptr),
    dataPtr_volFrac(nullptr),
    dataPtr_massFrac(nullptr),
    dataPtr_cloud_theta(nullptr),
    dataPtr_cloud_alpha(nullptr),
    dataPtr_cloud_rhoEff(nullptr),
    nParTotal_(1),  // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    massTotal_(1),   // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    pastIter_origId(0),
    pastIter_numPar(0),
    pastIter_numDensity(0),
    pastIter_conc(0),
    pastIter_volFrac(0),
    pastIter_massFrac(0)
{
    // define the pointers to the calculation data right from the start
    // then later just need to reset the variables each iteration of calculations
    initializeDataPointers();
}


template<class CloudType>
Foam::myParCalcs<CloudType>::myParCalcs
(
    const myParCalcs<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    dataPtr_numPar(nullptr),
    dataPtr_numDensity(nullptr),
    dataPtr_conc(nullptr),
    dataPtr_volFrac(nullptr),
    dataPtr_massFrac(nullptr),
    dataPtr_cloud_theta(nullptr),
    dataPtr_cloud_alpha(nullptr),
    dataPtr_cloud_rhoEff(nullptr),
    nParTotal_(1),  // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    massTotal_(1),   // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    pastIter_origId(0),
    pastIter_numPar(0),
    pastIter_numDensity(0),
    pastIter_conc(0),
    pastIter_volFrac(0),
    pastIter_massFrac(0)
{
    // define the pointers to the calculation data right from the start
    // then later just need to reset the variables each iteration of calculations
    initializeDataPointers();
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
    resetVars();
    
    // calculate the total nPar and mass information
    calc_totalVals();
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
    //// do the summation calculations, handling particles/parcels going multiple time iterations before switching particles/parcels, for each calculated variable
    ////  to handle particles going multiple time iterations before switching particles/parcels, past values are stored as temporary values each iteration
    ////  so only the values of the final time a particle/parcel has moved get added to the summation total.
    ////  So values are only temporary except during the final timestep a particle/parcel moves, in which they are added to the summation total.
    ////  The final timestep is used instead of the first timestep as the write times are always at the end of a set of timesteps. Technically, the values only really matter
    ////  at the write times, and this allows the values to match up at the right times. It does mean extra (double at least) temporary storage and some quirky coding though,
    ////  especially with the quirkiness of the particle/parcel iteration loops, which is annoying, but at least it is mostly straightforward.
    
    // get the current particle origId
    scalar currentIter_origId = p.origId();
    
    // detect whether the particle/parcel iteration loop started over without a call to preEvolve()
    // and if the particle/parcel iteration loop started over, need to reset the vars before any more 
    // calculations in prep for the current new set of particle/parcel iterations
    //  this check SHOULD work, so long as the particle/parcel iteration loop doesn't start over in the MIDDLE
    //  of the particle/parcel iteration loop instead of at the beginning. For sure the initial particle/parcel injections
    //  seem to be following this methology, still need to test it to see if it still holds when releasing particles/parcels
    //  over multiple timesteps, but it is possible that when particles/parcels need to switch processors, the list
    //  starts over again. If that happens in a way that causes problems, it doesn't seem like an easy thing to fix.
    //   Turns out that this works even when releasing particles/parcels over multiple timesteps, but not for the reasons I 
    //   originally expected. At the start of each simulation timestep till the injections are done, postMove() is called for
    //   the particles/parcels to be injected, and those particles/parcels are appended to the list of particles/parcels.
    //   So the origId does NOT start at zero for simulation timesteps where injections occur except for with the very first
    //   injection. The only reason my summation method ends up working, is because preEvolve() is always called just before
    //   the new simulation timestep, so just before each injection, so the summation variables start out as zero. Then, after
    //   the injection is done, the full particle/parcel loop is started over at zero and because the pastIter_origId is the
    //   value of the very last injected particle, this reset is called.
    //   So it still worked, but not in the way I originally intended.
    //   It turns out the method always drops the last one or two particles/parcels, and if there is only one dt 
    //   for the first iteration, the first particle/parcel is also dropped. Would need access to iteration information
    //   in a smarter way to fix this problem. At least now it isn't getting double counted.
    if ( currentIter_origId == 0 && pastIter_origId != 0 )
    {
        resetVars();
    }
    
    
    //// calc the current summation values for each calculated variable
    
    // calc the current summation value for number of particles
    scalar currentIter_numPar = p.nParticle();
    
    // calc the current summation value for number density
    scalar currentIter_numDensity = p.nParticle();
    
    // calc the current summation value for concentration
    scalar currentIter_conc = p.nParticle()*p.mass();
    
    // calc the current summation value for volume fraction
    scalar currentIter_volFrac = p.nParticle()*p.volume();
    
    // calc the current summation value for mass fraction
    scalar currentIter_massFrac = p.nParticle()*p.mass();
    
    //// if the current particle/parcel origId is different from the last iteration particle/parcel origId,
    //// add the past iteration summation values to the data storage values
    // this SHOULD avoid being off by one for the first iteration, because pastIter_origId always starts at zero and the first currentIter_origId should also be zero
    // and this SHOULD be restarted each new loop over the particles/parcels. If not, then pastIter_oridId would start at the last particle value when currentIter_origId
    // restarts at zero, so there would be an extra counting, an extra counting for each iteration past the first one of a series.
    // From my understanding, the code does a full simulation time for a given particle/parcel, so preEvolve() SHOULD be called resetting the value of pastIter_origId as desired
    // just be aware that this could potentially be not so straightforward as it looks.
    // Turns out preEvolve() wasn't getting called at the start of each particle/parcel loop, I fixed it with a call to reset just above. So long as the check method holds,
    // this should work.
    // If you wanted the first timestep rather than the last timestep of information for a given particle/parcel, would just add current iteration values to the data storage
    // rather than the past iteration values, so there wouldn't need to be so much storage and the coding of this could be condensed.
    //std::cout << "currentIter_origId = \"" << currentIter_origId << "\", pastIter_origId = \"" << pastIter_origId << "\", currentIter_numPar = \"" << currentIter_numPar << "\", pastIter_numPar = \"" << pastIter_numPar << "\"" << std::endl;
    if ( currentIter_origId != pastIter_origId )
    {
        // do the summation calculations for number of particles
        volScalarField& numPar = dataPtr_numPar();
        numPar[p.cell()] += pastIter_numPar;
        //std::cout << "numPar[" << p.cell() << "] = \"" << numPar[p.cell()] << "\"" << std::endl;
        
        // do the summation calculations for number density
        volScalarField& numDensity = dataPtr_numDensity();
        numDensity[p.cell()] += pastIter_numDensity;
        
        // do the summation calculations for concentration
        volScalarField& conc = dataPtr_conc();
        conc[p.cell()] += pastIter_conc;
        
        // do the summation calculations for volme fraction
        volScalarField& volFrac = dataPtr_volFrac();
        volFrac[p.cell()] += pastIter_volFrac;
        
        // do the summation calculations for mass fraction
        volScalarField& massFrac = dataPtr_massFrac();
        massFrac[p.cell()] += pastIter_massFrac;
    }
    
    
    //// for all iterations, set the current summation values to the past iteration summation values, to get ready for the next iteration
    
    // reset the origId in prep for the next iteration
    pastIter_origId = currentIter_origId;
    
    // reset the past iteration summation value for number of particles
    pastIter_numPar = currentIter_numPar;
    
    // reset the past iteration summation value for number density
    pastIter_numDensity = currentIter_numDensity;
    
    // reset the past iteration summation value for concentration
    pastIter_conc = currentIter_conc;
    
    // reset the past iteration summation value for volume fraction
    pastIter_volFrac = currentIter_volFrac;
    
    // reset the past iteration summation value for mass fraction
    pastIter_massFrac = currentIter_massFrac;
}

template<class CloudType>
void Foam::myParCalcs<CloudType>::postEvolve()
{
    // do the final calculation for each calculated variable
    
    const fvMesh& mesh = this->owner().mesh();
    
    // no final calculation for number of particles, it is already calculated with the summation
    
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
    theta.primitiveFieldRef() = this->owner().theta();
    
    // get the values of the cloud alpha, the cloud particle mass fraction field
    volScalarField& alpha = dataPtr_cloud_alpha();
    alpha.primitiveFieldRef() = this->owner().alpha();
    
    // get the values of the cloud rhoEff, the cloud particle effective density field (looks more like a concentration to me)
    volScalarField& rhoEff = dataPtr_cloud_rhoEff();
    rhoEff.primitiveFieldRef() = this->owner().rhoEff();
    

    CloudFunctionObject<CloudType>::postEvolve();
}

// ************************************************************************* //
