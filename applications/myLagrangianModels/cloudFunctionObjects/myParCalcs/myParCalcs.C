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
    // for now, just set nParTotal_ to 1000, volTotal_ to 1000 m^3, and massTotal_ to 1000 kg.
    //nParTotal_ = 1000;
    //volTotal_ = 1000;
    //massTotal_ = 1000;
    
    // turns out that while the injector information exists at constructor time, the injector information isn't filled until AFTER constructor
    // time, so this actually needs called after constructor time but before the calcs. Unfortunately that isn't possible, the closest
    // thing is to call it every time even though these values only need calculated once. I've placed it into preEvolve() with a boolean
    // that forces it to only be called once at the first call to preEvolve() (the first simulation timestep) to get around this problem.
    
    
    // seems like the best method will be to first figure out the start and end simulation times and the simulation timestep, then loop
    // over each simulation timestep summing up the number of particles/parcels to be released over the whole simulation. If injectors have slight
    // randomizations, the value might be a bit off, but it should be close enough. Then, using the total number of particles/parcels to be released, 
    // sample a diameter for each particle/parcel from the particle/parcel size distribution and use that to calculate the volume for that diameter, sample
    // the density for each given particle/parcel diameter and multiply that by the volume for that diameter to calculate the mass for that diameter, 
    // add each of the individual particle/parcel volumes and particle/parcel masses together to get a volume and a mass total.
    //  As the same random number generator is shared by all distributions and random number requests for the cloud, the sample of diameters will not be the same as the
    // sample that will be used when the particles/parcels are finally injected even if a copy of the random number generator is made with the 
    // position of the random number at the same spot. This is because there are other calls to the random number generator in between particle/parcel
    // injections, so the path of the position of the random number in the random number generator will not be the same. Now, sampling a second set of particle/parcel
    // diameters separate from the injection could also throw off the position of the random number generator, but it seems unlikely else why would the same 
    // number generator be used over all kinds of separate random number samplers? I was worried it would throw off a normal distribution somewhere. I think that 
    // so long as the number of particles/random number extractions from the random number generator for a given desired sampling is large before moving on to a new 
    // desired sampling, the distributions should still remain relatively unbiased.
    // 
    // turns out that my method is not good, simulations can have different start and end times than the duration of the particle/parcel release!!!
    // so the value would not be consistent. Turns out the right way to get the simulation start time, end time, and timestep are as follows though:
    // 
    //   const fvMesh& mesh = this->owner().mesh();
    //   scalar simulation_startTime = mesh.time().startTime().scalar();
    //   scalar simulation_endTime = mesh.time().endTime().scalar();
    //   scalar simulation_timeStep = mesh_.time().deltaTValue(); // technically this would be the simulation timestep for the current simulation time at the call of this function, so for this case, this should be the first simulation timestep if it were to be changing over the simulation
    // 
    // The plan to loop over each particle/parcel getting a sample of diameters and density and calculating the volume and mass of a given particle/parcel
    // for each sampled diameter and summing it up as a total volume and total mass does still sound good, but the method to get nTotalPars needs done differently
    // Hmm, could just get the parcelsPerSecond and duration for each injection and multiply them together, would be off a bit if the parcelsToInject() was slightly
    // different from the parcelsPerSecond not dividing evenly or randomness added to the given injection, oh crap parcelsPerSecond isn't guaranteed to be inside
    // each injection, for example cellZoneInjection uses numberDensity instead. Use an if statement using parcelsToInject() for that case? Or use parcelsToInject()
    // from 0 to duration? Hmm, the 0 to duration wouldn't work because of the weird and inconsistent if statements determining whether to inject or not for a given 
    // injection. I guess I should just use the old method, which is to just use parcelsToInject() for each injection, using the startTime() and endTime() returned 
    // for each injection. This method seems like it would fail in that the injections that only release over a single time like to set timeEnd() to the same value
    // as timeStart(), timeEnd() - timeStart() would then be zero, causing problems in the inconsistent if statements. It does look like the if statements for these
    // particular single time cases handle it to still output the parcelsToInject correctly though. But to test this, I've only studied 4 of the many injection models
    // it could become a problem in the future, watch to make sure nParTotal doesn't become something weird.
    // 
    // Dang it! I had carefully setup the plan to get the start and end times for a given injection, using them in parcelsToInject() to get nTotal, then looping
    // over each particle/parcel to get a sample of the diameters to calculate the individual and total volume and mass, but I quickly discovered that the size
    // distribution is private with no accessor function within all injection types that have it, and while all injection types I have studied have a size distribution,
    // not all injection types are guaranteed to have a size distribution. To try to get around this problem, I first tried thinking of ways to generate a 
    // temporary list of particles/parcels, but the problem with that is that each time a particle/parcel is generated, even if it is deleted at the end, it adds 
    // to the origID count, throwing off the origID count! So that won't work without writing a separate class that is a duplicate of particle.H type stuff, 
    // so NOT WORTH IT. No way to make a separate second list of origID without writing a new class/container.
    // My second attempt was to try generating the size distribution for myself, so replicate the injection type constructor methods for generating the
    // container for the size distribution, then deleting it after finishing with it. So just generate a second temporary version of the size distribution from the
    // dictionary held by each injection, which DOES have a public accessor function. The problem with this method is that it would fail for injection types
    // that do not have size distributions. Could do a, generate if exists, if not exists, do something else, for those injection models, but then there exists
    // NO method in those injection models for getting a particle/parcel diameter without generating a particle/parcel list, which as discussed is NOT a good idea.
    // I guess you could try to make like a fake distribution model, one that replicates how the diameters behave for those injection models that normally don't have
    // an injector list. But more likely, those injection models would just not be qualifiable for using this myParCalcs function on them. Luckily, all the injection 
    // types that I've used so far, even manual injection, seem to have a private sizeDistribution within them. Yay! So the generating a second temporary size 
    // distribution from dictionary should work as a method.
    // anyhow, it seems like the dictionary making of a temporary distribution model for each injection should work fine for now. OH CRAP, the regular cone injection
    // does NOT allow me to do parcelsToInject() this way, it spews out particles/parcels good if simulations start at zero, but outputs nothing if simulations start
    // later. BLARG!!! THIS MEANS myParCalcs CANNOT APPLY TO THE STANDARD coneInjection TYPE WITHOUT MODIFICATION ANYWAYS NO MATTER WHAT!!! BLARGG!!!
    // I guess this means screw it on trying to work with the standard provided injection types in OpenFOAM, if myParCalcs is to work for a given injection type
    // I just have to always use my own versions of the injection types, modified to the standards required by myParCalcs. This wouldn't be required IF THE STANDARD
    // INJECTION TYPES PROVIDED BY OPENFOAM HAPPENED TO BE CONSISTENT. So basically I am forced to always use modified versions of the standard OpenFOAM injection types
    // that are modified to exist by more consistent standards than the standard OpenFOAM standards. BLARG! But the nice thing about this, is that I can now make the
    // calculation of the total variables become much more consistent and less exhausting within myParCalcs, because I now have better inputs to the dang thing.
    // almost seems like what I am really after, is a myInjectionModel class that is more consistent with what I require for myParCalcs.
    // 
    // I tried requiring all injection types used with myParCalcs to have an accessor function to a sizeDistribution called sizeDistribution().
    // and the parcelsToInject() functions also required to be consistent and work for each startTime() and endTime() returned from the injector type.
    // But it turned out that accessing public NON-virtual functions from the injection types appears to require a definition of a type def or a type name
    // for the given injection, you have to get more specific on the call to the function when it is not virtual, and that method doesn't appear to exist 
    // without a bunch of extra modification. So instead, I focused on generating a second temporary size distribution from the dictionary referenced by the given 
    // injection type, as a virtual function for getting the dict held by a given injection model DID exist without modification to the code.
    // It also turns out that this required less modification of the given OpenFOAM provided injection types, only coneInjection needed updated to be
    // a myConeInjection type, as now only OpenFOAM provided injection types that have problems with parcelsToInject() are the problem.
    // 
    // After I got this running, and I was trying to verify the values to see if they made sense, and it turned out that the size distribution is SLIGHTLY
    // different for the injected diameters, and the temporary ones used to sample the particle/parcel sizes before injection in these calculations. So I 
    // added some extra output to the log to give the user an idea of what the values are, so they can see how they vary from the actual ones used in the
    // injections. Most of these outputs make sense, except for the mean total diameters and total diameters individually to the third power. It turns out that
    // you can turn a total into an average by dividing it by the number of particles/parcels sampled, or multiplying top and bottom of it by the number of particles/parcels
    // sampled. But the totals in an average are not always equal. Diameter total does NOT equal diameter to the third power total, and diameter to the third power
    // is the proper average diameter value to use when doing an average total volume (volume = pi/6*d^3). So (d_2 + d_3 + d_5)/(d_0 + d_1 + d_2 + d_3 + d_4 + d_5) does
    // NOT equal (d_2^3 + d_3^3 + d_5^3)/(d_0^3 + d_1^3 + d_2^3 + d_3^3 + d_4^3 + d_5^3). Yet if you multiply each term by pi/6*rho, so long as rho is constant, those
    // constants just cancel out, hence why the mass fraction and the volume fraction end up being equal. Some quirky stuff, the important part is that if only outputting
    // a single value for diameter information, d^3 is more useful as a quantity than d, so long as the totalling calculations are done correctly.
    // I would also output a total/averaged rho for the particles/parcels, but for now it is always a single constant value. But in the future if it is required
    // not to be a const value, probably will still work as is but I would need to add that calculation to the log output.
    //
    // so I'm good to go, I have specific things to call and use that are consistent each time (and if are not consistent, will be modified within the 
    // injection type itself to make it consistent instead of doing modifications here, yay!!!).
    
    
    std::cout << "running myParCalcs::calc_totalVals()" << std::endl;
    
    
    // initialize or get variables needed for the calculations
    nParTotal_ = 0.0;   // start nParTotal_ at zero and add the total nPar to be released over all simulation timesteps from each injector
    volTotal_ = 0.0;    // start volTotal_ at zero and add the total volume to be released over all simulation timesteps from each injector
    massTotal_ = 0.0;   // start massTotal_ at zero and add the total mass to be released over all simulation timesteps from each injector
    
    InjectionModelList<CloudType>& injectorList = this->owner().injectors();    // get the reference to the list of injectors
    
    // some side variables useful for debugging and log file information
    // total of particle/parcel diameters and total of particle/parcel diameters individually to the third power (for volume equation)
    // would use average for both, but that will be calculated at the end from these quantities, average is just a given total divided by nParTotal_
    scalar d_total = 0.0;    // start it at a value of zero
    scalar d_pow3_total = 0.0;  // start it at a value of zero
    
    
    // now calculate the nParTotal for each injection
    // sample the diameters and particle densities for the nParTotal for each given injection to calculate the volumes 
    // and masses for each particle/parcel of a given injection
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
        // just be aware that it could be a problem, I only studied four different injection models, this sample may not hold over all injection models
        // TURNS OUT THAT THE ORIGINAL coneInjection HAS A DIFFERENT PROBLEM. It inputs the correct timeStart() and timeEnd(), but it calculates parcelsToInject()
        // by subtracting the value from the total already found to be injected, WHICH IS ZERO FOR SIMULATIONS STARTED AND STOPPED OR STARTED AT A LATER TIME THAN 0!!!
        // I got around this by writing a new coneInjection, called myConeInjection, that modifies the parcelsToInject() function to do a more sensible calculation.
        scalar nTotal = currentInjectionModel.parcelsToInject(currentInjectionTimeStart,currentInjectionTimeEnd);
        if ( nTotal == 0.0 )
        {
            std::cout << "    !!! Warning !!! myParCalcs::calc_totalVals() function: nParTotal_ for injection type number \"" << i << "\" found to be 0 !!! This happens if an injection type does not output the right number of particles/parcels for calls to parcelsToInject() using the timeStart() and timeEnd() called from the injection. Check the parcelsToInject() code within injection type number \"" << i << "\" to make sure it is coded up correctly!!!"
                  << std::endl;
        }
        
        // increment the total number of particles/parcels to be released using the number of particles/parcels to be released over the simulation for the current injector
        nParTotal_ += nTotal;
        
        
        // need to get the sizeDistribution for the current injector for sampling particle/parcel diameters
        // create it from the dictionary of the current injection
        autoPtr<Foam::distributionModel> currentSizeDistribution
        (
            distributionModel::New( currentInjectionModel.coeffDict().subDict("sizeDistribution"), this->owner().rndGen() )
        );
        
        // now loop over each particle/parcel for the current injection and get a sample of the diameter and the density of the particle/parcel
        // to calculate the volume and mass for each sampled particle/parcel, summing up the total to be the total volume and total mass
        // for the given injection
        scalar volTotal = 0.0;  // start the current injection total off as zero
        scalar massTotal = 0.0; // start the current injection total off as zero
        for ( int i = 0; i < nTotal; i++ )
        {
            // this ends up sampling the size distribution twice for each particle/parcel diameter, once here, and once later when they are finally injected
            // but so long as there are enough particles/parcels to be released for each simulation, it shouldn't bias the sizeDistribution too much
            // besides, the random number generator used in a given size distribution is shared over all functions and objects within the lagrangian code
            // that need to sample a random number for anything, they all use the same random number generator so all would be biasing it together
            scalar sampledDiameter = currentSizeDistribution->sample();
            
            // calculate the current particle/parcel volume
            scalar sampledVolume = pi/6.0*pow3(sampledDiameter);  // 4/3*pi*r^3 == 4/3*pi*d^3/2^3
            
            // get the current particle/parcel density
            // KinematicCloud::setParcelThermoProperties() shows how individual particles/parcels receive their density
            // looks like density must be constant for all particles/parcels, as it is specified as a single value from within
            // constantProperties. Treat it as if it is changing each sample just in case.
            scalar sampledRho = this->owner().constProps().rho0();
            
            // calculate the current particle/parcel mass
            scalar sampledMass = sampledRho*sampledVolume;
            
            
            // now add the current particle/parcel volume and mass to the totals for the current injection
            volTotal += sampledVolume;
            massTotal += sampledMass;
            
            // also add the current particle/parcel diameter and diameter to the third power to the totals of these quantities
            // could also add the current particle/parcel rho, but it is const for all particles/parcels for now. Would just need to
            // add it here if that ever changed.
            d_total += sampledDiameter;
            d_pow3_total += pow3(sampledDiameter);
        }
        
        // increment the overall total volume and total mass of particles/parcels to be released over the simulation
        volTotal_ += volTotal;
        massTotal_ += massTotal;
    }
    
    // output the found information to the log file
    // calculate the mean/average form of the total particle/parcel diameters. Output the total instead of the mean for the
    //  total particle/parcel diameters individually to the third power, demicals look better with totals for this selection
    //  it is easy for the user to convert from a total particle/parcel diameters to a mean/average, average is just a given 
    //  total particle/parcel diameters divided by nParTotal_
    scalar d_avg_total = d_total/nParTotal_;
    std::cout << "myParCalcs::calc_totalVals() calculated values are \n\t nParTotal = " << nParTotal_ << "\n\t volTotal = " << volTotal_ << "\n\t massTotal = " << massTotal_
              << "\n\t mean of all particle/parcel diameters = " << d_avg_total << ", total of all particle/parcel diameters individually to the third power = " << d_pow3_total
              << std::endl;
    
    // if nParTotal_, volTotal_, or massTotal_ are still 0, need to warn, and set them to 1, to avoid a divide by 0 error later
    // I thought about setting the value to 1000 instead, so you don't get weird mass fractions greater than 1, but since mass fractions 
    // greater than one acts as a secondary warning that something is going on, I think the value of 1 is better.
    if ( nParTotal_ == 0 )
    {
        std::cout << "    !!! Warning !!! myParCalcs::calc_totalVals() function: nParTotal_ found to be 0, setting it to 1 to avoid divide by 0 problems in calculations"
                  << std::endl;
        nParTotal_ = 1;
    }
    if ( volTotal_ == 0 )
    {
        std::cout << "    !!! Warning !!! myParCalcs::calc_totalVals() function: volTotal_ found to be 0, setting it to 1 to avoid divide by 0 problems in calculations"
                  << std::endl;
        volTotal_ = 1;
    }
    if ( massTotal_ == 0 )
    {
        std::cout << "    !!! Warning !!! myParCalcs::calc_totalVals() function: massTotal_ found to be 0, setting it to 1 to avoid divide by 0 problems in calculations"
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
    dataPtr_volPar->write();
    dataPtr_massPar->write();
    dataPtr_numDensity->write();
    dataPtr_volFrac->write();
    dataPtr_massFrac->write();
    dataPtr_cellFillSpace->write();
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
    dataPtr_volPar.reset
    (
        new volScalarField
        (
            IOobject
            (
                "volPar",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    dataPtr_massPar.reset
    (
        new volScalarField
        (
            IOobject
            (
                "massPar",
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
    dataPtr_cellFillSpace.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellFillSpace",
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
    dataPtr_volPar->primitiveFieldRef() = 0.0;
    dataPtr_massPar->primitiveFieldRef() = 0.0;
    dataPtr_numDensity->primitiveFieldRef() = 0.0;
    dataPtr_volFrac->primitiveFieldRef() = 0.0;
    dataPtr_massFrac->primitiveFieldRef() = 0.0;
    dataPtr_cellFillSpace->primitiveFieldRef() = 0.0;
    
    // reset the pastIter variables
    pastIter_origId = 0.0;
    pastIter_cellIdx = 0.0;
    pastIter_numPar = 0.0;
    pastIter_volPar = 0.0;
    pastIter_massPar = 0.0;
    pastIter_numDensity = 0.0;
    pastIter_volFrac = 0.0;
    pastIter_massFrac = 0.0;
    pastIter_cellFillSpace = 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from dictionary
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
    dataPtr_volPar(nullptr),
    dataPtr_massPar(nullptr),
    dataPtr_numDensity(nullptr),
    dataPtr_volFrac(nullptr),
    dataPtr_massFrac(nullptr),
    dataPtr_cellFillSpace(nullptr),
    calcTotals_(true),  // start it out as true, gets set to false when the totals get calculated
    nParTotal_(1),  // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    volTotal_(1),  // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    massTotal_(1),   // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    pastIter_origId(0),
    pastIter_cellIdx(0),
    pastIter_numPar(0),
    pastIter_volPar(0),
    pastIter_massPar(0),
    pastIter_numDensity(0),
    pastIter_volFrac(0),
    pastIter_massFrac(0),
    pastIter_cellFillSpace(0)
{
    // define the pointers to the calculation data right from the start
    // then later just need to reset the variables each iteration of calculations
    initializeDataPointers();
}

//- Construct copy
template<class CloudType>
Foam::myParCalcs<CloudType>::myParCalcs
(
    const myParCalcs<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    dataPtr_numPar(nullptr),
    dataPtr_volPar(nullptr),
    dataPtr_massPar(nullptr),
    dataPtr_numDensity(nullptr),
    dataPtr_volFrac(nullptr),
    dataPtr_massFrac(nullptr),
    dataPtr_cellFillSpace(nullptr),
    calcTotals_(true),  // start it out as true, gets set to false when the totals get calculated
    nParTotal_(1),  // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    volTotal_(1),  // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    massTotal_(1),   // default value, a best choice doesn't necessarily exist because dividing by zero or a number smaller than 1 blows up, so I chose a value of 1 to be on the safe side
    pastIter_origId(0),
    pastIter_cellIdx(0),
    pastIter_numPar(0),
    pastIter_volPar(0),
    pastIter_massPar(0),
    pastIter_numDensity(0),
    pastIter_volFrac(0),
    pastIter_massFrac(0),
    pastIter_cellFillSpace(0)
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
    
    // calculate the total nPar, volume, and mass information
    // set to only be calculated ONCE, the first time preEvolve() is called, even though preEvolve() is called at every iteration
    // so it acts like a late call to the constructor of the values
    if ( calcTotals_ == true )
    {
        calc_totalVals();
        calcTotals_ = false;
    }
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
    //// do the summation calculations, handling particles/parcels going multiple time iterations before switching particles/parcels.
    ////
    //// for a given simulation timestep, postMove() appears to be called first to inject all the particles/parcels to be injected for that
    //// simulation timestep, going from nPar_t-1 to nPar_t, then second for all particle/parcel timesteps of a given particle/parcel still
    //// in the domain. So going from 0 to the simulation timestep broken up by various dt for a given particle/parcel, once in order for each of
    //// the particles/parcels. So a for 0 to simulation timestep loop within a for 0 to nPar_t loop.
    //// Looking at the code, these for loops are actually while loops, where it appears possible to kick it out and start over the particle/parcel loop
    //// out of order when hitting faces or switching processors. If this is the case, this algorythm is not correct and it is required to seek the 
    //// specific iteration variables from the code that is calling postMove() indirectly, probably by using this->owner() to get ahold of the 
    //// particle/parcel timestep fraction. But so far, it seems to be working as desired, so hopefully it is good enough as is. Actually, it might still 
    //// kick out of the loop, but without starting over the particle/parcel list or the timestep fraction, so like a pause to store or do extra weird things
    //// that are harder to handle, it is tough to tell. But if this is the case, then it would be resuming at the same spot in the for loops, making this
    //// current method work better than expected. I guess just be aware that I wasn't sure which of these behaviors was correct when writing this code.
    //// 
    //// To handle particles/parcels going multiple time iterations before switching particles/parcels, past values are stored as temporary values each iteration
    //// so only the values of the final time a particle/parcel has moved get added to the summation total. So values are only temporary except during the 
    //// final timestep a particle/parcel moves, in which they are added to the summation total one iteration later.
    //// The final timestep is used instead of the first timestep as the write times are always at the end of a set of timesteps. Technically, the values only really matter
    //// at the write times, and this allows the values to match up at the right times. It does mean extra (double at least) temporary storage and some quirky coding though,
    //// especially with the quirkiness of the particle/parcel iteration loops, which is annoying, but at least it is mostly straightforward.
    ////
    //// All values are reset whenever the particle/parcel list is detected to start over, but the method to detect when it starts over is a bit faulty and calls the reset
    //// for all timesteps in which the parID is 0, not just the first iteration after the particle/parcel loop starts over. Luckily, the when to sum method is such that 
    //// these values do not matter, the final values of parID of 0 are stored as the pastIter_ values the final timestep iteration over parID of 0, and the when to sum 
    //// then detects a change in parID for the next timestep and stores the final value. So the right value is still stored despite this extra resetting, the extra resetting
    //// shouldn't be a problem. However, the when to sum algorythm is off by one, leading to the final timestep of the final particle/parcel pastIter_ values not being added
    //// to the sum during postMove(). It can be shown that all final timestep particle/parcel pastIter_ values are added to the sum during postMove() except for the final 
    //// timestep particle/parcel pastIter_ values of the final particle/parcel, which is always skipped. Luckily, the pastIter_ values stored and left in storage at the end
    //// of this final call to postMove() happen to be the desired final timestep particle/parcel pastIter_ values of the final particle/parcel, so the final summation just
    //// needs done within postEvolve(). But this means the code is even more quirky for this algorythm, almost makes me think it would be better to seek out the timestep
    //// fraction information so that postMove() can finally know where in the iteration loops it is at. But it does work, at least so far, so hopefully it is good enough.
    //// if it ever stops working, I would just skip trying to do it this way and go immediately to seeking the timestep fraction or other iteration information so that postMove()
    //// can know exactly where it is within the loops, it was too much work to figure out this algorythm.
    
    
    // get the current particle origId
    scalar currentIter_origId = p.origId();
    
    
    // detect whether the particle/parcel iteration loop started over without a call to preEvolve()
    // and if the particle/parcel iteration loop started over, need to reset the vars before any more 
    // calculations in prep for the current new set of particle/parcel iterations.
    // note that this gets called for EVERY timestep of the first particle/parcel instead of just the first timestep, which seems bad,
    // but because the final timestep values of the first particle/parcel are the pastIter_ values for the first timestep iteration
    // of the second particle/parcel, and that value is immediately stored, resetting the values for each and every timestep of the 
    // first particle/parcel ends up not being a problem after all.
    // part of why this if statement works is because the particle/parcel loop always goes from particle/parcel 0 to nParticles/nParcels 
    // before starting over, so it completes a full loop in a specific order before starting the next loop. It is possible that this behavior is altered when
    // particles/parcels need to switch processors, so on the edges of how the mesh is broken up by decomposePar. If this behavior does get altered
    // and the values at that time of the simulation actually end up being important, I see no way to improve the algorythm without 
    // accessing the timestep fraction variables within the particle/parcel class that is calling this function (via this->owner() somehow).
    if ( currentIter_origId == 0 )
    {
        resetVars();
    }
    
    
    //// if the current particle/parcel origId is different from the last iteration particle/parcel origId,
    //// add the past iteration summation values to the data storage values
    // Note that this stores the final pastIter_ value for each and every particle/parcel timestep EXCEPT 
    // the last one. But at the end of each and every particle/parcel loop, it is also always at the end of 
    // the timestep loop for the last particle/parcel, so the final pastIter_ values are the pastIter_ values
    // of the final particle/parcel at the final timestep loop for each and every call to postEvolve().
    // Adding that summation value in postEvolve() completes a single summation of the final timestep pastIter_ values for 
    // each and every particle/parcel, despite this postMove() function only knowing current and past iteration values
    // but not able to know when it is at the final timestep iteration without accessing timestep fraction from some kind 
    // of nasty this->owner() call.
    //std::cout << "currentIter_origId = \"" << currentIter_origId << "\", pastIter_origId = \"" << pastIter_origId << "\", pastIter_numPar = \"" << pastIter_numPar << "\"" << std::endl;
    if ( currentIter_origId != pastIter_origId )
    {
        // get the pointers to each of the datasets that will be modified by the summation calculations
        volScalarField& numPar = dataPtr_numPar();
        volScalarField& volPar = dataPtr_volPar();
        volScalarField& massPar = dataPtr_massPar();
        volScalarField& numDensity = dataPtr_numDensity();
        volScalarField& volFrac = dataPtr_volFrac();
        volScalarField& massFrac = dataPtr_massFrac();
        volScalarField& cellFillSpace = dataPtr_cellFillSpace();
        
        // do the summation calculations
        numPar[pastIter_cellIdx] += pastIter_numPar;
        //std::cout << "numPar[" << pastIter_cellIdx << "] = \"" << numPar[pastIter_cellIdx] << "\"" << std::endl;
        volPar[pastIter_cellIdx] += pastIter_volPar;
        massPar[pastIter_cellIdx] += pastIter_massPar;
        numDensity[pastIter_cellIdx] += pastIter_numDensity;
        volFrac[pastIter_cellIdx] += pastIter_volFrac;
        massFrac[pastIter_cellIdx] += pastIter_massFrac;
        cellFillSpace[pastIter_cellIdx] += pastIter_cellFillSpace;
    }
    
    
    //// calculate the summation values for this current iteration and immediately store them within the past iteration temporary storage
    //// in preparation for the next iteration. So combine the calculation of current and the resetting of current values to the past values
    //// all in one step. LHS will be pastIter_ values, RHS will be currentIter_ calculations or values.
    
    // set the pastIter_ origId to the current value in prep for the next iteration
    pastIter_origId = currentIter_origId;
    
    // get the current cellIdx and set it to the pastIter_ value in prep for the next iteration
    pastIter_cellIdx = p.cell();
    
    // calc the summation value for number of particles and set it to the pastIter_ value in prep for the next iteration
    pastIter_numPar = p.nParticle();
    
    // calc the summation value for volume of particles and set it to the pastIter_ value in prep for the next iteration
    pastIter_volPar = p.nParticle()*p.volume();
    
    // calc the summation value for mass of particles and set it to the pastIter_ value in prep for the next iteration
    pastIter_massPar = p.nParticle()*p.mass();
    
    // calc the summation value for number density and set it to the pastIter_ value in prep for the next iteration
    pastIter_numDensity = p.nParticle();
    
    // calc the summation value for volume fraction and set it to the pastIter_ value in prep for the next iteration
    pastIter_volFrac = p.nParticle()*p.volume();
    
    // calc the summation value for mass fraction and set it to the pastIter_ value in prep for the next iteration
    pastIter_massFrac = p.nParticle()*p.mass();
    
    // calc the summation value for the filled cell space and set it to the pastIter_ value in prep for the next iteration
    pastIter_cellFillSpace = p.nParticle()*p.volume();
    
}

template<class CloudType>
void Foam::myParCalcs<CloudType>::postEvolve()
{
    // do the final particle/parcel summation (only necessary because of quirkiness of the calls to these functions)
    // followed by the final calculations for each calculated variable.
    
    
    // first, get the pointers to each of the datasets that will be modified/updated with the calculations
    volScalarField& numPar = dataPtr_numPar();
    volScalarField& volPar = dataPtr_volPar();
    volScalarField& massPar = dataPtr_massPar();
    volScalarField& numDensity = dataPtr_numDensity();
    volScalarField& volFrac = dataPtr_volFrac();
    volScalarField& massFrac = dataPtr_massFrac();
    volScalarField& cellFillSpace = dataPtr_cellFillSpace();
    
    
    // now do the final particle/parcel summation
    // 
    // because postMove() can only know currentIter_ and pastIter_ value information
    // and what was going on the last iteration, but without knowing what will happen at the next iteration
    // or whether it is even at the final iteration; postEvolve() needs to do the summation calculation
    // of the final timestep of the last particle/parcel summation for each calculated variable before doing the 
    // final calculations. Luckily, even though postMove() can't know when it ends, the final call to postMove() before calls to
    // postEvolve() should always be for the final timestep of this last particle/parcel. So calls to postEvolve()
    // SHOULD always start with the pastIter_ variables filled with the pastIter_ values of the final timestep 
    // of the last particle/parcel summation.
    // 
    numPar[pastIter_cellIdx] += pastIter_numPar;
    //std::cout << "numPar[" << pastIter_cellIdx << "] = \"" << numPar[pastIter_cellIdx] << "\"" << std::endl;
    volPar[pastIter_cellIdx] += pastIter_volPar;
    massPar[pastIter_cellIdx] += pastIter_massPar;
    numDensity[pastIter_cellIdx] += pastIter_numDensity;
    volFrac[pastIter_cellIdx] += pastIter_volFrac;
    massFrac[pastIter_cellIdx] += pastIter_massFrac;
    cellFillSpace[pastIter_cellIdx] += pastIter_cellFillSpace;
    
    
    // NOW do the final calculation for each calculated variable
    
    const fvMesh& mesh = this->owner().mesh();
    
    
    // no final calculation for number of particles, it is already calculated with the summation
    
    // no final calculation for volume of particles, it is already calculated with the summation
    
    // no final calculation for mass of particles, it is already calculated with the summation
    
    // do the final calculation for number density
    numDensity.primitiveFieldRef() /= nParTotal_;
    
    // do the final calculation for volume fraction
    volFrac.primitiveFieldRef() /= volTotal_;
    
    // do the final calculation for mass fraction
    massFrac.primitiveFieldRef() /= massTotal_;
    
    // do the final calculation for the filled cell volume
    cellFillSpace.primitiveFieldRef() /= mesh.V();
    

    CloudFunctionObject<CloudType>::postEvolve();
}

// ************************************************************************* //
