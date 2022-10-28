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

#include "myCuboidInjection.H"
#include "TimeFunction1.H"
#include "distributionModel.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::myCuboidInjection<CloudType>::setVelocityType()
{
    const word velocityType =
        this->coeffDict().template lookup<word>
        (
            "velocityType"
        );  // notice this disallows any default value

    if (velocityType == "constantVelocity")
    {
        velocityType_ = vtConstantVelocity;
        
        U0_ = this->coeffDict().lookup("U0");
    }
    else if (velocityType == "cellCenterValues")
    {
        velocityType_ = vtCellCenterValues;
    }
    else
    {
        FatalErrorInFunction
            << "velocityType must be either 'constantVelocity' or 'cellCenterValues'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myCuboidInjection<CloudType>::myCuboidInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    duration_(this->coeffDict().template lookup<scalar>("duration")),
    parcelsPerSecond_
    (
        this->coeffDict().template lookup<scalar>("parcelsPerSecond")
    ),
    flowRateProfile_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()
        )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
    velocityType_(vtConstantVelocity),  // setVelocityType() sets this later, disallowing a default value. Set this default value just to make the initializer list happy
    U0_(0,0,0),  // gets filled during setVelocityType() if it is actually used
    xMin_(this->coeffDict().template lookup<scalar>("xMin")),
    xMax_(this->coeffDict().template lookup<scalar>("xMax")),
    yMin_(this->coeffDict().template lookup<scalar>("yMin")),
    yMax_(this->coeffDict().template lookup<scalar>("yMax")),
    zMin_(this->coeffDict().template lookup<scalar>("zMin")),
    zMax_(this->coeffDict().template lookup<scalar>("zMax"))
{
    duration_ = owner.db().time().userTimeToTime(duration_);
    
    setVelocityType();
    
    // Set total volume/mass to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_);
    
    updateMesh();
}


template<class CloudType>
Foam::myCuboidInjection<CloudType>::myCuboidInjection
(
    const myCuboidInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    velocityType_(im.velocityType_),
    U0_(im.U0_),
    xMin_(im.xMin_),
    xMax_(im.xMax_),
    yMin_(im.yMin_),
    yMax_(im.yMax_),
    zMin_(im.zMin_),
    zMax_(im.zMax_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myCuboidInjection<CloudType>::~myCuboidInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::myCuboidInjection<CloudType>::updateMesh()
{
    // don't need to do anything fancy to set/calculate additional geometry info for doing release locations later
    // the vals are already set without calculation at constructor time for this injector
    // hmm, but it would be wise to check the input geometry information values to make sure that they make sense
    // make sure each min is less than each max, make sure each pair of min and max values is within the domain
    // not sure that there is a way to check to make sure each pair of min and max values is not inside a solid surface
    //  probably would be wise to allow it just in case a simulation starts out with particles already stuck on a wall
    //  or something like that
    
    if ( xMin_ > xMax_ )
    {
        std::cout << "!!! myCuboidInjection Error !!!  xMin is greater than xMax !!!" << std::endl;
    }
    if ( yMin_ > yMax_ )
    {
        std::cout << "!!! myCuboidInjection Error !!!  yMin is greater than yMax !!!" << std::endl;
    }
    if ( zMin_ > zMax_ )
    {
        std::cout << "!!! myCuboidInjection Error !!!  zMin is greater than zMax !!!" << std::endl;
    }
    
    // actually, if pars are outside the domain, it should just die or just allow it and few particles would be seen
    // would be nice to see a warning, but it would take some time to figure out how to get the domain edge information
    // so I'm going to skip that check
}


template<class CloudType>
Foam::scalar Foam::myCuboidInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::myCuboidInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ( time0 >= 0.0 && time0 < duration_ )
    {
        return floor(parcelsPerSecond_*(time1 - time0));
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::myCuboidInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ( time0 >= 0.0 && time0 < duration_ )
    {
        return flowRateProfile_.integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::myCuboidInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    Random& rndGen = this->owner().rndGen();
    
    scalar xPos = rndGen.scalarAB(xMin_,xMax_);  // draw a random value between xMin and xMax
    scalar yPos = rndGen.scalarAB(yMin_,yMax_);  // draw a random value between yMin and yMax
    scalar zPos = rndGen.scalarAB(zMin_,zMax_);  // draw a random value between zMin and zMax
    
    
    position[0] = xPos;
    position[1] = yPos;
    position[2] = zPos;
    
    // now use the position to set the rest of the cell information for this particle/parcel
    this->findCellAtPosition
    (
        cellOwner,
        tetFacei,
        tetPti,
        position,
        false
    );
}


template<class CloudType>
void Foam::myCuboidInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    if ( velocityType_ == vtConstantVelocity )
    {
        parcel.U() = U0_;
    }
    else // if ( velocityType_ == vtCellCenterValues )
    {
        // fill the parcel initial velocity with the eulerian grid cell center velocities
        // The eulerian velocities are stored in the owner (cloud) as U(). The parcel.cell() 
        //  is a reference to cellOwner set during setPositionAndCell(), it is a celli value.
        parcel.U() = this->owner().U()[parcel.cell()];
    }

    // set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::myCuboidInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::myCuboidInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
