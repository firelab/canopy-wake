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

#include "myPatchInjection.H"
#include "TimeFunction1.H"
#include "distributionModel.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::myPatchInjection<CloudType>::setVelocityType()
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
Foam::myPatchInjection<CloudType>::myPatchInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    patchInjectionBase(owner.mesh(), this->coeffDict().lookup("patchName")),
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
    U0_(0,0,0)  // gets filled during setVelocityType() if it is actually used
{
    duration_ = owner.db().time().userTimeToTime(duration_);
    
    setVelocityType();
    
    patchInjectionBase::updateMesh(owner.mesh());

    // Set total volume/mass to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_);
}


template<class CloudType>
Foam::myPatchInjection<CloudType>::myPatchInjection
(
    const myPatchInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    patchInjectionBase(im),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    velocityType_(im.velocityType_),
    U0_(im.U0_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myPatchInjection<CloudType>::~myPatchInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::myPatchInjection<CloudType>::updateMesh()
{
    patchInjectionBase::updateMesh(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::myPatchInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::myPatchInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar nParcels = (time1 - time0)*parcelsPerSecond_;

        Random& rnd = this->owner().rndGen();

        label nParcelsToInject = floor(nParcels);

        // Inject an additional parcel with a probability based on the
        // remainder after the floor function
        if
        (
            nParcelsToInject > 0
         && (nParcels - scalar(nParcelsToInject) > rnd.globalScalar01())
        )
        {
            ++nParcelsToInject;
        }

        return nParcelsToInject;
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::myPatchInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return flowRateProfile_.integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::myPatchInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    patchInjectionBase::setPositionAndCell
    (
        this->owner().mesh(),
        this->owner().rndGen(),
        position,
        cellOwner,
        tetFacei,
        tetPti
    );
}


template<class CloudType>
void Foam::myPatchInjection<CloudType>::setProperties
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
bool Foam::myPatchInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::myPatchInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
