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

#include "myConeInjection.H"
#include "TimeFunction1.H"
#include "Constant.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::myConeInjection<CloudType>::setInjectionMethod()
{
    const word injectionMethod =
        this->coeffDict().template lookupOrDefault<word>
        (
            "injectionMethod",
            word::null
        );

    if (injectionMethod == "point" || injectionMethod == word::null)
    {
        injectionMethod_ = imPoint;

        updateMesh();
    }
    else if (injectionMethod == "disc")
    {
        injectionMethod_ = imDisc;

        this->coeffDict().lookup("dInner") >> dInner_;
        this->coeffDict().lookup("dOuter") >> dOuter_;
    }
    else
    {
        FatalErrorInFunction
            << "injectionMethod must be either 'point' or 'disc'"
            << exit(FatalError);
    }
}


template<class CloudType>
void Foam::myConeInjection<CloudType>::setFlowType()
{
    const word flowType =
        this->coeffDict().template lookupOrDefault<word>
        (
            "flowType",
            word::null
        );

    if (flowType == "constantVelocity" || flowType == word::null)
    {
        flowType_ = ftConstantVelocity;

        Umag_.reset(this->coeffDict());
    }
    else if (flowType == "normalDistribution")
    {
        flowType_ = ftNormalDistribution;
        
        Umag_.reset(this->coeffDict());
        this->coeffDict().lookup("sigma_u") >> sigma_u_;        
    }
    else if (flowType == "pressureDrivenVelocity")
    {
        flowType_ = ftPressureDrivenVelocity;

        Pinj_.reset(this->coeffDict());
    }
    else if (flowType == "flowRateAndDischarge")
    {
        flowType_ = ftFlowRateAndDischarge;

        this->coeffDict().lookup("dInner") >> dInner_;
        this->coeffDict().lookup("dOuter") >> dOuter_;

        Cd_.reset(this->coeffDict());
    }
    else
    {
        FatalErrorInFunction
            << "flowType must be either 'constantVelocity', 'normalDistribution', "
            << "'pressureDrivenVelocity', or 'flowRateAndDischarge'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myConeInjection<CloudType>::myConeInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    injectionMethod_(imPoint),
    flowType_(ftConstantVelocity),
    position_
    (
        TimeFunction1<vector>
        (
            owner.db().time(),
            "position",
            this->coeffDict()
        )
    ),
    positionIsConstant_(isA<Function1s::Constant<vector>>(position_)),
    direction_
    (
        TimeFunction1<vector>
        (
            owner.db().time(),
            "direction",
            this->coeffDict()
        )
    ),
    injectorCell_(-1),
    injectorTetFace_(-1),
    injectorTetPt_(-1),
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
    thetaInner_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "thetaInner",
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "thetaOuter",
            this->coeffDict()
        )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"), owner.rndGen()
        )
    ),
    dInner_(vGreat),
    dOuter_(vGreat),
    Umag_(owner.db().time(), "Umag"),
    sigma_u_(0),    // should behave like a const velocity with this default value
    Cd_(owner.db().time(), "Cd"),
    Pinj_(owner.db().time(), "Pinj")
{
    duration_ = owner.db().time().userTimeToTime(duration_);

    setInjectionMethod();

    setFlowType();

    // Set total volume to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0, duration_);

    updateMesh();
}


template<class CloudType>
Foam::myConeInjection<CloudType>::myConeInjection
(
    const myConeInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    injectionMethod_(im.injectionMethod_),
    flowType_(im.flowType_),
    position_(im.position_),
    positionIsConstant_(im.positionIsConstant_),
    direction_(im.direction_),
    injectorCell_(im.injectorCell_),
    injectorTetFace_(im.injectorTetFace_),
    injectorTetPt_(im.injectorTetPt_),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    thetaInner_(im.thetaInner_),
    thetaOuter_(im.thetaOuter_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    dInner_(im.dInner_),
    dOuter_(im.dOuter_),
    Umag_(im.Umag_),
    sigma_u_(im.sigma_u_),
    Cd_(im.Cd_),
    Pinj_(im.Pinj_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myConeInjection<CloudType>::~myConeInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// turns out this reference would only be used by InjectionModel.H,
// would have to put this into InjectionModel.H as a virtual function for it to work correctly
// to be accessed by whatever is calling the separate injection types using InjectionModel.H
// so I made it, this seems to be the right method, but it didn't work as expected, so I'm leaving
// it here commented out in case it becomes useful in the future, don't want to lose the work 
// I did to figure out this format.
//template<class CloudType>
//inline Foam::distributionModel&
//Foam::myConeInjection<CloudType>::sizeDistribution()
//{
//    return sizeDistribution_;
//}


template<class CloudType>
void Foam::myConeInjection<CloudType>::updateMesh()
{
    if (injectionMethod_ == imPoint && positionIsConstant_)
    {
        vector position = position_.value(0);
        this->findCellAtPosition
        (
            injectorCell_,
            injectorTetFace_,
            injectorTetPt_,
            position
        );
    }
}


template<class CloudType>
Foam::scalar Foam::myConeInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::myConeInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (time0 >= 0 && time0 < duration_)
    {
        // Standard calculation
        return floor(parcelsPerSecond_*(time1 - time0));

        //// Modified calculation to make numbers exact
        //// causes me trouble with my cloudFunctionObjects
        //return floor(parcelsPerSecond_*time1 - this->parcelsAddedTotal());
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::myConeInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (time0 >= 0 && time0 < duration_)
    {
        return flowRateProfile_.integrate(time0, time1);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
void Foam::myConeInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    Random& rndGen = this->owner().rndGen();

    const scalar t = time - this->SOI_;

    switch (injectionMethod_)
    {
        case imPoint:
        {
            position = position_.value(t);
            if (positionIsConstant_)
            {
                cellOwner = injectorCell_;
                tetFacei = injectorTetFace_;
                tetPti = injectorTetPt_;
            }
            else
            {
                this->findCellAtPosition
                (
                    cellOwner,
                    tetFacei,
                    tetPti,
                    position,
                    false
                );
            }
            break;
        }
        case imDisc:
        {
            const scalar beta = twoPi*rndGen.globalScalar01();  // uniformly distributed random angle between 0 and 2 pi, random angle on a full circle
            const scalar frac = rndGen.globalScalar01();    // ends up acting like a uniformly distributed random distance between dInner_ and dOuter_
            const vector n = normalised(direction_.value(t));   // normalized (unit) vector of the input direction (angle of injection release)
            const vector t1 = normalised(perpendicular(n)); // normalized (unit) vector perpendicular to the input direction (angle of injection release)
            const vector t2 = normalised(n ^ t1);   // normalized (unit) vector that is the cross product of the normalized (unit) vector perpendicular to the input direction and the normalized (unit) input direction vector. This means that t2 is the normal to an imaginary plane formed between t1 and n/direction. It's almost like t1 is 90 deg from n/direction and t2 is 90 deg from t1 and n/direction but out of the plane. It seems to define a coord system of i,j,k x,y,z but rotated from the standard coord system by the angle/vector of n/direction. In my case, with n/direction as 45 deg release into the domain, n/direction is like the vector of release, t1 and t2 form the plane perpendicular to the release, with t1 going up and down in the x and z direction perpendicular to n, and t2 going to the left and right in the y direction.
            const vector tanVec = t1*cos(beta) + t2*sin(beta);  // the vector function for the edge of a circle, in the plane of t1 and t2, so the plane that is normal to input direction/n. So a uniformly distributed random location along the edge of a circle that is perpendicular/normal to the release angle.
            const scalar d = sqrt((1 - frac)*sqr(dInner_) + frac*sqr(dOuter_)); // uniformly distributed random spacing between dInner_ and dOuter_
            position = position_.value(t) + d/2*tanVec; // uniformly distributed position on a circle/disk perpendicular to the release direction. d acts as the uniformly distributed location between the origin of the circle/disk to the edge of the circle/disk (if dInner is not zero, from the inner edge of the donut er circle/disk to the outer edge dOuter), tanVec acts as the uniformly distributed direction to choose going out from the origin, the uniformly distributed position along the edge of the circle/disk. tanVec is setup to keep the circle/disk on a plane perpendicular/normal to the release direction, the input position determines where to place the origin of the circle/disk on that plane, dInner determines the inner edge of the circle/disk, dOuter determines the outer edge of the circle/disk.
            this->findCellAtPosition
            (
                cellOwner,
                tetFacei,
                tetPti,
                position,
                false
            );
            break;
        }
        default:
        {
            break;
        }
    }
}


template<class CloudType>
void Foam::myConeInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    Random& rndGen = this->owner().rndGen();

    const scalar t = time - this->SOI_;

    // Get the angle from the axis and the vector perpendicular from the axis.
    // If injecting at a point, then these are calculated from two new random
    // numbers. If a disc, then these calculations have already been done in
    // setPositionAndCell, so the angle and vector can be reverse engineered
    // from the position.
    scalar theta = vGreat;
    vector tanVec = vector::max;
    switch (injectionMethod_)
    {
        case imPoint:
        {
            const scalar beta = twoPi*rndGen.scalar01();  // uniformly distributed random angle between 0 and 2 pi, random angle on a full circle
            const scalar frac = rndGen.scalar01();    // ends up acting like a uniformly distributed random distance between thetaInner_ and thetaOuter_, separate from random numbers used for generating position stuff. A value between 0 and 1.
            const vector n = normalised(direction_.value(t));   // normalized (unit) vector of the input direction (angle of injection release)
            const vector t1 = normalised(perpendicular(n)); // normalized (unit) vector perpendicular to the input direction (angle of injection release)
            const vector t2 = normalised(n ^ t1);   // normalized (unit) vector that is the cross product of the normalized (unit) vector perpendicular to the input direction and the normalized (unit) input direction vector. This means that t2 is the normal to an imaginary plane formed between t1 and n/direction. It's almost like t1 is 90 deg from n/direction and t2 is 90 deg from t1 and n/direction but out of the plane. It seems to define a coord system of i,j,k x,y,z but rotated from the standard coord system by the angle/vector of n/direction. In my case, with n/direction as 45 deg release into the domain, n/direction is like the vector of release, t1 and t2 form the plane perpendicular to the release, with t1 going up and down in the x and z direction perpendicular to n, and t2 going to the left and right in the y direction.
            tanVec = t1*cos(beta) + t2*sin(beta);  // the vector function for the edge of a circle, in the plane of t1 and t2, so the plane that is normal to input direction/n. So a uniformly distributed random location along the edge of a circle that is perpendicular/normal to the release angle. But this time the circle is used as if it is farther away from the origin, a distance of the magnitude of direction*sin(theta), with the size of this circle changing proportionally by cos(theta), so the change in circle location/size determines the angle outward from direction. Technically theta varies tanVec and direction proportionally in such a way that theta is the cone angle.
            // up to here, same calcs repeated as for a disk in setPositionAndCell()
            theta =
                degToRad
                (
                    sqrt
                    (
                        (1 - frac)*sqr(thetaInner_.value(t))
                        + frac*sqr(thetaOuter_.value(t))
                    )
                );   // uniformly distributed random angle of value between thetaInner_ and thetaOuter_. Note that it forces negative thetaInner and thetaOuter to be treated as the same values but positive. Also note that this is NOT an angle difference, which would make it independent of location, thetaInner_ and thetaOuter_ could possibly not wrap around direction. It is also weird because this value acts more like 0 to 10 or 5 to 10, like a subsection of a sphere.
            break;
        }
        case imDisc:
        {
            const scalar r = mag(parcel.position() - position_.value(t));   // the distance between the point on the disk and the input position, so a distance along the disk
            const scalar frac = (2*r - dInner_)/(dOuter_ - dInner_);    // the fraction of distance along the disk for the current sampled point. Notice that this is NOT a separate random number from the random number distribution used for the original positions, so the cone angle ends up being distributed evenly as the random distribution of the positions are uniformly distributed
            tanVec = normalised(parcel.position() - position_.value(t));    // the vector function for the edge of a circle, in the plane of t1 and t2, so the plane that is normal to the input direction/n. But calculated from the already previously calculated positions. So it still holds the same randomized position as was already input into it before.
            theta =
                degToRad
                (
                    (1 - frac)*thetaInner_.value(t)
                    + frac*thetaOuter_.value(t)
                );  // an angle between thetaInner_ and thetaOuter_, but it is indirectly uniformly distributed, as it just uses the same uniformly distributed position information, it isn't built with a separate random number even though it probably should be
            break;
        }
        default:
        {
            break;
        }
    }

    // The direction of injection
    const vector dirVec =
        normalised
        (
            cos(theta)*normalised(direction_.value(t))
          + sin(theta)*tanVec
        );  // this is a WEIRD equation. The best I can tell, direction is the vector going out from the point, tanVec is the edge of the circle out at the end of the vector going out from the point, and the cos(theta) sin(theta) modify the magnitude of the direction vector and the circle size/position at the end of the vector proportionally. So it comes out as a cone, where the angle of the cone is determined by theta. It's like modifying the angle in the z direction, and proportionally modifying it in the y direction, the x direction changing size to go with it.

    // Set the velocity
    switch (flowType_)
    {
        // the final velocity is calculated the same each time, but these cases are various ways to calculate Umag_, or the initial speed of each vector.
        case ftConstantVelocity:
        {
            parcel.U() = Umag_.value(t)*dirVec;
            break;
        }
        case ftNormalDistribution:
        {
            // modify it with the normal distribution value
            // note, currently doesn't have any limits on the normal distribution, also no guards against negative values for U() if Umag_ is too small for the distribution
            parcel.U() = ( Umag_.value(t) + sigma_u_*rndGen.scalarNormal() )*dirVec;
            break;
        }
        case ftPressureDrivenVelocity:
        {
            const scalar pAmbient = this->owner().pAmbient();
            const scalar rho = parcel.rho();
            const scalar Umag = ::sqrt(2*(Pinj_.value(t) - pAmbient)/rho);
            parcel.U() = Umag*dirVec;
            break;
        }
        case ftFlowRateAndDischarge:
        {
            const scalar A = 0.25*pi*(sqr(dOuter_) - sqr(dInner_));
            const scalar massFlowRate =
                this->massTotal()*flowRateProfile_.value(t)/this->volumeTotal();    // note that volumeTotal_ and flowRateProfile_ are setup in such a way as to depend on inputs that aren't always consistent. nPar can be COMPLETELY separate from massTotal_, flowRateProfile_, and volumeTotal_, so the user would need to be careful to make sure all these get modified together correctly or there will be consistency issues.
            const scalar Umag =
                massFlowRate/(parcel.rho()*Cd_.value(t)*A);
            parcel.U() = Umag*dirVec;
            break;
        }
        default:
        {
            break;
        }
    }

    // Set the particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::myConeInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::myConeInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
