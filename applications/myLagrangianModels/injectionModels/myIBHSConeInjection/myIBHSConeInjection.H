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
    Foam::myIBHSConeInjection

Description
    This injector injects particles in a number of cones. The user specifies a
    position and a direction to inject at, and two angles to inject between.
    Optionally, this injector can introduce particles over a disc, instead of
    at a point, in which case inner and outer diameters of the disc are also
    specified.

    The velocity is specified either as constant, or it is calculated from a
    normal distribution. Note that the normal distribution of velocity won't work
    well for very low velocities, as it doesn't have a method to truncate the distribution
    or to stop negative values from switching the direction of the distribution.
    
    I have modified this to now do a uniform distribution or a normal distribution of theta,
    the half cone angle. The user specifies thetaDistributionType to be "uniform", or "normal",
    where "uniform" is uses a uniform distribution between thetaInner and thetaOuter, and
    "normal" does a normal distribution using a new variable sigma_theta instead of using 
    thetaInner and thetaOuter.

Usage
    \table
    Property        | Description                    | Required | Default
    position        | The injection position         | yes      |
    direction       | The injection direction        | yes      |
    thetaInner      | The inner cone angle           | yes      |           ignored if thetaDistributionType is "normal"
    thetaOuter      | The outer cone angle           | yes      |           ignored if thetaDistributionType is "normal"
    injectionMethod | Inject at a point or on a disc | no       | point
    dInner          | The inner disc diameter        |\\
                                               if disc or flowRateAndDischarge |
    dInner          | The outer disc diameter        |\\
                                               if disc or flowRateAndDischarge |
    flowType        | Inject with constantVelocity \\
                               or normalDistribution | no       | constantVelocity
    Umag            | The injection velocity         | if constantVelocity or normalDistribution |
    sigma_u         | The standard deviation of velocity | if normalDistribution |
    
    thetaDistributionType | type of distribution to use for theta | yes |
    sigma_theta     | the variance of the half cone angle theta | if thetaDistributionType is "normal" |
    \endtable

    Example specification:

    \verbatim
    injectionModels
    {
        model1
        {
            type            myIBHSConeInjection;

            // Times
            SOI             0;
            duration        1;

            // Quantities
            massTotal       0; // <-- not used with these settings
            parcelBasisType fixed;
            parcelsPerSecond 1000000;
            flowRateProfile constant 1;
            nParticle       1;

            // Sizes
            sizeDistribution
            {
                type        fixedValue;
                fixedValueDistribution
                {
                    value   0.0025;
                }
            }

            // Geometry
            position        (-0.15 -0.1 0);
            direction       (1 0 0);
            thetaInner      0;  // not used if thetaDistributionType is "normal"
            thetaOuter      45;  // not used if thetaDistributionType is "normal"

            // - Inject at a point
            injectionMethod point;

            //// - Or, inject over a disc:
            //injectionMethod disc;
            //dInner          0;
            //dOuter          0.05;

            // Velocity

            // - Inject with constant velocity
            flowType        constantVelocity;
            Umag            1;

            //// - Or, inject with a normal distribution of velocity speed
            ////   just uses a constantVelocity input, followed by a standard deviation of that velocity
            //flowType      normalDistribution;
            //Umag          1;  // takes 1 m/s as the mean
            //sigma_u       2;  // uses 2 m/s as the standard distribution. Assumes Umag large enough to avoid shifting into the negative velocities as no limit on the standard distribution
            
            // Injection half angle theta distribution
            
            // - Use a uniform distribution
            thetaDistributionType       uniform;
            
            //// - Or, use a normal distribution, requires specifying the variance of the injection half cone angle theta,
            //// causes thetaInner and thetaOuter to be ignored
            //thetaDistributionType      normal;
            //sigma_theta                    10;
        }
    }
    \endverbatim

SourceFiles
    myIBHSConeInjection.C

\*---------------------------------------------------------------------------*/

#ifndef myIBHSConeInjection_H
#define myIBHSConeInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "TimeFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class myIBHSConeInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class myIBHSConeInjection
:
    public InjectionModel<CloudType>
{
public:

    //- Injection method enumeration
    enum injectionMethod
    {
        imPoint,
        imDisc
    };

    //- Flow type enumeration
    enum flowType
    {
        ftConstantVelocity,
        ftNormalDistribution
    };
    
    // half cone angle theta distribution type enumeration
    enum thetaDistributionType
    {
        uniform,
        normal
    };


private:

    // Private Data

        //- Point/disc injection method
        injectionMethod injectionMethod_;

        //- Flow type
        flowType flowType_;
        
        //- half cone angle theta distribution type
        thetaDistributionType thetaDistributionType_;

        //- Position of the injector
        const TimeFunction1<vector> position_;

        //- Is the position constant?
        const bool positionIsConstant_;

        //- Centreline direction in which to inject
        const TimeFunction1<vector> direction_;

        //- Cell label corresponding to the injector position
        label injectorCell_;

        //- Tet-face label corresponding to the injector position
        label injectorTetFace_;

        //- Tet-point label corresponding to the injector position
        label injectorTetPt_;

        //- Injection duration [s]
        scalar duration_;

        //- Number of parcels to introduce per second
        const label parcelsPerSecond_;

        //- Flow rate profile relative to SOI []
        const TimeFunction1<scalar> flowRateProfile_;

        //- Inner half-cone angle relative to SOI [deg]
        const TimeFunction1<scalar> thetaInner_;

        //- Outer half-cone angle relative to SOI [deg]
        const TimeFunction1<scalar> thetaOuter_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;


        // Disc geometry

            //- The inner disc diameter [m]
            scalar dInner_;

            //- The outer disc diameter [m]
            scalar dOuter_;


        // Velocity model coefficients

            //- Parcel velocity [m/s]
            TimeFunction1<scalar> Umag_;

            //- Standard deviation of velocity [m/s]
            // set when flowType is normalDistribution
            scalar sigma_u_;
        
        // theta distribution type variables
            
            //- the variance of the half cone angle theta
            // set when thetaDistributionType is set to normal
            // when this is used, thetaInner and thetaOuter are ignored
            scalar sigma_theta_;


    // Private Member Functions

        //- Set the injection type
        void setInjectionMethod();

        //- Set the injection flow type
        void setFlowType();
        
        //- set the theta distribution type
        void setThetaDistributionType();


public:

    //- Runtime type information
    TypeName("myIBHSConeInjection");


    // Constructors

        //- Construct from dictionary
        myIBHSConeInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        myIBHSConeInjection(const myIBHSConeInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new myIBHSConeInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~myIBHSConeInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);



        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myIBHSConeInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
