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
    Foam::myConeInjection

Description
    This injector injects particles in a number of cones. The user specifies a
    position and a direction to inject at, and two angles to inject between.
    Optionally, this injector can introduce particles over a disc, instead of
    at a point, in which case inner and outer diameters of the disc are also
    specified.

    The velocity is specified either as constant, or it is calculated from an
    injection pressure, or it is calculated from the injector mass flow rate
    and a discharge coefficient; i.e.:

    Constant velocity:
    \f[
    U = U_{constant}
    \f]

    Pressure driven velocity:
    \f[
    U = \sqrt{2(p_{injection} - p)/\rho}
    \f]

    Flow rate and discharge:
    \f[
    U = \dot{m}/(\rho A C_{discharge})
    \f]
    
    I have modified this function to make it work for myParCalcs. Specifically, I
    made parcelsToInject() output the right number of particles/parcels when inputting
    timeStart() and timeEnd(), even if simulations start out at different times
    than the first time. I was also making an accessor function for the sizeDistribution
    to allow myParCalcs to be able to generate a list of particle/parcel diameters,
    but it turns out that stopped being required because it required a ton more stuff
    to be accessible. So I left that addition, but commented out, just in case it could
    be useful for the future.
    
    I also modified it to have an additional flow type, a normal distribution, to allow 
    modifying Umag by a normal distribution.

Usage
    \table
    Property        | Description                    | Required | Default
    position        | The injection position         | yes      |
    direction       | The injection direction        | yes      |
    thetaInner      | The inner cone angle           | yes      |
    thetaOuter      | The outer cone angle           | yes      |
    injectionMethod | Inject at a point or on a disc | no       | point
    dInner          | The inner disc diameter        |\\
                                               if disc or flowRateAndDischarge |
    dInner          | The outer disc diameter        |\\
                                               if disc or flowRateAndDischarge |
    flowType        | Inject with constantVelocity, pressureDrivenVelocity \\
                                 or flowRateAndDischarge | no | constantVelocity
    Umag            | The injection velocity         | if constantVelocity |
    sigma_u         | The standard deviation of velocity | if normalDistribution |
    Pinj            | The injection pressure         |\\
                                                     if pressureDrivenVelocity |
    Cd              | The discharge coefficient      | if flowRateAndDischarge |
    \endtable

    Example specification:

    \verbatim
    injectionModels
    {
        model1
        {
            type            myConeInjection;

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
            thetaInner      0;
            thetaOuter      45;

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

            //// - Or, inject with flow rate and discharge coefficient
            ////   This also requires massTotal, dInner and dOuter
            //flowType        flowRateAndDischarge;
            //Cd              0.9;

            //// - Or, inject at a pressure
            //flowType        pressureDrivenVelocity;
            //Pinj            10e5;
        }
    }
    \endverbatim

SourceFiles
    myConeInjection.C

\*---------------------------------------------------------------------------*/

#ifndef myConeInjection_H
#define myConeInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "TimeFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class myConeInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class myConeInjection
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
        ftNormalDistribution,
        ftPressureDrivenVelocity,
        ftFlowRateAndDischarge
    };


private:

    // Private Data

        //- Point/disc injection method
        injectionMethod injectionMethod_;

        //- Flow type
        flowType flowType_;

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
            scalar sigma_u_;

            //- Discharge coefficient []
            TimeFunction1<scalar> Cd_;

            //- Injection pressure [Pa]
            TimeFunction1<scalar> Pinj_;


    // Private Member Functions

        //- Set the injection type
        void setInjectionMethod();

        //- Set the injection flow type
        void setFlowType();


public:

    //- Runtime type information
    TypeName("myConeInjection");


    // Constructors

        //- Construct from dictionary
        myConeInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        myConeInjection(const myConeInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new myConeInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~myConeInjection();


    // Member Functions
        
        // turns out this reference would only be used by InjectionModel.H,
        // would have to put this into InjectionModel.H as a virtual function for it to work correctly
        // to be accessed by whatever is calling the separate injection types using InjectionModel.H
        // so I made it, this seems to be the right method, but it didn't work as expected, so I'm leaving
        // it here commented out in case it becomes useful in the future, don't want to lose the work 
        // I did to figure out this format.
        ////- Return reference to the sizeDistribution
        //inline distributionModel& sizeDistribution();

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
    #include "myConeInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
