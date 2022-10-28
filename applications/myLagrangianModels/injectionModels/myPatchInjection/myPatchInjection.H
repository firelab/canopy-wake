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
    Foam::myPatchInjection

Description
    Patch injection.

    User specifies:
      - Total mass to inject
      - Name of patch
      - Injection duration
      - Injection volume flow rate
      - Initial parcel velocity method, if constant, an initial parcel velocity

    Properties:
      - Parcel diameters obtained by distribution model
      - Parcels injected randomly across the patch

SourceFiles
    myPatchInjection.C

\*---------------------------------------------------------------------------*/

#ifndef myPatchInjection_H
#define myPatchInjection_H

#include "InjectionModel.H"
#include "patchInjectionBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class TimeFunction1;

class distributionModel;

/*---------------------------------------------------------------------------*\
                       Class myPatchInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class myPatchInjection
:
    public InjectionModel<CloudType>,
    public patchInjectionBase
{
    
    //- Velocity type enumeration
    enum velocityType
    {
        vtConstantVelocity,
        vtCellCenterValues
    };
    
    // Private Data

        //- Injection duration [s]
        scalar duration_;

        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Flow rate profile relative to SOI []
        const TimeFunction1<scalar> flowRateProfile_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;
        
        
        //- Velocity type
        velocityType velocityType_;
        
        // Velocity type coefficients
            
            //- Initial parcel velocity [m/s] for when velocityType is constant velocity (single velocity value for all parcels)
            // dang, using a const value means it HAS to be assigned in the initializer list, but I'm trying to only construct this AFTER
            // the initializer list ONLY if needed. So I had to drop the const from the type.
            // Looks like this is consistent with coneInjection though, so that at least is good. Wish there were a better way
            vector U0_;
        
    // Private Member Functions
    
        //- Set the injection velocity type
        void setVelocityType();


public:

    //- Runtime type information
    TypeName("myPatchInjection");


    // Constructors

        //- Construct from dictionary
        myPatchInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        myPatchInjection(const myPatchInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new myPatchInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~myPatchInjection();


    // Member Functions

        //- Inherit updateMesh from patchInjectionBase
        using patchInjectionBase::updateMesh;

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Inherit setPositionAndCell from patchInjectionBase
            using patchInjectionBase::setPositionAndCell;

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
    #include "myPatchInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
