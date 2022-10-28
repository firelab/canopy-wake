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
    Foam::myCuboidInjection

Description
    Cuboid injection. Does uniform distribution between each set of values. If a min and max are the same value for 1 dimension, 
    this collapses down to be a rectangle. If for 2 dimensions, this collapses down to be a line. If for 3 dimensions, this 
    collapses down to be a point, which makes no sense to use but is possible. This is technically an nPar/volume method, but
    it preserves nPar/area if looking at it from one side. I originally thought that I needed to define a plane, so which dimension
    to use first for the distributions, but it turns out to not matter what dimension is calculated first.
    If eventually it is desired to change the distribution type for a given dimension, or to rotate the cuboid, will likely have to
    update this code or write a new version of code.

    User specifies:
      - Total mass to inject
      - Injection duration
      - Injection volume flow rate
      - Initial parcel velocity method, if constant, an initial parcel velocity
      - xMin, xMax, yMin, yMax, zMin, zMax

    Properties:
      - Parcel diameters obtained by distribution model
      - Parcels injected randomly across the patch
    
    example:
        
        (should I convert the format to a box like in toposet? I guess this still works for now)
        
        // do a uniform distribution in a cuboid between xMin, xMax, yMin, yMax, zMin, and zMax
        xMin = 35;
        xMax = 45;
        yMin = 50;
        yMax = 150;
        zMin = 0;
        zMax = 30;
        
        // do a uniform distribution along a yz plane between yMin, yMax, zMin, zMax, at position xMin/xMax
        xMin = 40;  // notice that xMin and xMax are equal, no other min and max are equal, collapsing the cuboid down to be a rectangle
        xMax = 40;  // notice that xMin and xMax are equal, no other min and max are equal, collapsing the cuboid down to be a rectangle
        yMin = 50;
        yMax = 150;
        zMin = 0;
        zMax = 30;
        

SourceFiles
    myCuboidInjection.C

\*---------------------------------------------------------------------------*/

#ifndef myCuboidInjection_H
#define myCuboidInjection_H

#include "InjectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class TimeFunction1;

class distributionModel;

/*---------------------------------------------------------------------------*\
                       Class myCuboidInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class myCuboidInjection
:
    public InjectionModel<CloudType>
{
public:
    
    //- Velocity type enumeration
    enum velocityType
    {
        vtConstantVelocity,
        vtCellCenterValues
    };
    
    
private:

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
        
        
        // release area coefficients
            
            //- minimum distribution location in the x dimension
            const scalar xMin_;
            
            //- maximum distribution location in the x dimension
            const scalar xMax_;
            
            //- minimum distribution location in the y dimension
            const scalar yMin_;
            
            //- maximum distribution location in the y dimension
            const scalar yMax_;
            
            //- minimum distribution location in the z dimension
            const scalar zMin_;
            
            //- maximum distribution location in the z dimension
            const scalar zMax_;
        
    
    // Private Member Functions
    
        //- Set the injection velocity type
        void setVelocityType();


public:

    //- Runtime type information
    TypeName("myCuboidInjection");


    // Constructors

        //- Construct from dictionary
        myCuboidInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        myCuboidInjection(const myCuboidInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new myCuboidInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~myCuboidInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

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
    #include "myCuboidInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
