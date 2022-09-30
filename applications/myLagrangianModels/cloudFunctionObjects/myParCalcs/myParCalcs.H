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

Class
    Foam::myParCalcs

Description
    Creates particle calculation fields on the carrier/eulerian phase
    Also outputs a few of the cloud values that exist in the code but that are usually ignored
    
    Summation calculations are performed during postMove() where vals are reset to zero each full set of iterations
    by preEvolve(). Final single time calculations are done during postEvolve(). Values are written at calls 
    to write() which are done at each simulation write time.
    
    Calculated values include number density, concentration, volume fraction, and mass fraction.
    The internal cloud function versions of these variables are also called to output the cloud version set of these values.
    Time averaged and normalized versions of these values are also possible, but I decided not to implement them.
    It is also possible to do a per species or per injection form of each of these variables, I wanted to do this but in the end
    I decided it would be best to start simple with doing it over all the particles together 
     as total properties instead of per species properties.
     Even that was a bit tough without accessing the individual injector information directly.
    
    Need to make sure the injection function mass and volume calculations make sense 
    or all the outputs of this function will be skewed.
    

SourceFiles
    myParCalcs.C

\*---------------------------------------------------------------------------*/

#ifndef myParCalcs_H
#define myParCalcs_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class myParCalcs Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class myParCalcs
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;

        
        // calculated val data pointers
        
            //- number density, a concentration in terms of number of particles per cell volume instead of particle mass per cell volume
            autoPtr<volScalarField> dataPtr_numDensity;
            
            //- concentration, a concentration in terms of particle mass per cell volume instead of number of particles per cell volume
            autoPtr<volScalarField> dataPtr_conc;
            
            //- volume fraction, particle volume divided by cell volume
            autoPtr<volScalarField> dataPtr_volFrac;
            
            //- mass fraction, particle mass divided by the total mass of all particles to be injected
            autoPtr<volScalarField> dataPtr_massFrac;
        
        // cloud data pointers
            
            //- variable returned from calls to cloud theta() function, the cloud particle volume fraction field
            autoPtr<volScalarField> dataPtr_cloud_theta;
            
            //- variable returned from calls to cloud alpha() function, the cloud particle mass fraction field
            autoPtr<volScalarField> dataPtr_cloud_alpha;
            
            //- variable returned from calls to cloud rhoEff() function, the cloud particle effective density field (looks more like a concentration to me)
            autoPtr<volScalarField> dataPtr_cloud_rhoEff;
        
        // mass calculation variables
            
            //- the calculated total number of particles/parcels to be released by all injections over the time of the entire simulation
            scalar nParTotal_;
            
            //- the calculated total mass to be released by all injections over the time of the entire simulation
            scalar massTotal_;
        
    // Private Member Functions
        
        // mass calculation functions
        
            //- function for calculating the nParsTotal and massTotal from all the individual injectors.
            // This is kind of a complex function, probably start by just setting the desired vals to a default value
            void calc_totalVals();
            
            
protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();
        
        //- function for initializing data for data pointers, called at constructor time
        // this function holds information useful to understanding where write() will write data
        void initializeDataPointers();


public:

    //- Runtime type information
    TypeName("myParCalcs");


    // Constructors

        //- Construct from dictionary
        myParCalcs
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        myParCalcs(const myParCalcs<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new myParCalcs<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~myParCalcs();


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-move hook
            virtual void postMove
            (
                parcelType& p,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );
            
            //- Post-evolve hook
            virtual void postEvolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myParCalcs.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
