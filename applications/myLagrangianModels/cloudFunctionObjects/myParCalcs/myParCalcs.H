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
    by preEvolve(). Vals are also reset to zero during postMove() if the list of particles/parcels is detected to 
    start over before a new call to preEvolve(), which happens at least once during the first time step because of 
    injections. Final single time calculations are done during postEvolve(). Values are written at calls to write() 
    which are done at each simulation write time.
    Because postMove() is called multiple times for a single particle/parcel before moving on to the next particle/parcel,
    the summation calculations are done in such a way that they only do the summations over the final timestep 
    a given particle/parcel moves.
    
    Calculated values include number density, concentration, volume fraction, and mass fraction.
    The internal cloud function versions of these variables are also called to output the cloud version set of these values.
    Time averaged and normalized versions of these values are also possible, but I decided not to implement them.
    It is also possible to do a per species or per injection form of each of these variables, I wanted to do this but in the end
    I decided it would be best to start simple with doing it over all the particles together 
     as total properties instead of as per species properties.
     Even that was a bit tough without accessing the individual injector information directly.
    
    Need to make sure the MassTotal_ variable is specified correctly for EACH separate injection 
    for the mass fraction output to make sense.
    
    
    Turns out that the calculations of massTotal_ and nParTotal_ is complicated, the injector information exists in the inherited
    cloud at constructor time, but the data for the injections are not filled until later. Means the timing of some calculations
    must change to compensate for some of the containers being empty at constructor time.
    

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
        
            //- particle count, the number of particles counted for a given eulerian grid cell
            autoPtr<volScalarField> dataPtr_numPar;
            
            //- number density, a concentration in terms of number of particles per cell volume instead of particle mass per cell volume
            autoPtr<volScalarField> dataPtr_numDensity;
            
            //- concentration, a concentration in terms of particle mass per cell volume instead of number of particles per cell volume
            autoPtr<volScalarField> dataPtr_conc;
            
            //- volume fraction, particle volume divided by cell volume
            autoPtr<volScalarField> dataPtr_volFrac;
            
            //- mass fraction, particle mass divided by the total mass of all particles to be injected
            // be aware, that unlike the variables used in the calculation of the numDensity, conc, and volFrac, the massFrac
            // depends on the individual injector values of massTotal_, which may not always be specified correctly by the user.
            // The massTotal_ variable for a given injector is independent of the rest of the injector information if the given 
            // injector specifies parcelBasisType_ as the value "fixed"
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
        
        // postMove() particle/parcel iteration variables and iteration temporary values
            
            //- the origId from the last iteration, used to tell when switching from one particle/parcel to another
            scalar pastIter_origId;
            
            //- temporary numPar summation value from the last iteration, added to the storage of the actual value when switching from one particle/parcel to another
            scalar pastIter_numPar;
            
            //- temporary numDensity summation value from the last iteration, added to the storage of the actual value when switching from one particle/parcel to another
            scalar pastIter_numDensity;
            
            //- temporary conc summation value from the last iteration, added to the storage of the actual value when switching from one particle/parcel to another
            scalar pastIter_conc;
            
            //- temporary volFrac summation value from the last iteration, added to the storage of the actual value when switching from one particle/parcel to another
            scalar pastIter_volFrac;
            
            //- temporary massFrac summation value from the last iteration, added to the storage of the actual value when switching from one particle/parcel to another
            scalar pastIter_massFrac;
        
    // Private Member Functions
        
        // mass calculation functions
        
            //- function for calculating the nParsTotal and massTotal from all the individual injectors.
            // This is kind of a complex function, probably start by just setting the desired vals to a default value
            void calc_totalVals();
        
        // other functions
            
            //- function for resetting the values of everything at new iteration loops
            // called during preEvolve() each time it is called, called during postMove() each time it is detected the 
            // iteration loop starts over before finishing postMove()
            void resetVars();
            
            
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
