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

Application
    calcParPositionsFromCoords

Description
    For the chosen input times, converts the lagrangian barycentric coordinates 
    positions file values generated from running a kinematicCloud, to a 
    lagrangian x, y, z positions file called par_pos.

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "timeSelector.H"
#include "fvMesh.H"

//#include "Cloud.H"    // must be defined within passiveParticleCloud.H, doesn't seem to actually be needed
#include "passiveParticleCloud.H"


using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // create fields stuff
    
    const word cloudName = "kinematicCloud";

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "Scanning times to convert lagrangian barycentric coordinate positions to par_pos values for cloud " << cloudName
        << nl << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        
        Info<< "    Reading particle/parcel information into passiveParticleCloud" << endl;
        passiveParticleCloud myCloud(mesh, cloudName);
        Info<< "    Read information for " << returnReduce(myCloud.size(), sumOp<label>())
            << " particles/parcels" << endl;
        
        
        // make the storage for the output particle position information
        // got the ideas for the following lines of code from 
        //  /src/lagrangian/intermediate/submodels/CloudFunctionObjects/RelativeVelocity.C
        IOField<vector> par_pos
        (
            myCloud.fieldIOobject("par_pos", IOobject::NO_READ),
            myCloud.size()
        );
        
        Info<< "    calculating par_pos values" << endl;
        label i = 0;
        forAllConstIter(passiveParticleCloud, myCloud, iter)
        {
            par_pos[i] = iter().position();
            ++ i;
        }
        
        Info<< "    writing par_pos values" << endl;
        par_pos.write();
    }

    return 0;
}


// ************************************************************************* //
