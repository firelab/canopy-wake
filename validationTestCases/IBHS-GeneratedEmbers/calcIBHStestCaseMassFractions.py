#!/usr/bin/env python


"""
Utility for calculating and plotting the mass fractions for the 18 pans from the IBHS test case for a given case.

The user needs to comment out all case names except for the one of interest, and make sure the injection_yMid, nParsTotal, and par_rho are selected correctly for the case.

The script then reads in the lagrangian particle information, the origId, par_pos, and d, for the latest time. Then the script compares each particle location to the min and max edges of each pan to determine into which pan, if any, a given particle has landed. If a given particle lands in a particular pan, a particle count for that pan is incremented and the mass of the given particle is calculated from the particle diameter, the input particle density, and the volume of a particle equation, and the mass is then added to the summation of particle mass landed in that pan. After all pan information has been counted for all particles for each pan, the mass fractions are then calculated using just the total number of particles and amount of mass of particles that actually landed in the pans, the totals ignore all particles that did not land within any pans. The script then finally outputs the plot of the mass fractions to the /zz_pictures/particlePlots/ folder within the chosen case as zz_panMassFractions.png.


NOTE THAT TO RUN THIS PROPERLY, THE lagrangian positions FILE in the latest time needs to be converted to par_pos. Do this by running calcParPositionsFromCoords -latestTime or calcParPositionsFromCoords -time 3.9,4 in the directories of the cases of interest.
"""


##import os   ## used for "if not os.path.exists(filename): return None" in past scripts. Also for "times = map(int, os.listdir('../postProcessing/singleGraph'))" in other past scripts.
##import numpy as np  ## used for "x = np.zeros((npts,1))", "xx, zz = np.meshgrid(x_1D,z_1D)", "uxx = np.zeros((nz,nx))", "uMagg[k,i] = np.sqrt( (uxx[k,i] ** 2) + (uyy[k,i] ** 2) + (uzz[k,i] ** 2) )" and all kinds of other things and this style stuff in past scripts. Also for "U = np.genfromtxt(filename, dtype=np.float)" and "p_k_epsilon = np.genfromtxt(filename, dtype=np.float)" to read in super simple text data in other past scripts.
##from matplotlib import pyplot as plt    ## used for plotting

import os

import numpy as np
from matplotlib import pyplot as plt




### choose the case name
casename = './b_parSim'

### set the injection middle position location
### I was originally going to read in this value, but it looks like it will be a pain, so just set it here
### this means it needs changed if new cases have different geometry
injection_yMid = 4.572

### set the total number of particles released over the simulation
### I was originally going to read in this value, but it looks like it will be a pain, so just set it here
### this has to be manually hunted down from the simulation log files. Changes each time
### the injection number of parcels to release changes. Take care to remember to update this variable
#nParsTotal = 190000    # not actually used. Which is a good thing, because it turns out that this value changes for per time releases

### have the user specify the particle density
### if the particles have varying density, this no longer works
### I was lazy and didn't want to read it in
par_rho = 300   # kg/m^3


### function for using in a call to filter() from within a call to map()
### or for equivalent comprehension or generator expression
### ideas for this found at https://realpython.com/python-map-function/
### so it is designed to return true if it is to be a kept value
### and return false if it is to be an ignored value
def filterFunc_isFloat(number):
    try:
        test = float(number)
        return True
    except ValueError:
        return False


### function for reading in the most useful particle information for a given timestep
### if the data needs recycled, like in a for loop calling this function multiple times,
### that is up to the user of the function to handle the data cleanup
### useful info on reading files: https://stackoverflow.com/questions/31570237/reading-data-from-non-csv-files last lines, https://www.pythontutorial.net/python-basics/python-read-text-file/
### NOTICE THAT I'VE RENAMED THIS TO ORIGINAL, BECAUSE THIS TRIES TO USE THE positions FILE WHICH IS JUST coordinates
### AND IS NOT CORRECT. WOULD HAVE JUST COMMENTED THIS OUT, BUT DOING A NEW COPY OF THE FUNCTION WAS EASIER.
### I WANTED THIS KEPT FOR REFERENCE, I LEARNED A TON OF TECHNIQUES FOR READING IN WEIRD DATA, TRYING TO DEAL WITH
### THE ORIGINAL COORDINATES DATA.
def readParData_original(timestepString,casename):
    
    parDataFolder = casename + "/" + timestepString + "/lagrangian/kinematicCloud/"
    #parDataFolder = casename + "/3.9/lagrangian/kinematicCloud/"
    
    origID_file = parDataFolder + "origId"
    positions_file = parDataFolder + "positions"
    diameters_file = parDataFolder + "d"
    
    ### read nPars from origID file
    ### assume it is the same value for all the other files for this timestep (because it should be)
    ### notice that origID has the nPars at line 18, same for d, but positions starts it at line 19
    f = open(origID_file,'r')
    nPar = int(f.readlines()[18])   # notice it is just line 18, instead of calling the variable lines, I'm just immediately converting to int and using the desired variable as the storage
    f.close()
    #print(nPar)
    
    
    ### read the origID file
    f = open(origID_file,'r')
    lines = f.readlines()[20:20+nPar]   # notice I'm able to read all the lines in directly, still need to process them a bit though
    f.close()
    par_origID = [ int(x) for x in lines ]  # I still think that map() is cleaner
    #print(len(par_origID))
    #print(par_origID[0:10])
    #print(par_origID[-1])
    
    ### read the positions file
    ### this will require more processing
    f = open(positions_file,'r')
    lines = f.readlines()[19:19+nPar]   # notice I'm able to read all the lines in directly, still need to process them a bit though
    f.close()
    #print(len(lines))
    #print(lines[0:10])
    #print(lines[-1])
    ### now unpack the position data
    ### Looking at the /openfoam8/src/lagrangian/basic/particle/particleIO.C writePosition() file, the output is a barycentric coordinates_ variable followed by the particle celli, tetFacei, and tetPti. Looking at barycentric in /openfoam8/src/OpenFOAM/primitives/Barycentric/ files, it looks like a barycentric variable is just a vector storage with an extra 4th redundant component. Looking at the /openfoam8/src/lagrangian/basic/particle/particle.C correctAfterInteractionListReferral() function (other functions could be useful for this too, just this one is quickest and easiest to see), it looks like coordinates_.b() is x, coordinates_.c() is y, and coordinates_.d() is z, making coordinates_.a() the redundant component. Have to track a lot in this set of particle files, but it looks like coordinates_ is always generated with a 1 in the redundant .a() position, but in this correctAfterInteractionListReferral() function, the coordinate += means the .a() component changes as well, I did not see the .a() component used in any way that made it seem to be anything but the redundant coordinate.
    ### so it looks like just the 2nd, 3rd, and 4th variables from the (a b c d) celli tetFacei tetPti matter, so b c d.
    ### looks like if I replace all the "(" and ")" with nothing, and replace the endline characters "\n" with nothing, all the data stays in the same position but now just separated
    ### by spaces, and it looks like this still keeps these desired variables in the 2nd, 3rd, and 4th locations
    lines = [ x.replace("(","") for x in lines ]
    #print(lines[0:10])
    lines = [ x.replace(")","") for x in lines ]
    #print(lines[0:10])
    lines = [ x.replace("\n","") for x in lines ]
    #print(lines[0:10])
    lines = [ x.split(" ") for x in lines ]
    #print(lines[0:10][0:10])
    ### turns out that extracting columns of 2D lists is annoying and has to be done with a map style function
    par_x_string = [row[1] for row in lines]
    par_y_string = [row[2] for row in lines]
    par_z_string = [row[3] for row in lines]
    ### need another map style function to convert them to floats
    par_x = [ float(x) for x in par_x_string ]
    par_y = [ float(x) for x in par_y_string ]
    par_z = [ float(x) for x in par_z_string ]
    #print(par_x[0:10])
    #print(par_x[-1])
    #print(par_y[0:10])
    #print(par_y[-1])
    #print(par_z[0:10])
    #print(par_z[-1])
    ### TURNS OUT THAT WHILE THIS IS THE CORRECT WAY TO READ THIS STUFF IN, THE coordinates IS NOT WHAT I THOUGHT. These are not positions, but the whole
    ### set of values of the positions file is barycentric coordinates. You have to do an inner cross product between the original barycentric coordinate
    ### and a transform of those original values ( a Tensor cross vector ), where you do the transformation by looking up mesh position info from the mesh
    ### using the celli, tetFacei, and tetPti information to lookup the positions within the full mesh.
    ### As I didn't want to have to deal with reading in the full mesh, let alone doing cross product stuff, in this script, I just decided to convert the 
    ### values directly within the OpenFOAM code using an application, just output the position information as a new file, which then is simpler to read as well.
    
    ### read the diameters file
    f = open(diameters_file,'r')
    lines = f.readlines()[20:20+nPar]   # notice I'm able to read all the lines in directly, still need to process them a bit though
    f.close()
    par_d = [ float(x) for x in lines ]  # I still think that map() is cleaner
    #print(len(par_d))
    #print(par_d[0:10])
    #print(par_d[-1])
    
    return ( nPar,  par_origID,  par_x, par_y, par_z,  par_d )

### function for reading in the most useful particle information for a given timestep
### if the data needs recycled, like in a for loop calling this function multiple times,
### that is up to the user of the function to handle the data cleanup
### useful info on reading files: https://stackoverflow.com/questions/31570237/reading-data-from-non-csv-files last lines, https://www.pythontutorial.net/python-basics/python-read-text-file/
def readParData(timestepString,casename):
    
    parDataFolder = casename + "/" + timestepString + "/lagrangian/kinematicCloud/"
    #parDataFolder = casename + "/3.9/lagrangian/kinematicCloud/"
    
    origID_file = parDataFolder + "origId"
    positions_file = parDataFolder + "par_pos"
    diameters_file = parDataFolder + "d"
    
    ### read nPars from origID file
    ### assume it is the same value for all the other files for this timestep (because it should be)
    ### notice that origID has the nPars at line 18, same for d, but positions starts it at line 19
    f = open(origID_file,'r')
    nPar = int(f.readlines()[18])   # notice it is just line 18, instead of calling the variable lines, I'm just immediately converting to int and using the desired variable as the storage
    f.close()
    #print(nPar)
    
    
    ### read the origID file
    f = open(origID_file,'r')
    lines = f.readlines()[20:20+nPar]   # notice I'm able to read all the lines in directly, still need to process them a bit though
    f.close()
    par_origID = [ int(x) for x in lines ]  # I still think that map() is cleaner
    #print(len(par_origID))
    #print(par_origID[0:10])
    #print(par_origID[-1])
    
    ### read the par_pos file (NOT the positions file, that holds the positions in terms of barycentric coordinates)
    ### this will require more processing
    f = open(positions_file,'r')
    lines = f.readlines()[20:20+nPar]   # notice I'm able to read all the lines in directly, still need to process them a bit though
    f.close()
    #print(len(lines))
    #print(lines[0:10])
    #print(lines[-1])
    ### now unpack the position data
    ### looks like if I replace all the "(" and ")" with nothing, and replace the endline characters "\n" with nothing, all the data stays in the same position but now just separated
    ### by spaces, and it looks like this still keeps these desired variables in the 1st, 2nd, and 3rd locations
    lines = [ x.replace("(","") for x in lines ]
    #print(lines[0:10])
    lines = [ x.replace(")","") for x in lines ]
    #print(lines[0:10])
    lines = [ x.replace("\n","") for x in lines ]
    #print(lines[0:10])
    lines = [ x.split(" ") for x in lines ]
    #print(lines[0:10][0:10])
    ### turns out that extracting columns of 2D lists is annoying and has to be done with a map style function
    par_x_string = [row[0] for row in lines]
    par_y_string = [row[1] for row in lines]
    par_z_string = [row[2] for row in lines]
    ### need another map style function to convert them to floats
    par_x = [ float(x) for x in par_x_string ]
    par_y = [ float(x) for x in par_y_string ]
    par_z = [ float(x) for x in par_z_string ]
    #print(par_x[0:10])
    #print(par_x[-1])
    #print(par_y[0:10])
    #print(par_y[-1])
    #print(par_z[0:10])
    #print(par_z[-1])
    
    ### read the diameters file
    f = open(diameters_file,'r')
    lines = f.readlines()[20:20+nPar]   # notice I'm able to read all the lines in directly, still need to process them a bit though
    f.close()
    par_d = [ float(x) for x in lines ]  # I still think that map() is cleaner
    #print(len(par_d))
    #print(par_d[0:10])
    #print(par_d[-1])
    
    return ( nPar,  par_origID,  par_x, par_y, par_z,  par_d )


### setup the pan box locations
### centers of pans go from 40 ft to 74 ft every 2 ft for a total of 18 pans
### pans are 25.5" x 17.5" x 1", with 17.5" in the downwind direction, the direction of every 2 ft
### downwind direction is x in my simulations, pans go in a straight line directly downwind of the plume injection so y location is constant for each pan
### so output needs to be a list of pan x information, but single values for the pan y information
### convert everything to meters/SI units as quick as possible, output the information in terms of meters 
###  even though pan information was originally in inches and feet
### input is the injection y centerline value in meters, to be able to get the right y cell centers for the pans
###  would use the y center of the domain, but plumes could potentially be off to the side a bit (if this function
###  ever gets used for old cell zone injections)
def calc_panLocations(injection_yMid):
    
    inch_to_meter = 0.0254;
    ft_to_meter = 0.3048;
    
    nPans = 18
    pan_xSize = 17.5*inch_to_meter
    pan_ySize = 25.5*inch_to_meter
    ###pan_zSize = 1*inch_to_meter  ### not used
    
    firstPanXcenterLocation = 40*ft_to_meter
    xDistBetweenPanCenters = 2*ft_to_meter
    
    
    pan_halfXsize = pan_xSize/2
    pan_halfYsize = pan_ySize/2
    
    
    pan_xcenters = np.zeros((nPans,1))
    pan_xmins = np.zeros((nPans,1))
    pan_xmaxs = np.zeros((nPans,1))
    for i in xrange(nPans):
        pan_xcenters[i] = firstPanXcenterLocation + xDistBetweenPanCenters*i
        pan_xmins[i] = pan_xcenters[i] - pan_halfXsize
        pan_xmaxs[i] = pan_xcenters[i] + pan_halfXsize
    
    ### pan y information are single values
    ### doesn't change location in the y direction, pans go in a straight line in the x direction
    pans_ycenter = injection_yMid
    pans_ymin = pans_ycenter - pan_halfYsize
    pans_ymax = pans_ycenter + pan_halfYsize
    
    return ( nPans,  pan_xcenters, pan_xmins, pan_xmaxs, pans_ycenter, pans_ymin, pans_ymax )



### calculate the pan mass fractions
### start by reading in the particle information for the latest time directory
### then sum up the number of particles and mass per pan, calculating mass from the diameter and density of each particle
### in a given pan.
### then sum up the total number of particles and mass per pan to get the total number of particles and mass in all pans
### (ignores number of particles and mass outside the pans), and use that total number of particles and mass in all pans
### to calculate the mass fraction
### then delete the particle data, finished using it and it takes up extra space, doesn't need to hand around
def calcPanMassFractions( latestTime_asString, casename,  par_rho,  nPans,  pan_xmins, pan_xmaxs, pans_ymin, pans_ymax ):
    
    
    ### Use the particle information from the final time directory to bin up the particle information per pan
    ### also using the particle information per pan to calculate total values over all pans
    nPar,  par_origID,  par_x, par_y, par_z,  par_d = readParData( latestTime_asString, casename )
    
    
    ### initialize per pan sum totals to zero, to be incremented when particles are found to be in a given pan
    pan_nParTotals = np.zeros((nPans,1))
    pan_massTotals = np.zeros((nPans,1))
    
    ### for each particle, determine whether they are in any of the pans
    ### when found to be in a given pan, calculate values and add to the sums for that pan
    for parIdx in xrange(nPar):
        
        panIdx = 0
        while panIdx < nPans:
            
            #if parIdx == 0 or parIdx == 1 or parIdx == 2:
            #    print(parIdx, panIdx)
            
            ### use the >=, <= form when you want to count particles that land on the edges of the pan,
            ### use the >, < form when you want to count particles ignoring ones that land on the edges of the pan
            if par_x[parIdx] >= pan_xmins[panIdx] and par_x[parIdx] <= pan_xmaxs[panIdx] and par_y[parIdx] >= pans_ymin and par_y[parIdx] <= pans_ymax:
            #if par_x[parIdx] > pan_xmins[panIdx] and par_x[parIdx] < pan_xmaxs[panIdx] and par_y[parIdx] > pans_ymin and par_y[parIdx] < pans_ymax:
                
                ### add 1 to the particle count for the current pan
                pan_nParTotals[panIdx] += 1
                
                ### calculate the mass of the current particle using the current particle diameter and current particle rho, 
                ### adding it to the mass total for the current pan
                ### formula for mass is mass = rho*V = rho*4/3*pi*r^3 = rho*pi/6*d^3
                pan_massTotals[panIdx] += par_rho*np.pi/6*pow(par_d[parIdx],3)
                
                ### set panIdx to nPans to kick it out of the pan loop, but still moving on to the next particle in the particle loop
                ### I would use a for loop with break, but I believe that would kick it out of both loops, so I tried this method instead
                panIdx = nPans
            
            ### increment the panIdx. Notice that it gets set to a final value of nPans if the if statement ever gets entered
            panIdx += 1
            
    
    ### now calculate the total number of particles and total mass of particles across all pans
    #print(pan_xmins)
    #print(pan_xmaxs)
    #print(pan_nParTotals)
    #print(pan_massTotals)
    allPan_nParTotal = sum(pan_nParTotals)
    allPan_massTotal = sum(pan_massTotals)
    #print(allPan_nParTotal)
    #print(allPan_massTotal)
    
    
    ### now calculate the mass fractions for each pan using the mass of particles per pan and the total number of particles across all pans
    #pan_numDensities = map(lambda x: x/allPan_nParTotal, pan_nParTotals )
    #pan_massFractions = map(lambda x: x/allPan_massTotal, pan_massTotals )
    pan_numDensities = np.array(tuple([ x/allPan_nParTotal for x in pan_nParTotals ]))     #np.stack() also works to try to preserve the np.array() 2D array, but it is supposed to be slower, and an issue with certain versions of numpy
    pan_massFractions = np.array(tuple([ x/allPan_massTotal for x in pan_massTotals ]))
    #print(pan_numDensities)
    #print(pan_massFractions)
    #sum_pan_numDensities = sum(pan_numDensities)
    #sum_pan_massFractions = sum(pan_massFractions)
    #print(sum_pan_numDensities)     # should output a 1 or there is a problem
    #print(sum_pan_massFractions)    # should output a 1 or there is a problem
    
    
    return ( pan_nParTotals, pan_massTotals,  pan_numDensities, pan_massFractions,  allPan_nParTotal, allPan_massTotal )




def main():
    
    ### first, get the sorted list of time directories, including the final time
    times = sorted(   map(  float, filter( filterFunc_isFloat, os.listdir(casename) )  )   )     # adjusted to account for the non-int non-time directories, to only get the time directories
    #times = [ float(x) for x in filter(filterFunc_isFloat,os.listdir(casename)) ]    # attempt at equivalent comprehension expression form
    #times = ( float(x) for x in filter(filterFunc_isFloat,os.listdir(casename)) )    # attempt at equivalent generator expression form
    latestTime = max(times)
    #print(times)
    #print(latestTime)
    #print(os.listdir(casename))
    
    ### convert the list of time directories back to string
    times_asString = map(str,times)
    #times_asString = [ str(x) for x in times ]  # equivalent comprehension expression form
    #times_asString = ( str(x) for x in times )  # equivalent generator expression form
    #print(times_asString)
    latestTime_asString = str(latestTime)
    #print(latestTime_asString)
    
    ### get rid of trailing .0's on strings that happens when converting them to floats then back to strings
    #times_asString = map(lambda x: x.rstrip('0'), times_asString)   # I originally did this as a generator expression, but then had trouble with passing the variable into the function, it claimed to be a list. So I went back to this non-list type to see how that would go. Turns out that it seems to be happy either way, passing in a single variable not as a list was the trick. So it was separate from this variable, the problem could come back later of using non-list map vs list generator expression, I guess it means I might have to learn how to undo lists as that seems to be the more common expression for things and it is starting to be easier/clearer for me to do as well.
    #times_asString = map(lambda x: x.rstrip('.'), times_asString)
    times_asString = [ x.rstrip('0') for x in times_asString ]  # rstrip vs strip, strip does stripping to both outer edges of the string while rstrip just does stripping to the right outer edge
    times_asString = [ x.rstrip('.') for x in times_asString ]
    latestTime_asString = latestTime_asString.rstrip('0')   # oops!!! apparently the generator list and map do weird things if the thing to iterate over is just a single value, I got a list of 3 values from doing this particular line of code as a map or a generator expression instead of as a single variable expression as I am currently doing here
    latestTime_asString = latestTime_asString.rstrip('.')
    #print(times_asString)
    #print(latestTime_asString)
    
    ### get decimals turned into 'o' form of the times
    times_asString_oForm = map( lambda s: s.replace(".","o"),times_asString )
    #print(times_asString_oForm)
    times_asString_oForm_list = list(times_asString_oForm)  # supposedly map returns values not as a list but as something else, where len() isn't supposed to work on it
    #print(times_asString_oForm)
    nTimes = len(times_asString_oForm_list) # tried it, and apparently even though map shouldn't be a list, so len() shouldn't work on it, len() ALSO works on the original map output without turning it into a list. Guess this converting to list was redundant and not necessary
    #print(nTimes)
    latestTime_asString_oForm = times_asString_oForm[nTimes-1]  # equivalent without using nTimes is to just do listVals[-1] to get the last val of the list
    #print(latestTime_asString_oForm)
    
    ### finished playing with https://realpython.com/python-map-function/ methods
    
    
    
    ### I decided to have the user select the injection_yMid or domain_yMid, essentially selecting the pans_ycenter location in the domain
    ### rather than reading in the domain or the injection information. If it were not an input, this would be the time to read in the 
    ### domain and injection information
    
    
    
    ### calculate the pan sizes, the list of box information to be used for finding the locations of particles for summing up particles into boxes
    nPans,  pan_xcenters, pan_xmins, pan_xmaxs, pans_ycenter, pans_ymin, pans_ymax = calc_panLocations(injection_yMid)
    #print(nPans)
    #print(pan_xcenters)
    #print(pan_xmins)
    #print(pan_xmaxs)
    #print(pans_ycenter)
    #print(pans_ymin)
    #print(pans_ymax)
    
    
    
    ### This function is just for getting the mass fractions based off of particles that actually land in the bins
    ### so don't need to calculate total mass for all particles released during the simulation, just the total mass
    ### of particles found within the bins. If it needed to be total mass, the more complex function to calculate said
    ### values would need to be implemented here.
    ### I think the user should still specify the total number of particles to be released, as found manually from the log
    ### files. If a mass were needed, would have to do diameters each time particles were released and it gets more complicated
    ### and values can't just be input, I would need to go through counting up particle information at each simulation time.
    
    
    
    ### Use the particle information from the final time directory to bin up the particle information per pan
    ### also using the particle information per pan to calculate total values over all pans
    ### using these values to calculate the mass fractions per pan
    ### most of these things are handled/done within calcPanMassFractions, which calls other functions as needed
    ### to do these calculations
    pan_nParTotals, pan_massTotals,  pan_numDensities, pan_massFractions,  allPan_nParTotal, allPan_massTotal = calcPanMassFractions( latestTime_asString, casename,  par_rho,  nPans,  pan_xmins, pan_xmaxs, pans_ymin, pans_ymax )
    
    
    ### now plot the nParTotals and massTotals information
    plt.rc('figure', figsize=(18, 10))  # set the figure size (in inches)
    fig1, (ax1, ax2) = plt.subplots(1, 2)
    
    ax1.plot(pan_xcenters, pan_nParTotals, 'D', label='nParTotals', clip_on=False)
    ax1.minorticks_on()
    ax1.legend(loc='best', shadow=True)
    ax1.set_xlabel('distance downwind, x [m]')
    ax1.set_ylabel('nPars [-]')
    ax1.set_title('pan nParTotals')
    ax1.grid(True)
    
    ax2.plot(pan_xcenters, pan_massTotals, 'D', label='massTotals', clip_on=False)
    ax2.minorticks_on()
    ax2.legend(loc='best', shadow=True)
    ax2.set_xlabel('distance downwind, x [m]')
    ax2.set_ylabel('mass [kg]')
    ax2.set_title('pan massTotals')
    ax2.grid(b=True)
    
    fig1.tight_layout(rect=[0.02,0.02,0.98, 0.98])
    
    plt.savefig(casename + '/zz_pictures/zz_panTotals.png')
    
    plt.show()
    
    
    ### now plot the mass fraction information
    plt.rc('figure', figsize=(18, 10))  # set the figure size (in inches)
    fig2, (ax1, ax2) = plt.subplots(1, 2)
    
    ax1.plot(pan_xcenters, pan_numDensities, 'D', label='numberDensities', clip_on=False)
    ax1.minorticks_on()
    ax1.legend(loc='best', shadow=True)
    ax1.set_xlabel('distance downwind, x [m]')
    ax1.set_ylabel('numDensity [-]')
    ax1.set_title('pan number densities')
    ax1.grid(True)
    
    ax2.plot(pan_xcenters, pan_massFractions, 'D', label='massFractions', clip_on=False)
    ax2.minorticks_on()
    ax2.legend(loc='best', shadow=True)
    ax2.set_xlabel('distance downwind, x [m]')
    ax2.set_ylabel('massFraction [-]')
    ax2.set_title('pan mass fractions')
    ax2.grid(b=True)
    
    fig2.tight_layout(rect=[0.02,0.02,0.98, 0.98])
    
    plt.savefig(casename + '/zz_pictures/zz_panMassFractions.png')
    
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    main()
