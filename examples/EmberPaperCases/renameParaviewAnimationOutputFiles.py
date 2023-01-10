#!/usr/bin/env python


"""
Utility for renaming paraview animation output files found in a given folder of a given case to the correct timestep names.

The user needs to input the path to the case and the name of the folder where the paraview animation output files are located.


The script requires that no paraview animation output files have been deleted yet for a given time, that it is the complete list for all the timestep folders found within the case. It also requires that only the paraview animation output for a single animation are found in the desired folder, no other files or even folders can be within the paraview animation output folder than the paraview animation output files. It also requires that the folder for the paraview animation output files be named what the animation output pictures should be named before the timestep is added on, example: /sideView_d/ turns files into sideView_d_timeXXX.png.

NOTE this script assumes that the plot output, the paraview animation output folder, is located as /casename/zz_pictures/. The user can also specify an intermediate path that goes any number of levels deeper than /zz_pictures/ to get to the animation output folder, the default is to assume no intermediate path.
"""


import os



### try to setup a set of arguments to use as input
### example usages of this script:
###     python renameParaviewAnimationOutputFiles.py /aa_originalConeInjection_singleTimeRelease/ /sideView_d/ --intermediatePath /particlePlots/
###     python renameParaviewAnimationOutputFiles.py aa_originalConeInjection_singleTimeRelease/ /sideView_d --intermediatePath /particlePlots
###     python renameParaviewAnimationOutputFiles.py aa_originalConeInjection_singleTimeRelease sideView_d --intermediatePath particlePlots
###     python renameParaviewAnimationOutputFiles.py /aa_originalConeInjection_singleTimeRelease/ /sideView_d/
###     python renameParaviewAnimationOutputFiles.py /aa_originalConeInjection_singleTimeRelease sideView_d/
###     python renameParaviewAnimationOutputFiles.py aa_originalConeInjection_singleTimeRelease sideView_d
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("casename", help="path to the case in which to run the script", type=str)
parser.add_argument("foldername", help="folder name of the paraview animation output files to be renamed, is also the name to be used for the output files. Renamed output files will take the form of foldername_timeXXXX.png", type=str)
parser.add_argument("--intermediatePath", help="optional intermediate path to the paraview animation output folder, to be used if the folder is at levels deeper than /casename/zz_pictures/", type=str) # the - or -- in front of the argument makes it an optional argument, no - or -- makes it a required argument
args = parser.parse_args()

#print args.casename     # print the user input casename from cmd line.
#print sys.argv  # print all the sys argument passed from cmd line including the program name.

casename = args.casename
foldername = args.foldername
if args.intermediatePath == None:
    intermediatePath = ''
else:
    intermediatePath = args.intermediatePath
#print('casename: \"' + casename + '\"')
#print('foldername: \"' + foldername + '\"')
#print('intermediatePath: \"' + intermediatePath + '\"')



### process the input paths, like strip any extra / symbols off the ends
### so that when combining paths, adding / symbols is consistent and doesn't
### result in extra / symbols in any created paths
casename = casename.strip('/')
foldername = foldername.strip('/')
intermediatePath = intermediatePath.strip('/')
#print('casename: \"' + casename + '\"')
#print('foldername: \"' + foldername + '\"')
#print('intermediatePath: \"' + intermediatePath + '\"')


### now setup the path to the paraview animation output
if intermediatePath == "":
    animationFileFolder = './' + casename + '/zz_pictures/' + foldername + '/'
else:
    animationFileFolder = './' + casename + '/zz_pictures/' + intermediatePath + '/' + foldername + '/'
print('animationFileFolder = \"' + animationFileFolder + '\"')


### would check the casename and animationFileFolder NOW to see if they exist and to warn
### if they are bad values, but it appears that good errors are already thrown at the 
### calls to os to get the files/folders from within those directories
### so I won't worry about it.



### function for using in a call to filter() from within a call to map()
### or for equivalent comprehension or generator expression or list comprehension
### ideas for this found at https://realpython.com/python-map-function/
### so it is designed to return true if it is to be a kept value
### and return false if it is to be an ignored value
def filterFunc_isFloat(number):
    try:
        test = float(number)
        return True
    except ValueError:
        return False


def main():
    
    
    ### get the list of input time directory information
    
    ### first, get the sorted list of time directories, including the final time
    times = sorted([ float(x) for x in filter(filterFunc_isFloat,os.listdir(casename)) ])    # adjusted to account for the non-int non-time directories, to only get the time directories
    nTimes = len(times)
    latestTime = max(times)
    #print(nTimes)
    #print(times)
    #print(latestTime)
    #print(os.listdir(casename))
    
    ### convert the list of time directories back to string
    times_asString = [ str(x) for x in times ]
    #print(times_asString)
    latestTime_asString = str(latestTime)
    #print(latestTime_asString)
    
    ### get rid of trailing .0's on strings that happens when converting them to floats then back to strings
    ### got this idea from https://stackoverflow.com/questions/2440692/formatting-floats-without-trailing-zeros
    times_asString = [ x.rstrip('0') for x in times_asString ]  # rstrip vs strip, strip does stripping to both outer edges of the string while rstrip just does stripping to the right outer edge
    times_asString = [ x.rstrip('.') for x in times_asString ]
    latestTime_asString = latestTime_asString.rstrip('0')   # oops!!! apparently the generator list and map do weird things if the thing to iterate over is just a single value, I got a list of 3 values from doing this particular line of code as a map or a generator expression instead of as a single variable expression as I am currently doing here
    latestTime_asString = latestTime_asString.rstrip('.')
    #print(times_asString)
    #print(latestTime_asString)
    
    ### get decimals turned into 'o' form of the times
    times_asString_oForm = [ x.replace(".","o") for x in times_asString ]
    #print(times_asString_oForm)
    latestTime_asString_oForm = times_asString_oForm[nTimes-1]  # equivalent without using nTimes is to just do listVals[-1] to get the last val of the list
    #print(latestTime_asString_oForm)
    
    ### finished getting the list of input time directory information
    
    
    ### find the animation output filenames
    ### might have to adjust this if there are other folders within the folder, there shouldn't be any other folders within the folder if the user follows the directions of the script though.
    oldAnimationFileNames = sorted( os.listdir(animationFileFolder) )  # need to sort it because it is in a weird order. See https://stackoverflow.com/questions/4813061/non-alphanumeric-list-order-from-os-listdir if this sorting method ever breaks down
    #print(oldAnimationFileNames)
    
    
    ### now rename all the animation output files
    oldAnimationFiles = [ animationFileFolder + x for x in oldAnimationFileNames ]
    #print(oldAnimationFiles)
    ####newAnimationFileNames = [ "sideView_active_yClip" + "_time" + x + ".png" for x in times_asString_oForm ]
    newAnimationFileNames = [ foldername + "_time" + x + ".png" for x in times_asString_oForm ]
    #print(newAnimationFileNames)
    newAnimationFiles = [ animationFileFolder + x for x in newAnimationFileNames ]
    #print(newAnimationFiles)
    for x in xrange(nTimes):
        os.rename(oldAnimationFiles[x],newAnimationFiles[x])
    
    print('finished renaming animation files')
    
    
if __name__ == '__main__':
    main()
