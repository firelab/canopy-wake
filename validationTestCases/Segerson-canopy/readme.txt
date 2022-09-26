

This case is essentially an attempt to replicate the canopy case described in Segersson Dec2017 "A tutorial to urban wind flow using OpenFOAM", which is the canopy case described in Irvine and Gardiner and Hill Sept1997 "The Evolution of Turbulence Across a Forest Edge". The case was originally generated just using information from the Segersson paper, back before we had figured out the proper setup for the log profile Equilibrium Boundary Conditions for our cases when we were still using a domain average inlet wind. So it was originally an attempt to replicate Segersson, but with slightly different boundary conditions, numerics, our own set of models, and the case modified to be 3D instead of 2D so with slightly differing mesh resolutions. Later, it was found that Segersson had coded up the tutorial here: http://www.tfd.chalmers.se/~hani/kurser/OS_CFD/#YEAR_2017 which included some data to compare against and some files for plotting the profiles described in the paper.

The case has been modified to use the final log profile Equilibrium Boundary Conditions, solvers, and organization used in the final cases of the ember project. It still uses the original attempt mesh, so the 3D mesh compared to the Segersson 2D mesh, with the slightly differing mesh resolutions, but now using all the new solvers and methods for a log profile. In addition, the profile plotting scripts from the Segersson example code have been borrowed and modified for this case, which turns out to just be a slight update in the path for the new case organization, and adjusting the locations in the verticalProfiles sampling file. There was also a slight modification to get the inlet and outlet profiles right, for some odd reason these values were wrong in the original script.

This case offers the choice between a single LAD value at all heights, using the max LAD value provided in the Segersson paper, or to use an interpolated LAD profile, depending on the choice of setFields that is to be used. Calculations for converting the listed LAD values from the Segersson paper, from face centered values at the estimated heights to the desired cell center locations of the mesh for this case, are found in the extra interpolateLADtoCellCenterVals.xlsx file. Note the starting face centered locations are estimated assuming that Segersson used a non-graded mesh with dz 0.75, but looking at the Segersson code it looks like the mesh was technically graded, so this is not a perfect conversion from one set of LAD profile values to the other, but it should be close enough.

The case can also be run to test horizontal homogeneity by running with no porous media present. To do this, just use the same set of commands as are used for a normal run, but skipping the setFields command. This can be helpful for checking the Equilibrium Boundary Condition profiles. Note that this is NOT the method used in any of the standard cases of the ember project, a WindNinja no porous media set of models and solvers was normally used for non-porous simulations, with all commands calling myDragSimpleFoam being done as simpleFoam, and slight changes in the /0/ time directory files, the /constant/momentumTransport file, and other small changes in the /system/ folder. Should come out to be a similar result to do it this way with no calls to simpleFoam though.


So the only real difference from the Segersson coded version should be the mesh resolution, with 3D vs 2D, and the solvers and methods a bit, and slight differences in the choice of LAD profile from interpolation discrepancies. But one other difference is that the Equilibrium BCs set the inlet k value from the input velocity values, then the average inlet k value is used to calculate the inlet epsilon value from a specific formula given in the OpenFOAM documentation. But looking at the Segersson code, it looks like they used some set inlet values for k and epsilon that came out to be different from the ones used in this version of the case, that could be used instead of the standard calculated values. These inlet k and epsilon values may have come from measurements in the Irvine paper, not sure. So an option is to substitute the Segersson code values for inlet k and epsilon into the /0/ time directories of this case. The Segersson code value for inlet k is 1.3 m^2/s^2, the Segersson code value for inlet epsilon is 0.01 m^2/s^3. See the Segersson code to compare other differences between the two organizations of this canopy case.


No particle simulations involved with this test case, just the wind sims.





To run the case, use the following guide.


standard set of commands when in parallel:

  blockMesh 2>&1 | tee blockMesh.log
  checkMesh 2>&1 | tee checkMesh.log
  setFields 2>&1 | tee setFields.log
  decomposePar 2>&1 | tee decomposePar.log
  mpirun -np 3 myDragSimpleFoam -parallel 2>&1 | tee myDragSimpleFoam.log
  reconstructPar 2>&1 | tee reconstructPar.log

At the end after simulations, need to get the residuals:
  foamMonitor -l postProcessing/residuals/0/residuals.dat

At the end after simulations, need to do the sampling post processing for the python plots. 
Notice that I purposefully did NOT include the line "functions { #includeFunc verticalProfiles };" in controlDict,
as I just wanted data from the latest time, saves space on file since SS means only the last time matters for
plotting unless debugging. So I prefer post processing over run time processing for the plot sampling, to avoid
grabbing excess data. If you want data for all the timesteps, just drop the -latestTime argument.
command for without object cases
  myDragSimpleFoam -postProcess -func verticalProfiles -latestTime

To run the python scripts for the inlet outlet and mast vertical profile plots, in zz_plotFunctions run
  python plot_inlet_outlet_profiles.py
  python plot_vertical_profiles.py
If you are testing horizontal homogeneity, by running the solver without setFields, can still use the same 
commands as they are without replacing calls to myDragSimpleFoam with calls to simpleFoam.















