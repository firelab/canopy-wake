

This case is essentially an attempt to replicate the Equilibrium Boundary Conditions test of figure 2 in the WagenbrennerEtAl_2019 paper titled "Development and Evaluation of a ReynoldsAveraged Navier–Stokes Solver in WindNinja for Operational Wildland Fire Applications". In the original test cases for the ember project, this was the test case from /1stTestWindNinjaLibso/ where we were first developing the WindNinjaLibso library for OpenFOAM 8, to be used before adding porous media to the simulations. During that process, there were lots of side attempts, first with the inlet BCs, then adding in a new wall epsilon function to replicate the epsilon wall function from OpenFOAM 2.2.0 used by the original WindNinja because we thought the base epsilon wall function might be different enough in OpenFOAM 8 (might not be, but we decided it was different enough to replicate the OpenFOAM 2.2.0 version of it into OpenFOAM 8), then adding in a new wall nut function to replicate the one used in OpenFOAM 2.2.0 (this one was much less likely to be different, but at the time we decided it was different enough so we went ahead and started using the replicated one), then going through various iterations adding and updating a k-epsilon model that was supposed to replicate the k-epsilon model used in OpenFOAM 2.2.0 as we thought maybe this had changed quite a bit as well between versions of OpenFOAM and we knew the OpenFOAM 2.2.0 version was already vetted properly by WindNinja (turns out that it probably wasn't that different, but at the time we decided it was different enough so we went ahead and started using the replicated one). During the tests to add in and update all these libraries, the first attempt had major problems because the first cell height was accidentally set to zero which led to us discovering that if z is less than z0 that there is backwards flow at the inlet. Then we tried using an updated version of applyInit but as there were problems trying to get the velocity to propagate downwind it was decided to set the initial velocity field to the free stream velocity to help make it easier to see whether changes in the velocity field were from velocity propogating rather than just being stuck at the start values. Things went much smoother after that, as instead of trying to start out with WindNinja stuff and figure out what was breaking, we started out from a known working example with domain average winds instead of the log profile and tried to track down what was different to get it to work with WindNinja stuff. The resulting case (after the first few git updates to show the various most important versions of the code during this process of upgrading the code) is the end result of all this hard work to try to get WindNinja log profile equilibrium boundary conditions to work in OpenFOAM 8, with the organization of the files updated to the latest style (not really much different, got it worked out leading up to this particular case) and the wind speeds updated to more closely match that of the WagenbrennerEtAl_2019 paper.

It should be noted that for this and all other cases using the WindNinjaLibso, the BCs should match the original WindNinja BCs pretty well, though maybe the version of the BCs could be walked back a step to use the built in ones of OpenFOAM 8 in a few cases, but the numerics (fvSchemes, fvSolution) used for all these simulations are NOT the same as those used in WindNinja, they ended up coming from an attempt to switch from the OpenFOAM 8 windAroundBuildings tutorial numerics to the WindNinja numerics, with the end result being that a mix of the two sets of numerics were used for this test case and for all the test cases in this particular porous media project.


No particle simulations involved with this test case, just the wind sims.





To run the case, use the following guide.


standard set of commands when in parallel:

  blockMesh 2>&1 | tee blockMesh.log
  checkMesh 2>&1 | tee checkMesh.log
  decomposePar 2>&1 | tee decomposePar.log
  mpirun -np 3 simpleFoam -parallel 2>&1 | tee simpleFoam.log
  reconstructPar 2>&1 | tee reconstructPar.log

At the end after simulations, need to get the residuals:
  foamMonitor -l postProcessing/residuals/0/residuals.dat

At the end after simulations, need to do the sampling post processing for the python plots. 
Notice that I purposefully did NOT include the line "functions { #includeFunc singleGraph };" in controlDict,
as I just wanted data from the latest time, saves space on file since SS means only the last time matters for
plotting unless debugging. So I prefer post processing over run time processing for the plot sampling, to avoid
grabbing excess data. If you want data for all the timesteps, just drop the -latestTime argument.
command for without object cases
  simpleFoam -postProcess -func singleGraph -latestTime

To run the python script for the vertical profile plots, inlet vs outlet and some downwind profiles, in the folder with plot_vertical_profiles.py, run
  python plot_verticalProfiles.py
  python plot_inletOutletProfiles.py
You may need to adjust the singleGraph and python plot files for separate test cases.















