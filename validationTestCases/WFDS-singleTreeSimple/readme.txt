

This case is essentially an attempt to replicate the with and without trunk tree cases described in the masters thesis Mueller 2012 "LES Modeling of Flow through Vegetation with Applications to Wildland Fires". The case was originally generated back before we had figured out the proper setup for the log profile Equilibrium Boundary Conditions for our cases when we were still using a domain average inlet wind. So it was originally an attempt to replicate Mueller, but with slightly different boundary conditions, numerics, our own set of models, and with slightly differing mesh resolutions.

The case has been modified to use the final log profile Equilibrium Boundary Conditions, solvers, and organization used in the final cases of the ember project. It still uses the original attempt mesh, which had slightly differing mesh resolutions from Mueller, but now using all the new solvers and methods for a log profile.

No z0 was given for the case in the Mueller paper. So we assumed it to be 0.01, that of grass. All other values should be the same as were described in the paper.

One other major difference is that Mueller appears to have used a more accurate geometry for the tree than the blocks used in this case, converted into some form of mesh resolution block case in the code. Methods for this type of conversion have not yet been developed for this solver, and so a block tree of a single LAD value was used that divided evenly with the mesh, other than that the trunk had to be offset one cell to the side because it didn't divide quite evenly into the mesh. The blocks for the canopy top were estimated in size to get what seemed to make sense for the tree shape.

The script has three choices for Cd*LAD, 0.1, 1, and 5 as was done in the paper. Mueler may have used a slightly different Cd, and therefore a slightly different LAD, but as the solver is configured such that the two vary together, it works to just pick values of Cd*LAD that work to get the desired value of both combined, a choice of 0.5 for Cd was used to make the math a bit easier though technically the choice does not matter so long as the two values together come out to be 0.1, 1, or 5 (or whatever value the user might want to choose). To use the value of choice, just copy and paste the desired setFields file renaming it to setFields before running the setFields command.

The case can also be run to test horizontal homogeneity by running with no porous media present. To do this, just use the same set of commands as are used for a normal run, but skipping the setFields command. This can be helpful for checking the Equilibrium Boundary Condition profiles. Note that this is NOT the method used in any of the standard cases of the ember project, a WindNinja no porous media set of models and solvers was normally used for non-porous simulations, with all commands calling myDragSimpleFoam being done as simpleFoam, and slight changes in the /0/ time directory files, the /constant/momentumTransport file, and other small changes in the /system/ folder. Should come out to be a similar result to do it this way with no calls to simpleFoam though.

Plot scripts have not been developed for the case yet, though they may be developed for it in the future.



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














