

This case is essentially an attempt to replicate the ember transport model simulations described in the 2011 Thesis work of Hillary Harris 2011 "Analysis and Parameterization of the Flight of Ember Generation Experiments", which were an attempt to replicate ember generation and transport experiments performed at the Institute for Business Home and Safety (IBHS) facility humongous more than garage sized wind tunnel. Some notable differences in this replication is that the Harris 2011 methods appear to have just been 2D Lagrangian Particle Transport equations, this replication is in 3D. Also Harris stated that the measured 10 minute average turbulence intensity (I_x = sigma_ux/umean, I_z = sigma_uz/umean) was a really small value of 0.014, estimating the velocity fluctuations to be ufluct = umean*sqrt( 1 + I_u^2 ) and setting I_u to be I_x this came out to be uFluct ~= 1.0001*umean and so Harris 2011 ignored the velocity fluctuations, this replication does not ignore the velocity fluctuations but calculates them using the Lagrangian particle turbulence model described in equation 5 solved as equation 11 from Bailey 2017 "Numerical Considerations for Lagrangian Stochastic Dispersion Models: Eliminating Rogue Trajectories, and the Importance of Numerical Accuracy", even though the velocity fluctuations ended up coming out to be quite small.

In the simulations described in Harris 2011, embers of a lognormal distribution of mass with mean µ of -3.2811 and variance σ of 0.94289 grams were injected into a 20 mph wind at 6 ft, 11 ft, 21 ft heights with mean launch velocities of 6 m/s, 7.5 m/s, and 9 m/s respectively and mean vertical angles of 45 degrees. The launch velocities had a standard normal distribution between (+) or (-) 2 m/s from the mean launch velocity, the vertical launch angles were normally distributed with a variance of 10 degrees. A total of 190,000 particles were launched with 80,000 launched from the 6 ft height, with 60,000 launched from the 11 ft height, and with 50,000 launched from the 21 ft height. The resulting ember depositions were then compared against the experimental data, where the embers were launched from large firebrand generators, with embers caught in 18 water filled pans of size 25.5" x 17.5" and 1" placed from 12.2 m to 22.56 m straight downwind of the ember generators, pans placed such that their longest 25.5" dimension was normal to the wind direction and with 2 ft (~0.61 m) between each pan center.

Replicating the simulations described by Harris 2011 required additional information. Harris 2011 described the plumes coming from each ember generator to go straight downwind, extending out to the sides approximately 10 degrees to either side from the centerline. Using the pictures provided by Harris 2011 of the initial cones of firebrands produced by the firebrand generators, an initial shape for a cone injection was estimated. Harris 2011 described the pipe exits of the firebrand generators to be 6 inches in diameter, but the pictures showed embers to be exiting in a region that appeared to be only the top 1/3rd of the pipe exits. So, a cone injection diameter of 1/6 ft was used, with pipe exits assumed to be placed with the bottom of the pipe exits at the injection heights used by Harris, (injection goes from 6 2/6 ft to 6 3/6 ft, 11 2/6 ft to 11 3/6 ft, 21 2/6 ft to 12 3/6 ft), the cone injection extending into the domain 0.5 ft to allow enough space for the injections. In addition, no density was given in the paper, so an ember density of 300 kg/m3 was assumed (Finney 2004 "FARSITE: Fire Area Simulator--Model Development and Evaluation" gives that value for charred wood cylinders). Using this density, and assuming spherical particles, and taking the min and maximum weighted embers of 0.01 grams to 0.19 grams, this gave an approximate range of 1 mm to 30 mm for the lognormal distribution of ember sizes. This range of sizes also fits well as a good rough approximation with the sizes of embers shown in photographs by Harris 2011 of a sample of charred embers.

Some additional difficulties were experienced when trying to replicate the results. Initially, it was thought that the mean and variance of the lognormal distribution was in terms of mm, not grams. This caused the lognormal mean and variance given by Harris to cut off ember sizes that were below about 2 mm, disallowing use of the full desired distribution of sizes. To get around this, the lognormal distribution was estimated using a normal distribution with mean 0.01 and variance 0.015, with limits of 1 mm to 30 mm, which gave a distribution similar to the initial lognormal distribution seen when using the mean and variance as mm instead of grams, other than that the distribution was a little less steep to the left of the mean and steeper to the right of the mean. It IS possible that directly converting the mean and variance provided by Harris from grams to mm could allow these mean and variance values to now work correctly to actually include the 1 mm sized particles when doing a lognormal distribution, but we ended up just sticking to our original figured out normal distribution. The simulation was also performed without limits to the normal distributions of launch speed and horizontal launch angle, the horizontal launch angle varying the same as the vertical launch angle which appeared to not be true when analyzing the results.

The results still don't fully match the IBHS data or the simulation data described in Harris 2011, but the results are quite close. In trying to do the final replication methods for the case, simulations using all the same values but choosing between using no launch velocity variance and no launch angle variance, launch velocity variance only, angle velocity variance only, and both launch velocity and angle variances, each with single time and 30 second time releases, it was found that using both launch velocity and angle variances was the best, but only if the variances were kept as large as described in Harris 2011, that the variance in launch velocity was required for the larger particles to even reach far enough downwind, and that the launch variance controlled how much of a spread in a given direction the particles experienced. It was surmised that the departure from the best case results from the results described in Harris 2011 was because the launch velocity variance was not limited to (+) or (-) 2 m/s but was allowed to go way larger and smaller than these limits, that it caused a large portion of embers to land just barely on the front edge of the pans that were normally distributed more upwind and downwind. It also appears to be that the launch angle variance needs to be decoupled, that the horizontal and vertical launch angle variances needed to be treated separately rather than assuming the horizontal launch angle variance matched the vertical launch angle variance. This would mean that the vertical launch angle variance is correct, that it is a normal distribution with no limits using the same velocity variance value, but the horizontal launch angle variance should be limited to 10 deg on each side of the plume center, still using the similar value for the angle variance. It was surmised that to further improve the simulation results to not just match the simulation results described in Harris 2011, but to more closely match the ember generation experimental data from the IBHS experiment described in Harris, that the ember shapes need to be more cylindrical. There may be some small improvements that could still be made for the particle turbulence model as well, but as the turbulence intensity of the IBHS facility was so small, improving the particle properties would likely be more important. Whatever is missing, it seems related to giving the embers just a little more lift than is currently done, a greater portion of the larger particles seems to have been able to get further downwind in the experimental data than is done with the current modeling methods.



An additional folder with a matlab script used to figure out the lognormal vs normal size distributions was included, called /zz_lognormalDistributionScriptStuff/, which includes various pictures of the lognormal vs normal distributions. Keep in mind that these pictures were made before it was realized that the mean and variance values described for the lognormal distribution in Harris 2011 were actually in terms of grams and not mm. Pictures were named with line numbers in an attempt to help track down which code lines of the matlab script were used to generate the given picture and its corresponding size distributions. In practice, you can set the corresponding min and max limit values, and the corresponding mean and variance values, for the distributions in the kinematicCloudProperties, then run the particle simulations, then copy and paste a given diameters file from the particle simulations to excel to sort and plot the resulting diameters to see how they match the distributions described in the matlab script, the resulting distributions come out slightly different from the matlab script result, in particular you can see if the distribution somehow cuts off or truncates off certain sizes from the distribution that might be expected when doing calculations to select the desired values.





This case was saved such that to run the particle simulations, you need to first run the wind simulation, then copy and paste the result files and mesh from the wind simulation to the particle simulation. Normally, this would mean copying and pasting the mesh from the latest time directory that has mesh files from the wind sim to the /0/ time directory of the particle sim, the mesh from /constant/ of the wind sim to the /constant/ of the particle sim, the /system/ mesh generation files from the wind sim to the /system/ folder of the particle sim, all the result files from the /latestTime/ time directory of the wind sim to the /0/ directory of the particle sim (make sure to not copy over the /uniform/ folder in this copy), the momentumTransport and transportProperties files from the /constant/ directory of the wind sim to the /constant/ directory of the particle sim copying and pasting the full "rho" entry from the transportProperties file to a new line of the transportProperties file renaming "rho" to be "rhoInf", copying the fvSchemes, fvSolution, and decomposeParDict files from the /system/ directory of the wind sim to the /system/ folder of the particle sim then modifying the fvSchemes file to be a ddtScheme of Euler instead of SteadyState, copying and pasting and then doing all necessary adjustments to change to the new simulation type and simulation times of the controlDict file from the /system/ folder of the wind sim to the /system/ folder of the particle sim, and adding the "g" and "kinematicCloudProperties" to the /constant/ folder of the particle simulation. But many of these files are already copied over, added, and modified to be what is required to run the particle simulation so that all that is left to do is to copy and paste the mesh files from the /constant/ directory of the wind sim to the /constant directory of the particle sim (in this case, this is the /polyMesh/ folder only, other cases have even more mesh files) and to copy and paste the final result /latestTime/ time directory files from the wind sim to the /0/ time directory of the particle simulation.

For particle simulations, the general practice is to run them once for the very first timestep, the end time and write timestep set to the simulation timestep, the log file of this one timestep simulation saved separately from the log file of the rest of the particle simulation, so that the initial injection of particles can be inspected. Then the write timestep is set to something reasonable that is frequently enough to see what is going on with the simulation, but not so small that you get an insane amount of data (which is very easy to do), the end time set to some large time, the simulation cut before the final end time at the write timestep just after the number of particles still in the domain matches the number of particles that have stuck, which is when the last particle has either left the domain or gone inactive. Unfortunately, the Lagrangian particle code has not been set to carry over the active and sticking information and so while you can potentially split the simulation into write timesteps to try to see what is going on more particularly, if you ever let the simulation stop and start again then you have to manually keep track of the number of stuck particles to know when the simulation is done, hence why separate log files for each particle simulation run is very helpful.

The current configuration of the kinematicCloudProperties file in the particle simulation /constant/ directory is set to the best configuration found so far for replicating the Harris 2011 IBHS data. But there were other attempts made before reaching that configuration, the kinematicCloudProperties file has options to switch some of the values with commented values. In particular, a standard cone injection with just standard deviation in initial particle velocity can be used by changing each injection type from myIBHSConeInjection to myConeInjection, the standard distributions of initial velocity and angles can also be changed. There is also an option to switch from the normal distribution attempting to mimic the lognormal distribution given by Harris 2011, to the lognormal distribution that was used before discovering that the set lognormal mean and variance values didn't allow quite small enough particles as desired for the input particle size distribution.

There is also an option to re-enable the myParCalcs helper function, but it was found to be defective when running the particle simulations in parallel, there seems to still be something off with the algorythm used in it to detect which particles are which and when a new loop has started.

The current kinematicCloudProperties file is set to be a 30 second release. Single time releases end up having almost the exact same, if not the exact same, result as the same simulation released over time, and are less expensive to simulate. The 30 second releases end up being samples of the full injection done in a single time release of the same configuration, just splitting up the single time release to be a release over time. The two different release types are left in the /constant/ directory of the particle simulation case, where you can just switch between the two types of releases by  copying and pasting the desired kinematicCloudProperties file renaming it to kinematicCloudProperties before running the particle simulation.




The case COULD have been more particular in making sure the wind simulation BCs matched a more traditional wind tunnel simulation, these wind simulations were done using domain average inlet winds, but mixed with ideas for doing external ABL flows WITHOUT being log profile winds, so inletOutlet BCs rather than walled BCs. But as the resulting wind field was a single uniform wind of the same value as the inlet wind, it was assumed that this configuration of BCs for the wind sim was good enough. Future work could be to try out more correct wind tunnel style BCs, especially to see if this alters the turbulence fields.





To run the case, use the following guide.


standard set of commands for wind simulations in parallel:

  blockMesh 2>&1 | tee blockMesh.log
  checkMesh 2>&1 | tee checkMesh.log
  decomposePar 2>&1 | tee decomposePar.log
  mpirun -np 3 simpleFoam -parallel 2>&1 | tee simpleFoam.log
  reconstructPar 2>&1 | tee reconstructPar.log

At the end after simulations, need to get the residuals:
  foamMonitor -l postProcessing/residuals/0/residuals.dat




For particle simulations, originally they were done in parallel using decomposePar, but later it was discovered that the particle part of the simulations was still in series, the decomposePar was just for how the particle solver dealt with the Eulerian wind and turbulence fields. And it turned out to be faster to just run the particle simulations in series on a single processor. The parallel commands are given here just in case:
  
  decomposePar 2>&1 | tee decomposePar.log
  
For all cases with the above commands already done, the following commands are used to run the particle simulations:
  
  mpirun -np 3 particleFoam -parallel 2>&1 | tee particleFoam.log
  reconstructPar 2>&1 | tee reconstructPar.log
  
To run particle simulations the correct way, NOT in parallel, do the following:
  
  particleFoam 2>&1 | tee particleFoam.log



After particle simulations are done, the mass fractions of the embers that have landed in the pan regions can be generated by running the calcIBHStestCaseMassFractions.py function. From the command line within the particle simulation directory run the following to generate the particle position information in format useable by the script:
  
  calcParPositionsFromCoords -latestTime
  
Then go to the directory with the script and run the python script:
  
  python calcIBHStestCaseMassFractions.py
  









