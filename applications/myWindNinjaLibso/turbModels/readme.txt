

For this, I took my old myDragkEpsilon case, and used that organization, but copying and pasting the OpenFOAM 8 kEpsilon model.
I then commented out stuff, deleted stuff, and added in stuff from the OpenFOAM 2.2.0 incompressible version, which is what is used
by WindNinja. The goal was to try to replicate what is done by WindNinja more correctly, since it seemed there are extra 
terms and constants and stuff in the OpenFOAM 8 version of the kEpsilon model that weren't in the OpenFOAM 2.2.0 incompressible
version used by WindNinja. There were also some differences in the different formulations of the variables given that go into the 
equation terms (G, R(), devReff() style functions that are now sigma() and devTau() style functions).















Okay now the instructions of how to compile it. After compiling and running normal OpenFOAM 8, do the following:

Need to go to this /myWindNinjaLibso/turbModels/ directory and run

    wmakeLnInclude -u myWindNinjakEpsilon

then in the overall make folder /myWindNinjaLibso/ run 

    wmake


If need to reclean, in the overall make folder /myWindNinjaLibso/, run
    
    wclean

then delete the created lnInclude in the /myWindNinjakEpsilon/ folder and repeat the process.


Seems like if you just run wmake in the overall make directory /myWindNinjaLibso/, it does the linker process for you.
Also seems like deleting the lnInclude dir doesn't seem to matter much between wmake and wclean, so technically you 
can probably just do wmake, but if you want to be sure, do the full stuff above with the linking and deleting of the 
linked folder between each major instance of wclean and wmake.













