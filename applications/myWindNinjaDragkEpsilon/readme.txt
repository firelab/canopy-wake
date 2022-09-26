

For this, I took the myWindNinjakEpsilon model from the myWindNinjaLibso, then organized it using my old myDragkEpsilon model, 
adding the extra terms to make it be the same underlying kEpsilon model as the myWindNinjakEpsilon model, but with the drag
terms from the myDragkEpsilon. This is meant to be used alongside the myWindNinjaLibso log profile inlet and wall functions.










Okay now the instructions of how to compile it. After compiling and running normal OpenFOAM 8, do the following:

Need to go to this directory and run

    wmakeLnInclude -u myWindNinjaDragkEpsilon

then run 

    wmake


If need to reclean, run
    
    wclean

then delete the created lnInclude in the /myWindNinjaDragkEpsilon/ folder and repeat the process.


Seems like if you just run wmake in the overall make directory /myWindNinjaLibso/, it does the linker process for you.
Also seems like deleting the lnInclude dir doesn't seem to matter much between wmake and wclean, so technically you 
can probably just do wmake, but if you want to be sure, do the full stuff above with the linking and deleting of the 
linked folder between each major instance of wclean and wmake.













