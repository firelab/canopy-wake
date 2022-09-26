This will be a drag coefficient solver, written from scratch from simpleFoam, using myBuoyantBoussinesqSimpleFoam and porousSimpleFoam as references.

This is for openfoam 8, an openfoam.org distribution of OpenFOAM. In the end, I actually just copy and pasted simpleFoam and added the sources. So I edited createFields.H and UEqn.H and that's it.

There was a bit of trouble getting it to work for higher values of Cd*LAD. I ended up switching the source term from being specified as explicit to semi-implicit to help.



Okay now the instructions of how to compile it. After compiling and running normal OpenFOAM 8, do the following:

Need to go to this directy and run

    wmake


If need to reclean, run
    
    wclean







Update:

I found a better way to formulate the fvm::SuSp wrapper in the Segersson 2017 paper. The source in U equation is always implicit because there is a u in the eqn, so make it fvm::Sp() every time.
The source in the k equation has one part implicit and one part explicit cause the k is only in one part of the source, Segersson had a trick to do one part implicit ( fvm::Sp() ) and one part explicit (no wrapper or fvm::Su() ) for the k equation. The source in the epsilon equation is all implicit ( fvm::Sp() ) because there is an epsilon in each term.

So apparently the traditional meaning of implicit and explicit is that if the variable appears on both sides of the equation it is implicit. But for OpenFOAM, if the term is always negative, then it should be treated as always implicit with the fvm::Sp() wrapper. If it is always positive then it should be treated as always explicit with the fvm::Su() wrapper (or no wrapper). And if it could be positive or negative, then it should be treated as semi-implicit and use the fvm::SuSp() wrapper (tries to switch back and forth). The reasoning for this has to do with how terms are added or subtracted from the diagonal, making negative terms implicit and positive terms explicit improves the stability of the solver.

Now Segersson changed things up and broke up the k equation into the two separate terms, where one terms is always positive (explicit) and one term is always negative (implicit). But for some odd reason he did not do such a thing for the epsilon equation. And so when I used this version of the solver, it died on me for the simpler test cases, but worked for Segersson's test case. I investigated and found that the epsilon equation should be treated similar to the k equation, breaking up each term and finding one that is explicit and one that is implicit. I found that whether you use the wrapper fvm::Su() for the explicit terms or no, the result is the same.

So I've updated the solver with these methods and found it improves the solution and the stability.











