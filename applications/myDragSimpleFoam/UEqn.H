    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevSigma(U)
     ==
	  - fvm::Sp(   Cd*LAD*mag(U) , U   )		// line added. Used to have -Cd*LAD*mag(U)*U, but the solver blew up at higher Cd*LAD vals. Makes sense to make the term always implicit since the term is always negative (assuming U doesn't ever go negative from numerical error).
      + fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }


