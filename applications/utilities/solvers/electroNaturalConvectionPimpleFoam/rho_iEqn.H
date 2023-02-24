{
    alphat = turbulence->nut()/Sct;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff_c("alphaEff_c", Dmol + alphat);

    fvScalarMatrix rho_iEqn
    (
        fvm::ddt(rho_i)
      + fvm::div(phi, rho_i)
      - fvm::laplacian(alphaEff_c, rho_i)
     ==
        fvOptions(rho_i)
    );

    rho_iEqn.relax();

    fvOptions.constrain(rho_iEqn);

    rho_iEqn.solve();

    fvOptions.correct(rho_i);
   
    Info<< "Min/max C" << ' ' << min(rho_i).value() << ' ' << max(rho_i).value() << endl;
}
