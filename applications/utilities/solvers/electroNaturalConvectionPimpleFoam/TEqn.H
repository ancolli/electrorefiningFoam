{
    alphat_T = turbulence->nut()/Prt; 
    alphat.correctBoundaryConditions();

    volScalarField alphaEff_T("alphaEff_T", turbulence->nu()/Pr + alphat_T);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff_T, T)
     ==

        (fvc::grad(fi) & fvc::grad(fi))*keff/rho_cp // fuente
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);
    
    Info<< "Min/max T" << ' ' << min(T).value() << ' ' << max(T).value() << endl;
}
