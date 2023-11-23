{
    scalar Sch = 0.1;
    volScalarField D_("diff", turbulence->nu()/Sch);

    fvScalarMatrix Gch4Eqn
    (
        fvm::ddt(Gch4)
      + fvm::div(phi, Gch4)
      - fvm::laplacian(D_, Gch4)
     ==
        fvOptions(Gch4)
    );

    fvOptions.constrain(Gch4Eqn);

    Gch4Eqn.solve();

    fvOptions.correct(Gch4);
}