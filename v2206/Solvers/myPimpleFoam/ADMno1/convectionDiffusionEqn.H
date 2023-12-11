{
    // TODO: make it an input from ./constant/ADMpara file
    scalar SC = 1e-2; // Schmidt Number
    volScalarField Deff = turbulence->nu() / SC;
    fvScalarMatrix Sch4Eqn
    (
        fvm::ddt(Sch4)
      + fvm::div(phi, Sch4)
      - fvm::laplacian(Deff, Sch4)
     ==
        fvOptions(Sch4)
    );

    fvOptions.constrain(Sch4Eqn);

    Sch4Eqn.solve();

    fvOptions.correct(Sch4);
}