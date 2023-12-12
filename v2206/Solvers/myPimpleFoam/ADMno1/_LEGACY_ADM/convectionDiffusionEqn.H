{

    // TODO: make it an input from ./constant/ADMpara file
    scalar SC = 0.1; // Schmidt Number
    volScalarField Deff = turbulence->nu() / SC;

    // define convection-diffusion equations for soluables

    fvScalarMatrix SsuEqn
    (
        fvm::ddt(Ssu)
      + fvm::div(phi, Ssu)
      - fvm::laplacian(Deff, Ssu)
     ==
        fvOptions(Ssu)
    );

    fvScalarMatrix SaaEqn
    (
        fvm::ddt(Saa)
      + fvm::div(phi, Saa)
      - fvm::laplacian(Deff, Saa)
     ==
        fvOptions(Saa)
    );

    fvScalarMatrix SfaEqn
    (
        fvm::ddt(Sfa)
      + fvm::div(phi, Sfa)
      - fvm::laplacian(Deff, Sfa)
     ==
        fvOptions(Sfa)
    );

    fvScalarMatrix SvaEqn
    (
        fvm::ddt(Sva)
      + fvm::div(phi, Sva)
      - fvm::laplacian(Deff, Sva)
     ==
        fvOptions(Sva)
    );

    fvScalarMatrix SbuEqn
    (
        fvm::ddt(Sbu)
      + fvm::div(phi, Sbu)
      - fvm::laplacian(Deff, Sbu)
     ==
        fvOptions(Sbu)
    );

    fvScalarMatrix SproEqn
    (
        fvm::ddt(Spro)
      + fvm::div(phi, Spro)
      - fvm::laplacian(Deff, Spro)
     ==
        fvOptions(Spro)
    );

    fvScalarMatrix SacEqn
    (
        fvm::ddt(Sac)
      + fvm::div(phi, Sac)
      - fvm::laplacian(Deff, Sac)
     ==
        fvOptions(Sac)
    );

    fvScalarMatrix Sh2Eqn
    (
        fvm::ddt(Sh2)
      + fvm::div(phi, Sh2)
      - fvm::laplacian(Deff, Sh2)
     ==
        fvOptions(Sh2)
    );

    fvScalarMatrix Sch4Eqn
    (
        fvm::ddt(Sch4)
      + fvm::div(phi, Sch4)
      - fvm::laplacian(Deff, Sch4)
     ==
        fvOptions(Sch4)
    );

    fvScalarMatrix SICEqn
    (
        fvm::ddt(SIC)
      + fvm::div(phi, SIC)
      - fvm::laplacian(Deff, SIC)
     ==
        fvOptions(SIC)
    );

    fvScalarMatrix SINEqn
    (
        fvm::ddt(SIN)
      + fvm::div(phi, SIN)
      - fvm::laplacian(Deff, SIN)
     ==
        fvOptions(SIN)
    );

    fvScalarMatrix SIEqn
    (
        fvm::ddt(SI)
      + fvm::div(phi, SI)
      - fvm::laplacian(Deff, SI)
     ==
        fvOptions(SI)
    );

    // define convection-diffusion equations for gaseous

    fvScalarMatrix Gh2Eqn
    (
        fvm::ddt(Gh2)
      + fvm::div(phi, Gh2)
      - fvm::laplacian(Deff, Gh2)
     ==
        fvOptions(Gh2)
    );

    fvScalarMatrix Gch4Eqn
    (
        fvm::ddt(Gch4)
      + fvm::div(phi, Gch4)
      - fvm::laplacian(Deff, Gch4)
     ==
        fvOptions(Gch4)
    );

    fvScalarMatrix Gco2Eqn
    (
        fvm::ddt(Gco2)
      + fvm::div(phi, Gco2)
      - fvm::laplacian(Deff, Gco2)
     ==
        fvOptions(Gco2)
    );

    // solve equations for soluables

    fvOptions.constrain(SsuEqn);
    fvOptions.constrain(SaaEqn);
    fvOptions.constrain(SfaEqn);
    fvOptions.constrain(SvaEqn);
    fvOptions.constrain(SbuEqn);
    fvOptions.constrain(SproEqn);
    fvOptions.constrain(SacEqn);
    fvOptions.constrain(Sh2Eqn);
    fvOptions.constrain(Sch4Eqn);
    fvOptions.constrain(SICEqn);
    fvOptions.constrain(SINEqn);
    fvOptions.constrain(SIEqn);
    
    SsuEqn.solve();
    SaaEqn.solve();
    SfaEqn.solve();
    SvaEqn.solve();
    SbuEqn.solve();
    SproEqn.solve();
    SacEqn.solve();
    Sh2Eqn.solve();
    Sch4Eqn.solve();
    SICEqn.solve();
    SINEqn.solve();
    SIEqn.solve();

    fvOptions.correct(Ssu);
    fvOptions.correct(Saa);
    fvOptions.correct(Sfa);
    fvOptions.correct(Sva);
    fvOptions.correct(Sbu);
    fvOptions.correct(Spro);
    fvOptions.correct(Sac);
    fvOptions.correct(Sh2);
    fvOptions.correct(Sch4);
    fvOptions.correct(SIC);
    fvOptions.correct(SIN);
    fvOptions.correct(SI);

    // solve equations for gaseous

    fvOptions.constrain(Gh2Eqn);
    fvOptions.constrain(Gch4Eqn);
    fvOptions.constrain(Gco2Eqn);

    Gh2Eqn.solve();
    Gch4Eqn.solve();
    Gco2Eqn.solve();

    fvOptions.correct(Gh2);
    fvOptions.correct(Gch4);
    fvOptions.correct(Gco2);

    // define convection-diffusion equations for sluges

}