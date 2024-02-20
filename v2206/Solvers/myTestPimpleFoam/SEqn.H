volScalarField Deff = turbulence->nu();

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
fields.add(Sh2);

tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Sh2)")
    )
);

Info << endl;
Info << "phi: " << phi.dimensions() << endl;
Info << "ddt: " << fvm::ddt(Sh2)->dimensions() << endl;
Info << "div: " << fvm::div(phi, Sh2)->dimensions() << endl;
Info << "mvC: " << mvConvection->fvmDiv(phi, Sh2)->dimensions() << endl;
Info << "lap: " << fvm::laplacian(Deff, Sh2)->dimensions() << endl;
Info << "R(): " << tSu.dimensions() << endl;
Info << "Sh2: " << Sh2.dimensions() << endl;

{
    // fvScalarMatrix Sh2Eqn
    // (
    //     fvm::ddt(rho, Sh2)
    //   + fvm::div(phi, Sh2)
    //   - fvm::laplacian(Deff*rho, Sh2)
    //  ==
    //     fvOptions(rho, Sh2)
    //   + tSu
    // );

    fvScalarMatrix Sh2Eqn
    (
        fvm::ddt(Sh2)
      + fvm::div(phi, Sh2)
      - fvm::laplacian(Deff, Sh2)
     ==
        fvOptions(Sh2)
      + tSu
    );

    Sh2Eqn.relax();
    fvOptions.constrain(Sh2Eqn);
    Sh2Eqn.solve(mesh.solver("Sh2"));
    fvOptions.correct(Sh2);
}
