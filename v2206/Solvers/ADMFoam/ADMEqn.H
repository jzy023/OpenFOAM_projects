Info << "ADM convection-diffusion calculation" << endl;

PtrList<volScalarField>& YPtrs = reaction->Y();
PtrList<volScalarField>& GPtrs = reaction->G();

volScalarField Deff = turbulence->nu() / reaction->Sc();

{
    reaction->clear();
    reaction->correct(phi, T);

    // Liquid Phase ============================================================ //
    forAll(YPtrs, i)
    {
        if(i != 7)
        {
            volScalarField& Yi = YPtrs[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(Yi) == reaction->R(i)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);
        }
    }

    // Pseudo Gas Phase ======================================================== //
    forAll(GPtrs, i)
    {
        volScalarField& Gi = GPtrs[i];

        fvScalarMatrix GiEqn
        (
            fvm::ddt(Gi) == reaction->RG(i)
        );

        GiEqn.relax();

        fvOptions.constrain(GiEqn);

        GiEqn.solve(mesh.solver("Yi"));

        fvOptions.correct(Gi);
    }


}