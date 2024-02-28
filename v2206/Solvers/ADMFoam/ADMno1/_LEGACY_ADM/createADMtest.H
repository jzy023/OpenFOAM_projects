dimensionedScalar TopDefault("TopDefault", dimensionSet(0,0,0,1,0,0,0), 308.15); 
Info<< "Reading field Operating Temperature" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    TopDefault
);



Info<< "Reading field Soluables" << endl;
// GeometricField<admSoluables<scalar>, fvPatchField, volMesh> Sol
volADMSField Sol
(
    IOobject
    (
        "Sol",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

