dimensionedScalar TopDefault("TopDefault", dimensionSet(0,0,0,1,0,0,0), 308.15); 
Info<< "Reading field Operating Temperature" << endl;
volScalarField Top
(
    IOobject
    (
        "Top",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    TopDefault
);

Info<< "Reading field Soluables su" << endl;
volScalarField Ssu
(
    IOobject
    (
        "Ssu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables aa" << endl;
volScalarField Saa
(
    IOobject
    (
        "Saa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables fa" << endl;
volScalarField Sfa
(
    IOobject
    (
        "Sfa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables va" << endl;
volScalarField Sva
(
    IOobject
    (
        "Sva",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables bu" << endl;
volScalarField Sbu
(
    IOobject
    (
        "Sbu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables pro" << endl;
volScalarField Spro
(
    IOobject
    (
        "Spro",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables ac" << endl;
volScalarField Sac
(
    IOobject
    (
        "Sac",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables h2" << endl;
volScalarField Sh2
(
    IOobject
    (
        "Sh2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables ch4" << endl;
volScalarField Sch4
(
    IOobject
    (
        "Sch4",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables Inorganic Carbon" << endl;
volScalarField SIC
(
    IOobject
    (
        "SIC",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables Inorganic Nitrogen" << endl;
volScalarField SIN
(
    IOobject
    (
        "SIN",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Soluables Inerts" << endl;
volScalarField SI
(
    IOobject
    (
        "SI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Gaseous h2" << endl;
volScalarField Gh2
(
    IOobject
    (
        "Gh2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Gaseous ch4" << endl;
volScalarField Gch4
(
    IOobject
    (
        "Gch4",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Gaseous co2" << endl;
volScalarField Gco2
(
    IOobject
    (
        "Gco2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);