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

// making ptr list for all the components
// PtrList<volScalarField>& cmptPtrs(28); // 27 + Sco2 
PtrList<volScalarField> compositionPtrs(15);
// soluables
compositionPtrs.set(0, &Ssu);
compositionPtrs.set(1, &Saa);
compositionPtrs.set(2, &Sfa);
compositionPtrs.set(3, &Sva);
compositionPtrs.set(4, &Sbu);
compositionPtrs.set(5, &Spro);
compositionPtrs.set(6, &Sac);
compositionPtrs.set(7, &Sh2);
compositionPtrs.set(8, &Sch4);
compositionPtrs.set(9, &SIC);
compositionPtrs.set(10, &SIN);
compositionPtrs.set(11, &SI);
// gaseous
compositionPtrs.set(12, &Gh2);
compositionPtrs.set(13, &Gch4);
compositionPtrs.set(14, &Gco2);
// TODO: particulate
// TODO: for multiphase as well