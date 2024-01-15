#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"
// #include "createADMFields.H"

//- Make it a ptr with New() method
//- to read parameters from constant/admno1Properties
//- runTimeSelection mechanism -> Selector -> New()
//- ref>>
//       autoPtr<adModel> reaction
//       (
//           CombustionModel<psiReactionThermo>::New(thermo, turbulence())
//       );
adModel reaction(mesh, 1);

// label runMode = 1;
// autoPtr<adModel> reaction(adModel::New(mesh, runMode));
// auto* reaction = new adModel(mesh, runMode);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
