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


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// >>> TEST ------------------------------------------------------ //

volScalarField rho
(
    IOobject
    (
        "rho",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        dimDensity, 
        One
    )
);

volScalarField Sh2
(
    IOobject
    (
        "Sh2",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField dSh2
(
    IOobject
    (
        "dSh2",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        dimMass/dimVolume/dimTime, 
        0.01
    )
);

// tmp<fvScalarMatrix> tSu
// (
//     new fvScalarMatrix
//     (
//         Sh2,
//         dimMass/dimTime
//     )
// );
// fvScalarMatrix& Su = tSu.ref();
// Su += dSh2; // <<< causing dimension issue!!!!

// tmp<fvScalarMatrix> tSu
// (
//     new fvScalarMatrix
//     (
//         Sh2,
//         dimVolume/dimTime
//     )
// );
// fvScalarMatrix& Su = tSu.ref();
// Su += (dSh2/rho); // <<< causing dimension issue!!!!

fvScalarMatrix tSu
(
    Sh2,
    dimVolume/dimTime
    
);
// scalarField s(tSu.source().size(), 1e-2);
// tSu.source() += s;
tSu += (dSh2/rho); // <<< causing dimension issue!!!!

// --------------------------------------------------------------- //

#include "createMRF.H"
#include "createFvOptions.H"
