

// TODO: make them in a seperate ADM class
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

wordList namesSoluable = {"Ssu", "Saa", "Sfa", "Sva", "Sbu", "Spro",
                          "Sac", "Sh2", "Sch4", "SIC", "SIN", "SI", "Sco2"};

wordList namesGaseous = {"Gh2", "Gch4", "Gco2"};

wordList namesParticulate = {"Xc", "Xch", "Xpr", "Xli", "Xsu", "Xaa",
                              "Xfa", "Xc4", "Xpro", "Xac", "Xh2", "XI"};

wordList namesElectrolyte = {"Scat", "San", "ShP", "SohN"};

// making ptr list for all the components
// PtrList<volScalarField>& cmptPtrs(28); // 27 + Sco2 
PtrList<volScalarField> YPtrs(16);
label iNames = 0;

//-  Read soluables

Info<< "Reading ADMno1 initial concentrations for soluables" << endl;

forAll(namesSoluable, i)
{
    YPtrs.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                namesSoluable[i], // IOobject::groupName(namesSoluable[i]),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

iNames += namesSoluable.size();

//-  Read gaseuoses

Info<< "Reading ADMno1 initial concentrations for gaseuoses" << endl;

forAll(namesGaseous, i)
{
    YPtrs.set
    (
        i + iNames,
        new volScalarField
        (
            IOobject
            (
                namesGaseous[i],
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

// iNames += sizeof(namesGaseous);

// //-  Read particulates

// Info<< "Reading ADMno1 initial concentrations for particulates" << endl;

// forAll(namesParticulate, i)
// {
//     YPtrs.set
//     (
//         i + iNames,
//         new volScalarField
//         (
//             IOobject
//             (
//                 namesParticulate[i],
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh
//         )
//     );
// }

// iNames += sizeof(namesParticulate);

// //-  Read electrolytes

// Info<< "Reading ADMno1 initial concentrations for electrolytes" << endl;

// forAll(namesElectrolyte, i)
// {
//     YPtrs.set
//     (
//         i + iNames,
//         new volScalarField
//         (
//             IOobject
//             (
//                 namesElectrolyte[i],
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh
//         )
//     );
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


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