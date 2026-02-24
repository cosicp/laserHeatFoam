Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading diffusivity k\n" << endl;

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity, Zero)
);

// if (!DT.headerOk())
// {
//     IOdictionary transportProperties
//     (
//         IOobject
//         (
//             "transportProperties",
//             runTime.constant(),
//             mesh,
//             IOobject::MUST_READ_IF_MODIFIED,
//             IOobject::NO_WRITE
//         )
//     );
//     DT = dimensionedScalar("DT", dimViscosity, transportProperties);
// }

#include "createFvOptions.H"

// volScalarField rho(
//     IOobject(
//         "rho",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::NO_WRITE),
//     alphaM *((gT * rhoL) + ((scalar(1) - gT) * rhoS)) + (scalar(1) - alphaM) * rhoG,
//     zeroGradientFvPatchScalarField::typeName);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);