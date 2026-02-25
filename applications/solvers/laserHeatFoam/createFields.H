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

volScalarField solidificationTime
(
    IOobject
    (
        "solidificationTime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
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

IOdictionary LaserProperties(
    IOobject(
        "LaserProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

// scalar P (readScalar(LaserProperties.lookup("power")));
scalar A (readScalar(LaserProperties.lookup("absorptivity")));
scalar rb(readScalar(LaserProperties.lookup("laserRadius")));
scalar d (readScalar(LaserProperties.lookup("absorptionDepth")));

interpolationTable<scalar> laserPowerSeries
(
    LaserProperties.subDict("timeVsLaserPower")
);

interpolationTable<vector> laserPositionSeries
(
    LaserProperties.subDict("timeVsLaserPosition")
);

IOdictionary thermoDict
(
    IOobject
    (
        "thermophysicalProperties",   // dictionary name
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// 
// heSolidThermo<basicThermo, pureSolid> solidThermo(mesh, "solid");

// rho = solidThermo.rho();
// cp  = solidThermo.Cp();
// k   = solidThermo.kappa();
// scalar L   = solidThermo.L();
// scalar Tm  = solidThermo.Tm();
// volScalarField H = solidThermo.H(T); // sensible + latent

// Material properties
// const scalar rho = 4100.0;     // kg/m³
// const scalar Cp  = 570.0;      // J/kg/K
// const scalar k   = 8.0;        // W/m/K
// const scalar L   = 3.6e5;      // J/kg latent heat
// const scalar Tm  = 1878.0;     // K
// const scalar dT  = 100.0;      // smoothing interval for latent heat

// // Get patch ID by name
// label patchID = mesh.boundaryMesh().findPatchID("topWall");

// // Access the fvPatch
// // const fvPatch& topPatch = mesh.boundary()[patchID];
// const labelList& faceCells = mesh.boundary()[patchID].faceCells();

// // 1️⃣ Get the patch ID
// // label patchID = mesh.boundaryMesh().findPatchID("topWall");

// // 2️⃣ Get the patch field as fvPatchScalarField (non-const)
// fvPatchScalarField& topPatch =
//     refCast<fvPatchScalarField>(T.boundaryField()[patchID]);

// // Radiation cooling
// scalar epsilon = 0.3;
// scalar sigma = 5.67e-8;
// scalar h = 20.0;   // W/m²/K
// scalar Ta = 300.0; // ambient