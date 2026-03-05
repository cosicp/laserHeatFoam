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

volScalarField meltHistory
(
    IOobject
    (
        "meltHistory",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

#include "createFvOptions.H"

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

volScalarField cpEff
(
    IOobject("cpEff", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    cp
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
[[maybe_unused]] scalar A (readScalar(LaserProperties.lookup("absorptivity")));
[[maybe_unused]] scalar rb(readScalar(LaserProperties.lookup("laserRadius")));
[[maybe_unused]] scalar d (readScalar(LaserProperties.lookup("absorptionDepth")));

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

scalar Ts (readScalar(thermoDict.lookup("Ts")));
scalar Tl (readScalar(thermoDict.lookup("Tl")));
dimensionedScalar Tinf("Tinf", dimTemperature, 300);   // K
dimensionedScalar L("L", dimEnergy/dimMass, readScalar(thermoDict.lookup("L")));   // J/kg

const dictionary& controlDict = runTime.controlDict();

Switch adjustTimeStep
(
    controlDict.lookupOrDefault<Switch>("adjustTimeStep", false)
);

scalar maxDeltaT
(
    controlDict.lookupOrDefault<scalar>("maxDeltaT", GREAT)
);

scalar maxDiTarget
(
    controlDict.lookupOrDefault<scalar>("maxDi", 0.2)
);

volScalarField fL
(
    IOobject("fL", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField fLold
(
    IOobject("fLold", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    fL
);

volScalarField betaLoss
(
    IOobject("betaLoss", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimPower/dimVolume/dimTemperature, 0.0) // W/m3/K
);
