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

// k is computed from phase fractions every time step.
// Dimension: thermal conductivity [W/(m·K)] = [kg m s^-3 K^-1]
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("k", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0)
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

word laserModel(LaserProperties.lookupOrDefault<word>("laserModel", "cylindrical"));

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

dimensionedScalar Ts("Ts", dimTemperature, readScalar(thermoDict.lookup("Ts")));
dimensionedScalar Tl("Tl", dimTemperature, readScalar(thermoDict.lookup("Tl")));
dimensionedScalar Tinf("Tinf", dimTemperature, 300);   // K
dimensionedScalar L("L", dimEnergy/dimMass, readScalar(thermoDict.lookup("L")));   // J/kg

// Phase-dependent conductivities (Proell 2020, eq. 5 / Proell 2023, eq. 5)
// k = rp*kp + rm*km + rs*ks
// kp: powder (constant - highly porous, low conductivity)
// ks, km: temperature-dependent solid/melt conductivities (updated each timestep)
const dimensionedScalar kp("kp", dimensionSet(1, 1, -3, -1, 0, 0, 0),
    readScalar(thermoDict.lookup("kp")));

volScalarField ks
(
    IOobject("ks", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("ks", dimensionSet(1, 1, -3, -1, 0, 0, 0),
        readScalar(thermoDict.lookup("ks")))
);

volScalarField km
(
    IOobject("km", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("km", dimensionSet(1, 1, -3, -1, 0, 0, 0),
        readScalar(thermoDict.lookup("km")))
);

// Consolidated fraction rc (Proell 2020, eq. 23):
//   rc(0) = 1  for initially solid/substrate cells
//   rc(0) = 0  for initially powder cells
// rc can only increase (irreversible powder->melt transition).
// Initialize from file if present; default = 0 (all powder).
volScalarField rc
(
    IOobject
    (
        "rc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

const dictionary& controlDict = runTime.controlDict();

[[maybe_unused]] Switch adjustTimeStep
(
    controlDict.lookupOrDefault<Switch>("adjustTimeStep", false)
);

[[maybe_unused]] scalar maxDeltaT
(
    controlDict.lookupOrDefault<scalar>("maxDeltaT", GREAT)
);

[[maybe_unused]] scalar maxDiTarget
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

volScalarField rp
(
    IOobject("rp", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField rm
(
    IOobject("rm", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField rs
(
    IOobject("rs", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
