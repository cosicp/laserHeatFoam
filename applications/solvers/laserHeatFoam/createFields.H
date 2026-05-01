Info<< "Reading field T\n" << endl;

// Primary solution field.
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

// Effective thermal conductivity, rebuilt from phase fractions each time step.
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

// Cell marker for material that has reached the liquidus at least once.
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

// Persistent per-cell layer index on the base mesh. Filled by
// multiLayer.H from the BFS-derived per-layer cell sets and then mapped
// automatically by dynamicRefineFvMesh when the base mesh is refined
// (children inherit their parent's layerID, which is exactly what we
// want — a child of a layer-3 cell is still in layer 3).
//
// Cells with layerID < 0 are never activated (e.g. cells beyond
// nLayers*cellsPerLayer in the thickness direction).
volScalarField layerID
(
    IOobject
    (
        "layerID",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("layerID", dimless, -1.0)
);

// Refinement indicator: an interface-peaked function of a smoothed
// meltHistory. Computed each step by updateRefineIndicator.H and read
// by dynamicRefineFvMesh to decide where to refine/unrefine.
volScalarField refineIndicator
(
    IOobject
    (
        "refineIndicator",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

#include "createFvOptions.H"

// Thermophysical fields read from the case and updated during the run.
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
    IOobject
    (
        "cpEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    cp
);

// Volumetric laser source [W/m3].
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

// Laser model selection, scan path, and power history.
IOdictionary LaserProperties
(
    IOobject
    (
        "LaserProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word laserModel
(
    LaserProperties.lookupOrDefault<word>("laserModel", "cylindrical")
);

[[maybe_unused]] const scalar substrateProbeDepth
(
    LaserProperties.lookupOrDefault<scalar>
    (
        "substrateProbeDepth",
        LaserProperties.lookupOrDefault<scalar>("absorptionDepth", 0.0)
    )
);

interpolationTable<scalar> laserPowerSeries
(
    LaserProperties.subDict("timeVsLaserPower")
);

interpolationTable<vector> laserPositionSeries
(
    LaserProperties.subDict("timeVsLaserPosition")
);

// Material properties and coefficients.
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
dimensionedScalar L("L", dimEnergy/dimMass, readScalar(thermoDict.lookup("L")));

// Temperature-dependent property model coefficients.
[[maybe_unused]] const scalar T1val =
    thermoDict.lookupOrDefault<scalar>("T1", 1268.0);
[[maybe_unused]] const scalar Tlval = Tl.value();

// Powder porosity (Mohammadkamal eq. 1).
[[maybe_unused]] const scalar phiVal =
    thermoDict.lookupOrDefault<scalar>("phi", 0.35);

// Density coefficients.
[[maybe_unused]] const scalar rho0s =
    thermoDict.lookupOrDefault<scalar>("rho0s", 4557.892);
[[maybe_unused]] const scalar rho1s =
    thermoDict.lookupOrDefault<scalar>("rho1s", -0.154);
[[maybe_unused]] const scalar rho0m =
    thermoDict.lookupOrDefault<scalar>("rho0m", 5227.640);
[[maybe_unused]] const scalar rho1m =
    thermoDict.lookupOrDefault<scalar>("rho1m", -0.680);

// Specific heat coefficients.
[[maybe_unused]] const scalar cp0a =
    thermoDict.lookupOrDefault<scalar>("cp0a", 483.04);
[[maybe_unused]] const scalar cp1a =
    thermoDict.lookupOrDefault<scalar>("cp1a", 0.215);
[[maybe_unused]] const scalar cp0b =
    thermoDict.lookupOrDefault<scalar>("cp0b", 412.70);
[[maybe_unused]] const scalar cp1b =
    thermoDict.lookupOrDefault<scalar>("cp1b", 0.18);
[[maybe_unused]] const scalar cpL =
    thermoDict.lookupOrDefault<scalar>("cpL", 831.0);

// Solid/melt conductivity coefficients.
[[maybe_unused]] const scalar k0a =
    thermoDict.lookupOrDefault<scalar>("k0a", 1.250);
[[maybe_unused]] const scalar k1a =
    thermoDict.lookupOrDefault<scalar>("k1a", 0.015);
[[maybe_unused]] const scalar k0b =
    thermoDict.lookupOrDefault<scalar>("k0b", 3.15);
[[maybe_unused]] const scalar k1b =
    thermoDict.lookupOrDefault<scalar>("k1b", 0.012);
[[maybe_unused]] const scalar k0m =
    thermoDict.lookupOrDefault<scalar>("k0m", -12.752);
[[maybe_unused]] const scalar k1m =
    thermoDict.lookupOrDefault<scalar>("k1m", 0.024);

// Phase-dependent conductivity model:
//   k = rp*kp + rm*km + rs*ks
// where kp is a constant powder conductivity and ks/km are updated with T.
const dimensionedScalar kp
(
    "kp",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    readScalar(thermoDict.lookup("kp"))
);

volScalarField ks
(
    IOobject
    (
        "ks",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "ks",
        dimensionSet(1, 1, -3, -1, 0, 0, 0),
        readScalar(thermoDict.lookup("ks"))
    )
);

volScalarField km
(
    IOobject
    (
        "km",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "km",
        dimensionSet(1, 1, -3, -1, 0, 0, 0),
        readScalar(thermoDict.lookup("km"))
    )
);

// Consolidated fraction rc:
//   rc = 1  for consolidated solid/substrate
//   rc = 0  for unconsolidated powder
// rc only increases, giving an irreversible powder-to-melt transition.
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

// Explicit diffusion-number control. These are marked maybe_unused because
// postProcess.H also includes createFields.H without entering the time loop.
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

// Phase fraction fields.
volScalarField fL
(
    IOobject
    (
        "fL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField fLold
(
    IOobject
    (
        "fLold",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fL
);

volScalarField rp
(
    IOobject
    (
        "rp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField rm
(
    IOobject
    (
        "rm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField rs
(
    IOobject
    (
        "rs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
