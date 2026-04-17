// Explicit diffusion-number estimate used to limit the time step when
// adjustTimeStep is enabled in controlDict.
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    k/(rho*cpEff)
);

tmp<surfaceScalarField> tAlphaF = fvc::interpolate(alpha);

// mesh.deltaCoeffs() gives an inverse distance scale at faces.
tmp<surfaceScalarField> tDiF =
    tAlphaF * runTime.deltaTValue() * sqr(mesh.deltaCoeffs());

const scalar maxDiNow = gMax(tDiF().internalField());

if (adjustTimeStep && maxDiNow > VSMALL)
{
    // Keep the next step below the requested diffusion-number target.
    const scalar dtNew =
        min(maxDeltaT, runTime.deltaTValue()*(maxDiTarget/maxDiNow));

    runTime.setDeltaT(dtNew);
}

Info<< "maxDi = " << maxDiNow
    << "  deltaT = " << runTime.deltaTValue() << nl << endl;
