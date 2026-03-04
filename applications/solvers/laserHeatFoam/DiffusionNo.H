// Di ~ alpha*dt/dx^2, alpha = k/(rho*cpEff)
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

tmp<surfaceScalarField> tDiF =
    tAlphaF * runTime.deltaTValue() * sqr(mesh.deltaCoeffs());

scalar maxDiNow = gMax(tDiF().internalField());

if (adjustTimeStep && maxDiNow > VSMALL)
{
    scalar dtNew = min(maxDeltaT, runTime.deltaTValue()*(maxDiTarget/maxDiNow));
    runTime.setDeltaT(dtNew);
}

Info<< "maxDi = " << maxDiNow
    << "  deltaT = " << runTime.deltaTValue() << nl << endl;