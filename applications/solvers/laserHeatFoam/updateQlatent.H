fLold = fL;

forAll(mesh.cells(), cellI)
{
    scalar Tc = T[cellI];

    if (Tc <= Ts)      fL[cellI] = 0.0;
    else if (Tc >= Tl) fL[cellI] = 1.0;
    else               fL[cellI] = (Tc - Ts)/(Tl - Ts);
}
fL.correctBoundaryConditions();

volScalarField dfLdt
(
    IOobject("dfLdt", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0.0)
);

dfLdt = (fL - fLold)/runTime.deltaT();
dfLdt.correctBoundaryConditions();

volScalarField Qlatent
(
    IOobject("Qlatent", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimPower/dimVolume, 0.0)
);

Qlatent = -rho*L*dfLdt;
Qlatent.correctBoundaryConditions();