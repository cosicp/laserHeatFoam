// updateRadConvVol.H
scalar zSurf = laserPos.z();                 // depth
const scalar hConv  = 25.0;                 // W/m2/K
const scalar epsRad = 0.5;                  // emissivity
const scalar sigma  = 5.670374419e-8;       // Stefan-Boltzmann
scalar hLoss = 30e-6;   // smear thickness

radiationLoss = dimensionedScalar("zero", dimPower/dimVolume/dimTemperature, 0.0);

forAll(mesh.cells(), cellI)
{
    scalar depth = zSurf - mesh.C()[cellI].z();

    if (depth >= 0.0 && depth <= hLoss)
    {
        scalar Tc = T[cellI];
        // linearized radiation coefficient [W/m2/K]
        scalar hRad = 4.0*epsRad*sigma*Foam::pow3(Tc);
        radiationLoss[cellI] = hRad/hLoss; // W/m3/K
    }
}

radiationLoss.correctBoundaryConditions();