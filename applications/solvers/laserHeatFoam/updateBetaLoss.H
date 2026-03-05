// updateRadConvVol.H
const scalar hConv  = 19;                   // W/m2/K
const scalar epsRad = 0.5;                  // emissivity
const scalar sigma  = 5.670374419e-8;       // Stefan-Boltzmann

// scalar zSurf = z_layer;     // change to your actual top surface z
// scalar hLoss = 20e-6;   // smear thickness

betaLoss = dimensionedScalar("zero", dimPower/dimVolume/dimTemperature, 0.0);

forAll(mesh.cells(), cellI)
{
    scalar depth = z_layer - mesh.C()[cellI].z();

    if (depth >= 0.0 && depth <= d)
    {
        scalar Tc = T[cellI];
        // linearized radiation coefficient [W/m2/K]
        scalar hRad = 4.0*epsRad*sigma*Foam::pow3(Tc);

        betaLoss[cellI] = (hConv + hRad)/d; // W/m3/K
    }
}
betaLoss.correctBoundaryConditions();