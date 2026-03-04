scalar Tinf  = 300.0;

// updateConvectionVol.H
convectionLoss = dimensionedScalar("zero", dimPower/dimVolume/dimTemperature, 0.0);
forAll(mesh.cells(), cellI)
{
    scalar depth = zSurf - mesh.C()[cellI].z();

    if (depth >= 0.0 && depth <= hLoss)
    {
        convectionLoss[cellI] = hConv/hLoss; // W/m3/K
    }
}
convectionLoss.correctBoundaryConditions();