// Volumetric approximation of a surface convection sink over hLoss.
convectionLoss = dimensionedScalar("zero", dimPower/dimVolume/dimTemperature, 0.0);

forAll(mesh.cells(), cellI)
{
    const scalar depth = zSurf - mesh.C()[cellI].z();

    if (depth >= 0.0 && depth <= hLoss)
    {
        convectionLoss[cellI] = hConv/hLoss;
    }
}

convectionLoss.correctBoundaryConditions();
