// --- LPBF parameters ---
// --- Compute current laser position ---

Info << "Updating laser\n" << endl;
vector laserPos = laserPositionSeries(runTime.time().timeOutputValue());
scalar laserPower = laserPowerSeries(runTime.time().timeOutputValue());
scalar z_layer = laserPos.z();                 // depth
Info << "z_layer: " << z_layer << endl;

Q = dimensionedScalar("Q", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0);

// Info << Q<< endl;
// Loop over all 

forAll(mesh.cells(), cellI)
{
    // const cell& c = mesh.cells()[cellI]; // 'c' is the current cell
    // Info << "here1"<<endl;
    vector cellPos = mesh.C()[cellI];   // cell center
    vector rVec = cellPos - laserPos;
    scalar r = Foam::sqrt(rVec.x()*rVec.x() + rVec.y()*rVec.y());
    // Gaussian volumetric heat source
    if(cellPos.z() < z_layer)
     Q[cellI] = (2.0 * A * laserPower) / (M_PI * Foam::sqr(rb) * d)
               * Foam::exp(-2.0 * Foam::sqr(r)/Foam::sqr(rb))
               * Foam::exp(-(z_layer - cellPos.z())/d);
}

// Info << "Q After update\n" << endl;
// Info << Q<< endl;
