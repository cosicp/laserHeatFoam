// --- LPBF parameters ---
// --- Compute current laser position ---

Info << "Updating laser\n" << endl;
vector laserPos = laserPositionSeries(runTime.time().timeOutputValue());
scalar laserPower = laserPowerSeries(runTime.time().timeOutputValue());

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
    scalar z = rVec.z();                 // depth
    // Gaussian volumetric heat source
     Q[cellI] = (2.0 * A * laserPower) / (M_PI * Foam::sqr(rb) * d)
               * Foam::exp(-2.0 * Foam::sqr(r)/Foam::sqr(rb))
               * Foam::exp(-(200e-6-z)/d);
}
// Info << "Q After update\n" << endl;
// Info << Q<< endl;
