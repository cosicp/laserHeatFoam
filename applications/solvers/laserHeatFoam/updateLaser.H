// --- LPBF parameters ---
scalar P  = 200.0;        // W
scalar A  = 0.4;          // absorptivity
scalar rb = 40e-6;        // beam radius [m]
scalar d  = 30e-6;        // absorption depth [m]

vector v = vector(0.01, 0.01, 0.0); // laser velocity [m/s]
vector pos0 = vector(0.0, 0.0, 0.0); // initial laser position

// --- Compute current laser position ---
vector laserPos = pos0 + v * runTime.timeOutputValue();

// scalar deltaT = 0.001;
// scalar Q;

volScalarField Q (0);
// Loop over all cells
forAll(T.internalField(), celli)
{
    vector cellPos = mesh.C()[celli];   // cell center
    vector rVec = cellPos - laserPos;
    
    // Radial distance in xy-plane
    scalar r = Foam::sqrt(rVec.x()*rVec.x() + rVec.y()*rVec.y());
    scalar z = rVec.z();                 // depth

    // Gaussian volumetric heat source
     Q[celli] = (2.0 * A * P) / (M_PI * Foam::sqrt(rb) * d)
               * Foam::exp(-2.0 * Foam::sqrt(r)/Foam::sqrt(rb))
               * Foam::exp(-z/d);
}