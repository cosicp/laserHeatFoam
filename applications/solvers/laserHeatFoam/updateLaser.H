// --- LPBF parameters ---
// --- Compute current laser position ---

Info << "Updating laser\n" << endl;

vector laserPos = laserPositionSeries(runTime.time().timeOutputValue());
scalar laserPower = laserPowerSeries(runTime.time().timeOutputValue());
scalar z_layer = laserPos.z(); // depth

Q = dimensionedScalar("Q", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0);

forAll(mesh.cells(), cellI)
{
    vector cellPos = mesh.C()[cellI]; // cell center
    vector rVec = cellPos - laserPos;
    scalar r = Foam::sqrt(rVec.x() * rVec.x() + rVec.y() * rVec.y());
    scalar Qr = (2.0 * A * laserPower) / (M_PI * Foam::sqr(rb) * d)
                * Foam::exp(-2.0 * Foam::sqr(r) / Foam::sqr(rb));

    if (laserModel == "exponential")
    {
        // Exponential decay in z below the surface
        Q[cellI] = (cellPos.z() < z_layer)
            ? Qr * Foam::exp(-(z_layer - cellPos.z()) / d)
            : 0.0;
    }
    else // "cylindrical" (default): uniform within depth window d
    {
        Q[cellI] = (cellPos.z() < z_layer && cellPos.z() > z_layer - d)
            ? Qr
            : 0.0;
    }
}

// Info << "Q After update\n" << endl;
// Info << Q<< endl;

scalar Pdep = 0.0;

forAll(Q, cellI)
{
    Pdep += Q[cellI]*mesh.V()[cellI];
}

reduce(Pdep, sumOp<scalar>());
Info<< "Deposited power = " << Pdep << " W" << nl << endl;
