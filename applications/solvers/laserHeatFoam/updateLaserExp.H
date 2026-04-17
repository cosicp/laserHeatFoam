Info<< "Updating laser\n" << endl;

// Legacy exponential Gaussian source fragment. The unified updateLaser.H
// handles this model in normal solver runs.
const scalar timeValue = runTime.time().timeOutputValue();
const vector laserPos = laserPositionSeries(timeValue);
const scalar laserPower = laserPowerSeries(timeValue);
const scalar laserSurfaceZ = laserPos.z();

const scalar A = readScalar(LaserProperties.lookup("absorptivity"));
const scalar rb = readScalar(LaserProperties.lookup("laserRadius"));
const scalar d = readScalar(LaserProperties.lookup("absorptionDepth"));
const scalar pi = constant::mathematical::pi;

Q = dimensionedScalar("Q", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0);

// Gaussian radial profile with exponential attenuation below the surface.
forAll(mesh.cells(), cellI)
{
    const vector cellPos = mesh.C()[cellI];
    const vector rVec = cellPos - laserPos;
    const scalar rSqr = Foam::sqr(rVec.x()) + Foam::sqr(rVec.y());
    const scalar depth = laserSurfaceZ - cellPos.z();

    if (depth > 0.0)
    {
        Q[cellI] =          
            (2.0*A*laserPower)/(pi*Foam::sqr(rb)*d)
           *Foam::exp(-2.0*rSqr/Foam::sqr(rb))
           *Foam::exp(-depth/d);
    }
    else
    {
        Q[cellI] = 0.0;
    }
}

scalar Pdep = 0.0;

// Integrated source power, reduced across processors.
forAll(Q, cellI)
{
    Pdep += Q[cellI]*mesh.V()[cellI];
}

reduce(Pdep, sumOp<scalar>());

Info<< "Deposited power = " << Pdep << " W" << nl << endl;
