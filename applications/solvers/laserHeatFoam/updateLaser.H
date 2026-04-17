Info<< "Updating laser\n" << endl;

// Interpolate scan data to the current solver time.
const scalar timeValue = runTime.time().timeOutputValue();
const vector laserPos = laserPositionSeries(timeValue);
const scalar laserPower = laserPowerSeries(timeValue);
const scalar laserSurfaceZ = laserPos.z();

// Common laser input. AP is the absorbed laser power.
const scalar A = readScalar(LaserProperties.lookup("absorptivity"));
const scalar AP = A*laserPower;

// Gaussian beam and depth-scale parameters.
const scalar rb = LaserProperties.lookupOrDefault<scalar>("laserRadius", 0.0);
const scalar d = LaserProperties.lookupOrDefault<scalar>("absorptionDepth", 0.0);
const scalar axialDecayLength =
    LaserProperties.lookupOrDefault<scalar>("axialDecayLength", d);
const scalar absorptionCoeff =
    LaserProperties.lookupOrDefault<scalar>("absorptionCoefficient", 0.0);
const scalar surfaceFraction =
    LaserProperties.lookupOrDefault<scalar>("surfaceFraction", 0.5);

// Ellipsoidal and Goldak model dimensions.
const scalar ax = LaserProperties.lookupOrDefault<scalar>("ax", rb);
const scalar ay = LaserProperties.lookupOrDefault<scalar>("ay", rb);
const scalar az = LaserProperties.lookupOrDefault<scalar>("az", d);

const scalar goldakA = LaserProperties.lookupOrDefault<scalar>("goldakA", rb);
const scalar goldakB = LaserProperties.lookupOrDefault<scalar>("goldakB", rb);
const scalar goldakC1 = LaserProperties.lookupOrDefault<scalar>("goldakC1", d);
const scalar goldakC2 = LaserProperties.lookupOrDefault<scalar>("goldakC2", d);
const scalar goldakFrontFraction =
    LaserProperties.lookupOrDefault<scalar>("goldakFrontFraction", 1.0);
const scalar goldakRearFraction =
    LaserProperties.lookupOrDefault<scalar>("goldakRearFraction", 1.0);

// Conical/keyhole model dimensions.
const scalar rTop = LaserProperties.lookupOrDefault<scalar>("rTop", rb);
const scalar rBot = LaserProperties.lookupOrDefault<scalar>("rBot", rb);
const scalar H = LaserProperties.lookupOrDefault<scalar>("keyholeDepth", d);
const scalar cylindricalHeight =
    LaserProperties.lookupOrDefault<scalar>("cylindricalHeight", H);

// Direction used to split the Goldak front and rear ellipsoids.
vector scanDirection =
    LaserProperties.lookupOrDefault<vector>("scanDirection", vector(1, 0, 0));
const scalar scanDirectionMag = mag(scanDirection);
if (scanDirectionMag > SMALL)
{
    scanDirection /= scanDirectionMag;
}
else
{
    scanDirection = vector(1, 0, 0);
}

const scalar pi = constant::mathematical::pi;
const scalar piPowThreeHalves = pi*Foam::sqrt(pi);

// Keep model-specific parameter validation close to the model equations.
const auto requirePositive =
    [&](const scalar value, const word& parameterName) -> scalar
    {
        if (value <= SMALL)
        {
            FatalErrorInFunction
                << "Laser model '" << laserModel << "' requires '"
                << parameterName << "' > 0, but got " << value
                << exit(FatalError);
        }

        return value;
    };

// Estimate local vertical cell size for surface-flux-style source terms.
const auto computeCellThickness =
    [&](const label cellI) -> scalar
    {
        const cell& c = mesh.cells()[cellI];
        scalar zMin = GREAT;
        scalar zMax = -GREAT;

        forAll(c, faceI)
        {
            const face& f = mesh.faces()[c[faceI]];
            forAll(f, pointI)
            {
                const scalar z = mesh.points()[f[pointI]].z();
                zMin = min(zMin, z);
                zMax = max(zMax, z);
            }
        }

        return max(zMax - zMin, SMALL);
    };

// Linear cone radius, measured downward from the laser surface.
const auto linearRadiusAtDepth =
    [&](const scalar depth) -> scalar
    {
        const scalar safeH = requirePositive(H, "keyholeDepth");
        return rTop - (rTop - rBot)*depth/safeH;
    };

Q = dimensionedScalar("Q", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0);

// Fail early for misspelled model names instead of silently using no source.
if
(
    laserModel != "gaussianSurfaceFlux"
 && laserModel != "surfaceGaussian"
 && laserModel != "cylindrical"
 && laserModel != "cylindricalGaussian"
 && laserModel != "exponential"
 && laserModel != "ellipsoidal"
 && laserModel != "goldakDoubleEllipsoid"
 && laserModel != "goldak"
 && laserModel != "conical"
 && laserModel != "conicalExponential"
 && laserModel != "cylindricalExponential"
 && laserModel != "beerLambert"
 && laserModel != "hybrid"
)
{
    FatalErrorInFunction
        << "Unknown laserModel '" << laserModel << "'."
        << " Supported models are: gaussianSurfaceFlux, cylindrical,"
        << " exponential, ellipsoidal, goldakDoubleEllipsoid,"
        << " conical, conicalExponential, cylindricalExponential,"
        << " beerLambert, hybrid"
        << exit(FatalError);
}

forAll(mesh.cells(), cellI)
{
    const vector cellPos = mesh.C()[cellI];
    const vector rel = cellPos - laserPos;
    const scalar rSqr = Foam::sqr(rel.x()) + Foam::sqr(rel.y());
    const scalar r = Foam::sqrt(rSqr);
    const scalar depth = laserSurfaceZ - cellPos.z();

    scalar sourceValue = 0.0;

    // Surface Gaussian flux converted into a volumetric source over one cell.
    if (laserModel == "gaussianSurfaceFlux" || laserModel == "surfaceGaussian")
    {
        const scalar safeRb = requirePositive(rb, "laserRadius");
        const scalar dz = computeCellThickness(cellI);

        if (depth >= 0.0 && depth < dz)
        {
            sourceValue =
                (2.0*AP)/(pi*Foam::sqr(safeRb)*dz)
               *Foam::exp(-2.0*rSqr/Foam::sqr(safeRb));
        }
    }
    else if (laserModel == "cylindrical" || laserModel == "cylindricalGaussian")
    {
        // Gaussian radial profile distributed over a finite absorption depth.
        const scalar safeRb = requirePositive(rb, "laserRadius");
        const scalar safeD = requirePositive(d, "absorptionDepth");

        if (depth >= 0.0 && depth <= safeD)
        {
            sourceValue =
                (2.0*AP)/(pi*Foam::sqr(safeRb)*safeD)
               *Foam::exp(-2.0*rSqr/Foam::sqr(safeRb));
        }
    }
    else if (laserModel == "exponential")
    {
        // Gaussian beam with exponential attenuation below the surface.
        const scalar safeRb = requirePositive(rb, "laserRadius");
        const scalar safeD = requirePositive(d, "absorptionDepth");

        if (depth >= 0.0)
        {
            sourceValue =
                (2.0*AP)/(pi*Foam::sqr(safeRb)*safeD)
               *Foam::exp(-2.0*rSqr/Foam::sqr(safeRb))
               *Foam::exp(-depth/safeD);
        }
    }
    else if (laserModel == "ellipsoidal")
    {
        // Single ellipsoidal volumetric source.
        const scalar safeAx = requirePositive(ax, "ax");
        const scalar safeAy = requirePositive(ay, "ay");
        const scalar safeAz = requirePositive(az, "az");

        if (depth >= 0.0)
        {
            const scalar Q0 = (2.0*AP)/(piPowThreeHalves*safeAx*safeAy*safeAz);
            sourceValue =
                Q0
               *Foam::exp
                (
                    -Foam::sqr(rel.x())/Foam::sqr(safeAx)
                    -Foam::sqr(rel.y())/Foam::sqr(safeAy)
                    -Foam::sqr(depth)/Foam::sqr(safeAz)
                );
        }
    }
    else if (laserModel == "goldakDoubleEllipsoid" || laserModel == "goldak")
    {
        // Goldak double ellipsoid, split by the scan direction.
        const scalar safeA = requirePositive(goldakA, "goldakA");
        const scalar safeB = requirePositive(goldakB, "goldakB");
        const scalar safeC1 = requirePositive(goldakC1, "goldakC1");
        const scalar safeC2 = requirePositive(goldakC2, "goldakC2");

        if (depth >= 0.0)
        {
            const scalar xi = rel & scanDirection;
            const vector transverse = rel - xi*scanDirection;
            const scalar etaSqr = magSqr(transverse);

            if (xi >= 0.0)
            {
                const scalar frontPrefactor =
                    (6.0*Foam::sqrt(3.0)*goldakFrontFraction*AP)
                   /(piPowThreeHalves*safeA*safeB*safeC1);

                sourceValue =
                    frontPrefactor
                   *Foam::exp
                    (
                        -3.0*Foam::sqr(xi)/Foam::sqr(safeA)
                        -3.0*etaSqr/Foam::sqr(safeB)
                        -3.0*Foam::sqr(depth)/Foam::sqr(safeC1)
                    );
            }
            else
            {
                const scalar rearPrefactor =
                    (6.0*Foam::sqrt(3.0)*goldakRearFraction*AP)
                   /(piPowThreeHalves*safeA*safeB*safeC2);

                sourceValue =
                    rearPrefactor
                   *Foam::exp
                    (
                        -3.0*Foam::sqr(xi)/Foam::sqr(safeA)
                        -3.0*etaSqr/Foam::sqr(safeB)
                        -3.0*Foam::sqr(depth)/Foam::sqr(safeC2)
                    );
            }
        }
    }
    else if (laserModel == "conical")
    {
        // Conical source with Gaussian radial weighting at each depth.
        requirePositive(rTop, "rTop");
        requirePositive(rBot, "rBot");
        const scalar safeH = requirePositive(H, "keyholeDepth");
        const scalar denom =
            requirePositive
            (
                Foam::sqr(rTop) + rTop*rBot + Foam::sqr(rBot),
                "rTop^2 + rTop*rBot + rBot^2"
            );

        if (depth >= 0.0 && depth <= safeH)
        {
            const scalar localRadius =
                max(linearRadiusAtDepth(depth), SMALL);
            const scalar Q0 = (9.0*AP)/(pi*safeH*denom);

            sourceValue =
                Q0*Foam::exp(-3.0*rSqr/Foam::sqr(localRadius));
        }
    }
    else if (laserModel == "conicalExponential")
    {
        // Conical source with an additional axial exponential decay.
        requirePositive(rTop, "rTop");
        requirePositive(rBot, "rBot");
        const scalar safeH = requirePositive(H, "keyholeDepth");
        const scalar safeLambda =
            requirePositive(axialDecayLength, "axialDecayLength");
        const scalar denom =
            requirePositive
            (
                Foam::sqr(rTop) + rTop*rBot + Foam::sqr(rBot),
                "rTop^2 + rTop*rBot + rBot^2"
            );

        if (depth >= 0.0 && depth <= safeH)
        {
            const scalar localRadius =
                max(linearRadiusAtDepth(depth), SMALL);
            const scalar Q0 = (9.0*AP)/(pi*safeH*denom);

            sourceValue =
                Q0
               *Foam::exp(-3.0*rSqr/Foam::sqr(localRadius))
               *Foam::exp(-depth/safeLambda);
        }
    }
    else if (laserModel == "cylindricalExponential")
    {
        // Top-hat radial profile with exponential axial attenuation.
        const scalar safeRb = requirePositive(rb, "laserRadius");
        const scalar safeD = requirePositive(d, "absorptionDepth");
        const scalar safeH =
            requirePositive
            (
                cylindricalHeight > SMALL ? cylindricalHeight : safeD,
                "cylindricalHeight"
            );

        if (r <= safeRb && depth >= 0.0 && depth <= safeH)
        {
            sourceValue =
                (AP)/(pi*Foam::sqr(safeRb)*safeD)
               *Foam::exp(-depth/safeD);
        }
    }
    else if (laserModel == "beerLambert")
    {
        // Beer-Lambert absorption using alpha*I(z) as the volumetric source.
        const scalar safeRb = requirePositive(rb, "laserRadius");
        const scalar safeAlpha =
            requirePositive(absorptionCoeff, "absorptionCoefficient");

        if (depth >= 0.0)
        {
            const scalar I0 = (2.0*AP)/(pi*Foam::sqr(safeRb));
            sourceValue =
                safeAlpha*I0
               *Foam::exp(-2.0*rSqr/Foam::sqr(safeRb))
               *Foam::exp(-safeAlpha*depth);
        }
    }
    else if (laserModel == "hybrid")
    {
        // Blend a surface Gaussian contribution with a volumetric tail.
        const scalar safeRb = requirePositive(rb, "laserRadius");
        const scalar safeD = requirePositive(d, "absorptionDepth");
        const scalar clippedSurfaceFraction =
            min(max(surfaceFraction, scalar(0)), scalar(1));
        const scalar dz = computeCellThickness(cellI);

        if (depth >= 0.0)
        {
            const scalar gaussianRadial =
                Foam::exp(-2.0*rSqr/Foam::sqr(safeRb));

            if (depth < dz)
            {
                sourceValue +=
                    clippedSurfaceFraction
                   *(2.0*AP)/(pi*Foam::sqr(safeRb)*dz)
                   *gaussianRadial;
            }

            sourceValue +=
                (1.0 - clippedSurfaceFraction)
               *(2.0*AP)/(pi*Foam::sqr(safeRb)*safeD)
               *gaussianRadial
               *Foam::exp(-depth/safeD);
        }
    }

    Q[cellI] = sourceValue;
}

scalar Pdep = 0.0;

// Integrated source power, reduced across processors.
forAll(Q, cellI)
{
    Pdep += Q[cellI]*mesh.V()[cellI];
}

reduce(Pdep, sumOp<scalar>());

Info<< "Laser model = " << laserModel << nl
    << "Deposited power = " << Pdep << " W";

if (AP > SMALL)
{
    Info<< " (" << 100.0*Pdep/AP << "% of absorbed power)";
}

Info<< nl << endl;
