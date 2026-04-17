# Laser Model Tutorials

These cases are small variations of the existing titanium single-track setup and are meant to show how to select and parameterize the supported `laserModel` options in `constant/LaserProperties`.

Cases added:

- `laserModel_surfaceGaussian`
- `laserModel_cylindrical`
- `laserModel_exponential`
- `laserModel_ellipsoidal`
- `laserModel_goldak`
- `laserModel_conical`
- `laserModel_conicalExponential`
- `laserModel_cylindricalExponential`
- `laserModel_beerLambert`
- `laserModel_hybrid`

Each case reuses the same mesh, material properties, and scan path:

- scan line from `(120e-6 100e-6 300e-6)` to `(120e-6 700e-6 300e-6)`
- constant 195 W laser power during the track
- identical thermal and boundary-condition setup

Only `constant/LaserProperties` changes between cases. The parameter values are intended as working examples for the solver input format, not calibrated process parameters.

Run a case with:

```bash
cd tutorials/laserModel_exponential
./Allrun
```

Key model-specific inputs:

- `surfaceGaussian`: `laserRadius`
- `cylindrical`: `laserRadius`, `absorptionDepth`
- `exponential`: `laserRadius`, `absorptionDepth`
- `ellipsoidal`: `ax`, `ay`, `az`
- `goldak`: `goldakA`, `goldakB`, `goldakC1`, `goldakC2`, `goldakFrontFraction`, `goldakRearFraction`, `scanDirection`
- `conical`: `rTop`, `rBot`, `keyholeDepth`
- `conicalExponential`: `rTop`, `rBot`, `keyholeDepth`, `axialDecayLength`
- `cylindricalExponential`: `laserRadius`, `absorptionDepth`, `cylindricalHeight`
- `beerLambert`: `laserRadius`, `absorptionCoefficient`
- `hybrid`: `laserRadius`, `absorptionDepth`, `surfaceFraction`

Aliases supported by the solver but not duplicated here:

- `gaussianSurfaceFlux` behaves like `surfaceGaussian`
- `cylindricalGaussian` behaves like `cylindrical`
- `goldakDoubleEllipsoid` behaves like `goldak`

