// updateThermo.H
// Update temperature-dependent thermophysical properties each timestep.
//
// Piecewise polynomial fits for Ti-6Al-4V from:
//   Mohammadkamal & Caiazzo (2025), Table 1  [Paper 1]
//
// rho [kg/m3]:
//   T < Tl  : rho0s + rho1s * T
//   T >= Tl : rho0m + rho1m * T
//
// cp [J/kg/K]:
//   T < T1        : cp0a + cp1a * T
//   T1 <= T < Tl  : cp0b + cp1b * T
//   T >= Tl       : cpL
//
// k solid/melt [W/m/K]:
//   T < T1        : k0a + k1a * T
//   T1 <= T < Tl  : k0b + k1b * T
//   T >= Tl       : k0m + k1m * T
//
// kp (powder) remains constant.

{
    // Breakpoint and liquidus temperatures
    const scalar T1val  = thermoDict.lookupOrDefault<scalar>("T1",    1268.0);
    const scalar Tlval  = Tl.value();   // liquidus [K]

    // rho coefficients
    const scalar rho0s = thermoDict.lookupOrDefault<scalar>("rho0s",  4557.892);
    const scalar rho1s = thermoDict.lookupOrDefault<scalar>("rho1s",    -0.154);
    const scalar rho0m = thermoDict.lookupOrDefault<scalar>("rho0m",  5227.640);
    const scalar rho1m = thermoDict.lookupOrDefault<scalar>("rho1m",    -0.680);

    // cp coefficients
    const scalar cp0a  = thermoDict.lookupOrDefault<scalar>("cp0a",   483.04);
    const scalar cp1a  = thermoDict.lookupOrDefault<scalar>("cp1a",     0.215);
    const scalar cp0b  = thermoDict.lookupOrDefault<scalar>("cp0b",   412.70);
    const scalar cp1b  = thermoDict.lookupOrDefault<scalar>("cp1b",     0.18);
    const scalar cpL   = thermoDict.lookupOrDefault<scalar>("cpL",    831.0);

    // k (solid/melt) coefficients
    const scalar k0a   = thermoDict.lookupOrDefault<scalar>("k0a",     1.250);
    const scalar k1a   = thermoDict.lookupOrDefault<scalar>("k1a",     0.015);
    const scalar k0b   = thermoDict.lookupOrDefault<scalar>("k0b",     3.15);
    const scalar k1b   = thermoDict.lookupOrDefault<scalar>("k1b",     0.012);
    const scalar k0m   = thermoDict.lookupOrDefault<scalar>("k0m",   -12.752);
    const scalar k1m   = thermoDict.lookupOrDefault<scalar>("k1m",     0.024);

    forAll(T, cellI)
    {
        const scalar Tc = T[cellI];

        // Density
        rho[cellI] = (Tc < Tlval)
            ? rho0s + rho1s*Tc
            : rho0m + rho1m*Tc;

        // Specific heat
        cp[cellI] = (Tc < T1val)  ? cp0a + cp1a*Tc
                  : (Tc < Tlval)  ? cp0b + cp1b*Tc
                  : cpL;

        // Solid-phase conductivity (used weighted by rs in mushy zone)
        ks[cellI] = (Tc < T1val)  ? k0a + k1a*Tc
                  : (Tc < Tlval)  ? k0b + k1b*Tc
                  : k0m + k1m*Tc;

        // Melt-phase conductivity (used weighted by rm above solidus)
        km[cellI] = (Tc < T1val)  ? k0a + k1a*Tc
                  : (Tc < Tlval)  ? k0b + k1b*Tc
                  : k0m + k1m*Tc;
    }
}
