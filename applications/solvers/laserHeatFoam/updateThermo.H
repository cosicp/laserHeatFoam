// Update temperature-dependent thermophysical properties.
//
// Piecewise polynomial fits for Ti-6Al-4V from:
//   Mohammadkamal and Caiazzo (2025), Table 1
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

forAll(T, cellI)
{
    const scalar Tc = T[cellI];

    // Density is evaluated for solid/melt material, then reduced in cells
    // that still contain unconsolidated powder.
    const scalar rhoSolid = (Tc < Tlval)
        ? rho0s + rho1s*Tc
        : rho0m + rho1m*Tc;
    rho[cellI] = rhoSolid*(1.0 - phiVal*rp[cellI]);

    // Specific heat.
    cp[cellI] = (Tc < T1val)  ? cp0a + cp1a*Tc
              : (Tc < Tlval)  ? cp0b + cp1b*Tc
              : cpL;

    // Solid and melt conductivities are stored separately so the phase
    // fraction update can form k = rp*kp + rm*km + rs*ks.
    ks[cellI] = (Tc < T1val)  ? k0a + k1a*Tc
              : (Tc < Tlval)  ? k0b + k1b*Tc
              : k0m + k1m*Tc;

    km[cellI] = (Tc < T1val)  ? k0a + k1a*Tc
              : (Tc < Tlval)  ? k0b + k1b*Tc
              : k0m + k1m*Tc;
}
