// Apparent heat capacity latent heat model. The liquid fraction ramps
// linearly through the mushy interval and contributes L*d(fL)/dT to cpEff.
const dimensionedScalar dTmush
(
    "dTmush",
    dimTemperature,
    max(SMALL, (Tl - Ts).value())
);

// Limit the phase fraction to [0, 1] so overshoots in T do not feed back into
// the latent-heat source term.
fL =
    min
    (
        max
        (
            (T - Ts)/dTmush,
            dimensionedScalar("zero", dimless, 0)
        ),
        dimensionedScalar("one", dimless, 1)
    );

const volScalarField mushyMask(pos(T - Ts)*pos(Tl - T));
cpEff = cp + mushyMask*(L/dTmush);
