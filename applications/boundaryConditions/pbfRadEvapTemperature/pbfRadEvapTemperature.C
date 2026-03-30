#include "pbfRadEvapTemperature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(pbfRadEvapTemperatureFvPatchScalarField, 0);

addToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    pbfRadEvapTemperatureFvPatchScalarField
);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

pbfRadEvapTemperatureFvPatchScalarField::
pbfRadEvapTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    kappaName_("kappaEff"),
    Tinf_(300.0),
    hConv_(0.0),
    epsilon_(0.7),
    sigmaSB_(constant::physicoChemical::sigma.value()),
    Tv_(3000.0),
    Tmax_(4000.0),
    enableEvap_(true),
    CP_(54000.0),
    CT_(50000.0),
    CM_(0.001),
    hv_(6.0e6),
    Th0_(663.0),
    cpEvap_(500.0)
{
    refValue()    = scalarField(p.size(), Tinf_);
    refGrad()     = scalarField(p.size(), 0.0);
    valueFraction() = scalarField(p.size(), 0.0);
}


pbfRadEvapTemperatureFvPatchScalarField::
pbfRadEvapTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    kappaName_(dict.lookupOrDefault<word>("kappa", "kappaEff")),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    hConv_(dict.lookupOrDefault<scalar>("hConv", 0.0)),
    epsilon_(readScalar(dict.lookup("epsilon"))),
    sigmaSB_
    (
        dict.lookupOrDefault<scalar>
        (
            "sigmaSB",
            constant::physicoChemical::sigma.value()
        )
    ),
    Tv_(readScalar(dict.lookup("Tv"))),
    Tmax_(readScalar(dict.lookup("Tmax"))),
    enableEvap_(dict.lookupOrDefault<Switch>("enableEvap", true)),
    CP_(readScalar(dict.lookup("CP"))),
    CT_(readScalar(dict.lookup("CT"))),
    CM_(readScalar(dict.lookup("CM"))),
    hv_(readScalar(dict.lookup("hv"))),
    Th0_(readScalar(dict.lookup("Th0"))),
    cpEvap_(readScalar(dict.lookup("cpEvap")))
{
    refValue()      = scalarField(p.size(), Tinf_);
    refGrad()       = scalarField(p.size(), 0.0);
    valueFraction() = scalarField(p.size(), 0.0);

    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );

    Info<< "*** pbfRadEvapTemperature BC loaded successfully ***" << endl;
}


pbfRadEvapTemperatureFvPatchScalarField::
pbfRadEvapTemperatureFvPatchScalarField
(
    const pbfRadEvapTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    kappaName_(ptf.kappaName_),
    Tinf_(ptf.Tinf_),
    hConv_(ptf.hConv_),
    epsilon_(ptf.epsilon_),
    sigmaSB_(ptf.sigmaSB_),
    Tv_(ptf.Tv_),
    Tmax_(ptf.Tmax_),
    enableEvap_(ptf.enableEvap_),
    CP_(ptf.CP_),
    CT_(ptf.CT_),
    CM_(ptf.CM_),
    hv_(ptf.hv_),
    Th0_(ptf.Th0_),
    cpEvap_(ptf.cpEvap_)
{}


pbfRadEvapTemperatureFvPatchScalarField::
pbfRadEvapTemperatureFvPatchScalarField
(
    const pbfRadEvapTemperatureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    kappaName_(ptf.kappaName_),
    Tinf_(ptf.Tinf_),
    hConv_(ptf.hConv_),
    epsilon_(ptf.epsilon_),
    sigmaSB_(ptf.sigmaSB_),
    Tv_(ptf.Tv_),
    Tmax_(ptf.Tmax_),
    enableEvap_(ptf.enableEvap_),
    CP_(ptf.CP_),
    CT_(ptf.CT_),
    CM_(ptf.CM_),
    hv_(ptf.hv_),
    Th0_(ptf.Th0_),
    cpEvap_(ptf.cpEvap_)
{}


pbfRadEvapTemperatureFvPatchScalarField::
pbfRadEvapTemperatureFvPatchScalarField
(
    const pbfRadEvapTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    kappaName_(ptf.kappaName_),
    Tinf_(ptf.Tinf_),
    hConv_(ptf.hConv_),
    epsilon_(ptf.epsilon_),
    sigmaSB_(ptf.sigmaSB_),
    Tv_(ptf.Tv_),
    Tmax_(ptf.Tmax_),
    enableEvap_(ptf.enableEvap_),
    CP_(ptf.CP_),
    CT_(ptf.CT_),
    CM_(ptf.CM_),
    hv_(ptf.hv_),
    Th0_(ptf.Th0_),
    cpEvap_(ptf.cpEvap_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pbfRadEvapTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>
        (
            internalField().name()
        );

    const fvPatchScalarField& kappap =
        patch().lookupPatchField<volScalarField, scalar>(kappaName_);

    scalarField& rGrad = refGrad();
    scalarField& vFrac = valueFraction();
    scalarField& rVal  = refValue();

    forAll(Tp, faceI)
    {
        const scalar Tw = Tp[faceI];
        const scalar Tc = min(Tw, Tmax_);

        // Radiation
        const scalar qRad = epsilon_*sigmaSB_*(pow4(Tw) - pow4(Tinf_));

        // Convection
        const scalar qConv = hConv_*(Tw - Tinf_);

        // Evaporation
        scalar qEvap = 0.0;
        if (enableEvap_ && Tc > Tv_)
        {
            const scalar mdot =
                0.82*CP_
               *exp(-CT_*(1.0/Tc - 1.0/Tv_))
               *sqrt(CM_/Tc);

            qEvap = mdot*(hv_ + cpEvap_*(Tc - Th0_));
        }

        const scalar qTot = qRad + qConv + qEvap;

        // Pure gradient mode: valueFraction = 0 => -k dT/dn = qTot
        rGrad[faceI] = -qTot/max(kappap[faceI], SMALL);
        rVal[faceI]  = Tinf_;
        vFrac[faceI] = 0.0;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void pbfRadEvapTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("kappa",      kappaName_);
    os.writeEntry("Tinf",       Tinf_);
    os.writeEntry("hConv",      hConv_);
    os.writeEntry("epsilon",    epsilon_);
    os.writeEntry("sigmaSB",    sigmaSB_);
    os.writeEntry("Tv",         Tv_);
    os.writeEntry("Tmax",       Tmax_);
    os.writeEntry("enableEvap", enableEvap_);
    os.writeEntry("CP",         CP_);
    os.writeEntry("CT",         CT_);
    os.writeEntry("CM",         CM_);
    os.writeEntry("hv",         hv_);
    os.writeEntry("Th0",        Th0_);
    os.writeEntry("cpEvap",     cpEvap_);

    writeEntry("value", os);
}

} // End namespace Foam
