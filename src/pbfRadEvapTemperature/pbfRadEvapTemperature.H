#ifndef pbfRadEvapTemperatureFvPatchScalarField_H
#define pbfRadEvapTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "Switch.H"

namespace Foam
{

class pbfRadEvapTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    //- Name of conductivity field
    word kappaName_;

    //- Ambient temperature [K], read from thermophysicalProperties
    scalar Tinf_;

    //- Optional convection coefficient [W/m2/K], read from thermophysicalProperties
    scalar hConv_;

    //- Emissivity [-], read from thermophysicalProperties
    scalar epsilon_;

    //- Stefan-Boltzmann constant [W/m2/K4]
    scalar sigmaSB_;

    //- Boiling/evaporation onset temperature [K], read from thermophysicalProperties
    scalar Tv_;

    //- Temperature clip for evaporation law [K], read from thermophysicalProperties
    scalar Tmax_;

    //- Enable evaporation heat loss term
    Switch enableEvap_;

    //- Pre-exponential pressure factor [Pa]
    scalar CP_;

    //- Exponential temperature factor [K]
    scalar CT_;

    //- Molecular mass factor [kg/mol]
    scalar CM_;

    //- Latent heat of evaporation [J/kg]
    scalar hv_;

    //- Enthalpy reference temperature [K]
    scalar Th0_;

    //- Specific heat inside qEvap enthalpy term [J/kg/K]
    scalar cpEvap_;


public:

    TypeName("pbfRadEvapTemperature");

    // Constructors

    pbfRadEvapTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    pbfRadEvapTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    pbfRadEvapTemperatureFvPatchScalarField
    (
        const pbfRadEvapTemperatureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    pbfRadEvapTemperatureFvPatchScalarField
    (
        const pbfRadEvapTemperatureFvPatchScalarField&
    );

    pbfRadEvapTemperatureFvPatchScalarField
    (
        const pbfRadEvapTemperatureFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new pbfRadEvapTemperatureFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new pbfRadEvapTemperatureFvPatchScalarField(*this, iF)
        );
    }

    // Member functions

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
