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

    //- Ambient temperature
    scalar Tinf_;

    //- Optional convection coefficient [W/m2/K]
    scalar hConv_;

    //- Emissivity
    scalar epsilon_;

    //- Stefan-Boltzmann constant
    scalar sigmaSB_;

    //- Boiling temperature
    scalar Tv_;

    //- Temperature clip for evaporation law
    scalar Tmax_;

    //- Enable evaporation term
    Switch enableEvap_;

    //- Recoil pressure factor
    scalar CP_;

    //- Recoil pressure temperature factor
    scalar CT_;

    //- Heat loss temperature factor
    scalar CM_;

    //- Latent heat of evaporation
    scalar hv_;

    //- Enthalpy reference temperature
    scalar Th0_;

    //- Specific heat used inside qEvap enthalpy term
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