/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Author
    Tarang Garg
    gargtarang@gmail.com

Class
    Foam::compressible::ebeamTurbulentHeatFluxFvPatchScalarField

Description
    Heat Boundary condition to specify electron beam heating source or laser
    beam heating source. Evaporative heat loss and radiation condition can
    also be specified by specifying relevant parameters. Multiple heating
    sources can be specified using multiple dictionaries.

    The thermal conductivity \c kappa can either be retrieved from various
    possible sources, as detailed in the class temperatureCoupledBase.

Usage
    \table
        Property     | Description             | Required    | Default value
        heatSource   | 'circle' or 'ellipse' or 'strip'| yes |
        heatProfile  | 'flat' or 'gaussian | yes | 'flat'
        Q	     | Heat Absorbed [W] | yes |
        radius       | radius of circular beam [m] (sigma in case of gaussian) | yes |
        semiMajorAxis| SemiMajorAxis X-direction of ellipse beam [m] (sigmaX in case of gaussian) | yes |
        semiMinorAxis| SemiMinorAxis Z-direction of ellipse beam [m] (sigmaY in case of gaussian) | yes |
        length	     | length X-direction of strip beam [m] | yes |
        sigmaX       | sigmaX in X-direction of strip [m] | no | 1e-9
        breadth	     | breadth Z-direction of strip beam [m] | yes |
        sigmaZ       | sigmaX in Z-direction of strip [m] | no | 1e-9
        Qr           | name of the radiative flux field | no | none
        value        | initial temperature value | no | calculated
        gradient     | initial gradient value | no | 0.0
        kappaMethod  | inherited from temperatureCoupledBase | inherited |
        kappa        | inherited from temperatureCoupledBase | inherited |
        amplitude    | amplitude of oscillation [vector] | no | (0 0 0)  |
        frequency    | frequency of oscillation [Hz] [scalar] | no | 0.0 |
        evaporation  | Switch for evaporation | optional |off
        A            | A for log10(p) = -A/T+B |optional|1e+10
        B            | B for log10(p) = -A/T+B |optional|0
        massNumber   | mass number of evaporated material | 1
        vapourisationEnthalpy | enthalpy of vaporisation[J/kg] |0
        units        | units of pressure in above equation |Pa
    \endtable

    Note: If needed, both 'value' and 'gradient' must be defined to be used.

    Example usage:
    \verbatim
    hotSurface
    {
        type            compressible::ebeamTurbulentHeatFlux;
        sources         3;
        source1
        {
            heatSource      strip;
            heatProfile     flat;
            Q               10000;
            centre          (-0.015 0 0);
            length          0.01;
            breadth         0.01;
            amplitude       (0.02 0 0);
            frequency       50;
        }
        source2
        {
            heatSource      circle;
            heatProfile     flat;
            Q               15000;
            radius          0.005;
            centre          (0 0 0);
            amplitude       (0.01 0 0);
            frequency       0;
        }
        source3
        {
            heatSource      strip;
            heatProfile     flat;
            Q               10000;
            centre          (0.015 0 0);
            length          0.01;
            breadth         0.01;
            amplitude       (-0.02 0 0);
            frequency       25;
        }
        kappaMethod     fluidThermo;
        kappa           none;
        alphaAni        Anialpha;
        evaporation     on;
        A               16450;
        B               12.35;
        massNumber      27;
        vapourisationEnthalpy   1.15e7;
        units           mmHg;
        value           uniform 298;
    }

    \endverbatim

See also
    Foam::temperatureCoupledBase

SourceFiles
    ebeamTurbulentHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef ebeamTurbulentHeatFluxFvPatchScalarFields_H
#define ebeamTurbulentHeatFluxFvPatchScalarFields_H

#include "fixedGradientFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
      Class ebeamTurbulentHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ebeamTurbulentHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    // Data types

        //- Enumeration listing the possible heat source input modes
        enum heatSourceType
        {
            hsCircle,
            hsEllipse,
            hsStrip
        };
        
        //- Enumeration listing the possible heat flux profile modes
        enum heatFluxProfileType
        {
            hfpGaussian,
            hfpFlat
        };
        
        //- Enumeration listing the possible units name for evaporation
        enum unitsType
        {
            unPa,
            unatm,
            unmmHg,
            uncmHg,
            unbar,
            unkPa,
            unMPa
        };


private:

    // Private data

        //- Heat source type names
        static const NamedEnum<heatSourceType, 3> heatSourceTypeNames_;
        
        static const NamedEnum<heatFluxProfileType, 2> heatFluxProfileTypeNames_;
        
        static const NamedEnum<unitsType,7> unitsName_;

        //- Number of sources
        label nSources_;

        //- Heat source type
        List<heatSourceType> heatSource_;

        //- Source dictionaries
        List<dictionary> sourceDicts_;

        //- Name of source dict
        List<word> sourceDictName_;
        
        //- Heat Flux Profile i.e. gaussian or flat
        List<heatFluxProfileType> heatFluxProfile_;

        //- Heat power [W] or flux [W/m2]
        scalarField q_;
        
        //- Absorbed Beam Power [W]
        List<scalar> Q_;

        //- Beam Radius for Circular Beam
        List<scalar> radius_;

		//- Semi Major Axis in X-direction for Ellipse Beam
        List<scalar> semiMajorAxis_;

		//- Semi Minor Axis in Z-direction for Ellipse Beam
        List<scalar> semiMinorAxis_;

        //- Length of Strip Beam in X-direction
        List<scalar> length_;
		
		//- Breadth of Strip Beam in Z-direction
        List<scalar> breadth_;
		
		//- Centre of Beam
        List<vector> centre_;

        //- Name of radiative in flux field
        word QrName_;
        
        //- Amplitude of Oscillation
        List<vector> amplitude_;
        
        //- Frequency of Oscillation
        List<scalar> frequency_;

        //- Emissivity
        scalar emissivity_;
        
        //- Ambient Temperature for radiation exchange
        scalar Tamb_;

        //- Ask if Evaporation is ON
        Switch evaporation_;

        //- A term in log10(p) = -A+B/T
        scalar  A_;

        //- B term in log10(p) = -A+B/T
        scalar  B_;

        //- Mass Number for evaporated material
        scalar MassNumber_;

        //- Latent Heat of Vapourisation for evaporated material
        scalar  vapourisationEnthalpy_;

        //- Unit of Evaporation equation i.e. Pa, mmHg, atm, Pa etc.
        unitsType units;


public:

    //- Runtime type information
    TypeName("compressible::ebeamTurbulentHeatFlux")


    // Constructors

        //- Construct from patch and internal field
        ebeamTurbulentHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ebeamTurbulentHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  ebeamTurbulentHeatFluxFvPatchScalarField onto
        //  a new patch
        ebeamTurbulentHeatFluxFvPatchScalarField
        (
            const ebeamTurbulentHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ebeamTurbulentHeatFluxFvPatchScalarField
        (
            const ebeamTurbulentHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new ebeamTurbulentHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ebeamTurbulentHeatFluxFvPatchScalarField
        (
            const ebeamTurbulentHeatFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new ebeamTurbulentHeatFluxFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
            
            //- Check if face within beam boundary
            virtual void checkBoundary();
            
            //- Calculate evaporation loss
            virtual scalarField evaporationLoss();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
