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

Class
    Foam::marangoniGradientConditionFvPatchVectorField

Group
    grpTurbulenceBoundaryConditions

Description
    Set a marangoni condition as nuEff dU/dn = (dSigDT/rho) dT/ds.

Usage
    \table
        Property	| Description		| Required	| Default Value
        dSigDT		| temperature coefficient of Surface Tension | yes | 0.0
        rho			| density of the fluid (constant)	| yes | 1.0
    \endtable

Example usage:
    \verbatim
    myPatch
    {
        type    marangoniCondition;
        dSigDT  -1.0e-4;
        rho     1000;
    }
    \endverbatim

See also
    Foam::fixedShearStressFvPatchVectorField

SourceFiles
    marangoniGradientConditionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef marangoniGradientConditionFvPatchVectorField_H
#define marangoniGradientConditionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class marangoniGradientConditionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class marangoniGradientConditionFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Constant temperature coefficient of Surface Tension
        const scalar dSigDT_;
        
        //- Constant Density of fluid
        const scalar rho_;


public:

    //- Runtime type information
    TypeName("marangoniCondition");


    // Constructors

        //- Construct from patch and internal field
        marangoniGradientConditionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        marangoniGradientConditionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        marangoniGradientConditionFvPatchVectorField
        (
            const marangoniGradientConditionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        marangoniGradientConditionFvPatchVectorField
        (
            const marangoniGradientConditionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new marangoniGradientConditionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        marangoniGradientConditionFvPatchVectorField
        (
            const marangoniGradientConditionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new marangoniGradientConditionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
    
        //- Return gradient at boundary
        virtual tmp<vectorField> snGrad() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
        
        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::blocking
        );
        
        //- Return face-gradient transform diagonal
        virtual tmp<vectorField> snGradTransformDiag() const;
            
        //- Write
        virtual void write(Ostream&) const;
        
        //- Calculating the refGrad and refValue
        virtual void getValues();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
