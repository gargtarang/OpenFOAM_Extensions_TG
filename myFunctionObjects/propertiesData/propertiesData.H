/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::functionObjects::propertiesData

Group
    grpFieldFunctionObjects

Description
    This function object stores and writes thermophysical data.
    This functionObject is functional but some code cleanup is 
    required. Use with care.

    Fields are stored as copies of the original, with the prefix
    <functionObjectName>_, e.g.:

    \verbatim
    propertiesData1_Cp
    \endverbatim

    Example of function object specification:
    \verbatim
    propertiesData1
    {
        type        propertiesData;
        libs ("libfieldFunctionObjects.so");
        ...
        fields
        (
            Cp
            Mu
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description                 | Required | Default value
        type         | type name: processorField   | yes      |
        fields       | fields to store (see below) | yes      |
    \endtable

    Where \c fields can include:
    \plaintable
        Cp          | Constant pressure heat capacity
        kappa       | Thermal Conductivity
        mu          | Dynamic Viscosity
        alpha       | Thermal Diffusivity
    \endplaintable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    propertiesData.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_propertiesData_H
#define functionObjects_propertiesData_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"
#include "NamedEnum.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class propertiesData Declaration
\*---------------------------------------------------------------------------*/

class propertiesData
:
    public fvMeshFunctionObject
{
public:

    enum compressibleField
    {
        cfMu,
        cfKappa,
        cfalpha,
        cfCp
    };
    static const NamedEnum<compressibleField, 4> compressibleFieldNames_;
    static const word modelName;


protected:

    // Protected data

        //- Fields to load
        wordHashSet fieldSet_;


    // Protected Member Functions

        //- Return true if compressible turbulence model is identified
        bool compressible();

        //- Process the turbulence field
        template<class Type>
        void processField
        (
            const word& fieldName,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
        );

        //- Return omega calculated from k and epsilon
        template<class Model>
        tmp<volScalarField> omega(const Model& model) const;


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        propertiesData(const propertiesData&);

        //- Disallow default bitwise assignment
        void operator=(const propertiesData&);


public:

    //- Runtime type information
    TypeName("propertiesData");


    // Constructors

        //- Construct from Time and dictionary
        propertiesData
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~propertiesData();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate turbulence fields
        virtual bool execute();

        //- Do nothing.
        //  The turbulence fields are registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "propertiesDataTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
