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

\*---------------------------------------------------------------------------*/

#include "propertiesData.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(propertiesData, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        propertiesData,
        dictionary
    );
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::propertiesData::compressibleField,
    4
>::names[] =
{
    "mu",
    "kappa",
    "alpha",
    "Cp"
};

const Foam::NamedEnum
<
    Foam::functionObjects::propertiesData::compressibleField,
    4
> Foam::functionObjects::propertiesData::compressibleFieldNames_;

const Foam::word Foam::functionObjects::propertiesData::modelName
(
    Foam::turbulenceModel::propertiesName
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::propertiesData::compressible()
{
    if (obr_.foundObject<basicThermo>(basicThermo::dictName))
    {
        return true;
    }
    else
    {
        FatalErrorInFunction
            << "Turbulence model not found in database, deactivating"
            << exit(FatalError);
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::propertiesData::propertiesData
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::propertiesData::~propertiesData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::propertiesData::read(const dictionary& dict)
{
    if (dict.found("field"))
    {
        fieldSet_.insert(word(dict.lookup("field")));
    }
    else
    {
        fieldSet_.insert(wordList(dict.lookup("fields")));
    }

    Info<< type() << " " << name() << ": ";
    if (fieldSet_.size())
    {
        Info<< "storing fields:" << nl;
        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            Info<< "    " << name() << '_' << iter.key() << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no fields requested to be stored" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::propertiesData::execute()
{
    bool comp = compressible();

    if (comp)
    {
        const fluidThermo& model =
            obr_.lookupObject<fluidThermo>("thermophysicalProperties");

        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            const word& f = iter.key();
            switch (compressibleFieldNames_[f])
            {
                case cfMu:
                {
                    processField<scalar>(f, model.mu());
                    break;
                }
                case cfKappa:
                {
                    processField<scalar>(f, model.kappa());
                    break;
                }
                case cfalpha:
                {
                    processField<scalar>(f, model.alpha());
                    break;
                }
                case cfCp:
                {
                    processField<scalar>(f, model.Cp());
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }
    return true;
}


bool Foam::functionObjects::propertiesData::write()
{
    forAllConstIter(wordHashSet, fieldSet_, iter)
    {
        const word fieldName = name() + '_' + iter.key();
        writeObject(fieldName);
    }

    return true;
}


// ************************************************************************* //
