/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "myWallHeatFlux.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "emptyPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(myWallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, myWallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::myWallHeatFlux::writeFileHeader(const label i)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    
    // Add headers to output data
    writeHeader(file(), "Wall heat-flux");
    writeCommented(file(), "Time");
    forAll(pbm, patchi)
    {
        if(isA<wallPolyPatch>(pbm[patchi]))
            writeTabbed(file(),pbm[patchi].name());
    }

    writeTabbed(file(), "heatBalance");
    file() << endl;
}


void Foam::functionObjects::myWallHeatFlux::calcHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he,
    volScalarField& myWallHeatFlux
)
{
    surfaceScalarField heatFlux
    (
        fvc::interpolate(alpha)*fvc::snGrad(he)
    );

    volScalarField::Boundary& myWallHeatFluxBf =
        myWallHeatFlux.boundaryFieldRef();

    const surfaceScalarField::Boundary& heatFluxBf =
        heatFlux.boundaryField();

    forAll(myWallHeatFluxBf, patchi)
    {
        myWallHeatFluxBf[patchi] = heatFluxBf[patchi];
    }

    if (foundObject<volScalarField>("qr"))
    {
        const volScalarField& qr = lookupObject<volScalarField>("qr");

        const volScalarField::Boundary& radHeatFluxBf =
            qr.boundaryField();

        forAll(myWallHeatFluxBf, patchi)
        {
            myWallHeatFluxBf[patchi] += radHeatFluxBf[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::myWallHeatFlux::myWallHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    volScalarField* wallHeatFluxPtr
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimMass/pow3(dimTime), 0)
        )
    );

    mesh_.objectRegistry::store(wallHeatFluxPtr);

    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::myWallHeatFlux::~myWallHeatFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::myWallHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    writeObject_ = dict.lookupOrDefault("writeObject",false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::myWallHeatFlux::execute()
{
    volScalarField& myWallHeatFlux = lookupObjectRef<volScalarField>(type());

    if
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        calcHeatFlux
        (
            turbModel.alphaEff(),
            turbModel.transport().he(),
            myWallHeatFlux
        );
    }
    else if (foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            lookupObject<solidThermo>(solidThermo::dictName);

        calcHeatFlux(thermo.alpha(), thermo.he(), myWallHeatFlux);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::myWallHeatFlux::write()
{
    Log << type() << " " << name() << " write:" << nl
        <<"    writing field " << name() << endl;

    if(writeObject_)
    {
        writeLocalObjects::write();
    }

    logFiles::write();

    const volScalarField& myWallHeatFlux =
        obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();
    
    scalar sumHeat = 0;
    if (Pstream::master())
    {
        file() << mesh_.time().value() << token::TAB;
    }
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = myWallHeatFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);
        
        sumHeat += integralHfp;
        if (Pstream::master())
        {
            file() << integralHfp << token::TAB;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
    }
    if (Pstream::master())
    {
        file() << sumHeat << endl;
    }
    Log << "    net heat balance = " << sumHeat << endl;
    Log << endl;

    return true;
}


// ************************************************************************* //
