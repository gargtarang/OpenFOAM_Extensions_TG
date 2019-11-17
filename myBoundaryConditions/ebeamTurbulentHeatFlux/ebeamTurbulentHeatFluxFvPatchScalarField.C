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

\*---------------------------------------------------------------------------*/

#include "ebeamTurbulentHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        Foam::compressible::
        ebeamTurbulentHeatFluxFvPatchScalarField::heatSourceType,
        3
    >::names[] =
    {
        "circle",
        "ellipse",
        "strip"
    };
    
    template<>
    const char* NamedEnum
    <
        Foam::compressible::
        ebeamTurbulentHeatFluxFvPatchScalarField::heatFluxProfileType,
        2
    >::names[] =
    {
        "gaussian",
        "flat"
    };
    
    template<>
    const char* NamedEnum
    <
        Foam::compressible::
        ebeamTurbulentHeatFluxFvPatchScalarField::unitsType,
        7
    >::names[] =
    {
        "Pa",
        "atm",
        "mmHg",
        "cmHg",
        "bar",
        "kPa",
        "MPa"
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const NamedEnum
<
    ebeamTurbulentHeatFluxFvPatchScalarField::heatSourceType,
    3
> ebeamTurbulentHeatFluxFvPatchScalarField::heatSourceTypeNames_;

const NamedEnum
<
    ebeamTurbulentHeatFluxFvPatchScalarField::heatFluxProfileType,
    2
> ebeamTurbulentHeatFluxFvPatchScalarField::heatFluxProfileTypeNames_;

const NamedEnum
<
    ebeamTurbulentHeatFluxFvPatchScalarField::unitsType,
    7
> ebeamTurbulentHeatFluxFvPatchScalarField::unitsName_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ebeamTurbulentHeatFluxFvPatchScalarField::
ebeamTurbulentHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    nSources_(0),
//    heatSource_(hsStrip),
    sourceDicts_(0,NULL),
    sourceDictName_(0,"undefined"),
//    heatFluxProfile_(hfpFlat),
    q_(p.size(), 0.0),
//    Q_(0.0),
    QrName_("undefinedQr"),
    emissivity_(0.0),
    Tamb_(0.0),
    evaporation_(false),
    A_(0.0),
    B_(0.0),
    MassNumber_(1.0),
    vapourisationEnthalpy_(0.0),
    units(unPa)
{}


ebeamTurbulentHeatFluxFvPatchScalarField::
ebeamTurbulentHeatFluxFvPatchScalarField
(
    const ebeamTurbulentHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    nSources_(ptf.nSources_),
//    heatSource_(ptf.heatSource_),
    sourceDicts_(ptf.sourceDicts_),
    sourceDictName_(ptf.sourceDictName_),
//    heatFluxProfile_(ptf.heatFluxProfile_),
    q_(ptf.q_, mapper),
//    Q_(ptf.Q_),
//    radius_(ptf.radius_),
//    semiMajorAxis_(ptf.semiMajorAxis_),
//    semiMinorAxis_(ptf.semiMinorAxis_),
//    length_(ptf.length_),
//    breadth_(ptf.breadth_),
//    centre_(ptf.centre_),
    QrName_(ptf.QrName_),
//    amplitude_(ptf.amplitude_),
//    frequency_(ptf.frequency_),
    emissivity_(ptf.emissivity_),
    Tamb_(ptf.Tamb_),
    evaporation_(ptf.evaporation_),
    A_(ptf.A_),
    B_(ptf.B_),
    MassNumber_(ptf.MassNumber_),
    vapourisationEnthalpy_(ptf.vapourisationEnthalpy_),
    units(ptf.units)
{}


ebeamTurbulentHeatFluxFvPatchScalarField::
ebeamTurbulentHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    nSources_(dict.lookupOrDefault<label>("sources",0)),
    QrName_(dict.lookupOrDefault<word>("Qr", "none")),
    emissivity_(dict.lookupOrDefault<scalar>("emissivity",Zero)),
    Tamb_(dict.lookupOrDefault<scalar>("Tambient",Zero)),
    evaporation_(dict.lookupOrDefault<Switch>("evaporation",false))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        // Still reading so cannot yet evaluate. Make up a value.
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }

    if(nSources_ < 1)
    {
        FatalErrorInFunction
            << "Atleast 1 source is required"
            << nl << exit(FatalError);
    }

    // Here we are searching for any subdictionaries present.
    // dict.toc() gives names of all words present on the left hand side
    forAll(dict.toc(),c)
    {
        if(dict.isDict(dict.toc()[c]))      //checking if any word is a sub-dict.
        {
            //- Subdict found -> new Ebeam source
            sourceDictName_.append(dict.toc()[c]);
        }
    }
    //- We have now names of all Ebeam source

    for (label s = 0; s < nSources_; s++)
    {
        sourceDicts_.append(dict.subDict(sourceDictName_[s]));
        const dictionary& sourceDict = sourceDicts_[s];
        heatSource_.append(heatSourceTypeNames_.read(sourceDict.lookup("heatSource")));
        heatFluxProfile_.append(heatFluxProfileTypeNames_.read(sourceDict.lookup("heatProfile")));
        Q_.append(sourceDict.lookupOrDefault<scalar>("Q",Zero));
        centre_.append(sourceDict.lookupOrDefault<vector>("centre",Zero));
        amplitude_.append(sourceDict.lookupOrDefault<vector>("amplitude",Zero));
        frequency_.append(sourceDict.lookupOrDefault<scalar>("frequency",Zero));
        switch(heatSource_[s])
        {
            case hsCircle:
            {
                radius_.append(sourceDict.lookupOrDefault<scalar>("radius",1e+15));
                break;
            }
            case hsEllipse:
            {
                semiMajorAxis_.append(sourceDict.lookupOrDefault<scalar>("semiMajorAxis",1e+15));
                semiMinorAxis_.append(sourceDict.lookupOrDefault<scalar>("semiMinorAxis",1e+15));
                break;
            }
            case hsStrip:
            {
                length_.append(sourceDict.lookupOrDefault<scalar>("length",1e+15));
                breadth_.append(sourceDict.lookupOrDefault<scalar>("breadth",1e+15));
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown heat source type. Valid types are: "
                    << heatSourceTypeNames_ << nl << exit(FatalError);
                break;
            }
        }
    }
    if(evaporation_)
    {
        dict.lookup("A") >> A_;
        dict.lookup("B") >> B_;
        dict.lookup("massNumber") >> MassNumber_;
        dict.lookup("vapourisationEnthalpy") >> vapourisationEnthalpy_;
        units = unitsName_.read(dict.lookup("units"));
        volScalarField* evaporationMassPtr
        (
             new volScalarField
             (
                IOobject
                (
                    "evaporationMassFlux",
                    this->internalField().mesh().time().timeName(),
                    this->internalField().mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->internalField().mesh(),
                dimensionedScalar("0",dimMass/(dimArea*dimTime),0)
              )
        );
        this->internalField().mesh().objectRegistry::store(evaporationMassPtr);
    }
    q_.setSize(p.size());
    checkBoundary();
}


ebeamTurbulentHeatFluxFvPatchScalarField::
ebeamTurbulentHeatFluxFvPatchScalarField
(
    const ebeamTurbulentHeatFluxFvPatchScalarField& ethfpsf
)
:
    fixedGradientFvPatchScalarField(ethfpsf),
    temperatureCoupledBase(patch(), ethfpsf),
    nSources_(ethfpsf.nSources_),
//    heatSource_(ethfpsf.heatSource_),
//    sourceDicts_(ethfpsf.sourceDicts_),
//    sourceDictName_(ethfpsf.sourceDictName_),
//    heatFluxProfile_(ethfpsf.heatFluxProfile_),
    q_(ethfpsf.q_),
//    Q_(ethfpsf.Q_),
//    radius_(ethfpsf.radius_),
//    semiMajorAxis_(ethfpsf.semiMajorAxis_),
//    semiMinorAxis_(ethfpsf.semiMinorAxis_),
//    length_(ethfpsf.length_),
//    breadth_(ethfpsf.breadth_),
//    centre_(ethfpsf.centre_),
    QrName_(ethfpsf.QrName_),
//    amplitude_(ethfpsf.amplitude_),
//    frequency_(ethfpsf.frequency_),
    emissivity_(ethfpsf.emissivity_),
    Tamb_(ethfpsf.Tamb_),
    evaporation_(ethfpsf.evaporation_),
    A_(ethfpsf.A_),
    B_(ethfpsf.B_),
    MassNumber_(ethfpsf.MassNumber_),
    vapourisationEnthalpy_(ethfpsf.vapourisationEnthalpy_),
    units(ethfpsf.units)
{}


ebeamTurbulentHeatFluxFvPatchScalarField::
ebeamTurbulentHeatFluxFvPatchScalarField
(
    const ebeamTurbulentHeatFluxFvPatchScalarField& ethfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ethfpsf, iF),
    temperatureCoupledBase(patch(), ethfpsf),
    nSources_(ethfpsf.nSources_),
    //heatSource_(ethfpsf.heatSource_),
    sourceDicts_(ethfpsf.sourceDicts_),
    sourceDictName_(ethfpsf.sourceDictName_),
    //heatFluxProfile_(ethfpsf.heatFluxProfile_),
    q_(ethfpsf.q_),
    //Q_(ethfpsf.Q_),
    //radius_(ethfpsf.radius_),
//    semiMajorAxis_(ethfpsf.semiMajorAxis_),
//    semiMinorAxis_(ethfpsf.semiMinorAxis_),
//    length_(ethfpsf.length_),
//    breadth_(ethfpsf.breadth_),
//    centre_(ethfpsf.centre_),
    QrName_(ethfpsf.QrName_),
//    amplitude_(ethfpsf.amplitude_),
//    frequency_(ethfpsf.frequency_),
    emissivity_(ethfpsf.emissivity_),
    Tamb_(ethfpsf.Tamb_),
    evaporation_(ethfpsf.evaporation_),
    A_(ethfpsf.A_),
    B_(ethfpsf.B_),
    MassNumber_(ethfpsf.MassNumber_),
    vapourisationEnthalpy_(ethfpsf.vapourisationEnthalpy_),
    units(ethfpsf.units)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ebeamTurbulentHeatFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    q_.autoMap(m);
}


void ebeamTurbulentHeatFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const ebeamTurbulentHeatFluxFvPatchScalarField& ethfpsf =
        refCast<const ebeamTurbulentHeatFluxFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(ethfpsf.q_, addr);
}


void ebeamTurbulentHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tp = *this;

    scalarField qr(this->size(), 0.0);

    //- Qr is negative going into the domain
    if (QrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField, scalar>(QrName_);
    }
    
    checkBoundary();
    scalarField Tamb(this->size(),Tamb_);
    qr = qr - emissivity_ * 5.6704e-8 * (pow4(Tp) - pow4(Tamb));
    gradient() = (q_ + qr - (((evaporation_)?evaporationLoss():scalarField(this->size(),0.0))))/kappa(Tp);
    fixedGradientFvPatchScalarField::updateCoeffs();
}

scalarField ebeamTurbulentHeatFluxFvPatchScalarField::evaporationLoss()
{
    const scalarField& Ts = *this;
    scalarField evapLoss(this->size(),0.0);
    scalarField pressure(this->size(),0.0);
    scalarField temp_(-A_/(Ts) + B_);
    pressure = pow(10,temp_);
    switch(units)
    {
        case unPa:
            pressure = pow(10,temp_); break;
        case unatm:
            pressure = pow(10,temp_) * 101325; break;
        case unmmHg:
            pressure = pow(10,temp_) * 101325/760; break;
        case uncmHg:
            pressure = pow(10,temp_) * 101325/76; break;
        case unbar:
            pressure = pow(10,temp_) * 100000; break;
        case unkPa:
            pressure = pow(10,temp_) * 1000; break;
        case unMPa:
            pressure = pow(10,temp_) * 1000000; break;
        default:
        {
            FatalErrorInFunction
                << "Unknown units name. Valid units name are: "
                << unitsName_ << nl << exit(FatalError); break;
        }
    }
    
    scalarField massFlux(pressure * sqrt(8*0.001*MassNumber_/(M_PI*8.314*Ts)));
//    scalarField evaporationMass = this->patch().magSf()*massFlux;
    volScalarField& evaporationMassLoss =
                this->internalField().mesh().lookupObjectRef<volScalarField>("evaporationMassFlux");
    evaporationMassLoss.boundaryFieldRef()[this->patch().index()] = massFlux;
    evapLoss = massFlux * (2*8.314*Ts/(0.001*MassNumber_) + vapourisationEnthalpy_);
    //Info << "mass Flux" << massFlux <<nl;
    return evapLoss;
}

void ebeamTurbulentHeatFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("sources") << nSources_<<token::END_STATEMENT<<nl;
    forAll(sourceDicts_,d)
        os.writeKeyword(sourceDictName_[d]) << sourceDicts_[d];

    if(emissivity_!=0.0)
    {
        os.writeKeyword("emissivity") << emissivity_ 
                << token::END_STATEMENT << nl;
        os.writeKeyword("Tambient") << Tamb_
                << token::END_STATEMENT << nl;
    }

    if(evaporation_)
    {
        os.writeKeyword("evaporation") << evaporation_
                << token::END_STATEMENT << nl;
        os.writeKeyword("A") << A_
                << token::END_STATEMENT << nl;
        os.writeKeyword("B") << B_
                << token::END_STATEMENT << nl;
        os.writeKeyword("massNumber") << MassNumber_
                << token::END_STATEMENT << nl;
        os.writeKeyword("vapourisationEnthalpy") << vapourisationEnthalpy_
                << token::END_STATEMENT << nl;
        os.writeKeyword("units") << unitsName_[units]
                << token::END_STATEMENT << nl;
    }
    os.writeKeyword("Qr")<< QrName_ << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
    writeEntry("value", os);
}

void ebeamTurbulentHeatFluxFvPatchScalarField::checkBoundary()
{
//    Info<< "radius= "<< radius_<<nl
//        << "length= "<<length_<<nl
//        << "breadth= "<<breadth_<<nl;
    scalar area_beam = 0.0;
    scalar flux_ = 0.0;
    vector temp = Zero;
    q_=0.0;
    label cCtr=0;
    label eCtr=0;
    label sCtr=0;
    const scalar& time = this->db().time().value();
    //- Patch
    const fvPatch& p = this->patch();

    //- Face Centres
    const List<point>& cf = p.Cf();
    forAll(sourceDicts_,d)
    {
        temp[0] = centre_[d][0]+amplitude_[d][0]*sin(2*M_PI*frequency_[d]*time);
        temp[1] = centre_[d][1]+amplitude_[d][1]*sin(2*M_PI*frequency_[d]*time);
        temp[2] = centre_[d][2]+amplitude_[d][2]*sin(2*M_PI*frequency_[d]*time);

        switch (heatSource_[d])
        {
            case hsCircle:
            {
                area_beam = M_PI*radius_[cCtr]*radius_[cCtr];
                flux_ = Q_[d]/area_beam;
                forAll(cf, c)
                {
                    const scalar xCf = cf[c][0];
                    const scalar zCf = cf[c][2];
                    scalar residue =
                        pow(((temp[0]-xCf)/radius_[cCtr]),2)+pow(((temp[2]-zCf)/radius_[cCtr]),2);
                    switch (heatFluxProfile_[d])
                    {
                        case hfpGaussian:
                        {
                            q_[c] = 0.5*flux_*exp(-0.5*residue);
                            break;
                        }

                        case hfpFlat:
                        {
                            if(residue<=1.0)
                            {
                                q_[c] = flux_;
                            }
                            break;
                        }
                    
                        default:
                        {
                            FatalErrorInFunction
                                << "Unknown heat flux profile type. Valid types are: "
                                << heatFluxProfileTypeNames_ << nl << exit(FatalError);
                        }
                    }
                }
                cCtr++;
                break;
            }
            case hsEllipse:
            {
                area_beam = M_PI*semiMajorAxis_[eCtr]*semiMinorAxis_[eCtr];
                flux_ = Q_[d]/area_beam;
                forAll(cf, c)
                {
                    const scalar xCf = cf[c][0];
                    const scalar zCf = cf[c][2];
                    scalar residue =
                        pow(((xCf-temp[0])/semiMajorAxis_[eCtr]),2) +
                        pow(((zCf-temp[2])/semiMinorAxis_[eCtr]),2);
                    switch (heatFluxProfile_[d])
                    {
                        case hfpGaussian:
                        {
                            q_[c] = 0.5*flux_*exp(-1*residue/2);
                            break;
                        }
                    
                        case hfpFlat:
                        {
                            if((residue/radius_[cCtr])<=1.0)
                            {
                                q_[c] = flux_;
                            }
                            break;
                        }
                    
                        default:
                        {
                            FatalErrorInFunction
                                << "Unknown heat flux profile type. Valid types are: "
                                << heatFluxProfileTypeNames_ << nl << exit(FatalError);
                        }
                    }
                }
                eCtr++;
                break;
            }
            case hsStrip:
            {
                area_beam = length_[sCtr]*breadth_[sCtr];
                flux_ = Q_[d]/area_beam;
                forAll(cf, c)
                {
                    const scalar xCf = cf[c][0];
                    const scalar zCf = cf[c][2];
                    scalar residue1 = 2*fabs(xCf-temp[0])/length_[sCtr];
                    scalar residue2 = 2*fabs(zCf-temp[2])/breadth_[sCtr];
                    if(max(residue1,residue2)<=1.0)
                    {
                        q_[c] = flux_;
                    }
                }
                sCtr++;
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown heat source type. Valid types are: "
                    << heatSourceTypeNames_ << nl << exit(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    ebeamTurbulentHeatFluxFvPatchScalarField
)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
