/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.
\*---------------------------------------------------------------------------*/

#include "temperatureTable.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(temperatureTable, 0);
    addToRunTimeSelectionTable
    (
        viscosityModel,
        temperatureTable,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::scalar Foam::viscosityModels::temperatureTable::nuValue
(
    const scalar T
) const
{
    if (nuPtr_)
    {
        return nuPtr_->value(T);
    }

    return muPtr_->value(T)/rhoPtr_->value(T);
}


void Foam::viscosityModels::temperatureTable::updateNu()
{
    const auto* TPtr = mesh_.findObject<volScalarField>(TName_);

    if (!TPtr)
    {
        // Temperature not registered yet (the transport model is constructed
        // before the primal solver allocates T). Fall back to the reference.
        nu_ = dimensionedScalar(nu_.name(), dimViscosity, nuValue(TRef_));
        return;
    }

    const volScalarField& T = *TPtr;

    scalarField& nui = nu_.primitiveFieldRef();
    const scalarField& Ti = T.primitiveField();
    forAll(nui, celli)
    {
        nui[celli] = nuValue(Ti[celli]);
    }

    volScalarField::Boundary& bnu = nu_.boundaryFieldRef();
    forAll(bnu, patchi)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchi];
        fvPatchScalarField& nup = bnu[patchi];

        forAll(nup, facei)
        {
            nup[facei] = nuValue(Tp[facei]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::viscosityModels::temperatureTable::temperatureTable
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    coeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    muPtr_(nullptr),
    rhoPtr_(nullptr),
    nuPtr_(nullptr),
    TName_(coeffs_.getOrDefault<word>("field", "T")),
    TRef_(coeffs_.get<scalar>("TRef")),
    mesh_(U.mesh()),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar(name, dimViscosity, Zero)
    )
{
    read(viscosityProperties);

    nu_ = dimensionedScalar(name, dimViscosity, nuValue(TRef_));

    // Picks up T if it already exists; otherwise the first correct() will
    updateNu();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

bool Foam::viscosityModels::temperatureTable::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    coeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    TName_ = coeffs_.getOrDefault<word>("field", "T");
    TRef_ = coeffs_.get<scalar>("TRef");

    muPtr_.reset(nullptr);
    rhoPtr_.reset(nullptr);
    nuPtr_.reset(nullptr);

    if (coeffs_.found("nuTable"))
    {
        nuPtr_ = Function1<scalar>::New("nuTable", coeffs_, &mesh_);
    }
    else
    {
        if (!coeffs_.found("muTable") || !coeffs_.found("rhoTable"))
        {
            FatalIOErrorInFunction(coeffs_)
                << "temperatureTable needs either nuTable [m2/s], or both "
                << "muTable [Pa s] and rhoTable [kg/m3]."
                << exit(FatalIOError);
        }

        muPtr_ = Function1<scalar>::New("muTable", coeffs_, &mesh_);
        rhoPtr_ = Function1<scalar>::New("rhoTable", coeffs_, &mesh_);
    }

    return true;
}


// ************************************************************************* //
