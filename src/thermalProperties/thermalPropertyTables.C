/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.
\*---------------------------------------------------------------------------*/

#include "thermalPropertyTables.H"
#include "fvMesh.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::thermalPropertyTables::thermalPropertyTables()
:
    DFluid_(Zero),
    DSolid_(Zero),
    rhoCpRef_(Zero),
    TRef_(Zero),
    rhoPtr_(nullptr),
    cpPtr_(nullptr),
    kFluidPtr_(nullptr),
    DSolidPtr_(nullptr)
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::uniform
(
    const volScalarField& T,
    const word& name,
    const scalar value,
    const dimensionSet& dims
) const
{
    const fvMesh& mesh = T.mesh();

    return tmp<volScalarField>::New
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(name, dims, value)
    );
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::evaluate
(
    const Function1<scalar>& f,
    const volScalarField& T,
    const word& name,
    const scalar scale,
    const dimensionSet& dims
) const
{
    tmp<volScalarField> tfld(uniform(T, name, Zero, dims));
    volScalarField& fld = tfld.ref();

    scalarField& fi = fld.primitiveFieldRef();
    const scalarField& Ti = T.primitiveField();
    forAll(fi, celli)
    {
        fi[celli] = scale*f.value(Ti[celli]);
    }

    // Boundary faces are evaluated at the boundary temperature, not left at
    // the internal or a constant value: on a patch where solid meets the
    // domain boundary (a heated wall, say) the wrong face diffusivity would
    // otherwise be used in the Laplacian flux.
    volScalarField::Boundary& bfld = fld.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchi];
        fvPatchScalarField& fp = bfld[patchi];

        forAll(fp, facei)
        {
            fp[facei] = scale*f.value(Tp[facei]);
        }
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::thermalPropertyTables::read
(
    const dictionary& thermalDict,
    const fvMesh& mesh
)
{
    rhoPtr_.reset(nullptr);
    cpPtr_.reset(nullptr);
    kFluidPtr_.reset(nullptr);
    DSolidPtr_.reset(nullptr);

    if (thermalDict.found("rhoTable") != thermalDict.found("cpTable"))
    {
        FatalIOErrorInFunction(thermalDict)
            << "rhoTable and cpTable must be given together: the equation is "
            << "scaled by the volumetric heat capacity rho*cp, so one without "
            << "the other is not meaningful."
            << exit(FatalIOError);
    }

    if (thermalDict.found("rhoTable"))
    {
        rhoPtr_ = Function1<scalar>::New("rhoTable", thermalDict, &mesh);
        cpPtr_ = Function1<scalar>::New("cpTable", thermalDict, &mesh);
    }
    if (thermalDict.found("kFluidTable"))
    {
        kFluidPtr_ = Function1<scalar>::New("kFluidTable", thermalDict, &mesh);
    }
    if (thermalDict.found("DSolidTable"))
    {
        DSolidPtr_ = Function1<scalar>::New("DSolidTable", thermalDict, &mesh);
    }

    // Constants: required only where the corresponding table is absent
    DFluid_ =
        kFluidPtr_
      ? thermalDict.getOrDefault<scalar>("DFluid", Zero)
      : thermalDict.get<scalar>("DFluid");

    DSolid_ =
        DSolidPtr_
      ? thermalDict.getOrDefault<scalar>("DSolid", Zero)
      : thermalDict.get<scalar>("DSolid");

    // Reference volumetric heat capacity
    if (thermalDict.found("rhoCpRef"))
    {
        rhoCpRef_ = thermalDict.get<scalar>("rhoCpRef");
        TRef_ = thermalDict.getOrDefault<scalar>("TRef", Zero);
    }
    else if (rhoPtr_)
    {
        TRef_ = thermalDict.get<scalar>("TRef");
        rhoCpRef_ = rhoPtr_->value(TRef_)*cpPtr_->value(TRef_);
    }
    else if (kFluidPtr_)
    {
        FatalIOErrorInFunction(thermalDict)
            << "kFluidTable gives k [W/m/K], which is converted to a "
            << "diffusivity by (rho cp)_ref. Supply either rhoCpRef, or "
            << "rhoTable + cpTable + TRef."
            << exit(FatalIOError);
    }

    if (rhoCpRef_ < SMALL && (rhoPtr_ || kFluidPtr_))
    {
        FatalIOErrorInFunction(thermalDict)
            << "Non-positive reference volumetric heat capacity "
            << rhoCpRef_ << exit(FatalIOError);
    }

    if (rhoPtr_ || kFluidPtr_ || DSolidPtr_)
    {
        Info<< "thermalTopO: temperature-dependent properties active"
            << " [fluid k: " << Switch(bool(kFluidPtr_))
            << ", fluid rho*cp: " << Switch(bool(rhoPtr_))
            << ", solid D: " << Switch(bool(DSolidPtr_)) << "]" << nl;

        if (rhoPtr_ || kFluidPtr_)
        {
            Info<< "    (rho cp)_ref = " << rhoCpRef_ << " J/m3/K";
            if (rhoPtr_)
            {
                Info<< " at TRef = " << TRef_ << " K";
            }
            Info<< endl;
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::DFluidField
(
    const volScalarField& T,
    const word& name
) const
{
    if (kFluidPtr_)
    {
        return evaluate(*kFluidPtr_, T, name, scalar(1)/rhoCpRef_, dimViscosity);
    }

    return uniform(T, name, DFluid_, dimViscosity);
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::DSolidField
(
    const volScalarField& T,
    const word& name
) const
{
    if (DSolidPtr_)
    {
        return evaluate(*DSolidPtr_, T, name, scalar(1), dimViscosity);
    }

    return uniform(T, name, DSolid_, dimViscosity);
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::CField
(
    const volScalarField& T,
    const word& name
) const
{
    if (!rhoPtr_)
    {
        // Constant properties: C == 1, and the formulation collapses to the
        // constant-property one exactly.
        return uniform(T, name, scalar(1), dimless);
    }

    tmp<volScalarField> tC(uniform(T, name, scalar(1), dimless));
    volScalarField& C = tC.ref();

    const scalar rrhoCpRef = scalar(1)/rhoCpRef_;

    scalarField& Ci = C.primitiveFieldRef();
    const scalarField& Ti = T.primitiveField();
    forAll(Ci, celli)
    {
        const scalar Tc = Ti[celli];
        Ci[celli] = rrhoCpRef*rhoPtr_->value(Tc)*cpPtr_->value(Tc);
    }

    volScalarField::Boundary& bC = C.boundaryFieldRef();
    forAll(bC, patchi)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchi];
        fvPatchScalarField& Cp = bC[patchi];

        forAll(Cp, facei)
        {
            const scalar Tf = Tp[facei];
            Cp[facei] = rrhoCpRef*rhoPtr_->value(Tf)*cpPtr_->value(Tf);
        }
    }

    return tC;
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::alphaField
(
    const volScalarField& T,
    const word& name
) const
{
    tmp<volScalarField> tAlpha(DFluidField(T, name));

    if (rhoPtr_)
    {
        tAlpha.ref() /= CField(T, "alphaRhoCpNorm")();
    }

    return tAlpha;
}


// ************************************************************************* //
