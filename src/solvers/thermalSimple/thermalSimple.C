/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.
\*---------------------------------------------------------------------------*/

#include "thermalSimple.H"
#include "topOVariablesBase.H"
#include "fvOptions.H"
#include "fvm.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalSimple, 0);
    addToRunTimeSelectionTable
    (
        incompressiblePrimalSolver,
        thermalSimple,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::thermalSimple::allocateThermalFields()
{
    TPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "T",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::thermalSimple::kIndicator() const
{
    auto tindicator =
        tmp<volScalarField>::New
        (
            IOobject
            (
                "kIndicator",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
    volScalarField& indicator = tindicator.ref();

    if (mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        const topOVariablesBase& vars =
            mesh_.lookupObject<topOVariablesBase>("topOVars");

        DimensionedField<scalar, volMesh> interpolant
        (
            IOobject
            (
                "kInterpolant",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );

        // unit multiplier: interpolant becomes Ik(beta) in [0,1]
        vars.sourceTerm(interpolant, kInterpolation_(), scalar(1), "beta");

        indicator.primitiveFieldRef() = interpolant.field();
        indicator.correctBoundaryConditions();
    }

    return tindicator;
}


Foam::tmp<Foam::volScalarField> Foam::thermalSimple::DEff() const
{
    const volScalarField Ik(kIndicator());

    const autoPtr<incompressible::turbulenceModel>& turbulence =
        incoVars_.turbulence();

    auto tDEff =
        tmp<volScalarField>::New
        (
            IOobject
            (
                "DEff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("DFluid", dimViscosity, DFluid_)
          + dimensionedScalar("dD", dimViscosity, DSolid_ - DFluid_)*Ik
          + (scalar(1) - Ik)*turbulence->nut()/Prt_
        );

    return tDEff;
}


void Foam::thermalSimple::solveTEqn()
{
    const surfaceScalarField& phi = incoVars_.phiInst();
    volScalarField& T = TPtr_.ref();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    const volScalarField DEff(this->DEff());

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian(DEff, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();
    fvOptions.constrain(TEqn);
    TEqn.solve().initialResidual();
    fvOptions.correct(T);

    Info<< "Min/max T: " << gMin(T.primitiveField()) << ", "
        << gMax(T.primitiveField()) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::thermalSimple::thermalSimple
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& solverName
)
:
    simple(mesh, managerType, dict, solverName),
    TPtr_(nullptr),
    DFluid_(Zero),
    DSolid_(Zero),
    Prt_(0.85),
    kInterpolation_(nullptr)
{
    const dictionary& thermalDict = dict.subDict("thermal");

    DFluid_ = thermalDict.get<scalar>("DFluid");
    DSolid_ = thermalDict.get<scalar>("DSolid");
    Prt_ = thermalDict.getOrDefault<scalar>("Prt", 0.85);
    kInterpolation_ =
        topOInterpolationFunction::New
        (
            mesh_,
            thermalDict.subDict("kInterpolation")
        );

    allocateThermalFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

bool Foam::thermalSimple::readDict(const dictionary& dict)
{
    if (simple::readDict(dict))
    {
        const dictionary& thermalDict = dict.subDict("thermal");
        DFluid_ = thermalDict.get<scalar>("DFluid");
        DSolid_ = thermalDict.get<scalar>("DSolid");
        Prt_ = thermalDict.getOrDefault<scalar>("Prt", 0.85);

        return true;
    }

    return false;
}


void Foam::thermalSimple::mainIter()
{
    simple::mainIter();
    solveTEqn();
}


// ************************************************************************* //
