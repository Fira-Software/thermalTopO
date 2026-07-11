/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "thermalAdjointSimple.H"
#include "topOVariablesBase.H"
#include "objectiveIncompressible.H"
#include "fixedGradientFvPatchFields.H"
#include "sensitivityTopO.H"
#include "fvm.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalAdjointSimple, 0);
    addToRunTimeSelectionTable
    (
        incompressibleAdjointSolver,
        thermalAdjointSimple,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::volScalarField& Foam::thermalAdjointSimple::TRef() const
{
    if (!mesh_.foundObject<volScalarField>("T"))
    {
        FatalErrorInFunction
            << "No registered temperature field 'T'. "
            << "thermalAdjointSimple requires the thermalSimple primal solver."
            << exit(FatalError);
    }

    return mesh_.lookupObject<volScalarField>("T");
}


Foam::tmp<Foam::volScalarField> Foam::thermalAdjointSimple::kIndicator() const
{
    auto tindicator =
        tmp<volScalarField>::New
        (
            IOobject
            (
                "kIndicatorAdj",
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
                "kInterpolantAdj",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );

        vars.sourceTerm(interpolant, kInterpolation_(), scalar(1), "beta");

        indicator.primitiveFieldRef() = interpolant.field();
        indicator.correctBoundaryConditions();
    }

    return tindicator;
}


Foam::tmp<Foam::volScalarField> Foam::thermalAdjointSimple::DEff() const
{
    const volScalarField Ik(kIndicator());

    const autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();

    auto tDs =
        tmp<volScalarField>::New
        (
            IOobject
            (
                "DSolidFieldAdj",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("DSolid", dimViscosity, DSolid_)
        );
    if (DSolidTablePtr_)
    {
        const volScalarField& T = TRef();
        scalarField& Ds = tDs.ref().primitiveFieldRef();
        forAll(Ds, celli)
        {
            Ds[celli] = DSolidTablePtr_->value(T[celli]);
        }
        tDs.ref().correctBoundaryConditions();
    }

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "DEffAdj",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("DFluid", dimViscosity, DFluid_)
      + (tDs() - dimensionedScalar(dimViscosity, DFluid_))*Ik
      + (scalar(1) - Ik)*turbulence->nut()/Prt_
    );
}


void Foam::thermalAdjointSimple::updateTaBCs(const volScalarField& DEffField)
{
    volScalarField& Ta = TaPtr_.ref();

    PtrList<objective>& functions =
        objectiveManager_.getObjectiveFunctions();

    forAll(mesh_.boundary(), pI)
    {
        fvPatchScalarField& Tap = Ta.boundaryFieldRef()[pI];

        if (isA<fixedGradientFvPatchScalarField>(Tap))
        {
            const fvPatch& patch = mesh_.boundary()[pI];
            scalarField flux(patch.size(), Zero);

            // Objective contributions: D dTa/dn = - sum w bdJdT
            for (objective& func : functions)
            {
                objectiveIncompressible& funcI =
                    refCast<objectiveIncompressible>(func);

                if (funcI.hasBoundarydJdT())
                {
                    flux -= func.weight()*funcI.boundarydJdT(pI);
                }
            }

            const scalarField DEffp
            (
                DEffField.boundaryField()[pI].patchInternalField()
            );

            refCast<fixedGradientFvPatchScalarField>(Tap).gradient() =
                flux/DEffp;
        }
    }
}


void Foam::thermalAdjointSimple::solveTaEqn()
{
    const surfaceScalarField& phi = primalVars_.phi();
    volScalarField& Ta = TaPtr_.ref();

    const volScalarField DEff(this->DEff());

    updateTaBCs(DEff);

    // Volume objective sources
    volScalarField source
    (
        IOobject
        (
            "TaSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(TaPtr_().dimensions()/dimTime, Zero)
    );

    PtrList<objective>& functions =
        objectiveManager_.getObjectiveFunctions();
    for (objective& func : functions)
    {
        objectiveIncompressible& funcI =
            refCast<objectiveIncompressible>(func);

        if (funcI.hasdJdT())
        {
            source.primitiveFieldRef() +=
                func.weight()*funcI.dJdT().primitiveField();
        }
    }

    fvScalarMatrix TaEqn
    (
        fvm::div(-phi, Ta)
      - fvm::laplacian(DEff, Ta)
    );

    Foam::solve(TaEqn == -source);

    Info<< "Min/max Ta: " << gMin(Ta.primitiveField()) << ", "
        << gMax(Ta.primitiveField()) << endl;
}


void Foam::thermalAdjointSimple::addMomentumSource(fvVectorMatrix& matrix)
{
    adjointSimple::addMomentumSource(matrix);

    const volScalarField& T = TRef();
    const volScalarField& Ta = TaPtr_();

    matrix += couplingSign_*(Ta*fvc::grad(T))();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::thermalAdjointSimple::thermalAdjointSimple
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName,
    const word& solverName
)
:
    adjointSimple(mesh, managerType, dict, primalSolverName, solverName),
    TaPtr_(nullptr),
    DFluid_(Zero),
    DSolid_(Zero),
    DSolidTablePtr_(nullptr),
    Prt_(0.85),
    couplingSign_(1),
    thermalSensScale_(1),
    kInterpolation_(nullptr)
{
    const dictionary& thermalDict = dict.subDict("thermal");

    DFluid_ = thermalDict.get<scalar>("DFluid");
    DSolid_ = thermalDict.get<scalar>("DSolid");
    Prt_ = thermalDict.getOrDefault<scalar>("Prt", 0.85);
    if (thermalDict.found("DSolidTable"))
    {
        DSolidTablePtr_ =
            Function1<scalar>::New("DSolidTable", thermalDict, &mesh_);
    }
    couplingSign_ = thermalDict.getOrDefault<scalar>("couplingSign", 1);
    thermalSensScale_ =
        thermalDict.getOrDefault<scalar>("thermalSensScale", 1);
    kInterpolation_ =
        topOInterpolationFunction::New
        (
            mesh_,
            thermalDict.subDict("kInterpolation")
        );

    TaPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Ta",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

bool Foam::thermalAdjointSimple::readDict(const dictionary& dict)
{
    if (adjointSimple::readDict(dict))
    {
        const dictionary& thermalDict = dict.subDict("thermal");
        DFluid_ = thermalDict.get<scalar>("DFluid");
        DSolid_ = thermalDict.get<scalar>("DSolid");
        Prt_ = thermalDict.getOrDefault<scalar>("Prt", 0.85);
        couplingSign_ =
            thermalDict.getOrDefault<scalar>("couplingSign", 1);
        thermalSensScale_ =
            thermalDict.getOrDefault<scalar>("thermalSensScale", 1);

        return true;
    }

    return false;
}


void Foam::thermalAdjointSimple::mainIter()
{
    solveTaEqn();
    adjointSimple::mainIter();
}


void Foam::thermalAdjointSimple::topOSensMultiplier
(
    scalarField& betaMult,
    const word& designVariablesName,
    const scalar dt
)
{
    // Brinkman momentum contribution (upstream machinery)
    adjointSimple::topOSensMultiplier(betaMult, designVariablesName, dt);

    // Conductivity-interpolation contribution:
    //   thermalSensScale * Ik'(beta) (DSolid - DFluid - nut/Prt)
    //   (grad(T) . grad(Ta)) dt
    if (mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        const volScalarField& T = TRef();
        const volScalarField& Ta = TaPtr_();
        const autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();

        const volVectorField gradT(fvc::grad(T));
        const volVectorField gradTa(fvc::grad(Ta));

        scalarField DsField(mesh_.nCells(), DSolid_);
        if (DSolidTablePtr_)
        {
            forAll(DsField, celli)
            {
                DsField[celli] = DSolidTablePtr_->value(T[celli]);
            }
        }
        scalarField thermSens
        (
            thermalSensScale_
           *(gradT.primitiveField() & gradTa.primitiveField())
           *(
                DsField - DFluid_
              - turbulence->nut()().primitiveField()/Prt_
            )
           *dt
        );

        const topOVariablesBase& vars =
            mesh_.lookupObject<topOVariablesBase>("topOVars");
        vars.sourceTermSensitivities
        (
            thermSens,
            kInterpolation_(),
            scalar(1),
            designVariablesName,
            "beta"
        );

        betaMult += thermSens;
    }
}


// ************************************************************************* //
