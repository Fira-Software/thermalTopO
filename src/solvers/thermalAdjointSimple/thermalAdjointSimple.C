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
    const volScalarField& T = TRef();

    const autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();

    // Same property state as the primal, by construction: both solvers build
    // their diffusivities from a thermalPropertyTables read from the same
    // "thermal" sub-dictionary.
    const volScalarField Df(props_.DFluidField(T, "DFluidFieldAdj"));
    const volScalarField Ds(props_.DSolidField(T, "DSolidFieldAdj"));

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
        Df
      + (Ds - Df)*Ik
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
    const volScalarField& T = TRef();

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

    // Adjoint convection. The transpose of the primal term C(T) u.grad(T) is
    // the conservative divergence -div(C u Ta), so the adjoint carries the
    // C-weighted flux -C_f phi_f.
    //
    // This flux is NOT solenoidal (div(C u) != 0 where C varies), so it must
    // use an unbounded convection scheme: a "bounded" scheme subtracts
    // Sp(div(C u), Ta), which would quietly replace the conservative operator
    // by the non-conservative C u.grad(Ta) and break adjoint consistency.
    // Hence its own scheme key, div(-phiC,Ta), which a variable-property case
    // must declare; the constant-property path keeps div(-phi,Ta) unchanged,
    // where div(u) = 0 makes bounded and conservative forms identical.
    tmp<surfaceScalarField> tNegPhiC;
    word divScheme;

    if (props_.variableRhoCp())
    {
        const volScalarField C(props_.CField(T, "rhoCpNorm"));

        tNegPhiC = -fvc::interpolate(C)*phi;
        divScheme = "div(-phiC,Ta)";
    }
    else
    {
        tNegPhiC = -phi;
        divScheme = "div(-phi,Ta)";
    }

    fvScalarMatrix TaEqn
    (
        fvm::div(tNegPhiC(), Ta, divScheme)
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

    // d/du of the primal convective term C(T) u.grad(T), contracted with Ta
    if (props_.variableRhoCp())
    {
        const volScalarField C(props_.CField(T, "rhoCpNorm"));

        matrix += couplingSign_*(C*Ta*fvc::grad(T))();
    }
    else
    {
        matrix += couplingSign_*(Ta*fvc::grad(T))();
    }
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
    props_(),
    Prt_(0.85),
    couplingSign_(1),
    thermalSensScale_(1),
    kInterpolation_(nullptr)
{
    const dictionary& thermalDict = dict.subDict("thermal");

    props_.read(thermalDict, mesh_);
    Prt_ = thermalDict.getOrDefault<scalar>("Prt", 0.85);
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
        props_.read(thermalDict, mesh_);
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
    //   thermalSensScale * Ik'(beta) (DSolid(T) - DFluid(T) - nut/Prt)
    //   (grad(T) . grad(Ta)) dt
    //
    // beta enters DEff only through Ik, so dDEff/dIk = Ds - Df - nut/Prt with
    // the properties evaluated at the frozen primal T. C(T) carries no beta
    // dependence: it multiplies convection, which the Brinkman term drives to
    // zero in solidified cells.
    if (mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        const volScalarField& T = TRef();
        const volScalarField& Ta = TaPtr_();
        const autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();

        const volVectorField gradT(fvc::grad(T));
        const volVectorField gradTa(fvc::grad(Ta));

        const volScalarField Df(props_.DFluidField(T, "DFluidFieldSens"));
        const volScalarField Ds(props_.DSolidField(T, "DSolidFieldSens"));

        scalarField thermSens
        (
            thermalSensScale_
           *(gradT.primitiveField() & gradTa.primitiveField())
           *(
                Ds.primitiveField() - Df.primitiveField()
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
