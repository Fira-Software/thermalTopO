/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "thermalAdjointSimple.H"
#include "topOVariablesBase.H"
#include "objectiveIncompressible.H"
#include "fixedGradientFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "sensitivityTopO.H"
#include "topOZones.H"
#include "bitSet.H"
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
    // extrapolatedCalculated: see the note in thermalSimple::kIndicator().
    // The adjoint must see the same DEff as the primal, boundaries included.
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
            dimensionedScalar(dimless, Zero),
            extrapolatedCalculatedFvPatchScalarField::typeName
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

            // Face diffusivity, not the adjacent cell value: the adjoint BC
            // D_face*snGrad(Ta) = -dJ/dT must use the SAME face coefficient
            // the primal Laplacian uses. patchInternalField() would be
            // inconsistent with it, and the discrepancy is real now that the
            // conductivity indicator is extrapolated to boundary faces.
            const scalarField DEffp(DEffField.boundaryField()[pI]);

            refCast<fixedGradientFvPatchScalarField>(Tap).gradient() =
                flux/max(DEffp, scalarField(DEffp.size(), SMALL));
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
        // rebuild: otherwise RAMP continuation through a dictionary re-read
        // would be silently ignored
        kInterpolation_ =
            topOInterpolationFunction::New
            (
                mesh_,
                thermalDict.subDict("kInterpolation")
            );

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
    //   thermalSensScale * Ik'(beta) * (DSolid(T) - DFluid(T) - nut/Prt)
    //                    * dJ/dDEff_P * dt
    //
    // beta enters DEff only through Ik, so dDEff/dIk = Ds - Df - nut/Prt with
    // the properties evaluated at the frozen primal T. C(T) carries no beta
    // dependence: it multiplies convection, which the Brinkman term drives to
    // zero in solidified cells.
    //
    // dJ/dDEff_P is computed FACE-BASED, i.e. by differentiating the discrete
    // operator that the primal actually solves, not its continuous analogue.
    // fvm::laplacian(DEff, T) assembles, for internal face f (owner P,
    // neighbour N),
    //
    //     R_P -= DEff_f * gamma_f * (T_N - T_P),   gamma_f = |Sf|*deltaCoeff
    //     R_N += DEff_f * gamma_f * (T_N - T_P)
    //
    // so, with L = J + sum_P Ta_P R_P,
    //
    //     dJ/dDEff_f = gamma_f * (T_N - T_P) * (Ta_N - Ta_P)
    //
    // and DEff_f = w_f*DEff_P + (1-w_f)*DEff_N (linear interpolation, which is
    // what "Gauss linear" means for the diffusivity). Distributing dJ/dDEff_f
    // onto the two cells by those same weights gives dJ/dDEff_P.
    //
    // In a smooth region this equals V_P*(grad(T).grad(Ta)) exactly, so it is a
    // drop-in for the continuous form used before. It differs precisely where
    // the continuous form is wrong: across a face with a large conductivity
    // jump, which is every fluid/solid interface in a topology optimisation
    // once the projection re-sharpens beta.
    //
    // It is retained because it matches the finite-volume operator rather than
    // its continuous analogue, and is necessary for hard material-interface
    // cases. It is NOT, on its own, sufficient to verify the
    // examples/coOptimiseChannelAndBody stress test: that case also exercises
    // the filter/projection/fixed-zone chain, and switching to the face-based
    // form left its adjoint/FD discrepancy unchanged. See that example's
    // README for the open issue.
    if (mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        const volScalarField& T = TRef();
        const volScalarField& Ta = TaPtr_();
        const autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();

        const volScalarField Df(props_.DFluidField(T, "DFluidFieldSens"));
        const volScalarField Ds(props_.DSolidField(T, "DSolidFieldSens"));

        // dJ/dDEff_P, accumulated face by face
        scalarField dJdDEff(mesh_.nCells(), Zero);

        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const surfaceScalarField& w = mesh_.weights();
        const surfaceScalarField& deltaCoeffs = mesh_.nonOrthDeltaCoeffs();
        const surfaceScalarField& magSf = mesh_.magSf();

        const scalarField& Ti = T.primitiveField();
        const scalarField& Tai = Ta.primitiveField();

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];

            const scalar dJdDf =
                magSf[facei]*deltaCoeffs[facei]
               *(Ti[N] - Ti[P])*(Tai[N] - Tai[P]);

            dJdDEff[P] += w[facei]*dJdDf;
            dJdDEff[N] += (scalar(1) - w[facei])*dJdDf;
        }

        // Boundary faces. The diffusive flux is DEff_b*|Sf|*snGrad(T), so
        // dJ/dDEff_b = -Ta_P*|Sf|*snGrad(T), and DEff_b depends on the
        // adjacent cell alone (the conductivity indicator is extrapolated).
        // Zero on adiabatic/zeroGradient patches; non-zero on a heated wall,
        // where it is the beta-dependence of the boundary heat input.
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];

            if (patch.coupled() || patch.size() == 0)
            {
                continue;
            }

            const labelUList& fc = patch.faceCells();
            const scalarField snGradT(T.boundaryField()[patchi].snGrad());
            const scalarField& magSfb = patch.magSf();

            forAll(fc, i)
            {
                dJdDEff[fc[i]] -= Tai[fc[i]]*magSfb[i]*snGradT[i];
            }
        }

        // Per unit volume, to keep the convention of the Brinkman term (and to
        // reduce to grad(T).grad(Ta) in the smooth limit)
        scalarField thermSens
        (
            thermalSensScale_
           *dJdDEff/mesh_.V()
           *(
                Ds.primitiveField() - Df.primitiveField()
              - turbulence->nut()().primitiveField()/Prt_
            )
           *dt
        );

        const topOVariablesBase& vars =
            mesh_.lookupObject<topOVariablesBase>("topOVars");

        // Mask cells the optimiser may not move (pinned solid/fluid zones,
        // anything outside the design space, and the cells next to inlets and
        // outlets). Upstream restricts the *update* to the active design
        // variables anyway, but this term is ours: zero it explicitly rather
        // than rely on that, so a stray sensitivity can never leak into a
        // fixed zone and so the written topOSens field is honest about where
        // the design can actually change.
        const topOZones& zones = vars.getTopOZones();

        auto applyMask = [&](scalarField& sens)
        {
            auto zeroZone = [&](const label zoneID)
            {
                if (zoneID != -1)
                {
                    for (const label celli : mesh_.cellZones()[zoneID])
                    {
                        sens[celli] = Zero;
                    }
                }
            };

            for (const label zoneID : zones.fixedPorousZoneIDs())
            {
                zeroZone(zoneID);
            }
            for (const label zoneID : zones.fixedZeroPorousZoneIDs())
            {
                zeroZone(zoneID);
            }
            for (const label celli : zones.IOCells())
            {
                sens[celli] = Zero;
            }

            // If an explicit design space is given, keep only those cells
            if (!zones.adjointPorousZoneIDs().empty())
            {
                bitSet isDesign(mesh_.nCells(), false);
                for (const label zoneID : zones.adjointPorousZoneIDs())
                {
                    for (const label celli : mesh_.cellZones()[zoneID])
                    {
                        isDesign.set(celli);
                    }
                }
                forAll(sens, celli)
                {
                    if (!isDesign.test(celli))
                    {
                        sens[celli] = Zero;
                    }
                }
            }
        };

        // Mask the physical sensitivity BEFORE the chain ...
        applyMask(thermSens);

        vars.sourceTermSensitivities
        (
            thermSens,
            kInterpolation_(),
            scalar(1),
            designVariablesName,
            "beta"
        );

        // ... and again AFTER it. sourceTermSensitivities only multiplies by
        // the interpolation derivative, so this is currently a no-op, but the
        // regularisation transpose applied downstream spreads values across
        // cells: masking on both sides of the chain keeps the guarantee true
        // regardless of what the chain does.
        applyMask(thermSens);

        betaMult += thermSens;
    }
}


// ************************************************************************* //
