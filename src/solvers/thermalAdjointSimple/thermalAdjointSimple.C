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
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvOptions.H"
#include "inletOutletFvPatchFields.H"
#include "fixedFluxExtrapolatedPressureFvPatchScalarField.H"
#include "linear.H"
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

    if (exactScalarThermalTranspose_)
    {
        volScalarField exactSource(source);
        scalarField& exactSourceI = exactSource.primitiveFieldRef();
        const scalarField& V = mesh_.V();

        for (objective& func : functions)
        {
            objectiveIncompressible& funcI =
                refCast<objectiveIncompressible>(func);

            if (!funcI.hasBoundarydJdT())
            {
                continue;
            }

            forAll(mesh_.boundary(), patchi)
            {
                if (mesh_.boundary()[patchi].type() == "empty")
                {
                    continue;
                }

                const fvPatch& patch = mesh_.boundary()[patchi];
                const scalarField& bdJdT = funcI.boundarydJdT(patchi);
                const fvsPatchScalarField& magSfp =
                    mesh_.magSf().boundaryField()[patchi];
                const labelUList& faceCells = patch.faceCells();

                tmp<scalarField> tValueInternalCoeffs =
                    T.boundaryField()[patchi].valueInternalCoeffs
                    (
                        patch.weights()
                    );
                const scalarField& valueInternalCoeffs =
                    tValueInternalCoeffs();

                forAll(bdJdT, facei)
                {
                    const label celli = faceCells[facei];
                    exactSourceI[celli] +=
                        func.weight()
                       *bdJdT[facei]
                       *magSfp[facei]
                       *valueInternalCoeffs[facei]
                       /V[celli];
                }
            }
        }

        volScalarField::Boundary& Tab = Ta.boundaryFieldRef();
        forAll(Tab, patchi)
        {
            if (isA<fixedGradientFvPatchScalarField>(Tab[patchi]))
            {
                refCast<fixedGradientFvPatchScalarField>
                (
                    Tab[patchi]
                ).gradient() = scalarField(Tab[patchi].size(), Zero);
            }
        }
        Ta.correctBoundaryConditions();

        volScalarField& TMutable = const_cast<volScalarField&>(T);
        tmp<fvScalarMatrix> tPrimalTEqn;

        if (props_.variableRhoCp())
        {
            const volScalarField C(props_.CField(T, "rhoCpNormTaExact"));

            tPrimalTEqn =
            (
                C.internalField()*fvm::div(phi, TMutable)
              - fvm::laplacian(DEff, TMutable)
            );
        }
        else
        {
            tPrimalTEqn =
            (
                fvm::div(phi, TMutable)
              - fvm::laplacian(DEff, TMutable)
            );
        }

        fvScalarMatrix& primalTEqn = tPrimalTEqn.ref();

        tmp<surfaceScalarField> tNegPhiCExact;
        word exactDivScheme;

        if (props_.variableRhoCp())
        {
            const volScalarField C(props_.CField(T, "rhoCpNormTaExactAdj"));

            tNegPhiCExact = -fvc::interpolate(C)*phi;
            exactDivScheme = "div(-phiC,Ta)";
        }
        else
        {
            tNegPhiCExact = -phi;
            exactDivScheme = "div(-phi,Ta)";
        }

        fvScalarMatrix TaEqn
        (
            fvm::div(tNegPhiCExact(), Ta, exactDivScheme)
          - fvm::laplacian(DEff, Ta)
        );

        TaEqn.diag() = primalTEqn.diag();
        TaEqn.upper() = primalTEqn.lower();
        TaEqn.lower() = primalTEqn.upper();
        TaEqn.internalCoeffs() = primalTEqn.internalCoeffs();
        TaEqn.source() = Zero;

        FieldField<Field, scalar>& bCoeffs = TaEqn.boundaryCoeffs();
        forAll(bCoeffs, patchi)
        {
            bCoeffs[patchi] = Zero;
        }

        Foam::solve(TaEqn == -exactSource);

        checkScalarThermalTranspose(DEff, source);

        Info<< "Min/max Ta: " << gMin(Ta.primitiveField()) << ", "
            << gMax(Ta.primitiveField()) << endl;

        return;
    }

    updateTaBCs(DEff);

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

    checkScalarThermalTranspose(DEff, source);

    Info<< "Min/max Ta: " << gMin(Ta.primitiveField()) << ", "
        << gMax(Ta.primitiveField()) << endl;
}


void Foam::thermalAdjointSimple::checkScalarThermalTranspose
(
    const volScalarField& DEffField,
    const volScalarField& source
)
{
    if (!checkScalarThermalTranspose_)
    {
        return;
    }

    const surfaceScalarField& phi = primalVars_.phi();
    const volScalarField& T = TRef();
    volScalarField& Ta = TaPtr_.ref();
    volScalarField& TMutable = const_cast<volScalarField&>(T);

    volScalarField dT
    (
        IOobject
        (
            "ATCTdTScalarTranspose",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T
    );
    dT = dimensionedScalar(T.dimensions(), Zero);

    scalarField& dTI = dT.primitiveFieldRef();
    forAll(dTI, celli)
    {
        dTI[celli] =
            scalar((53*celli + 29) % 127)/scalar(127) - scalar(0.5);
    }

    volScalarField::Boundary& dTb = dT.boundaryFieldRef();
    forAll(dTb, patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        if (isA<fixedGradientFvPatchScalarField>(dTb[patchi]))
        {
            refCast<fixedGradientFvPatchScalarField>(dTb[patchi]).gradient() =
                scalarField(dTb[patchi].size(), Zero);
        }

        if (dTb[patchi].fixesValue())
        {
            dTb[patchi] == scalar(0);
        }
    }

    dT.correctBoundaryConditions();

    forAll(dTb, patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() != "empty"
         && dTb[patchi].fixesValue()
        )
        {
            dTb[patchi] == scalar(0);
        }
    }

    tmp<fvScalarMatrix> tPrimalTEqn;

    if (props_.variableRhoCp())
    {
        const volScalarField C(props_.CField(T, "rhoCpNormScalarCheck"));

        tPrimalTEqn =
        (
            C.internalField()*fvm::div(phi, TMutable)
          - fvm::laplacian(DEffField, TMutable)
        );
    }
    else
    {
        tPrimalTEqn =
        (
            fvm::div(phi, TMutable)
          - fvm::laplacian(DEffField, TMutable)
        );
    }

    fvScalarMatrix& primalTEqn = tPrimalTEqn.ref();

    tmp<surfaceScalarField> tNegPhiC;
    word divScheme;

    if (props_.variableRhoCp())
    {
        const volScalarField C(props_.CField(T, "rhoCpNormScalarAdjCheck"));

        tNegPhiC = -fvc::interpolate(C)*phi;
        divScheme = "div(-phiC,Ta)";
    }
    else
    {
        tNegPhiC = -phi;
        divScheme = "div(-phi,Ta)";
    }

    fvScalarMatrix currentTaEqn
    (
        fvm::div(tNegPhiC(), Ta, divScheme)
      - fvm::laplacian(DEffField, Ta)
    );

    auto zeroKnownTerms = [](fvScalarMatrix& matrix)
    {
        matrix.source() = Zero;

        FieldField<Field, scalar>& bCoeffs = matrix.boundaryCoeffs();
        forAll(bCoeffs, patchi)
        {
            bCoeffs[patchi] = Zero;
        }
    };

    zeroKnownTerms(primalTEqn);
    zeroKnownTerms(currentTaEqn);

    auto integratedAction =
        [this](fvScalarMatrix& matrix, const volScalarField& field)
        {
            scalarField diag(matrix.diag());

            const FieldField<Field, scalar>& intCoeffs =
                matrix.internalCoeffs();
            forAll(intCoeffs, patchi)
            {
                const labelUList& faceCells =
                    mesh_.boundary()[patchi].faceCells();

                forAll(intCoeffs[patchi], facei)
                {
                    diag[faceCells[facei]] += intCoeffs[patchi][facei];
                }
            }

            scalarField action(diag*field.primitiveField());

            const scalarField& upper = matrix.upper();
            const scalarField& lower = matrix.lower();
            const labelUList& own = mesh_.owner();
            const labelUList& nei = mesh_.neighbour();

            forAll(own, facei)
            {
                action[own[facei]] += upper[facei]*field[nei[facei]];
                action[nei[facei]] += lower[facei]*field[own[facei]];
            }

            return action;
        };

    const scalarField primalActionIntegrated =
        integratedAction(primalTEqn, dT);
    const scalarField currentAdjointActionIntegrated =
        integratedAction(currentTaEqn, Ta);

    scalar objectiveVolume = Zero;
    scalar primalTerm = Zero;
    scalar currentAdjointTerm = Zero;

    const scalarField& V = mesh_.V();

    forAll(dTI, celli)
    {
        objectiveVolume += V[celli]*source[celli]*dT[celli];
        primalTerm += Ta[celli]*primalActionIntegrated[celli];
        currentAdjointTerm +=
            dT[celli]*currentAdjointActionIntegrated[celli];
    }

    scalar objectiveBoundary = Zero;
    PtrList<objective>& functions =
        objectiveManager_.getObjectiveFunctions();

    for (objective& func : functions)
    {
        objectiveIncompressible& funcI =
            refCast<objectiveIncompressible>(func);

        if (!funcI.hasBoundarydJdT())
        {
            continue;
        }

        forAll(mesh_.boundary(), patchi)
        {
            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            const scalarField& bdJdT = funcI.boundarydJdT(patchi);
            const fvPatchScalarField& dTp = dT.boundaryField()[patchi];
            const fvsPatchScalarField& magSfp =
                mesh_.magSf().boundaryField()[patchi];

            forAll(bdJdT, facei)
            {
                objectiveBoundary +=
                    func.weight()*bdJdT[facei]*dTp[facei]*magSfp[facei];
            }
        }
    }

    const scalar objectiveTerm = objectiveVolume + objectiveBoundary;

    boolList nearBoundary(mesh_.nCells(), false);
    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].size()
         && !mesh_.boundary()[patchi].coupled()
        )
        {
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(faceCells, facei)
            {
                nearBoundary[faceCells[facei]] = true;
            }
        }
    }

    scalar primalTermInterior = Zero;
    scalar currentAdjointTermInterior = Zero;
    scalar interiorScale = VSMALL;

    forAll(dTI, celli)
    {
        if (nearBoundary[celli])
        {
            continue;
        }

        const scalar pContrib =
            Ta[celli]*primalActionIntegrated[celli];
        const scalar bContrib =
            dT[celli]*currentAdjointActionIntegrated[celli];

        primalTermInterior += pContrib;
        currentAdjointTermInterior += bContrib;
        interiorScale += mag(pContrib) + mag(bContrib);
    }

    scalarField exactTransposeIntegrated(mesh_.nCells(), Zero);
    {
        scalarField diag(primalTEqn.diag());

        const FieldField<Field, scalar>& intCoeffs =
            primalTEqn.internalCoeffs();
        forAll(intCoeffs, patchi)
        {
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(intCoeffs[patchi], facei)
            {
                diag[faceCells[facei]] += intCoeffs[patchi][facei];
            }
        }

        const scalarField& upper = primalTEqn.upper();
        const scalarField& lower = primalTEqn.lower();
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();

        forAll(diag, celli)
        {
            exactTransposeIntegrated[celli] = diag[celli]*Ta[celli];
        }

        forAll(own, facei)
        {
            exactTransposeIntegrated[own[facei]] +=
                lower[facei]*Ta[nei[facei]];
            exactTransposeIntegrated[nei[facei]] +=
                upper[facei]*Ta[own[facei]];
        }
    }

    scalar exactTransposeTerm = Zero;
    scalar exactTransposeDiff = Zero;
    scalar exactTransposeScale = VSMALL;
    scalar exactTransposeTermInterior = Zero;
    scalar exactTransposeDiffInterior = Zero;
    scalar exactTransposeScaleInterior = VSMALL;

    forAll(dTI, celli)
    {
        const scalar exactContrib = dT[celli]*exactTransposeIntegrated[celli];
        const scalar currentContrib =
            dT[celli]*currentAdjointActionIntegrated[celli];

        exactTransposeTerm += exactContrib;
        exactTransposeDiff += mag(exactContrib - currentContrib);
        exactTransposeScale += mag(exactContrib);

        if (nearBoundary[celli])
        {
            continue;
        }

        exactTransposeTermInterior += exactContrib;
        exactTransposeDiffInterior += mag(exactContrib - currentContrib);
        exactTransposeScaleInterior += mag(exactContrib);
    }

    const scalar objPrimalScale =
        max(max(mag(objectiveTerm), mag(primalTerm)), VSMALL);
    const scalar objCurrentScale =
        max(max(mag(objectiveTerm), mag(currentAdjointTerm)), VSMALL);
    const scalar operatorScale =
        max(max(mag(primalTerm), mag(currentAdjointTerm)), VSMALL);
    const scalar exactInteriorScale =
        max
        (
            max(mag(exactTransposeTermInterior), mag(primalTermInterior)),
            VSMALL
        );
    const scalar exactScale =
        max(max(mag(exactTransposeTerm), mag(primalTerm)), VSMALL);

    Info<< "ATC-T scalar thermal transpose check: objective = "
        << objectiveTerm
        << " (volume " << objectiveVolume
        << ", boundary " << objectiveBoundary << ")"
        << ", TaT_A_dT = " << primalTerm
        << ", dTT_B_Ta = " << currentAdjointTerm
        << ", rel(objective + TaT_A_dT) = "
        << mag(objectiveTerm + primalTerm)/objPrimalScale
        << ", rel(objective + dTT_B_Ta) = "
        << mag(objectiveTerm + currentAdjointTerm)/objCurrentScale
        << ", rel(B_vs_Atranspose_all) = "
        << mag(currentAdjointTerm - primalTerm)/operatorScale
        << ", exactAtranspose = " << exactTransposeTerm
        << ", rel(exactCoeff_vs_primalAction_all) = "
        << mag(exactTransposeTerm - primalTerm)/exactScale
        << ", rel(B_vs_exactCoeff_L1_all) = "
        << exactTransposeDiff/max(exactTransposeScale, VSMALL)
        << endl;

    Info<< "ATC-T scalar thermal interior transpose check: TaT_A_dT = "
        << primalTermInterior
        << ", dTT_B_Ta = " << currentAdjointTermInterior
        << ", exactAtranspose = " << exactTransposeTermInterior
        << ", rel(B_vs_Atranspose) = "
        << mag(currentAdjointTermInterior - primalTermInterior)
          /max(interiorScale, VSMALL)
        << ", rel(exactCoeff_vs_primalAction) = "
        << mag(exactTransposeTermInterior - primalTermInterior)
          /exactInteriorScale
        << ", rel(B_vs_exactCoeff_L1) = "
        << exactTransposeDiffInterior
          /max(exactTransposeScaleInterior, VSMALL)
        << endl;
}


void Foam::thermalAdjointSimple::validateExactThermalCouplingOptions() const
{
    const bool needsExactScalar =
        couplingForm_ == "exactFluxTranspose"
     || checkFixedPointTangentSensitivity_
     || checkFixedPointTangentAgainstFD_
     || usePredictorReverseMomentumSens_
     || checkFixedPointMapAdjoint_
     || checkFullStateMapTranspose_;

    if (needsExactScalar && !exactScalarThermalTranspose_)
    {
        FatalErrorInFunction
            << "The exact thermal-flow coupling diagnostics require "
            << "exactScalarThermalTranspose true. The legacy "
            << "continuous/reversed-flux Ta equation is not the exact "
            << "discrete scalar adjoint, so thermalFluxSensitivity() cannot "
            << "be presented as an exact flux derivative with the current "
            << "settings. Enable exactScalarThermalTranspose, or disable "
            << "couplingForm exactFluxTranspose and the fixed-point/exact "
            << "thermal coupling diagnostics."
            << exit(FatalError);
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::thermalAdjointSimple::thermalFluxSensitivity() const
{
    const volScalarField& T = TRef();
    const volScalarField& Ta = TaPtr_();
    const surfaceScalarField& phi = primalVars_.phi();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    auto tGPhi =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "ATCTgPhi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(T.dimensions()*Ta.dimensions(), Zero)
        );
    surfaceScalarField& gPhi = tGPhi.ref();
    scalarField& gPhiI = gPhi.primitiveFieldRef();

    const scalar chi = boundedConvection_ ? 1 : 0;

    OStringStream divSchemeText;
    divSchemeText << mesh_.divScheme("div(phi,T)");
    const string divScheme(divSchemeText.str());
    const bool hasBounded =
        divScheme.find("bounded") != string::npos;
    const bool hasGauss =
        divScheme.find("Gauss") != string::npos;
    const bool hasUpwind =
        divScheme.find("upwind") != string::npos;
    const bool hasLinearUpwind =
        divScheme.find("linearUpwind") != string::npos;
    const bool supportedScheme =
        hasGauss
     && hasUpwind
     && !hasLinearUpwind
     && (hasBounded == boundedConvection_);

    if (!supportedScheme)
    {
        FatalErrorInFunction
            << "thermalFluxSensitivity() currently implements the exact "
            << "face derivative only for Gauss upwind, with the "
            << "boundedConvection switch matching the fvSchemes entry. "
            << "div(phi,T) scheme is " << mesh_.divScheme("div(phi,T)")
            << ", boundedConvection = " << boundedConvection_
            << ". Derive the scheme-specific flux derivative before using "
            << "this exact flux-transpose path with another scalar "
            << "convection scheme."
            << exit(FatalError);
    }

    tmp<volScalarField> tC;
    if (props_.variableRhoCp())
    {
        tC = props_.CField(T, "rhoCpNormATC");
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];

        const scalar CP = props_.variableRhoCp() ? tC()[P] : scalar(1);
        const scalar CN = props_.variableRhoCp() ? tC()[N] : scalar(1);

        const scalar AP = CP*Ta[P];
        const scalar AN = CN*Ta[N];

        const scalar Tf = (phi[facei] >= 0) ? T[P] : T[N];

        gPhiI[facei] = Tf*(AP - AN) - chi*(AP*T[P] - AN*T[N]);
    }

    surfaceScalarField::Boundary& gPhib = gPhi.boundaryFieldRef();
    forAll(gPhib, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (patch.type() == "empty")
        {
            continue;
        }

        fvsPatchScalarField& gPatch = gPhib[patchi];
        const fvPatchScalarField& TPatch = T.boundaryField()[patchi];
        const labelUList& faceCells = patch.faceCells();

        forAll(gPatch, facei)
        {
            const label celli = faceCells[facei];
            const scalar CP =
                props_.variableRhoCp() ? tC()[celli] : scalar(1);

            gPatch[facei] =
                CP*Ta[celli]*(TPatch[facei] - chi*T[celli]);
        }
    }

    return tGPhi;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::thermalAdjointSimple::momentumFluxSensitivity() const
{
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = getAdjointVars().Ua();
    const surfaceScalarField& phi = primalVars_.phi();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    auto tHPhi =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "ATCMomhPhi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(U.dimensions()*Ua.dimensions(), Zero)
        );
    surfaceScalarField& hPhi = tHPhi.ref();
    scalarField& hPhiI = hPhi.primitiveFieldRef();

    const scalar chi = boundedConvection_ ? 1 : 0;

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const vector Uf = (phi[facei] >= 0) ? U[P] : U[N];

        hPhiI[facei] =
            (Uf & (Ua[P] - Ua[N]))
          - chi*((Ua[P] & U[P]) - (Ua[N] & U[N]));
    }

    surfaceScalarField::Boundary& hPatches = hPhi.boundaryFieldRef();
    const surfaceScalarField::Boundary& phiPatches = phi.boundaryField();
    const volVectorField::Boundary& UPatches = U.boundaryField();

    forAll(hPatches, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (patch.type() == "empty" || patch.coupled())
        {
            continue;
        }

        fvsPatchScalarField& hPatch = hPatches[patchi];
        const fvsPatchScalarField& phiPatch = phiPatches[patchi];
        const fvPatchVectorField& UPatch = UPatches[patchi];
        const labelUList& fc = patch.faceCells();

        forAll(hPatch, facei)
        {
            const label celli = fc[facei];
            const vector Uf =
                (phiPatch[facei] >= 0) ? U[celli] : UPatch[facei];

            hPatch[facei] = Ua[celli] & (Uf - chi*U[celli]);
        }
    }

    return tHPhi;
}


Foam::tmp<Foam::volVectorField>
Foam::thermalAdjointSimple::exactMomentumATCSource()
{
    tmp<surfaceScalarField> tHPhi = momentumFluxSensitivity();
    return projectedFluxMomentumSource(tHPhi(), "ATCMom");
}


Foam::tmp<Foam::volVectorField>
Foam::thermalAdjointSimple::exactConvectionTransposeSource()
{
    const volVectorField& Ua = getAdjointVars().Ua();
    const volVectorField& U = primalVars_.U();
    const surfaceScalarField& phi = primalVars_.phi();

    tmp<volVectorField> tSource =
        tmp<volVectorField>::New
        (
            IOobject
            (
                "ATCExactConvTSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimLength/sqr(dimTime), Zero)
        );
    volVectorField& source = tSource.ref();
    vectorField& S = source.primitiveFieldRef();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& V = mesh_.V();
    const scalar chi = boundedConvection_ ? 1 : 0;

    vectorField exactIntegrated(mesh_.nCells(), vector::zero);

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar phif = phi[facei];

        if (phif >= 0)
        {
            exactIntegrated[P] += phif*(scalar(1) - chi)*Ua[P] - phif*Ua[N];
            exactIntegrated[N] += chi*phif*Ua[N];
        }
        else
        {
            exactIntegrated[P] += -chi*phif*Ua[P];
            exactIntegrated[N] += phif*Ua[P] - phif*(scalar(1) - chi)*Ua[N];
        }
    }

    const surfaceScalarField::Boundary& phib = phi.boundaryField();
    const volVectorField::Boundary& Ub = U.boundaryField();

    forAll(phib, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (patch.type() == "empty" || patch.coupled())
        {
            continue;
        }

        tmp<vectorField> tValueInternalCoeffs =
            Ub[patchi].valueInternalCoeffs(patch.weights());
        const vectorField& valueInternalCoeffs = tValueInternalCoeffs();
        const fvsPatchScalarField& phiPatch = phib[patchi];
        const labelUList& faceCells = patch.faceCells();

        forAll(phiPatch, facei)
        {
            const label celli = faceCells[facei];
            const scalar phif = phiPatch[facei];
            const vector coeff =
                (phif >= 0)
              ? vector::one
              : valueInternalCoeffs[facei];

            exactIntegrated[celli] +=
                phif*cmptMultiply(coeff - chi*vector::one, Ua[celli]);
        }
    }

    tmp<fvVectorMatrix> tBaseConv(fvm::div(-phi, const_cast<volVectorField&>(Ua)));
    fvVectorMatrix& baseConv = tBaseConv.ref();
    baseConv.boundaryManipulate
    (
        const_cast<volVectorField&>(Ua).boundaryFieldRef()
    );
    baseConv.source() = vector::zero;
    FieldField<Field, vector>& bCoeffs = baseConv.boundaryCoeffs();
    forAll(bCoeffs, patchi)
    {
        bCoeffs[patchi] = vector::zero;
    }

    tmp<volVectorField> tBaseAction =
        baseConv & static_cast<const DimensionedField<vector, volMesh>&>(Ua);
    const volVectorField& baseAction = tBaseAction();

    forAll(S, celli)
    {
        S[celli] = exactIntegrated[celli]/V[celli] - baseAction[celli];
    }

    return tSource;
}


Foam::tmp<Foam::volVectorField>
Foam::thermalAdjointSimple::projectedFluxMomentumSource
(
    const surfaceScalarField& gPhi,
    const word& namePrefix
)
{
    volScalarField& p = primalVars_.p();
    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
    fvOptions.constrain(UEqn);

    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA(namePrefix + "phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> trAtU(rAU);

    if (solverControl_().consistent())
    {
        trAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(trAtU() - rAU)*fvc::snGrad(p)*mesh_.magSf();
        HbyA -= (rAU - trAtU())*fvc::grad(p);
    }

    const volScalarField& rAtU = trAtU();

    constrainPressure(p, U, phiHbyA, rAtU);

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(solverControl_().pRefCell(), solverControl_().pRefValue());

    volScalarField lambda
    (
        IOobject
        (
            "paATCTProjection",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    lambda = dimensionedScalar(lambda.dimensions(), Zero);

    fvScalarMatrix lambdaEqn(fvm::laplacian(rAtU, lambda));
    scalarField& lambdaSource = lambdaEqn.source();
    lambdaSource = Zero;

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& lower = pEqn.lower();
    const scalarField& upper = pEqn.upper();
    const scalarField& gPhiI = gPhi.primitiveField();
    const FieldField<Field, scalar>& pInternalCoeffs =
        pEqn.internalCoeffs();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];

        lambdaSource[P] += -lower[facei]*gPhiI[facei];
        lambdaSource[N] +=  upper[facei]*gPhiI[facei];
    }

    forAll(pInternalCoeffs, patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        const scalarField& intCoeffs = pInternalCoeffs[patchi];
        const fvsPatchScalarField& gPhip = gPhi.boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(intCoeffs, facei)
        {
            lambdaSource[faceCells[facei]] += intCoeffs[facei]*gPhip[facei];
        }
    }

    lambdaEqn.setReference(solverControl_().pRefCell(), Zero);
    dictionary projectionSolver(lambdaEqn.solverDict("p"));
    projectionSolver.set("relTol", scalar(0));
    projectionSolver.set("tolerance", scalar(1e-12));
    lambdaEqn.solve(projectionSolver);

    surfaceScalarField gF
    (
        IOobject
        (
            namePrefix + "gF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gPhi
    );

    scalarField& gFI = gF.primitiveFieldRef();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];

        gFI[facei] -= lambda[P] - lambda[N];
    }

    surfaceScalarField::Boundary& gFb = gF.boundaryFieldRef();

    forAll(gFb, patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        fvsPatchScalarField& gFpatch = gFb[patchi];
        const fvsPatchScalarField& gPhip = gPhi.boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(gFpatch, facei)
        {
            gFpatch[facei] = gPhip[facei] - lambda[faceCells[facei]];
        }
    }

    auto tSrc =
        tmp<volVectorField>::New
        (
            IOobject
            (
                namePrefix + "source",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimLength/sqr(dimTime), Zero)
        );
    volVectorField& source = tSrc.ref();
    vectorField& S = source.primitiveFieldRef();

    const surfaceScalarField& w = mesh_.weights();
    const surfaceVectorField& Sf = mesh_.Sf();
    const scalarField& V = mesh_.V();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const vector Gf(gFI[facei]*Sf[facei]);

        S[P] += w[facei]*Gf/V[P];
        S[N] += (scalar(1) - w[facei])*Gf/V[N];
    }

    const surfaceVectorField::Boundary& Sfb = Sf.boundaryField();
    const volVectorField::Boundary& Ub = U.boundaryField();

    forAll(Ub, patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        tmp<vectorField> tValueInternalCoeffs =
            Ub[patchi].valueInternalCoeffs
            (
                mesh_.boundary()[patchi].weights()
            );
        const vectorField& valueInternalCoeffs = tValueInternalCoeffs();
        const fvsPatchScalarField& gFpatch = gF.boundaryField()[patchi];
        const fvsPatchVectorField& SfPatch = Sfb[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(gFpatch, facei)
        {
            const label celli = faceCells[facei];
            S[celli] +=
                gFpatch[facei]
               *cmptMultiply(valueInternalCoeffs[facei], SfPatch[facei])
               /V[celli];
        }
    }

    checkProjectedFluxTranspose(gPhi, gF, source, rAU, rAtU, pEqn);

    return tSrc;
}


Foam::tmp<Foam::volVectorField>
Foam::thermalAdjointSimple::computePredictorY
(
    const volVectorField& Uprobe,
    const word& namePrefix,
    const bool doRelax,
    const bool doFvOptionsConstrain
)
{
    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    const volVectorField Uold(U);

    U.primitiveFieldRef() = Uprobe.primitiveField();
    forAll(U.boundaryFieldRef(), patchi)
    {
        U.boundaryFieldRef()[patchi] == Uprobe.boundaryField()[patchi];
    }
    U.correctBoundaryConditions();

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    if (doRelax)
    {
        UEqn.relax();
    }

    if (doFvOptionsConstrain)
    {
        fvOptions.constrain(UEqn);
    }

    volScalarField rAU(1.0/UEqn.A());
    tmp<volVectorField> tY(rAU*UEqn.H());
    tY.ref().rename(namePrefix + "Y");

    U.primitiveFieldRef() = Uold.primitiveField();
    forAll(U.boundaryFieldRef(), patchi)
    {
        U.boundaryFieldRef()[patchi] == Uold.boundaryField()[patchi];
    }
    U.correctBoundaryConditions();

    return tY;
}


void Foam::thermalAdjointSimple::computePredictorYPieces
(
    const volVectorField& Uprobe,
    const word& namePrefix,
    volVectorField& Ytotal,
    volVectorField& Yldu,
    volVectorField& Ysource,
    const bool doRelax,
    const bool doFvOptionsConstrain,
    scalarField* rAUI,
    volVectorField* YsourceRaw,
    volVectorField* YsourceRelax,
    volVectorField* YsourceConstrain
)
{
    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    const volVectorField Uold(U);

    U.primitiveFieldRef() = Uprobe.primitiveField();
    forAll(U.boundaryFieldRef(), patchi)
    {
        U.boundaryFieldRef()[patchi] == Uprobe.boundaryField()[patchi];
    }
    U.correctBoundaryConditions();

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    const vectorField sourceRaw(UEqn.source());
    vectorField sourceAfterRelax(sourceRaw);

    if (doRelax)
    {
        UEqn.relax();
        sourceAfterRelax = UEqn.source();
    }

    vectorField sourceAfterConstrain(sourceAfterRelax);

    if (doFvOptionsConstrain)
    {
        fvOptions.constrain(UEqn);
        sourceAfterConstrain = UEqn.source();
    }

    volScalarField rAU(1.0/UEqn.A());
    if (rAUI)
    {
        *rAUI = rAU.primitiveField();
    }

    Ytotal = rAU*UEqn.H();
    Ytotal.rename(namePrefix + "Ytotal");

    Ysource.primitiveFieldRef() =
        rAU.primitiveField()*UEqn.source()/mesh_.V();
    Ysource.correctBoundaryConditions();
    Ysource.rename(namePrefix + "Ysource");

    if (YsourceRaw)
    {
        YsourceRaw->primitiveFieldRef() =
            rAU.primitiveField()*sourceRaw/mesh_.V();
        YsourceRaw->correctBoundaryConditions();
        YsourceRaw->rename(namePrefix + "YsourceRaw");
    }

    if (YsourceRelax)
    {
        YsourceRelax->primitiveFieldRef() =
            rAU.primitiveField()*(sourceAfterRelax - sourceRaw)/mesh_.V();
        YsourceRelax->correctBoundaryConditions();
        YsourceRelax->rename(namePrefix + "YsourceRelax");
    }

    if (YsourceConstrain)
    {
        YsourceConstrain->primitiveFieldRef() =
            rAU.primitiveField()
           *(sourceAfterConstrain - sourceAfterRelax)/mesh_.V();
        YsourceConstrain->correctBoundaryConditions();
        YsourceConstrain->rename(namePrefix + "YsourceConstrain");
    }

    Yldu = Ytotal;
    Yldu.primitiveFieldRef() -= Ysource.primitiveField();
    Yldu.correctBoundaryConditions();
    Yldu.rename(namePrefix + "Yldu");

    U.primitiveFieldRef() = Uold.primitiveField();
    forAll(U.boundaryFieldRef(), patchi)
    {
        U.boundaryFieldRef()[patchi] == Uold.boundaryField()[patchi];
    }
    U.correctBoundaryConditions();
}


void Foam::thermalAdjointSimple::checkPredictorMapTranspose
(
    const vectorField& barYInt,
    const word& checkName,
    const bool doRelax,
    const bool doFvOptionsConstrain
)
{
    const volVectorField& U = primalVars_.U();

    const label nTarget = mesh_.nCells() < 12 ? mesh_.nCells() : 12;
    labelList selected(nTarget, -1);
    scalarField selectedMag(nTarget, -GREAT);

    forAll(barYInt, celli)
    {
        const scalar m = mag(barYInt[celli]);

        for (label slot = 0; slot < nTarget; ++slot)
        {
            if (m > selectedMag[slot])
            {
                for (label shift = nTarget - 1; shift > slot; --shift)
                {
                    selected[shift] = selected[shift - 1];
                    selectedMag[shift] = selectedMag[shift - 1];
                }

                selected[slot] = celli;
                selectedMag[slot] = m;
                break;
            }
        }
    }

    label nSelected = 0;
    forAll(selected, i)
    {
        if (selected[i] != -1)
        {
            ++nSelected;
        }
    }

    if (!nSelected)
    {
        Info<< "ATC-T predictor FD transpose check (" << checkName
            << "): skipped, no selected cells" << endl;
        return;
    }

    volVectorField dU
    (
        IOobject
        (
            "ATCTdUPredictorFDTranspose",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    dU = dimensionedVector(dU.dimensions(), vector::zero);

    for (label i = 0; i < nSelected; ++i)
    {
        const label celli = selected[i];
        dU.primitiveFieldRef()[celli] =
            vector
            (
                scalar((53*celli + 23) % 127)/scalar(127) - scalar(0.5),
                scalar((59*celli + 29) % 131)/scalar(131) - scalar(0.5),
                scalar((61*celli + 31) % 137)/scalar(137) - scalar(0.5)
            );
    }
    dU.correctBoundaryConditions();

    forAll(dU.boundaryFieldRef(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() != "empty"
         && !U.boundaryField()[patchi].assignable()
        )
        {
            dU.boundaryFieldRef()[patchi] == vector::zero;
        }
    }

    const scalar maxU =
        max(gMax(mag(U.primitiveField())), scalar(1));
    const scalar maxdU =
        max(gMax(mag(dU.primitiveField())), SMALL);
    const scalar eps = scalar(1e-6)*maxU/maxdU;

    volVectorField Uplus
    (
        IOobject
        (
            "ATCTUplusPredictorFDTranspose",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField Uminus
    (
        IOobject
        (
            "ATCTUminusPredictorFDTranspose",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );

    Uplus.primitiveFieldRef() =
        U.primitiveField() + eps*dU.primitiveField();
    Uminus.primitiveFieldRef() =
        U.primitiveField() - eps*dU.primitiveField();

    forAll(Uplus.boundaryFieldRef(), patchi)
    {
        Uplus.boundaryFieldRef()[patchi] ==
            U.boundaryField()[patchi] + eps*dU.boundaryField()[patchi];
        Uminus.boundaryFieldRef()[patchi] ==
            U.boundaryField()[patchi] - eps*dU.boundaryField()[patchi];
    }
    Uplus.correctBoundaryConditions();
    Uminus.correctBoundaryConditions();

    scalarField rAUBase(mesh_.nCells(), Zero);
    scalarField rAUPlus(mesh_.nCells(), Zero);
    scalarField rAUMinus(mesh_.nCells(), Zero);

    volVectorField YtotalBase
    (
        IOobject
        (
            word("ATCTPredBaseTotal") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField YlduBase
    (
        IOobject
        (
            word("ATCTPredBaseLdu") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField YsourceBase
    (
        IOobject
        (
            word("ATCTPredBaseSource") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );

    volVectorField YtotalPlus
    (
        IOobject
        (
            word("ATCTPredPlusTotal") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField YlduPlus
    (
        IOobject
        (
            word("ATCTPredPlusLdu") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField YsourcePlus
    (
        IOobject
        (
            word("ATCTPredPlusSource") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField YtotalMinus
    (
        IOobject
        (
            word("ATCTPredMinusTotal") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField YlduMinus
    (
        IOobject
        (
            word("ATCTPredMinusLdu") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField YsourceMinus
    (
        IOobject
        (
            word("ATCTPredMinusSource") + checkName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField YsourceRawPlus(YsourcePlus);
    volVectorField YsourceRelaxPlus(YsourcePlus);
    volVectorField YsourceConstrainPlus(YsourcePlus);
    volVectorField YsourceRawMinus(YsourceMinus);
    volVectorField YsourceRelaxMinus(YsourceMinus);
    volVectorField YsourceConstrainMinus(YsourceMinus);

    computePredictorYPieces
    (
        U,
        word("ATCTPredBase") + checkName,
        YtotalBase,
        YlduBase,
        YsourceBase,
        doRelax,
        doFvOptionsConstrain,
        &rAUBase
    );

    computePredictorYPieces
    (
        Uplus,
        word("ATCTPredPlus") + checkName,
        YtotalPlus,
        YlduPlus,
        YsourcePlus,
        doRelax,
        doFvOptionsConstrain,
        &rAUPlus,
        &YsourceRawPlus,
        &YsourceRelaxPlus,
        &YsourceConstrainPlus
    );
    computePredictorYPieces
    (
        Uminus,
        word("ATCTPredMinus") + checkName,
        YtotalMinus,
        YlduMinus,
        YsourceMinus,
        doRelax,
        doFvOptionsConstrain,
        &rAUMinus,
        &YsourceRawMinus,
        &YsourceRelaxMinus,
        &YsourceConstrainMinus
    );

    vectorField dY
    (
        (YtotalPlus.primitiveField() - YtotalMinus.primitiveField())/(2*eps)
    );
    vectorField dYldu
    (
        (YlduPlus.primitiveField() - YlduMinus.primitiveField())/(2*eps)
    );
    vectorField dYsource
    (
        (
            YsourcePlus.primitiveField()
          - YsourceMinus.primitiveField()
        )/(2*eps)
    );
    vectorField dYsourceRaw
    (
        (
            YsourceRawPlus.primitiveField()
          - YsourceRawMinus.primitiveField()
        )/(2*eps)
    );
    vectorField dYsourceRelax
    (
        (
            YsourceRelaxPlus.primitiveField()
          - YsourceRelaxMinus.primitiveField()
        )/(2*eps)
    );
    vectorField dYsourceConstrain
    (
        (
            YsourceConstrainPlus.primitiveField()
          - YsourceConstrainMinus.primitiveField()
        )/(2*eps)
    );

    scalar lhs = Zero;
    scalar lhsLdu = Zero;
    scalar lhsSource = Zero;
    scalar lhsSourceRaw = Zero;
    scalar lhsSourceRelax = Zero;
    scalar lhsSourceConstrain = Zero;
    scalar lhsLduFixedA = Zero;
    scalar lhsSourceFixedA = Zero;
    scalar lhsDiag = Zero;
    forAll(dY, celli)
    {
        const scalar rPlus =
            mag(rAUPlus[celli]) > VSMALL
          ? rAUBase[celli]/rAUPlus[celli]
          : Zero;
        const scalar rMinus =
            mag(rAUMinus[celli]) > VSMALL
          ? rAUBase[celli]/rAUMinus[celli]
          : Zero;

        const vector dYlduFixedA =
            (rPlus*YlduPlus[celli] - rMinus*YlduMinus[celli])/(2*eps);
        const vector dYsourceFixedA =
            (
                rPlus*YsourcePlus[celli]
              - rMinus*YsourceMinus[celli]
            )/(2*eps);
        const vector dYFixedA = dYlduFixedA + dYsourceFixedA;

        lhs += barYInt[celli] & dY[celli];
        lhsLdu += barYInt[celli] & dYldu[celli];
        lhsSource += barYInt[celli] & dYsource[celli];
        lhsSourceRaw += barYInt[celli] & dYsourceRaw[celli];
        lhsSourceRelax += barYInt[celli] & dYsourceRelax[celli];
        lhsSourceConstrain += barYInt[celli] & dYsourceConstrain[celli];
        lhsLduFixedA += barYInt[celli] & dYlduFixedA;
        lhsSourceFixedA += barYInt[celli] & dYsourceFixedA;
        lhsDiag += barYInt[celli] & (dY[celli] - dYFixedA);
    }

    vectorField barUFD(mesh_.nCells(), vector::zero);
    const scalar epsComp = scalar(1e-6)*maxU;

    for (label i = 0; i < nSelected; ++i)
    {
        const label cellj = selected[i];

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            volVectorField UcompPlus
            (
                IOobject
                (
                    "ATCTUCompPlusPredictorFDTranspose",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                U
            );
            volVectorField UcompMinus
            (
                IOobject
                (
                    "ATCTUCompMinusPredictorFDTranspose",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                U
            );

            UcompPlus.primitiveFieldRef()[cellj][cmpt] += epsComp;
            UcompMinus.primitiveFieldRef()[cellj][cmpt] -= epsComp;
            UcompPlus.correctBoundaryConditions();
            UcompMinus.correctBoundaryConditions();

            tmp<volVectorField> tYcompPlus =
                computePredictorY
                (
                    UcompPlus,
                    word("ATCTPredCompPlus") + checkName,
                    doRelax,
                    doFvOptionsConstrain
                );
            tmp<volVectorField> tYcompMinus =
                computePredictorY
                (
                    UcompMinus,
                    word("ATCTPredCompMinus") + checkName,
                    doRelax,
                    doFvOptionsConstrain
                );

            scalar contraction = Zero;
            forAll(barYInt, celli)
            {
                contraction +=
                    barYInt[celli]
                  & (
                        tYcompPlus().primitiveField()[celli]
                      - tYcompMinus().primitiveField()[celli]
                    )/(2*epsComp);
            }

            barUFD[cellj][cmpt] = contraction;
        }
    }

    scalar rhs = Zero;
    for (label i = 0; i < nSelected; ++i)
    {
        const label celli = selected[i];
        rhs += barUFD[celli] & dU[celli];
    }

    tmp<vectorField> tBarUAnalytic =
        reversePredictorY(barYInt, doRelax, doFvOptionsConstrain);
    const vectorField& barUAnalytic = tBarUAnalytic();

    scalar rhsAnalytic = Zero;
    for (label i = 0; i < nSelected; ++i)
    {
        const label celli = selected[i];
        rhsAnalytic += barUAnalytic[celli] & dU[celli];
    }

    const scalar scale = max(max(mag(lhs), mag(rhs)), VSMALL);

    Info<< "ATC-T predictor FD transpose check (" << checkName
        << ", relax=" << doRelax
        << ", fvOptions=" << doFvOptionsConstrain
        << "): lhs = " << lhs
        << ", rhs = " << rhs
        << ", rel = " << mag(lhs - rhs)/scale
        << ", rhsAnalytic = " << rhsAnalytic
        << ", relAnalytic = " << mag(lhs - rhsAnalytic)/scale
        << ", lhsLdu = " << lhsLdu
        << ", lhsSource = " << lhsSource
        << ", sourceFrac = "
        << mag(lhsSource)/max(mag(lhs), VSMALL)
        << ", lhsLduFixedA = " << lhsLduFixedA
        << ", lhsSourceFixedA = " << lhsSourceFixedA
        << ", lhsDiag = " << lhsDiag
        << ", diagFrac = " << mag(lhsDiag)/max(mag(lhs), VSMALL)
        << ", lhsSourceRaw = " << lhsSourceRaw
        << ", lhsSourceRelax = " << lhsSourceRelax
        << ", lhsSourceConstrain = " << lhsSourceConstrain
        << ", eps = " << eps
        << ", epsComp = " << epsComp
        << ", nSelected = " << nSelected
        << endl;
}


Foam::tmp<Foam::vectorField> Foam::thermalAdjointSimple::reversePredictorY
(
    const vectorField& barYInt,
    const bool doRelax,
    const bool doFvOptionsConstrain
)
{
    PredictorReverseSeed seed =
        reversePredictorYFullState
        (
            barYInt,
            doRelax,
            doFvOptionsConstrain
        );

    auto tBarU = tmp<vectorField>::New(mesh_.nCells(), vector::zero);
    tBarU.ref() = seed.barUold;

    return tBarU;
}


Foam::vectorField Foam::thermalAdjointSimple::reverseExplicitStokesViscousSource
(
    const vectorField& barSourceIntegrated
)
{
    volVectorField& U = primalVars_.U();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();

    if (turbulence->type() != "Stokes")
    {
        FatalErrorInFunction
            << "Exact explicit viscous predictor transpose is implemented "
            << "only for the supported Stokes laminar model. Runtime model is "
            << turbulence->type() << "."
            << exit(FatalError);
    }

    Info<< "ATC-T explicit Stokes source reverse uses "
        << "$FOAM_SRC/TurbulenceModels/turbulenceModels/"
        << "linearViscousStress/linearViscousStress.C::divDevRhoReff, "
        << "-fvc::div(nuEff*dev2(T(fvc::grad(U)))) with frozen nuEff."
        << endl;

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].coupled() && mesh_.boundary()[patchi].size())
        {
            FatalErrorInFunction
                << "Exact explicit viscous predictor transpose currently "
                << "supports the serial non-coupled diagnostic case only. "
                << "Coupled patch: " << mesh_.boundary()[patchi].name()
                << exit(FatalError);
        }
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Sf = mesh_.Sf();
    const surfaceScalarField& weights = mesh_.weights();
    const scalarField& V = mesh_.V();

    tmp<volScalarField> tNuEff = turbulence->nuEff();
    const volScalarField& nuEff = tNuEff();

    tensorField barTau(mesh_.nCells(), tensor::zero);
    List<tensorField> barTauBoundary(mesh_.boundary().size());
    forAll(barTauBoundary, patchi)
    {
        barTauBoundary[patchi].setSize
        (
            mesh_.boundary()[patchi].size(),
            tensor::zero
        );
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar w = weights[facei];

        const vector barFaceFlux =
            barSourceIntegrated[P] - barSourceIntegrated[N];

        const tensor barTauFace = Sf[facei]*barFaceFlux;

        barTau[P] += w*barTauFace;
        barTau[N] += (scalar(1) - w)*barTauFace;
    }

    // Boundary faces of fvc::div(tau) depend on the tau patch field, not on
    // the owner-cell tau by interpolation.  The exact patch-gradient reverse
    // is handled separately for the supported scratch scope; do not route this
    // seed through owner-cell tau as if it were an internal interpolation.
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (patch.type() == "empty")
        {
            continue;
        }

        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const labelUList& faceCells = patch.faceCells();

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            barTauBoundary[patchi][facei] +=
                Sfp[facei]*barSourceIntegrated[celli];
        }
    }

    tensorField barGrad(mesh_.nCells(), tensor::zero);
    forAll(barGrad, celli)
    {
        const tensor C = nuEff[celli]*barTau[celli];
        barGrad[celli] += T(C) - (scalar(2)/scalar(3))*tr(C)*tensor::I;
    }

    vectorField barU(mesh_.nCells(), vector::zero);

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (patch.type() == "empty")
        {
            continue;
        }

        const vectorField n
        (
            mesh_.Sf().boundaryField()[patchi]
          / mesh_.magSf().boundaryField()[patchi]
        );
        const scalarField& nuEffp = nuEff.boundaryField()[patchi];
        const labelUList& faceCells = patch.faceCells();
        tmp<vectorField> tGradientInternalCoeffs =
            U.boundaryField()[patchi].gradientInternalCoeffs();
        const vectorField& gradientInternalCoeffs =
            tGradientInternalCoeffs();

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            const tensor C = nuEffp[facei]*barTauBoundary[patchi][facei];
            const tensor barGradCorrected =
                T(C) - (scalar(2)/scalar(3))*tr(C)*tensor::I;

            const vector barSnGrad = n[facei] & barGradCorrected;
            const tensor barGradUncorrected =
                barGradCorrected - n[facei]*barSnGrad;

            barGrad[celli] += barGradUncorrected;

            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                barU[celli][cmpt] +=
                    gradientInternalCoeffs[facei][cmpt]
                   *barSnGrad[cmpt];
            }
        }
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar w = weights[facei];

        const vector barUf =
            (Sf[facei] & barGrad[P])/V[P]
          - (Sf[facei] & barGrad[N])/V[N];

        barU[P] += w*barUf;
        barU[N] += (scalar(1) - w)*barUf;
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (patch.type() == "empty")
        {
            continue;
        }

        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const labelUList& faceCells = patch.faceCells();
        tmp<vectorField> tValueInternalCoeffs =
            U.boundaryField()[patchi].valueInternalCoeffs
            (
                patch.weights()
            );
        const vectorField& valueInternalCoeffs = tValueInternalCoeffs();

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            const vector barUb = (Sfp[facei] & barGrad[celli])/V[celli];

            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                barU[celli][cmpt] +=
                    valueInternalCoeffs[facei][cmpt]*barUb[cmpt];
            }
        }
    }

    return barU;
}


Foam::thermalAdjointSimple::PredictorReverseSeed
Foam::thermalAdjointSimple::reversePredictorYFullState
(
    const vectorField& barYInt,
    const bool doRelax,
    const bool doFvOptionsConstrain
)
{
    scalarField zeroD(mesh_.nCells(), Zero);
    scalarField zeroH1(mesh_.nCells(), Zero);

    return reversePredictorYFullState
    (
        barYInt,
        zeroD,
        zeroH1,
        doRelax,
        doFvOptionsConstrain
    );
}


Foam::thermalAdjointSimple::FinalDiagonalPhiSeed
Foam::thermalAdjointSimple::reverseExtraFinalDiagonalToPhi
(
    const scalarField& extraBarDfinal,
    const scalar alphaU
) const
{
    if (extraBarDfinal.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "extraBarDfinal size mismatch: "
            << extraBarDfinal.size() << " vs " << mesh_.nCells()
            << exit(FatalError);
    }

    const volVectorField& U = primalVars_.U();
    const surfaceScalarField& phi = primalVars_.phi();
    const scalar chi = boundedConvection_ ? scalar(1) : scalar(0);

    FinalDiagonalPhiSeed result(mesh_);

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar phif = phi[facei];
        scalar dDPhiP = Zero;
        scalar dDPhiN = Zero;

        if (phif >= scalar(0))
        {
            dDPhiP = scalar(1) - chi;
            dDPhiN = chi;
        }
        else
        {
            dDPhiP = -chi;
            dDPhiN = -scalar(1) + chi;
        }

        result.internal[facei] +=
            (extraBarDfinal[P]*dDPhiP + extraBarDfinal[N]*dDPhiN)/alphaU;
    }

    scalar maxMagPhi = gMax(mag(phi.primitiveField()));
    forAll(phi.boundaryField(), patchi)
    {
        maxMagPhi = max(maxMagPhi, gMax(mag(phi.boundaryField()[patchi])));
    }
    const scalar nearZeroPhi = scalar(1e-12)*max(maxMagPhi, scalar(1));

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        if (mesh_.boundary()[patchi].coupled())
        {
            FatalErrorInFunction
                << "reverseExtraFinalDiagonalToPhi supports the current "
                << "serial non-coupled diagnostic case only. Coupled patch "
                << mesh_.boundary()[patchi].name()
                << " is active."
                << exit(FatalError);
        }

        const fvsPatchScalarField& phip = phi.boundaryField()[patchi];
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();
        scalarField& barPhip = result.boundary[patchi];

        tmp<vectorField> tValueInternalCoeffs =
            Up.valueInternalCoeffs(mesh_.boundary()[patchi].weights());
        const vectorField& valueInternalCoeffs = tValueInternalCoeffs();

        typename pTraits<vector>::labelType validVectorComponents
        (
            mesh_.validComponents<vector>()
        );

        forAll(phip, facei)
        {
            const label celli = faceCells[facei];
            const scalar phif = phip[facei];

            if
            (
                Up.type() == "inletOutlet"
             && phif > nearZeroPhi
             && boundedConvection_
            )
            {
                for (label cmpt = 0; cmpt < vector::nComponents; ++cmpt)
                {
                    if (component(validVectorComponents, cmpt) == -1)
                    {
                        continue;
                    }

                    if
                    (
                        mag(valueInternalCoeffs[facei][cmpt] - scalar(1))
                      > scalar(1e-10)
                    )
                    {
                        FatalErrorInFunction
                            << "Unsupported inletOutlet final-diagonal "
                            << "coefficient on patch "
                            << mesh_.boundary()[patchi].name()
                            << " face " << facei
                            << ": valueInternalCoeffs = "
                            << valueInternalCoeffs[facei]
                            << exit(FatalError);
                    }
                }

                // OpenFOAM relax() subtracts cmptMin(internalCoeffs) after
                // the max/alpha operation.  For this outlet branch the
                // face-basis diagnostics prove d(maxSelected)/dphi = 0 and
                // d(cmptMin(internalCoeffs))/dphi = 1 on valid components.
                barPhip[facei] -= extraBarDfinal[celli];
            }
            else
            {
                if (mag(phif) <= nearZeroPhi)
                {
                    continue;
                }

                const scalar dDPhi =
                    (phif >= scalar(0)) ? (scalar(1) - chi) : -chi;
                barPhip[facei] += extraBarDfinal[celli]*dDPhi/alphaU;
            }
        }
    }

    return result;
}


Foam::scalarField
Foam::thermalAdjointSimple::reversePressureRAtUExact
(
    const volScalarField& rAtUBase,
    const volScalarField& pSolveBase,
    const surfaceScalarField& barF2
)
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "reversePressureRAtUExact currently supports only serial "
            << "diagnostic runs." << exit(FatalError);
    }
    if (rAtUBase.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "rAtUBase size mismatch: " << rAtUBase.size()
            << " vs " << mesh_.nCells() << exit(FatalError);
    }
    if (pSolveBase.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "pSolveBase size mismatch: " << pSolveBase.size()
            << " vs " << mesh_.nCells() << exit(FatalError);
    }
    if (barF2.primitiveField().size() != mesh_.nInternalFaces())
    {
        FatalErrorInFunction
            << "barF2 internal size mismatch: "
            << barF2.primitiveField().size()
            << " vs " << mesh_.nInternalFaces() << exit(FatalError);
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (patch.type() == "empty")
        {
            continue;
        }
        if (patch.coupled())
        {
            FatalErrorInFunction
                << "reversePressureRAtUExact supports only non-coupled "
                << "patches. Patch " << patch.name() << " is coupled."
                << exit(FatalError);
        }

        const fvPatchScalarField& rPatch = rAtUBase.boundaryField()[patchi];
        if (rPatch.type() != "calculated")
        {
            FatalErrorInFunction
                << "reversePressureRAtUExact only supports the validated "
                << "calculated rAtU boundary branch. Patch " << patch.name()
                << " has rAtU type " << rPatch.type() << exit(FatalError);
        }
    }

    surfaceScalarField gammaFaceUnit
    (
        IOobject
        (
            "ATCTPressureRAtUExactGammaUnit",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(rAtUBase)
    );
    gammaFaceUnit.primitiveFieldRef() = scalar(1);
    forAll(gammaFaceUnit.boundaryFieldRef(), patchi)
    {
        gammaFaceUnit.boundaryFieldRef()[patchi] == scalar(1);
    }

    fvScalarMatrix unitEqn(fvm::laplacian(gammaFaceUnit, pSolveBase));
    unitEqn.setReference(solverControl_().pRefCell(), solverControl_().pRefValue());

    surfaceScalarField unitFlux
    (
        IOobject
        (
            "ATCTPressureRAtUExactUnitFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        unitEqn.flux()
    );

    scalarField barrAtU(mesh_.nCells(), Zero);
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceScalarField& w = mesh_.weights();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar barGamma =
           -barF2.primitiveField()[facei]*unitFlux.primitiveField()[facei];

        barrAtU[P] += w[facei]*barGamma;
        barrAtU[N] += (scalar(1) - w[facei])*barGamma;
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        // In the validated active branch, the surface-coefficient field built
        // by linearInterpolate(rAtU) has zero boundary tangent with respect to
        // the stored cell-centred rAtU state.  Boundary barGamma is therefore
        // deliberately not mapped to owner cells here.
    }

    return barrAtU;
}


Foam::thermalAdjointSimple::PredictorReverseSeed
Foam::thermalAdjointSimple::reversePredictorYFullState
(
    const vectorField& barYInt,
    const scalarField& extraBarDfinal,
    const scalarField& extraBarUEqnH1Coeff,
    const bool doRelax,
    const bool doFvOptionsConstrain
)
{
    List<vectorField> zeroBoundary(mesh_.boundary().size());
    forAll(zeroBoundary, patchi)
    {
        zeroBoundary[patchi].setSize
        (
            mesh_.boundary()[patchi].size(),
            vector::zero
        );
    }

    return reversePredictorYFullState
    (
        barYInt,
        zeroBoundary,
        extraBarDfinal,
        extraBarUEqnH1Coeff,
        doRelax,
        doFvOptionsConstrain
    );
}


Foam::thermalAdjointSimple::PredictorReverseSeed
Foam::thermalAdjointSimple::reversePredictorYFullState
(
    const vectorField& barYInt,
    const List<vectorField>& barYBoundary,
    const scalarField& extraBarDfinal,
    const scalarField& extraBarUEqnH1Coeff,
    const bool doRelax,
    const bool doFvOptionsConstrain
)
{
    volVectorField& U = primalVars_.U();
    const surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    vectorField barYWork(barYInt);
    if (barYBoundary.size() != mesh_.boundary().size())
    {
        FatalErrorInFunction
            << "barYBoundary list size mismatch: "
            << barYBoundary.size() << " vs "
            << mesh_.boundary().size()
            << exit(FatalError);
    }
    forAll(barYBoundary, patchi)
    {
        if (barYBoundary[patchi].size() != mesh_.boundary()[patchi].size())
        {
            FatalErrorInFunction
                << "barYBoundary patch size mismatch for patch "
                << mesh_.boundary()[patchi].name() << ": "
                << barYBoundary[patchi].size() << " vs "
                << mesh_.boundary()[patchi].size()
                << exit(FatalError);
        }
    }

    typename pTraits<vector>::labelType validVectorComponents
    (
        mesh_.validComponents<vector>()
    );
    forAll(barYWork, celli)
    {
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            if (component(validVectorComponents, cmpt) == -1)
            {
                barYWork[celli][cmpt] = Zero;
            }
        }
    }

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    const scalarField diagRaw(UEqn.diag());
    const scalarField lowerRaw(UEqn.lower());
    const scalarField upperRaw(UEqn.upper());
    const FieldField<Field, vector> internalCoeffsRaw(UEqn.internalCoeffs());
    const vectorField sourceRaw(UEqn.source());
    scalarField diagRelax(diagRaw);
    vectorField sourceRelax(sourceRaw);

    if (doRelax)
    {
        UEqn.relax();
        diagRelax = UEqn.diag();
        sourceRelax = UEqn.source();
    }

    scalarField diagFinal(diagRelax);
    vectorField sourceFinal(sourceRelax);

    if (doFvOptionsConstrain)
    {
        fvOptions.constrain(UEqn);
        diagFinal = UEqn.diag();
        sourceFinal = UEqn.source();
    }

    volScalarField rAU(1.0/UEqn.A());
    volVectorField Y
    (
        IOobject
        (
            "ATCTPredictorReverseY",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rAU*UEqn.H()
    );

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        const fvPatchVectorField& Yp = Y.boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        if (Yp.type() == fvPatchFieldBase::extrapolatedCalculatedType())
        {
            forAll(faceCells, facei)
            {
                barYWork[faceCells[facei]] +=
                    barYBoundary[patchi][facei];
            }
        }
        else if (Yp.type() == "calculated")
        {
            FatalErrorInFunction
                << "Unsupported calculated Y patch without an "
                << "extrapolatedCalculated value rule on patch "
                << mesh_.boundary()[patchi].name()
                << exit(FatalError);
        }
        else
        {
            tmp<vectorField> tValueInternalCoeffs =
                Yp.valueInternalCoeffs(mesh_.boundary()[patchi].weights());
            const vectorField& valueInternalCoeffs =
                tValueInternalCoeffs();

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                barYWork[celli] += cmptMultiply
                (
                    valueInternalCoeffs[facei],
                    barYBoundary[patchi][facei]
                );
            }
        }
    }

    const vectorField residual(UEqn.residual());
    const scalarField& V = mesh_.V();

    scalar relPlusDivV = Zero;
    scalar relMinusDivV = Zero;
    scalar relPlus = Zero;
    scalar relMinus = Zero;
    scalar scale = VSMALL;

    forAll(barYWork, celli)
    {
        const vector& Ui = U[celli];
        const vector& Yi = Y[celli];
        const vector& Ri = residual[celli];
        const scalar ri = rAU[celli];
        const scalar Vi = V[celli];

        const vector cPlusDivV = Ui + ri*Ri/Vi;
        const vector cMinusDivV = Ui - ri*Ri/Vi;
        const vector cPlus = Ui + ri*Ri;
        const vector cMinus = Ui - ri*Ri;

        relPlusDivV += mag(Yi - cPlusDivV);
        relMinusDivV += mag(Yi - cMinusDivV);
        relPlus += mag(Yi - cPlus);
        relMinus += mag(Yi - cMinus);
        scale += mag(Yi);
    }

    const scalar predictorIdentityRel =
        min(min(relPlusDivV, relMinusDivV), min(relPlus, relMinus))/scale;

    static bool printedPredictorIdentity = false;
    if (!printedPredictorIdentity)
    {
        Info<< "ATC-T predictor residual identity: "
            << "rel(U+rAU*residual/V)=" << relPlusDivV/scale
            << ", rel(U-rAU*residual/V)=" << relMinusDivV/scale
            << ", rel(U+rAU*residual)=" << relPlus/scale
            << ", rel(U-rAU*residual)=" << relMinus/scale
            << endl;
        printedPredictorIdentity = true;
    }

    if (predictorIdentityRel > scalar(1e-10))
    {
        FatalErrorInFunction
            << "Cannot identify OpenFOAM predictor residual convention for "
            << "Y = rAU*UEqn.H(). Best relative error = "
            << predictorIdentityRel
            << exit(FatalError);
    }

    PredictorReverseSeed result(mesh_);
    vectorField& barU = result.barUold;
    scalarField& barPhiI = result.barPhiInternal;

    vectorField barE(mesh_.nCells(), vector::zero);
    scalarField barD(mesh_.nCells(), Zero);

    forAll(barE, celli)
    {
        barE[celli] = rAU[celli]*barYWork[celli]/V[celli];
        barD[celli] =
            (barYWork[celli] & (U[celli] - Y[celli]))
           *rAU[celli]/V[celli];
        barU[celli] += barYWork[celli];
    }

    if (extraBarDfinal.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "extraBarDfinal size mismatch: "
            << extraBarDfinal.size() << " vs " << mesh_.nCells()
            << exit(FatalError);
    }
    if (extraBarUEqnH1Coeff.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "extraBarUEqnH1Coeff size mismatch: "
            << extraBarUEqnH1Coeff.size() << " vs " << mesh_.nCells()
            << exit(FatalError);
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    forAll(barU, celli)
    {
        barU[celli] -= diagRaw[celli]*barE[celli];
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];

        barU[P] -= lowerRaw[facei]*barE[N];
        barU[N] -= upperRaw[facei]*barE[P];
    }

    forAll(internalCoeffsRaw, patchi)
    {
        const vectorField& coeffs = internalCoeffsRaw[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(coeffs, facei)
        {
            const label celli = faceCells[facei];
            const vector& coeff = coeffs[facei];

            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                barU[celli][cmpt] -= coeff[cmpt]*barE[celli][cmpt];
            }
        }
    }

    const vectorField barUExplicitStokes =
        reverseExplicitStokesViscousSource(barE);
    barU += barUExplicitStokes;

    const scalar chi = boundedConvection_ ? scalar(1) : scalar(0);
    scalar maxMagPhi = gMax(mag(phi.primitiveField()));
    forAll(phi.boundaryField(), patchi)
    {
        maxMagPhi = max(maxMagPhi, gMax(mag(phi.boundaryField()[patchi])));
    }
    const scalar nearZeroPhi = scalar(1e-12)*max(maxMagPhi, scalar(1));
    label nNearZeroPhiInternal = 0;
    label nNearZeroPhiBoundary = 0;

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar phif = phi[facei];
        const vector& Uf = (phif >= scalar(0)) ? U[P] : U[N];

        if (mag(phif) <= nearZeroPhi)
        {
            ++nNearZeroPhiInternal;
        }

        barPhiI[facei] +=
            (Uf & (barE[N] - barE[P]))
          + chi*((barE[P] & U[P]) - (barE[N] & U[N]));
    }

    scalar alphaU = scalar(1);
    word URelaxName = U.select(mesh_.data().isFinalIteration());
    scalar relaxCoeff = scalar(0);
    if (mesh_.relaxEquation(URelaxName, relaxCoeff))
    {
        alphaU = relaxCoeff;
    }

    // The downstream SIMPLE coefficient seed is on UEqn.D()==UEqn.A()*V,
    // not on UEqn.diag().  Route it through the same D-scalar face
    // transpose as the intrinsic predictor diagonal seed.
    const scalarField barDFromY(barD);
    barD += extraBarDfinal;

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar phif = phi[facei];
        scalar dDPhiP = Zero;
        scalar dDPhiN = Zero;

        if (phif >= scalar(0))
        {
            dDPhiP = scalar(1) - chi;
            dDPhiN = chi;
        }
        else
        {
            dDPhiP = -chi;
            dDPhiN = -scalar(1) + chi;
        }

        barPhiI[facei] += (barD[P]*dDPhiP + barD[N]*dDPhiN)/alphaU;

        // fvMatrix::H1() is lduMatrix::H1()/V for non-coupled patches.
        // For Gauss upwind, lower=-phi and upper=0 on P->N outflow,
        // while upper=phi and lower=0 on inflow.
        if (phif >= scalar(0))
        {
            barPhiI[facei] += extraBarUEqnH1Coeff[N]/V[N];
        }
        else
        {
            barPhiI[facei] -= extraBarUEqnH1Coeff[P]/V[P];
        }
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        const fvsPatchScalarField& phip = phi.boundaryField()[patchi];
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();
        scalarField& barPhip = result.barPhiBoundary[patchi];

        forAll(phip, facei)
        {
            const label celli = faceCells[facei];
            const scalar phif = phip[facei];
            const vector& Uf = (phif >= scalar(0)) ? U[celli] : Up[facei];

            if (mag(phif) <= nearZeroPhi)
            {
                ++nNearZeroPhiBoundary;
            }

            barPhip[facei] += barE[celli] & (-Uf + chi*U[celli]);

            const scalar dDPhi =
                (phif >= scalar(0)) ? (scalar(1) - chi) : -chi;
            barPhip[facei] += barD[celli]*dDPhi/alphaU;
        }
    }

    static bool printedPredictorPhiCounts = false;
    if (!printedPredictorPhiCounts)
    {
        reduce(nNearZeroPhiInternal, sumOp<label>());
        reduce(nNearZeroPhiBoundary, sumOp<label>());

        Info<< "ATC-T predictor momentum-flux branch counts: "
            << "nearZeroPhi = " << nearZeroPhi
            << ", nNearZeroPhiInternal = " << nNearZeroPhiInternal
            << ", nNearZeroPhiBoundary = " << nNearZeroPhiBoundary
            << endl;
        printedPredictorPhiCounts = true;
    }

    if (doRelax)
    {
        scalar relaxDiff = Zero;
        scalar relaxScale = VSMALL;

        forAll(barU, celli)
        {
            const vector relaxSourceExpected =
                (diagRelax[celli] - diagRaw[celli])*U[celli];
            const vector relaxSourceActual =
                sourceRelax[celli] - sourceRaw[celli];

            relaxDiff += mag(relaxSourceActual - relaxSourceExpected);
            relaxScale +=
                max(mag(relaxSourceActual), mag(relaxSourceExpected));
        }

        scalar constrainSourceMag = Zero;
        scalar constrainDiagMag = Zero;
        if (doFvOptionsConstrain)
        {
            constrainSourceMag = gSum(mag(sourceFinal - sourceRelax));
            constrainDiagMag = gSum(mag(diagFinal - diagRelax));
        }

        Info<< "ATC-T predictor relaxation-source consistency: rel = "
            << relaxDiff/relaxScale
            << ", sourceConstrainL1 = " << constrainSourceMag
            << ", diagConstrainL1 = " << constrainDiagMag
            << endl;
    }

    forAll(barU, celli)
    {
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            if (component(validVectorComponents, cmpt) == -1)
            {
                barU[celli][cmpt] = Zero;
            }
        }
    }

    return result;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::thermalAdjointSimple::reverseAdjustPhiExact
(
    const surfaceScalarField& phiBeforeAdjust,
    const surfaceScalarField& barPhiAfterAdjust,
    const volVectorField& U
)
{
    auto tBarPhiBefore = tmp<surfaceScalarField>::New
    (
        IOobject
        (
            "ATCTReverseAdjustPhiSeed",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        barPhiAfterAdjust
    );
    surfaceScalarField& barPhiBefore = tBarPhiBefore.ref();

    const volScalarField& p = primalVars_.p();

    static bool printedAdjustPhiReverseSource = false;
    if (!printedAdjustPhiReverseSource)
    {
        Info<< "ATC-T reverseAdjustPhiExact uses "
            << "$FOAM_SRC/finiteVolume/cfdTools/general/adjustPhi/"
            << "adjustPhi.C::adjustPhi. p.needReference() = "
            << p.needReference()
            << ". If false, OpenFOAM returns without changing phi."
            << endl;
        printedAdjustPhiReverseSource = true;
    }

    if (!p.needReference())
    {
        return tBarPhiBefore;
    }

    scalar massIn = Zero;
    scalar fixedMassOut = Zero;
    scalar adjustableMassOut = Zero;

    forAll(phiBeforeAdjust.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvsPatchScalarField& phip =
            phiBeforeAdjust.boundaryField()[patchi];

        if (phip.coupled())
        {
            continue;
        }

        const bool fixedPatch =
            Up.fixesValue() && !isA<inletOutletFvPatchVectorField>(Up);

        forAll(phip, facei)
        {
            if (phip[facei] < scalar(0))
            {
                massIn -= phip[facei];
            }
            else if (fixedPatch)
            {
                fixedMassOut += phip[facei];
            }
            else
            {
                adjustableMassOut += phip[facei];
            }
        }
    }

    const scalar totalFlux = VSMALL + sum(mag(phiBeforeAdjust)).value();
    reduce(massIn, sumOp<scalar>());
    reduce(fixedMassOut, sumOp<scalar>());
    reduce(adjustableMassOut, sumOp<scalar>());

    const scalar magAdjustableMassOut = mag(adjustableMassOut);
    if
    (
        magAdjustableMassOut <= VSMALL
     || magAdjustableMassOut/totalFlux <= SMALL
    )
    {
        return tBarPhiBefore;
    }

    const scalar massCorr =
        (massIn - fixedMassOut)/adjustableMassOut;

    scalar barMassCorr = Zero;
    forAll(phiBeforeAdjust.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvsPatchScalarField& phiPatch =
            phiBeforeAdjust.boundaryField()[patchi];
        const fvsPatchScalarField& seedPatch =
            barPhiAfterAdjust.boundaryField()[patchi];

        if (phiPatch.coupled())
        {
            continue;
        }

        const bool adjustablePatch =
            !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

        if (!adjustablePatch)
        {
            continue;
        }

        forAll(phiPatch, facei)
        {
            if (phiPatch[facei] > scalar(0))
            {
                barMassCorr += seedPatch[facei]*phiPatch[facei];
            }
        }
    }
    reduce(barMassCorr, sumOp<scalar>());

    forAll(phiBeforeAdjust.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvsPatchScalarField& phiPatch =
            phiBeforeAdjust.boundaryField()[patchi];
        fvsPatchScalarField& seedBeforePatch =
            barPhiBefore.boundaryFieldRef()[patchi];

        if (phiPatch.coupled())
        {
            continue;
        }

        const bool fixedPatch =
            Up.fixesValue() && !isA<inletOutletFvPatchVectorField>(Up);
        const bool adjustablePatch =
            !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

        forAll(phiPatch, facei)
        {
            if (phiPatch[facei] < scalar(0))
            {
                seedBeforePatch[facei] -=
                    barMassCorr/adjustableMassOut;
            }
            else if (fixedPatch)
            {
                seedBeforePatch[facei] -=
                    barMassCorr/adjustableMassOut;
            }
            else if (adjustablePatch)
            {
                seedBeforePatch[facei] =
                    massCorr*seedBeforePatch[facei]
                  - massCorr*barMassCorr/adjustableMassOut;
            }
        }
    }

    return tBarPhiBefore;
}


Foam::thermalAdjointSimple::SimpleMapSeed
Foam::thermalAdjointSimple::reverseOneSimpleMapSeedImpl
(
    const SimpleMapSeed& seedNew,
    const volScalarField& rAtUBase,
    SmoothPhiCoefficientTrace* trace
)
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceScalarField& w = mesh_.weights();
    const surfaceVectorField& Sf = mesh_.Sf();
    const scalarField& V = mesh_.V();
    const volScalarField& p = primalVars_.p();
    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();
    UEqn.relax();
    fvOptions.constrain(UEqn);

    volScalarField rAU(1.0/UEqn.A());
    tmp<volScalarField> trAtU(rAU);

    volScalarField pForFinalCorrection
    (
        IOobject
        (
            "ATCTReverseMapPForFinalCorrection",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );

    volVectorField HbyABeforeAdjust
    (
        constrainHbyA(rAU*UEqn.H(), U, pForFinalCorrection)
    );
    surfaceScalarField phiBeforeAdjust
    (
        IOobject
        (
            "ATCTReverseMapPhiBeforeAdjust",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(HbyABeforeAdjust)
    );

    surfaceScalarField phiForPressure
    (
        IOobject
        (
            "ATCTReverseMapPhiForPressure",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiBeforeAdjust
    );
    adjustPhi(phiForPressure, U, pForFinalCorrection);

    if (solverControl_().consistent())
    {
        trAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiForPressure +=
            fvc::interpolate(trAtU() - rAU)
           *fvc::snGrad(pForFinalCorrection)
           *mesh_.magSf();
    }

    const volScalarField& rAtUMap = trAtU();

    constrainPressure(pForFinalCorrection, U, phiForPressure, rAtUMap);

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtUMap, pForFinalCorrection) == fvc::div(phiForPressure)
    );
    pEqn.setReference
    (
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );

    vectorField source(mesh_.nCells(), vector::zero);
    forAll(source, celli)
    {
        source[celli] = seedNew.barU[celli]/V[celli];
    }

    scalarField barrAtU(mesh_.nCells(), Zero);
    scalarField barrAU(mesh_.nCells(), Zero);
    scalarField barq(mesh_.nCells(), Zero);
    scalarField barqFromH1(mesh_.nCells(), Zero);
    scalarField barqFromF2InternalFaces(mesh_.nCells(), Zero);
    scalarField barqFromF2BoundaryCurrent(mesh_.nCells(), Zero);
    scalarField barQFaceInternal(mesh_.nInternalFaces(), Zero);
    List<scalarField> barQFaceBoundary(mesh_.boundary().size());
    forAll(barQFaceBoundary, patchi)
    {
        barQFaceBoundary[patchi].setSize
        (
            mesh_.boundary()[patchi].size(),
            Zero
        );
    }
    scalarField barUEqnH1Coeff(mesh_.nCells(), Zero);
    scalarField barDfinal(mesh_.nCells(), Zero);

    tmp<volVectorField> tGradPFinal = fvc::grad(pForFinalCorrection);
    const volVectorField& gradPFinal = tGradPFinal();
    forAll(barrAtU, celli)
    {
        barrAtU[celli] -= seedNew.barU[celli] & gradPFinal[celli];
    }
    if (trace)
    {
        trace->barrAtUFromFinalCorrection = barrAtU;
    }

    scalarField pSeed(mesh_.nCells(), Zero);

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar SfSource =
            Sf[facei] & (rAtUMap[P]*source[P] - rAtUMap[N]*source[N]);

        pSeed[P] -= w[facei]*SfSource;
        pSeed[N] -= (scalar(1) - w[facei])*SfSource;
    }

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        const fvPatch& patch = mesh_.boundary()[patchi];
        const labelUList& faceCells = patch.faceCells();
        const fvsPatchVectorField& SfPatch = Sf.boundaryField()[patchi];

        tmp<scalarField> tValueInternalCoeffs =
            pForFinalCorrection.boundaryField()[patchi].valueInternalCoeffs
            (
                patch.weights()
            );
        const scalarField& valueInternalCoeffs = tValueInternalCoeffs();

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            pSeed[celli] -=
                valueInternalCoeffs[facei]
               *(SfPatch[facei] & (rAtUMap[celli]*source[celli]));
        }
    }

    pSeed += seedNew.barp;
    if (trace)
    {
        trace->barpNew = pSeed;
    }

    scalar pRelaxCoeff = scalar(1);
    word pRelaxName = p.name();
    if (p.mesh().data().isFinalIteration())
    {
        pRelaxName += "Final";
    }
    p.mesh().relaxField(pRelaxName, pRelaxCoeff);

    scalarField pSolveSeed(pSeed);
    pSolveSeed *= pRelaxCoeff;

    const scalarField& pLower = pEqn.lower();
    const scalarField& pUpper = pEqn.upper();
    const FieldField<Field, scalar>& pInternalCoeffs =
        pEqn.internalCoeffs();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar barPhi = seedNew.barPhiInternal[facei];

        pSolveSeed[P] += pLower[facei]*barPhi;
        pSolveSeed[N] -= pUpper[facei]*barPhi;
    }

    forAll(pInternalCoeffs, patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        const scalarField& intCoeffs = pInternalCoeffs[patchi];
        const scalarField& barPhip = seedNew.barPhiBoundary[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(intCoeffs, facei)
        {
            pSolveSeed[faceCells[facei]] -=
                intCoeffs[facei]*barPhip[facei];
        }
    }
    if (trace)
    {
        trace->barpSolve = pSolveSeed;
    }

    scalarField pOldSeed(pSeed);
    pOldSeed *= scalar(1) - pRelaxCoeff;

    volScalarField mu
    (
        IOobject
        (
            "paATCTReverseSimpleMapSeed",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    mu = dimensionedScalar(mu.dimensions(), Zero);

    fvScalarMatrix muEqn(fvm::laplacian(rAtUMap, mu));
    muEqn.source() = pSolveSeed;
    muEqn.setReference(solverControl_().pRefCell(), Zero);
    dictionary muSolver(muEqn.solverDict("p"));
    muSolver.set("relTol", scalar(0));
    muSolver.set("tolerance", scalar(1e-12));
    muEqn.solve(muSolver);

    surfaceScalarField predictorFSeedAfterAdjust
    (
        IOobject
        (
            "ATCTReverseMapFSeedAfterAdjust",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiBeforeAdjust
    );
    predictorFSeedAfterAdjust =
        dimensionedScalar(predictorFSeedAfterAdjust.dimensions(), Zero);

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        predictorFSeedAfterAdjust.primitiveFieldRef()[facei] =
            seedNew.barPhiInternal[facei] + mu[P] - mu[N];
    }

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        fvsPatchScalarField& seedp =
            predictorFSeedAfterAdjust.boundaryFieldRef()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(seedp, facei)
        {
            seedp[facei] =
                seedNew.barPhiBoundary[patchi][facei]
              + mu[faceCells[facei]];
        }
    }
    if (trace)
    {
        trace->barF2Internal =
            predictorFSeedAfterAdjust.primitiveField();
        trace->barF1Internal =
            predictorFSeedAfterAdjust.primitiveField();
        forAll(trace->barF2Boundary, patchi)
        {
            trace->barF2Boundary[patchi] =
                predictorFSeedAfterAdjust.boundaryField()[patchi];
            trace->barF1Boundary[patchi] =
                predictorFSeedAfterAdjust.boundaryField()[patchi];
        }
    }

    volScalarField pSolveForPressureCoeff
    (
        IOobject
        (
            "ATCTReverseMapPressureCoeffPSolve",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    pSolveForPressureCoeff.storePrevIter();
    constrainPressure(pSolveForPressureCoeff, U, phiForPressure, rAtUMap);
    mesh_.setFluxRequired(pSolveForPressureCoeff.name());

    fvScalarMatrix pCoeffEqn
    (
        fvm::laplacian(rAtUMap, pSolveForPressureCoeff)
     == fvc::div(phiForPressure)
    );
    pCoeffEqn.setReference
    (
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );
    dictionary pCoeffSolver(pCoeffEqn.solverDict("p"));
    pCoeffSolver.set("relTol", scalar(0));
    pCoeffSolver.set("tolerance", scalar(1e-12));
    pCoeffEqn.solve(pCoeffSolver);

    scalarField pressureBarrAtU =
        reversePressureRAtUExact
        (
            rAtUMap,
            pSolveForPressureCoeff,
            predictorFSeedAfterAdjust
        );
    barrAtU += pressureBarrAtU;
    if (trace)
    {
        trace->barrAtUFromPressureStage = pressureBarrAtU;
    }

    tmp<surfaceScalarField> tPredictorFSeedBeforeAdjust =
        reverseAdjustPhiExact
        (
            phiBeforeAdjust,
            predictorFSeedAfterAdjust,
            U
        );
    const surfaceScalarField& predictorFSeedBeforeAdjust =
        tPredictorFSeedBeforeAdjust();
    if (trace)
    {
        trace->barF0Internal =
            predictorFSeedBeforeAdjust.primitiveField();
        forAll(trace->barF0Boundary, patchi)
        {
            trace->barF0Boundary[patchi] =
                predictorFSeedBeforeAdjust.boundaryField()[patchi];
        }
    }

    VolVectorValueSeed barH0Direct(mesh_);
    barH0Direct.internal = seedNew.barU;

    VolVectorValueSeed barH0FromFlux(mesh_);

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const vector FSeedSf =
            predictorFSeedBeforeAdjust[facei]*Sf[facei];

        barH0FromFlux.internal[P] += w[facei]*FSeedSf;
        barH0FromFlux.internal[N] += (scalar(1) - w[facei])*FSeedSf;
    }

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        const fvsPatchScalarField& seedPatch =
            predictorFSeedBeforeAdjust.boundaryField()[patchi];
        const fvsPatchVectorField& SfPatch = Sf.boundaryField()[patchi];

        forAll(seedPatch, facei)
        {
            barH0FromFlux.boundary[patchi][facei] +=
                seedPatch[facei]*SfPatch[facei];
        }
    }

    VolVectorValueSeed barH0(mesh_);
    barH0.internal = barH0Direct.internal;
    barH0.internal += barH0FromFlux.internal;
    forAll(barH0.boundary, patchi)
    {
        barH0.boundary[patchi] = barH0Direct.boundary[patchi];
        barH0.boundary[patchi] += barH0FromFlux.boundary[patchi];
    }

    if (trace)
    {
        trace->barH0Internal = barH0.internal;
        trace->barH0DirectInternal = barH0Direct.internal;
        trace->barH0FromFluxInternal = barH0FromFlux.internal;
        trace->barH0TotalInternal = barH0.internal;
        forAll(trace->barH0Boundary, patchi)
        {
            trace->barH0Boundary[patchi] = barH0.boundary[patchi];
            trace->barH0DirectBoundary[patchi] =
                barH0Direct.boundary[patchi];
            trace->barH0FromFluxBoundary[patchi] =
                barH0FromFlux.boundary[patchi];
            trace->barH0TotalBoundary[patchi] = barH0.boundary[patchi];
        }
    }

    if (solverControl_().consistent())
    {
        const volScalarField q(rAtUMap - rAU);
        tmp<volVectorField> tGradPOld = fvc::grad(p);
        const volVectorField& gradPOld = tGradPOld();
        tmp<surfaceScalarField> tqf = fvc::interpolate(q);
        const surfaceScalarField& qf = tqf();
        const surfaceScalarField& magSf = mesh_.magSf();
        tmp<fv::snGradScheme<scalar>> tSnGradScheme =
            fv::snGradScheme<scalar>::New
            (
                mesh_,
                mesh_.snGradScheme("snGrad(" + p.name() + ')')
            );
        tmp<surfaceScalarField> tSnGradDeltaCoeffs =
            tSnGradScheme().deltaCoeffs(p);
        const surfaceScalarField& snGradDeltaCoeffs =
            tSnGradDeltaCoeffs();

        if (trace)
        {
            trace->barrAtUBeforeQ = barrAtU;
            trace->barrAUBeforeQ = barrAU;
        }

        forAll(barq, celli)
        {
            const scalar h1Seed = seedNew.barU[celli] & gradPOld[celli];
            barq[celli] += h1Seed;
            barqFromH1[celli] += h1Seed;
        }

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];
            const scalar eta =
                Sf[facei] & (q[P]*source[P] - q[N]*source[N]);

            pOldSeed[P] += w[facei]*eta;
            pOldSeed[N] += (scalar(1) - w[facei])*eta;

            const scalar nCoeff =
                qf[facei]*magSf[facei]*snGradDeltaCoeffs[facei];
            const scalar FSeed = predictorFSeedAfterAdjust[facei];
            const scalar qFaceSeed =
                FSeed*magSf[facei]*snGradDeltaCoeffs[facei]
               *(p[N] - p[P]);
            barQFaceInternal[facei] = qFaceSeed;

            pOldSeed[P] -= nCoeff*FSeed;
            pOldSeed[N] += nCoeff*FSeed;

            const scalar ownerSeed = w[facei]*qFaceSeed;
            const scalar neighbourSeed = (scalar(1) - w[facei])*qFaceSeed;
            barq[P] += ownerSeed;
            barq[N] += neighbourSeed;
            barqFromF2InternalFaces[P] += ownerSeed;
            barqFromF2InternalFaces[N] += neighbourSeed;
        }

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();
            const fvsPatchScalarField& qfp =
                qf.boundaryField()[patchi];
            const fvsPatchScalarField& magSfp =
                magSf.boundaryField()[patchi];
            const fvsPatchScalarField& FSeedp =
                predictorFSeedAfterAdjust.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp =
                Sf.boundaryField()[patchi];

            tmp<scalarField> tpValueInternalCoeffs =
                p.boundaryField()[patchi].valueInternalCoeffs
                (
                    mesh_.boundary()[patchi].weights()
                );
            const scalarField& pValueInternalCoeffs =
                tpValueInternalCoeffs();

            tmp<scalarField> tpGradientInternalCoeffs =
                p.boundaryField()[patchi].gradientInternalCoeffs();
            const scalarField& pGradientInternalCoeffs =
                tpGradientInternalCoeffs();
            tmp<scalarField> tpGradientBoundaryCoeffs =
                p.boundaryField()[patchi].gradientBoundaryCoeffs();
            const scalarField& pGradientBoundaryCoeffs =
                tpGradientBoundaryCoeffs();

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];

                pOldSeed[celli] +=
                    pValueInternalCoeffs[facei]
                   *(Sfp[facei] & (q[celli]*source[celli]));

                pOldSeed[faceCells[facei]] +=
                    FSeedp[facei]
                   *qfp[facei]
                   *magSfp[facei]
                   *pGradientInternalCoeffs[facei];

                const scalar qFaceSeed =
                    FSeedp[facei]
                   *magSfp[facei]
                   *(
                        pGradientInternalCoeffs[facei]*p[celli]
                      + pGradientBoundaryCoeffs[facei]
                    );

                barQFaceBoundary[patchi][facei] = qFaceSeed;
                barq[celli] += qFaceSeed;
                barqFromF2BoundaryCurrent[celli] += qFaceSeed;
            }
        }

        barrAtU += barq;
        barrAU -= barq;

        if (trace)
        {
            trace->barq = barq;
            trace->barrAtUAfterQ = barrAtU;
            trace->barrAUAfterQ = barrAU;
            trace->barrAUBeforeRAtUReverse = barrAU;
            trace->barqFromH1 = barqFromH1;
            trace->barqFromF2InternalFaces = barqFromF2InternalFaces;
            trace->barqFromF2BoundaryCurrent = barqFromF2BoundaryCurrent;
            trace->barQFaceInternal = barQFaceInternal;
            forAll(trace->barQFaceBoundary, patchi)
            {
                trace->barQFaceBoundary[patchi] =
                    barQFaceBoundary[patchi];
            }
            trace->barrAtUFromQ = barq;
            trace->barrAUFromQ = -barq;
        }
    }

    forAll(barrAtU, celli)
    {
        barrAU[celli] +=
            barrAtU[celli]*sqr(rAtUMap[celli])/sqr(rAU[celli]);
        barUEqnH1Coeff[celli] += barrAtU[celli]*sqr(rAtUMap[celli]);
    }
    if (trace)
    {
        trace->rAUUsed = rAU.primitiveField();
        trace->rAtUUsed = rAtUMap.primitiveField();
        trace->barrAUAfterRAtUReverse = barrAU;
        trace->barH1Expected = barUEqnH1Coeff;
    }

    forAll(barDfinal, celli)
    {
        barDfinal[celli] +=
            -barrAU[celli]*sqr(rAU[celli])/V[celli];
    }
    if (trace)
    {
        trace->barDExpected = barDfinal;
        trace->barDfinal = barDfinal;
        trace->barDDownstream = barDfinal;
        trace->barDCombined = barDfinal;
        trace->barUEqnH1Coeff = barUEqnH1Coeff;
    }

    ConstrainHbyAReverseSeed hbyASeed(mesh_);
    hbyASeed.barYInternal = barH0.internal;

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        const fvPatchVectorField& Up =
            primalVars_.U().boundaryField()[patchi];
        const fvPatchScalarField& pp =
            primalVars_.p().boundaryField()[patchi];
        const bool overwritten =
            !Up.assignable()
         && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>(pp);

        if (overwritten)
        {
            hbyASeed.barUBoundary[patchi] += barH0.boundary[patchi];
        }
        else
        {
            hbyASeed.barYBoundary[patchi] += barH0.boundary[patchi];
        }
    }
    if (trace)
    {
        trace->barYInternal = hbyASeed.barYInternal;
        forAll(trace->barYBoundary, patchi)
        {
            trace->barYBoundary[patchi] =
                hbyASeed.barYBoundary[patchi];
            trace->barUBoundary[patchi] =
                hbyASeed.barUBoundary[patchi];
        }
    }

    vectorField barUFromBoundary(mesh_.nCells(), vector::zero);
    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        const fvPatchVectorField& Up =
            primalVars_.U().boundaryField()[patchi];
        tmp<vectorField> tValueInternalCoeffs =
            Up.valueInternalCoeffs(mesh_.boundary()[patchi].weights());
        const vectorField& valueInternalCoeffs = tValueInternalCoeffs();
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(faceCells, facei)
        {
            barUFromBoundary[faceCells[facei]] += cmptMultiply
            (
                valueInternalCoeffs[facei],
                hbyASeed.barUBoundary[patchi][facei]
            );
        }
    }
    if (trace)
    {
        trace->barUFromBoundary = barUFromBoundary;
    }

    PredictorReverseSeed predictorOld =
        reversePredictorYFullState
        (
            hbyASeed.barYInternal,
            hbyASeed.barYBoundary,
            barDfinal,
            barUEqnH1Coeff,
            true,
            true
        );

    SimpleMapSeed seedOld(mesh_);
    seedOld.barU = predictorOld.barUold;
    seedOld.barU += barUFromBoundary;
    seedOld.barp = pOldSeed;
    seedOld.barPhiInternal = predictorOld.barPhiInternal;
    forAll(seedOld.barPhiBoundary, patchi)
    {
        seedOld.barPhiBoundary[patchi] = predictorOld.barPhiBoundary[patchi];
    }

    return seedOld;
}


Foam::thermalAdjointSimple::SimpleMapSeed
Foam::thermalAdjointSimple::reverseOneSimpleMapSeed
(
    const SimpleMapSeed& seedNew,
    const volScalarField& rAtUBase
)
{
    return reverseOneSimpleMapSeedImpl(seedNew, rAtUBase, nullptr);
}


Foam::tmp<Foam::vectorField>
Foam::thermalAdjointSimple::reverseSimpleMapVelocitySeed
(
    const vectorField& barUnewInt,
    const volScalarField& rAtU
)
{
    SimpleMapSeed seedNew(mesh_);
    seedNew.barU = barUnewInt;

    SimpleMapSeed seedOld = reverseOneSimpleMapSeed(seedNew, rAtU);

    auto tBarU = tmp<vectorField>::New(mesh_.nCells(), vector::zero);
    tBarU.ref() = seedOld.barU;

    return tBarU;
}


void Foam::thermalAdjointSimple::checkProjectedFluxTranspose
(
    const surfaceScalarField& gPhi,
    const surfaceScalarField& gF,
    const volVectorField& source,
    const volScalarField& rAU,
    const volScalarField& rAtU,
    const fvScalarMatrix& pEqn
)
{
    if (!checkProjectionTranspose_)
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    surfaceScalarField dF
    (
        IOobject
        (
            "ATCTdF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(primalVars_.phi().dimensions(), Zero)
    );

    forAll(dF.primitiveFieldRef(), facei)
    {
        const label k = (37*facei + 17) % 101;
        dF.primitiveFieldRef()[facei] = scalar(k)/scalar(101) - scalar(0.5);
    }

    volScalarField dp
    (
        IOobject
        (
            "paATCTDelta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.p()
    );
    dp = dimensionedScalar(dp.dimensions(), Zero);

    fvScalarMatrix dpEqn(fvm::laplacian(rAtU, dp) == fvc::div(dF));
    dpEqn.setReference(solverControl_().pRefCell(), Zero);
    dictionary projectionSolver(dpEqn.solverDict("p"));
    projectionSolver.set("relTol", scalar(0));
    projectionSolver.set("tolerance", scalar(1e-12));
    dpEqn.solve(projectionSolver);

    const scalarField& lower = pEqn.lower();
    const scalarField& upper = pEqn.upper();

    scalar lhs = Zero;
    scalar rhs = Zero;

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];

        const scalar Cdp = upper[facei]*dp[N] - lower[facei]*dp[P];
        lhs += gPhi[facei]*(dF[facei] - Cdp);
        rhs += gF[facei]*dF[facei];
    }

    const scalar projectionScale = max(max(mag(lhs), mag(rhs)), VSMALL);

    Info<< "ATC-T projection transpose check: lhs = " << lhs
        << ", rhs = " << rhs
        << ", rel = " << mag(lhs - rhs)/projectionScale << endl;

    volVectorField dH
    (
        IOobject
        (
            "ATCTdH",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(primalVars_.U().dimensions(), Zero)
    );

    vectorField& dHI = dH.primitiveFieldRef();
    forAll(dHI, celli)
    {
        const scalar x = scalar((19*celli + 3) % 97)/scalar(97) - scalar(0.5);
        const scalar y = scalar((23*celli + 5) % 89)/scalar(89) - scalar(0.5);
        const scalar z = scalar((29*celli + 7) % 83)/scalar(83) - scalar(0.5);
        dHI[celli] = vector(x, y, z);
    }

    const surfaceScalarField& w = mesh_.weights();
    const surfaceVectorField& Sf = mesh_.Sf();
    const scalarField& V = mesh_.V();

    scalar lhsH = Zero;
    scalar rhsH = Zero;
    vectorField HbyAGrad(mesh_.nCells(), Zero);

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];

        const vector dHf = w[facei]*dH[P] + (scalar(1) - w[facei])*dH[N];
        const vector Gf = gF[facei]*Sf[facei];
        lhsH += gF[facei]*(Sf[facei] & dHf);
        HbyAGrad[P] += w[facei]*Gf;
        HbyAGrad[N] += (scalar(1) - w[facei])*Gf;
    }

    forAll(dHI, celli)
    {
        rhsH += HbyAGrad[celli] & dH[celli];
    }

    const scalar sourceScale = max(max(mag(lhsH), mag(rhsH)), VSMALL);

    Info<< "ATC-T internal HbyA transpose check: lhs = " << lhsH
        << ", rhs = " << rhsH
        << ", rel = " << mag(lhsH - rhsH)/sourceScale << endl;

    volScalarField dpProbe
    (
        IOobject
        (
            "paATCTProbe",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.p()
    );
    dpProbe = dimensionedScalar(dpProbe.dimensions(), Zero);

    scalarField& dpProbeI = dpProbe.primitiveFieldRef();
    forAll(dpProbeI, celli)
    {
        dpProbeI[celli] =
            scalar((31*celli + 11) % 103)/scalar(103) - scalar(0.5);
    }
    dpProbe.correctBoundaryConditions();

    mesh_.setFluxRequired(dpProbe.name());
    fvScalarMatrix dpProbeEqn(fvm::laplacian(rAtU, dpProbe));
    surfaceScalarField CdpFull(dpProbeEqn.flux());

    volVectorField dU
    (
        IOobject
        (
            "ATCTdU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.U()
    );
    dU = dimensionedVector(dU.dimensions(), Zero);

    vectorField& dUI = dU.primitiveFieldRef();
    forAll(dUI, celli)
    {
        const scalar x = scalar((41*celli + 13) % 107)/scalar(107) - scalar(0.5);
        const scalar y = scalar((43*celli + 17) % 109)/scalar(109) - scalar(0.5);
        const scalar z = scalar((47*celli + 19) % 113)/scalar(113) - scalar(0.5);
        dUI[celli] = vector(x, y, z);
    }
    dU.correctBoundaryConditions();

    volVectorField::Boundary& dUb = dU.boundaryFieldRef();

    forAll(dUb, patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() != "empty"
         && !primalVars_.U().boundaryField()[patchi].assignable()
        )
        {
            dUb[patchi] == vector::zero;
        }
    }

    surfaceScalarField BdUFull(fvc::flux(dU));

    scalar pressureNullFull = Zero;
    scalar pressureNullFullScale = VSMALL;
    scalar lhsUFull = Zero;
    scalar rhsU = Zero;

    forAll(own, facei)
    {
        pressureNullFull += gF[facei]*CdpFull[facei];
        pressureNullFullScale += mag(gF[facei]*CdpFull[facei]);
        lhsUFull += gF[facei]*(BdUFull[facei] + CdpFull[facei]);
    }

    forAll(dUI, celli)
    {
        rhsU += V[celli]*(source[celli] & dU[celli]);
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        const fvsPatchScalarField& gFpatch = gF.boundaryField()[patchi];
        const fvsPatchScalarField& CdpPatch =
            CdpFull.boundaryField()[patchi];
        const fvsPatchScalarField& BdUPatch =
            BdUFull.boundaryField()[patchi];

        forAll(gFpatch, facei)
        {
            pressureNullFull += gFpatch[facei]*CdpPatch[facei];
            pressureNullFullScale += mag(gFpatch[facei]*CdpPatch[facei]);
            lhsUFull +=
                gFpatch[facei]*(BdUPatch[facei] + CdpPatch[facei]);
        }
    }

    const scalar correctedVelocityFullScale =
        max(max(mag(lhsUFull), mag(rhsU)), VSMALL);

    Info<< "ATC-T full pressure-null check: value = " << pressureNullFull
        << ", rel = " << mag(pressureNullFull)/pressureNullFullScale
        << endl;

    Info<< "ATC-T full corrected-velocity check: lhs = " << lhsUFull
        << ", rhs = " << rhsU
        << ", rel = "
        << mag(lhsUFull - rhsU)/correctedVelocityFullScale << endl;

    volVectorField dHFinal
    (
        IOobject
        (
            "ATCTdHFinalCorrection",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.U()
    );
    dHFinal = dimensionedVector(dHFinal.dimensions(), Zero);

    vectorField& dHFinalI = dHFinal.primitiveFieldRef();
    forAll(dHFinalI, celli)
    {
        const scalar x = scalar((61*celli + 23) % 137)/scalar(137) - scalar(0.5);
        const scalar y = scalar((67*celli + 29) % 139)/scalar(139) - scalar(0.5);
        const scalar z = scalar((71*celli + 31) % 149)/scalar(149) - scalar(0.5);
        dHFinalI[celli] = vector(x, y, z);
    }
    dHFinal.correctBoundaryConditions();

    volScalarField dpFinal
    (
        IOobject
        (
            "ATCTdpFinalCorrection",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.p()
    );
    dpFinal = dimensionedScalar(dpFinal.dimensions(), Zero);

    scalarField& dpFinalI = dpFinal.primitiveFieldRef();
    forAll(dpFinalI, celli)
    {
        dpFinalI[celli] =
            scalar((73*celli + 37) % 151)/scalar(151) - scalar(0.5);
    }
    dpFinal.correctBoundaryConditions();

    const volVectorField gradDpFinal(fvc::grad(dpFinal));

    scalar finalCorrectionLhs = Zero;
    scalar finalCorrectionH = Zero;
    scalar finalCorrectionP = Zero;

    forAll(dHFinalI, celli)
    {
        const vector dUFinal =
            dHFinal[celli] - rAtU[celli]*gradDpFinal[celli];

        finalCorrectionLhs += V[celli]*(source[celli] & dUFinal);
        finalCorrectionH += V[celli]*(source[celli] & dHFinal[celli]);
        finalCorrectionP -=
            V[celli]*(source[celli] & (rAtU[celli]*gradDpFinal[celli]));
    }

    const scalar finalCorrectionScale =
        max
        (
            mag(finalCorrectionLhs)
          + mag(finalCorrectionH)
          + mag(finalCorrectionP),
            VSMALL
        );

    Info<< "ATC-T final velocity-correction transpose seed check: lhs = "
        << finalCorrectionLhs
        << ", Hseed = " << finalCorrectionH
        << ", pSeedTerm = " << finalCorrectionP
        << ", rel = "
        << mag(finalCorrectionLhs - finalCorrectionH - finalCorrectionP)
          /finalCorrectionScale
        << endl;

    scalarField pSeed(mesh_.nCells(), Zero);
    const volScalarField& p = primalVars_.p();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar SfSource =
            mesh_.Sf()[facei] & (rAtU[P]*source[P] - rAtU[N]*source[N]);

        pSeed[P] -= w[facei]*SfSource;
        pSeed[N] -= (scalar(1) - w[facei])*SfSource;
    }

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        const fvPatch& patch = mesh_.boundary()[patchi];
        const labelUList& faceCells = patch.faceCells();
        const fvsPatchVectorField& SfPatch =
            mesh_.Sf().boundaryField()[patchi];

        tmp<scalarField> tValueInternalCoeffs =
            p.boundaryField()[patchi].valueInternalCoeffs(patch.weights());
        const scalarField& valueInternalCoeffs = tValueInternalCoeffs();

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            pSeed[celli] -=
                valueInternalCoeffs[facei]
               *(SfPatch[facei] & (rAtU[celli]*source[celli]));
        }
    }

    volScalarField mu
    (
        IOobject
        (
            "paATCTFinalCorrectionSeed",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    mu = dimensionedScalar(mu.dimensions(), Zero);

    fvScalarMatrix muEqn(fvm::laplacian(rAtU, mu));
    muEqn.source() = pSeed;
    muEqn.setReference(solverControl_().pRefCell(), Zero);
    dictionary muSolver(muEqn.solverDict("p"));
    muSolver.set("relTol", scalar(0));
    muSolver.set("tolerance", scalar(1e-12));
    muEqn.solve(muSolver);

    surfaceScalarField dF2
    (
        IOobject
        (
            "ATCTdF2FinalCorrection",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(primalVars_.phi().dimensions(), Zero)
    );

    forAll(dF2.primitiveFieldRef(), facei)
    {
        dF2.primitiveFieldRef()[facei] =
            scalar((79*facei + 41) % 157)/scalar(157) - scalar(0.5);
    }

    volScalarField dpFromF2
    (
        IOobject
        (
            "paATCTFinalCorrectionDelta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    dpFromF2 = dimensionedScalar(dpFromF2.dimensions(), Zero);

    fvScalarMatrix dpFromF2Eqn
    (
        fvm::laplacian(rAtU, dpFromF2) == fvc::div(dF2)
    );
    dpFromF2Eqn.setReference(solverControl_().pRefCell(), Zero);
    dictionary dpFromF2Solver(dpFromF2Eqn.solverDict("p"));
    dpFromF2Solver.set("relTol", scalar(0));
    dpFromF2Solver.set("tolerance", scalar(1e-12));
    dpFromF2Eqn.solve(dpFromF2Solver);

    const volVectorField gradDpFromF2(fvc::grad(dpFromF2));

    scalar pressureSeedSolveLhs = Zero;
    scalar pressureSeedSolveRhs = Zero;
    scalar combinedLhs = Zero;
    scalar combinedRhs = Zero;

    forAll(pSeed, celli)
    {
        pressureSeedSolveLhs += pSeed[celli]*dpFromF2[celli];
        combinedLhs +=
            V[celli]
           *(
                source[celli]
              & (
                    dHFinal[celli]
                  - rAtU[celli]*gradDpFromF2[celli]
                )
            );
        combinedRhs += V[celli]*(source[celli] & dHFinal[celli]);
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar FSeed = mu[P] - mu[N];

        pressureSeedSolveRhs += FSeed*dF2[facei];
        combinedRhs += FSeed*dF2[facei];
    }

    const scalar pressureSeedSolveScale =
        max(max(mag(pressureSeedSolveLhs), mag(pressureSeedSolveRhs)), VSMALL);
    const scalar combinedScale =
        max(max(mag(combinedLhs), mag(combinedRhs)), VSMALL);

    Info<< "ATC-T final-correction pressure-solve transpose check: lhs = "
        << pressureSeedSolveLhs
        << ", rhs = " << pressureSeedSolveRhs
        << ", rel = "
        << mag(pressureSeedSolveLhs - pressureSeedSolveRhs)
          /pressureSeedSolveScale
        << endl;

    Info<< "ATC-T final correction plus pressure solve check: lhs = "
        << combinedLhs
        << ", rhs = " << combinedRhs
        << ", rel = "
        << mag(combinedLhs - combinedRhs)/combinedScale
        << endl;

    surfaceScalarField dF2Boundary
    (
        IOobject
        (
            "ATCTdF2BoundaryPressure",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dF2
    );

    label nBoundaryFluxPerturb = 0;

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        fvsPatchScalarField& dF2p =
            dF2Boundary.boundaryFieldRef()[patchi];

        forAll(dF2p, facei)
        {
            dF2p[facei] =
                scalar((197*(facei + 1) + 151*(patchi + 1)) % 281)
               /scalar(281) - scalar(0.5);
            ++nBoundaryFluxPerturb;
        }
    }

    volVectorField zeroU
    (
        IOobject
        (
            "ATCTzeroUBoundaryPressure",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.U()
    );
    zeroU = dimensionedVector(zeroU.dimensions(), vector::zero);

    volScalarField dpBoundary
    (
        IOobject
        (
            "paATCTBoundaryPressureDelta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    dpBoundary = dimensionedScalar(dpBoundary.dimensions(), Zero);

    constrainPressure(dpBoundary, zeroU, dF2Boundary, rAtU);

    fvScalarMatrix dpBoundaryEqn
    (
        fvm::laplacian(rAtU, dpBoundary) == fvc::div(dF2Boundary)
    );
    dpBoundaryEqn.setReference(solverControl_().pRefCell(), Zero);
    dictionary dpBoundarySolver(dpBoundaryEqn.solverDict("p"));
    dpBoundarySolver.set("relTol", scalar(0));
    dpBoundarySolver.set("tolerance", scalar(1e-12));
    dpBoundaryEqn.solve(dpBoundarySolver);

    scalar pressureBoundaryLhs = Zero;
    scalar pressureBoundaryRhsInternal = Zero;
    scalar pressureBoundaryRhsWithBoundary = Zero;

    forAll(pSeed, celli)
    {
        pressureBoundaryLhs += pSeed[celli]*dpBoundary[celli];
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar FSeed = mu[P] - mu[N];

        pressureBoundaryRhsInternal += FSeed*dF2Boundary[facei];
    }

    pressureBoundaryRhsWithBoundary = pressureBoundaryRhsInternal;

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        const fvsPatchScalarField& dF2p =
            dF2Boundary.boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(dF2p, facei)
        {
            pressureBoundaryRhsWithBoundary +=
                mu[faceCells[facei]]*dF2p[facei];
        }
    }

    const scalar pressureBoundaryScale =
        max
        (
            max(mag(pressureBoundaryLhs), mag(pressureBoundaryRhsInternal)),
            VSMALL
        );
    const scalar pressureBoundaryWithPatchScale =
        max
        (
            max
            (
                mag(pressureBoundaryLhs),
                mag(pressureBoundaryRhsWithBoundary)
            ),
            VSMALL
        );

    Info<< "ATC-T pressure boundary transpose check: lhs = "
        << pressureBoundaryLhs
        << ", rhsInternal = " << pressureBoundaryRhsInternal
        << ", relInternal = "
        << mag(pressureBoundaryLhs - pressureBoundaryRhsInternal)
          /pressureBoundaryScale
        << ", rhsWithBoundary = " << pressureBoundaryRhsWithBoundary
        << ", relWithBoundary = "
        << mag(pressureBoundaryLhs - pressureBoundaryRhsWithBoundary)
          /pressureBoundaryWithPatchScale
        << ", nBoundaryFluxPerturb = " << nBoundaryFluxPerturb
        << endl;

    volScalarField q
    (
        IOobject
        (
            "ATCTconsistentQ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rAtU - rAU
    );

    volScalarField dpOld
    (
        IOobject
        (
            "paATCTConsistentOldPressure",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    dpOld = dimensionedScalar(dpOld.dimensions(), Zero);

    scalarField& dpOldI = dpOld.primitiveFieldRef();
    forAll(dpOldI, celli)
    {
        dpOldI[celli] =
            scalar((83*celli + 43) % 163)/scalar(163) - scalar(0.5);
    }
    dpOld.correctBoundaryConditions();

    vectorField dH0(mesh_.nCells(), vector::zero);
    forAll(dH0, celli)
    {
        const scalar x = scalar((89*celli + 47) % 167)/scalar(167) - scalar(0.5);
        const scalar y = scalar((97*celli + 53) % 173)/scalar(173) - scalar(0.5);
        const scalar z = scalar((101*celli + 59) % 179)/scalar(179) - scalar(0.5);
        dH0[celli] = vector(x, y, z);
    }

    scalarField dF1(mesh_.nInternalFaces(), Zero);
    forAll(dF1, facei)
    {
        dF1[facei] =
            scalar((103*facei + 61) % 181)/scalar(181) - scalar(0.5);
    }

    tmp<surfaceScalarField> tqf = fvc::interpolate(q);
    const surfaceScalarField& qf = tqf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceScalarField& deltaCoeffs = mesh_.deltaCoeffs();

    vectorField dHq(mesh_.nCells(), vector::zero);
    scalarField dN(mesh_.nInternalFaces(), Zero);
    scalarField pOldSeed(mesh_.nCells(), Zero);

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar wp = w[facei];
        const scalar wn = scalar(1) - wp;
        const scalar dpFace = wp*dpOld[P] + wn*dpOld[N];

        dHq[P] += q[P]*Sf[facei]*dpFace/V[P];
        dHq[N] -= q[N]*Sf[facei]*dpFace/V[N];

        const scalar eta =
            Sf[facei] & (q[P]*source[P] - q[N]*source[N]);

        pOldSeed[P] += wp*eta;
        pOldSeed[N] += wn*eta;

        const scalar nCoeff =
            qf[facei]*magSf[facei]*deltaCoeffs[facei];
        const scalar FSeed = mu[P] - mu[N];

        dN[facei] = nCoeff*(dpOld[N] - dpOld[P]);
        pOldSeed[P] -= nCoeff*FSeed;
        pOldSeed[N] += nCoeff*FSeed;
    }

    scalar consistentLhs = Zero;
    scalar consistentRhs = Zero;

    forAll(dH0, celli)
    {
        consistentLhs +=
            V[celli]*(source[celli] & (dH0[celli] + dHq[celli]));
        consistentRhs += V[celli]*(source[celli] & dH0[celli]);
        consistentRhs += pOldSeed[celli]*dpOld[celli];
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar FSeed = mu[P] - mu[N];

        consistentLhs += FSeed*(dF1[facei] + dN[facei]);
        consistentRhs += FSeed*dF1[facei];
    }

    const scalar consistentScale =
        max(max(mag(consistentLhs), mag(consistentRhs)), VSMALL);

    Info<< "ATC-T internal consistent-SIMPLE transpose check: lhs = "
        << consistentLhs
        << ", rhs = " << consistentRhs
        << ", rel = "
        << mag(consistentLhs - consistentRhs)/consistentScale
        << ", maxMag(q) = " << gMax(mag(q.primitiveField()))
        << endl;

    surfaceScalarField phiBeforeAdjust
    (
        IOobject
        (
            "ATCTphiBeforeAdjustProbe",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.phi()
    );

    scalar massIn = Zero;
    scalar fixedMassOut = Zero;
    scalar adjustableMassOut = Zero;
    scalar totalFlux = VSMALL + sum(mag(phiBeforeAdjust)).value();
    label nIn = 0;
    label nFixedOut = 0;
    label nAdjustableOut = 0;

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatchVectorField& Up =
            primalVars_.U().boundaryField()[patchi];
        const fvsPatchScalarField& phip =
            phiBeforeAdjust.boundaryField()[patchi];

        if (phip.coupled())
        {
            continue;
        }

        const bool fixedPatch =
            Up.fixesValue() && !isA<inletOutletFvPatchVectorField>(Up);

        forAll(phip, facei)
        {
            if (phip[facei] < 0)
            {
                massIn -= phip[facei];
                ++nIn;
            }
            else if (fixedPatch)
            {
                fixedMassOut += phip[facei];
                ++nFixedOut;
            }
            else
            {
                adjustableMassOut += phip[facei];
                ++nAdjustableOut;
            }
        }
    }

    scalar massCorr = scalar(1);
    bool adjustActive = false;

    if
    (
        mag(adjustableMassOut) > VSMALL
     && mag(adjustableMassOut)/totalFlux > SMALL
    )
    {
        massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        adjustActive = true;
    }

    surfaceScalarField phiAfterAdjust(phiBeforeAdjust);
    surfaceScalarField dPhiBeforeAdjust(phiBeforeAdjust);
    surfaceScalarField adjSeedAfterAdjust(phiBeforeAdjust);
    surfaceScalarField adjSeedBeforeAdjust(phiBeforeAdjust);

    dPhiBeforeAdjust = dimensionedScalar(dPhiBeforeAdjust.dimensions(), Zero);
    adjSeedAfterAdjust =
        dimensionedScalar(adjSeedAfterAdjust.dimensions(), Zero);
    adjSeedBeforeAdjust =
        dimensionedScalar(adjSeedBeforeAdjust.dimensions(), Zero);

    scalar barMassCorr = Zero;

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatchVectorField& Up =
            primalVars_.U().boundaryField()[patchi];
        fvsPatchScalarField& phiAfterPatch =
            phiAfterAdjust.boundaryFieldRef()[patchi];
        fvsPatchScalarField& dPhiPatch =
            dPhiBeforeAdjust.boundaryFieldRef()[patchi];
        fvsPatchScalarField& seedAfterPatch =
            adjSeedAfterAdjust.boundaryFieldRef()[patchi];
        fvsPatchScalarField& seedBeforePatch =
            adjSeedBeforeAdjust.boundaryFieldRef()[patchi];
        const fvsPatchScalarField& phiBeforePatch =
            phiBeforeAdjust.boundaryField()[patchi];

        if (phiBeforePatch.coupled())
        {
            continue;
        }

        const bool adjustablePatch =
            !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

        forAll(phiBeforePatch, facei)
        {
            dPhiPatch[facei] =
                scalar((107*(facei + 1) + 67*(patchi + 1)) % 191)
               /scalar(191) - scalar(0.5);

            seedAfterPatch[facei] =
                scalar((109*(facei + 1) + 71*(patchi + 1)) % 193)
               /scalar(193) - scalar(0.5);

            if
            (
                adjustActive
             && adjustablePatch
             && phiBeforePatch[facei] > 0
            )
            {
                phiAfterPatch[facei] = massCorr*phiBeforePatch[facei];
                barMassCorr += seedAfterPatch[facei]*phiBeforePatch[facei];
            }
            else
            {
                phiAfterPatch[facei] = phiBeforePatch[facei];
            }

            seedBeforePatch[facei] = seedAfterPatch[facei];
        }
    }

    if (adjustActive)
    {
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatchVectorField& Up =
                primalVars_.U().boundaryField()[patchi];
            const fvsPatchScalarField& phiBeforePatch =
                phiBeforeAdjust.boundaryField()[patchi];
            fvsPatchScalarField& seedBeforePatch =
                adjSeedBeforeAdjust.boundaryFieldRef()[patchi];

            if (phiBeforePatch.coupled())
            {
                continue;
            }

            const bool fixedPatch =
                Up.fixesValue()
             && !isA<inletOutletFvPatchVectorField>(Up);
            const bool adjustablePatch =
                !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

            forAll(phiBeforePatch, facei)
            {
                if (phiBeforePatch[facei] < 0)
                {
                    seedBeforePatch[facei] -=
                        barMassCorr/adjustableMassOut;
                }
                else if (fixedPatch)
                {
                    seedBeforePatch[facei] -=
                        barMassCorr/adjustableMassOut;
                }
                else if (adjustablePatch)
                {
                    seedBeforePatch[facei] =
                        massCorr*seedBeforePatch[facei]
                      - massCorr*barMassCorr/adjustableMassOut;
                }
            }
        }
    }

    scalar adjustPhiLhs = Zero;
    scalar adjustPhiRhs = Zero;
    scalar dMassIn = Zero;
    scalar dFixedMassOut = Zero;
    scalar dAdjustableMassOut = Zero;

    forAll(mesh_.boundary(), patchi)
    {
        const fvsPatchScalarField& phiBeforePatch =
            phiBeforeAdjust.boundaryField()[patchi];
        const fvsPatchScalarField& dPhiPatch =
            dPhiBeforeAdjust.boundaryField()[patchi];

        if (phiBeforePatch.coupled())
        {
            continue;
        }

        const fvPatchVectorField& Up =
            primalVars_.U().boundaryField()[patchi];
        const bool adjustablePatch =
            !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

        forAll(phiBeforePatch, facei)
        {
            if (phiBeforePatch[facei] < 0)
            {
                dMassIn -= dPhiPatch[facei];
            }
            else if (adjustablePatch)
            {
                dAdjustableMassOut += dPhiPatch[facei];
            }
            else
            {
                dFixedMassOut += dPhiPatch[facei];
            }
        }
    }

    scalar dMassCorr = Zero;
    if (adjustActive)
    {
        dMassCorr =
            (
                dMassIn
              - dFixedMassOut
              - massCorr*dAdjustableMassOut
            )/adjustableMassOut;
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvsPatchScalarField& phiBeforePatch =
            phiBeforeAdjust.boundaryField()[patchi];
        const fvsPatchScalarField& dPhiPatch =
            dPhiBeforeAdjust.boundaryField()[patchi];
        const fvsPatchScalarField& seedAfterPatch =
            adjSeedAfterAdjust.boundaryField()[patchi];
        const fvsPatchScalarField& seedBeforePatch =
            adjSeedBeforeAdjust.boundaryField()[patchi];

        if (phiBeforePatch.coupled())
        {
            continue;
        }

        const fvPatchVectorField& Up =
            primalVars_.U().boundaryField()[patchi];
        const bool adjustablePatch =
            !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

        forAll(phiBeforePatch, facei)
        {
            scalar dPhiAfter = dPhiPatch[facei];

            if
            (
                adjustActive
             && adjustablePatch
             && phiBeforePatch[facei] > 0
            )
            {
                dPhiAfter =
                    massCorr*dPhiPatch[facei]
                  + phiBeforePatch[facei]*dMassCorr;
            }

            adjustPhiLhs += seedAfterPatch[facei]*dPhiAfter;
            adjustPhiRhs += seedBeforePatch[facei]*dPhiPatch[facei];
        }
    }

    const scalar adjustPhiScale =
        max(max(mag(adjustPhiLhs), mag(adjustPhiRhs)), VSMALL);

    Info<< "ATC-T adjustPhi transpose check: lhs = " << adjustPhiLhs
        << ", rhs = " << adjustPhiRhs
        << ", rel = " << mag(adjustPhiLhs - adjustPhiRhs)/adjustPhiScale
        << ", active = " << adjustActive
        << ", massCorr = " << massCorr
        << ", nIn = " << nIn
        << ", nFixedOut = " << nFixedOut
        << ", nAdjustableOut = " << nAdjustableOut
        << endl;

    volVectorField dYHbyA
    (
        IOobject
        (
            "ATCTdYConstrainHbyA",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.U()
    );
    dYHbyA = dimensionedVector(dYHbyA.dimensions(), vector::zero);

    volVectorField dUHbyA
    (
        IOobject
        (
            "ATCTdUConstrainHbyA",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.U()
    );
    dUHbyA = dimensionedVector(dUHbyA.dimensions(), vector::zero);

    volVectorField HSeed
    (
        IOobject
        (
            "ATCTHSeedConstrainHbyA",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.U()
    );
    HSeed = dimensionedVector(HSeed.dimensions(), vector::zero);

    volVectorField dHConstrained(dYHbyA);
    volVectorField YSeed(dYHbyA);
    volVectorField USeed(dUHbyA);

    dHConstrained =
        dimensionedVector(dHConstrained.dimensions(), vector::zero);
    YSeed = dimensionedVector(YSeed.dimensions(), vector::zero);
    USeed = dimensionedVector(USeed.dimensions(), vector::zero);

    forAll(dYHbyA.primitiveFieldRef(), celli)
    {
        dYHbyA.primitiveFieldRef()[celli] =
            vector
            (
                scalar((113*celli + 73) % 197)/scalar(197) - scalar(0.5),
                scalar((127*celli + 79) % 199)/scalar(199) - scalar(0.5),
                scalar((131*celli + 83) % 211)/scalar(211) - scalar(0.5)
            );
        HSeed.primitiveFieldRef()[celli] =
            vector
            (
                scalar((137*celli + 89) % 223)/scalar(223) - scalar(0.5),
                scalar((139*celli + 97) % 227)/scalar(227) - scalar(0.5),
                scalar((149*celli + 101) % 229)/scalar(229) - scalar(0.5)
            );
        dHConstrained.primitiveFieldRef()[celli] = dYHbyA[celli];
        YSeed.primitiveFieldRef()[celli] = HSeed[celli];
    }

    label nHbyAOverwritten = 0;
    label nHbyAAssignable = 0;

    forAll(mesh_.boundary(), patchi)
    {
        fvPatchVectorField& dYp = dYHbyA.boundaryFieldRef()[patchi];
        fvPatchVectorField& dUp = dUHbyA.boundaryFieldRef()[patchi];
        fvPatchVectorField& seedHp = HSeed.boundaryFieldRef()[patchi];
        fvPatchVectorField& dHp = dHConstrained.boundaryFieldRef()[patchi];
        fvPatchVectorField& seedYp = YSeed.boundaryFieldRef()[patchi];
        fvPatchVectorField& seedUp = USeed.boundaryFieldRef()[patchi];

        const fvPatchVectorField& Up =
            primalVars_.U().boundaryField()[patchi];
        const fvPatchScalarField& pp =
            primalVars_.p().boundaryField()[patchi];

        const bool overwritten =
            !Up.assignable()
         && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>(pp);

        forAll(dYp, facei)
        {
            dYp[facei] =
                vector
                (
                    scalar((151*(facei + 1) + 103*(patchi + 1)) % 233)
                   /scalar(233) - scalar(0.5),
                    scalar((157*(facei + 1) + 107*(patchi + 1)) % 239)
                   /scalar(239) - scalar(0.5),
                    scalar((163*(facei + 1) + 109*(patchi + 1)) % 241)
                   /scalar(241) - scalar(0.5)
                );
            dUp[facei] =
                vector
                (
                    scalar((167*(facei + 1) + 113*(patchi + 1)) % 251)
                   /scalar(251) - scalar(0.5),
                    scalar((173*(facei + 1) + 127*(patchi + 1)) % 257)
                   /scalar(257) - scalar(0.5),
                    scalar((179*(facei + 1) + 131*(patchi + 1)) % 263)
                   /scalar(263) - scalar(0.5)
                );
            seedHp[facei] =
                vector
                (
                    scalar((181*(facei + 1) + 137*(patchi + 1)) % 269)
                   /scalar(269) - scalar(0.5),
                    scalar((191*(facei + 1) + 139*(patchi + 1)) % 271)
                   /scalar(271) - scalar(0.5),
                    scalar((193*(facei + 1) + 149*(patchi + 1)) % 277)
                   /scalar(277) - scalar(0.5)
                );

            if (overwritten)
            {
                dHp[facei] = dUp[facei];
                seedYp[facei] = vector::zero;
                seedUp[facei] = seedHp[facei];
                ++nHbyAOverwritten;
            }
            else
            {
                dHp[facei] = dYp[facei];
                seedYp[facei] = seedHp[facei];
                seedUp[facei] = vector::zero;
                ++nHbyAAssignable;
            }
        }
    }

    scalar constrainHbyALhs = Zero;
    scalar constrainHbyARhs = Zero;

    forAll(dYHbyA, celli)
    {
        constrainHbyALhs += HSeed[celli] & dHConstrained[celli];
        constrainHbyARhs += YSeed[celli] & dYHbyA[celli];
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatchVectorField& dYp = dYHbyA.boundaryField()[patchi];
        const fvPatchVectorField& dUp = dUHbyA.boundaryField()[patchi];
        const fvPatchVectorField& seedHp = HSeed.boundaryField()[patchi];
        const fvPatchVectorField& dHp =
            dHConstrained.boundaryField()[patchi];
        const fvPatchVectorField& seedYp = YSeed.boundaryField()[patchi];
        const fvPatchVectorField& seedUp = USeed.boundaryField()[patchi];

        forAll(dYp, facei)
        {
            constrainHbyALhs += seedHp[facei] & dHp[facei];
            constrainHbyARhs +=
                (seedYp[facei] & dYp[facei])
              + (seedUp[facei] & dUp[facei]);
        }
    }

    const scalar constrainHbyAScale =
        max(max(mag(constrainHbyALhs), mag(constrainHbyARhs)), VSMALL);

    Info<< "ATC-T constrainHbyA transpose check: lhs = "
        << constrainHbyALhs
        << ", rhs = " << constrainHbyARhs
        << ", rel = "
        << mag(constrainHbyALhs - constrainHbyARhs)/constrainHbyAScale
        << ", overwrittenFaces = " << nHbyAOverwritten
        << ", assignableFaces = " << nHbyAAssignable
        << endl;

    volVectorField& UForPredictor = primalVars_.U();
    surfaceScalarField& phiForPredictor = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    const volVectorField UBase(UForPredictor);

    volVectorField dUPredictor
    (
        IOobject
        (
            "ATCTdUPredictorMap",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        UBase
    );
    dUPredictor = dimensionedVector(dUPredictor.dimensions(), vector::zero);

    forAll(dUPredictor.primitiveFieldRef(), celli)
    {
        dUPredictor.primitiveFieldRef()[celli] =
            vector
            (
                scalar((199*celli + 151) % 283)/scalar(283) - scalar(0.5),
                scalar((211*celli + 157) % 293)/scalar(293) - scalar(0.5),
                scalar((223*celli + 163) % 307)/scalar(307) - scalar(0.5)
            );
    }

    forAll(dUPredictor.boundaryFieldRef(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() != "empty"
         && !mesh_.boundary()[patchi].coupled()
        )
        {
            dUPredictor.boundaryFieldRef()[patchi] == vector::zero;
        }
    }

    const scalar predictorEps = scalar(1e-6);

    auto setPredictorU = [&](const scalar coeff)
    {
        UForPredictor = UBase;
        UForPredictor.primitiveFieldRef() =
            UBase.primitiveField()
          + coeff*dUPredictor.primitiveField();

        forAll(UForPredictor.boundaryFieldRef(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() != "empty"
             && !mesh_.boundary()[patchi].coupled()
            )
            {
                UForPredictor.boundaryFieldRef()[patchi] ==
                    UBase.boundaryField()[patchi];
            }
        }
    };

    auto predictorY = [&](const word& fieldName)
    {
        tmp<fvVectorMatrix> tPredUEqn
        (
            fvm::div(phiForPredictor, UForPredictor)
          + turbulence->divDevReff(UForPredictor)
         ==
            fvOptions(UForPredictor)
        );
        fvVectorMatrix& predUEqn = tPredUEqn.ref();

        predUEqn.relax();
        fvOptions.constrain(predUEqn);

        volScalarField predRAU(1.0/predUEqn.A());
        tmp<volVectorField> tY(predRAU*predUEqn.H());
        tY.ref().rename(fieldName);
        return tY;
    };

    setPredictorU(predictorEps);
    tmp<volVectorField> tYPlus = predictorY("ATCTYPlusPredictorMap");
    vectorField YPlus(tYPlus().primitiveField());

    setPredictorU(-predictorEps);
    tmp<volVectorField> tYMinus = predictorY("ATCTYMinusPredictorMap");
    vectorField YMinus(tYMinus().primitiveField());

    setPredictorU(Zero);

    tmp<fvVectorMatrix> tFrozenUEqn
    (
        fvm::div(phiForPredictor, UForPredictor)
      + turbulence->divDevReff(UForPredictor)
     ==
        fvOptions(UForPredictor)
    );
    fvVectorMatrix& frozenUEqn = tFrozenUEqn.ref();

    const scalarField frozenDiag0(frozenUEqn.diag());
    frozenUEqn.relax();
    const scalarField frozenDiagRelaxed(frozenUEqn.diag());
    fvOptions.constrain(frozenUEqn);

    volScalarField frozenRAU(1.0/frozenUEqn.A());

    setPredictorU(predictorEps);
    tmp<volVectorField> tYFrozenPlus(frozenRAU*frozenUEqn.H());
    vectorField YFrozenPlus(tYFrozenPlus().primitiveField());

    setPredictorU(-predictorEps);
    tmp<volVectorField> tYFrozenMinus(frozenRAU*frozenUEqn.H());
    vectorField YFrozenMinus(tYFrozenMinus().primitiveField());

    UForPredictor = UBase;

    scalar predictorFDSourceSeed = Zero;
    scalar predictorFrozenSourceSeed = Zero;
    scalar predictorFrozenRelaxSourceSeed = Zero;
    scalar predictorFDRandomSeed = Zero;
    scalar predictorFrozenRandomSeed = Zero;
    scalar predictorFrozenRelaxRandomSeed = Zero;
    scalar predictorDiffL1 = Zero;
    scalar predictorRelaxDiffL1 = Zero;
    scalar predictorFDL1 = VSMALL;

    forAll(YPlus, celli)
    {
        const vector dYFD =
            (YPlus[celli] - YMinus[celli])/(2*predictorEps);
        const vector dYFrozen =
            (YFrozenPlus[celli] - YFrozenMinus[celli])/(2*predictorEps);
        const vector dYFrozenRelax =
            dYFrozen
          + frozenRAU[celli]
           *(frozenDiagRelaxed[celli] - frozenDiag0[celli])
           *dUPredictor[celli]/V[celli];

        const vector sourceSeed = V[celli]*source[celli];
        const vector randomSeed
        (
            scalar((227*celli + 167) % 311)/scalar(311) - scalar(0.5),
            scalar((229*celli + 173) % 313)/scalar(313) - scalar(0.5),
            scalar((233*celli + 179) % 317)/scalar(317) - scalar(0.5)
        );

        predictorFDSourceSeed += sourceSeed & dYFD;
        predictorFrozenSourceSeed += sourceSeed & dYFrozen;
        predictorFrozenRelaxSourceSeed += sourceSeed & dYFrozenRelax;
        predictorFDRandomSeed += randomSeed & dYFD;
        predictorFrozenRandomSeed += randomSeed & dYFrozen;
        predictorFrozenRelaxRandomSeed += randomSeed & dYFrozenRelax;
        predictorDiffL1 += mag(dYFD - dYFrozen);
        predictorRelaxDiffL1 += mag(dYFD - dYFrozenRelax);
        predictorFDL1 += mag(dYFD);
    }

    const scalar predictorSourceScale =
        max
        (
            max
            (
                max
                (
                    mag(predictorFDSourceSeed),
                    mag(predictorFrozenSourceSeed)
                ),
                mag(predictorFrozenRelaxSourceSeed)
            ),
            VSMALL
        );
    const scalar predictorRandomScale =
        max
        (
            max
            (
                max
                (
                    mag(predictorFDRandomSeed),
                    mag(predictorFrozenRandomSeed)
                ),
                mag(predictorFrozenRelaxRandomSeed)
            ),
            VSMALL
        );

    Info<< "ATC-T predictor map action check: sourceSeedFD = "
        << predictorFDSourceSeed
        << ", sourceSeedFrozen = " << predictorFrozenSourceSeed
        << ", sourceSeedFrozenRelax = " << predictorFrozenRelaxSourceSeed
        << ", relSource = "
        << mag(predictorFDSourceSeed - predictorFrozenSourceSeed)
          /predictorSourceScale
        << ", relSourceRelax = "
        << mag(predictorFDSourceSeed - predictorFrozenRelaxSourceSeed)
          /predictorSourceScale
        << ", randomSeedFD = " << predictorFDRandomSeed
        << ", randomSeedFrozen = " << predictorFrozenRandomSeed
        << ", randomSeedFrozenRelax = " << predictorFrozenRelaxRandomSeed
        << ", relRandom = "
        << mag(predictorFDRandomSeed - predictorFrozenRandomSeed)
          /predictorRandomScale
        << ", relRandomRelax = "
        << mag(predictorFDRandomSeed - predictorFrozenRelaxRandomSeed)
          /predictorRandomScale
        << ", relL1(dYFD-dYFrozen) = " << predictorDiffL1/predictorFDL1
        << ", relL1(dYFD-dYFrozenRelax) = "
        << predictorRelaxDiffL1/predictorFDL1
        << endl;

    vectorField predictorSourceSeed(mesh_.nCells(), vector::zero);
    vectorField predictorRandomSeed(mesh_.nCells(), vector::zero);
    vectorField predictorReverseSeed(mesh_.nCells(), vector::zero);

    forAll(predictorSourceSeed, celli)
    {
        predictorSourceSeed[celli] = V[celli]*source[celli];
        predictorReverseSeed[celli] = predictorSourceSeed[celli];
        predictorRandomSeed[celli] =
            vector
            (
                scalar((227*celli + 167) % 311)/scalar(311) - scalar(0.5),
                scalar((229*celli + 173) % 313)/scalar(313) - scalar(0.5),
                scalar((233*celli + 179) % 317)/scalar(317) - scalar(0.5)
            );
    }

    surfaceScalarField predictorFSeedAfterAdjust
    (
        IOobject
        (
            "ATCTPredictorFSeedAfterAdjust",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiBeforeAdjust
    );
    predictorFSeedAfterAdjust =
        dimensionedScalar(predictorFSeedAfterAdjust.dimensions(), Zero);

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        predictorFSeedAfterAdjust.primitiveFieldRef()[facei] = mu[P] - mu[N];
    }

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        fvsPatchScalarField& seedp =
            predictorFSeedAfterAdjust.boundaryFieldRef()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(seedp, facei)
        {
            seedp[facei] = mu[faceCells[facei]];
        }
    }

    surfaceScalarField predictorFSeedBeforeAdjust
    (
        predictorFSeedAfterAdjust
    );

    if (adjustActive)
    {
        scalar predictorBarMassCorr = Zero;

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatchVectorField& Up =
                primalVars_.U().boundaryField()[patchi];
            const fvsPatchScalarField& phiBeforePatch =
                phiBeforeAdjust.boundaryField()[patchi];
            const fvsPatchScalarField& seedAfterPatch =
                predictorFSeedAfterAdjust.boundaryField()[patchi];

            if (phiBeforePatch.coupled())
            {
                continue;
            }

            const bool adjustablePatch =
                !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

            forAll(phiBeforePatch, facei)
            {
                if
                (
                    adjustablePatch
                 && phiBeforePatch[facei] > 0
                )
                {
                    predictorBarMassCorr +=
                        seedAfterPatch[facei]*phiBeforePatch[facei];
                }
            }
        }

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatchVectorField& Up =
                primalVars_.U().boundaryField()[patchi];
            const fvsPatchScalarField& phiBeforePatch =
                phiBeforeAdjust.boundaryField()[patchi];
            fvsPatchScalarField& seedBeforePatch =
                predictorFSeedBeforeAdjust.boundaryFieldRef()[patchi];

            if (phiBeforePatch.coupled())
            {
                continue;
            }

            const bool fixedPatch =
                Up.fixesValue()
             && !isA<inletOutletFvPatchVectorField>(Up);
            const bool adjustablePatch =
                !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

            forAll(phiBeforePatch, facei)
            {
                if (phiBeforePatch[facei] < 0)
                {
                    seedBeforePatch[facei] -=
                        predictorBarMassCorr/adjustableMassOut;
                }
                else if (fixedPatch)
                {
                    seedBeforePatch[facei] -=
                        predictorBarMassCorr/adjustableMassOut;
                }
                else if (adjustablePatch)
                {
                    seedBeforePatch[facei] =
                        massCorr*seedBeforePatch[facei]
                      - massCorr*predictorBarMassCorr/adjustableMassOut;
                }
            }
        }
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const vector FSeedSf =
            predictorFSeedBeforeAdjust[facei]*Sf[facei];

        predictorReverseSeed[P] += w[facei]*FSeedSf;
        predictorReverseSeed[N] += (scalar(1) - w[facei])*FSeedSf;
    }

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        const fvPatchVectorField& Up =
            primalVars_.U().boundaryField()[patchi];
        const fvPatchScalarField& pp =
            primalVars_.p().boundaryField()[patchi];
        const bool overwritten =
            !Up.assignable()
         && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>(pp);

        if (overwritten)
        {
            continue;
        }

        tmp<vectorField> tValueInternalCoeffs =
            Up.valueInternalCoeffs(mesh_.boundary()[patchi].weights());
        const vectorField& valueInternalCoeffs = tValueInternalCoeffs();
        const fvsPatchScalarField& seedPatch =
            predictorFSeedBeforeAdjust.boundaryField()[patchi];
        const fvsPatchVectorField& SfPatch =
            Sf.boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        forAll(seedPatch, facei)
        {
            predictorReverseSeed[faceCells[facei]] +=
                seedPatch[facei]
               *cmptMultiply(valueInternalCoeffs[facei], SfPatch[facei]);
        }
    }

    checkPredictorMapTranspose
    (
        predictorSourceSeed,
        "sourceSeed_P0",
        false,
        false
    );
    checkPredictorMapTranspose
    (
        predictorSourceSeed,
        "sourceSeed_P1_relax",
        true,
        false
    );
    checkPredictorMapTranspose
    (
        predictorSourceSeed,
        "sourceSeed_P2_relaxFvOptions",
        true,
        true
    );
    checkPredictorMapTranspose
    (
        predictorRandomSeed,
        "randomSeed_P2_relaxFvOptions",
        true,
        true
    );
    checkPredictorMapTranspose
    (
        predictorReverseSeed,
        "reverseSeed_P2_relaxFvOptions",
        true,
        true
    );

    if (checkPressureMapTranspose_)
    {
        checkPressureMapTranspose(predictorSourceSeed, rAtU);
    }

    if (checkStateMapTranspose_)
    {
        checkStateMapTranspose(predictorSourceSeed, rAtU);
    }

    if (usePredictorReverseMomentumSens_)
    {
        tmp<vectorField> tUbarExactInt =
            reversePredictorY(predictorReverseSeed, true, true);
        const vectorField& UbarExactInt = tUbarExactInt();

        const incompressibleAdjointVars& adjointVars = getAdjointVars();
        const volVectorField& Ua = adjointVars.Ua();

        volVectorField UbarExact
        (
            IOobject
            (
                "ATCTUbarExactPredictorReverse",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Ua
        );

        forAll(UbarExact.primitiveFieldRef(), celli)
        {
            UbarExact.primitiveFieldRef()[celli] =
                UbarExactInt[celli]/V[celli];
        }
        UbarExact.correctBoundaryConditions();
        UbarExact.write();

        scalar dotUa = Zero;
        scalar normExact = VSMALL;
        scalar normUa = VSMALL;
        scalar diffNorm = Zero;

        forAll(UbarExact, celli)
        {
            dotUa += V[celli]*(UbarExact[celli] & Ua[celli]);
            normExact += V[celli]*magSqr(UbarExact[celli]);
            normUa += V[celli]*magSqr(Ua[celli]);
            diffNorm += V[celli]*magSqr(UbarExact[celli] - Ua[celli]);
        }

        Info<< "ATC-T predictor reverse momentum seed written: "
            << UbarExact.name()
            << ", cos(UbarExact,Ua) = "
            << dotUa/sqrt(normExact*normUa)
            << ", relL2(UbarExact-Ua) = "
            << sqrt(diffNorm)/sqrt(normExact)
            << ". Diagnostic only; betaMult is unchanged."
            << endl;
    }

    if (checkOneStepBetaMapSensitivity_)
    {
        tmp<vectorField> tUbarExactInt =
            reversePredictorY(predictorReverseSeed, true, true);

        checkOneStepBetaMapSensitivity
        (
            predictorSourceSeed,
            tUbarExactInt(),
            "beta",
            scalar(1)
        );
    }

    if (checkFixedPointMapAdjoint_)
    {
        checkFixedPointMapAdjointSensitivity
        (
            predictorSourceSeed,
            rAtU,
            "beta",
            scalar(1)
        );
    }

    if (checkFixedPointTangentSensitivity_)
    {
        checkFixedPointTangentSensitivity
        (
            predictorSourceSeed,
            "beta",
            scalar(1)
        );
    }
}


void Foam::thermalAdjointSimple::checkFlowBlockTranspose
(
    const volVectorField& source
)
{
    if (!checkProjectionTranspose_)
    {
        return;
    }

    Info<< "WARNING: ATC-T flow-block transpose check is deprecated. "
        << "It is a hybrid residual/SIMPLE-map diagnostic, not an exact "
        << "transpose pass/fail gate. Use the predictor FD transpose checks "
        << "for the current target." << endl;

    const incompressibleAdjointVars& adjointVars = getAdjointVars();
    const volVectorField& Ua = adjointVars.Ua();
    const volScalarField& pa = adjointVars.pa();

    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
    fvOptions.constrain(UEqn);

    volVectorField dU
    (
        IOobject
        (
            "ATCTdUFlowBlock",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.U()
    );
    dU = dimensionedVector(dU.dimensions(), Zero);

    vectorField& dUI = dU.primitiveFieldRef();
    forAll(dUI, celli)
    {
        const scalar x = scalar((53*celli + 23) % 127)/scalar(127) - scalar(0.5);
        const scalar y = scalar((59*celli + 29) % 131)/scalar(131) - scalar(0.5);
        const scalar z = scalar((61*celli + 31) % 137)/scalar(137) - scalar(0.5);
        dUI[celli] = vector(x, y, z);
    }

    dU.correctBoundaryConditions();

    volVectorField::Boundary& dUb = dU.boundaryFieldRef();
    forAll(dUb, patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() != "empty"
         && !mesh_.boundary()[patchi].coupled()
        )
        {
            dUb[patchi] == vector::zero;
        }
    }

    volScalarField dp
    (
        IOobject
        (
            "ATCTdpFlowBlock",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.p()
    );
    dp = dimensionedScalar(dp.dimensions(), Zero);

    scalarField& dpI = dp.primitiveFieldRef();
    forAll(dpI, celli)
    {
        dpI[celli] =
            scalar((67*celli + 37) % 139)/scalar(139) - scalar(0.5);
    }

    dp.correctBoundaryConditions();

    volScalarField::Boundary& dpb = dp.boundaryFieldRef();
    forAll(dpb, patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() != "empty"
         && !mesh_.boundary()[patchi].coupled()
        )
        {
            dpb[patchi] == scalar(0);
        }
    }

    auto integratedVectorAction =
        [this](fvVectorMatrix& matrix, const volVectorField& field)
        {
            vectorField action
            (
                matrix.diag()*field.primitiveField()
            );

            const FieldField<Field, vector>& intCoeffs =
                matrix.internalCoeffs();
            forAll(intCoeffs, patchi)
            {
                const labelUList& faceCells =
                    mesh_.boundary()[patchi].faceCells();

                forAll(intCoeffs[patchi], facei)
                {
                    const label celli = faceCells[facei];
                    action[celli] += cmptMultiply
                    (
                        intCoeffs[patchi][facei],
                        field[celli]
                    );
                }
            }

            const scalarField& upper = matrix.upper();
            const scalarField& lower = matrix.lower();
            const labelUList& own = mesh_.owner();
            const labelUList& nei = mesh_.neighbour();

            forAll(own, facei)
            {
                action[own[facei]] += upper[facei]*field[nei[facei]];
                action[nei[facei]] += lower[facei]*field[own[facei]];
            }

            return action;
        };

    tmp<fvVectorMatrix> tRawUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& rawUEqn = tRawUEqn.ref();

    const vectorField rawAction(integratedVectorAction(rawUEqn, dU));
    const scalarField rawDiag(rawUEqn.diag());

    rawUEqn.relax();
    const scalarField relaxDiagDelta(rawUEqn.diag() - rawDiag);
    fvOptions.constrain(rawUEqn);

    const vectorField relaxedAction(integratedVectorAction(rawUEqn, dU));
    vectorField relaxSourceDerivative(mesh_.nCells(), vector::zero);
    forAll(relaxSourceDerivative, celli)
    {
        relaxSourceDerivative[celli] = relaxDiagDelta[celli]*dU[celli];
    }

    // Apply only the homogeneous tangent matrix to the perturbation.
    UEqn.source() = vector::zero;
    FieldField<Field, vector>& bCoeffs = UEqn.boundaryCoeffs();
    forAll(bCoeffs, patchi)
    {
        bCoeffs[patchi] = vector::zero;
    }

    tmp<volVectorField> tMdU =
        UEqn & static_cast<const DimensionedField<vector, volMesh>&>(dU);

    volVectorField deltaRU
    (
        IOobject
        (
            "ATCTdeltaRU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tMdU() + fvc::grad(dp)
    );

    volVectorField gradDp
    (
        IOobject
        (
            "ATCTgradDpFlowBlock",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(dp)
    );

    surfaceScalarField dPhiU("ATCTdPhiU", fvc::flux(dU));
    volScalarField deltaRc("ATCTdeltaRc", fvc::div(dPhiU));

    volVectorField deltaRUViscous
    (
        IOobject
        (
            "ATCTdeltaRUViscous",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       -fvc::div
        (
            turbulence->nuEff()()*dev2(T(fvc::grad(dU))),
            "div((nuEff*dev2(T(grad(U)))))"
        )
    );

    volVectorField deltaRUFlux
    (
        IOobject
        (
            "ATCTdeltaRUFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(dPhiU, U, "div(phi,U)")
    );

    const scalarField& V = mesh_.V();

    tmp<surfaceScalarField> tHPhi = momentumFluxSensitivity();
    tmp<volVectorField> tExactATC =
        projectedFluxMomentumSource(tHPhi(), "ATCMomCheck");
    tmp<volVectorField> tExactConvT = exactConvectionTransposeSource();

    volVectorField standardATC
    (
        IOobject
        (
            "ATCTstandardATC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(U, "gradUATC") & Ua
    );

    scalar objTerm = Zero;
    scalar momTerm = Zero;
    scalar momViscTerm = Zero;
    scalar momFluxTerm = Zero;
    scalar momStandardATCTerm = Zero;
    scalar momExactATCTerm = Zero;
    scalar momExactConvTTerm = Zero;
    scalar momFluxFaceTerm = Zero;
    scalar contTerm = Zero;
    scalar momRawMatrixTerm = Zero;
    scalar momRelaxedMatrixTerm = Zero;
    scalar momRelaxedWithSourceDerivativeTerm = Zero;

    forAll(V, celli)
    {
        objTerm +=
            V[celli]*couplingSign_*(source[celli] & dU[celli]);
        momTerm += V[celli]*(Ua[celli] & deltaRU[celli]);
        momViscTerm +=
            V[celli]*(Ua[celli] & (deltaRU[celli] + deltaRUViscous[celli]));
        momFluxTerm +=
            V[celli]
           *(
                Ua[celli]
              & (
                    deltaRU[celli]
                  + deltaRUViscous[celli]
                  + deltaRUFlux[celli]
                )
            );
        momStandardATCTerm +=
            V[celli]*(standardATC[celli] & dU[celli]);
        momExactATCTerm +=
            V[celli]*(tExactATC()[celli] & dU[celli]);
        momExactConvTTerm +=
            V[celli]*(tExactConvT()[celli] & dU[celli]);
        contTerm += V[celli]*pa[celli]*deltaRc[celli];
        momRawMatrixTerm +=
            (Ua[celli] & rawAction[celli])
          + V[celli]*(Ua[celli] & gradDp[celli]);
        momRelaxedMatrixTerm +=
            (Ua[celli] & relaxedAction[celli])
          + V[celli]*(Ua[celli] & gradDp[celli]);
        momRelaxedWithSourceDerivativeTerm +=
            (
                Ua[celli]
              & (
                    relaxedAction[celli]
                  - relaxSourceDerivative[celli]
                )
            )
          + V[celli]*(Ua[celli] & gradDp[celli]);
    }

    forAll(dPhiU.primitiveField(), facei)
    {
        momFluxFaceTerm += tHPhi()[facei]*dPhiU[facei];
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        const fvsPatchScalarField& hPatch = tHPhi().boundaryField()[patchi];
        const fvsPatchScalarField& dPhiPatch = dPhiU.boundaryField()[patchi];

        forAll(hPatch, facei)
        {
            momFluxFaceTerm += hPatch[facei]*dPhiPatch[facei];
        }
    }

    auto reportFlowBlock =
        [&](const word& name, const scalar mom)
        {
            const scalar blockScale =
                max(mag(objTerm) + mag(mom) + mag(contTerm), VSMALL);

            Info<< "ATC-T flow-block transpose check (" << name
                << "): obj = " << objTerm
                << ", mom = " << mom
                << ", cont = " << contTerm
                << ", rel(obj+mom+cont) = "
                << mag(objTerm + mom + contTerm)/blockScale
                << ", rel(obj+mom-cont) = "
                << mag(objTerm + mom - contTerm)/blockScale
                << ", rel(obj-mom+cont) = "
                << mag(objTerm - mom + contTerm)/blockScale
                << ", rel(obj-mom-cont) = "
                << mag(objTerm - mom - contTerm)/blockScale
                << endl;
        };

    reportFlowBlock("frozenMatrix", momTerm);
    reportFlowBlock("withExplicitViscous", momViscTerm);
    reportFlowBlock("withFluxConvection", momFluxTerm);
    reportFlowBlock("withStandardATC", momViscTerm + momStandardATCTerm);
    reportFlowBlock("withExactMomentumATC", momViscTerm + momExactATCTerm);
    reportFlowBlock("withExactConvT", momTerm + momExactConvTTerm);
    reportFlowBlock
    (
        "withExactConvTAndMomentumATC",
        momTerm + momExactConvTTerm + momExactATCTerm
    );
    reportFlowBlock("rawExplicitMatrix", momRawMatrixTerm);
    reportFlowBlock("relaxedExplicitMatrix", momRelaxedMatrixTerm);
    reportFlowBlock
    (
        "relaxedExplicitMatrixWithRelaxSourceDerivative",
        momRelaxedWithSourceDerivativeTerm
    );

    Info<< "ATC-T flow-block pieces: obj = " << objTerm
        << ", momFrozen = " << momTerm
        << ", momWithViscous = " << momViscTerm
        << ", momWithFluxConvection = " << momFluxTerm
        << ", standardATC = " << momStandardATCTerm
        << ", exactMomentumATC = " << momExactATCTerm
        << ", exactConvT = " << momExactConvTTerm
        << ", rawFluxFace = " << momFluxFaceTerm
        << ", cont = " << contTerm
        << ", rawExplicitMatrix = " << momRawMatrixTerm
        << ", relaxedExplicitMatrix = " << momRelaxedMatrixTerm
        << ", relaxedExplicitMatrixWithRelaxSourceDerivative = "
        << momRelaxedWithSourceDerivativeTerm
        << endl;

    volScalarField da
    (
        IOobject
        (
            "ATCTdeltaA",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    scalarField& daI = da.primitiveFieldRef();
    forAll(daI, celli)
    {
        daI[celli] =
            scalar((71*celli + 41) % 149)/scalar(149) - scalar(0.5);
    }

    tmp<fvVectorMatrix> tSpEqn(fvm::Sp(da, U));
    fvVectorMatrix& spEqn = tSpEqn.ref();
    tmp<volVectorField> tSpU =
        spEqn & static_cast<const DimensionedField<vector, volMesh>&>(U);

    scalar spAssembled = Zero;
    scalar spDirect = Zero;

    forAll(V, celli)
    {
        spAssembled += V[celli]*(Ua[celli] & tSpU()[celli]);
        spDirect += V[celli]*da[celli]*(U[celli] & Ua[celli]);
    }

    const scalar spScale =
        max(max(mag(spAssembled), mag(spDirect)), VSMALL);

    Info<< "ATC-T Brinkman Sp residual check: assembled = "
        << spAssembled
        << ", direct = " << spDirect
        << ", rel = " << mag(spAssembled - spDirect)/spScale
        << endl;
}


Foam::tmp<Foam::vectorField>
Foam::thermalAdjointSimple::primalMomentumResidualLHS()
{
    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
    fvOptions.constrain(UEqn);

    tmp<vectorField> tResidual = UEqn.residual();
    auto tLHS = tmp<vectorField>::New(tResidual());
    tLHS.ref() *= -1;

    return tLHS;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::thermalAdjointSimple::primalProjectedFluxAtFrozenState
(
    const word& namePrefix
)
{
    volScalarField& p = primalVars_.p();
    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
    fvOptions.constrain(UEqn);

    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA(namePrefix + "phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> trAtU(rAU);

    if (solverControl_().consistent())
    {
        trAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(trAtU() - rAU)*fvc::snGrad(p)*mesh_.magSf();
        HbyA -= (rAU - trAtU())*fvc::grad(p);
    }

    const volScalarField& rAtU = trAtU();
    constrainPressure(p, U, phiHbyA, rAtU);

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU, p) == fvc::div(phiHbyA)
    );
    pEqn.setReference(solverControl_().pRefCell(), solverControl_().pRefValue());

    auto tProjectedPhi =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                namePrefix + "phiProjected",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phiHbyA
        );
    tProjectedPhi.ref() -= pEqn.flux();

    return tProjectedPhi;
}


Foam::tmp<Foam::volVectorField>
Foam::thermalAdjointSimple::primalSimpleMapUAtFrozenState
(
    const word& namePrefix
)
{
    const volScalarField& p = primalVars_.p();
    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
    fvOptions.constrain(UEqn);

    volScalarField rAU(1.0/UEqn.A());

    volScalarField pMap
    (
        IOobject
        (
            namePrefix + "pMap",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    pMap.storePrevIter();

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, pMap));
    surfaceScalarField phiHbyA(namePrefix + "phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, U, pMap);

    tmp<volScalarField> trAtU(rAU);

    if (solverControl_().consistent())
    {
        trAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(trAtU() - rAU)*fvc::snGrad(pMap)*mesh_.magSf();
        HbyA -= (rAU - trAtU())*fvc::grad(pMap);
    }

    const volScalarField& rAtU = trAtU();
    constrainPressure(pMap, U, phiHbyA, rAtU);
    mesh_.setFluxRequired(pMap.name());

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU, pMap) == fvc::div(phiHbyA)
    );
    pEqn.setReference
    (
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );
    dictionary pSolver(pEqn.solverDict("p"));
    pSolver.set("relTol", scalar(0));
    pSolver.set("tolerance", scalar(1e-12));
    pEqn.solve(pSolver);

    scalar pRelaxCoeff = scalar(1);
    word pRelaxName = p.name();
    if (p.mesh().data().isFinalIteration())
    {
        pRelaxName += "Final";
    }
    p.mesh().relaxField(pRelaxName, pRelaxCoeff);
    pMap.relax(pRelaxCoeff);

    auto tUnew =
        tmp<volVectorField>::New
        (
            IOobject
            (
                namePrefix + "Unew",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            HbyA - rAtU*fvc::grad(pMap)
    );

    tUnew.ref().correctBoundaryConditions();
    fvOptions.correct(tUnew.ref());

    return tUnew;
}


Foam::thermalAdjointSimple::SimpleMapState
Foam::thermalAdjointSimple::primalSimpleMapStateAtFrozenState
(
    const word& namePrefix
)
{
    const volScalarField& p = primalVars_.p();
    volVectorField& U = primalVars_.U();
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
    fvOptions.constrain(UEqn);

    volScalarField rAU(1.0/UEqn.A());

    volScalarField pMap
    (
        IOobject
        (
            namePrefix + "pMap",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    pMap.storePrevIter();

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, pMap));
    surfaceScalarField phiHbyA(namePrefix + "phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, U, pMap);

    tmp<volScalarField> trAtU(rAU);

    if (solverControl_().consistent())
    {
        trAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(trAtU() - rAU)*fvc::snGrad(pMap)*mesh_.magSf();
        HbyA -= (rAU - trAtU())*fvc::grad(pMap);
    }

    const volScalarField& rAtU = trAtU();
    constrainPressure(pMap, U, phiHbyA, rAtU);
    mesh_.setFluxRequired(pMap.name());

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU, pMap) == fvc::div(phiHbyA)
    );
    pEqn.setReference
    (
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );
    dictionary pSolver(pEqn.solverDict("p"));
    pSolver.set("relTol", scalar(0));
    pSolver.set("tolerance", scalar(1e-12));
    pEqn.solve(pSolver);

    surfaceScalarField phiNew
    (
        IOobject
        (
            namePrefix + "phiNew",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHbyA - pEqn.flux()
    );

    scalar pRelaxCoeff = scalar(1);
    word pRelaxName = p.name();
    if (p.mesh().data().isFinalIteration())
    {
        pRelaxName += "Final";
    }
    p.mesh().relaxField(pRelaxName, pRelaxCoeff);
    pMap.relax(pRelaxCoeff);

    volVectorField Unew
    (
        IOobject
        (
            namePrefix + "Unew",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        HbyA - rAtU*fvc::grad(pMap)
    );
    Unew.correctBoundaryConditions();
    fvOptions.correct(Unew);

    SimpleMapState state(mesh_);
    state.U = Unew.primitiveField();
    state.p = pMap.primitiveField();
    state.phiInternal = phiNew.primitiveField();
    forAll(state.phiBoundary, patchi)
    {
        state.phiBoundary[patchi] = phiNew.boundaryField()[patchi];
    }

    return state;
}


Foam::tmp<Foam::volScalarField>
Foam::thermalAdjointSimple::computePressureMap
(
    const surfaceScalarField& F2,
    const volVectorField& Uold,
    const volScalarField& pold,
    const volScalarField& rAtU,
    const word& namePrefix
)
{
    auto tpWork =
        tmp<volScalarField>::New
        (
            IOobject
            (
                namePrefix + "pWork",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pold
        );
    volScalarField& pWork = tpWork.ref();

    surfaceScalarField FWork
    (
        IOobject
        (
            namePrefix + "FWork",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F2
    );

    constrainPressure(pWork, Uold, FWork, rAtU);

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU, pWork) == fvc::div(FWork)
    );

    pEqn.setReference
    (
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );

    dictionary pSolver(pEqn.solverDict("p"));
    pSolver.set("relTol", scalar(0));
    pSolver.set("tolerance", scalar(1e-12));
    pEqn.solve(pSolver);

    return tpWork;
}


void Foam::thermalAdjointSimple::checkOneStepBetaMapSensitivity
(
    const vectorField& barUnewInt,
    const vectorField& UbarExactInt,
    const word& designVariablesName,
    const scalar dt
)
{
    if (!checkOneStepBetaMapSensitivity_)
    {
        return;
    }

    if (!mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        Info<< "ATC-T one-step beta-map check: no topOVars object" << endl;
        return;
    }

    const topOVariablesBase& vars =
        mesh_.lookupObject<topOVariablesBase>("topOVars");
    const topOZones& zones = vars.getTopOZones();
    const scalarField& V = mesh_.V();
    const incompressibleAdjointVars& adjointVars = getAdjointVars();
    const volVectorField& Ua = adjointVars.Ua();
    volScalarField& beta = const_cast<volScalarField&>(vars.beta());

    bitSet isDesign(mesh_.nCells(), false);
    if (zones.adjointPorousZoneIDs().empty())
    {
        isDesign.fill(true);
    }
    else
    {
        for (const label zoneID : zones.adjointPorousZoneIDs())
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.set(celli);
            }
        }
    }

    auto unsetZone = [&](const label zoneID)
    {
        if (zoneID != -1)
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.unset(celli);
            }
        }
    };

    for (const label zoneID : zones.fixedPorousZoneIDs())
    {
        unsetZone(zoneID);
    }
    for (const label zoneID : zones.fixedZeroPorousZoneIDs())
    {
        unsetZone(zoneID);
    }
    for (const label celli : zones.IOCells())
    {
        isDesign.unset(celli);
    }

    const label nReport = 15;
    labelList reportCells(nReport, -1);

    auto addReportCell = [&](const label celli)
    {
        if (celli < 0 || celli >= mesh_.nCells() || !isDesign.test(celli))
        {
            return;
        }

        forAll(reportCells, i)
        {
            if (reportCells[i] == celli)
            {
                return;
            }
        }

        forAll(reportCells, i)
        {
            if (reportCells[i] == -1)
            {
                reportCells[i] = celli;
                return;
            }
        }
    };

    addReportCell(707);
    addReportCell(708);
    addReportCell(709);

    const label firstRankedSlot = 3;

    forAll(barUnewInt, celli)
    {
        if (!isDesign.test(celli))
        {
            continue;
        }

        const scalar metric = mag(barUnewInt[celli]);
        for (label slot = firstRankedSlot; slot < reportCells.size(); ++slot)
        {
            if (reportCells[slot] == -1)
            {
                reportCells[slot] = celli;
                break;
            }

            if (metric > mag(barUnewInt[reportCells[slot]]))
            {
                for (label shift = reportCells.size() - 1; shift > slot; --shift)
                {
                    reportCells[shift] = reportCells[shift - 1];
                }
                reportCells[slot] = celli;
                break;
            }
        }
    }

    scalarField conductivitySens(mesh_.nCells(), Zero);
    {
        const volScalarField& T = TRef();
        const volScalarField& Ta = TaPtr_();
        const autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();
        const volScalarField Df(props_.DFluidField(T, "DFluidFieldMapCheck"));
        const volScalarField Ds(props_.DSolidField(T, "DSolidFieldMapCheck"));

        scalarField dJdDEff(mesh_.nCells(), Zero);
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const surfaceScalarField& w = mesh_.weights();
        const surfaceScalarField& deltaCoeffs = mesh_.nonOrthDeltaCoeffs();
        const surfaceScalarField& magSf = mesh_.magSf();

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];

            const scalar dJdDf =
                magSf[facei]*deltaCoeffs[facei]
               *(T[N] - T[P])*(Ta[N] - Ta[P]);

            dJdDEff[P] += w[facei]*dJdDf;
            dJdDEff[N] += (scalar(1) - w[facei])*dJdDf;
        }

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
                dJdDEff[fc[i]] -= Ta[fc[i]]*magSfb[i]*snGradT[i];
            }
        }

        conductivitySens =
            thermalSensScale_
           *dJdDEff/V
           *(
                Ds.primitiveField() - Df.primitiveField()
              - turbulence->nut()().primitiveField()/Prt_
            )
           *dt;

        vars.sourceTermSensitivities
        (
            conductivitySens,
            kInterpolation_(),
            scalar(1),
            designVariablesName,
            "beta"
        );
    }

    label nChecked = 0;
    label nSkipped = 0;

    Info<< "ATC-T one-step beta-map samples: "
        << "cell beta eps mapFD/dt mapFDVol/dt residualExactVol/dt "
        << "residualStockVol/dt conductivityVol/dt totalMapPlusCond/dt "
        << "ratioExact/map ratioStock/map"
        << endl;

    forAll(reportCells, i)
    {
        const label cellj = reportCells[i];
        if (cellj == -1 || !isDesign.test(cellj))
        {
            continue;
        }

        const scalar oldBeta = beta[cellj];
        const scalar room = min(oldBeta, scalar(1) - oldBeta);
        const scalar eps =
            min(oneStepBetaMapFDEps_, scalar(0.25)*max(SMALL, room));

        if (eps <= SMALL)
        {
            ++nSkipped;
            continue;
        }

        beta.primitiveFieldRef()[cellj] = oldBeta + eps;
        beta.correctBoundaryConditions();
        tmp<volVectorField> tUPlus =
            primalSimpleMapUAtFrozenState("ATCMapPlus");
        vectorField UPlus(tUPlus().primitiveField());
        tmp<vectorField> tRPlus = primalMomentumResidualLHS();
        vectorField RPlus(tRPlus());

        beta.primitiveFieldRef()[cellj] = oldBeta - eps;
        beta.correctBoundaryConditions();
        tmp<volVectorField> tUMinus =
            primalSimpleMapUAtFrozenState("ATCMapMinus");
        vectorField UMinus(tUMinus().primitiveField());
        tmp<vectorField> tRMinus = primalMomentumResidualLHS();
        vectorField RMinus(tRMinus());

        beta.primitiveFieldRef()[cellj] = oldBeta;
        beta.correctBoundaryConditions();

        scalar mapFD = Zero;
        scalar residualExact = Zero;
        scalar residualStock = Zero;

        forAll(V, celli)
        {
            const vector dUnew =
                (UPlus[celli] - UMinus[celli])/(2*eps);
            const vector dR =
                (RPlus[celli] - RMinus[celli])/(2*eps);
            const vector UbarExact = UbarExactInt[celli]/V[celli];

            mapFD += barUnewInt[celli] & dUnew;
            residualExact += UbarExact & dR;
            residualStock += Ua[celli] & dR;
        }

        const scalar mapFDVol = mapFD/V[cellj];
        const scalar residualExactVol = residualExact/V[cellj];
        const scalar residualStockVol = residualStock/V[cellj];
        const scalar condVol = conductivitySens[cellj];
        const scalar totalVol = mapFDVol + condVol;

        const scalar signedMapScale = mag(mapFD) > VSMALL ? mapFD : VSMALL;

        Info<< "ATC-T one-step beta-map sample: "
            << cellj << " "
            << oldBeta << " "
            << eps << " "
            << mapFD/dt << " "
            << mapFDVol/dt << " "
            << residualExactVol/dt << " "
            << residualStockVol/dt << " "
            << condVol/dt << " "
            << totalVol/dt << " "
            << residualExact/signedMapScale << " "
            << residualStock/signedMapScale
            << endl;

        ++nChecked;
    }

    beta.correctBoundaryConditions();

    Info<< "ATC-T one-step beta-map check: nChecked = " << nChecked
        << ", nSkipped = " << nSkipped
        << ", eps = " << oneStepBetaMapFDEps_
        << ". Diagnostic only; betaMult is unchanged."
        << endl;
}


void Foam::thermalAdjointSimple::checkStateMapTranspose
(
    const vectorField& barUnewInt,
    const volScalarField& rAtU
)
{
    if (!checkStateMapTranspose_)
    {
        return;
    }

    volVectorField& U = primalVars_.U();
    volScalarField& p =
        const_cast<volScalarField&>(primalVars_.p());

    const volVectorField UBase
    (
        IOobject
        (
            "ATCTStateMapUBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    const volScalarField pBase
    (
        IOobject
        (
            "ATCTStateMappBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );

    volVectorField dUold
    (
        IOobject
        (
            "ATCTStateMapdUold",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    dUold = dimensionedVector(dUold.dimensions(), vector::zero);

    forAll(dUold.primitiveFieldRef(), celli)
    {
        dUold.primitiveFieldRef()[celli] =
            vector
            (
                scalar((257*celli + 191) % 337)/scalar(337) - scalar(0.5),
                scalar((263*celli + 193) % 347)/scalar(347) - scalar(0.5),
                scalar((269*celli + 197) % 349)/scalar(349) - scalar(0.5)
            );
    }

    volScalarField dpold
    (
        IOobject
        (
            "ATCTStateMapdpold",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    dpold = dimensionedScalar(dpold.dimensions(), Zero);

    forAll(dpold.primitiveFieldRef(), celli)
    {
        dpold.primitiveFieldRef()[celli] =
            scalar((271*celli + 199) % 353)/scalar(353) - scalar(0.5);
    }

    const bool pressureNeedsReference = p.needReference();
    const label pRefCell = solverControl_().pRefCell();

    auto removeMean = [](scalarField& fld)
    {
        const label n =
            returnReduce(fld.size(), sumOp<label>());

        if (n > 0)
        {
            fld -= gSum(fld)/scalar(n);
        }
    };

    auto projectPressureDirection = [&](scalarField& dp)
    {
        if (!pressureNeedsReference)
        {
            return;
        }

        if (pRefCell >= 0 && pRefCell < mesh_.nCells())
        {
            dp[pRefCell] = Zero;
        }
        else
        {
            removeMean(dp);
        }
    };

    projectPressureDirection(dpold.primitiveFieldRef());

    forAll(dUold.boundaryFieldRef(), patchi)
    {
        if (!mesh_.boundary()[patchi].coupled())
        {
            dUold.boundaryFieldRef()[patchi] == vector::zero;
        }
    }
    forAll(dpold.boundaryFieldRef(), patchi)
    {
        if (!mesh_.boundary()[patchi].coupled())
        {
            dpold.boundaryFieldRef()[patchi] == Zero;
        }
    }

    scalar maxdU = VSMALL;
    scalar maxdp = VSMALL;
    forAll(dUold, celli)
    {
        maxdU = max(maxdU, mag(dUold[celli]));
        maxdp = max(maxdp, mag(dpold[celli]));
    }

    const scalar eps = scalar(1e-6);

    auto enforcePressureBoundaryState = [&]()
    {
        bool hasUpdateablePressureSnGrad = false;

        forAll(p.boundaryField(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            if
            (
                isA<updateablePatchTypes::updateableSnGrad>
                (
                    p.boundaryField()[patchi]
                )
            )
            {
                hasUpdateablePressureSnGrad = true;
            }
        }

        forAll(p.boundaryFieldRef(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            if (pBase.boundaryField()[patchi].fixesValue())
            {
                p.boundaryFieldRef()[patchi] ==
                    pBase.boundaryField()[patchi];
            }
        }

        if (!hasUpdateablePressureSnGrad)
        {
            p.correctBoundaryConditions();
        }
    };

    auto restoreState = [&]()
    {
        U.primitiveFieldRef() = UBase.primitiveField();
        p.primitiveFieldRef() = pBase.primitiveField();

        forAll(U.boundaryFieldRef(), patchi)
        {
            U.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
        }
        U.correctBoundaryConditions();

        forAll(p.boundaryFieldRef(), patchi)
        {
            p.boundaryFieldRef()[patchi] == pBase.boundaryField()[patchi];
        }
        enforcePressureBoundaryState();
    };

    auto setStatePerturbationComponents =
        [&](const scalar coeff, const bool perturbU, const bool perturbP)
    {
        U.primitiveFieldRef() = UBase.primitiveField();
        p.primitiveFieldRef() = pBase.primitiveField();

        if (perturbU)
        {
            U.primitiveFieldRef() += coeff*eps*dUold.primitiveField();
        }
        if (perturbP)
        {
            p.primitiveFieldRef() += coeff*eps*dpold.primitiveField();
        }

        forAll(U.boundaryFieldRef(), patchi)
        {
            U.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
        }
        U.correctBoundaryConditions();

        forAll(p.boundaryFieldRef(), patchi)
        {
            p.boundaryFieldRef()[patchi] == pBase.boundaryField()[patchi];
        }
        enforcePressureBoundaryState();
    };

    auto setStatePerturbation = [&](const scalar coeff)
    {
        setStatePerturbationComponents(coeff, true, true);
    };

    setStatePerturbation(scalar(1));
    SimpleMapState statePlus =
        primalSimpleMapStateAtFrozenState("ATCStateMapPlus");

    setStatePerturbation(scalar(-1));
    SimpleMapState stateMinus =
        primalSimpleMapStateAtFrozenState("ATCStateMapMinus");

    restoreState();

    volVectorField H1Plus
    (
        IOobject
        (
            "ATCTStateMapH1Plus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField H1Minus
    (
        IOobject
        (
            "ATCTStateMapH1Minus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField H0Plus
    (
        IOobject
        (
            "ATCTStateMapH0Plus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField H0Minus
    (
        IOobject
        (
            "ATCTStateMapH0Minus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    surfaceScalarField F2Plus
    (
        IOobject
        (
            "ATCTStateMapF2Plus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.phi()
    );
    surfaceScalarField F2Minus
    (
        IOobject
        (
            "ATCTStateMapF2Minus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.phi()
    );
    surfaceScalarField F1Plus
    (
        IOobject
        (
            "ATCTStateMapF1Plus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.phi()
    );
    surfaceScalarField F1Minus
    (
        IOobject
        (
            "ATCTStateMapF1Minus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primalVars_.phi()
    );
    volVectorField UrawPlus
    (
        IOobject
        (
            "ATCTStateMapUrawPlus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField UrawMinus
    (
        IOobject
        (
            "ATCTStateMapUrawMinus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField UbcPlus
    (
        IOobject
        (
            "ATCTStateMapUbcPlus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField UbcMinus
    (
        IOobject
        (
            "ATCTStateMapUbcMinus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField UstagePlus
    (
        IOobject
        (
            "ATCTStateMapUstagePlus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField UstageMinus
    (
        IOobject
        (
            "ATCTStateMapUstageMinus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volScalarField pStagePlus
    (
        IOobject
        (
            "ATCTStateMappStagePlus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    volScalarField pStageMinus
    (
        IOobject
        (
            "ATCTStateMappStageMinus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    volScalarField rAtUStagePlus
    (
        IOobject
        (
            "ATCTStateMaprAtUStagePlus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rAtU
    );
    volScalarField rAtUStageMinus
    (
        IOobject
        (
            "ATCTStateMaprAtUStageMinus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rAtU
    );
    volVectorField UrawBaseStage
    (
        IOobject
        (
            "ATCTStateMapUrawBaseStage",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField UbcBaseStage
    (
        IOobject
        (
            "ATCTStateMapUbcBaseStage",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volVectorField UfinalBaseStage
    (
        IOobject
        (
            "ATCTStateMapUfinalBaseStage",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    volScalarField pStageBase
    (
        IOobject
        (
            "ATCTStateMappStageBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    volScalarField rAtUStageBase
    (
        IOobject
        (
            "ATCTStateMaprAtUStageBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rAtU
    );

    auto computePrePressureMap =
        [&]
        (
            const word& namePrefix,
            volVectorField& H0,
            surfaceScalarField& F1,
            volVectorField& H1,
            surfaceScalarField& F2
        )
    {
        autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();
        fv::options& fvOptions(fv::options::New(this->mesh_));

        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(primalVars_.phi(), U)
          + turbulence->divDevReff(U)
         ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();

        UEqn.relax();
        fvOptions.constrain(UEqn);

        volScalarField rAU(1.0/UEqn.A());
        volScalarField pMap
        (
            IOobject
            (
                namePrefix + "pMap",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p
        );

        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, pMap));
        surfaceScalarField phiHbyA(namePrefix + "phiHbyA", fvc::flux(HbyA));
        adjustPhi(phiHbyA, U, pMap);

        H0 = HbyA;
        F1 = phiHbyA;

        if (solverControl_().consistent())
        {
            tmp<volScalarField> trAtU =
                1.0/(1.0/rAU - UEqn.H1());
            phiHbyA +=
                fvc::interpolate(trAtU() - rAU)
               *fvc::snGrad(pMap)
               *mesh_.magSf();
            HbyA -= (rAU - trAtU())*fvc::grad(pMap);
        }

        H1 = HbyA;
        F2 = phiHbyA;
    };

    auto computeFullMapStages =
        [&]
        (
            const word& namePrefix,
            volVectorField& Uraw,
            volVectorField& Ubc,
            volVectorField& Ustage,
            volScalarField& pStage,
            volScalarField& rAtUStage
        )
    {
        autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();
        fv::options& fvOptions(fv::options::New(this->mesh_));

        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(primalVars_.phi(), U)
          + turbulence->divDevReff(U)
         ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();

        UEqn.relax();
        fvOptions.constrain(UEqn);

        volScalarField rAU(1.0/UEqn.A());
        volScalarField pMap
        (
            IOobject
            (
                namePrefix + "pMap",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p
        );
        pMap.storePrevIter();

        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, pMap));
        surfaceScalarField phiHbyA(namePrefix + "phiHbyA", fvc::flux(HbyA));
        adjustPhi(phiHbyA, U, pMap);

        tmp<volScalarField> trAtU(rAU);

        if (solverControl_().consistent())
        {
            trAtU = 1.0/(1.0/rAU - UEqn.H1());
            phiHbyA +=
                fvc::interpolate(trAtU() - rAU)
               *fvc::snGrad(pMap)
               *mesh_.magSf();
            HbyA -= (rAU - trAtU())*fvc::grad(pMap);
        }

        const volScalarField& rAtULocal = trAtU();
        rAtUStage = rAtULocal;
        constrainPressure(pMap, U, phiHbyA, rAtULocal);

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtULocal, pMap) == fvc::div(phiHbyA)
        );
        pEqn.setReference
        (
            solverControl_().pRefCell(),
            solverControl_().pRefValue()
        );
        dictionary pSolver(pEqn.solverDict("p"));
        pSolver.set("relTol", scalar(0));
        pSolver.set("tolerance", scalar(1e-12));
        pEqn.solve(pSolver);

        scalar pRelaxCoeff = scalar(1);
        word pRelaxName = p.name();
        if (p.mesh().data().isFinalIteration())
        {
            pRelaxName += "Final";
        }
        p.mesh().relaxField(pRelaxName, pRelaxCoeff);
        pMap.relax(pRelaxCoeff);

        Uraw = HbyA - rAtULocal*fvc::grad(pMap);
        Ubc = Uraw;
        Ubc.correctBoundaryConditions();
        Ustage = Ubc;
        fvOptions.correct(Ustage);
        pStage = pMap;
        rAtUStage.correctBoundaryConditions();
    };

    auto computeConsistentQ =
        [&](const word& namePrefix)
    {
        autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();
        fv::options& fvOptions(fv::options::New(this->mesh_));

        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(primalVars_.phi(), U)
          + turbulence->divDevReff(U)
         ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();

        UEqn.relax();
        fvOptions.constrain(UEqn);

        volScalarField rAU(1.0/UEqn.A());
        tmp<volScalarField> trAtU =
            1.0/(1.0/rAU - UEqn.H1());

        return tmp<volScalarField>::New
        (
            IOobject
            (
                namePrefix + "q",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            trAtU() - rAU
        );
    };

    setStatePerturbation(scalar(1));
    computePrePressureMap
    (
        "ATCStateMapPrePlus",
        H0Plus,
        F1Plus,
        H1Plus,
        F2Plus
    );

    setStatePerturbation(scalar(-1));
    computePrePressureMap
    (
        "ATCStateMapPreMinus",
        H0Minus,
        F1Minus,
        H1Minus,
        F2Minus
    );

    setStatePerturbation(scalar(1));
    computeFullMapStages
    (
        "ATCStateMapStagePlus",
        UrawPlus,
        UbcPlus,
        UstagePlus,
        pStagePlus,
        rAtUStagePlus
    );

    setStatePerturbation(scalar(-1));
    computeFullMapStages
    (
        "ATCStateMapStageMinus",
        UrawMinus,
        UbcMinus,
        UstageMinus,
        pStageMinus,
        rAtUStageMinus
    );

    setStatePerturbation(scalar(0));
    computeFullMapStages
    (
        "ATCStateMapStageBase",
        UrawBaseStage,
        UbcBaseStage,
        UfinalBaseStage,
        pStageBase,
        rAtUStageBase
    );

    restoreState();

    vectorField dUnew(mesh_.nCells(), vector::zero);
    scalarField dpnew(mesh_.nCells(), Zero);
    vectorField dUrawStage(mesh_.nCells(), vector::zero);
    vectorField dUbcStage(mesh_.nCells(), vector::zero);
    vectorField dUstage(mesh_.nCells(), vector::zero);
    scalarField dpstage(mesh_.nCells(), Zero);
    vectorField dH0Map(mesh_.nCells(), vector::zero);
    vectorField dH1Map(mesh_.nCells(), vector::zero);
    surfaceScalarField dF1Map
    (
        IOobject
        (
            "ATCTStateMapdF1Map",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F1Plus
    );
    surfaceScalarField dF2Map
    (
        IOobject
        (
            "ATCTStateMapdF2Map",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F2Plus
    );

    forAll(dUnew, celli)
    {
        dUnew[celli] =
            (statePlus.U[celli] - stateMinus.U[celli])/(2*eps);
        dpnew[celli] =
            (statePlus.p[celli] - stateMinus.p[celli])/(2*eps);
        dUrawStage[celli] =
            (UrawPlus[celli] - UrawMinus[celli])/(2*eps);
        dUbcStage[celli] =
            (UbcPlus[celli] - UbcMinus[celli])/(2*eps);
        dUstage[celli] =
            (UstagePlus[celli] - UstageMinus[celli])/(2*eps);
        dpstage[celli] =
            (pStagePlus[celli] - pStageMinus[celli])/(2*eps);
        dH0Map[celli] =
            (H0Plus[celli] - H0Minus[celli])/(2*eps);
        dH1Map[celli] =
            (H1Plus[celli] - H1Minus[celli])/(2*eps);
    }
    dF1Map.primitiveFieldRef() =
        (F1Plus.primitiveField() - F1Minus.primitiveField())/(2*eps);
    dF2Map.primitiveFieldRef() =
        (F2Plus.primitiveField() - F2Minus.primitiveField())/(2*eps);
    forAll(dF1Map.boundaryFieldRef(), patchi)
    {
        dF1Map.boundaryFieldRef()[patchi] ==
            (
                F1Plus.boundaryField()[patchi]
              - F1Minus.boundaryField()[patchi]
            )/(2*eps);
    }
    forAll(dF2Map.boundaryFieldRef(), patchi)
    {
        dF2Map.boundaryFieldRef()[patchi] ==
            (
                F2Plus.boundaryField()[patchi]
              - F2Minus.boundaryField()[patchi]
            )/(2*eps);
    }

    scalarField dqMap(mesh_.nCells(), Zero);
    if (solverControl_().consistent())
    {
        setStatePerturbation(scalar(1));
        tmp<volScalarField> tqPlus =
            computeConsistentQ("ATCStateMapQPlus");

        setStatePerturbation(scalar(-1));
        tmp<volScalarField> tqMinus =
            computeConsistentQ("ATCStateMapQMinus");

        restoreState();

        dqMap =
            (tqPlus().primitiveField() - tqMinus().primitiveField())
           /(2*eps);
    }

    vectorField dH0MapU(mesh_.nCells(), vector::zero);
    vectorField dH1MapU(mesh_.nCells(), vector::zero);
    vectorField dH0MapP(mesh_.nCells(), vector::zero);
    vectorField dH1MapP(mesh_.nCells(), vector::zero);

    surfaceScalarField dF1MapU
    (
        IOobject
        (
            "ATCTStateMapdF1MapU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dF1Map
    );
    surfaceScalarField dF2MapU
    (
        IOobject
        (
            "ATCTStateMapdF2MapU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dF2Map
    );
    surfaceScalarField dF1MapP
    (
        IOobject
        (
            "ATCTStateMapdF1MapP",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dF1Map
    );
    surfaceScalarField dF2MapP
    (
        IOobject
        (
            "ATCTStateMapdF2MapP",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dF2Map
    );

    scalarField dqMapU(mesh_.nCells(), Zero);
    scalarField dqMapP(mesh_.nCells(), Zero);

    auto computePrePressureDelta =
        [&]
        (
            const word& suffix,
            const bool perturbU,
            const bool perturbP,
            vectorField& dH0,
            surfaceScalarField& dF1,
            vectorField& dH1,
            surfaceScalarField& dF2,
            scalarField& dq
        )
    {
        setStatePerturbationComponents(scalar(1), perturbU, perturbP);
        computePrePressureMap
        (
            word("ATCStateMap") + suffix + "Plus",
            H0Plus,
            F1Plus,
            H1Plus,
            F2Plus
        );

        scalarField qPlus(mesh_.nCells(), Zero);
        if (solverControl_().consistent())
        {
            tmp<volScalarField> tqPlus =
                computeConsistentQ(word("ATCStateMap") + suffix + "QPlus");
            qPlus = tqPlus().primitiveField();
        }

        setStatePerturbationComponents(scalar(-1), perturbU, perturbP);
        computePrePressureMap
        (
            word("ATCStateMap") + suffix + "Minus",
            H0Minus,
            F1Minus,
            H1Minus,
            F2Minus
        );

        scalarField qMinus(mesh_.nCells(), Zero);
        if (solverControl_().consistent())
        {
            tmp<volScalarField> tqMinus =
                computeConsistentQ(word("ATCStateMap") + suffix + "QMinus");
            qMinus = tqMinus().primitiveField();
        }

        restoreState();

        forAll(dH0, celli)
        {
            dH0[celli] = (H0Plus[celli] - H0Minus[celli])/(2*eps);
            dH1[celli] = (H1Plus[celli] - H1Minus[celli])/(2*eps);
        }

        dF1.primitiveFieldRef() =
            (F1Plus.primitiveField() - F1Minus.primitiveField())/(2*eps);
        dF2.primitiveFieldRef() =
            (F2Plus.primitiveField() - F2Minus.primitiveField())/(2*eps);

        forAll(dF1.boundaryFieldRef(), patchi)
        {
            dF1.boundaryFieldRef()[patchi] ==
                (
                    F1Plus.boundaryField()[patchi]
                  - F1Minus.boundaryField()[patchi]
                )/(2*eps);
        }
        forAll(dF2.boundaryFieldRef(), patchi)
        {
            dF2.boundaryFieldRef()[patchi] ==
                (
                    F2Plus.boundaryField()[patchi]
                  - F2Minus.boundaryField()[patchi]
                )/(2*eps);
        }

        if (solverControl_().consistent())
        {
            dq = (qPlus - qMinus)/(2*eps);
        }
        else
        {
            dq = Zero;
        }
    };

    computePrePressureDelta
    (
        "UOnly",
        true,
        false,
        dH0MapU,
        dF1MapU,
        dH1MapU,
        dF2MapU,
        dqMapU
    );

    computePrePressureDelta
    (
        "POnly",
        false,
        true,
        dH0MapP,
        dF1MapP,
        dH1MapP,
        dF2MapP,
        dqMapP
    );

    auto projectPressureSeed = [&](SimpleMapSeed& seed)
    {
        if (!pressureNeedsReference)
        {
            return;
        }

        if (pRefCell >= 0 && pRefCell < mesh_.nCells())
        {
            seed.barp[pRefCell] = Zero;
        }
        else
        {
            removeMean(seed.barp);
        }
    };

    projectPressureDirection(dpnew);
    projectPressureDirection(dpstage);

    volScalarField dpStageInternalBC
    (
        IOobject
        (
            "ATCTStateMapdpStageInternalBC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pStageBase
    );
    dpStageInternalBC.primitiveFieldRef() = dpstage;
    dpStageInternalBC.correctBoundaryConditions();

    volVectorField dUrawInternalBC
    (
        IOobject
        (
            "ATCTStateMapdUrawInternalBC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    dUrawInternalBC.primitiveFieldRef() = dH1Map;
    dUrawInternalBC -= rAtU*fvc::grad(dpStageInternalBC);

    surfaceScalarField pFaceStagePlus
    (
        IOobject
        (
            "ATCTStateMapPFaceStagePlus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(pStagePlus)
    );
    surfaceScalarField pFaceStageMinus
    (
        IOobject
        (
            "ATCTStateMapPFaceStageMinus",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(pStageMinus)
    );
    surfaceScalarField pFaceStageInternalBC
    (
        IOobject
        (
            "ATCTStateMapPFaceStageInternalBC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(dpStageInternalBC)
    );

    surfaceScalarField dpFaceStage
    (
        IOobject
        (
            "ATCTStateMapDPFaceStage",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pFaceStagePlus
    );
    dpFaceStage =
        (pFaceStagePlus - pFaceStageMinus)/(2*eps);

    volVectorField dGradPStageActual
    (
        IOobject
        (
            "ATCTStateMapDGradPStageActual",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(pStagePlus)
    );
    dGradPStageActual -= fvc::grad(pStageMinus);
    dGradPStageActual *= scalar(1)/(2*eps);

    volVectorField dGradPStageInternalBC
    (
        IOobject
        (
            "ATCTStateMapDGradPStageInternalBC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(dpStageInternalBC)
    );

    auto seedNorms = [&](const SimpleMapSeed& seed, scalar& Unorm, scalar& pnorm)
    {
        Unorm = VSMALL;
        pnorm = VSMALL;
        forAll(seed.barU, celli)
        {
            Unorm += magSqr(seed.barU[celli]);
            pnorm += sqr(seed.barp[celli]);
        }
        Unorm = sqrt(Unorm);
        pnorm = sqrt(pnorm);
    };

    auto runCheck = [&](const word& checkName, const SimpleMapSeed& rawSeed)
    {
        SimpleMapSeed seedNew(rawSeed);
        projectPressureSeed(seedNew);

        SimpleMapSeed seedOld = reverseOneSimpleMapSeed(seedNew, rAtU);
        projectPressureSeed(seedOld);

        scalar knownBoundaryLhs = Zero;
        scalar knownBoundaryRhsF2 = Zero;
        scalar knownBoundaryRhsU = Zero;
        scalar knownBoundaryRhsRAtU = Zero;
        scalar knownBoundaryRhsOwner = Zero;
        scalar knownBoundaryFaceLhs = Zero;
        scalar knownBoundaryGradDelta = Zero;

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            const bool updateablePressureSnGrad =
                isA<updateablePatchTypes::updateableSnGrad>
                (
                    pStageBase.boundaryField()[patchi]
                );

            const fvPatch& patch = mesh_.boundary()[patchi];
            const labelUList& faceCells = patch.faceCells();
            const fvsPatchVectorField& SfPatch =
                mesh_.Sf().boundaryField()[patchi];
            const fvsPatchScalarField& magSfPatch =
                mesh_.magSf().boundaryField()[patchi];
            const fvsPatchScalarField& dF2Patch =
                dF2Map.boundaryField()[patchi];
            const fvsPatchScalarField& dpFacePatch =
                dpFaceStage.boundaryField()[patchi];
            const fvsPatchScalarField& pFaceInternalPatch =
                pFaceStageInternalBC.boundaryField()[patchi];
            const fvsPatchScalarField& F2PlusPatch =
                F2Plus.boundaryField()[patchi];
            const fvsPatchScalarField& F2MinusPatch =
                F2Minus.boundaryField()[patchi];
            const fvPatchScalarField& pPlusPatch =
                pStagePlus.boundaryField()[patchi];
            const fvPatchScalarField& pMinusPatch =
                pStageMinus.boundaryField()[patchi];
            const fvPatchVectorField& UBasePatch =
                UBase.boundaryField()[patchi];
            const fvPatchVectorField& dUPatch =
                dUold.boundaryField()[patchi];
            const fvPatchScalarField& rAtUBasePatch =
                rAtUStageBase.boundaryField()[patchi];
            const fvPatchScalarField& rAtUPlusPatch =
                rAtUStagePlus.boundaryField()[patchi];
            const fvPatchScalarField& rAtUMinusPatch =
                rAtUStageMinus.boundaryField()[patchi];
            const scalarField& deltaCoeffs = patch.deltaCoeffs();

            tmp<scalarField> tpValueInternalCoeffs =
                pStageBase.boundaryField()[patchi].valueInternalCoeffs
                (
                    patch.weights()
                );
            const scalarField& pValueInternalCoeffs =
                tpValueInternalCoeffs();

            scalar patchLhs = Zero;
            scalar patchRhsF2 = Zero;
            scalar patchRhsU = Zero;
            scalar patchRhsRAtU = Zero;
            scalar patchRhsOwner = Zero;
            scalar patchFaceLhs = Zero;
            scalar patchDpBoundaryActual = Zero;
            scalar patchDpBoundaryOwner = Zero;
            scalar patchDpBoundaryKnown = Zero;

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                const scalar rBoundary = rAtUBasePatch[facei];
                const scalar denom =
                    deltaCoeffs[facei]*magSfPatch[facei]*rBoundary;

                if (mag(denom) <= VSMALL)
                {
                    continue;
                }

                const vector sourceCell =
                    seedNew.barU[celli]/mesh_.V()[celli];
                const scalar barPb =
                   -(SfPatch[facei]
                   & (rAtUStageBase[celli]*sourceCell));

                const scalar dpBoundaryActual =
                    (pPlusPatch[facei] - pMinusPatch[facei])/(2*eps);
                const scalar dpBoundaryKnown =
                    dpBoundaryActual
                  - pValueInternalCoeffs[facei]*dpstage[celli];
                const scalar dpBoundaryOwner =
                    pValueInternalCoeffs[facei]*dpstage[celli];
                const scalar dpFaceKnown =
                    dpFacePatch[facei] - pFaceInternalPatch[facei];

                const scalar F2Base =
                    scalar(0.5)
                   *(F2PlusPatch[facei] + F2MinusPatch[facei]);
                const scalar dRAtUBoundary =
                    (rAtUPlusPatch[facei] - rAtUMinusPatch[facei])
                   /(2*eps);

                const scalar barF2Known = barPb/denom;
                const vector barUbKnown = -barPb*SfPatch[facei]/denom;
                const scalar barRAtUKnown =
                   -barPb
                   *(F2Base - (SfPatch[facei] & UBasePatch[facei]))
                   /(deltaCoeffs[facei]*magSfPatch[facei]*sqr(rBoundary));

                patchLhs += barPb*dpBoundaryKnown;
                patchFaceLhs += barPb*dpFaceKnown;
                patchDpBoundaryActual += barPb*dpBoundaryActual;
                patchDpBoundaryOwner += barPb*dpBoundaryOwner;
                patchDpBoundaryKnown += barPb*dpBoundaryKnown;

                if (updateablePressureSnGrad)
                {
                    patchRhsF2 += barF2Known*dF2Patch[facei];
                    patchRhsU += barUbKnown & dUPatch[facei];
                    patchRhsRAtU += barRAtUKnown*dRAtUBoundary;
                }
            }

            knownBoundaryLhs += patchLhs;
            knownBoundaryRhsF2 += patchRhsF2;
            knownBoundaryRhsU += patchRhsU;
            knownBoundaryRhsRAtU += patchRhsRAtU;
            knownBoundaryRhsOwner += patchRhsOwner;
            knownBoundaryFaceLhs += patchFaceLhs;
            const scalar patchRhs =
                patchRhsF2 + patchRhsU + patchRhsRAtU + patchRhsOwner;

            const scalar patchScale =
                max(max(mag(patchLhs), mag(patchRhs)), VSMALL);

            Info<< "ATC-T pressure-known-boundary transpose check: "
                << checkName
                << ", patch = " << patch.name()
                << ", pType = " << pStageBase.boundaryField()[patchi].type()
                << ", updateable = " << updateablePressureSnGrad
                << ", lhsKnown = " << patchLhs
                << ", rhsF2 = " << patchRhsF2
                << ", rhsU = " << patchRhsU
                << ", rhsRAtU = " << patchRhsRAtU
                << ", rhsOwner = " << patchRhsOwner
                << ", lhsFaceKnown = " << patchFaceLhs
                << ", dpBoundaryActualSeeded = "
                << patchDpBoundaryActual
                << ", dpBoundaryOwnerSeeded = "
                << patchDpBoundaryOwner
                << ", dpBoundaryKnownSeeded = "
                << patchDpBoundaryKnown
                << ", rel = " << mag(patchLhs - patchRhs)/patchScale
                << endl;
        }

        const scalar knownBoundaryRhs =
            knownBoundaryRhsF2
          + knownBoundaryRhsU
          + knownBoundaryRhsRAtU
          + knownBoundaryRhsOwner;
        const scalar knownBoundaryScale =
            max(max(mag(knownBoundaryLhs), mag(knownBoundaryRhs)), VSMALL);
        const scalar knownBoundaryFaceScale =
            max
            (
                max(mag(knownBoundaryFaceLhs), mag(knownBoundaryGradDelta)),
                VSMALL
            );

        scalar lhsU = Zero;
        scalar lhsP = Zero;
        scalar lhsUrawStage = Zero;
        scalar lhsUbcStage = Zero;
        scalar lhsUstage = Zero;
        scalar lhsUrawInternalBC = Zero;
        scalar lhsPstage = Zero;
        scalar rhsU = Zero;
        scalar rhsP = Zero;
        scalar barpSum = Zero;

        forAll(dUnew, celli)
        {
            lhsU += seedNew.barU[celli] & dUnew[celli];
            lhsP += seedNew.barp[celli]*dpnew[celli];
            lhsUrawStage += seedNew.barU[celli] & dUrawStage[celli];
            lhsUbcStage += seedNew.barU[celli] & dUbcStage[celli];
            lhsUstage += seedNew.barU[celli] & dUstage[celli];
            lhsUrawInternalBC +=
                seedNew.barU[celli] & dUrawInternalBC[celli];
            lhsPstage += seedNew.barp[celli]*dpstage[celli];
            rhsU += seedOld.barU[celli] & dUold[celli];
            rhsP += seedOld.barp[celli]*dpold[celli];
            barpSum += seedNew.barp[celli];
            knownBoundaryGradDelta +=
                seedNew.barU[celli]
              & (
                    -rAtU[celli]
                   *(
                        dGradPStageActual[celli]
                      - dGradPStageInternalBC[celli]
                    )
                );
        }

        const scalar lhs = lhsU + lhsP;
        const scalar lhsRawStage = lhsUrawStage + lhsPstage;
        const scalar lhsBcStage = lhsUbcStage + lhsPstage;
        const scalar lhsFinalStage = lhsUstage + lhsPstage;
        const scalar lhsRawInternalBC = lhsUrawInternalBC + lhsPstage;
        const scalar rhs = rhsU + rhsP;
        const scalar scale = max(max(mag(lhs), mag(rhs)), VSMALL);
        const scalar tapeScale =
            max(max(mag(lhs), mag(lhsFinalStage)), VSMALL);

        scalar seedUNorm = Zero;
        scalar seedPNorm = Zero;
        seedNorms(seedNew, seedUNorm, seedPNorm);

        const scalar pRefSeed =
            (pRefCell >= 0 && pRefCell < mesh_.nCells())
          ? seedNew.barp[pRefCell]
          : Zero;

        Info<< "ATC-T state-map transpose check: "
            << checkName
            << " lhsU = " << lhsU
            << ", lhsP = " << lhsP
            << ", lhsRawStage = " << lhsRawStage
            << ", lhsRawInternalBC = " << lhsRawInternalBC
            << ", pressureBoundaryKnownDelta = "
            << lhsRawStage - lhsRawInternalBC
            << ", knownBoundaryLhs = " << knownBoundaryLhs
            << ", knownBoundaryRhs = " << knownBoundaryRhs
            << ", knownBoundaryRhsF2 = " << knownBoundaryRhsF2
            << ", knownBoundaryRhsU = " << knownBoundaryRhsU
            << ", knownBoundaryRhsRAtU = " << knownBoundaryRhsRAtU
            << ", knownBoundaryRhsOwner = " << knownBoundaryRhsOwner
            << ", knownBoundaryFaceLhs = " << knownBoundaryFaceLhs
            << ", knownBoundaryGradDelta = " << knownBoundaryGradDelta
            << ", knownBoundaryFaceRel = "
            << mag(knownBoundaryFaceLhs - knownBoundaryGradDelta)
              /knownBoundaryFaceScale
            << ", knownBoundaryRel = "
            << mag(knownBoundaryLhs - knownBoundaryRhs)
              /knownBoundaryScale
            << ", lhsBcStage = " << lhsBcStage
            << ", lhsFinalStage = " << lhsFinalStage
            << ", postBcDelta = " << lhsBcStage - lhsRawStage
            << ", postFvOptionsDelta = " << lhsFinalStage - lhsBcStage
            << ", finalTapeRel = " << mag(lhs - lhsFinalStage)/tapeScale
            << ", rhsU = " << rhsU
            << ", rhsP = " << rhsP
            << ", lhs = " << lhs
            << ", rhs = " << rhs
            << ", rel = " << mag(lhs - rhs)/scale
            << ", seedUNorm = " << seedUNorm
            << ", seedPNorm = " << seedPNorm
            << ", sumBarp = " << barpSum
            << ", pRefSeed = " << pRefSeed
            << endl;
    };

    SimpleMapSeed randomSeed(mesh_);
    forAll(randomSeed.barU, celli)
    {
        randomSeed.barU[celli] =
            vector
            (
                scalar((277*celli + 211) % 359)/scalar(359) - scalar(0.5),
                scalar((281*celli + 223) % 367)/scalar(367) - scalar(0.5),
                scalar((283*celli + 227) % 373)/scalar(373) - scalar(0.5)
            );
        randomSeed.barp[celli] =
            scalar((293*celli + 229) % 379)/scalar(379) - scalar(0.5);
    }
    projectPressureSeed(randomSeed);

    SimpleMapSeed initialSeed(mesh_);
    initialSeed.barU = barUnewInt;
    projectPressureSeed(initialSeed);

    SimpleMapSeed seedAfterOne = reverseOneSimpleMapSeed(initialSeed, rAtU);
    projectPressureSeed(seedAfterOne);

    SimpleMapSeed seedAfterFour(seedAfterOne);
    for (label step = 1; step < 4; ++step)
    {
        seedAfterFour = reverseOneSimpleMapSeed(seedAfterFour, rAtU);
        projectPressureSeed(seedAfterFour);
    }

    auto pressureDot = [&](const scalarField& a, const scalarField& b)
    {
        scalar value = Zero;
        forAll(a, celli)
        {
            value += a[celli]*b[celli];
        }
        return value;
    };

    auto vectorDot = [&](const vectorField& a, const vectorField& b)
    {
        scalar value = Zero;
        forAll(a, celli)
        {
            value += a[celli] & b[celli];
        }
        return value;
    };

    auto faceDot =
        [&](const surfaceScalarField& a, const surfaceScalarField& b)
    {
        scalar value = Zero;
        forAll(a.primitiveField(), facei)
        {
            value += a[facei]*b[facei];
        }

        forAll(mesh_.boundary(), patchi)
        {
            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            const fvsPatchScalarField& ap = a.boundaryField()[patchi];
            const fvsPatchScalarField& bp = b.boundaryField()[patchi];
            forAll(ap, facei)
            {
                value += ap[facei]*bp[facei];
            }
        }

        return value;
    };

    volVectorField dH1
    (
        IOobject
        (
            "ATCTStateMapdH1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    dH1 = dimensionedVector(dH1.dimensions(), vector::zero);

    forAll(dH1.primitiveFieldRef(), celli)
    {
        dH1.primitiveFieldRef()[celli] =
            vector
            (
                scalar((359*celli + 281) % 431)/scalar(431) - scalar(0.5),
                scalar((367*celli + 283) % 433)/scalar(433) - scalar(0.5),
                scalar((373*celli + 293) % 439)/scalar(439) - scalar(0.5)
            );
    }
    forAll(dH1.boundaryFieldRef(), patchi)
    {
        if (!mesh_.boundary()[patchi].coupled())
        {
            dH1.boundaryFieldRef()[patchi] == vector::zero;
        }
    }

    volScalarField dpSolve
    (
        IOobject
        (
            "ATCTStateMapdpSolve",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    dpSolve = dimensionedScalar(dpSolve.dimensions(), Zero);

    volScalarField dpRelax
    (
        IOobject
        (
            "ATCTStateMapdpRelax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    dpRelax = dimensionedScalar(dpRelax.dimensions(), Zero);

    forAll(dpSolve.primitiveFieldRef(), celli)
    {
        dpSolve.primitiveFieldRef()[celli] =
            scalar((379*celli + 307) % 443)/scalar(443) - scalar(0.5);
        dpRelax.primitiveFieldRef()[celli] =
            scalar((383*celli + 311) % 449)/scalar(449) - scalar(0.5);
    }
    projectPressureDirection(dpSolve.primitiveFieldRef());
    projectPressureDirection(dpRelax.primitiveFieldRef());
    forAll(dpSolve.boundaryFieldRef(), patchi)
    {
        if (!mesh_.boundary()[patchi].coupled())
        {
            dpSolve.boundaryFieldRef()[patchi] == Zero;
            dpRelax.boundaryFieldRef()[patchi] == Zero;
        }
    }
    dpRelax.correctBoundaryConditions();

    scalar pRelaxCoeff = scalar(1);
    word pRelaxName = p.name();
    if (p.mesh().data().isFinalIteration())
    {
        pRelaxName += "Final";
    }
    p.mesh().relaxField(pRelaxName, pRelaxCoeff);

    scalarField dpRelaxFromSolve(dpSolve.primitiveField());
    dpRelaxFromSolve *= pRelaxCoeff;
    dpRelaxFromSolve += (scalar(1) - pRelaxCoeff)*dpold.primitiveField();
    projectPressureDirection(dpRelaxFromSolve);

    auto finalCorrectionPressureSeed =
        [&](const vectorField& barU)
    {
        vectorField source(mesh_.nCells(), vector::zero);
        forAll(source, celli)
        {
            source[celli] = barU[celli]/mesh_.V()[celli];
        }

        scalarField barp(mesh_.nCells(), Zero);
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const surfaceVectorField& Sf = mesh_.Sf();
        const surfaceScalarField& w = mesh_.weights();

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];
            const scalar SfSource =
                Sf[facei] & (rAtU[P]*source[P] - rAtU[N]*source[N]);

            barp[P] -= w[facei]*SfSource;
            barp[N] -= (scalar(1) - w[facei])*SfSource;
        }

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            const fvPatch& patch = mesh_.boundary()[patchi];
            const labelUList& faceCells = patch.faceCells();
            const fvsPatchVectorField& SfPatch =
                Sf.boundaryField()[patchi];

            tmp<scalarField> tValueInternalCoeffs =
                pStageBase.boundaryField()[patchi].valueInternalCoeffs
                (
                    patch.weights()
                );
            const scalarField& valueInternalCoeffs =
                tValueInternalCoeffs();

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                barp[celli] -=
                    valueInternalCoeffs[facei]
                   *(SfPatch[facei] & (rAtU[celli]*source[celli]));
            }
        }

        return barp;
    };

    auto pressureSolveReverseF2 =
        [&](const scalarField& rawBarpSolve, surfaceScalarField& barF2)
    {
        scalarField barpSolve(rawBarpSolve);
        projectPressureDirection(barpSolve);

        volScalarField mu
        (
            IOobject
            (
                "paATCTPrePressureSeed",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p
        );
        mu = dimensionedScalar(mu.dimensions(), Zero);

        fvScalarMatrix muEqn(fvm::laplacian(rAtU, mu));
        muEqn.source() = barpSolve;
        muEqn.setReference(solverControl_().pRefCell(), Zero);
        dictionary muSolver(muEqn.solverDict("p"));
        muSolver.set("relTol", scalar(0));
        muSolver.set("tolerance", scalar(1e-12));
        muEqn.solve(muSolver);

        barF2 = dimensionedScalar(barF2.dimensions(), Zero);

        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];
            barF2.primitiveFieldRef()[facei] = mu[P] - mu[N];
        }

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            fvsPatchScalarField& barF2p =
                barF2.boundaryFieldRef()[patchi];
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(barF2p, facei)
            {
                barF2p[facei] = mu[faceCells[facei]];
            }
        }
    };

    auto consistentOldPressureSeed =
        [&]
        (
            const vectorField& barH1,
            const surfaceScalarField& barF2,
            scalarField& pOldSeed,
            scalarField& qSeed,
            scalarField& pOldSeedH,
            scalarField& pOldSeedF,
            scalarField& pOldSeedBoundaryH,
            scalarField& pOldSeedBoundaryF
        )
    {
        pOldSeed = Zero;
        qSeed = Zero;
        pOldSeedH = Zero;
        pOldSeedF = Zero;
        pOldSeedBoundaryH = Zero;
        pOldSeedBoundaryF = Zero;

        if (!solverControl_().consistent())
        {
            return;
        }

        autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();
        fv::options& fvOptions(fv::options::New(this->mesh_));

        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(primalVars_.phi(), U)
          + turbulence->divDevReff(U)
         ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();

        UEqn.relax();
        fvOptions.constrain(UEqn);

        volScalarField rAU(1.0/UEqn.A());
        tmp<volScalarField> trAtU =
            1.0/(1.0/rAU - UEqn.H1());
        const volScalarField q(trAtU() - rAU);
        const volVectorField gradP(fvc::grad(p));
        const surfaceScalarField snGradP(fvc::snGrad(p));

        tmp<surfaceScalarField> tqf = fvc::interpolate(q);
        const surfaceScalarField& qf = tqf();
        const surfaceScalarField& magSf = mesh_.magSf();
        tmp<fv::snGradScheme<scalar>> tSnGradScheme =
            fv::snGradScheme<scalar>::New
            (
                mesh_,
                mesh_.snGradScheme("snGrad(" + p.name() + ')')
            );
        tmp<surfaceScalarField> tSnGradDeltaCoeffs =
            tSnGradScheme().deltaCoeffs(p);
        const surfaceScalarField& snGradDeltaCoeffs =
            tSnGradDeltaCoeffs();
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const surfaceVectorField& Sf = mesh_.Sf();
        const surfaceScalarField& w = mesh_.weights();
        tmp<surfaceVectorField> tNonOrthCorrectionVectors =
            mesh_.nonOrthCorrectionVectors();
        const surfaceVectorField& nonOrthCorrectionVectors =
            tNonOrthCorrectionVectors();

        vectorField sourceH(mesh_.nCells(), vector::zero);
        vectorField gradPSeed(mesh_.nCells(), vector::zero);
        forAll(sourceH, celli)
        {
            sourceH[celli] = barH1[celli]/mesh_.V()[celli];
            qSeed[celli] += barH1[celli] & gradP[celli];
        }

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];

            const scalar eta =
                Sf[facei] & (q[P]*sourceH[P] - q[N]*sourceH[N]);

            pOldSeedH[P] += w[facei]*eta;
            pOldSeedH[N] += (scalar(1) - w[facei])*eta;

            const scalar nCoeff =
                qf[facei]*magSf[facei]*snGradDeltaCoeffs[facei];

            pOldSeedF[P] -= nCoeff*barF2[facei];
            pOldSeedF[N] += nCoeff*barF2[facei];

            const scalar qFaceSeed =
                barF2[facei]*snGradP[facei]*magSf[facei];

            qSeed[P] += w[facei]*qFaceSeed;
            qSeed[N] += (scalar(1) - w[facei])*qFaceSeed;

            const vector corrGradSeed =
                barF2[facei]
               *qf[facei]
               *magSf[facei]
               *nonOrthCorrectionVectors[facei];

            gradPSeed[P] += w[facei]*corrGradSeed;
            gradPSeed[N] += (scalar(1) - w[facei])*corrGradSeed;
        }

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];

            const scalar eta =
                Sf[facei]
              & (
                    gradPSeed[P]/mesh_.V()[P]
                  - gradPSeed[N]/mesh_.V()[N]
                );

            pOldSeedF[P] += w[facei]*eta;
            pOldSeedF[N] += (scalar(1) - w[facei])*eta;
        }

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            const fvPatch& patch = mesh_.boundary()[patchi];
            const labelUList& faceCells = patch.faceCells();
            const fvsPatchScalarField& barF2p =
                barF2.boundaryField()[patchi];
            const fvsPatchScalarField& snGradPp =
                snGradP.boundaryField()[patchi];
            const fvsPatchScalarField& magSfp =
                magSf.boundaryField()[patchi];
            const fvsPatchScalarField& qfp =
                qf.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp =
                Sf.boundaryField()[patchi];

            tmp<scalarField> tpValueInternalCoeffs =
                p.boundaryField()[patchi].valueInternalCoeffs
                (
                    patch.weights()
                );
            const scalarField& pValueInternalCoeffs =
                tpValueInternalCoeffs();

            tmp<scalarField> tpGradientInternalCoeffs =
                p.boundaryField()[patchi].gradientInternalCoeffs();
            const scalarField& pGradientInternalCoeffs =
                tpGradientInternalCoeffs();

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];

                pOldSeedBoundaryH[celli] +=
                    pValueInternalCoeffs[facei]
                   *(Sfp[facei] & (q[celli]*sourceH[celli]));

                pOldSeedBoundaryF[celli] +=
                    barF2p[facei]
                   *qfp[facei]
                   *magSfp[facei]
                   *pGradientInternalCoeffs[facei];
            }

            if
            (
                isA<calculatedFvPatchScalarField>
                (
                    q.boundaryField()[patchi]
                )
            )
            {
                continue;
            }

            tmp<scalarField> tValueInternalCoeffs =
                q.boundaryField()[patchi].valueInternalCoeffs
                (
                    patch.weights()
                );
            const scalarField& valueInternalCoeffs =
                tValueInternalCoeffs();

            forAll(faceCells, facei)
            {
                const scalar qFaceSeed =
                    barF2p[facei]*snGradPp[facei]*magSfp[facei];

                qSeed[faceCells[facei]] +=
                    valueInternalCoeffs[facei]*qFaceSeed;
            }
        }

        pOldSeedH += pOldSeedBoundaryH;
        pOldSeedF += pOldSeedBoundaryF;
        pOldSeed = pOldSeedH;
        pOldSeed += pOldSeedF;
    };

    auto runPrePressureCheck =
        [&](const word& checkName, const SimpleMapSeed& rawSeed)
    {
        SimpleMapSeed seedNew(rawSeed);
        projectPressureSeed(seedNew);

        scalarField barpRelax =
            finalCorrectionPressureSeed(seedNew.barU);
        barpRelax += seedNew.barp;
        projectPressureDirection(barpRelax);

        scalarField barpSolve(barpRelax);
        barpSolve *= pRelaxCoeff;

        scalarField barpOldDirect(barpRelax);
        barpOldDirect *= scalar(1) - pRelaxCoeff;

        scalarField dpSolveFromStage(dpstage);
        dpSolveFromStage -=
            (scalar(1) - pRelaxCoeff)*dpold.primitiveField();
        if (mag(pRelaxCoeff) > VSMALL)
        {
            dpSolveFromStage /= pRelaxCoeff;
        }
        projectPressureDirection(dpSolveFromStage);

        surfaceScalarField barF2
        (
            IOobject
            (
                "ATCTStateMapbarF2PrePressure",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dF2Map
        );
        pressureSolveReverseF2(barpSolve, barF2);

        scalarField pOldConsistent(mesh_.nCells(), Zero);
        scalarField qConsistentSeed(mesh_.nCells(), Zero);
        scalarField pOldConsistentH(mesh_.nCells(), Zero);
        scalarField pOldConsistentF(mesh_.nCells(), Zero);
        scalarField pOldConsistentBoundaryH(mesh_.nCells(), Zero);
        scalarField pOldConsistentBoundaryF(mesh_.nCells(), Zero);
        consistentOldPressureSeed
        (
            seedNew.barU,
            barF2,
            pOldConsistent,
            qConsistentSeed,
            pOldConsistentH,
            pOldConsistentF,
            pOldConsistentBoundaryH,
            pOldConsistentBoundaryF
        );

        const scalar consistentLhs =
            vectorDot(seedNew.barU, dH1Map)
          + faceDot(barF2, dF2Map);
        const scalar consistentRhs =
            vectorDot(seedNew.barU, dH0Map)
          + faceDot(barF2, dF1Map)
          + pressureDot(pOldConsistent, dpold.primitiveField());
        const scalar consistentQCoeff =
            pressureDot(qConsistentSeed, dqMap);
        const scalar consistentRhsWithQ =
            consistentRhs + consistentQCoeff;
        const scalar consistentLhsU =
            vectorDot(seedNew.barU, dH1MapU)
          + faceDot(barF2, dF2MapU);
        const scalar consistentRhsU =
            vectorDot(seedNew.barU, dH0MapU)
          + faceDot(barF2, dF1MapU);
        const scalar consistentQCoeffU =
            pressureDot(qConsistentSeed, dqMapU);
        const scalar consistentRhsUWithQ =
            consistentRhsU + consistentQCoeffU;
        const scalar consistentLhsP =
            vectorDot(seedNew.barU, dH1MapP)
          + faceDot(barF2, dF2MapP);
        const scalar consistentHActualP =
            vectorDot(seedNew.barU, dH1MapP - dH0MapP);
        const scalar consistentFActualP =
            faceDot(barF2, dF2MapP - dF1MapP);
        const scalar consistentHTransposeP =
            pressureDot(pOldConsistentH, dpold.primitiveField());
        const scalar consistentFTransposeP =
            pressureDot(pOldConsistentF, dpold.primitiveField());
        const scalar consistentBoundaryHTransposeP =
            pressureDot(pOldConsistentBoundaryH, dpold.primitiveField());
        const scalar consistentBoundaryFTransposeP =
            pressureDot(pOldConsistentBoundaryF, dpold.primitiveField());
        const scalar consistentRhsP =
            vectorDot(seedNew.barU, dH0MapP)
          + faceDot(barF2, dF1MapP)
          + pressureDot(pOldConsistent, dpold.primitiveField());
        const scalar consistentQCoeffP =
            pressureDot(qConsistentSeed, dqMapP);
        const scalar consistentRhsPWithQ =
            consistentRhsP + consistentQCoeffP;
        const scalar consistentCoeffGap = consistentLhs - consistentRhs;
        const scalar consistentScale =
            max(max(mag(consistentLhs), mag(consistentRhs)), VSMALL);
        const scalar consistentScaleWithQ =
            max(max(mag(consistentLhs), mag(consistentRhsWithQ)), VSMALL);
        const scalar consistentScaleU =
            max(max(mag(consistentLhsU), mag(consistentRhsU)), VSMALL);
        const scalar consistentScaleUWithQ =
            max(max(mag(consistentLhsU), mag(consistentRhsUWithQ)), VSMALL);
        const scalar consistentScaleP =
            max(max(mag(consistentLhsP), mag(consistentRhsP)), VSMALL);
        const scalar consistentScalePWithQ =
            max(max(mag(consistentLhsP), mag(consistentRhsPWithQ)), VSMALL);

        scalar qSeedL1 = Zero;
        scalar qSeedMax = Zero;
        scalar dqL1 = Zero;
        scalar dqMax = Zero;
        scalar dqUL1 = Zero;
        scalar dqUMax = Zero;
        scalar dqPL1 = Zero;
        scalar dqPMax = Zero;
        forAll(qConsistentSeed, celli)
        {
            qSeedL1 += mag(qConsistentSeed[celli]);
            qSeedMax = max(qSeedMax, mag(qConsistentSeed[celli]));
            dqL1 += mag(dqMap[celli]);
            dqMax = max(dqMax, mag(dqMap[celli]));
            dqUL1 += mag(dqMapU[celli]);
            dqUMax = max(dqUMax, mag(dqMapU[celli]));
            dqPL1 += mag(dqMapP[celli]);
            dqPMax = max(dqPMax, mag(dqMapP[celli]));
        }

        SimpleMapSeed seedOld = reverseOneSimpleMapSeed(seedNew, rAtU);
        projectPressureSeed(seedOld);

        const scalar lhsH = vectorDot(seedNew.barU, dH1Map);
        const scalar lhsF = faceDot(barF2, dF2Map);
        const scalar lhsPressureDirect =
            pressureDot(barpSolve, dpSolveFromStage);
        const scalar lhsPoldDirect =
            pressureDot(barpOldDirect, dpold.primitiveField());
        const scalar lhsComposedDirect =
            lhsH + lhsPressureDirect + lhsPoldDirect;
        const scalar lhs = lhsH + lhsF + lhsPoldDirect;

        const scalar rhsU = vectorDot(seedOld.barU, dUold.primitiveField());
        const scalar rhsP = pressureDot(seedOld.barp, dpold.primitiveField());
        const scalar rhs = rhsU + rhsP;

        const scalar scale = max(max(mag(lhs), mag(rhs)), VSMALL);

        Info<< "ATC-T pre-pressure transpose check: "
            << checkName
            << " lhsH = " << lhsH
            << ", lhsF = " << lhsF
            << ", lhsPressureDirect = " << lhsPressureDirect
            << ", pressureGap = " << lhsPressureDirect - lhsF
            << ", lhsComposedDirect = " << lhsComposedDirect
            << ", lhsPoldDirect = " << lhsPoldDirect
            << ", rhsU = " << rhsU
            << ", rhsP = " << rhsP
            << ", lhs = " << lhs
            << ", rhs = " << rhs
            << ", rel = " << mag(lhs - rhs)/scale
            << ", consistentRel = "
            << mag(consistentLhs - consistentRhs)/consistentScale
            << ", consistentRelWithQ = "
            << mag(consistentLhs - consistentRhsWithQ)
              /consistentScaleWithQ
            << " (consistentLhs = " << consistentLhs
            << ", consistentRhs = " << consistentRhs
            << ", consistentQCoeff = " << consistentQCoeff
            << ", consistentRhsWithQ = " << consistentRhsWithQ
            << ", consistentCoeffGap = " << consistentCoeffGap
            << ", consistentUOnlyRel = "
            << mag(consistentLhsU - consistentRhsU)/consistentScaleU
            << ", consistentUOnlyRelWithQ = "
            << mag(consistentLhsU - consistentRhsUWithQ)
              /consistentScaleUWithQ
            << ", consistentUOnlyQCoeff = " << consistentQCoeffU
            << ", consistentPOnlyRel = "
            << mag(consistentLhsP - consistentRhsP)/consistentScaleP
            << ", consistentPOnlyRelWithQ = "
            << mag(consistentLhsP - consistentRhsPWithQ)
              /consistentScalePWithQ
            << ", consistentPOnlyQCoeff = " << consistentQCoeffP
            << ", consistentHActualP = " << consistentHActualP
            << ", consistentHTransposeP = " << consistentHTransposeP
            << ", consistentFActualP = " << consistentFActualP
            << ", consistentFTransposeP = " << consistentFTransposeP
            << ", consistentBoundaryHTransposeP = "
            << consistentBoundaryHTransposeP
            << ", consistentBoundaryFTransposeP = "
            << consistentBoundaryFTransposeP
            << ", qSeedL1 = " << qSeedL1
            << ", qSeedMax = " << qSeedMax
            << ", dqL1 = " << dqL1
            << ", dqMax = " << dqMax
            << ", dqUL1 = " << dqUL1
            << ", dqUMax = " << dqUMax
            << ", dqPL1 = " << dqPL1
            << ", dqPMax = " << dqPMax
            << ")"
            << endl;
    };

    auto runPostPressureCheck =
        [&](const word& checkName, const SimpleMapSeed& rawSeed)
    {
        SimpleMapSeed seedNew(rawSeed);
        projectPressureSeed(seedNew);

        scalarField barpRelax =
            finalCorrectionPressureSeed(seedNew.barU);
        barpRelax += seedNew.barp;
        projectPressureDirection(barpRelax);

        scalarField barpSolve(barpRelax);
        barpSolve *= pRelaxCoeff;

        scalarField barpOld(barpRelax);
        barpOld *= scalar(1) - pRelaxCoeff;

        volVectorField dUraw
        (
            IOobject
            (
                "ATCTStateMapdUraw",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dH1 - rAtU*fvc::grad(dpRelax)
        );

        const scalar lhsFinal =
            vectorDot(seedNew.barU, dUraw.primitiveField())
          + pressureDot(seedNew.barp, dpRelax.primitiveField());
        const scalar rhsFinal =
            vectorDot(seedNew.barU, dH1.primitiveField())
          + pressureDot(barpRelax, dpRelax.primitiveField());

        const scalar lhsActualFinal =
            vectorDot(seedNew.barU, dUrawStage)
          + pressureDot(seedNew.barp, dpstage);
        const scalar rhsActualFinal =
            vectorDot(seedNew.barU, dH1Map)
          + pressureDot(barpRelax, dpstage);

        const scalar lhsRelax =
            pressureDot(barpRelax, dpRelaxFromSolve);
        const scalar rhsRelax =
            pressureDot(barpSolve, dpSolve.primitiveField())
          + pressureDot(barpOld, dpold.primitiveField());

        volVectorField dUrawCombined
        (
            IOobject
            (
                "ATCTStateMapdUrawCombined",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dH1
        );
        volScalarField dpRelaxCombined(dpRelax);
        dpRelaxCombined.primitiveFieldRef() = dpRelaxFromSolve;
        dpRelaxCombined.correctBoundaryConditions();
        dUrawCombined -= rAtU*fvc::grad(dpRelaxCombined);

        const scalar lhsCombined =
            vectorDot(seedNew.barU, dUrawCombined.primitiveField())
          + pressureDot(seedNew.barp, dpRelaxFromSolve);
        const scalar rhsCombined =
            vectorDot(seedNew.barU, dH1.primitiveField())
          + pressureDot(barpSolve, dpSolve.primitiveField())
          + pressureDot(barpOld, dpold.primitiveField());

        const scalar scaleFinal =
            max(max(mag(lhsFinal), mag(rhsFinal)), VSMALL);
        const scalar scaleActualFinal =
            max(max(mag(lhsActualFinal), mag(rhsActualFinal)), VSMALL);
        const scalar scaleRelax =
            max(max(mag(lhsRelax), mag(rhsRelax)), VSMALL);
        const scalar scaleCombined =
            max(max(mag(lhsCombined), mag(rhsCombined)), VSMALL);

        Info<< "ATC-T post-pressure transpose check: "
            << checkName
            << " finalRel = " << mag(lhsFinal - rhsFinal)/scaleFinal
            << " (lhs = " << lhsFinal
            << ", rhs = " << rhsFinal
            << "), actualFinalRel = "
            << mag(lhsActualFinal - rhsActualFinal)/scaleActualFinal
            << " (lhs = " << lhsActualFinal
            << ", rhs = " << rhsActualFinal
            << "), pRelaxRel = " << mag(lhsRelax - rhsRelax)/scaleRelax
            << " (lhs = " << lhsRelax
            << ", rhs = " << rhsRelax
            << "), combinedRel = "
            << mag(lhsCombined - rhsCombined)/scaleCombined
            << " (lhs = " << lhsCombined
            << ", rhs = " << rhsCombined
            << "), alphaP = " << pRelaxCoeff
            << endl;
    };

    Info<< "ATC-T state-map transpose perturbation: eps = " << eps
        << ", maxdU = " << maxdU
        << ", maxdp = " << maxdp
        << ", needReference = " << pressureNeedsReference
        << ", pRefCell = " << pRefCell
        << ", pressureGauge = "
        <<
        (
            !pressureNeedsReference
          ? "fixed"
          : (pRefCell >= 0 && pRefCell < mesh_.nCells() ? "reference" : "mean")
        )
        << ". Internal perturbations; non-coupled boundary values held fixed."
        << endl;

    runCheck("randomSeed", randomSeed);
    runCheck("initialThermalSeed", initialSeed);
    runCheck("seedAfterOneReverse", seedAfterOne);
    runCheck("seedAfterFourReverse", seedAfterFour);

    runPrePressureCheck("initialThermalSeed", initialSeed);
    runPrePressureCheck("seedAfterOneReverse", seedAfterOne);
    runPrePressureCheck("seedAfterFourReverse", seedAfterFour);

    runPostPressureCheck("initialThermalSeed", initialSeed);
    runPostPressureCheck("seedAfterOneReverse", seedAfterOne);
    runPostPressureCheck("seedAfterFourReverse", seedAfterFour);

    restoreState();
}


void Foam::thermalAdjointSimple::checkFullStateMapTranspose()
{
    if (!checkFullStateMapTranspose_)
    {
        return;
    }

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "checkFullStateMapTranspose currently supports the validated "
            << "serial diagnostic case only."
            << exit(FatalError);
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].coupled() && mesh_.boundary()[patchi].size())
        {
            FatalErrorInFunction
                << "checkFullStateMapTranspose does not yet support coupled "
                << "non-empty patch " << mesh_.boundary()[patchi].name()
                << "."
                << exit(FatalError);
        }
    }

    volVectorField& U = primalVars_.U();
    volScalarField& p = const_cast<volScalarField&>(primalVars_.p());
    surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();

    Info<< "ATC-T full-state map transpose diagnostic: turbulence model "
        << turbulence->type()
        << ", turbulence coefficients treated as frozen by this map."
        << endl;

    const volVectorField UBase
    (
        IOobject
        (
            "ATCTFullStateUBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    const volScalarField pBase
    (
        IOobject
        (
            "ATCTFullStatepBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    const surfaceScalarField phiBase
    (
        IOobject
        (
            "ATCTFullStatePhiBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    const bool pressureNeedsReference = p.needReference();
    const label pRefCell = solverControl_().pRefCell();

    auto removeMean = [](scalarField& fld)
    {
        const label n = returnReduce(fld.size(), sumOp<label>());
        if (n > 0)
        {
            fld -= gSum(fld)/scalar(n);
        }
    };

    auto projectPressureDirection = [&](scalarField& fld)
    {
        if (!pressureNeedsReference)
        {
            return;
        }

        if (pRefCell >= 0 && pRefCell < mesh_.nCells())
        {
            fld[pRefCell] = Zero;
        }
        else
        {
            removeMean(fld);
        }
    };

    auto enforcePressureBoundaryState = [&]()
    {
        bool hasUpdateablePressureSnGrad = false;

        forAll(p.boundaryField(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            if
            (
                isA<updateablePatchTypes::updateableSnGrad>
                (
                    p.boundaryField()[patchi]
                )
            )
            {
                hasUpdateablePressureSnGrad = true;
            }
        }

        forAll(p.boundaryFieldRef(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            if (pBase.boundaryField()[patchi].fixesValue())
            {
                p.boundaryFieldRef()[patchi] ==
                    pBase.boundaryField()[patchi];
            }
        }

        if (!hasUpdateablePressureSnGrad)
        {
            p.correctBoundaryConditions();
        }
    };

    auto restoreState = [&]()
    {
        U.primitiveFieldRef() = UBase.primitiveField();
        p.primitiveFieldRef() = pBase.primitiveField();
        phi.primitiveFieldRef() = phiBase.primitiveField();

        forAll(U.boundaryFieldRef(), patchi)
        {
            U.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
        }
        U.correctBoundaryConditions();

        forAll(p.boundaryFieldRef(), patchi)
        {
            p.boundaryFieldRef()[patchi] == pBase.boundaryField()[patchi];
        }
        enforcePressureBoundaryState();

        forAll(phi.boundaryFieldRef(), patchi)
        {
            phi.boundaryFieldRef()[patchi] == phiBase.boundaryField()[patchi];
        }
    };

    auto checkStateSizes =
        [&](const SimpleMapState& state, const word& where)
    {
        if (state.U.size() != mesh_.nCells())
        {
            FatalErrorInFunction
                << "SimpleMapState U size mismatch in " << where
                << exit(FatalError);
        }
        if (state.p.size() != mesh_.nCells())
        {
            FatalErrorInFunction
                << "SimpleMapState p size mismatch in " << where
                << exit(FatalError);
        }
        if (state.phiInternal.size() != mesh_.nInternalFaces())
        {
            FatalErrorInFunction
                << "SimpleMapState phiInternal size mismatch in " << where
                << exit(FatalError);
        }
        if (state.phiBoundary.size() != mesh_.boundary().size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiBoundary list size mismatch in "
                << where << exit(FatalError);
        }
        forAll(state.phiBoundary, patchi)
        {
            if (state.phiBoundary[patchi].size() != mesh_.boundary()[patchi].size())
            {
                FatalErrorInFunction
                    << "SimpleMapState phiBoundary patch size mismatch in "
                    << where << " for patch "
                    << mesh_.boundary()[patchi].name()
                    << exit(FatalError);
            }
        }
    };

    auto checkSeedSizes =
        [&](const SimpleMapSeed& seed, const word& where)
    {
        if (seed.barU.size() != mesh_.nCells())
        {
            FatalErrorInFunction
                << "SimpleMapSeed barU size mismatch in " << where
                << exit(FatalError);
        }
        if (seed.barp.size() != mesh_.nCells())
        {
            FatalErrorInFunction
                << "SimpleMapSeed barp size mismatch in " << where
                << exit(FatalError);
        }
        if (seed.barPhiInternal.size() != mesh_.nInternalFaces())
        {
            FatalErrorInFunction
                << "SimpleMapSeed barPhiInternal size mismatch in " << where
                << exit(FatalError);
        }
        if (seed.barPhiBoundary.size() != mesh_.boundary().size())
        {
            FatalErrorInFunction
                << "SimpleMapSeed barPhiBoundary list size mismatch in "
                << where << exit(FatalError);
        }
        forAll(seed.barPhiBoundary, patchi)
        {
            if (seed.barPhiBoundary[patchi].size() != mesh_.boundary()[patchi].size())
            {
                FatalErrorInFunction
                    << "SimpleMapSeed barPhiBoundary patch size mismatch in "
                    << where << " for patch "
                    << mesh_.boundary()[patchi].name()
                    << exit(FatalError);
            }
        }
    };

    auto projectPressureSeed = [&](SimpleMapSeed& seed)
    {
        if (!pressureNeedsReference)
        {
            return;
        }

        if (pRefCell >= 0 && pRefCell < mesh_.nCells())
        {
            seed.barp[pRefCell] = Zero;
        }
        else
        {
            removeMean(seed.barp);
        }
    };

    auto baseStateFromFields = [&]()
    {
        SimpleMapState state(mesh_);
        state.U = UBase.primitiveField();
        state.p = pBase.primitiveField();
        state.phiInternal = phiBase.primitiveField();
        forAll(state.phiBoundary, patchi)
        {
            state.phiBoundary[patchi] = phiBase.boundaryField()[patchi];
        }
        return state;
    };

    auto installState = [&](const SimpleMapState& state)
    {
        checkStateSizes(state, "installState");

        U.primitiveFieldRef() = state.U;
        p.primitiveFieldRef() = state.p;
        phi.primitiveFieldRef() = state.phiInternal;

        forAll(U.boundaryFieldRef(), patchi)
        {
            U.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
        }
        U.correctBoundaryConditions();

        forAll(p.boundaryFieldRef(), patchi)
        {
            p.boundaryFieldRef()[patchi] == pBase.boundaryField()[patchi];
        }
        enforcePressureBoundaryState();

        forAll(phi.boundaryFieldRef(), patchi)
        {
            phi.boundaryFieldRef()[patchi] == state.phiBoundary[patchi];
        }
    };

    auto stateAxpy =
        [&](SimpleMapState& y, const scalar a, const SimpleMapState& x)
    {
        checkStateSizes(y, "stateAxpy y");
        checkStateSizes(x, "stateAxpy x");
        y.U += a*x.U;
        y.p += a*x.p;
        y.phiInternal += a*x.phiInternal;
        forAll(y.phiBoundary, patchi)
        {
            y.phiBoundary[patchi] += a*x.phiBoundary[patchi];
        }
    };

    auto stateScale = [&](SimpleMapState& y, const scalar a)
    {
        checkStateSizes(y, "stateScale");
        y.U *= a;
        y.p *= a;
        y.phiInternal *= a;
        forAll(y.phiBoundary, patchi)
        {
            y.phiBoundary[patchi] *= a;
        }
    };

    auto projectStateDirection = [&](SimpleMapState& state)
    {
        projectPressureDirection(state.p);
    };

    auto mapStateDifference =
        [&]
        (
            const SimpleMapState& plus,
            const SimpleMapState& minus,
            const scalar eps
        )
    {
        SimpleMapState diff(plus);
        stateAxpy(diff, scalar(-1), minus);
        stateScale(diff, scalar(1)/(2*eps));
        projectStateDirection(diff);
        return diff;
    };

    auto seedStateBlocks =
        [&]
        (
            const SimpleMapSeed& seed,
            const SimpleMapState& state,
            scalar& uPart,
            scalar& pPart,
            scalar& phiInternalPart,
            scalar& phiBoundaryPart
        )
    {
        checkSeedSizes(seed, "seedStateBlocks seed");
        checkStateSizes(state, "seedStateBlocks state");

        uPart = Zero;
        pPart = Zero;
        phiInternalPart = Zero;
        phiBoundaryPart = Zero;

        forAll(seed.barU, celli)
        {
            uPart += seed.barU[celli] & state.U[celli];
            pPart += seed.barp[celli]*state.p[celli];
        }

        forAll(seed.barPhiInternal, facei)
        {
            phiInternalPart +=
                seed.barPhiInternal[facei]*state.phiInternal[facei];
        }

        forAll(seed.barPhiBoundary, patchi)
        {
            forAll(seed.barPhiBoundary[patchi], facei)
            {
                phiBoundaryPart +=
                    seed.barPhiBoundary[patchi][facei]
                   *state.phiBoundary[patchi][facei];
            }
        }

        reduce(uPart, sumOp<scalar>());
        reduce(pPart, sumOp<scalar>());
        reduce(phiInternalPart, sumOp<scalar>());
        reduce(phiBoundaryPart, sumOp<scalar>());
    };

    auto makePerturbedMap =
        [&]
        (
            const SimpleMapState& base,
            const SimpleMapState& direction,
            const scalar eps,
            const scalar sign,
            const word& namePrefix
        )
    {
        SimpleMapState probe(base);
        stateAxpy(probe, sign*eps, direction);
        installState(probe);
        SimpleMapState mapped = primalSimpleMapStateAtFrozenState(namePrefix);
        restoreState();
        return mapped;
    };

    auto applyMx =
        [&]
        (
            const SimpleMapState& base,
            const SimpleMapState& direction,
            const scalar eps,
            const word& namePrefix
        )
    {
        SimpleMapState plus =
            makePerturbedMap(base, direction, eps, scalar(1), namePrefix + "Plus");
        SimpleMapState minus =
            makePerturbedMap(base, direction, eps, scalar(-1), namePrefix + "Minus");
        return mapStateDifference(plus, minus, eps);
    };

    auto computeRAtU = [&]()
    {
        fv::options& fvOptions(fv::options::New(this->mesh_));
        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, U)
          + turbulence->divDevReff(U)
         ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();
        UEqn.relax();
        fvOptions.constrain(UEqn);

        volScalarField rAU(1.0/UEqn.A());
        tmp<volScalarField> trAtU(rAU);
        if (solverControl_().consistent())
        {
            tmp<volVectorField> tYForH1(rAU*UEqn.H());
            trAtU = 1.0/(1.0/rAU - UEqn.H1());
        }

        return tmp<volScalarField>::New
        (
            IOobject
            (
                "ATCTFullStateRAtU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            trAtU()
        );
    };

    restoreState();
    SimpleMapState baseState = baseStateFromFields();
    tmp<volScalarField> trAtUBase = computeRAtU();
    const volScalarField& rAtUBase = trAtUBase();

    typename pTraits<vector>::labelType validVectorComponents
    (
        mesh_.validComponents<vector>()
    );
    auto projectValidVectorComponents = [&](vectorField& vf)
    {
        forAll(vf, celli)
        {
            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                if (component(validVectorComponents, cmpt) == -1)
                {
                    vf[celli][cmpt] = Zero;
                }
            }
        }
    };

    SimpleMapState direction(mesh_);
    forAll(direction.U, celli)
    {
        direction.U[celli] =
            vector
            (
                scalar((257*celli + 191) % 337)/scalar(337) - scalar(0.5),
                scalar((263*celli + 193) % 347)/scalar(347) - scalar(0.5),
                scalar((269*celli + 197) % 349)/scalar(349) - scalar(0.5)
            );
        direction.p[celli] =
            scalar((271*celli + 199) % 353)/scalar(353) - scalar(0.5);
    }
    projectValidVectorComponents(direction.U);
    projectPressureDirection(direction.p);

    forAll(direction.phiInternal, facei)
    {
        direction.phiInternal[facei] =
            scalar((277*facei + 211) % 359)/scalar(359) - scalar(0.5);
    }

    forAll(direction.phiBoundary, patchi)
    {
        forAll(direction.phiBoundary[patchi], facei)
        {
            direction.phiBoundary[patchi][facei] =
                scalar((281*(facei + 1) + 223*(patchi + 1)) % 367)
               /scalar(367) - scalar(0.5);
        }
    }

    SimpleMapState directionU(direction);
    directionU.p = Zero;
    directionU.phiInternal = Zero;
    forAll(directionU.phiBoundary, patchi)
    {
        directionU.phiBoundary[patchi] = Zero;
    }

    SimpleMapState directionP(direction);
    directionP.U = vector::zero;
    directionP.phiInternal = Zero;
    forAll(directionP.phiBoundary, patchi)
    {
        directionP.phiBoundary[patchi] = Zero;
    }

    SimpleMapState directionPhiI(direction);
    directionPhiI.U = vector::zero;
    directionPhiI.p = Zero;
    forAll(directionPhiI.phiBoundary, patchi)
    {
        directionPhiI.phiBoundary[patchi] = Zero;
    }

    SimpleMapState directionPhiB(direction);
    directionPhiB.U = vector::zero;
    directionPhiB.p = Zero;
    directionPhiB.phiInternal = Zero;

    SimpleMapState directionPhiISmooth(directionPhiI);
    directionPhiISmooth.phiInternal = Zero;
    const scalar maxMagPhiBase =
        max(gMax(mag(phiBase.primitiveField())), scalar(1));
    const scalar nearZeroPhiBase = scalar(1e-12)*maxMagPhiBase;
    forAll(directionPhiISmooth.phiInternal, facei)
    {
        const scalar phif = phiBase.primitiveField()[facei];
        if (mag(phif) > nearZeroPhiBase)
        {
            const scalar rnd =
                scalar((277*facei + 211) % 359)/scalar(359)
              - scalar(0.5);
            directionPhiISmooth.phiInternal[facei] = rnd*mag(phif);
        }
    }

    SimpleMapState directionPhiBSmooth(directionPhiB);
    forAll(directionPhiBSmooth.phiBoundary, patchi)
    {
        directionPhiBSmooth.phiBoundary[patchi] = Zero;
        const fvsPatchScalarField& basep = phiBase.boundaryField()[patchi];
        const scalar patchMaxPhi = max(gMax(mag(basep)), scalar(1));
        const scalar nearZeroPatchPhi = scalar(1e-12)*patchMaxPhi;

        forAll(directionPhiBSmooth.phiBoundary[patchi], facei)
        {
            const scalar phif = basep[facei];
            if (mag(phif) > nearZeroPatchPhi)
            {
                const scalar rnd =
                    scalar((281*(facei + 1) + 223*(patchi + 1)) % 367)
                   /scalar(367) - scalar(0.5);
                directionPhiBSmooth.phiBoundary[patchi][facei] =
                    rnd*mag(phif);
            }
        }
    }

    SimpleMapState directionSmooth(directionU);
    directionSmooth.phiInternal = directionPhiISmooth.phiInternal;
    forAll(directionSmooth.phiBoundary, patchi)
    {
        directionSmooth.phiBoundary[patchi] =
            directionPhiBSmooth.phiBoundary[patchi];
    }

    SimpleMapSeed randomSeed(mesh_);
    forAll(randomSeed.barU, celli)
    {
        randomSeed.barU[celli] =
            vector
            (
                scalar((283*celli + 227) % 373)/scalar(373) - scalar(0.5),
                scalar((293*celli + 229) % 379)/scalar(379) - scalar(0.5),
                scalar((307*celli + 233) % 383)/scalar(383) - scalar(0.5)
            );
        randomSeed.barp[celli] =
            scalar((311*celli + 239) % 389)/scalar(389) - scalar(0.5);
    }
    projectValidVectorComponents(randomSeed.barU);
    forAll(randomSeed.barPhiInternal, facei)
    {
        randomSeed.barPhiInternal[facei] =
            scalar((313*facei + 241) % 397)/scalar(397) - scalar(0.5);
    }
    forAll(randomSeed.barPhiBoundary, patchi)
    {
        forAll(randomSeed.barPhiBoundary[patchi], facei)
        {
            randomSeed.barPhiBoundary[patchi][facei] =
                scalar((317*(facei + 1) + 251*(patchi + 1)) % 401)
               /scalar(401) - scalar(0.5);
        }
    }
    projectPressureSeed(randomSeed);

    SimpleMapSeed thermalSeed(mesh_);
    tmp<surfaceScalarField> tGPhi = thermalFluxSensitivity();
    const surfaceScalarField& gPhi = tGPhi();
    thermalSeed.barPhiInternal = couplingSign_*gPhi.primitiveField();
    forAll(thermalSeed.barPhiBoundary, patchi)
    {
        thermalSeed.barPhiBoundary[patchi] =
            couplingSign_*gPhi.boundaryField()[patchi];
    }

    SimpleMapSeed seedAfterOne =
        reverseOneSimpleMapSeed(thermalSeed, rAtUBase);
    projectPressureSeed(seedAfterOne);

    SimpleMapSeed seedAfterFour(seedAfterOne);
    for (label step = 1; step < 4; ++step)
    {
        seedAfterFour = reverseOneSimpleMapSeed(seedAfterFour, rAtUBase);
        projectPressureSeed(seedAfterFour);
    }

    scalarList epsList(5, Zero);
    epsList[0] = 1e-5;
    epsList[1] = 3e-6;
    epsList[2] = 1e-6;
    epsList[3] = 3e-7;
    epsList[4] = 1e-7;

    struct SimpleMapTape
    {
        vectorField Y;
        scalarField Dfinal;
        scalarField Ddiag;
        scalarField DboundaryTotal;
        List<scalarField> DboundaryByPatch;
        scalarField Dscalar;
        scalarField Ainternal;
        scalarField AinternalTimesV;
        scalarField rAU;
        scalarField UEqnH1Coeff;
        scalarField rAtU;
        List<scalarField> rAtUBoundary;
        scalarField q;
        List<scalarField> qBoundary;
        scalarField qFaceInternal;
        List<scalarField> qFaceBoundary;

        List<vectorField> UoldBoundary;
        List<vectorField> YBoundary;
        vectorField H0;
        List<vectorField> H0Boundary;
        vectorField H1;

        scalarField F0Internal;
        List<scalarField> F0Boundary;
        scalarField F1Internal;
        List<scalarField> F1Boundary;
        scalarField F2Internal;
        List<scalarField> F2Boundary;

        scalarField pBoundaryAfterConstrainInternal;
        List<scalarField> pBoundaryAfterConstrainBoundary;
        List<scalarField> pBoundaryAfterConstrainKnownBoundary;
        scalarField pSolve;
        List<scalarField> pSolveBoundary;
        List<scalarField> pSolveKnownBoundary;
        scalarField pressureFluxInternal;
        List<scalarField> pressureFluxBoundary;
        scalarField phiNewInternal;
        List<scalarField> phiNewBoundary;

        scalarField pNew;
        List<scalarField> pNewBoundary;
        List<scalarField> pNewKnownBoundary;
        vectorField Uraw;
        vectorField UafterBC;
        vectorField Unew;

        explicit SimpleMapTape(const fvMesh& mesh)
        :
            Y(mesh.nCells(), vector::zero),
            Dfinal(mesh.nCells(), Zero),
            Ddiag(mesh.nCells(), Zero),
            DboundaryTotal(mesh.nCells(), Zero),
            DboundaryByPatch(mesh.boundary().size()),
            Dscalar(mesh.nCells(), Zero),
            Ainternal(mesh.nCells(), Zero),
            AinternalTimesV(mesh.nCells(), Zero),
            rAU(mesh.nCells(), Zero),
            UEqnH1Coeff(mesh.nCells(), Zero),
            rAtU(mesh.nCells(), Zero),
            rAtUBoundary(mesh.boundary().size()),
            q(mesh.nCells(), Zero),
            qBoundary(mesh.boundary().size()),
            qFaceInternal(mesh.nInternalFaces(), Zero),
            qFaceBoundary(mesh.boundary().size()),
            UoldBoundary(mesh.boundary().size()),
            YBoundary(mesh.boundary().size()),
            H0(mesh.nCells(), vector::zero),
            H0Boundary(mesh.boundary().size()),
            H1(mesh.nCells(), vector::zero),
            F0Internal(mesh.nInternalFaces(), Zero),
            F0Boundary(mesh.boundary().size()),
            F1Internal(mesh.nInternalFaces(), Zero),
            F1Boundary(mesh.boundary().size()),
            F2Internal(mesh.nInternalFaces(), Zero),
            F2Boundary(mesh.boundary().size()),
            pBoundaryAfterConstrainInternal(mesh.nCells(), Zero),
            pBoundaryAfterConstrainBoundary(mesh.boundary().size()),
            pBoundaryAfterConstrainKnownBoundary(mesh.boundary().size()),
            pSolve(mesh.nCells(), Zero),
            pSolveBoundary(mesh.boundary().size()),
            pSolveKnownBoundary(mesh.boundary().size()),
            pressureFluxInternal(mesh.nInternalFaces(), Zero),
            pressureFluxBoundary(mesh.boundary().size()),
            phiNewInternal(mesh.nInternalFaces(), Zero),
            phiNewBoundary(mesh.boundary().size()),
            pNew(mesh.nCells(), Zero),
            pNewBoundary(mesh.boundary().size()),
            pNewKnownBoundary(mesh.boundary().size()),
            Uraw(mesh.nCells(), vector::zero),
            UafterBC(mesh.nCells(), vector::zero),
            Unew(mesh.nCells(), vector::zero)
        {
            forAll(F0Boundary, patchi)
            {
                const label nFaces = mesh.boundary()[patchi].size();
                DboundaryByPatch[patchi].setSize(nFaces, Zero);
                UoldBoundary[patchi].setSize(nFaces, vector::zero);
                rAtUBoundary[patchi].setSize(nFaces, Zero);
                qBoundary[patchi].setSize(nFaces, Zero);
                qFaceBoundary[patchi].setSize(nFaces, Zero);
                YBoundary[patchi].setSize(nFaces, vector::zero);
                H0Boundary[patchi].setSize(nFaces, vector::zero);
                F0Boundary[patchi].setSize(nFaces, Zero);
                F1Boundary[patchi].setSize(nFaces, Zero);
                F2Boundary[patchi].setSize(nFaces, Zero);
                pBoundaryAfterConstrainBoundary[patchi].setSize
                (
                    nFaces,
                    Zero
                );
                pBoundaryAfterConstrainKnownBoundary[patchi].setSize
                (
                    nFaces,
                    Zero
                );
                pSolveBoundary[patchi].setSize(nFaces, Zero);
                pSolveKnownBoundary[patchi].setSize(nFaces, Zero);
                pressureFluxBoundary[patchi].setSize(nFaces, Zero);
                phiNewBoundary[patchi].setSize(nFaces, Zero);
                pNewBoundary[patchi].setSize(nFaces, Zero);
                pNewKnownBoundary[patchi].setSize(nFaces, Zero);
            }
        }
    };

    struct PressureFluxProbeTape
    {
        scalarField totalInternal;
        List<scalarField> totalBoundary;
        scalarField coefficientInternal;
        List<scalarField> coefficientBoundary;
        scalarField correctionInternal;
        List<scalarField> correctionBoundary;

        explicit PressureFluxProbeTape(const fvMesh& mesh)
        :
            totalInternal(mesh.nInternalFaces(), Zero),
            totalBoundary(mesh.boundary().size()),
            coefficientInternal(mesh.nInternalFaces(), Zero),
            coefficientBoundary(mesh.boundary().size()),
            correctionInternal(mesh.nInternalFaces(), Zero),
            correctionBoundary(mesh.boundary().size())
        {
            forAll(totalBoundary, patchi)
            {
                const label nFaces = mesh.boundary()[patchi].size();
                totalBoundary[patchi].setSize(nFaces, Zero);
                coefficientBoundary[patchi].setSize(nFaces, Zero);
                correctionBoundary[patchi].setSize(nFaces, Zero);
            }
        }
    };

    auto copySurfaceToTape =
        [&](const surfaceScalarField& sf, scalarField& internal, List<scalarField>& boundary)
    {
        internal = sf.primitiveField();
        forAll(boundary, patchi)
        {
            boundary[patchi] = sf.boundaryField()[patchi];
        }
    };

    auto copyBoundaryPToTape =
        [&](const volScalarField& pf, SimpleMapTape& tape)
    {
        tape.pBoundaryAfterConstrainInternal = pf.primitiveField();
        forAll(tape.pBoundaryAfterConstrainBoundary, patchi)
        {
            tape.pBoundaryAfterConstrainBoundary[patchi] =
                pf.boundaryField()[patchi];
        }
    };

    auto copyKnownBoundaryContribution =
        [&](const volScalarField& pf, List<scalarField>& knownBoundary)
    {
        forAll(knownBoundary, patchi)
        {
            knownBoundary[patchi] = Zero;

            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            const fvPatch& patch = mesh_.boundary()[patchi];
            const labelUList& faceCells = patch.faceCells();

            tmp<scalarField> tValueInternalCoeffs =
                pf.boundaryField()[patchi].valueInternalCoeffs
                (
                    patch.weights()
                );
            const scalarField& valueInternalCoeffs =
                tValueInternalCoeffs();

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                knownBoundary[patchi][facei] =
                    pf.boundaryField()[patchi][facei]
                  - valueInternalCoeffs[facei]*pf[celli];
            }
        }
    };

    auto buildTapeAtState =
        [&](const SimpleMapState& state, const word& namePrefix)
    {
        installState(state);

        SimpleMapTape tape(mesh_);
        fv::options& localFvOptions(fv::options::New(this->mesh_));

        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, U)
          + turbulence->divDevReff(U)
         ==
            localFvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();
        UEqn.relax();
        localFvOptions.constrain(UEqn);

        forAll(tape.UoldBoundary, patchi)
        {
            tape.UoldBoundary[patchi] = U.boundaryField()[patchi];
        }

        volScalarField rAUField(1.0/UEqn.A());
        tmp<scalarField> tDScalar = UEqn.D();
        const scalarField& DScalar = tDScalar();
        tmp<volScalarField> tAField = UEqn.A();
        const volScalarField& AField = tAField();
        volVectorField YField(rAUField*UEqn.H());
        volScalarField H1Coeff(UEqn.H1());

        tape.rAU = rAUField.primitiveField();
        tape.Dfinal = UEqn.diag();
        tape.Ddiag = UEqn.diag();
        tape.DboundaryTotal = Zero;
        forAll(tape.DboundaryByPatch, patchi)
        {
            tape.DboundaryByPatch[patchi] = Zero;

            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            const vectorField& patchCoeffs = UEqn.internalCoeffs()[patchi];
            const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

            forAll(patchCoeffs, facei)
            {
                const label celli = faceCells[facei];
                const scalar contribution = cmptAv(patchCoeffs[facei]);
                tape.DboundaryByPatch[patchi][facei] = contribution;
                tape.DboundaryTotal[celli] += contribution;
            }
        }
        tape.Dscalar = DScalar;
        tape.Ainternal = AField.primitiveField();
        forAll(tape.AinternalTimesV, celli)
        {
            tape.AinternalTimesV[celli] =
                tape.Ainternal[celli]*mesh_.V()[celli];
        }
        tape.UEqnH1Coeff = H1Coeff.primitiveField();
        tape.Y = YField.primitiveField();
        forAll(tape.YBoundary, patchi)
        {
            tape.YBoundary[patchi] = YField.boundaryField()[patchi];
        }

        volScalarField pMap
        (
            IOobject
            (
                namePrefix + "TapeP",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p
        );
        pMap.storePrevIter();

        volVectorField H0Field(constrainHbyA(YField, U, pMap));
        tape.H0 = H0Field.primitiveField();
        forAll(tape.H0Boundary, patchi)
        {
            tape.H0Boundary[patchi] = H0Field.boundaryField()[patchi];
        }

        surfaceScalarField F0(namePrefix + "TapeF0", fvc::flux(H0Field));
        copySurfaceToTape(F0, tape.F0Internal, tape.F0Boundary);

        surfaceScalarField F1
        (
            IOobject
            (
                namePrefix + "TapeF1",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            F0
        );
        adjustPhi(F1, U, pMap);
        copySurfaceToTape(F1, tape.F1Internal, tape.F1Boundary);

        tmp<volScalarField> trAtU(rAUField);
        volVectorField H1Velocity(H0Field);
        surfaceScalarField F2
        (
            IOobject
            (
                namePrefix + "TapeF2",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            F1
        );

        if (solverControl_().consistent())
        {
            trAtU = 1.0/(1.0/rAUField - H1Coeff);
            volScalarField qField(trAtU() - rAUField);
            tape.q = qField.primitiveField();
            forAll(tape.qBoundary, patchi)
            {
                tape.qBoundary[patchi] = qField.boundaryField()[patchi];
            }
            surfaceScalarField qFace(fvc::interpolate(qField));
            copySurfaceToTape(qFace, tape.qFaceInternal, tape.qFaceBoundary);
            F2 += qFace*fvc::snGrad(pMap)*mesh_.magSf();
            H1Velocity -= (rAUField - trAtU())*fvc::grad(pMap);
        }

        tape.rAtU = trAtU().primitiveField();
        forAll(tape.rAtUBoundary, patchi)
        {
            tape.rAtUBoundary[patchi] =
                trAtU().boundaryField()[patchi];
        }
        tape.H1 = H1Velocity.primitiveField();
        copySurfaceToTape(F2, tape.F2Internal, tape.F2Boundary);

        const volScalarField& rAtUTape = trAtU();
        constrainPressure(pMap, U, F2, rAtUTape);
        copyBoundaryPToTape(pMap, tape);
        copyKnownBoundaryContribution
        (
            pMap,
            tape.pBoundaryAfterConstrainKnownBoundary
        );
        mesh_.setFluxRequired(pMap.name());

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtUTape, pMap) == fvc::div(F2)
        );
        pEqn.setReference
        (
            solverControl_().pRefCell(),
            solverControl_().pRefValue()
        );
        dictionary pSolver(pEqn.solverDict("p"));
        pSolver.set("relTol", scalar(0));
        pSolver.set("tolerance", scalar(1e-12));
        pEqn.solve(pSolver);

        tape.pSolve = pMap.primitiveField();
        forAll(tape.pSolveBoundary, patchi)
        {
            tape.pSolveBoundary[patchi] = pMap.boundaryField()[patchi];
        }
        copyKnownBoundaryContribution(pMap, tape.pSolveKnownBoundary);

        surfaceScalarField pFlux
        (
            IOobject
            (
                namePrefix + "TapePressureFlux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pEqn.flux()
        );
        copySurfaceToTape
        (
            pFlux,
            tape.pressureFluxInternal,
            tape.pressureFluxBoundary
        );

        surfaceScalarField phiNew
        (
            IOobject
            (
                namePrefix + "TapePhiNew",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            F2 - pFlux
        );
        copySurfaceToTape(phiNew, tape.phiNewInternal, tape.phiNewBoundary);

        scalar pRelaxCoeff = scalar(1);
        word pRelaxName = p.name();
        if (p.mesh().data().isFinalIteration())
        {
            pRelaxName += "Final";
        }
        p.mesh().relaxField(pRelaxName, pRelaxCoeff);
        pMap.relax(pRelaxCoeff);
        tape.pNew = pMap.primitiveField();
        forAll(tape.pNewBoundary, patchi)
        {
            tape.pNewBoundary[patchi] = pMap.boundaryField()[patchi];
        }
        copyKnownBoundaryContribution(pMap, tape.pNewKnownBoundary);

        volVectorField Uraw
        (
            IOobject
            (
                namePrefix + "TapeUraw",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            H1Velocity - rAtUTape*fvc::grad(pMap)
        );
        tape.Uraw = Uraw.primitiveField();

        volVectorField UafterBC(Uraw);
        UafterBC.correctBoundaryConditions();
        tape.UafterBC = UafterBC.primitiveField();

        localFvOptions.correct(UafterBC);
        tape.Unew = UafterBC.primitiveField();

        restoreState();
        return tape;
    };

    auto tapeDifference =
        [&](const SimpleMapTape& plus, const SimpleMapTape& minus, const scalar eps)
    {
        SimpleMapTape diff(mesh_);
        diff.Y = (plus.Y - minus.Y)/(2*eps);
        diff.Dfinal = (plus.Dfinal - minus.Dfinal)/(2*eps);
        diff.Ddiag = (plus.Ddiag - minus.Ddiag)/(2*eps);
        diff.DboundaryTotal =
            (plus.DboundaryTotal - minus.DboundaryTotal)/(2*eps);
        forAll(diff.DboundaryByPatch, patchi)
        {
            diff.DboundaryByPatch[patchi] =
                (
                    plus.DboundaryByPatch[patchi]
                  - minus.DboundaryByPatch[patchi]
                )/(2*eps);
        }
        diff.Dscalar = (plus.Dscalar - minus.Dscalar)/(2*eps);
        diff.Ainternal = (plus.Ainternal - minus.Ainternal)/(2*eps);
        diff.AinternalTimesV =
            (plus.AinternalTimesV - minus.AinternalTimesV)/(2*eps);
        diff.rAU = (plus.rAU - minus.rAU)/(2*eps);
        diff.UEqnH1Coeff =
            (plus.UEqnH1Coeff - minus.UEqnH1Coeff)/(2*eps);
        diff.rAtU = (plus.rAtU - minus.rAtU)/(2*eps);
        diff.q = (plus.q - minus.q)/(2*eps);
        diff.qFaceInternal =
            (plus.qFaceInternal - minus.qFaceInternal)/(2*eps);
        diff.H0 = (plus.H0 - minus.H0)/(2*eps);
        diff.H1 = (plus.H1 - minus.H1)/(2*eps);
        diff.F0Internal = (plus.F0Internal - minus.F0Internal)/(2*eps);
        diff.F1Internal = (plus.F1Internal - minus.F1Internal)/(2*eps);
        diff.F2Internal = (plus.F2Internal - minus.F2Internal)/(2*eps);
        diff.pBoundaryAfterConstrainInternal =
            (
                plus.pBoundaryAfterConstrainInternal
              - minus.pBoundaryAfterConstrainInternal
            )/(2*eps);
        diff.pSolve = (plus.pSolve - minus.pSolve)/(2*eps);
        diff.pressureFluxInternal =
            (plus.pressureFluxInternal - minus.pressureFluxInternal)/(2*eps);
        diff.phiNewInternal =
            (plus.phiNewInternal - minus.phiNewInternal)/(2*eps);
        diff.pNew = (plus.pNew - minus.pNew)/(2*eps);
        diff.Uraw = (plus.Uraw - minus.Uraw)/(2*eps);
        diff.UafterBC = (plus.UafterBC - minus.UafterBC)/(2*eps);
        diff.Unew = (plus.Unew - minus.Unew)/(2*eps);

        forAll(diff.F0Boundary, patchi)
        {
            diff.UoldBoundary[patchi] =
                (plus.UoldBoundary[patchi] - minus.UoldBoundary[patchi])
               /(2*eps);
            diff.rAtUBoundary[patchi] =
                (plus.rAtUBoundary[patchi] - minus.rAtUBoundary[patchi])
               /(2*eps);
            diff.qBoundary[patchi] =
                (plus.qBoundary[patchi] - minus.qBoundary[patchi])
               /(2*eps);
            diff.qFaceBoundary[patchi] =
                (plus.qFaceBoundary[patchi] - minus.qFaceBoundary[patchi])
               /(2*eps);
            diff.YBoundary[patchi] =
                (plus.YBoundary[patchi] - minus.YBoundary[patchi])
               /(2*eps);
            diff.H0Boundary[patchi] =
                (plus.H0Boundary[patchi] - minus.H0Boundary[patchi])
               /(2*eps);
            diff.F0Boundary[patchi] =
                (plus.F0Boundary[patchi] - minus.F0Boundary[patchi])
               /(2*eps);
            diff.F1Boundary[patchi] =
                (plus.F1Boundary[patchi] - minus.F1Boundary[patchi])
               /(2*eps);
            diff.F2Boundary[patchi] =
                (plus.F2Boundary[patchi] - minus.F2Boundary[patchi])
               /(2*eps);
            diff.pBoundaryAfterConstrainBoundary[patchi] =
                (
                    plus.pBoundaryAfterConstrainBoundary[patchi]
                  - minus.pBoundaryAfterConstrainBoundary[patchi]
                )/(2*eps);
            diff.pBoundaryAfterConstrainKnownBoundary[patchi] =
                (
                    plus.pBoundaryAfterConstrainKnownBoundary[patchi]
                  - minus.pBoundaryAfterConstrainKnownBoundary[patchi]
                )/(2*eps);
            diff.pressureFluxBoundary[patchi] =
                (
                    plus.pressureFluxBoundary[patchi]
                  - minus.pressureFluxBoundary[patchi]
                )/(2*eps);
            diff.pSolveBoundary[patchi] =
                (
                    plus.pSolveBoundary[patchi]
                  - minus.pSolveBoundary[patchi]
                )/(2*eps);
            diff.pSolveKnownBoundary[patchi] =
                (
                    plus.pSolveKnownBoundary[patchi]
                  - minus.pSolveKnownBoundary[patchi]
                )/(2*eps);
            diff.phiNewBoundary[patchi] =
                (
                    plus.phiNewBoundary[patchi]
                  - minus.phiNewBoundary[patchi]
                )/(2*eps);
            diff.pNewBoundary[patchi] =
                (
                    plus.pNewBoundary[patchi]
                  - minus.pNewBoundary[patchi]
                )/(2*eps);
            diff.pNewKnownBoundary[patchi] =
                (
                    plus.pNewKnownBoundary[patchi]
                  - minus.pNewKnownBoundary[patchi]
                )/(2*eps);
        }

        return diff;
    };

    auto tapeScalarNorm = [](const scalarField& fld)
    {
        return gSum(mag(fld));
    };

    auto tapeVectorNorm = [](const vectorField& fld)
    {
        return gSum(mag(fld));
    };

    auto tapeBoundaryNorm = [](const List<scalarField>& b)
    {
        scalar n = Zero;
        forAll(b, patchi)
        {
            n += gSum(mag(b[patchi]));
        }
        reduce(n, sumOp<scalar>());
        return n;
    };

    auto tapeBoundaryVectorNorm = [](const List<vectorField>& b)
    {
        scalar n = Zero;
        forAll(b, patchi)
        {
            n += gSum(mag(b[patchi]));
        }
        reduce(n, sumOp<scalar>());
        return n;
    };

    auto makeTapeDirection =
        [&]
        (
            const SimpleMapState& base,
            const SimpleMapState& dir,
            const scalar eps,
            const word& namePrefix
        )
    {
        SimpleMapState plus(base);
        SimpleMapState minus(base);
        stateAxpy(plus, eps, dir);
        stateAxpy(minus, -eps, dir);

        SimpleMapTape tapePlus =
            buildTapeAtState(plus, namePrefix + "Plus");
        SimpleMapTape tapeMinus =
            buildTapeAtState(minus, namePrefix + "Minus");

        return tapeDifference(tapePlus, tapeMinus, eps);
    };

    vectorField randomBarY(mesh_.nCells(), vector::zero);
    forAll(randomBarY, celli)
    {
        randomBarY[celli] =
            vector
            (
                scalar((331*celli + 257) % 409)/scalar(409) - scalar(0.5),
                scalar((337*celli + 263) % 419)/scalar(419) - scalar(0.5),
                scalar((347*celli + 269) % 421)/scalar(421) - scalar(0.5)
            );
    }
    projectValidVectorComponents(randomBarY);

    auto computePredictorYForState =
        [&]
        (
            const SimpleMapState& probe,
            const word& namePrefix
        )
    {
        installState(probe);

        fv::options& fvOptions(fv::options::New(this->mesh_));
        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, U)
          + turbulence->divDevReff(U)
         ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();
        UEqn.relax();
        fvOptions.constrain(UEqn);

        volScalarField rAU(1.0/UEqn.A());
        volVectorField Y
        (
            IOobject
            (
                namePrefix + "Y",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rAU*UEqn.H()
        );

        vectorField result(Y.primitiveField());
        restoreState();
        return result;
    };

    auto computePredictorResidualForState =
        [&]
        (
            const SimpleMapState& probe,
            const word& namePrefix
        )
    {
        installState(probe);

        fv::options& fvOptions(fv::options::New(this->mesh_));
        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, U)
          + turbulence->divDevReff(U)
         ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();
        UEqn.relax();
        fvOptions.constrain(UEqn);

        vectorField result(UEqn.residual());
        restoreState();
        return result;
    };

    auto predictorDirection =
        [&]
        (
            const SimpleMapState& dir,
            const scalar eps,
            const word& namePrefix
        )
    {
        SimpleMapState plus(baseState);
        SimpleMapState minus(baseState);
        stateAxpy(plus, eps, dir);
        stateAxpy(minus, -eps, dir);

        vectorField Yplus =
            computePredictorYForState(plus, namePrefix + "Plus");
        vectorField Yminus =
            computePredictorYForState(minus, namePrefix + "Minus");
        return (Yplus - Yminus)/(2*eps);
    };

    auto countMomentumSignCrossings =
        [&](const SimpleMapState& dir, const scalar eps)
    {
        label nCross = 0;

        forAll(phiBase.primitiveField(), facei)
        {
            const scalar base = phiBase.primitiveField()[facei];
            const scalar plus = base + eps*dir.phiInternal[facei];
            const scalar minus = base - eps*dir.phiInternal[facei];
            if ((plus > scalar(0) && minus < scalar(0))
             || (plus < scalar(0) && minus > scalar(0)))
            {
                ++nCross;
            }
        }

        forAll(phiBase.boundaryField(), patchi)
        {
            const fvsPatchScalarField& basep = phiBase.boundaryField()[patchi];
            forAll(basep, facei)
            {
                const scalar base = basep[facei];
                const scalar plus =
                    base + eps*dir.phiBoundary[patchi][facei];
                const scalar minus =
                    base - eps*dir.phiBoundary[patchi][facei];
                if ((plus > scalar(0) && minus < scalar(0))
                 || (plus < scalar(0) && minus > scalar(0)))
                {
                    ++nCross;
                }
            }
        }

        reduce(nCross, sumOp<label>());
        return nCross;
    };

    auto predictorSeedStateBlocks =
        [&]
        (
            const PredictorReverseSeed& seed,
            const SimpleMapState& state,
            scalar& uPart,
            scalar& phiInternalPart,
            scalar& phiBoundaryPart
        )
    {
        uPart = Zero;
        phiInternalPart = Zero;
        phiBoundaryPart = Zero;

        forAll(seed.barUold, celli)
        {
            uPart += seed.barUold[celli] & state.U[celli];
        }
        forAll(seed.barPhiInternal, facei)
        {
            phiInternalPart +=
                seed.barPhiInternal[facei]*state.phiInternal[facei];
        }
        forAll(seed.barPhiBoundary, patchi)
        {
            forAll(seed.barPhiBoundary[patchi], facei)
            {
                phiBoundaryPart +=
                    seed.barPhiBoundary[patchi][facei]
                   *state.phiBoundary[patchi][facei];
            }
        }

        reduce(uPart, sumOp<scalar>());
        reduce(phiInternalPart, sumOp<scalar>());
        reduce(phiBoundaryPart, sumOp<scalar>());
    };

    PredictorReverseSeed randomPredictorReverse =
        reversePredictorYFullState(randomBarY, true, true);

    {
        SimpleMapTape baseTape =
            buildTapeAtState(baseState, "ATCTFullStateBaseTape");
        installState(baseState);
        SimpleMapState mapBase =
            primalSimpleMapStateAtFrozenState("ATCTFullStateBaseMap");
        restoreState();

        scalar tapeDiff = Zero;
        scalar tapeScale = VSMALL;
        forAll(mapBase.U, celli)
        {
            tapeDiff += mag(baseTape.Unew[celli] - mapBase.U[celli]);
            tapeDiff += mag(baseTape.pNew[celli] - mapBase.p[celli]);
            tapeScale += max(mag(baseTape.Unew[celli]), mag(mapBase.U[celli]));
            tapeScale += max(mag(baseTape.pNew[celli]), mag(mapBase.p[celli]));
        }
        forAll(mapBase.phiInternal, facei)
        {
            tapeDiff +=
                mag(baseTape.phiNewInternal[facei] - mapBase.phiInternal[facei]);
            tapeScale +=
                max
                (
                    mag(baseTape.phiNewInternal[facei]),
                    mag(mapBase.phiInternal[facei])
                );
        }
        forAll(mapBase.phiBoundary, patchi)
        {
            forAll(mapBase.phiBoundary[patchi], facei)
            {
                tapeDiff +=
                    mag
                    (
                        baseTape.phiNewBoundary[patchi][facei]
                      - mapBase.phiBoundary[patchi][facei]
                    );
                tapeScale +=
                    max
                    (
                        mag(baseTape.phiNewBoundary[patchi][facei]),
                        mag(mapBase.phiBoundary[patchi][facei])
                    );
            }
        }
        reduce(tapeDiff, sumOp<scalar>());
        reduce(tapeScale, sumOp<scalar>());

        Info<< "ATC-T full-state tape reproduction rel "
            << tapeDiff/tapeScale << endl;
    }

    SimpleMapTape baseTapeForStages =
        buildTapeAtState(baseState, "ATCTFullStateBaseTapeStages");

    auto tapeOutputContraction =
        [&](const SimpleMapSeed& seed, const SimpleMapTape& dtape)
    {
        scalar part = Zero;
        forAll(seed.barU, celli)
        {
            part += seed.barU[celli] & dtape.Unew[celli];
            part += seed.barp[celli]*dtape.pNew[celli];
        }
        forAll(seed.barPhiInternal, facei)
        {
            part += seed.barPhiInternal[facei]*dtape.phiNewInternal[facei];
        }
        forAll(seed.barPhiBoundary, patchi)
        {
            forAll(seed.barPhiBoundary[patchi], facei)
            {
                part +=
                    seed.barPhiBoundary[patchi][facei]
                   *dtape.phiNewBoundary[patchi][facei];
            }
        }
        reduce(part, sumOp<scalar>());
        return part;
    };

    Info<< "ATC-T full-state U-only tape columns: eps "
        << "normDrAU normDH1Coeff normDRAtU normDQ "
        << "normDY normDH0 normDF0 normDF1 normDH1 normDF2 "
        << "normDPBoundaryInternal normDPBoundaryBoundary "
        << "normDPBoundaryKnown normDPSolve normDPSolveBoundary "
        << "normDPSolveKnown "
        << "normDPressureFlux normDPhiNew "
        << "normDPNew normDPNewBoundary normDPNewKnown "
        << "pRelaxAlpha pRelaxInternalMaxErr pRelaxBoundaryMaxErr "
        << "pRelaxInternalRel pRelaxBoundaryRel "
        << "normDUraw normDUafterBC normDUnew "
        << "randomOutputContraction thermalOutputContraction"
        << endl;

    scalar pRelaxCoeffTape = scalar(1);
    {
        word pRelaxName = p.name();
        if (p.mesh().data().isFinalIteration())
        {
            pRelaxName += "Final";
        }
        p.mesh().relaxField(pRelaxName, pRelaxCoeffTape);
    }

    auto tapeVectorDot = [&](const vectorField& a, const vectorField& b)
    {
        scalar value = Zero;
        forAll(a, celli)
        {
            value += a[celli] & b[celli];
        }
        reduce(value, sumOp<scalar>());
        return value;
    };

    auto tapeScalarDot = [&](const scalarField& a, const scalarField& b)
    {
        scalar value = Zero;
        forAll(a, celli)
        {
            value += a[celli]*b[celli];
        }
        reduce(value, sumOp<scalar>());
        return value;
    };

    auto scalarFieldIdentityMetrics =
        [&]
        (
            const scalarField& a,
            const scalarField& b,
            scalar& relDiff,
            scalar& maxDiff
        )
    {
        scalar diff = Zero;
        scalar scale = VSMALL;
        maxDiff = Zero;
        forAll(a, celli)
        {
            const scalar localDiff = mag(a[celli] - b[celli]);
            diff += localDiff;
            scale += max(mag(a[celli]), mag(b[celli]));
            if (localDiff > maxDiff)
            {
                maxDiff = localDiff;
            }
        }
        reduce(diff, sumOp<scalar>());
        reduce(scale, sumOp<scalar>());
        reduce(maxDiff, maxOp<scalar>());
        relDiff = diff/scale;
    };

    {
        scalarField invA(mesh_.nCells(), Zero);
        scalarField VOverDdiag(mesh_.nCells(), Zero);
        scalarField VOverDscalar(mesh_.nCells(), Zero);
        forAll(invA, celli)
        {
            invA[celli] = scalar(1)/baseTapeForStages.Ainternal[celli];
            VOverDdiag[celli] =
                mesh_.V()[celli]/baseTapeForStages.Ddiag[celli];
            VOverDscalar[celli] =
                mesh_.V()[celli]/baseTapeForStages.Dscalar[celli];
        }

        scalar relDiagVsD = Zero, maxDiagVsD = Zero;
        scalar relDiagVsAV = Zero, maxDiagVsAV = Zero;
        scalar relDVsAV = Zero, maxDVsAV = Zero;
        scalar relRAUVsInvA = Zero, maxRAUVsInvA = Zero;
        scalar relRAUVsVDiag = Zero, maxRAUVsVDiag = Zero;
        scalar relRAUVsVD = Zero, maxRAUVsVD = Zero;
        scalar relDDecomposition = Zero, maxDDecomposition = Zero;

        scalarField DdiagPlusBoundary(baseTapeForStages.Ddiag);
        DdiagPlusBoundary += baseTapeForStages.DboundaryTotal;

        scalarFieldIdentityMetrics
        (
            baseTapeForStages.Ddiag,
            baseTapeForStages.Dscalar,
            relDiagVsD,
            maxDiagVsD
        );
        scalarFieldIdentityMetrics
        (
            baseTapeForStages.Ddiag,
            baseTapeForStages.AinternalTimesV,
            relDiagVsAV,
            maxDiagVsAV
        );
        scalarFieldIdentityMetrics
        (
            baseTapeForStages.Dscalar,
            baseTapeForStages.AinternalTimesV,
            relDVsAV,
            maxDVsAV
        );
        scalarFieldIdentityMetrics
        (
            baseTapeForStages.rAU,
            invA,
            relRAUVsInvA,
            maxRAUVsInvA
        );
        scalarFieldIdentityMetrics
        (
            baseTapeForStages.rAU,
            VOverDdiag,
            relRAUVsVDiag,
            maxRAUVsVDiag
        );
        scalarFieldIdentityMetrics
        (
            baseTapeForStages.rAU,
            VOverDscalar,
            relRAUVsVD,
            maxRAUVsVD
        );
        scalarFieldIdentityMetrics
        (
            DdiagPlusBoundary,
            baseTapeForStages.Dscalar,
            relDDecomposition,
            maxDDecomposition
        );

        Info<< "ATC-T coefficient base identity: "
            << "relDdiagDscalar " << relDiagVsD
            << " maxDdiagDscalar " << maxDiagVsD
            << " relDdiagATimesV " << relDiagVsAV
            << " maxDdiagATimesV " << maxDiagVsAV
            << " relDscalarATimesV " << relDVsAV
            << " maxDscalarATimesV " << maxDVsAV
            << " relRAUInvA " << relRAUVsInvA
            << " maxRAUInvA " << maxRAUVsInvA
            << " relRAUVOverDdiag " << relRAUVsVDiag
            << " maxRAUVOverDdiag " << maxRAUVsVDiag
            << " relRAUVOverDscalar " << relRAUVsVD
            << " maxRAUVOverDscalar " << maxRAUVsVD
            << " relDdiagPlusBoundaryDscalar " << relDDecomposition
            << " maxDdiagPlusBoundaryDscalar " << maxDDecomposition
            << endl;

        forAll(baseTapeForStages.DboundaryByPatch, patchi)
        {
            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            Info<< "ATC-T coefficient D boundary contribution: patch "
                << mesh_.boundary()[patchi].name()
                << " sum "
                << gSum(baseTapeForStages.DboundaryByPatch[patchi])
                << " l1 "
                << gSum(mag(baseTapeForStages.DboundaryByPatch[patchi]))
                << endl;
        }
    }

    auto tapePhiDot =
        [&](const SimpleMapSeed& seed, const SimpleMapTape& dtape)
    {
        scalar value = Zero;
        forAll(seed.barPhiInternal, facei)
        {
            value +=
                seed.barPhiInternal[facei]*dtape.phiNewInternal[facei];
        }
        forAll(seed.barPhiBoundary, patchi)
        {
            forAll(seed.barPhiBoundary[patchi], facei)
            {
                value +=
                    seed.barPhiBoundary[patchi][facei]
                   *dtape.phiNewBoundary[patchi][facei];
            }
        }
        reduce(value, sumOp<scalar>());
        return value;
    };

    auto tapeF2Dot =
        [&](const surfaceScalarField& seed, const SimpleMapTape& dtape)
    {
        scalar value = Zero;
        forAll(seed.primitiveField(), facei)
        {
            value += seed.primitiveField()[facei]*dtape.F2Internal[facei];
        }
        forAll(seed.boundaryField(), patchi)
        {
            const fvsPatchScalarField& seedp = seed.boundaryField()[patchi];
            forAll(seedp, facei)
            {
                value += seedp[facei]*dtape.F2Boundary[patchi][facei];
            }
        }
        reduce(value, sumOp<scalar>());
        return value;
    };

    auto tapeBoundaryVectorDot =
        [&](const List<vectorField>& a, const List<vectorField>& b)
    {
        scalar value = Zero;
        forAll(a, patchi)
        {
            forAll(a[patchi], facei)
            {
                value += a[patchi][facei] & b[patchi][facei];
            }
        }
        reduce(value, sumOp<scalar>());
        return value;
    };

    auto collapsedFinalPressureSeed =
        [&](const SimpleMapSeed& seed)
    {
        const labelUList& ownLocal = mesh_.owner();
        const labelUList& neiLocal = mesh_.neighbour();
        const surfaceVectorField& SfLocal = mesh_.Sf();
        const surfaceScalarField& wLocal = mesh_.weights();

        vectorField source(mesh_.nCells(), vector::zero);
        forAll(source, celli)
        {
            source[celli] = seed.barU[celli]/mesh_.V()[celli];
        }

        scalarField barp(seed.barp);
        forAll(ownLocal, facei)
        {
            const label P = ownLocal[facei];
            const label N = neiLocal[facei];
            const scalar SfSource =
                SfLocal[facei]
              & (
                    rAtUBase[P]*source[P]
                  - rAtUBase[N]*source[N]
                );

            barp[P] -= wLocal[facei]*SfSource;
            barp[N] -= (scalar(1) - wLocal[facei])*SfSource;
        }

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            const fvPatch& patch = mesh_.boundary()[patchi];
            const labelUList& faceCells = patch.faceCells();
            const fvsPatchVectorField& SfPatch =
                SfLocal.boundaryField()[patchi];

            tmp<scalarField> tValueInternalCoeffs =
                p.boundaryField()[patchi].valueInternalCoeffs
                (
                    patch.weights()
                );
            const scalarField& valueInternalCoeffs =
                tValueInternalCoeffs();

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                barp[celli] -=
                    valueInternalCoeffs[facei]
                   *(SfPatch[facei] & (rAtUBase[celli]*source[celli]));
            }
        }

        projectPressureDirection(barp);
        return barp;
    };

    auto pressureStageF2Seed =
        [&](const SimpleMapSeed& seed, const scalarField& barpSolve)
    {
        surfaceScalarField F2Base
        (
            IOobject
            (
                "ATCTFullStatePressureStageF2Base",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phiBase
        );
        F2Base.primitiveFieldRef() = baseTapeForStages.F2Internal;
        forAll(F2Base.boundaryFieldRef(), patchi)
        {
            F2Base.boundaryFieldRef()[patchi] ==
                baseTapeForStages.F2Boundary[patchi];
        }

        volScalarField pWork
        (
            IOobject
            (
                "ATCTFullStatePressureStagePWork",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p
        );
        constrainPressure(pWork, U, F2Base, rAtUBase);

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtUBase, pWork) == fvc::div(F2Base)
        );
        pEqn.setReference
        (
            solverControl_().pRefCell(),
            solverControl_().pRefValue()
        );

        scalarField pSolveSeed(barpSolve);
        const scalarField& pLower = pEqn.lower();
        const scalarField& pUpper = pEqn.upper();
        const FieldField<Field, scalar>& pInternalCoeffs =
            pEqn.internalCoeffs();
        const labelUList& ownLocal = mesh_.owner();
        const labelUList& neiLocal = mesh_.neighbour();

        forAll(ownLocal, facei)
        {
            const label P = ownLocal[facei];
            const label N = neiLocal[facei];
            const scalar barPhi = seed.barPhiInternal[facei];

            pSolveSeed[P] += pLower[facei]*barPhi;
            pSolveSeed[N] -= pUpper[facei]*barPhi;
        }

        forAll(pInternalCoeffs, patchi)
        {
            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            const scalarField& intCoeffs = pInternalCoeffs[patchi];
            const scalarField& barPhip = seed.barPhiBoundary[patchi];
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(intCoeffs, facei)
            {
                pSolveSeed[faceCells[facei]] -=
                    intCoeffs[facei]*barPhip[facei];
            }
        }

        volScalarField mu
        (
            IOobject
            (
                "ATCTFullStatePressureStageMu",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p
        );
        mu = dimensionedScalar(mu.dimensions(), Zero);

        fvScalarMatrix muEqn(fvm::laplacian(rAtUBase, mu));
        muEqn.source() = pSolveSeed;
        muEqn.setReference(solverControl_().pRefCell(), Zero);
        dictionary muSolver(muEqn.solverDict("p"));
        muSolver.set("relTol", scalar(0));
        muSolver.set("tolerance", scalar(1e-12));
        muEqn.solve(muSolver);

        surfaceScalarField barF2
        (
            IOobject
            (
                "ATCTFullStatePressureStageBarF2",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            F2Base
        );
        barF2 =
            dimensionedScalar(barF2.dimensions(), Zero);

        forAll(ownLocal, facei)
        {
            const label P = ownLocal[facei];
            const label N = neiLocal[facei];
            barF2.primitiveFieldRef()[facei] =
                seed.barPhiInternal[facei] + mu[P] - mu[N];
        }

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            fvsPatchScalarField& barF2p =
                barF2.boundaryFieldRef()[patchi];
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(barF2p, facei)
            {
                barF2p[facei] =
                    seed.barPhiBoundary[patchi][facei]
                  + mu[faceCells[facei]];
            }
        }

        return tmp<surfaceScalarField>::New(barF2);
    };

    auto printPostPressureStageDiagnostic =
        [&]
        (
            const word& seedName,
            const SimpleMapSeed& seed,
            const SimpleMapTape& dtape,
            const scalar eps
        )
    {
        const scalarField barpNew(collapsedFinalPressureSeed(seed));
        scalarField barpSolve(barpNew);
        barpSolve *= pRelaxCoeffTape;
        tmp<surfaceScalarField> tBarF2 =
            pressureStageF2Seed(seed, barpSolve);
        const surfaceScalarField& barF2 = tBarF2();

        const scalar lhs =
            tapeVectorDot(seed.barU, dtape.Unew)
          + tapeScalarDot(seed.barp, dtape.pNew)
          + tapePhiDot(seed, dtape);

        const scalar rhsFinal =
            tapeVectorDot(seed.barU, dtape.H1)
          + tapeScalarDot(barpNew, dtape.pNew)
          + tapePhiDot(seed, dtape);

        const scalar rhsRelax =
            tapeVectorDot(seed.barU, dtape.H1)
          + tapeScalarDot(barpSolve, dtape.pSolve)
          + tapePhiDot(seed, dtape);

        const scalar pressureLhs =
            tapeScalarDot(barpSolve, dtape.pSolve)
          + tapePhiDot(seed, dtape);
        const scalar pressureRhs = tapeF2Dot(barF2, dtape);

        const scalar finalScale =
            max(max(mag(lhs), mag(rhsFinal)), VSMALL);
        const scalar relaxScale =
            max(max(mag(rhsFinal), mag(rhsRelax)), VSMALL);
        const scalar pressureScale =
            max(max(mag(pressureLhs), mag(pressureRhs)), VSMALL);

        Info<< "ATC-T full-state U-only post-pressure stage: "
            << seedName
            << " eps " << eps
            << " lhsOutput " << lhs
            << " rhsAfterFinalCorrection " << rhsFinal
            << " finalRel " << mag(lhs - rhsFinal)/finalScale
            << " rhsAfterPressureRelax " << rhsRelax
            << " relaxRel " << mag(rhsFinal - rhsRelax)/relaxScale
            << " pressureLhs " << pressureLhs
            << " pressureRhs " << pressureRhs
            << " pressureRel "
            << mag(pressureLhs - pressureRhs)/pressureScale
            << endl;
    };

    auto printPrePressureBoundaryStageDiagnostic =
        [&]
        (
            const word& seedName,
            const SimpleMapSeed& seed,
            const SimpleMapTape& dtape,
            const SimpleMapState& directionState,
            const scalar eps
        )
    {
        const scalarField barpNew(collapsedFinalPressureSeed(seed));
        scalarField barpSolve(barpNew);
        barpSolve *= pRelaxCoeffTape;
        tmp<surfaceScalarField> tBarF1 =
            pressureStageF2Seed(seed, barpSolve);
        const surfaceScalarField& barF1 = tBarF1();

        surfaceScalarField F0Base
        (
            IOobject
            (
                "ATCTFullStateAdjustPhiF0Base",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phiBase
        );
        F0Base.primitiveFieldRef() = baseTapeForStages.F0Internal;
        forAll(F0Base.boundaryFieldRef(), patchi)
        {
            F0Base.boundaryFieldRef()[patchi] ==
                baseTapeForStages.F0Boundary[patchi];
        }

        tmp<surfaceScalarField> tBarF0 =
            reverseAdjustPhiExact(F0Base, barF1, U);
        const surfaceScalarField& barF0 = tBarF0();

        scalar adjustLhs = Zero;
        scalar adjustRhs = Zero;
        scalar adjustInternal = Zero;
        scalar adjustBoundary = Zero;
        label nInflow = 0;
        label nFixedOutflow = 0;
        label nAdjustableOutflow = 0;
        label nBoundarySignCrossings = 0;

        forAll(barF1.primitiveField(), facei)
        {
            adjustInternal +=
                barF1.primitiveField()[facei]*dtape.F1Internal[facei];
        }

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatchVectorField& Up =
                primalVars_.U().boundaryField()[patchi];
            const fvsPatchScalarField& F0p =
                F0Base.boundaryField()[patchi];
            const fvsPatchScalarField& barF1p =
                barF1.boundaryField()[patchi];
            const fvsPatchScalarField& barF0p =
                barF0.boundaryField()[patchi];
            const scalarField& dF0p = dtape.F0Boundary[patchi];
            const scalarField& dF1p = dtape.F1Boundary[patchi];

            if (F0p.coupled())
            {
                continue;
            }

            const bool fixedPatch =
                Up.fixesValue()
             && !isA<inletOutletFvPatchVectorField>(Up);
            const bool adjustablePatch =
                !Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up);

            forAll(F0p, facei)
            {
                adjustBoundary += barF1p[facei]*dF1p[facei];

                if (F0p[facei] < scalar(0))
                {
                    ++nInflow;
                }
                else if (fixedPatch)
                {
                    ++nFixedOutflow;
                }
                else if (adjustablePatch)
                {
                    ++nAdjustableOutflow;
                }

                if
                (
                    (F0p[facei] + eps*dF0p[facei])
                   *(F0p[facei] - eps*dF0p[facei])
                  < scalar(0)
                )
                {
                    ++nBoundarySignCrossings;
                }

                adjustRhs += barF0p[facei]*dF0p[facei];
            }
        }
        adjustLhs = adjustInternal + adjustBoundary;
        forAll(barF0.primitiveField(), facei)
        {
            adjustRhs +=
                barF0.primitiveField()[facei]*dtape.F0Internal[facei];
        }

        reduce(adjustInternal, sumOp<scalar>());
        reduce(adjustBoundary, sumOp<scalar>());
        reduce(adjustLhs, sumOp<scalar>());
        reduce(adjustRhs, sumOp<scalar>());
        reduce(nInflow, sumOp<label>());
        reduce(nFixedOutflow, sumOp<label>());
        reduce(nAdjustableOutflow, sumOp<label>());
        reduce(nBoundarySignCrossings, sumOp<label>());

        VolVectorValueSeed barH0(mesh_);
        forAll(barF0.primitiveField(), facei)
        {
            const label P = mesh_.owner()[facei];
            const label N = mesh_.neighbour()[facei];
            const vector FSeedSf =
                barF0.primitiveField()[facei]*mesh_.Sf()[facei];

            barH0.internal[P] += mesh_.weights()[facei]*FSeedSf;
            barH0.internal[N] +=
                (scalar(1) - mesh_.weights()[facei])*FSeedSf;
        }
        forAll(mesh_.boundary(), patchi)
        {
            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            const fvsPatchScalarField& barF0p =
                barF0.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp =
                mesh_.Sf().boundaryField()[patchi];

            forAll(barF0p, facei)
            {
                barH0.boundary[patchi][facei] +=
                    barF0p[facei]*Sfp[facei];
            }
        }

        ConstrainHbyAReverseSeed hbyASeed(mesh_);
        hbyASeed.barYInternal = barH0.internal;
        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            const fvPatchVectorField& Up =
                primalVars_.U().boundaryField()[patchi];
            const fvPatchScalarField& pp =
                primalVars_.p().boundaryField()[patchi];
            const bool overwritten =
                !Up.assignable()
             && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>(pp);

            if (overwritten)
            {
                hbyASeed.barUBoundary[patchi] += barH0.boundary[patchi];
            }
            else
            {
                hbyASeed.barYBoundary[patchi] += barH0.boundary[patchi];
            }
        }

        vectorField barUFromBoundary(mesh_.nCells(), vector::zero);
        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            const fvPatchVectorField& Up =
                primalVars_.U().boundaryField()[patchi];
            tmp<vectorField> tValueInternalCoeffs =
                Up.valueInternalCoeffs(mesh_.boundary()[patchi].weights());
            const vectorField& valueInternalCoeffs =
                tValueInternalCoeffs();
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(faceCells, facei)
            {
                barUFromBoundary[faceCells[facei]] += cmptMultiply
                (
                    valueInternalCoeffs[facei],
                    hbyASeed.barUBoundary[patchi][facei]
                );
            }
        }

        scalar fluxLhs = Zero;
        scalar fluxRhs = Zero;
        forAll(barF0.primitiveField(), facei)
        {
            fluxLhs += barF0.primitiveField()[facei]*dtape.F0Internal[facei];
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvsPatchScalarField& barF0p =
                barF0.boundaryField()[patchi];
            forAll(barF0p, facei)
            {
                fluxLhs += barF0p[facei]*dtape.F0Boundary[patchi][facei];
            }
        }
        reduce(fluxLhs, sumOp<scalar>());

        fluxRhs =
            tapeVectorDot(barH0.internal, dtape.H0)
          + tapeBoundaryVectorDot(barH0.boundary, dtape.H0Boundary);

        const scalar constrainLhs = fluxRhs;
        const scalar constrainRhs =
            tapeVectorDot(hbyASeed.barYInternal, dtape.Y)
          + tapeBoundaryVectorDot(hbyASeed.barYBoundary, dtape.YBoundary)
          + tapeBoundaryVectorDot(hbyASeed.barUBoundary, dtape.UoldBoundary);

        scalarField zeroD(mesh_.nCells(), Zero);
        scalarField zeroH1(mesh_.nCells(), Zero);
        PredictorReverseSeed predictorSeed =
            reversePredictorYFullState
            (
                hbyASeed.barYInternal,
                hbyASeed.barYBoundary,
                zeroD,
                zeroH1,
                true,
                true
            );
        predictorSeed.barUold += barUFromBoundary;

        const scalar predictorRhs =
            tapeVectorDot(predictorSeed.barUold, directionState.U);
        const scalar predictorLhs = constrainRhs;

        const scalar fluxScale =
            max(max(mag(fluxLhs), mag(fluxRhs)), VSMALL);
        const scalar adjustScale =
            max(max(mag(adjustLhs), mag(adjustRhs)), VSMALL);
        const scalar constrainScale =
            max(max(mag(constrainLhs), mag(constrainRhs)), VSMALL);
        const scalar predictorScale =
            max(max(mag(predictorLhs), mag(predictorRhs)), VSMALL);

        Info<< "ATC-T full-state U-only pre-pressure boundary stage: "
            << seedName
            << " eps " << eps
            << " adjustLhs " << adjustLhs
            << " adjustRhs " << adjustRhs
            << " adjustRel " << mag(adjustLhs - adjustRhs)/adjustScale
            << " adjustInternal " << adjustInternal
            << " adjustBoundary " << adjustBoundary
            << " adjustNeedReference " << primalVars_.p().needReference()
            << " nInflow " << nInflow
            << " nFixedOutflow " << nFixedOutflow
            << " nAdjustableOutflow " << nAdjustableOutflow
            << " nBoundarySignCrossings " << nBoundarySignCrossings
            << " fluxLhs " << fluxLhs
            << " fluxRhs " << fluxRhs
            << " fluxRel " << mag(fluxLhs - fluxRhs)/fluxScale
            << " constrainLhs " << constrainLhs
            << " constrainRhs " << constrainRhs
            << " constrainRel "
            << mag(constrainLhs - constrainRhs)/constrainScale
            << " predictorLhs " << predictorLhs
            << " predictorRhs " << predictorRhs
            << " predictorRel "
            << mag(predictorLhs - predictorRhs)/predictorScale
            << " barYBoundaryL1 "
            << tapeBoundaryVectorNorm(hbyASeed.barYBoundary)
            << " barUBoundaryL1 "
            << tapeBoundaryVectorNorm(hbyASeed.barUBoundary)
            << endl;
    };

    auto copySurfaceToPressureProbe =
        [&]
        (
            const surfaceScalarField& sf,
            scalarField& internal,
            List<scalarField>& boundary
        )
    {
        internal = sf.primitiveField();
        forAll(boundary, patchi)
        {
            boundary[patchi] = sf.boundaryField()[patchi];
        }
    };

    auto fillPressureFluxProbe =
        [&]
        (
            const fvScalarMatrix& pEqn,
            const volScalarField& pField,
            PressureFluxProbeTape& out,
            const word& namePrefix
        )
    {
        surfaceScalarField totalFlux
        (
            IOobject
            (
                namePrefix + "TotalPressureFlux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pEqn.flux()
        );
        copySurfaceToPressureProbe
        (
            totalFlux,
            out.totalInternal,
            out.totalBoundary
        );

        out.coefficientInternal = Zero;
        out.correctionInternal = Zero;
        forAll(out.coefficientBoundary, patchi)
        {
            out.coefficientBoundary[patchi] = Zero;
            out.correctionBoundary[patchi] = Zero;
        }

        const scalarField& lower = pEqn.lower();
        const scalarField& upper = pEqn.upper();
        const FieldField<Field, scalar>& internalCoeffs =
            pEqn.internalCoeffs();
        const FieldField<Field, scalar>& boundaryCoeffs =
            pEqn.boundaryCoeffs();

        forAll(mesh_.owner(), facei)
        {
            const label P = mesh_.owner()[facei];
            const label N = mesh_.neighbour()[facei];
            out.coefficientInternal[facei] =
                upper[facei]*pField[N] - lower[facei]*pField[P];
            out.correctionInternal[facei] =
                out.totalInternal[facei] - out.coefficientInternal[facei];
        }

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];
            if (patch.type() == "empty")
            {
                continue;
            }
            if (patch.coupled())
            {
                FatalErrorInFunction
                    << "Pressure-rAtU fixed-pressure flux diagnostic "
                    << "supports only non-coupled patches. Patch "
                    << patch.name() << " is coupled."
                    << exit(FatalError);
            }

            tmp<scalarField> tPatchInternal =
                pField.boundaryField()[patchi].patchInternalField();
            const scalarField& patchInternal = tPatchInternal();
            const scalarField& intCoeffs = internalCoeffs[patchi];
            const scalarField& bCoeffs = boundaryCoeffs[patchi];

            forAll(out.coefficientBoundary[patchi], facei)
            {
                out.coefficientBoundary[patchi][facei] =
                    intCoeffs[facei]*patchInternal[facei] - bCoeffs[facei];
                out.correctionBoundary[patchi][facei] =
                    out.totalBoundary[patchi][facei]
                  - out.coefficientBoundary[patchi][facei];
            }
        }
    };

    Info<< "ATC-T pressure-rAtU source paths: "
        << "fvMatrix::flux "
        << "/usr/lib/openfoam/openfoam2512/src/finiteVolume/fvMatrices/"
        << "fvMatrix/fvMatrix.C:1446-1511 "
        << "lduMatrix::faceH "
        << "/usr/lib/openfoam/openfoam2512/src/OpenFOAM/matrices/"
        << "lduMatrix/lduMatrix/lduMatrixTemplates.C:76-101 "
        << "gaussLaplacianScheme::fvmLaplacian "
        << "/usr/lib/openfoam/openfoam2512/src/finiteVolume/"
        << "finiteVolume/laplacianSchemes/gaussLaplacianScheme/"
        << "gaussLaplacianScheme.C:159-197 "
        << "constrainPressure "
        << "/usr/lib/openfoam/openfoam2512/src/finiteVolume/cfdTools/"
        << "general/constrainPressure/constrainPressure.C:38-78 "
        << "activeLaplacianScheme Gauss linear corrected "
        << "activeInterpolationScheme linear "
        << "activeSnGradScheme corrected "
        << "fluxFormula faceH plus boundary internalCoeffs-minus-"
        << "boundaryCoeffs plus faceFluxCorrectionPtr"
        << endl;

    surfaceScalarField pressureRAtUF2Base
    (
        IOobject
        (
            "ATCTPressureRAtUF2Base",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiBase
    );
    pressureRAtUF2Base.primitiveFieldRef() = baseTapeForStages.F2Internal;
    forAll(pressureRAtUF2Base.boundaryFieldRef(), patchi)
    {
        pressureRAtUF2Base.boundaryFieldRef()[patchi] ==
            baseTapeForStages.F2Boundary[patchi];
    }

    volScalarField pressureRAtUPSolveBase
    (
        IOobject
        (
            "ATCTPressureRAtUPSolveBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    pressureRAtUPSolveBase.storePrevIter();
    constrainPressure
    (
        pressureRAtUPSolveBase,
        U,
        pressureRAtUF2Base,
        rAtUBase
    );
    mesh_.setFluxRequired(pressureRAtUPSolveBase.name());

    fvScalarMatrix pressureRAtUBaseEqn
    (
        fvm::laplacian(rAtUBase, pressureRAtUPSolveBase)
     == fvc::div(pressureRAtUF2Base)
    );
    pressureRAtUBaseEqn.setReference
    (
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );
    dictionary pressureRAtUBaseSolver(pressureRAtUBaseEqn.solverDict("p"));
    pressureRAtUBaseSolver.set("relTol", scalar(0));
    pressureRAtUBaseSolver.set("tolerance", scalar(1e-12));
    pressureRAtUBaseEqn.solve(pressureRAtUBaseSolver);

    surfaceScalarField pressureRAtUFluxBase
    (
        IOobject
        (
            "ATCTPressureRAtUFluxBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pressureRAtUBaseEqn.flux()
    );

    scalar pressureBasePSolveDiff = Zero;
    scalar pressureBasePSolveScale = VSMALL;
    forAll(pressureRAtUPSolveBase, celli)
    {
        pressureBasePSolveDiff +=
            mag(pressureRAtUPSolveBase[celli] - baseTapeForStages.pSolve[celli]);
        pressureBasePSolveScale +=
            max
            (
                mag(pressureRAtUPSolveBase[celli]),
                mag(baseTapeForStages.pSolve[celli])
            );
    }
    forAll(pressureRAtUPSolveBase.boundaryField(), patchi)
    {
        const fvPatchScalarField& pPatch =
            pressureRAtUPSolveBase.boundaryField()[patchi];
        forAll(pPatch, facei)
        {
            pressureBasePSolveDiff +=
                mag(pPatch[facei] - baseTapeForStages.pSolveBoundary[patchi][facei]);
            pressureBasePSolveScale +=
                max
                (
                    mag(pPatch[facei]),
                    mag(baseTapeForStages.pSolveBoundary[patchi][facei])
                );
        }
    }
    reduce(pressureBasePSolveDiff, sumOp<scalar>());
    reduce(pressureBasePSolveScale, sumOp<scalar>());

    scalar pressureBaseFluxDiff = Zero;
    scalar pressureBaseFluxScale = VSMALL;
    forAll(pressureRAtUFluxBase.primitiveField(), facei)
    {
        pressureBaseFluxDiff +=
            mag
            (
                pressureRAtUFluxBase.primitiveField()[facei]
              - baseTapeForStages.pressureFluxInternal[facei]
            );
        pressureBaseFluxScale +=
            max
            (
                mag(pressureRAtUFluxBase.primitiveField()[facei]),
                mag(baseTapeForStages.pressureFluxInternal[facei])
            );
    }
    forAll(pressureRAtUFluxBase.boundaryField(), patchi)
    {
        const fvsPatchScalarField& fluxPatch =
            pressureRAtUFluxBase.boundaryField()[patchi];
        forAll(fluxPatch, facei)
        {
            pressureBaseFluxDiff +=
                mag
                (
                    fluxPatch[facei]
                  - baseTapeForStages.pressureFluxBoundary[patchi][facei]
                );
            pressureBaseFluxScale +=
                max
                (
                    mag(fluxPatch[facei]),
                    mag(baseTapeForStages.pressureFluxBoundary[patchi][facei])
                );
        }
    }
    reduce(pressureBaseFluxDiff, sumOp<scalar>());
    reduce(pressureBaseFluxScale, sumOp<scalar>());

    Info<< "ATC-T pressure-rAtU base tape reproduction: "
        << "pSolveRel " << pressureBasePSolveDiff/pressureBasePSolveScale
        << " pressureFluxRel "
        << pressureBaseFluxDiff/pressureBaseFluxScale
        << endl;

    auto pressureFluxProbeAtRAtUDirection =
        [&]
        (
            const scalarField& drAtU,
            const scalar h,
            const scalar sign,
            const bool reconstrainPressurePatch,
            const word& namePrefix
        )
    {
        volScalarField rAtUProbe
        (
            IOobject
            (
                namePrefix + "RAtUProbe",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rAtUBase
        );
        forAll(rAtUProbe, celli)
        {
            rAtUProbe.primitiveFieldRef()[celli] += sign*h*drAtU[celli];
        }
        rAtUProbe.correctBoundaryConditions();

        volScalarField pFixed
        (
            IOobject
            (
                namePrefix + "FixedPressure",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p
        );

        if (reconstrainPressurePatch)
        {
            pFixed.storePrevIter();
            constrainPressure(pFixed, U, pressureRAtUF2Base, rAtUProbe);
        }
        else
        {
            pFixed = pressureRAtUPSolveBase;
        }
        pFixed.primitiveFieldRef() =
            pressureRAtUPSolveBase.primitiveField();
        mesh_.setFluxRequired(pFixed.name());

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtUProbe, pFixed) == fvc::div(pressureRAtUF2Base)
        );
        pEqn.setReference
        (
            solverControl_().pRefCell(),
            solverControl_().pRefValue()
        );

        PressureFluxProbeTape out(mesh_);
        fillPressureFluxProbe(pEqn, pFixed, out, namePrefix);
        return out;
    };

    auto pressureFluxProbeDifference =
        [&](const PressureFluxProbeTape& plus, const PressureFluxProbeTape& minus, const scalar eps)
    {
        PressureFluxProbeTape diff(mesh_);
        diff.totalInternal =
            (plus.totalInternal - minus.totalInternal)/(2*eps);
        diff.coefficientInternal =
            (plus.coefficientInternal - minus.coefficientInternal)/(2*eps);
        diff.correctionInternal =
            (plus.correctionInternal - minus.correctionInternal)/(2*eps);

        forAll(diff.totalBoundary, patchi)
        {
            diff.totalBoundary[patchi] =
                (plus.totalBoundary[patchi] - minus.totalBoundary[patchi])
               /(2*eps);
            diff.coefficientBoundary[patchi] =
                (
                    plus.coefficientBoundary[patchi]
                  - minus.coefficientBoundary[patchi]
                )/(2*eps);
            diff.correctionBoundary[patchi] =
                (
                    plus.correctionBoundary[patchi]
                  - minus.correctionBoundary[patchi]
                )/(2*eps);
        }

        return diff;
    };

    auto pressureFluxSeedContraction =
        [&]
        (
            const scalarField& seedInternal,
            const List<scalarField>& seedBoundary,
            const scalarField& fluxInternal,
            const List<scalarField>& fluxBoundary
        )
    {
        scalar internal = Zero;
        scalar boundary = Zero;
        forAll(seedInternal, facei)
        {
            internal -= seedInternal[facei]*fluxInternal[facei];
        }
        forAll(seedBoundary, patchi)
        {
            forAll(seedBoundary[patchi], facei)
            {
                boundary -=
                    seedBoundary[patchi][facei]*fluxBoundary[patchi][facei];
            }
        }
        reduce(internal, sumOp<scalar>());
        reduce(boundary, sumOp<scalar>());
        return Pair<scalar>(internal, boundary);
    };

    auto makeLinearGammaFace =
        [&](const volScalarField& rField, const word& namePrefix)
    {
        surfaceScalarField gammaFace
        (
            IOobject
            (
                namePrefix + "GammaFace",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(rField)
        );

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];
            if (patch.size() && patch.coupled())
            {
                FatalErrorInFunction
                    << "Pressure-rAtU face-space diagnostic supports only "
                    << "non-coupled patches. Patch " << patch.name()
                    << " is coupled." << exit(FatalError);
            }
        }

        return gammaFace;
    };

    surfaceScalarField pressureRAtUGammaFaceBase =
        makeLinearGammaFace(rAtUBase, "ATCTPressureRAtUBase");

    fvScalarMatrix pressureRAtUFaceEqn
    (
        fvm::laplacian
        (
            pressureRAtUGammaFaceBase,
            pressureRAtUPSolveBase
        )
     == fvc::div(pressureRAtUF2Base)
    );
    pressureRAtUFaceEqn.setReference
    (
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );

    auto matrixCoeffRel =
        [](const scalarField& a, const scalarField& b)
    {
        scalar diff = Zero;
        scalar scale = VSMALL;
        forAll(a, i)
        {
            diff += mag(a[i] - b[i]);
            scale += max(mag(a[i]), mag(b[i]));
        }
        reduce(diff, sumOp<scalar>());
        reduce(scale, sumOp<scalar>());
        return diff/scale;
    };

    auto matrixBoundaryCoeffRel =
        []
        (
            const FieldField<Field, scalar>& a,
            const FieldField<Field, scalar>& b
        )
    {
        scalar diff = Zero;
        scalar scale = VSMALL;
        forAll(a, patchi)
        {
            forAll(a[patchi], facei)
            {
                diff += mag(a[patchi][facei] - b[patchi][facei]);
                scale += max(mag(a[patchi][facei]), mag(b[patchi][facei]));
            }
        }
        reduce(diff, sumOp<scalar>());
        reduce(scale, sumOp<scalar>());
        return diff/scale;
    };

    PressureFluxProbeTape pressureRAtUVolProbe(mesh_);
    PressureFluxProbeTape pressureRAtUFaceProbe(mesh_);
    fillPressureFluxProbe
    (
        pressureRAtUBaseEqn,
        pressureRAtUPSolveBase,
        pressureRAtUVolProbe,
        "ATCTPressureRAtUVolBase"
    );
    fillPressureFluxProbe
    (
        pressureRAtUFaceEqn,
        pressureRAtUPSolveBase,
        pressureRAtUFaceProbe,
        "ATCTPressureRAtUFaceBase"
    );

    scalar volFaceFluxDiff = Zero;
    scalar volFaceFluxScale = VSMALL;
    forAll(pressureRAtUVolProbe.totalInternal, facei)
    {
        volFaceFluxDiff +=
            mag
            (
                pressureRAtUVolProbe.totalInternal[facei]
              - pressureRAtUFaceProbe.totalInternal[facei]
            );
        volFaceFluxScale +=
            max
            (
                mag(pressureRAtUVolProbe.totalInternal[facei]),
                mag(pressureRAtUFaceProbe.totalInternal[facei])
            );
    }
    forAll(pressureRAtUVolProbe.totalBoundary, patchi)
    {
        scalar patchDiff = Zero;
        scalar patchScale = VSMALL;
        forAll(pressureRAtUVolProbe.totalBoundary[patchi], facei)
        {
            const scalar v = pressureRAtUVolProbe.totalBoundary[patchi][facei];
            const scalar f = pressureRAtUFaceProbe.totalBoundary[patchi][facei];
            patchDiff += mag(v - f);
            patchScale += max(mag(v), mag(f));
        }
        reduce(patchDiff, sumOp<scalar>());
        reduce(patchScale, sumOp<scalar>());
        Info<< "ATC-T pressure-rAtU vol-surface patch flux: "
            << mesh_.boundary()[patchi].name()
            << " rel " << patchDiff/patchScale
            << endl;

        volFaceFluxDiff += patchDiff;
        volFaceFluxScale += patchScale;
    }
    reduce(volFaceFluxDiff, sumOp<scalar>());
    reduce(volFaceFluxScale, sumOp<scalar>());

    Info<< "ATC-T pressure-rAtU vol-surface laplacian: "
        << "diagRel "
        << matrixCoeffRel(pressureRAtUBaseEqn.diag(), pressureRAtUFaceEqn.diag())
        << " upperRel "
        << matrixCoeffRel(pressureRAtUBaseEqn.upper(), pressureRAtUFaceEqn.upper())
        << " lowerRel "
        << matrixCoeffRel(pressureRAtUBaseEqn.lower(), pressureRAtUFaceEqn.lower())
        << " internalCoeffRel "
        << matrixBoundaryCoeffRel
           (
               pressureRAtUBaseEqn.internalCoeffs(),
               pressureRAtUFaceEqn.internalCoeffs()
           )
        << " boundaryCoeffRel "
        << matrixBoundaryCoeffRel
           (
               pressureRAtUBaseEqn.boundaryCoeffs(),
               pressureRAtUFaceEqn.boundaryCoeffs()
           )
        << " fluxRel " << volFaceFluxDiff/volFaceFluxScale
        << endl;

    surfaceScalarField pressureRAtUGammaFaceUnit
    (
        IOobject
        (
            "ATCTPressureRAtUGammaFaceUnit",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pressureRAtUGammaFaceBase
    );
    pressureRAtUGammaFaceUnit.primitiveFieldRef() = scalar(1);
    forAll(pressureRAtUGammaFaceUnit.boundaryFieldRef(), patchi)
    {
        pressureRAtUGammaFaceUnit.boundaryFieldRef()[patchi] == scalar(1);
    }

    fvScalarMatrix pressureRAtUUnitEqn
    (
        fvm::laplacian
        (
            pressureRAtUGammaFaceUnit,
            pressureRAtUPSolveBase
        )
     == fvc::div(pressureRAtUF2Base)
    );
    pressureRAtUUnitEqn.setReference
    (
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );
    surfaceScalarField pressureRAtUUnitFlux
    (
        IOobject
        (
            "ATCTPressureRAtUUnitFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pressureRAtUUnitEqn.flux()
    );

    auto pressureFluxFromGammaFace =
        [&]
        (
            const surfaceScalarField& gammaFace,
            const word& namePrefix
        )
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(gammaFace, pressureRAtUPSolveBase)
         == fvc::div(pressureRAtUF2Base)
        );
        pEqn.setReference
        (
            solverControl_().pRefCell(),
            solverControl_().pRefValue()
        );
        surfaceScalarField flux
        (
            IOobject
            (
                namePrefix + "Flux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pEqn.flux()
        );
        return flux;
    };

    auto gammaFaceAtRAtUDirection =
        [&]
        (
            const scalarField& drAtU,
            const scalar h,
            const scalar sign,
            const word& namePrefix
        )
    {
        volScalarField rAtUProbe
        (
            IOobject
            (
                namePrefix + "RAtUProbe",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rAtUBase
        );
        forAll(rAtUProbe, celli)
        {
            rAtUProbe.primitiveFieldRef()[celli] += sign*h*drAtU[celli];
        }
        rAtUProbe.correctBoundaryConditions();
        return makeLinearGammaFace(rAtUProbe, namePrefix);
    };

    auto validateUnitFaceDerivative =
        [&](const scalar eps)
    {
        scalar internalMaxRel = Zero;
        scalar boundaryMaxRel = Zero;
        scalar maxOffFace = Zero;
        label nInternalChecked = 0;
        label nBoundaryChecked = 0;

        labelList selectedInternal(0);
        if (mesh_.nInternalFaces())
        {
            label nSelected = 1;
            if (mesh_.nInternalFaces() > 2)
            {
                ++nSelected;
            }
            if (mesh_.nInternalFaces() > 1)
            {
                ++nSelected;
            }
            selectedInternal.setSize(nSelected);
            label selectedi = 0;
            selectedInternal[selectedi++] = 0;
            if (mesh_.nInternalFaces() > 2)
            {
                selectedInternal[selectedi++] = mesh_.nInternalFaces()/2;
            }
            if (mesh_.nInternalFaces() > 1)
            {
                selectedInternal[selectedi++] = mesh_.nInternalFaces() - 1;
            }
        }

        forAll(selectedInternal, i)
        {
            const label facei = selectedInternal[i];
            surfaceScalarField gammaPlus
            (
                IOobject
                (
                    "ATCTPressureRAtUUnitInternalPlus",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pressureRAtUGammaFaceBase
            );
            surfaceScalarField gammaMinus
            (
                IOobject
                (
                    "ATCTPressureRAtUUnitInternalMinus",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pressureRAtUGammaFaceBase
            );
            gammaPlus.primitiveFieldRef()[facei] += eps;
            gammaMinus.primitiveFieldRef()[facei] -= eps;
            const surfaceScalarField fluxPlus =
                pressureFluxFromGammaFace(gammaPlus, "ATCTUnitInternalPlus");
            const surfaceScalarField fluxMinus =
                pressureFluxFromGammaFace(gammaMinus, "ATCTUnitInternalMinus");

            forAll(fluxPlus.primitiveField(), f)
            {
                const scalar dFlux =
                    (fluxPlus.primitiveField()[f] - fluxMinus.primitiveField()[f])
                   /(2*eps);
                if (f == facei)
                {
                    const scalar unit = pressureRAtUUnitFlux.primitiveField()[f];
                    internalMaxRel =
                        max(internalMaxRel, mag(dFlux - unit)/(mag(unit) + VSMALL));
                }
                else
                {
                    maxOffFace = max(maxOffFace, mag(dFlux));
                }
            }
            forAll(fluxPlus.boundaryField(), patchi)
            {
                const fvsPatchScalarField& plusPatch =
                    fluxPlus.boundaryField()[patchi];
                const fvsPatchScalarField& minusPatch =
                    fluxMinus.boundaryField()[patchi];
                forAll(plusPatch, facej)
                {
                    maxOffFace =
                        max
                        (
                            maxOffFace,
                            mag((plusPatch[facej] - minusPatch[facej])/(2*eps))
                        );
                }
            }
            ++nInternalChecked;
        }

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];
            if (!patch.size() || patch.type() == "empty")
            {
                continue;
            }
            forAll(patch, facei)
            {
                surfaceScalarField gammaPlus
                (
                    IOobject
                    (
                        "ATCTPressureRAtUUnitBoundaryPlus",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    pressureRAtUGammaFaceBase
                );
                surfaceScalarField gammaMinus
                (
                    IOobject
                    (
                        "ATCTPressureRAtUUnitBoundaryMinus",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    pressureRAtUGammaFaceBase
                );
                gammaPlus.boundaryFieldRef()[patchi][facei] += eps;
                gammaMinus.boundaryFieldRef()[patchi][facei] -= eps;
                const surfaceScalarField fluxPlus =
                    pressureFluxFromGammaFace(gammaPlus, "ATCTUnitBoundaryPlus");
                const surfaceScalarField fluxMinus =
                    pressureFluxFromGammaFace(gammaMinus, "ATCTUnitBoundaryMinus");

                forAll(fluxPlus.primitiveField(), f)
                {
                    maxOffFace =
                        max
                        (
                            maxOffFace,
                            mag
                            (
                                (fluxPlus.primitiveField()[f]
                               - fluxMinus.primitiveField()[f])/(2*eps)
                            )
                        );
                }
                forAll(fluxPlus.boundaryField(), patchj)
                {
                    const fvsPatchScalarField& plusPatch =
                        fluxPlus.boundaryField()[patchj];
                    const fvsPatchScalarField& minusPatch =
                        fluxMinus.boundaryField()[patchj];
                    forAll(plusPatch, facej)
                    {
                        const scalar dFlux =
                            (plusPatch[facej] - minusPatch[facej])/(2*eps);
                        if (patchj == patchi && facej == facei)
                        {
                            const scalar unit =
                                pressureRAtUUnitFlux.boundaryField()[patchj][facej];
                            boundaryMaxRel =
                                max
                                (
                                    boundaryMaxRel,
                                    mag(dFlux - unit)/(mag(unit) + VSMALL)
                                );
                        }
                        else
                        {
                            maxOffFace = max(maxOffFace, mag(dFlux));
                        }
                    }
                }
                ++nBoundaryChecked;
            }
        }

        reduce(internalMaxRel, maxOp<scalar>());
        reduce(boundaryMaxRel, maxOp<scalar>());
        reduce(maxOffFace, maxOp<scalar>());
        reduce(nInternalChecked, sumOp<label>());
        reduce(nBoundaryChecked, sumOp<label>());

        Info<< "ATC-T pressure-rAtU unit face derivative: "
            << "eps " << eps
            << " nInternalChecked " << nInternalChecked
            << " nBoundaryChecked " << nBoundaryChecked
            << " internalMaxRel " << internalMaxRel
            << " boundaryMaxRel " << boundaryMaxRel
            << " maxOffFace " << maxOffFace
            << endl;
    };

    forAll(epsList, epsi)
    {
        validateUnitFaceDerivative(epsList[epsi]);
    }

    auto pressureStageAtRAtUDirection =
        [&]
        (
            const scalarField& drAtU,
            const scalar h,
            const scalar sign,
            const word& namePrefix
        )
    {
        volScalarField rAtUProbe
        (
            IOobject
            (
                namePrefix + "RAtUProbe",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rAtUBase
        );
        forAll(rAtUProbe, celli)
        {
            rAtUProbe.primitiveFieldRef()[celli] += sign*h*drAtU[celli];
        }
        rAtUProbe.correctBoundaryConditions();

        surfaceScalarField F2Base
        (
            IOobject
            (
                namePrefix + "F2Base",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phiBase
        );
        F2Base.primitiveFieldRef() = baseTapeForStages.F2Internal;
        forAll(F2Base.boundaryFieldRef(), patchi)
        {
            F2Base.boundaryFieldRef()[patchi] ==
                baseTapeForStages.F2Boundary[patchi];
        }

        volScalarField pWork
        (
            IOobject
            (
                namePrefix + "PressureStageP",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p
        );
        constrainPressure(pWork, U, F2Base, rAtUProbe);
        mesh_.setFluxRequired(pWork.name());

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtUProbe, pWork) == fvc::div(F2Base)
        );
        pEqn.setReference
        (
            solverControl_().pRefCell(),
            solverControl_().pRefValue()
        );
        dictionary pSolver(pEqn.solverDict("p"));
        pSolver.set("relTol", scalar(0));
        pSolver.set("tolerance", scalar(1e-12));
        pEqn.solve(pSolver);

        SimpleMapTape out(mesh_);
        out.pSolve = pWork.primitiveField();

        surfaceScalarField phiNew
        (
            IOobject
            (
                namePrefix + "PhiNew",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            F2Base - pEqn.flux()
        );
        out.phiNewInternal = phiNew.primitiveField();
        forAll(out.phiNewBoundary, patchi)
        {
            out.phiNewBoundary[patchi] = phiNew.boundaryField()[patchi];
        }

        return out;
    };

    forAll(epsList, epsi)
    {
        const scalar eps = epsList[epsi];
        SimpleMapTape dTapeU =
            makeTapeDirection(baseState, directionU, eps, "ATCTTapeUOnly");

        scalar pRelaxInternalMaxErr = Zero;
        scalar pRelaxInternalScale = VSMALL;
        forAll(dTapeU.pNew, celli)
        {
            const scalar expected = pRelaxCoeffTape*dTapeU.pSolve[celli];
            pRelaxInternalMaxErr =
                max(pRelaxInternalMaxErr, mag(dTapeU.pNew[celli] - expected));
            pRelaxInternalScale =
                max(pRelaxInternalScale, max(mag(dTapeU.pNew[celli]), mag(expected)));
        }
        reduce(pRelaxInternalMaxErr, maxOp<scalar>());
        reduce(pRelaxInternalScale, maxOp<scalar>());

        scalar pRelaxBoundaryMaxErr = Zero;
        scalar pRelaxBoundaryScale = VSMALL;
        forAll(dTapeU.pNewBoundary, patchi)
        {
            forAll(dTapeU.pNewBoundary[patchi], facei)
            {
                const scalar expected =
                    pRelaxCoeffTape*dTapeU.pSolveBoundary[patchi][facei];
                pRelaxBoundaryMaxErr =
                    max
                    (
                        pRelaxBoundaryMaxErr,
                        mag(dTapeU.pNewBoundary[patchi][facei] - expected)
                    );
                pRelaxBoundaryScale =
                    max
                    (
                        pRelaxBoundaryScale,
                        max
                        (
                            mag(dTapeU.pNewBoundary[patchi][facei]),
                            mag(expected)
                        )
                    );
            }
        }
        reduce(pRelaxBoundaryMaxErr, maxOp<scalar>());
        reduce(pRelaxBoundaryScale, maxOp<scalar>());

        Info<< "ATC-T full-state U-only tape: "
            << eps << " "
            << tapeScalarNorm(dTapeU.rAU) << " "
            << tapeScalarNorm(dTapeU.UEqnH1Coeff) << " "
            << tapeScalarNorm(dTapeU.rAtU) << " "
            << tapeScalarNorm(dTapeU.q) << " "
            << tapeVectorNorm(dTapeU.Y) << " "
            << tapeVectorNorm(dTapeU.H0) << " "
            << tapeScalarNorm(dTapeU.F0Internal)
               + tapeBoundaryNorm(dTapeU.F0Boundary) << " "
            << tapeScalarNorm(dTapeU.F1Internal)
               + tapeBoundaryNorm(dTapeU.F1Boundary) << " "
            << tapeVectorNorm(dTapeU.H1) << " "
            << tapeScalarNorm(dTapeU.F2Internal)
               + tapeBoundaryNorm(dTapeU.F2Boundary) << " "
            << tapeScalarNorm(dTapeU.pBoundaryAfterConstrainInternal) << " "
            << tapeBoundaryNorm(dTapeU.pBoundaryAfterConstrainBoundary) << " "
            << tapeBoundaryNorm(dTapeU.pBoundaryAfterConstrainKnownBoundary) << " "
            << tapeScalarNorm(dTapeU.pSolve) << " "
            << tapeBoundaryNorm(dTapeU.pSolveBoundary) << " "
            << tapeBoundaryNorm(dTapeU.pSolveKnownBoundary) << " "
            << tapeScalarNorm(dTapeU.pressureFluxInternal)
               + tapeBoundaryNorm(dTapeU.pressureFluxBoundary) << " "
            << tapeScalarNorm(dTapeU.phiNewInternal)
               + tapeBoundaryNorm(dTapeU.phiNewBoundary) << " "
            << tapeScalarNorm(dTapeU.pNew) << " "
            << tapeBoundaryNorm(dTapeU.pNewBoundary) << " "
            << tapeBoundaryNorm(dTapeU.pNewKnownBoundary) << " "
            << pRelaxCoeffTape << " "
            << pRelaxInternalMaxErr << " "
            << pRelaxBoundaryMaxErr << " "
            << pRelaxInternalMaxErr/pRelaxInternalScale << " "
            << pRelaxBoundaryMaxErr/pRelaxBoundaryScale << " "
            << tapeVectorNorm(dTapeU.Uraw) << " "
            << tapeVectorNorm(dTapeU.UafterBC) << " "
            << tapeVectorNorm(dTapeU.Unew) << " "
            << tapeOutputContraction(randomSeed, dTapeU) << " "
            << tapeOutputContraction(thermalSeed, dTapeU)
            << endl;

        printPostPressureStageDiagnostic
        (
            "randomSeed",
            randomSeed,
            dTapeU,
            eps
        );
        printPostPressureStageDiagnostic
        (
            "thermalSeed",
            thermalSeed,
            dTapeU,
            eps
        );
        printPrePressureBoundaryStageDiagnostic
        (
            "randomSeed",
            randomSeed,
            dTapeU,
            directionU,
            eps
        );
        printPrePressureBoundaryStageDiagnostic
        (
            "thermalSeed",
            thermalSeed,
            dTapeU,
            directionU,
            eps
        );
    }

    {
        fv::options& fvOptions(fv::options::New(this->mesh_));
        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, U)
          + turbulence->divDevReff(U)
         ==
            fvOptions(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();

        const scalarField diagRaw(UEqn.diag());
        const scalarField lowerRaw(UEqn.lower());
        const scalarField upperRaw(UEqn.upper());
        const FieldField<Field, vector> internalCoeffsRaw
        (
            UEqn.internalCoeffs()
        );
        fvVectorMatrix rawUEqnForAction(UEqn);

        UEqn.relax();
        fvOptions.constrain(UEqn);

        volScalarField rAU(1.0/UEqn.A());
        volVectorField Y
        (
            IOobject
            (
                "ATCTPredictorDiagnosticY",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rAU*UEqn.H()
        );

        vectorField barE(mesh_.nCells(), vector::zero);
        scalarField barD(mesh_.nCells(), Zero);
        forAll(barE, celli)
        {
            barE[celli] =
                rAU[celli]*randomBarY[celli]/mesh_.V()[celli];
            barD[celli] =
                (randomBarY[celli] & (U[celli] - Y[celli]))
               *rAU[celli]/mesh_.V()[celli];
        }

        const scalar eps = scalar(1e-6);
        SimpleMapState plus(baseState);
        SimpleMapState minus(baseState);
        stateAxpy(plus, eps, directionU);
        stateAxpy(minus, -eps, directionU);

        const vectorField Eplus =
            computePredictorResidualForState(plus, "ATCTResUPlus");
        const vectorField Eminus =
            computePredictorResidualForState(minus, "ATCTResUMinus");
        const vectorField dE = (Eplus - Eminus)/(2*eps);
        const vectorField dYPredictor =
            predictorDirection(directionU, eps, "ATCTPredictorIdentityU");

        scalar lhs = Zero;
        scalar dYIdentityDiff = Zero;
        scalar dYIdentityScale = VSMALL;
        scalar diagPart = Zero;
        scalar offDiagPart = Zero;
        scalar offDiagSwappedPart = Zero;
        scalar boundaryInternalPart = Zero;

        forAll(barE, celli)
        {
            lhs += barE[celli] & dE[celli];
            const vector dYIdentity =
                directionU.U[celli]
              + rAU[celli]*dE[celli]/mesh_.V()[celli];
            dYIdentityDiff += mag(dYPredictor[celli] - dYIdentity);
            dYIdentityScale +=
                max(mag(dYPredictor[celli]), mag(dYIdentity));
            diagPart -=
                barE[celli] & (diagRaw[celli]*directionU.U[celli]);
        }

        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];

            offDiagPart -=
                barE[P] & (lowerRaw[facei]*directionU.U[N]);
            offDiagPart -=
                barE[N] & (upperRaw[facei]*directionU.U[P]);

            offDiagSwappedPart -=
                barE[P] & (upperRaw[facei]*directionU.U[N]);
            offDiagSwappedPart -=
                barE[N] & (lowerRaw[facei]*directionU.U[P]);
        }

        forAll(internalCoeffsRaw, patchi)
        {
            const vectorField& coeffs = internalCoeffsRaw[patchi];
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(coeffs, facei)
            {
                const label celli = faceCells[facei];
                for (label cmpt = 0; cmpt < vector::nComponents; ++cmpt)
                {
                    boundaryInternalPart -=
                        barE[celli][cmpt]
                       *coeffs[facei][cmpt]
                       *directionU.U[celli][cmpt];
                }
            }
        }

        reduce(lhs, sumOp<scalar>());
        reduce(dYIdentityDiff, sumOp<scalar>());
        reduce(dYIdentityScale, sumOp<scalar>());
        reduce(diagPart, sumOp<scalar>());
        reduce(offDiagPart, sumOp<scalar>());
        reduce(offDiagSwappedPart, sumOp<scalar>());
        reduce(boundaryInternalPart, sumOp<scalar>());

        volVectorField dUAction
        (
            IOobject
            (
                "ATCTPredictorActionDU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U
        );
        dUAction.primitiveFieldRef() = directionU.U;
        dUAction.correctBoundaryConditions();

        rawUEqnForAction.source() = vector::zero;
        forAll(rawUEqnForAction.boundaryCoeffs(), patchi)
        {
            rawUEqnForAction.boundaryCoeffs()[patchi] = vector::zero;
        }

        tmp<volVectorField> tMatrixAction = rawUEqnForAction & dUAction;
        const vectorField& matrixAction = tMatrixAction().primitiveField();

        vectorField manualMatrix(mesh_.nCells(), vector::zero);
        forAll(manualMatrix, celli)
        {
            manualMatrix[celli] +=
                diagRaw[celli]*directionU.U[celli];
        }
        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];

            manualMatrix[P] += upperRaw[facei]*directionU.U[N];
            manualMatrix[N] += lowerRaw[facei]*directionU.U[P];
        }
        forAll(internalCoeffsRaw, patchi)
        {
            const vectorField& coeffs = internalCoeffsRaw[patchi];
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(coeffs, facei)
            {
                const label celli = faceCells[facei];
                for (label cmpt = 0; cmpt < vector::nComponents; ++cmpt)
                {
                    manualMatrix[celli][cmpt] +=
                        coeffs[facei][cmpt]
                       *directionU.U[celli][cmpt];
                }
            }
        }

        scalar matrixActionDiff = Zero;
        scalar matrixActionNegDiff = Zero;
        scalar matrixActionScale = VSMALL;
        forAll(manualMatrix, celli)
        {
            const vector matrixActionIntegrated =
                mesh_.V()[celli]*matrixAction[celli];

            matrixActionDiff +=
                mag(matrixActionIntegrated - manualMatrix[celli]);
            matrixActionNegDiff +=
                mag(matrixActionIntegrated + manualMatrix[celli]);
            matrixActionScale +=
                max(mag(matrixActionIntegrated), mag(manualMatrix[celli]));
        }
        reduce(matrixActionDiff, sumOp<scalar>());
        reduce(matrixActionNegDiff, sumOp<scalar>());
        reduce(matrixActionScale, sumOp<scalar>());

        const scalar analytic =
            diagPart + offDiagPart + boundaryInternalPart;
        const scalar analyticSwapped =
            diagPart + offDiagSwappedPart + boundaryInternalPart;
        const scalar remainder = lhs - analytic;
        const scalar remainderSwapped = lhs - analyticSwapped;

        Info<< "ATC-T predictor residual-U action check: "
            << "eps " << eps
            << " lhs " << lhs
            << " diag " << diagPart
            << " offDiagWrongOrientation " << offDiagPart
            << " offDiagOpenFOAM " << offDiagSwappedPart
            << " boundaryInternal " << boundaryInternalPart
            << " analyticWrongOrientation " << analytic
            << " wrongOrientationRemainder " << remainder
            << " wrongOrientationRelRemainder "
            << mag(remainder)/max(max(mag(lhs), mag(analytic)), VSMALL)
            << " analyticOpenFOAM " << analyticSwapped
            << " openFOAMRemainder " << remainderSwapped
            << " openFOAMRelRemainder "
            << mag(remainderSwapped)
              /max(max(mag(lhs), mag(analyticSwapped)), VSMALL)
            << " matrixActionRel "
            << matrixActionDiff/matrixActionScale
            << " matrixActionNegRel "
            << matrixActionNegDiff/matrixActionScale
            << " dYIdentityRel "
            << dYIdentityDiff/dYIdentityScale
            << endl;

        auto residualForStage =
            [&]
            (
                const SimpleMapState& state,
                const bool relaxMatrix,
                const bool constrainMatrix
            )
        {
            installState(state);

            fv::options& localFvOptions(fv::options::New(this->mesh_));
            tmp<fvVectorMatrix> tStageUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
             ==
                localFvOptions(U)
            );
            fvVectorMatrix& stageUEqn = tStageUEqn.ref();

            if (relaxMatrix)
            {
                stageUEqn.relax();
            }
            if (constrainMatrix)
            {
                localFvOptions.constrain(stageUEqn);
            }

            vectorField res(stageUEqn.residual());
            restoreState();
            return res;
        };

        auto diagForStage =
            [&]
            (
                const SimpleMapState& state,
                const bool relaxMatrix,
                const bool constrainMatrix,
                const bool useD
            )
        {
            installState(state);

            fv::options& localFvOptions(fv::options::New(this->mesh_));
            tmp<fvVectorMatrix> tStageUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
             ==
                localFvOptions(U)
            );
            fvVectorMatrix& stageUEqn = tStageUEqn.ref();

            if (relaxMatrix)
            {
                stageUEqn.relax();
            }
            if (constrainMatrix)
            {
                localFvOptions.constrain(stageUEqn);
            }

            scalarField diag(useD ? stageUEqn.D()() : stageUEqn.diag());
            restoreState();
            return diag;
        };

        auto residualStageContraction =
            [&]
            (
                const bool relaxMatrix,
                const bool constrainMatrix
            )
        {
            const vectorField plusResidual =
                residualForStage(plus, relaxMatrix, constrainMatrix);
            const vectorField minusResidual =
                residualForStage(minus, relaxMatrix, constrainMatrix);
            const vectorField dResidual =
                (plusResidual - minusResidual)/(2*eps);

            scalar contraction = Zero;
            forAll(barE, celli)
            {
                contraction += barE[celli] & dResidual[celli];
            }
            reduce(contraction, sumOp<scalar>());
            return contraction;
        };

        auto diagStageContraction =
            [&]
            (
            const bool relaxMatrix,
                const bool constrainMatrix,
                const bool useD
            )
        {
            const scalarField plusDiag =
                diagForStage(plus, relaxMatrix, constrainMatrix, useD);
            const scalarField minusDiag =
                diagForStage(minus, relaxMatrix, constrainMatrix, useD);
            const scalarField dDiag = (plusDiag - minusDiag)/(2*eps);

            scalar contraction = Zero;
            forAll(barD, celli)
            {
                contraction += barD[celli]*dDiag[celli];
            }
            reduce(contraction, sumOp<scalar>());
            return contraction;
        };

        const scalar barE_dEraw =
            residualStageContraction(false, false);
        const scalar barE_dErelaxed =
            residualStageContraction(true, false);
        const scalar barE_dEconstrained =
            residualStageContraction(true, true);

        const scalar barD_dDraw =
            diagStageContraction(false, false, false);
        const scalar barD_dDrelaxed =
            diagStageContraction(true, false, false);
        const scalar barD_dDfinal =
            diagStageContraction(true, true, false);
        const scalar barD_dDrawD =
            diagStageContraction(false, false, true);
        const scalar barD_dDrelaxedD =
            diagStageContraction(true, false, true);
        const scalar barD_dDfinalD =
            diagStageContraction(true, true, true);

        Info<< "ATC-T predictor residual-stage U derivative: "
            << "eps " << eps
            << " barE_dEraw " << barE_dEraw
            << " barE_dErelaxed " << barE_dErelaxed
            << " barE_dEconstrained " << barE_dEconstrained
            << " relaxedMinusRaw " << barE_dErelaxed - barE_dEraw
            << " constrainedMinusRelaxed "
            << barE_dEconstrained - barE_dErelaxed
            << " barD_dDraw " << barD_dDraw
            << " barD_dDrelaxed " << barD_dDrelaxed
            << " barD_dDfinal " << barD_dDfinal
            << " barD_dDrawUsingD " << barD_dDrawD
            << " barD_dDrelaxedUsingD " << barD_dDrelaxedD
            << " barD_dDfinalUsingD " << barD_dDfinalD
            << endl;

        auto sourcePlusBoundary =
            [&]
            (
                const fvVectorMatrix& matrix,
                scalarList* patchContractions,
                const vectorField* seed
            )
        {
            vectorField s(matrix.source());

            const FieldField<Field, vector>& bCoeffs =
                matrix.boundaryCoeffs();

            forAll(bCoeffs, patchi)
            {
                if (mesh_.boundary()[patchi].coupled())
                {
                    continue;
                }

                const labelUList& faceCells =
                    mesh_.boundary()[patchi].faceCells();

                forAll(bCoeffs[patchi], facei)
                {
                    const label celli = faceCells[facei];
                    s[celli] += bCoeffs[patchi][facei];

                    if (patchContractions && seed)
                    {
                        (*patchContractions)[patchi] +=
                            (*seed)[celli] & bCoeffs[patchi][facei];
                    }
                }
            }

            return s;
        };

        auto sourceLikeContraction =
            [&]
            (
                const word& termName,
                const SimpleMapState& state
            )
        {
            installState(state);

            tmp<fvVectorMatrix> tTerm;
            if (termName == "convection")
            {
                tTerm = fvm::div(phi, U);
            }
            else if (termName == "viscous")
            {
                tTerm = turbulence->divDevReff(U);
            }
            else
            {
                tTerm =
                (
                    fvm::div(phi, U)
                  + turbulence->divDevReff(U)
                 ==
                    fvOptions(U)
                );
            }

            vectorField sourceLike = sourcePlusBoundary
            (
                tTerm.ref(),
                nullptr,
                nullptr
            );
            restoreState();
            return sourceLike;
        };

        auto sourceLikeFDPart =
            [&]
            (
                const word& termName
            )
        {
            const vectorField plusSource =
                sourceLikeContraction(termName, plus);
            const vectorField minusSource =
                sourceLikeContraction(termName, minus);
            const vectorField dSource = (plusSource - minusSource)/(2*eps);

            scalar part = Zero;
            forAll(barE, celli)
            {
                part += barE[celli] & dSource[celli];
            }
            reduce(part, sumOp<scalar>());
            return part;
        };

        const scalar sourceLikeFull = sourceLikeFDPart("full");
        const scalar sourceLikeConvection = sourceLikeFDPart("convection");
        const scalar sourceLikeViscous = sourceLikeFDPart("viscous");

        const vectorField explicitReverse =
            reverseExplicitStokesViscousSource(barE);
        scalar explicitReverseContraction = Zero;
        scalar barYDotDU = Zero;
        scalar randomReverseUContraction = Zero;
        forAll(explicitReverse, celli)
        {
            explicitReverseContraction +=
                explicitReverse[celli] & directionU.U[celli];
            barYDotDU += randomBarY[celli] & directionU.U[celli];
            randomReverseUContraction +=
                randomPredictorReverse.barUold[celli]
              & directionU.U[celli];
        }
        reduce(explicitReverseContraction, sumOp<scalar>());
        reduce(barYDotDU, sumOp<scalar>());
        reduce(randomReverseUContraction, sumOp<scalar>());

        auto sourceBoundaryByPatch =
            [&](const SimpleMapState& state, scalarList& patchParts)
        {
            installState(state);

            tmp<fvVectorMatrix> tTerm = turbulence->divDevReff(U);
            patchParts.setSize(mesh_.boundary().size(), Zero);
            (void)sourcePlusBoundary(tTerm.ref(), &patchParts, &barE);
            restoreState();
        };

        scalarList viscousPatchPlus(mesh_.boundary().size(), Zero);
        scalarList viscousPatchMinus(mesh_.boundary().size(), Zero);
        sourceBoundaryByPatch(plus, viscousPatchPlus);
        sourceBoundaryByPatch(minus, viscousPatchMinus);

        Info<< "ATC-T predictor raw source-like U derivative: "
            << "eps " << eps
            << " full " << sourceLikeFull
            << " convection " << sourceLikeConvection
            << " viscous " << sourceLikeViscous
            << " explicitReverse " << explicitReverseContraction
            << " explicitReverseMinusViscous "
            << explicitReverseContraction - sourceLikeViscous
            << " barYDotDU " << barYDotDU
            << " expectedPredictorUOpenFOAM "
            << barYDotDU + analyticSwapped + explicitReverseContraction
            << " actualPredictorReverseU "
            << randomReverseUContraction
            << " fullMinusSplit "
            << sourceLikeFull - sourceLikeConvection - sourceLikeViscous
            << endl;

        forAll(viscousPatchPlus, patchi)
        {
            const scalar patchPart =
                (viscousPatchPlus[patchi] - viscousPatchMinus[patchi])
               /(2*eps);
            if (mag(patchPart) > scalar(1e-14))
            {
                Info<< "ATC-T predictor explicit viscous boundary derivative: "
                    << mesh_.boundary()[patchi].name()
                    << " " << patchPart
                    << endl;
            }
        }
    }

    Info<< "ATC-T predictor full-state transpose columns: seed direction eps "
        << "lhs rhsU rhsPhiInternal rhsPhiBoundary rhs rel "
        << "nMomentumPhiSignCrossings"
        << endl;

    forAll(epsList, epsi)
    {
        const scalar eps = epsList[epsi];
        wordList dirNames(7);
        dirNames[0] = "UOnly";
        dirNames[1] = "PhiInternalOnly";
        dirNames[2] = "PhiBoundaryOnly";
        dirNames[3] = "PhiInternalSmoothOnly";
        dirNames[4] = "PhiBoundarySmoothOnly";
        dirNames[5] = "combinedSmooth";
        dirNames[6] = "combined";

        List<const SimpleMapState*> dirs(7);
        dirs[0] = &directionU;
        dirs[1] = &directionPhiI;
        dirs[2] = &directionPhiB;
        dirs[3] = &directionPhiISmooth;
        dirs[4] = &directionPhiBSmooth;
        dirs[5] = &directionSmooth;
        dirs[6] = &direction;

        forAll(dirNames, diri)
        {
            const SimpleMapState& dir = *dirs[diri];
            vectorField dY =
                predictorDirection
                (
                    dir,
                    eps,
                    "ATCTPredictor" + dirNames[diri]
                );

            scalar lhs = Zero;
            forAll(dY, celli)
            {
                lhs += randomBarY[celli] & dY[celli];
            }
            reduce(lhs, sumOp<scalar>());

            scalar rhsU = Zero;
            scalar rhsPhiI = Zero;
            scalar rhsPhiB = Zero;
            predictorSeedStateBlocks
            (
                randomPredictorReverse,
                dir,
                rhsU,
                rhsPhiI,
                rhsPhiB
            );

            const scalar rhs = rhsU + rhsPhiI + rhsPhiB;
            const scalar rel =
                mag(lhs - rhs)/max(max(mag(lhs), mag(rhs)), VSMALL);

            Info<< "ATC-T predictor full-state transpose check: "
                << "randomBarY " << dirNames[diri] << " " << eps << " "
                << lhs << " " << rhsU << " " << rhsPhiI << " "
                << rhsPhiB << " " << rhs << " " << rel << " "
                << countMomentumSignCrossings(dir, eps)
                << endl;
        }
    }

    auto runSeedDirection =
        [&]
        (
            const word& seedName,
            const SimpleMapSeed& rawSeed,
            const word& directionName,
            const SimpleMapState& dir,
            const scalar eps
        )
    {
        SimpleMapSeed seedNew(rawSeed);
        projectPressureSeed(seedNew);
        SimpleMapState dNew =
            applyMx(baseState, dir, eps, "ATCTFullState" + seedName + directionName);
        SmoothPhiCoefficientTrace trace(mesh_);
        SimpleMapSeed seedOld =
            reverseOneSimpleMapSeedImpl(seedNew, rAtUBase, &trace);
        SimpleMapSeed seedOldWrapper =
            reverseOneSimpleMapSeed(seedNew, rAtUBase);

        scalar traceSeedDiff = Zero;
        forAll(seedOld.barU, celli)
        {
            traceSeedDiff =
                max(traceSeedDiff, mag(seedOld.barU[celli] - seedOldWrapper.barU[celli]));
            traceSeedDiff =
                max(traceSeedDiff, mag(seedOld.barp[celli] - seedOldWrapper.barp[celli]));
        }
        forAll(seedOld.barPhiInternal, facei)
        {
            traceSeedDiff =
                max
                (
                    traceSeedDiff,
                    mag
                    (
                        seedOld.barPhiInternal[facei]
                      - seedOldWrapper.barPhiInternal[facei]
                    )
                );
        }
        forAll(seedOld.barPhiBoundary, patchi)
        {
            forAll(seedOld.barPhiBoundary[patchi], facei)
            {
                traceSeedDiff =
                    max
                    (
                        traceSeedDiff,
                        mag
                        (
                            seedOld.barPhiBoundary[patchi][facei]
                          - seedOldWrapper.barPhiBoundary[patchi][facei]
                        )
                    );
            }
        }
        reduce(traceSeedDiff, maxOp<scalar>());

        projectPressureSeed(seedOld);

        scalar lhsU = Zero, lhsP = Zero, lhsPhiI = Zero, lhsPhiB = Zero;
        scalar rhsU = Zero, rhsP = Zero, rhsPhiI = Zero, rhsPhiB = Zero;
        seedStateBlocks(seedNew, dNew, lhsU, lhsP, lhsPhiI, lhsPhiB);
        seedStateBlocks(seedOld, dir, rhsU, rhsP, rhsPhiI, rhsPhiB);

        const scalar lhs = lhsU + lhsP + lhsPhiI + lhsPhiB;
        const scalar rhs = rhsU + rhsP + rhsPhiI + rhsPhiB;
        const scalar scale = max(max(mag(lhs), mag(rhs)), VSMALL);
        const scalar absoluteGap = mag(lhs - rhs);
        const scalar symmetricRelativeGap =
            absoluteGap/(mag(lhs) + mag(rhs) + VSMALL);
        const scalar contractionL1Scale =
            mag(lhsU) + mag(lhsP) + mag(lhsPhiI) + mag(lhsPhiB)
          + mag(rhsU) + mag(rhsP) + mag(rhsPhiI) + mag(rhsPhiB)
          + VSMALL;
        const scalar gapOverL1Scale = absoluteGap/contractionL1Scale;

        Info<< "ATC-T full-state map transpose check: "
            << seedName << " "
            << directionName << " "
            << eps << " "
            << lhsU << " "
            << lhsP << " "
            << lhsPhiI << " "
            << lhsPhiB << " "
            << rhsU << " "
            << rhsP << " "
            << rhsPhiI << " "
            << rhsPhiB << " "
            << lhs << " "
            << rhs << " "
            << mag(lhs - rhs)/scale << " "
            << absoluteGap << " "
            << symmetricRelativeGap << " "
            << contractionL1Scale << " "
            << gapOverL1Scale
            << endl;

        if
        (
            directionName == "PhiInternalSmoothOnly"
         || directionName == "PhiBoundarySmoothOnly"
        )
        {
            static bool printedH1Source = false;
            if (!printedH1Source)
            {
                Info<< "ATC-T UEqn.H1 source paths: "
                    << "$FOAM_SRC/finiteVolume/fvMatrices/fvMatrix/"
                    << "fvMatrix.C::fvMatrix<Type>::H1(), "
                    << "$FOAM_SRC/OpenFOAM/matrices/lduMatrix/lduMatrix/"
                    << "lduMatrixATmul.C::lduMatrix::H1()."
                    << endl;
                printedH1Source = true;
            }

            SimpleMapTape dTape =
                makeTapeDirection
                (
                    baseState,
                    dir,
                    eps,
                    "ATCTCoeffTrace" + seedName + directionName
                );

            const scalar fullSignedGap = lhs - rhs;

            const scalar dDNorm = tapeScalarNorm(dTape.Dscalar);
            const scalar dDdiagNorm = tapeScalarNorm(dTape.Ddiag);
            const scalar dDboundaryNorm =
                tapeScalarNorm(dTape.DboundaryTotal);
            const scalar dH1Norm = tapeScalarNorm(dTape.UEqnH1Coeff);
            const scalar drAUNorm = tapeScalarNorm(dTape.rAU);
            const scalar drAtUNorm = tapeScalarNorm(dTape.rAtU);
            const scalar dqNorm = tapeScalarNorm(dTape.q);
            const scalar dDMax = gMax(mag(dTape.Dscalar));
            const scalar dDdiagMax = gMax(mag(dTape.Ddiag));
            const scalar dDboundaryMax = gMax(mag(dTape.DboundaryTotal));
            const scalar dH1Max = gMax(mag(dTape.UEqnH1Coeff));
            const scalar drAUMax = gMax(mag(dTape.rAU));
            const scalar drAtUMax = gMax(mag(dTape.rAtU));
            const scalar dqMax = gMax(mag(dTape.q));

            List<vectorField> zeroYBoundary(mesh_.boundary().size());
            forAll(zeroYBoundary, patchi)
            {
                zeroYBoundary[patchi].setSize
                (
                    mesh_.boundary()[patchi].size(),
                    vector::zero
                );
            }
            vectorField zeroY(mesh_.nCells(), vector::zero);
            scalarField zeroCoeff(mesh_.nCells(), Zero);

            auto traceSurfaceDot =
                [&]
                (
                    const scalarField& aInternal,
                    const List<scalarField>& aBoundary,
                    const scalarField& bInternal,
                    const List<scalarField>& bBoundary
                )
            {
                scalar value = Zero;
                forAll(aInternal, facei)
                {
                    value += aInternal[facei]*bInternal[facei];
                }
                forAll(aBoundary, patchi)
                {
                    forAll(aBoundary[patchi], facei)
                    {
                        value +=
                            aBoundary[patchi][facei]
                           *bBoundary[patchi][facei];
                    }
                }
                reduce(value, sumOp<scalar>());
                return value;
            };

            auto traceGapOverL1 =
                [](const scalar a, const scalar b)
            {
                return mag(a - b)/(mag(a) + mag(b) + VSMALL);
            };

            auto printAdjacentStage =
                [&]
                (
                    const word& fromName,
                    const word& toName,
                    const scalar stageLhs,
                    const scalar stageRhs
                )
            {
                const scalar signedGap = stageLhs - stageRhs;
                Info<< "ATC-T smooth-phi adjacent trace: "
                    << seedName << " " << directionName << " eps " << eps
                    << " " << fromName << "->" << toName
                    << " lhs " << stageLhs
                    << " rhs " << stageRhs
                    << " signedGap " << signedGap
                    << " absoluteGap " << mag(signedGap)
                    << " gapOverL1 " << traceGapOverL1(stageLhs, stageRhs)
                    << endl;
            };

            const scalar dYInternalNorm = tapeVectorNorm(dTape.Y);
            const scalar dYBoundaryNorm =
                tapeBoundaryVectorNorm(dTape.YBoundary);
            const scalar barYBoundaryNorm =
                tapeBoundaryVectorNorm(trace.barYBoundary);

            scalar dFinalBoundaryActualTotal = Zero;
            scalar dFinalBoundaryCurrentTotal = Zero;
            scalar dFinalBoundaryCandidateTotal = Zero;
            scalar dFinalBoundaryCandidateGapL1 = Zero;

            if (directionName == "PhiBoundarySmoothOnly")
            {
                const scalar chi = boundedConvection_ ? scalar(1) : scalar(0);
                scalar alphaU = scalar(1);
                word URelaxName = primalVars_.U().select
                (
                    mesh_.data().isFinalIteration()
                );
                scalar relaxCoeff = scalar(0);
                if (mesh_.relaxEquation(URelaxName, relaxCoeff))
                {
                    alphaU = relaxCoeff;
                }

                struct MomentumRelaxDiagStages
                {
                    scalarField rawDiag;
                    scalarField boundaryMaxTerm;
                    scalarField boundaryMinTerm;
                    scalarField preMaxDiag;
                    scalarField boundaryAugmentedDiag;
                    scalarField sumOffDiag;
                    scalarField dominanceDiag;
                    scalarField maxSelectedTerm;
                    scalarField relaxFormulaDiag;
                    scalarField postMinDiag;
                    scalarField relaxedDiag;
                    scalarField finalDiag;
                    vectorField rawSource;
                    vectorField relaxedSource;
                    vectorField finalSource;
                    List<vectorField> internalCoeffs;

                    explicit MomentumRelaxDiagStages(const fvMesh& mesh)
                    :
                        rawDiag(mesh.nCells(), Zero),
                        boundaryMaxTerm(mesh.nCells(), Zero),
                        boundaryMinTerm(mesh.nCells(), Zero),
                        preMaxDiag(mesh.nCells(), Zero),
                        boundaryAugmentedDiag(mesh.nCells(), Zero),
                        sumOffDiag(mesh.nCells(), Zero),
                        dominanceDiag(mesh.nCells(), Zero),
                        maxSelectedTerm(mesh.nCells(), Zero),
                        relaxFormulaDiag(mesh.nCells(), Zero),
                        postMinDiag(mesh.nCells(), Zero),
                        relaxedDiag(mesh.nCells(), Zero),
                        finalDiag(mesh.nCells(), Zero),
                        rawSource(mesh.nCells(), vector::zero),
                        relaxedSource(mesh.nCells(), vector::zero),
                        finalSource(mesh.nCells(), vector::zero),
                        internalCoeffs(mesh.boundary().size())
                    {
                        forAll(internalCoeffs, patchi)
                        {
                            internalCoeffs[patchi].setSize
                            (
                                mesh.boundary()[patchi].size(),
                                vector::zero
                            );
                        }
                    }
                };

                struct MomentumRawTermStages
                {
                    scalarField diag;
                    List<vectorField> internalCoeffs;

                    explicit MomentumRawTermStages(const fvMesh& mesh)
                    :
                        diag(mesh.nCells(), Zero),
                        internalCoeffs(mesh.boundary().size())
                    {
                        forAll(internalCoeffs, patchi)
                        {
                            internalCoeffs[patchi].setSize
                            (
                                mesh.boundary()[patchi].size(),
                                vector::zero
                            );
                        }
                    }
                };

                auto momentumRelaxDiagStages =
                    [&](const SimpleMapState& state)
                {
                    installState(state);

                    fv::options& localFvOptions(fv::options::New(this->mesh_));
                    tmp<fvVectorMatrix> tStageUEqn
                    (
                        fvm::div(phi, U)
                      + turbulence->divDevReff(U)
                     ==
                        localFvOptions(U)
                    );
                    fvVectorMatrix& stageUEqn = tStageUEqn.ref();

                    MomentumRelaxDiagStages stages(mesh_);
                    stages.rawDiag = stageUEqn.diag();
                    stages.preMaxDiag = stages.rawDiag;
                    stages.rawSource = stageUEqn.source();

                    const scalarField& lower = stageUEqn.lower();
                    const scalarField& upper = stageUEqn.upper();
                    const labelUList& own = mesh_.owner();
                    const labelUList& nei = mesh_.neighbour();
                    forAll(own, facei)
                    {
                        stages.sumOffDiag[nei[facei]] += mag(lower[facei]);
                        stages.sumOffDiag[own[facei]] += mag(upper[facei]);
                    }

                    const FieldField<Field, vector>& iCoeffs =
                        stageUEqn.internalCoeffs();

                    forAll(mesh_.boundary(), patchi)
                    {
                        const fvPatchVectorField& Up =
                            U.boundaryField()[patchi];
                        if (!Up.size())
                        {
                            continue;
                        }

                        if (Up.coupled())
                        {
                            FatalErrorInFunction
                                << "The outlet-Dfinal diagnostic supports "
                                << "non-coupled patches only. Coupled patch "
                                << mesh_.boundary()[patchi].name()
                                << " is active."
                                << exit(FatalError);
                        }

                        const labelUList& faceCells =
                            mesh_.boundary()[patchi].faceCells();
                        forAll(faceCells, facei)
                        {
                            const vector& coeff = iCoeffs[patchi][facei];
                            stages.internalCoeffs[patchi][facei] = coeff;
                            stages.boundaryMaxTerm[faceCells[facei]] +=
                                cmptMax(cmptMag(coeff));
                            stages.boundaryMinTerm[faceCells[facei]] +=
                                cmptMin(coeff);
                        }
                    }

                    forAll(stages.rawDiag, celli)
                    {
                        stages.preMaxDiag[celli] +=
                            stages.boundaryMaxTerm[celli];
                        stages.boundaryAugmentedDiag[celli] =
                            stages.preMaxDiag[celli];
                        stages.dominanceDiag[celli] =
                            max
                            (
                                mag(stages.preMaxDiag[celli]),
                                stages.sumOffDiag[celli]
                            );
                        stages.maxSelectedTerm[celli] =
                            stages.dominanceDiag[celli];
                        stages.postMinDiag[celli] =
                            stages.maxSelectedTerm[celli]/alphaU
                          - stages.boundaryMinTerm[celli];
                        stages.relaxFormulaDiag[celli] =
                            stages.postMinDiag[celli];
                    }

                    stageUEqn.relax();
                    stages.relaxedDiag = stageUEqn.diag();
                    stages.relaxedSource = stageUEqn.source();

                    localFvOptions.constrain(stageUEqn);
                    stages.finalDiag = stageUEqn.diag();
                    stages.finalSource = stageUEqn.source();

                    restoreState();
                    return stages;
                };

                auto momentumRawTermStages =
                    [&](const SimpleMapState& state, const word& termName)
                {
                    installState(state);

                    MomentumRawTermStages stages(mesh_);

                    if (termName == "convection")
                    {
                        tmp<fvVectorMatrix> tTerm(fvm::div(phi, U));
                        fvVectorMatrix& termEqn = tTerm.ref();
                        stages.diag = termEqn.diag();
                        const FieldField<Field, vector>& iCoeffs =
                            termEqn.internalCoeffs();
                        forAll(iCoeffs, patchi)
                        {
                            stages.internalCoeffs[patchi] = iCoeffs[patchi];
                        }
                    }
                    else if (termName == "viscous")
                    {
                        tmp<fvVectorMatrix> tTerm
                        (
                            turbulence->divDevReff(U)
                        );
                        fvVectorMatrix& termEqn = tTerm.ref();
                        stages.diag = termEqn.diag();
                        const FieldField<Field, vector>& iCoeffs =
                            termEqn.internalCoeffs();
                        forAll(iCoeffs, patchi)
                        {
                            stages.internalCoeffs[patchi] = iCoeffs[patchi];
                        }
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Unsupported raw term diagnostic: "
                            << termName
                            << exit(FatalError);
                    }

                    restoreState();
                    return stages;
                };

                auto scalarStageContraction =
                    [&](const scalarField& field)
                {
                    scalar contraction = Zero;
                    forAll(trace.barDfinal, celli)
                    {
                        contraction += trace.barDfinal[celli]*field[celli];
                    }
                    reduce(contraction, sumOp<scalar>());
                    return contraction;
                };

                auto vectorStageNorm =
                    [](const vectorField& field)
                {
                    return gSum(mag(field));
                };

                static bool printedRelaxSource = false;
                if (!printedRelaxSource)
                {
                    Info<< "ATC-T OpenFOAM relax source: "
                        << "$FOAM_SRC/finiteVolume/fvMatrices/fvMatrix/"
                        << "fvMatrix.C::fvMatrix<Type>::relax(alpha), "
                        << "$FOAM_SRC/OpenFOAM/matrices/lduMatrix/lduMatrix/"
                        << "lduMatrixOperations.C::lduMatrix::sumMagOffDiag. "
                        << "Formula: non-coupled patches add "
                        << "cmptMax(cmptMag(internalCoeffs)) to D, "
                        << "D=max(mag(D),sumMagOffDiag), D/=alpha, "
                        << "then subtract cmptMin(internalCoeffs) and add "
                        << "(D-D0)*psi to source."
                        << endl;
                    printedRelaxSource = true;
                }

                typename pTraits<vector>::labelType validVectorComponents
                (
                    mesh_.validComponents<vector>()
                );

                forAll(mesh_.boundary(), patchi)
                {
                    if
                    (
                        mesh_.boundary()[patchi].type() == "empty"
                     || mesh_.boundary()[patchi].coupled()
                    )
                    {
                        continue;
                    }

                    SimpleMapState patchDir(mesh_);
                    patchDir.phiBoundary[patchi] =
                        dir.phiBoundary[patchi];
                    SimpleMapTape dPatchTape =
                        makeTapeDirection
                        (
                            baseState,
                            patchDir,
                            eps,
                            "ATCTBoundaryDfinal"
                          + seedName
                          + directionName
                          + mesh_.boundary()[patchi].name()
                        );

                    const fvPatchVectorField& Up =
                        primalVars_.U().boundaryField()[patchi];
                    const fvsPatchScalarField& phip =
                        primalVars_.phi().boundaryField()[patchi];
                    const labelUList& faceCells =
                        mesh_.boundary()[patchi].faceCells();

                    tmp<vectorField> tValueInternalCoeffs =
                        Up.valueInternalCoeffs
                        (
                            mesh_.boundary()[patchi].weights()
                        );
                    const vectorField& valueInternalCoeffs =
                        tValueInternalCoeffs();

                    scalar patchActual = Zero;
                    forAll(trace.barDfinal, celli)
                    {
                        patchActual +=
                            trace.barDfinal[celli]
                           *dPatchTape.Dfinal[celli];
                    }
                    reduce(patchActual, sumOp<scalar>());

                    SimpleMapState patchPlus(baseState);
                    stateAxpy(patchPlus, eps, patchDir);
                    SimpleMapState patchMinus(baseState);
                    stateAxpy(patchMinus, -eps, patchDir);
                    const MomentumRelaxDiagStages plusStages =
                        momentumRelaxDiagStages(patchPlus);
                    const MomentumRelaxDiagStages minusStages =
                        momentumRelaxDiagStages(patchMinus);

                    const scalarField dRawDiag =
                        (plusStages.rawDiag - minusStages.rawDiag)/(2*eps);
                    const scalarField dBoundaryMaxTerm =
                        (
                            plusStages.boundaryMaxTerm
                          - minusStages.boundaryMaxTerm
                        )/(2*eps);
                    const scalarField dBoundaryMinTerm =
                        (
                            plusStages.boundaryMinTerm
                          - minusStages.boundaryMinTerm
                        )/(2*eps);
                    const scalarField dPreMaxDiag =
                        (
                            plusStages.preMaxDiag
                          - minusStages.preMaxDiag
                        )/(2*eps);
                    const scalarField dBoundaryAugmentedDiag =
                        (
                            plusStages.boundaryAugmentedDiag
                          - minusStages.boundaryAugmentedDiag
                        )/(2*eps);
                    const scalarField dSumOffDiag =
                        (
                            plusStages.sumOffDiag
                          - minusStages.sumOffDiag
                        )/(2*eps);
                    const scalarField dDominanceDiag =
                        (
                            plusStages.dominanceDiag
                          - minusStages.dominanceDiag
                        )/(2*eps);
                    const scalarField dMaxSelectedTerm =
                        (
                            plusStages.maxSelectedTerm
                          - minusStages.maxSelectedTerm
                        )/(2*eps);
                    const scalarField dRelaxFormulaDiag =
                        (
                            plusStages.relaxFormulaDiag
                          - minusStages.relaxFormulaDiag
                        )/(2*eps);
                    const scalarField dPostMinDiag =
                        (
                            plusStages.postMinDiag
                          - minusStages.postMinDiag
                        )/(2*eps);
                    const scalarField dRelaxedDiag =
                        (
                            plusStages.relaxedDiag
                          - minusStages.relaxedDiag
                        )/(2*eps);
                    const scalarField dFinalDiag =
                        (
                            plusStages.finalDiag
                          - minusStages.finalDiag
                        )/(2*eps);
                    const vectorField dRawSource =
                        (
                            plusStages.rawSource
                          - minusStages.rawSource
                        )/(2*eps);
                    const vectorField dRelaxedSource =
                        (
                            plusStages.relaxedSource
                          - minusStages.relaxedSource
                        )/(2*eps);
                    const vectorField dFinalSource =
                        (
                            plusStages.finalSource
                          - minusStages.finalSource
                        )/(2*eps);

                    scalar patchCurrent = Zero;
                    scalar patchCandidate = Zero;
                    scalar coeffMin = GREAT;
                    scalar coeffMax = -GREAT;
                    scalar componentSpreadMax = Zero;
                    label nPositive = 0;
                    label nNegative = 0;

                    forAll(faceCells, facei)
                    {
                        const label celli = faceCells[facei];
                        const scalar phif = phip[facei];
                        if (phif >= scalar(0))
                        {
                            ++nPositive;
                        }
                        else
                        {
                            ++nNegative;
                        }

                        scalar valueCoeff = Zero;
                        bool haveCoeff = false;
                        scalar coeffCmptMin = GREAT;
                        scalar coeffCmptMax = -GREAT;
                        for (label cmpt = 0; cmpt < vector::nComponents; ++cmpt)
                        {
                            if (component(validVectorComponents, cmpt) == -1)
                            {
                                continue;
                            }

                            const scalar c = valueInternalCoeffs[facei][cmpt];
                            if (!haveCoeff)
                            {
                                valueCoeff = c;
                                haveCoeff = true;
                            }
                            coeffCmptMin = min(coeffCmptMin, c);
                            coeffCmptMax = max(coeffCmptMax, c);
                        }

                        if (!haveCoeff)
                        {
                            valueCoeff = Zero;
                            coeffCmptMin = Zero;
                            coeffCmptMax = Zero;
                        }

                        coeffMin = min(coeffMin, valueCoeff);
                        coeffMax = max(coeffMax, valueCoeff);
                        componentSpreadMax =
                            max(componentSpreadMax, mag(coeffCmptMax - coeffCmptMin));

                        const scalar currentCoeff =
                            (
                                (phif >= scalar(0))
                              ? (scalar(1) - chi)
                              : -chi
                            )/alphaU;
                        const scalar candidateCoeff =
                            (
                                (phif >= scalar(0))
                              ? (scalar(1) - chi)
                              : (valueCoeff - chi)
                            )/alphaU;

                        patchCurrent +=
                            trace.barDfinal[celli]
                           *currentCoeff
                           *dir.phiBoundary[patchi][facei];
                        patchCandidate +=
                            trace.barDfinal[celli]
                           *candidateCoeff
                           *dir.phiBoundary[patchi][facei];
                    }

                    reduce(patchCurrent, sumOp<scalar>());
                    reduce(patchCandidate, sumOp<scalar>());
                    reduce(nPositive, sumOp<label>());
                    reduce(nNegative, sumOp<label>());
                    reduce(coeffMin, minOp<scalar>());
                    reduce(coeffMax, maxOp<scalar>());
                    reduce(componentSpreadMax, maxOp<scalar>());

                    dFinalBoundaryActualTotal += patchActual;
                    dFinalBoundaryCurrentTotal += patchCurrent;
                    dFinalBoundaryCandidateTotal += patchCandidate;
                    dFinalBoundaryCandidateGapL1 +=
                        mag(patchActual - patchCandidate);

                    Info<< "ATC-T boundary Dfinal matrix diagnostic: "
                        << seedName << " " << directionName
                        << " eps " << eps
                        << " patch " << mesh_.boundary()[patchi].name()
                        << " UPatchType " << Up.type()
                        << " nPositive " << nPositive
                        << " nNegative " << nNegative
                        << " valueInternalCoeffMin " << coeffMin
                        << " valueInternalCoeffMax " << coeffMax
                        << " valueInternalCoeffComponentSpreadMax "
                        << componentSpreadMax
                        << " actualBarDdotDfinal " << patchActual
                        << " currentAnalytic " << patchCurrent
                        << " candidateAnalytic " << patchCandidate
                        << " currentGap " << patchActual - patchCurrent
                        << " candidateGap " << patchActual - patchCandidate
                        << " candidateRel "
                        << mag(patchActual - patchCandidate)
                          /(mag(patchActual) + mag(patchCandidate) + VSMALL)
                        << endl;

                    const scalar rawContraction =
                        scalarStageContraction(dRawDiag);
                    const scalar boundaryMaxContraction =
                        scalarStageContraction(dBoundaryMaxTerm);
                    const scalar preMaxContraction =
                        scalarStageContraction(dPreMaxDiag);
                    const scalar boundaryAugmentedContraction =
                        scalarStageContraction(dBoundaryAugmentedDiag);
                    const scalar sumOffContraction =
                        scalarStageContraction(dSumOffDiag);
                    const scalar dominanceContraction =
                        scalarStageContraction(dDominanceDiag);
                    const scalar maxSelectedContraction =
                        scalarStageContraction(dMaxSelectedTerm);
                    const scalar boundaryMinContraction =
                        scalarStageContraction(dBoundaryMinTerm);
                    const scalar postMinContraction =
                        scalarStageContraction(dPostMinDiag);
                    const scalar relaxFormulaContraction =
                        scalarStageContraction(dRelaxFormulaDiag);
                    const scalar relaxedContraction =
                        scalarStageContraction(dRelaxedDiag);
                    const scalar finalContraction =
                        scalarStageContraction(dFinalDiag);
                    const scalar relaxFormulaRel =
                        gSum(mag(plusStages.relaxFormulaDiag - plusStages.relaxedDiag))
                       /(
                            gSum(mag(plusStages.relaxFormulaDiag))
                          + gSum(mag(plusStages.relaxedDiag))
                          + VSMALL
                        );
                    const scalar postMinActualRel =
                        gSum(mag(plusStages.postMinDiag - plusStages.relaxedDiag))
                       /(
                            gSum(mag(plusStages.postMinDiag))
                          + gSum(mag(plusStages.relaxedDiag))
                          + VSMALL
                        );

                    Info<< "ATC-T outlet Dfinal relax-stage diagnostic: "
                        << seedName << " " << directionName
                        << " eps " << eps
                        << " patch " << mesh_.boundary()[patchi].name()
                        << " rawDiagContraction " << rawContraction
                        << " boundaryMaxContraction "
                        << boundaryMaxContraction
                        << " preMaxContraction " << preMaxContraction
                        << " boundaryAugmentedContraction "
                        << boundaryAugmentedContraction
                        << " sumOffContraction " << sumOffContraction
                        << " dominanceContraction " << dominanceContraction
                        << " maxSelectedContraction "
                        << maxSelectedContraction
                        << " boundaryMinContraction "
                        << boundaryMinContraction
                        << " postMinContraction " << postMinContraction
                        << " relaxFormulaContraction "
                        << relaxFormulaContraction
                        << " relaxedDiagContraction "
                        << relaxedContraction
                        << " finalDiagContraction " << finalContraction
                        << " fvOptionsDiagContraction "
                        << finalContraction - relaxedContraction
                        << " relaxFormulaVsActualRel "
                        << relaxFormulaRel
                        << " postMinVsActualRel "
                        << postMinActualRel
                        << " rawSourceDerivativeL1 "
                        << vectorStageNorm(dRawSource)
                        << " relaxedSourceDerivativeL1 "
                        << vectorStageNorm(dRelaxedSource)
                        << " finalSourceDerivativeL1 "
                        << vectorStageNorm(dFinalSource)
                        << endl;

                    if (mesh_.boundary()[patchi].name() == "outlet")
                    {
                        bool maxTieHarmless = true;
                        bool componentTieHarmless = true;
                        scalar maxTieWorstRel = Zero;
                        scalar cmptTieWorstDerivativeGap = Zero;

                        const scalar branchTol =
                            scalar(100)*SMALL
                           *max
                            (
                                gMax(mag(plusStages.dominanceDiag)),
                                scalar(1)
                            );
                        forAll(faceCells, facei)
                        {
                            const label celli = faceCells[facei];
                            const scalar diagCandidate =
                                plusStages.boundaryAugmentedDiag[celli];
                            const scalar sumOff =
                                plusStages.sumOffDiag[celli];
                            const scalar branchMargin =
                                mag(diagCandidate) - sumOff;
                            const word activeBranch =
                                (mag(diagCandidate) >= sumOff)
                              ? word("diag")
                              : word("sumOff");
                            if (mag(branchMargin) <= branchTol)
                            {
                                Info<< "ATC-T outlet Dfinal branch warning: "
                                    << seedName << " eps " << eps
                                    << " face " << facei
                                    << " owner " << celli
                                    << " branchMargin " << branchMargin
                                    << " branchTol " << branchTol
                                    << endl;
                            }

                            Info<< "ATC-T outlet Dfinal owner diagnostic: "
                                << seedName << " " << directionName
                                << " eps " << eps
                                << " patch " << mesh_.boundary()[patchi].name()
                                << " face " << facei
                                << " owner " << celli
                                << " phi " << phip[facei]
                                << " rawDiag " << plusStages.rawDiag[celli]
                                << " boundaryMaxTerm "
                                << plusStages.boundaryMaxTerm[celli]
                                << " boundaryMinTerm "
                                << plusStages.boundaryMinTerm[celli]
                                << " preMaxDiag "
                                << plusStages.preMaxDiag[celli]
                                << " boundaryAugmentedDiag "
                                << diagCandidate
                                << " sumOff " << sumOff
                                << " dominance "
                                << plusStages.dominanceDiag[celli]
                                << " relaxedDiag "
                                << plusStages.relaxedDiag[celli]
                                << " finalDiag "
                                << plusStages.finalDiag[celli]
                                << " dRawDiag " << dRawDiag[celli]
                                << " dBoundaryMaxTerm "
                                << dBoundaryMaxTerm[celli]
                                << " dPreMaxDiag " << dPreMaxDiag[celli]
                                << " dBoundaryAugmentedDiag "
                                << dBoundaryAugmentedDiag[celli]
                                << " dSumOff " << dSumOffDiag[celli]
                                << " dMaxSelectedTerm "
                                << dMaxSelectedTerm[celli]
                                << " dBoundaryMinTerm "
                                << dBoundaryMinTerm[celli]
                                << " dPostMinDiag "
                                << dPostMinDiag[celli]
                                << " dDominance " << dDominanceDiag[celli]
                                << " dRelaxFormula "
                                << dRelaxFormulaDiag[celli]
                                << " dRelaxedDiag " << dRelaxedDiag[celli]
                                << " dFinalDiag " << dFinalDiag[celli]
                                << " diagMinusSumOff " << branchMargin
                                << " activeRelaxBranch " << activeBranch
                                << " barD " << trace.barDfinal[celli]
                                << endl;
                        }

                        MomentumRelaxDiagStages baseStages =
                            momentumRelaxDiagStages(baseState);

                        static bool printedConvectionSources = false;
                        if (!printedConvectionSources)
                        {
                            Info<< "ATC-T outlet Dfinal matrix source paths: "
                                << "$FOAM_SRC/finiteVolume/finiteVolume/"
                                << "convectionSchemes/gaussConvectionScheme/"
                                << "gaussConvectionScheme.C::fvmDiv sets "
                                << "lower, upper, negSumDiag and boundary "
                                << "internalCoeffs; "
                                << "$FOAM_SRC/finiteVolume/finiteVolume/"
                                << "convectionSchemes/boundedConvectionScheme/"
                                << "boundedConvectionScheme.C::fvmDiv subtracts "
                                << "fvm::Sp(fvc::surfaceIntegrate(phi), U)."
                                << endl;
                            printedConvectionSources = true;
                        }

                        forAll(faceCells, basisFacei)
                        {
                            const scalar basisAmp = mag(phip[basisFacei]);
                            if (basisAmp <= VSMALL)
                            {
                                continue;
                            }

                            SimpleMapState basisDir(mesh_);
                            basisDir.phiBoundary[patchi][basisFacei] =
                                basisAmp;
                            SimpleMapState basisPlus(baseState);
                            stateAxpy(basisPlus, eps, basisDir);
                            SimpleMapState basisMinus(baseState);
                            stateAxpy(basisMinus, -eps, basisDir);

                            const MomentumRelaxDiagStages basisPlusStages =
                                momentumRelaxDiagStages(basisPlus);
                            const MomentumRelaxDiagStages basisMinusStages =
                                momentumRelaxDiagStages(basisMinus);

                            const MomentumRawTermStages convectionPlus =
                                momentumRawTermStages(basisPlus, "convection");
                            const MomentumRawTermStages convectionMinus =
                                momentumRawTermStages(basisMinus, "convection");
                            const MomentumRawTermStages viscousPlus =
                                momentumRawTermStages(basisPlus, "viscous");
                            const MomentumRawTermStages viscousMinus =
                                momentumRawTermStages(basisMinus, "viscous");

                            const label celli = faceCells[basisFacei];
                            const vector coeffBase =
                                baseStages.internalCoeffs[patchi][basisFacei];
                            const vector dCoeff =
                                (
                                    basisPlusStages.internalCoeffs[patchi][basisFacei]
                                  - basisMinusStages.internalCoeffs[patchi][basisFacei]
                                )/(2*eps);
                            const vector dConvCoeff =
                                (
                                    convectionPlus.internalCoeffs[patchi][basisFacei]
                                  - convectionMinus.internalCoeffs[patchi][basisFacei]
                                )/(2*eps);
                            const vector dViscCoeff =
                                (
                                    viscousPlus.internalCoeffs[patchi][basisFacei]
                                  - viscousMinus.internalCoeffs[patchi][basisFacei]
                                )/(2*eps);

                            const scalar bRaw =
                                (
                                    basisPlusStages.rawDiag[celli]
                                  - basisMinusStages.rawDiag[celli]
                                )/(2*eps);
                            const scalar bBoundaryMax =
                                (
                                    basisPlusStages.boundaryMaxTerm[celli]
                                  - basisMinusStages.boundaryMaxTerm[celli]
                                )/(2*eps);
                            const scalar bPreMax =
                                (
                                    basisPlusStages.preMaxDiag[celli]
                                  - basisMinusStages.preMaxDiag[celli]
                                )/(2*eps);
                            const scalar bSumOff =
                                (
                                    basisPlusStages.sumOffDiag[celli]
                                  - basisMinusStages.sumOffDiag[celli]
                                )/(2*eps);
                            const scalar bMaxSelected =
                                (
                                    basisPlusStages.maxSelectedTerm[celli]
                                  - basisMinusStages.maxSelectedTerm[celli]
                                )/(2*eps);
                            const scalar bBoundaryMin =
                                (
                                    basisPlusStages.boundaryMinTerm[celli]
                                  - basisMinusStages.boundaryMinTerm[celli]
                                )/(2*eps);
                            const scalar bPostMin =
                                (
                                    basisPlusStages.postMinDiag[celli]
                                  - basisMinusStages.postMinDiag[celli]
                                )/(2*eps);
                            const scalar bActualRelaxed =
                                (
                                    basisPlusStages.relaxedDiag[celli]
                                  - basisMinusStages.relaxedDiag[celli]
                                )/(2*eps);
                            const scalar bConvRawDiag =
                                (
                                    convectionPlus.diag[celli]
                                  - convectionMinus.diag[celli]
                                )/(2*eps);
                            const scalar bViscRawDiag =
                                (
                                    viscousPlus.diag[celli]
                                  - viscousMinus.diag[celli]
                                )/(2*eps);

                            const word plusSelected =
                                (
                                    mag(basisPlusStages.preMaxDiag[celli])
                                 >= basisPlusStages.sumOffDiag[celli]
                                )
                              ? word("preMax")
                              : word("sumOff");
                            const word minusSelected =
                                (
                                    mag(basisMinusStages.preMaxDiag[celli])
                                 >= basisMinusStages.sumOffDiag[celli]
                                )
                              ? word("preMax")
                              : word("sumOff");

                            const scalar maxTieRel =
                                mag(bPreMax - bSumOff)
                               /(mag(bPreMax) + mag(bSumOff) + VSMALL);
                            maxTieWorstRel =
                                max(maxTieWorstRel, maxTieRel);
                            if (maxTieRel > scalar(1e-10))
                            {
                                maxTieHarmless = false;
                            }

                            scalar minVal = GREAT;
                            scalar maxMagVal = -GREAT;
                            label nMin = 0;
                            label nMaxMag = 0;
                            scalar minDerivRef = Zero;
                            scalar maxMagDerivRef = Zero;
                            bool haveMinDeriv = false;
                            bool haveMaxMagDeriv = false;

                            for (label cmpt = 0; cmpt < vector::nComponents; ++cmpt)
                            {
                                if (component(validVectorComponents, cmpt) == -1)
                                {
                                    continue;
                                }
                                minVal = min(minVal, coeffBase[cmpt]);
                                maxMagVal = max(maxMagVal, mag(coeffBase[cmpt]));
                            }

                            const scalar componentTieTol =
                                scalar(100)*SMALL
                               *max(max(mag(minVal), mag(maxMagVal)), scalar(1));
                            for (label cmpt = 0; cmpt < vector::nComponents; ++cmpt)
                            {
                                if (component(validVectorComponents, cmpt) == -1)
                                {
                                    continue;
                                }

                                if (mag(coeffBase[cmpt] - minVal) <= componentTieTol)
                                {
                                    ++nMin;
                                    if (!haveMinDeriv)
                                    {
                                        minDerivRef = dCoeff[cmpt];
                                        haveMinDeriv = true;
                                    }
                                    else
                                    {
                                        const scalar gap =
                                            mag(dCoeff[cmpt] - minDerivRef);
                                        cmptTieWorstDerivativeGap =
                                            max(cmptTieWorstDerivativeGap, gap);
                                        if (gap > scalar(1e-10))
                                        {
                                            componentTieHarmless = false;
                                        }
                                    }
                                }

                                if (mag(mag(coeffBase[cmpt]) - maxMagVal) <= componentTieTol)
                                {
                                    ++nMaxMag;
                                    const scalar maxMagDeriv =
                                        sign(coeffBase[cmpt])*dCoeff[cmpt];
                                    if (!haveMaxMagDeriv)
                                    {
                                        maxMagDerivRef = maxMagDeriv;
                                        haveMaxMagDeriv = true;
                                    }
                                    else
                                    {
                                        const scalar gap =
                                            mag(maxMagDeriv - maxMagDerivRef);
                                        cmptTieWorstDerivativeGap =
                                            max(cmptTieWorstDerivativeGap, gap);
                                        if (gap > scalar(1e-10))
                                        {
                                            componentTieHarmless = false;
                                        }
                                    }
                                }
                            }

                            const scalar formulaGap =
                                bActualRelaxed
                              - (bMaxSelected/alphaU - bBoundaryMin);

                            Info<< "ATC-T outlet Dfinal face-basis diagnostic: "
                                << seedName << " eps " << eps
                                << " patch " << mesh_.boundary()[patchi].name()
                                << " face " << basisFacei
                                << " owner " << celli
                                << " basisAmp " << basisAmp
                                << " dRawDiag " << bRaw
                                << " dBoundaryMaxTerm " << bBoundaryMax
                                << " dPreMaxDiag " << bPreMax
                                << " dSumOffDiag " << bSumOff
                                << " dMaxSelectedTerm " << bMaxSelected
                                << " dBoundaryMinTerm " << bBoundaryMin
                                << " dPostMinDiag " << bPostMin
                                << " dActualRelaxedDiag " << bActualRelaxed
                                << " dMaxOverAlpha " << bMaxSelected/alphaU
                                << " minusDBoundaryMin " << -bBoundaryMin
                                << " formulaGap " << formulaGap
                                << " plusSelected " << plusSelected
                                << " minusSelected " << minusSelected
                                << " maxTieRel " << maxTieRel
                                << " coeffBase " << coeffBase
                                << " dCoeff " << dCoeff
                                << " validComponents "
                                << validVectorComponents
                                << " minVal " << minVal
                                << " nMinComponents " << nMin
                                << " maxMagVal " << maxMagVal
                                << " nMaxMagComponents " << nMaxMag
                                << " convectionRawDiagDerivative "
                                << bConvRawDiag
                                << " viscousRawDiagDerivative "
                                << bViscRawDiag
                                << " fullRawDiagDerivative " << bRaw
                                << " convectionInternalCoeffDerivative "
                                << dConvCoeff
                                << " viscousInternalCoeffDerivative "
                                << dViscCoeff
                                << endl;
                        }

                        Info<< "ATC-T outlet Dfinal differentiability decision: "
                            << seedName << " " << directionName
                            << " eps " << eps
                            << " patch " << mesh_.boundary()[patchi].name()
                            << " maxTieHarmless " << maxTieHarmless
                            << " maxTieWorstRel " << maxTieWorstRel
                            << " componentTieHarmless "
                            << componentTieHarmless
                            << " componentTieWorstDerivativeGap "
                            << cmptTieWorstDerivativeGap
                            << " uniqueClassicalDerivative "
                            << (maxTieHarmless && componentTieHarmless)
                            << endl;
                    }
                }

                reduce(dFinalBoundaryActualTotal, sumOp<scalar>());
                reduce(dFinalBoundaryCurrentTotal, sumOp<scalar>());
                reduce(dFinalBoundaryCandidateTotal, sumOp<scalar>());
                reduce(dFinalBoundaryCandidateGapL1, sumOp<scalar>());
            }

            forAll(mesh_.boundary(), patchi)
            {
                const scalar patchDYBoundaryL1 =
                    gSum(mag(dTape.YBoundary[patchi]));
                const scalar patchBarYBoundaryL1 =
                    gSum(mag(trace.barYBoundary[patchi]));
                const scalar patchYBoundaryDot =
                    tapeBoundaryVectorDot
                    (
                        List<vectorField>(1, trace.barYBoundary[patchi]),
                        List<vectorField>(1, dTape.YBoundary[patchi])
                    );

                if
                (
                    patchDYBoundaryL1 > SMALL
                 || patchBarYBoundaryL1 > SMALL
                 || mag(patchYBoundaryDot) > SMALL
                )
                {
                    Info<< "ATC-T smooth-phi actual-Y patch: "
                        << seedName << " " << directionName
                        << " eps " << eps
                        << " patch " << mesh_.boundary()[patchi].name()
                        << " UPatchType "
                        << primalVars_.U().boundaryField()[patchi].type()
                        << " dYBoundaryL1 " << patchDYBoundaryL1
                        << " barYBoundaryL1 " << patchBarYBoundaryL1
                        << " barYBoundaryDotDYBoundary "
                        << patchYBoundaryDot
                        << endl;
                }
            }

            PredictorReverseSeed dReverse =
                reversePredictorYFullState
                (
                    zeroY,
                    zeroYBoundary,
                    trace.barDfinal,
                    zeroCoeff,
                    true,
                    true
                );
            scalar dRhsU = Zero;
            scalar dRhsPhiI = Zero;
            scalar dRhsPhiB = Zero;
            predictorSeedStateBlocks(dReverse, dir, dRhsU, dRhsPhiI, dRhsPhiB);
            const scalar lhsDdiag =
                tapeScalarDot(trace.barDfinal, dTape.Ddiag);
            const scalar lhsDboundary =
                tapeScalarDot(trace.barDfinal, dTape.DboundaryTotal);
            const scalar lhsD =
                tapeScalarDot(trace.barDfinal, dTape.Dscalar);
            const scalar rhsD = dRhsU + dRhsPhiI + dRhsPhiB;
            const scalar gapD = lhsD - rhsD;
            const scalar scaleD =
                mag(lhsD) + mag(rhsD) + mag(dRhsU)
              + mag(dRhsPhiI) + mag(dRhsPhiB) + VSMALL;

            scalar alphaUForDiagCandidate = scalar(1);
            word URelaxNameForDiagCandidate =
                primalVars_.U().select(mesh_.data().isFinalIteration());
            scalar relaxCoeffForDiagCandidate = scalar(0);
            if
            (
                mesh_.relaxEquation
                (
                    URelaxNameForDiagCandidate,
                    relaxCoeffForDiagCandidate
                )
            )
            {
                alphaUForDiagCandidate = relaxCoeffForDiagCandidate;
            }
            const FinalDiagonalPhiSeed diagCandidate =
                reverseExtraFinalDiagonalToPhi
                (
                    trace.barDfinal,
                    alphaUForDiagCandidate
                );
            const scalar rhsDiagCandidate =
                traceSurfaceDot
                (
                    diagCandidate.internal,
                    diagCandidate.boundary,
                    dir.phiInternal,
                    dir.phiBoundary
                );

            PredictorReverseSeed h1Reverse =
                reversePredictorYFullState
                (
                    zeroY,
                    zeroYBoundary,
                    zeroCoeff,
                    trace.barUEqnH1Coeff,
                    true,
                    true
                );
            scalar h1RhsU = Zero;
            scalar h1RhsPhiI = Zero;
            scalar h1RhsPhiB = Zero;
            predictorSeedStateBlocks(h1Reverse, dir, h1RhsU, h1RhsPhiI, h1RhsPhiB);
            const scalar lhsH1 =
                tapeScalarDot(trace.barUEqnH1Coeff, dTape.UEqnH1Coeff);
            const scalar rhsH1 = h1RhsU + h1RhsPhiI + h1RhsPhiB;
            const scalar gapH1 = lhsH1 - rhsH1;
            const scalar scaleH1 =
                mag(lhsH1) + mag(rhsH1) + mag(h1RhsU)
              + mag(h1RhsPhiI) + mag(h1RhsPhiB) + VSMALL;

            PredictorReverseSeed yReverse =
                reversePredictorYFullState
                (
                    trace.barYInternal,
                    trace.barYBoundary,
                    zeroCoeff,
                    zeroCoeff,
                    true,
                    true
                );
            scalar yRhsU = Zero;
            scalar yRhsPhiI = Zero;
            scalar yRhsPhiB = Zero;
            predictorSeedStateBlocks(yReverse, dir, yRhsU, yRhsPhiI, yRhsPhiB);
            const scalar lhsYInternal =
                tapeVectorDot(trace.barYInternal, dTape.Y);
            const scalar lhsYBoundary =
                tapeBoundaryVectorDot(trace.barYBoundary, dTape.YBoundary);
            const scalar lhsY = lhsYInternal + lhsYBoundary;
            const scalar rhsY = yRhsU + yRhsPhiI + yRhsPhiB;
            const scalar gapY = lhsY - rhsY;
            const scalar scaleY =
                mag(lhsY) + mag(rhsY) + mag(yRhsU)
              + mag(yRhsPhiI) + mag(yRhsPhiB) + VSMALL;

            forAll(mesh_.boundary(), patchi)
            {
                scalar lhsYPatch = Zero;
                scalar rhsYPatch = Zero;
                forAll(trace.barYBoundary[patchi], facei)
                {
                    lhsYPatch +=
                        trace.barYBoundary[patchi][facei]
                      & dTape.YBoundary[patchi][facei];
                }
                forAll(yReverse.barPhiBoundary[patchi], facei)
                {
                    rhsYPatch +=
                        yReverse.barPhiBoundary[patchi][facei]
                       *dir.phiBoundary[patchi][facei];
                }
                reduce(lhsYPatch, sumOp<scalar>());
                reduce(rhsYPatch, sumOp<scalar>());

                if (mag(lhsYPatch) > SMALL || mag(rhsYPatch) > SMALL)
                {
                    Info<< "ATC-T smooth-phi actual-Y patch contribution: "
                        << seedName << " " << directionName
                        << " eps " << eps
                        << " patch " << mesh_.boundary()[patchi].name()
                        << " lhsYBoundaryPatch " << lhsYPatch
                        << " rhsYBoundaryPhiPatch " << rhsYPatch
                        << endl;
                }
            }

            PredictorReverseSeed combinedPredictorReverse =
                reversePredictorYFullState
                (
                    trace.barYInternal,
                    trace.barYBoundary,
                    trace.barDfinal,
                    trace.barUEqnH1Coeff,
                    true,
                    true
                );
            scalar combinedRhsU = Zero;
            scalar combinedRhsPhiI = Zero;
            scalar combinedRhsPhiB = Zero;
            predictorSeedStateBlocks
            (
                combinedPredictorReverse,
                dir,
                combinedRhsU,
                combinedRhsPhiI,
                combinedRhsPhiB
            );
            const scalar lhsPredictor = lhsY + lhsD + lhsH1;
            const scalar rhsPredictor =
                combinedRhsU + combinedRhsPhiI + combinedRhsPhiB;
            const scalar gapPredictorCombined =
                lhsPredictor - rhsPredictor;
            const scalar rhsPredictorLinear =
                rhsY + rhsD + rhsH1;
            const scalar predictorLinearityGap =
                rhsPredictor - rhsPredictorLinear;
            const scalar predictorLinearityRel =
                mag(predictorLinearityGap)
               /(
                    mag(rhsPredictor)
                  + mag(rhsPredictorLinear)
                  + VSMALL
                );

            const scalarField barpNew(collapsedFinalPressureSeed(seedNew));
            scalarField barpSolve(barpNew);
            barpSolve *= pRelaxCoeffTape;

            const scalar barF2TraceDotF2 =
                traceSurfaceDot
                (
                    trace.barF2Internal,
                    trace.barF2Boundary,
                    dTape.F2Internal,
                    dTape.F2Boundary
                );
            const scalar barF1TraceDotF1 =
                traceSurfaceDot
                (
                    trace.barF1Internal,
                    trace.barF1Boundary,
                    dTape.F1Internal,
                    dTape.F1Boundary
                );
            const scalar barF0TraceDotF0 =
                traceSurfaceDot
                (
                    trace.barF0Internal,
                    trace.barF0Boundary,
                    dTape.F0Internal,
                    dTape.F0Boundary
                );
            const scalar barH0DirectContraction =
                tapeVectorDot(trace.barH0DirectInternal, dTape.H0)
              + tapeBoundaryVectorDot
                (
                    trace.barH0DirectBoundary,
                    dTape.H0Boundary
                );
            const scalar barH0FromFluxContraction =
                tapeVectorDot(trace.barH0FromFluxInternal, dTape.H0)
              + tapeBoundaryVectorDot
                (
                    trace.barH0FromFluxBoundary,
                    dTape.H0Boundary
                );
            const scalar barH0TotalContraction =
                tapeVectorDot(trace.barH0TotalInternal, dTape.H0)
              + tapeBoundaryVectorDot
                (
                    trace.barH0TotalBoundary,
                    dTape.H0Boundary
                );
            const scalar barH0SplitGap =
                barH0TotalContraction
              - barH0DirectContraction
              - barH0FromFluxContraction;

            scalarField barqSplitCurrent(trace.barqFromH1);
            barqSplitCurrent += trace.barqFromF2InternalFaces;
            barqSplitCurrent += trace.barqFromF2BoundaryCurrent;
            scalar qSplitMaxDiff = Zero;
            forAll(trace.barq, celli)
            {
                qSplitMaxDiff =
                    max(qSplitMaxDiff, mag(trace.barq[celli] - barqSplitCurrent[celli]));
            }
            reduce(qSplitMaxDiff, maxOp<scalar>());

            scalarField barqZeroBoundary(trace.barqFromH1);
            barqZeroBoundary += trace.barqFromF2InternalFaces;

            const scalar lhsQH1 =
                tapeScalarDot(trace.barqFromH1, dTape.q);
            scalar lhsQInternal = Zero;
            forAll(trace.barQFaceInternal, facei)
            {
                lhsQInternal +=
                    trace.barQFaceInternal[facei]*dTape.qFaceInternal[facei];
            }
            reduce(lhsQInternal, sumOp<scalar>());

            scalar lhsQBoundary = Zero;
            scalar rhsQBoundaryCurrent = Zero;
            scalar qBoundarySeedL1 = Zero;
            scalar qBoundaryForwardL1 = Zero;
            forAll(mesh_.boundary(), patchi)
            {
                const labelUList& faceCells =
                    mesh_.boundary()[patchi].faceCells();
                forAll(trace.barQFaceBoundary[patchi], facei)
                {
                    const scalar seedFace =
                        trace.barQFaceBoundary[patchi][facei];
                    lhsQBoundary +=
                        seedFace*dTape.qFaceBoundary[patchi][facei];
                    rhsQBoundaryCurrent +=
                        seedFace*dTape.q[faceCells[facei]];
                    qBoundarySeedL1 += mag(seedFace);
                    qBoundaryForwardL1 +=
                        mag(dTape.qFaceBoundary[patchi][facei]);
                }
            }
            reduce(lhsQBoundary, sumOp<scalar>());
            reduce(rhsQBoundaryCurrent, sumOp<scalar>());
            reduce(qBoundarySeedL1, sumOp<scalar>());
            reduce(qBoundaryForwardL1, sumOp<scalar>());

            const scalar lhsQ =
                lhsQH1 + lhsQInternal + lhsQBoundary;
            const scalar rhsQCurrent =
                tapeScalarDot(barqSplitCurrent, dTape.q);
            const scalar rhsQZeroBoundary =
                tapeScalarDot(barqZeroBoundary, dTape.q);
            const scalar gapQBoundaryCurrent = lhsQ - rhsQCurrent;
            const scalar gapQBoundaryZero = lhsQ - rhsQZeroBoundary;
            const scalar qCurrentRel =
                mag(gapQBoundaryCurrent)
               /(mag(lhsQ) + mag(rhsQCurrent) + VSMALL);
            const scalar qZeroRel =
                mag(gapQBoundaryZero)
               /(mag(lhsQ) + mag(rhsQZeroBoundary) + VSMALL);

            volScalarField qBoundaryProbe
            (
                IOobject
                (
                    "ATCTQBoundaryProbe",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rAtUBase
            );
            qBoundaryProbe.primitiveFieldRef() = baseTapeForStages.q;
            forAll(qBoundaryProbe.boundaryFieldRef(), patchi)
            {
                qBoundaryProbe.boundaryFieldRef()[patchi] ==
                    baseTapeForStages.qBoundary[patchi];
            }
            surfaceScalarField qFaceBoundaryProbe
            (
                IOobject
                (
                    "ATCTQFaceBoundaryProbe",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(qBoundaryProbe)
            );

            forAll(mesh_.boundary(), patchi)
            {
                if
                (
                    mesh_.boundary()[patchi].type() == "empty"
                 || mesh_.boundary()[patchi].coupled()
                )
                {
                    continue;
                }

                const labelUList& faceCells =
                    mesh_.boundary()[patchi].faceCells();
                scalar dqOwnerL1 = Zero;
                scalar dqBoundaryL1 = Zero;
                scalar dqFaceBoundaryL1 = Zero;
                scalar maxDQBoundaryMinusOwner = Zero;
                scalar maxDQFaceMinusOwner = Zero;
                scalar maxDQFaceMinusBoundary = Zero;
                scalar patchForward = Zero;
                scalar patchCurrent = Zero;
                scalar patchZero = Zero;
                scalar patchSeedL1 = Zero;

                forAll(faceCells, facei)
                {
                    const label celli = faceCells[facei];
                    const scalar dqOwner = dTape.q[celli];
                    const scalar dqBoundary =
                        dTape.qBoundary[patchi][facei];
                    const scalar dqFaceBoundary =
                        dTape.qFaceBoundary[patchi][facei];
                    const scalar qFaceSeed =
                        trace.barQFaceBoundary[patchi][facei];

                    dqOwnerL1 += mag(dqOwner);
                    dqBoundaryL1 += mag(dqBoundary);
                    dqFaceBoundaryL1 += mag(dqFaceBoundary);
                    maxDQBoundaryMinusOwner =
                        max(maxDQBoundaryMinusOwner, mag(dqBoundary - dqOwner));
                    maxDQFaceMinusOwner =
                        max(maxDQFaceMinusOwner, mag(dqFaceBoundary - dqOwner));
                    maxDQFaceMinusBoundary =
                        max
                        (
                            maxDQFaceMinusBoundary,
                            mag(dqFaceBoundary - dqBoundary)
                        );
                    patchForward += qFaceSeed*dqFaceBoundary;
                    patchCurrent += qFaceSeed*dqOwner;
                    patchSeedL1 += mag(qFaceSeed);
                }

                reduce(dqOwnerL1, sumOp<scalar>());
                reduce(dqBoundaryL1, sumOp<scalar>());
                reduce(dqFaceBoundaryL1, sumOp<scalar>());
                reduce(maxDQBoundaryMinusOwner, maxOp<scalar>());
                reduce(maxDQFaceMinusOwner, maxOp<scalar>());
                reduce(maxDQFaceMinusBoundary, maxOp<scalar>());
                reduce(patchForward, sumOp<scalar>());
                reduce(patchCurrent, sumOp<scalar>());
                reduce(patchSeedL1, sumOp<scalar>());

                if
                (
                    dqOwnerL1 > SMALL
                 || dqBoundaryL1 > SMALL
                 || dqFaceBoundaryL1 > SMALL
                 || patchSeedL1 > SMALL
                 || mag(patchForward) > SMALL
                 || mag(patchCurrent) > SMALL
                )
                {
                    const fvPatchScalarField& qp =
                        qBoundaryProbe.boundaryField()[patchi];
                    const fvsPatchScalarField& qfp =
                        qFaceBoundaryProbe.boundaryField()[patchi];
                    Info<< "ATC-T consistent-q boundary tangent: "
                        << seedName << " " << directionName
                        << " eps " << eps
                        << " patch " << mesh_.boundary()[patchi].name()
                        << " qPatchType " << qp.type()
                        << " qFacePatchType " << qfp.type()
                        << " qAssignable " << qp.assignable()
                        << " qFixesValue " << qp.fixesValue()
                        << " dqOwnerL1 " << dqOwnerL1
                        << " dqBoundaryL1 " << dqBoundaryL1
                        << " dqFaceBoundaryL1 " << dqFaceBoundaryL1
                        << " maxDQBoundaryMinusOwner "
                        << maxDQBoundaryMinusOwner
                        << " maxDQFaceMinusOwner "
                        << maxDQFaceMinusOwner
                        << " maxDQFaceMinusBoundary "
                        << maxDQFaceMinusBoundary
                        << " barQFaceBoundaryDotDQFaceBoundary "
                        << patchForward
                        << " currentOwnerContraction "
                        << patchCurrent
                        << " zeroBoundaryContraction " << patchZero
                        << " patchSeedL1 " << patchSeedL1
                        << endl;
                }
            }

            if (directionName == "PhiBoundarySmoothOnly")
            {
                forAll(mesh_.boundary(), patchi)
                {
                    if
                    (
                        mesh_.boundary()[patchi].type() == "empty"
                     || mesh_.boundary()[patchi].coupled()
                    )
                    {
                        continue;
                    }

                    const labelUList& faceCells =
                        mesh_.boundary()[patchi].faceCells();
                    const fvsPatchScalarField& phiPatch =
                        primalVars_.phi().boundaryField()[patchi];
                    const fvPatchScalarField& qp =
                        qBoundaryProbe.boundaryField()[patchi];
                    const fvsPatchScalarField& qfp =
                        qFaceBoundaryProbe.boundaryField()[patchi];

                    forAll(faceCells, facei)
                    {
                        const scalar qFaceSeed =
                            trace.barQFaceBoundary[patchi][facei];
                        if (mag(qFaceSeed) <= SMALL)
                        {
                            continue;
                        }

                        const scalar basisAmp =
                            max(mag(phiPatch[facei]), scalar(1));
                        SimpleMapState qBasisDir(mesh_);
                        qBasisDir.phiBoundary[patchi][facei] = basisAmp;
                        SimpleMapTape qBasisTape =
                            makeTapeDirection
                            (
                                baseState,
                                qBasisDir,
                                eps,
                                "ATCTQBoundaryBasis"
                              + seedName
                              + directionName
                              + mesh_.boundary()[patchi].name()
                              + Foam::name(facei)
                            );

                        const label celli = faceCells[facei];
                        const scalar dqOwner = qBasisTape.q[celli];
                        const scalar dqBoundary =
                            qBasisTape.qBoundary[patchi][facei];
                        const scalar dqFaceBoundary =
                            qBasisTape.qFaceBoundary[patchi][facei];
                        const scalar exactForward =
                            qFaceSeed*dqFaceBoundary;
                        const scalar currentOwner =
                            qFaceSeed*dqOwner;

                        Info<< "ATC-T consistent-q boundary face-basis: "
                            << seedName << " " << directionName
                            << " eps " << eps
                            << " patch " << mesh_.boundary()[patchi].name()
                            << " face " << facei
                            << " owner " << celli
                            << " qPatchType " << qp.type()
                            << " qFacePatchType " << qfp.type()
                            << " basePhi " << phiPatch[facei]
                            << " basisAmp " << basisAmp
                            << " qFaceSeed " << qFaceSeed
                            << " dqOwner " << dqOwner
                            << " dqBoundary " << dqBoundary
                            << " dqFaceBoundary " << dqFaceBoundary
                            << " currentMappingContraction "
                            << currentOwner
                            << " zeroBoundaryContraction " << scalar(0)
                            << " exactForwardContraction "
                            << exactForward
                            << " currentGap "
                            << exactForward - currentOwner
                            << " zeroGap " << exactForward
                            << endl;
                    }
                }
            }

            const scalar lhsCons =
                tapeVectorDot(seedNew.barU, dTape.H1)
              + barF2TraceDotF2;
            const scalar rhsCons =
                barH0DirectContraction
              + tapeScalarDot(trace.barq, dTape.q);

            const scalar rhsConsTotal = rhsCons + barF1TraceDotF1;
            const scalar gapCons = lhsCons - rhsConsTotal;
            const scalar scaleCons =
                mag(lhsCons) + mag(rhsConsTotal) + VSMALL;

            const scalar lhsAdjust = barF1TraceDotF1;
            const scalar rhsAdjust = barF0TraceDotF0;
            const scalar gapAdjust = lhsAdjust - rhsAdjust;

            const scalar lhsFlux = rhsAdjust;
            const scalar rhsFlux = barH0FromFluxContraction;
            const scalar gapFlux = lhsFlux - rhsFlux;

            const scalar lhsConstrain = barH0TotalContraction;
            const scalar rhsConstrain =
                lhsY
              + tapeBoundaryVectorDot(trace.barUBoundary, dTape.UoldBoundary);
            const scalar gapConstrain = lhsConstrain - rhsConstrain;

            const scalar hPressure = eps;
            const SimpleMapTape pressurePlus =
                pressureStageAtRAtUDirection
                (
                    dTape.rAtU,
                    hPressure,
                    scalar(1),
                    "ATCTPressureCoeffPlus" + seedName + directionName
                );
            const SimpleMapTape pressureMinus =
                pressureStageAtRAtUDirection
                (
                    dTape.rAtU,
                    hPressure,
                    scalar(-1),
                    "ATCTPressureCoeffMinus" + seedName + directionName
                );
            SimpleMapTape dPressure =
                tapeDifference(pressurePlus, pressureMinus, hPressure);

            const scalar pSolveCoeff =
                tapeScalarDot(barpSolve, dPressure.pSolve);
            scalar phiCoeffI = Zero;
            forAll(seedNew.barPhiInternal, facei)
            {
                phiCoeffI +=
                    seedNew.barPhiInternal[facei]
                   *dPressure.phiNewInternal[facei];
            }
            scalar phiCoeffB = Zero;
            forAll(seedNew.barPhiBoundary, patchi)
            {
                forAll(seedNew.barPhiBoundary[patchi], facei)
                {
                    phiCoeffB +=
                        seedNew.barPhiBoundary[patchi][facei]
                       *dPressure.phiNewBoundary[patchi][facei];
                }
            }
            reduce(phiCoeffI, sumOp<scalar>());
            reduce(phiCoeffB, sumOp<scalar>());
            const scalar lhsPressureCoeff =
                pSolveCoeff + phiCoeffI + phiCoeffB;
            const scalar rhsPressureCoeff =
                tapeScalarDot(trace.barrAtUFromPressureStage, dTape.rAtU);
            const scalar gapPressure =
                lhsPressureCoeff - rhsPressureCoeff;
            const scalar scalePressure =
                mag(lhsPressureCoeff) + mag(rhsPressureCoeff) + VSMALL;

            const PressureFluxProbeTape frozenPatchFluxPlus =
                pressureFluxProbeAtRAtUDirection
                (
                    dTape.rAtU,
                    hPressure,
                    scalar(1),
                    false,
                    "ATCTPressureFixedPPlus" + seedName + directionName
                );
            const PressureFluxProbeTape frozenPatchFluxMinus =
                pressureFluxProbeAtRAtUDirection
                (
                    dTape.rAtU,
                    hPressure,
                    scalar(-1),
                    false,
                    "ATCTPressureFixedPMinus" + seedName + directionName
                );
            const PressureFluxProbeTape reconstrainedFluxPlus =
                pressureFluxProbeAtRAtUDirection
                (
                    dTape.rAtU,
                    hPressure,
                    scalar(1),
                    true,
                    "ATCTPressureReconstrainedPlus" + seedName + directionName
                );
            const PressureFluxProbeTape reconstrainedFluxMinus =
                pressureFluxProbeAtRAtUDirection
                (
                    dTape.rAtU,
                    hPressure,
                    scalar(-1),
                    true,
                    "ATCTPressureReconstrainedMinus" + seedName + directionName
                );
            PressureFluxProbeTape dFrozenPatchFlux =
                pressureFluxProbeDifference
                (
                    frozenPatchFluxPlus,
                    frozenPatchFluxMinus,
                    hPressure
                );
            PressureFluxProbeTape dReconstrainedFlux =
                pressureFluxProbeDifference
                (
                    reconstrainedFluxPlus,
                    reconstrainedFluxMinus,
                    hPressure
                );

            const Pair<scalar> frozenTotal =
                pressureFluxSeedContraction
                (
                    trace.barF2Internal,
                    trace.barF2Boundary,
                    dFrozenPatchFlux.totalInternal,
                    dFrozenPatchFlux.totalBoundary
                );
            const Pair<scalar> frozenCoefficient =
                pressureFluxSeedContraction
                (
                    trace.barF2Internal,
                    trace.barF2Boundary,
                    dFrozenPatchFlux.coefficientInternal,
                    dFrozenPatchFlux.coefficientBoundary
                );
            const Pair<scalar> frozenCorrection =
                pressureFluxSeedContraction
                (
                    trace.barF2Internal,
                    trace.barF2Boundary,
                    dFrozenPatchFlux.correctionInternal,
                    dFrozenPatchFlux.correctionBoundary
                );
            const Pair<scalar> reconstrainedTotal =
                pressureFluxSeedContraction
                (
                    trace.barF2Internal,
                    trace.barF2Boundary,
                    dReconstrainedFlux.totalInternal,
                    dReconstrainedFlux.totalBoundary
                );

            scalar fixedPFluxInternal = reconstrainedTotal.first();
            scalar fixedPFluxBoundary = reconstrainedTotal.second();
            const scalar fixedPFluxTotal =
                fixedPFluxInternal + fixedPFluxBoundary;
            const scalar pressureCoefficientOnly =
                frozenTotal.first() + frozenTotal.second();
            const scalar pressureBoundaryStateEffect =
                fixedPFluxTotal - pressureCoefficientOnly;
            const scalar pressureCoefficientMatrixFlux =
                frozenCoefficient.first() + frozenCoefficient.second();
            const scalar pressureFaceFluxCorrection =
                frozenCorrection.first() + frozenCorrection.second();
            const scalar pressureFluxDecompositionGap =
                pressureCoefficientOnly
              - pressureCoefficientMatrixFlux
              - pressureFaceFluxCorrection;
            const scalar fixedPIdentityGap =
                lhsPressureCoeff - fixedPFluxTotal;
            const scalar fixedPReverseGap =
                fixedPFluxTotal - rhsPressureCoeff;

            const surfaceScalarField gammaPlus =
                gammaFaceAtRAtUDirection
                (
                    dTape.rAtU,
                    hPressure,
                    scalar(1),
                    "ATCTPressureGammaPlus" + seedName + directionName
                );
            const surfaceScalarField gammaMinus =
                gammaFaceAtRAtUDirection
                (
                    dTape.rAtU,
                    hPressure,
                    scalar(-1),
                    "ATCTPressureGammaMinus" + seedName + directionName
                );

            scalarField dGammaInternal(mesh_.nInternalFaces(), Zero);
            forAll(dGammaInternal, facei)
            {
                dGammaInternal[facei] =
                    (
                        gammaPlus.primitiveField()[facei]
                      - gammaMinus.primitiveField()[facei]
                    )/(2*hPressure);
            }
            List<scalarField> dGammaBoundary(mesh_.boundary().size());
            forAll(dGammaBoundary, patchi)
            {
                dGammaBoundary[patchi].setSize
                (
                    mesh_.boundary()[patchi].size(),
                    Zero
                );
                forAll(dGammaBoundary[patchi], facei)
                {
                    dGammaBoundary[patchi][facei] =
                        (
                            gammaPlus.boundaryField()[patchi][facei]
                          - gammaMinus.boundaryField()[patchi][facei]
                        )/(2*hPressure);
                }
            }

            scalarField barGammaInternal(mesh_.nInternalFaces(), Zero);
            List<scalarField> barGammaBoundary(mesh_.boundary().size());
            forAll(barGammaBoundary, patchi)
            {
                barGammaBoundary[patchi].setSize
                (
                    mesh_.boundary()[patchi].size(),
                    Zero
                );
            }

            scalar faceRhsInternal = Zero;
            scalar faceRhsBoundary = Zero;
            scalar faceLhsInternal = frozenTotal.first();
            scalar faceLhsBoundary = frozenTotal.second();
            forAll(barGammaInternal, facei)
            {
                barGammaInternal[facei] =
                   -trace.barF2Internal[facei]
                   *pressureRAtUUnitFlux.primitiveField()[facei];
                faceRhsInternal += barGammaInternal[facei]*dGammaInternal[facei];
            }
            forAll(barGammaBoundary, patchi)
            {
                forAll(barGammaBoundary[patchi], facei)
                {
                    barGammaBoundary[patchi][facei] =
                       -trace.barF2Boundary[patchi][facei]
                       *pressureRAtUUnitFlux.boundaryField()[patchi][facei];
                    faceRhsBoundary +=
                        barGammaBoundary[patchi][facei]
                       *dGammaBoundary[patchi][facei];
                }
            }
            reduce(faceRhsInternal, sumOp<scalar>());
            reduce(faceRhsBoundary, sumOp<scalar>());

            scalar internalOrientADiff = Zero;
            scalar internalOrientBDiff = Zero;
            scalar internalOrientScale = VSMALL;
            const surfaceScalarField& w = mesh_.weights();
            forAll(mesh_.owner(), facei)
            {
                const label P = mesh_.owner()[facei];
                const label N = mesh_.neighbour()[facei];
                const scalar candA =
                    w[facei]*dTape.rAtU[P]
                  + (scalar(1) - w[facei])*dTape.rAtU[N];
                const scalar candB =
                    (scalar(1) - w[facei])*dTape.rAtU[P]
                  + w[facei]*dTape.rAtU[N];
                internalOrientADiff += mag(candA - dGammaInternal[facei]);
                internalOrientBDiff += mag(candB - dGammaInternal[facei]);
                internalOrientScale += mag(dGammaInternal[facei]);
            }
            reduce(internalOrientADiff, sumOp<scalar>());
            reduce(internalOrientBDiff, sumOp<scalar>());
            reduce(internalOrientScale, sumOp<scalar>());

            const bool useOrientationA =
                internalOrientADiff <= internalOrientBDiff;
            scalarField barrAtUInternalFaces(mesh_.nCells(), Zero);
            scalarField barrAtUBoundaryFaces(mesh_.nCells(), Zero);
            scalarField barrAtUDirectOwnerVariant(mesh_.nCells(), Zero);
            scalarField barrAtUOmitBoundaryVariant(mesh_.nCells(), Zero);

            forAll(mesh_.owner(), facei)
            {
                const label P = mesh_.owner()[facei];
                const label N = mesh_.neighbour()[facei];
                const scalar ownerW =
                    useOrientationA ? w[facei] : scalar(1) - w[facei];
                const scalar neiW =
                    useOrientationA ? scalar(1) - w[facei] : w[facei];
                barrAtUInternalFaces[P] += ownerW*barGammaInternal[facei];
                barrAtUInternalFaces[N] += neiW*barGammaInternal[facei];
                barrAtUOmitBoundaryVariant[P] += ownerW*barGammaInternal[facei];
                barrAtUOmitBoundaryVariant[N] += neiW*barGammaInternal[facei];
                barrAtUDirectOwnerVariant[P] += ownerW*barGammaInternal[facei];
                barrAtUDirectOwnerVariant[N] += neiW*barGammaInternal[facei];
            }

            scalar boundaryOwnerMapUnresolved = Zero;
            scalar boundaryPatchContraction = Zero;
            scalar boundaryOwnerContraction = Zero;
            scalar boundaryOmittedContraction = Zero;
            scalar directOwnerVariantContraction = Zero;

            forAll(mesh_.boundary(), patchi)
            {
                const fvPatch& patch = mesh_.boundary()[patchi];
                scalar drOwnerL1 = Zero;
                scalar drBoundaryL1 = Zero;
                scalar dGammaBoundaryL1 = Zero;
                scalar maxBoundaryOwnerDiff = Zero;
                scalar patchBarGammaDotDGamma = Zero;
                scalar patchOwnerDotDR = Zero;
                scalar patchDirectOwnerDotDR = Zero;

                const labelUList& faceCells = patch.faceCells();
                forAll(patch, facei)
                {
                    const label celli = faceCells[facei];
                    const scalar drOwner = dTape.rAtU[celli];
                    const scalar drBoundary = dTape.rAtUBoundary[patchi][facei];
                    const scalar dGamma = dGammaBoundary[patchi][facei];
                    drOwnerL1 += mag(drOwner);
                    drBoundaryL1 += mag(drBoundary);
                    dGammaBoundaryL1 += mag(dGamma);
                    maxBoundaryOwnerDiff =
                        max(maxBoundaryOwnerDiff, mag(dGamma - drOwner));
                    patchBarGammaDotDGamma +=
                        barGammaBoundary[patchi][facei]*dGamma;
                    barrAtUDirectOwnerVariant[celli] +=
                        barGammaBoundary[patchi][facei];
                    patchDirectOwnerDotDR +=
                        barGammaBoundary[patchi][facei]*drOwner;
                }

                const scalar ownerRel =
                    maxBoundaryOwnerDiff/(dGammaBoundaryL1 + drOwnerL1 + VSMALL);
                const bool zeroBoundaryTangent =
                    dGammaBoundaryL1 <= SMALL;
                const bool mapsToOwner =
                    !zeroBoundaryTangent && ownerRel < scalar(1e-10);
                if (mapsToOwner)
                {
                    forAll(patch, facei)
                    {
                        const label celli = faceCells[facei];
                        barrAtUBoundaryFaces[celli] +=
                            barGammaBoundary[patchi][facei];
                        patchOwnerDotDR +=
                            barGammaBoundary[patchi][facei]*dTape.rAtU[celli];
                    }
                }
                else if (dGammaBoundaryL1 > SMALL)
                {
                    boundaryOwnerMapUnresolved += dGammaBoundaryL1;
                }

                boundaryPatchContraction += patchBarGammaDotDGamma;
                boundaryOwnerContraction += patchOwnerDotDR;
                directOwnerVariantContraction += patchDirectOwnerDotDR;

                reduce(drOwnerL1, sumOp<scalar>());
                reduce(drBoundaryL1, sumOp<scalar>());
                reduce(dGammaBoundaryL1, sumOp<scalar>());
                reduce(maxBoundaryOwnerDiff, maxOp<scalar>());
                reduce(patchBarGammaDotDGamma, sumOp<scalar>());
                reduce(patchOwnerDotDR, sumOp<scalar>());
                reduce(patchDirectOwnerDotDR, sumOp<scalar>());

                if
                (
                    drOwnerL1 > SMALL
                 || drBoundaryL1 > SMALL
                 || dGammaBoundaryL1 > SMALL
                 || mag(patchBarGammaDotDGamma) > SMALL
                )
                {
                    const fvPatchScalarField& rPatch =
                        rAtUBase.boundaryField()[patchi];
                    Info<< "ATC-T pressure-rAtU boundary interpolation: "
                        << seedName << " " << directionName
                        << " eps " << eps
                        << " patch " << patch.name()
                        << " rAtUType " << rPatch.type()
                        << " assignable " << rPatch.assignable()
                        << " fixesValue " << rPatch.fixesValue()
                        << " drOwnerL1 " << drOwnerL1
                        << " drBoundaryL1 " << drBoundaryL1
                        << " dGammaBoundaryL1 " << dGammaBoundaryL1
                        << " maxDGammaMinusOwner "
                        << maxBoundaryOwnerDiff
                        << " ownerRel " << ownerRel
                        << " zeroBoundaryTangent " << zeroBoundaryTangent
                        << " mapsToOwner " << mapsToOwner
                        << " barGammaDotDGamma "
                        << patchBarGammaDotDGamma
                        << " ownerMappedContraction "
                        << patchOwnerDotDR
                        << " directOwnerContraction "
                        << patchDirectOwnerDotDR
                        << endl;
                }
            }
            reduce(boundaryOwnerMapUnresolved, sumOp<scalar>());
            reduce(boundaryPatchContraction, sumOp<scalar>());
            reduce(boundaryOwnerContraction, sumOp<scalar>());
            reduce(boundaryOmittedContraction, sumOp<scalar>());
            reduce(directOwnerVariantContraction, sumOp<scalar>());

            scalar internalCellContribution =
                tapeScalarDot(barrAtUInternalFaces, dTape.rAtU);
            scalar boundaryCellContribution =
                tapeScalarDot(barrAtUBoundaryFaces, dTape.rAtU);
            scalarField barrAtUTotalFaces =
                barrAtUInternalFaces + barrAtUBoundaryFaces;
            scalar totalCellContribution =
                tapeScalarDot(barrAtUTotalFaces, dTape.rAtU);
            scalar directOwnerCellContribution =
                tapeScalarDot(barrAtUDirectOwnerVariant, dTape.rAtU);
            scalar omittedBoundaryContribution =
                tapeScalarDot(barrAtUOmitBoundaryVariant, dTape.rAtU);
            const scalar faceSpaceLhs =
                faceLhsInternal + faceLhsBoundary;
            const scalar faceSpaceRhs =
                faceRhsInternal + faceRhsBoundary;
            const scalar faceSpaceScale =
                mag(faceSpaceLhs) + mag(faceSpaceRhs) + VSMALL;
            const scalar cellSpaceScale =
                mag(faceSpaceLhs) + mag(totalCellContribution) + VSMALL;

            forAll(mesh_.boundary(), patchi)
            {
                scalar patchFixedP = Zero;
                scalar patchFluxNorm = Zero;
                scalar patchRNorm = Zero;
                forAll(trace.barF2Boundary[patchi], facei)
                {
                    patchFixedP -=
                        trace.barF2Boundary[patchi][facei]
                       *dReconstrainedFlux.totalBoundary[patchi][facei];
                    patchFluxNorm +=
                        mag(dReconstrainedFlux.totalBoundary[patchi][facei]);
                    patchRNorm +=
                        mag(dTape.rAtU[mesh_.boundary()[patchi].faceCells()[facei]]);
                }
                reduce(patchFixedP, sumOp<scalar>());
                reduce(patchFluxNorm, sumOp<scalar>());
                reduce(patchRNorm, sumOp<scalar>());

                if
                (
                    mag(patchFixedP) > SMALL
                 || patchFluxNorm > SMALL
                 || patchRNorm > SMALL
                )
                {
                    const fvPatchScalarField& pp =
                        primalVars_.p().boundaryField()[patchi];
                    Info<< "ATC-T pressure-rAtU fixedP patch: "
                        << seedName << " " << directionName
                        << " eps " << eps
                        << " patch " << mesh_.boundary()[patchi].name()
                        << " pType " << pp.type()
                        << " updateableSnGrad "
                        << isA<updateablePatchTypes::updateableSnGrad>(pp)
                        << " fixesValue " << pp.fixesValue()
                        << " fixedPContraction " << patchFixedP
                        << " dPressureFluxFixedPNorm " << patchFluxNorm
                        << " drAtUOwnerBoundaryNorm " << patchRNorm
                        << endl;
                }
            }

            const scalar finalRAtUContribution =
                tapeScalarDot(trace.barrAtUFromFinalCorrection, dTape.rAtU);
            const scalar barUBoundaryContraction =
                tapeBoundaryVectorDot(trace.barUBoundary, dTape.UoldBoundary);
            const scalar barUFromBoundaryContraction =
                tapeVectorDot(trace.barUFromBoundary, dir.U);
            const scalar coefficientInputsCurrent = lhsD + lhsH1;
            const scalar coefficientInputsWithPressureFD =
                coefficientInputsCurrent + lhsPressureCoeff;

            auto printCoefficientIdentity =
                [&]
                (
                    const word& identityName,
                    const scalarField& predicted,
                    const scalarField& actual
                )
            {
                scalar signedDiff = Zero;
                scalar diffL1 = Zero;
                scalar scaleL1 = VSMALL;
                scalar maxAbsDiff = Zero;
                label maxCell = -1;

                forAll(predicted, celli)
                {
                    const scalar residual = predicted[celli] - actual[celli];
                    const scalar absResidual = mag(residual);
                    signedDiff += residual;
                    diffL1 += absResidual;
                    scaleL1 += max(mag(predicted[celli]), mag(actual[celli]));
                    if (absResidual > maxAbsDiff)
                    {
                        maxAbsDiff = absResidual;
                        maxCell = celli;
                    }
                }

                reduce(signedDiff, sumOp<scalar>());
                reduce(diffL1, sumOp<scalar>());
                reduce(scaleL1, sumOp<scalar>());
                reduce(maxAbsDiff, maxOp<scalar>());
                reduce(maxCell, maxOp<label>());

                Info<< "ATC-T coefficient forward identity: "
                    << seedName << " " << directionName
                    << " eps " << eps
                    << " identity " << identityName
                    << " signedDiff " << signedDiff
                    << " relativeL1 " << diffL1/scaleL1
                    << " maxAbsDiff " << maxAbsDiff
                    << " maxCell " << maxCell
                    << endl;
            };

            auto printCoefficientSeedCompare =
                [&]
                (
                    const word& seedLabel,
                    const scalarField& actualSeed,
                    const scalarField& expectedSeed
                )
            {
                scalar diffL1 = Zero;
                scalar scaleL1 = VSMALL;
                scalar maxAbsDiff = Zero;
                label maxCell = -1;

                forAll(actualSeed, celli)
                {
                    const scalar residual =
                        actualSeed[celli] - expectedSeed[celli];
                    const scalar absResidual = mag(residual);
                    diffL1 += absResidual;
                    scaleL1 +=
                        max(mag(actualSeed[celli]), mag(expectedSeed[celli]));
                    if (absResidual > maxAbsDiff)
                    {
                        maxAbsDiff = absResidual;
                        maxCell = celli;
                    }
                }

                reduce(diffL1, sumOp<scalar>());
                reduce(scaleL1, sumOp<scalar>());
                reduce(maxAbsDiff, maxOp<scalar>());
                reduce(maxCell, maxOp<label>());

                Info<< "ATC-T coefficient seed algebra: "
                    << seedName << " " << directionName
                    << " eps " << eps
                    << " seed " << seedLabel
                    << " relativeL1 " << diffL1/scaleL1
                    << " maxAbsDiff " << maxAbsDiff
                    << " maxCell " << maxCell
                    << endl;
            };

            scalarField drAUFromDiag(mesh_.nCells(), Zero);
            scalarField drAUFromD(mesh_.nCells(), Zero);
            scalarField drAUFromA(mesh_.nCells(), Zero);
            scalarField drAUFromATimesV(mesh_.nCells(), Zero);
            scalarField drAtUFromDiag(mesh_.nCells(), Zero);
            scalarField drAtUFromD(mesh_.nCells(), Zero);
            scalarField drAtUFromA(mesh_.nCells(), Zero);
            scalarField drAtUFromATimesV(mesh_.nCells(), Zero);
            scalarField dqFromDiag(mesh_.nCells(), Zero);
            scalarField dqFromD(mesh_.nCells(), Zero);
            scalarField dqFromA(mesh_.nCells(), Zero);
            scalarField dqFromATimesV(mesh_.nCells(), Zero);

            forAll(drAUFromDiag, celli)
            {
                const scalar rAUCell = baseTapeForStages.rAU[celli];
                const scalar rAtUCell = baseTapeForStages.rAtU[celli];
                const scalar VCell = mesh_.V()[celli];
                drAUFromDiag[celli] =
                    -sqr(rAUCell)/VCell*dTape.Ddiag[celli];
                drAUFromD[celli] =
                    -sqr(rAUCell)/VCell*dTape.Dscalar[celli];
                drAUFromA[celli] =
                    -sqr(rAUCell)*dTape.Ainternal[celli];
                drAUFromATimesV[celli] =
                    -sqr(rAUCell)/VCell*dTape.AinternalTimesV[celli];

                drAtUFromDiag[celli] =
                    sqr(rAtUCell)/sqr(rAUCell)*drAUFromDiag[celli]
                  + sqr(rAtUCell)*dTape.UEqnH1Coeff[celli];
                drAtUFromD[celli] =
                    sqr(rAtUCell)/sqr(rAUCell)*drAUFromD[celli]
                  + sqr(rAtUCell)*dTape.UEqnH1Coeff[celli];
                drAtUFromA[celli] =
                    sqr(rAtUCell)/sqr(rAUCell)*drAUFromA[celli]
                  + sqr(rAtUCell)*dTape.UEqnH1Coeff[celli];
                drAtUFromATimesV[celli] =
                    sqr(rAtUCell)/sqr(rAUCell)*drAUFromATimesV[celli]
                  + sqr(rAtUCell)*dTape.UEqnH1Coeff[celli];

                dqFromDiag[celli] =
                    drAtUFromDiag[celli] - drAUFromDiag[celli];
                dqFromD[celli] =
                    drAtUFromD[celli] - drAUFromD[celli];
                dqFromA[celli] =
                    drAtUFromA[celli] - drAUFromA[celli];
                dqFromATimesV[celli] =
                    drAtUFromATimesV[celli] - drAUFromATimesV[celli];
            }

            printCoefficientIdentity("drAUFromDdiag", drAUFromDiag, dTape.rAU);
            printCoefficientIdentity("drAUFromD", drAUFromD, dTape.rAU);
            printCoefficientIdentity("drAUFromA", drAUFromA, dTape.rAU);
            printCoefficientIdentity
            (
                "drAUFromATimesV",
                drAUFromATimesV,
                dTape.rAU
            );
            printCoefficientIdentity
            (
                "drAtUFromDdiag",
                drAtUFromDiag,
                dTape.rAtU
            );
            printCoefficientIdentity("drAtUFromD", drAtUFromD, dTape.rAtU);
            printCoefficientIdentity("drAtUFromA", drAtUFromA, dTape.rAtU);
            printCoefficientIdentity
            (
                "drAtUFromATimesV",
                drAtUFromATimesV,
                dTape.rAtU
            );
            printCoefficientIdentity("dqFromDdiag", dqFromDiag, dTape.q);
            printCoefficientIdentity("dqFromD", dqFromD, dTape.q);
            printCoefficientIdentity("dqFromA", dqFromA, dTape.q);
            printCoefficientIdentity("dqFromATimesV", dqFromATimesV, dTape.q);

            scalarField barR0Expected(trace.barrAtUFromFinalCorrection);
            barR0Expected += trace.barrAtUFromPressureStage;
            scalarField barRExpected(barR0Expected);
            barRExpected += trace.barq;
            scalarField barAExpected(mesh_.nCells(), Zero);
            scalarField barH1Expected(mesh_.nCells(), Zero);
            scalarField barDExpected(mesh_.nCells(), Zero);
            scalarField barAInternalExpected(mesh_.nCells(), Zero);
            forAll(barAExpected, celli)
            {
                const scalar rAUCell = trace.rAUUsed[celli];
                const scalar rAtUCell = trace.rAtUUsed[celli];
                const scalar VCell = mesh_.V()[celli];
                barAExpected[celli] =
                    -trace.barq[celli]
                  + barRExpected[celli]*sqr(rAtUCell)/sqr(rAUCell);
                barH1Expected[celli] =
                    barRExpected[celli]*sqr(rAtUCell);
                barDExpected[celli] =
                    -barAExpected[celli]*sqr(rAUCell)/VCell;
                barAInternalExpected[celli] =
                    -barAExpected[celli]*sqr(rAUCell);
            }

            printCoefficientSeedCompare
            (
                "barrAtUBeforeQ",
                trace.barrAtUBeforeQ,
                barR0Expected
            );
            printCoefficientSeedCompare
            (
                "barrAtUAfterQ",
                trace.barrAtUAfterQ,
                barRExpected
            );
            printCoefficientSeedCompare
            (
                "barrAUAfterRAtUReverse",
                trace.barrAUAfterRAtUReverse,
                barAExpected
            );
            printCoefficientSeedCompare
            (
                "barUEqnH1Coeff",
                trace.barUEqnH1Coeff,
                barH1Expected
            );
            printCoefficientSeedCompare
            (
                "barDfinal",
                trace.barDfinal,
                barDExpected
            );

            {
                scalar relRAUUsed = Zero, maxRAUUsed = Zero;
                scalar relRAtUUsed = Zero, maxRAtUUsed = Zero;
                scalarFieldIdentityMetrics
                (
                    trace.rAUUsed,
                    baseTapeForStages.rAU,
                    relRAUUsed,
                    maxRAUUsed
                );
                scalarFieldIdentityMetrics
                (
                    trace.rAtUUsed,
                    baseTapeForStages.rAtU,
                    relRAtUUsed,
                    maxRAtUUsed
                );
                Info<< "ATC-T coefficient reverse base comparison: "
                    << seedName << " " << directionName
                    << " eps " << eps
                    << " relRAUUsedVsTape " << relRAUUsed
                    << " maxRAUUsedVsTape " << maxRAUUsed
                    << " relRAtUUsedVsTape " << relRAtUUsed
                    << " maxRAtUUsedVsTape " << maxRAtUUsed
                    << endl;
            }

            auto coefficientSeedCase =
                [&]
                (
                    const word& caseName,
                    const scalarField& barR0Case,
                    const scalarField& barQCase,
                    scalar& lhsCoeff,
                    scalar& rhsDiag,
                    scalar& rhsD,
                    scalar& rhsA,
                    scalar& rhsATimesV
                )
            {
                scalarField barRCase(barR0Case);
                barRCase += barQCase;
                scalarField barACase(mesh_.nCells(), Zero);
                scalarField barH1Case(mesh_.nCells(), Zero);
                scalarField barDCase(mesh_.nCells(), Zero);
                scalarField barAInternalCase(mesh_.nCells(), Zero);

                forAll(barACase, celli)
                {
                    const scalar rAUCell = trace.rAUUsed[celli];
                    const scalar rAtUCell = trace.rAtUUsed[celli];
                    const scalar VCell = mesh_.V()[celli];
                    barACase[celli] =
                        -barQCase[celli]
                      + barRCase[celli]*sqr(rAtUCell)/sqr(rAUCell);
                    barH1Case[celli] =
                        barRCase[celli]*sqr(rAtUCell);
                    barDCase[celli] =
                        -barACase[celli]*sqr(rAUCell)/VCell;
                    barAInternalCase[celli] =
                        -barACase[celli]*sqr(rAUCell);
                }

                lhsCoeff =
                    tapeScalarDot(barR0Case, dTape.rAtU)
                  + tapeScalarDot(barQCase, dTape.q);
                rhsDiag =
                    tapeScalarDot(barDCase, dTape.Ddiag)
                  + tapeScalarDot(barH1Case, dTape.UEqnH1Coeff);
                rhsD =
                    tapeScalarDot(barDCase, dTape.Dscalar)
                  + tapeScalarDot(barH1Case, dTape.UEqnH1Coeff);
                rhsA =
                    tapeScalarDot(barAInternalCase, dTape.Ainternal)
                  + tapeScalarDot(barH1Case, dTape.UEqnH1Coeff);
                rhsATimesV =
                    tapeScalarDot(barDCase, dTape.AinternalTimesV)
                  + tapeScalarDot(barH1Case, dTape.UEqnH1Coeff);

                Info<< "ATC-T coefficient seed-source contraction: "
                    << seedName << " " << directionName
                    << " eps " << eps
                    << " case " << caseName
                    << " lhsCoeff " << lhsCoeff
                    << " rhsCoeffDdiag " << rhsDiag
                    << " gapDdiag " << lhsCoeff - rhsDiag
                    << " relDdiag " << traceGapOverL1(lhsCoeff, rhsDiag)
                    << " rhsCoeffD " << rhsD
                    << " gapD " << lhsCoeff - rhsD
                    << " relD " << traceGapOverL1(lhsCoeff, rhsD)
                    << " rhsCoeffA " << rhsA
                    << " gapA " << lhsCoeff - rhsA
                    << " relA " << traceGapOverL1(lhsCoeff, rhsA)
                    << " rhsCoeffATimesV " << rhsATimesV
                    << " gapATimesV " << lhsCoeff - rhsATimesV
                    << " relATimesV "
                    << traceGapOverL1(lhsCoeff, rhsATimesV)
                    << endl;
            };

            scalarField zeroCellCoeff(mesh_.nCells(), Zero);
            scalar lhsCoeffFinal = Zero, rhsCoeffFinalDiag = Zero;
            scalar rhsCoeffFinalD = Zero, rhsCoeffFinalA = Zero;
            scalar rhsCoeffFinalATimesV = Zero;
            coefficientSeedCase
            (
                "finalRAtUOnly",
                trace.barrAtUFromFinalCorrection,
                zeroCellCoeff,
                lhsCoeffFinal,
                rhsCoeffFinalDiag,
                rhsCoeffFinalD,
                rhsCoeffFinalA,
                rhsCoeffFinalATimesV
            );

            scalar lhsCoeffPressure = Zero, rhsCoeffPressureDiag = Zero;
            scalar rhsCoeffPressureD = Zero, rhsCoeffPressureA = Zero;
            scalar rhsCoeffPressureATimesV = Zero;
            coefficientSeedCase
            (
                "pressureRAtUOnly",
                trace.barrAtUFromPressureStage,
                zeroCellCoeff,
                lhsCoeffPressure,
                rhsCoeffPressureDiag,
                rhsCoeffPressureD,
                rhsCoeffPressureA,
                rhsCoeffPressureATimesV
            );

            scalar lhsCoeffQ = Zero, rhsCoeffQDiag = Zero;
            scalar rhsCoeffQD = Zero, rhsCoeffQA = Zero;
            scalar rhsCoeffQATimesV = Zero;
            coefficientSeedCase
            (
                "qOnly",
                zeroCellCoeff,
                trace.barq,
                lhsCoeffQ,
                rhsCoeffQDiag,
                rhsCoeffQD,
                rhsCoeffQA,
                rhsCoeffQATimesV
            );

            scalar lhsCoeffCombined = Zero, rhsCoeffCombinedDiag = Zero;
            scalar rhsCoeffCombinedD = Zero, rhsCoeffCombinedA = Zero;
            scalar rhsCoeffCombinedATimesV = Zero;
            coefficientSeedCase
            (
                "combinedProduction",
                barR0Expected,
                trace.barq,
                lhsCoeffCombined,
                rhsCoeffCombinedDiag,
                rhsCoeffCombinedD,
                rhsCoeffCombinedA,
                rhsCoeffCombinedATimesV
            );

            const scalar c0Output = lhs;
            const scalar c1AfterFinal =
                tapeVectorDot(seedNew.barU, dTape.H1)
              + tapeScalarDot(trace.barpNew, dTape.pNew)
              + tapePhiDot(seedNew, dTape)
              + finalRAtUContribution;
            const scalar c2AfterRelax =
                tapeVectorDot(seedNew.barU, dTape.H1)
              + tapeScalarDot(trace.barpSolve, dTape.pSolve)
              + tapePhiDot(seedNew, dTape)
              + finalRAtUContribution;
            const scalar c3AfterPressureFixed =
                tapeVectorDot(seedNew.barU, dTape.H1)
              + traceSurfaceDot
                (
                    trace.barF2Internal,
                    trace.barF2Boundary,
                    dTape.F2Internal,
                    dTape.F2Boundary
                )
              + finalRAtUContribution;
            const scalar c4AfterPressureCoeffFD =
                c3AfterPressureFixed + lhsPressureCoeff;
            const scalar c4ActualPressureReverse =
                c3AfterPressureFixed + rhsPressureCoeff;
            const scalar c5AfterConsistent =
                barH0DirectContraction
              + barF1TraceDotF1
              + tapeScalarDot(trace.barq, dTape.q)
              + finalRAtUContribution;
            const scalar c5AfterConsistentWithPressureFD =
                c5AfterConsistent + lhsPressureCoeff;
            const scalar c5AfterConsistentActual =
                c5AfterConsistent + rhsPressureCoeff;
            const scalar c6AfterCoeffAlgebra =
                barH0DirectContraction
              + barF0TraceDotF0
              + coefficientInputsCurrent;
            const scalar c6AfterCoeffAlgebraWithPressureFD =
                barH0DirectContraction
              + barF0TraceDotF0
              + coefficientInputsWithPressureFD;
            const scalar c7AfterFlux =
                barH0TotalContraction
              + coefficientInputsCurrent;
            const scalar c7AfterFluxWithPressureFD =
                barH0TotalContraction
              + coefficientInputsWithPressureFD;
            const scalar c7AfterPredictorInputs =
                lhsPredictor + barUBoundaryContraction;
            const scalar c8AfterBoundaryEval =
                lhsPredictor + barUFromBoundaryContraction;
            const scalar c9AfterPredictorReverse =
                rhsPredictor + barUFromBoundaryContraction;
            const scalar c10OldState = rhs;

            const scalar gapCoeffDiag =
                lhsCoeffCombined - rhsCoeffCombinedDiag;
            const scalar gapCoeffD =
                lhsCoeffCombined - rhsCoeffCombinedD;
            const scalar gapCoeffA =
                lhsCoeffCombined - rhsCoeffCombinedA;
            const scalar gapCoeffATimesV =
                lhsCoeffCombined - rhsCoeffCombinedATimesV;

            word dExactName("Ddiag");
            const scalarField* dDExactPtr = &dTape.Ddiag;
            scalar gapCoeffExact = gapCoeffDiag;
            if (mag(gapCoeffD) < mag(gapCoeffExact))
            {
                dExactName = "D";
                dDExactPtr = &dTape.Dscalar;
                gapCoeffExact = gapCoeffD;
            }
            if (mag(gapCoeffATimesV) < mag(gapCoeffExact))
            {
                dExactName = "ATimesV";
                dDExactPtr = &dTape.AinternalTimesV;
                gapCoeffExact = gapCoeffATimesV;
            }
            if (mag(gapCoeffA) < mag(gapCoeffExact))
            {
                dExactName = "A";
                dDExactPtr = &dTape.Ainternal;
                gapCoeffExact = gapCoeffA;
            }
            const scalarField& dDExact = *dDExactPtr;

            const scalar c5ActualMinusC6 =
                c5AfterConsistentActual - c6AfterCoeffAlgebra;
            const scalar c5c6ReconDiag = gapAdjust + gapCoeffDiag;
            const scalar c5c6ReconD = gapAdjust + gapCoeffD;
            const scalar c5c6ReconA = gapAdjust + gapCoeffA;
            const scalar c5c6ReconATimesV = gapAdjust + gapCoeffATimesV;
            const scalar c5c6ReconExact = gapAdjust + gapCoeffExact;

            Info<< "ATC-T coefficient C5C6 decomposition: "
                << seedName << " " << directionName
                << " eps " << eps
                << " c5ActualMinusC6 " << c5ActualMinusC6
                << " gapAdjust " << gapAdjust
                << " gapCoeffDdiag " << gapCoeffDiag
                << " reconDdiag " << c5c6ReconDiag
                << " reconResidualDdiag "
                << c5ActualMinusC6 - c5c6ReconDiag
                << " gapCoeffD " << gapCoeffD
                << " reconD " << c5c6ReconD
                << " reconResidualD "
                << c5ActualMinusC6 - c5c6ReconD
                << " gapCoeffA " << gapCoeffA
                << " reconA " << c5c6ReconA
                << " reconResidualA "
                << c5ActualMinusC6 - c5c6ReconA
                << " gapCoeffATimesV " << gapCoeffATimesV
                << " reconATimesV " << c5c6ReconATimesV
                << " reconResidualATimesV "
                << c5ActualMinusC6 - c5c6ReconATimesV
                << " selectedDExact " << dExactName
                << " selectedGapCoeff " << gapCoeffExact
                << " selectedRecon " << c5c6ReconExact
                << " selectedReconResidual "
                << c5ActualMinusC6 - c5c6ReconExact
                << " selectedReconRel "
                << mag(c5ActualMinusC6 - c5c6ReconExact)
                  /(mag(c5ActualMinusC6) + mag(c5c6ReconExact) + VSMALL)
                << endl;

            scalarField coeffResidualCell(mesh_.nCells(), Zero);
            scalar interiorCoeffResidual = Zero;
            List<scalar> patchCoeffResidual(mesh_.boundary().size(), Zero);
            boolList boundaryAdjacent(mesh_.nCells(), false);

            forAll(coeffResidualCell, celli)
            {
                scalar dTerm = barDExpected[celli]*dDExact[celli];
                if (dExactName == "A")
                {
                    dTerm = barAInternalExpected[celli]*dDExact[celli];
                }
                coeffResidualCell[celli] =
                    barR0Expected[celli]*dTape.rAtU[celli]
                  + trace.barq[celli]*dTape.q[celli]
                  - dTerm
                  - barH1Expected[celli]*dTape.UEqnH1Coeff[celli];
            }

            forAll(mesh_.boundary(), patchi)
            {
                const labelUList& faceCells =
                    mesh_.boundary()[patchi].faceCells();
                forAll(faceCells, facei)
                {
                    boundaryAdjacent[faceCells[facei]] = true;
                    patchCoeffResidual[patchi] +=
                        coeffResidualCell[faceCells[facei]];
                }
            }
            forAll(coeffResidualCell, celli)
            {
                if (!boundaryAdjacent[celli])
                {
                    interiorCoeffResidual += coeffResidualCell[celli];
                }
            }
            reduce(interiorCoeffResidual, sumOp<scalar>());
            forAll(patchCoeffResidual, patchi)
            {
                reduce(patchCoeffResidual[patchi], sumOp<scalar>());
                if (mag(patchCoeffResidual[patchi]) > SMALL)
                {
                    Info<< "ATC-T coefficient residual patch sum: "
                        << seedName << " " << directionName
                        << " eps " << eps
                        << " patch " << mesh_.boundary()[patchi].name()
                        << " selectedDExact " << dExactName
                        << " residualSum " << patchCoeffResidual[patchi]
                        << endl;
                }
            }
            Info<< "ATC-T coefficient residual interior sum: "
                << seedName << " " << directionName
                << " eps " << eps
                << " selectedDExact " << dExactName
                << " residualSum " << interiorCoeffResidual
                << endl;

            FixedList<label, 20> topCells(-1);
            FixedList<scalar, 20> topAbsResidual(Zero);
            forAll(coeffResidualCell, celli)
            {
                const scalar absResidual = mag(coeffResidualCell[celli]);
                forAll(topAbsResidual, sloti)
                {
                    if (absResidual > topAbsResidual[sloti])
                    {
                        for
                        (
                            label shifti = topAbsResidual.size() - 1;
                            shifti > sloti;
                            --shifti
                        )
                        {
                            topAbsResidual[shifti] =
                                topAbsResidual[shifti - 1];
                            topCells[shifti] = topCells[shifti - 1];
                        }
                        topAbsResidual[sloti] = absResidual;
                        topCells[sloti] = celli;
                        break;
                    }
                }
            }

            forAll(topCells, sloti)
            {
                const label celli = topCells[sloti];
                if (celli < 0 || topAbsResidual[sloti] <= SMALL)
                {
                    continue;
                }

                word adjacentPatchName("interior");
                forAll(mesh_.boundary(), patchi)
                {
                    const labelUList& faceCells =
                        mesh_.boundary()[patchi].faceCells();
                    forAll(faceCells, facei)
                    {
                        if (faceCells[facei] == celli)
                        {
                            adjacentPatchName = mesh_.boundary()[patchi].name();
                        }
                    }
                }

                Info<< "ATC-T coefficient largest cell residual: "
                    << seedName << " " << directionName
                    << " eps " << eps
                    << " rank " << sloti
                    << " cell " << celli
                    << " adjacentPatch " << adjacentPatchName
                    << " selectedDExact " << dExactName
                    << " Ddiag " << baseTapeForStages.Ddiag[celli]
                    << " Dscalar " << baseTapeForStages.Dscalar[celli]
                    << " AinternalTimesV "
                    << baseTapeForStages.AinternalTimesV[celli]
                    << " rAU " << baseTapeForStages.rAU[celli]
                    << " rAUUsed " << trace.rAUUsed[celli]
                    << " H1 " << baseTapeForStages.UEqnH1Coeff[celli]
                    << " rAtU " << baseTapeForStages.rAtU[celli]
                    << " rAtUUsed " << trace.rAtUUsed[celli]
                    << " q " << baseTapeForStages.q[celli]
                    << " dDexact " << dDExact[celli]
                    << " drAU " << dTape.rAU[celli]
                    << " dH1 " << dTape.UEqnH1Coeff[celli]
                    << " drAtU " << dTape.rAtU[celli]
                    << " dq " << dTape.q[celli]
                    << " barR0 " << barR0Expected[celli]
                    << " barQ " << trace.barq[celli]
                    << " barRExpected " << barRExpected[celli]
                    << " barAExpected " << barAExpected[celli]
                    << " barH1Expected " << barH1Expected[celli]
                    << " barDExpected " << barDExpected[celli]
                    << " residualCell " << coeffResidualCell[celli]
                    << endl;
            }

            const scalar isolatedGapSum =
                gapD + gapH1 + gapCons + gapPressure;
            const scalar reconstructionRel =
                mag(fullSignedGap - isolatedGapSum)
               /(mag(fullSignedGap) + mag(isolatedGapSum) + VSMALL);
            const scalar correctedGapSum =
                gapY
              + gapD
              + gapH1
              + gapAdjust
              + gapFlux
              + gapConstrain
              + gapCons
              + gapPressure;
            const scalar correctedReconstructionRel =
                mag(fullSignedGap - correctedGapSum)
               /(
                    mag(fullSignedGap)
                  + mag(gapY)
                  + mag(gapD)
                  + mag(gapH1)
                  + mag(gapAdjust)
                  + mag(gapFlux)
                  + mag(gapConstrain)
                  + mag(gapCons)
                  + mag(gapPressure)
                  + VSMALL
                );
            const scalar telescopingSum =
                (c0Output - c1AfterFinal)
              + (c1AfterFinal - c2AfterRelax)
              + (c2AfterRelax - c3AfterPressureFixed)
              + (c3AfterPressureFixed - c4AfterPressureCoeffFD)
              + (c4AfterPressureCoeffFD - c5AfterConsistentWithPressureFD)
              + (
                    c5AfterConsistentWithPressureFD
                  - c6AfterCoeffAlgebraWithPressureFD
                )
              + (
                    c6AfterCoeffAlgebraWithPressureFD
                  - c7AfterFluxWithPressureFD
                )
              + (c7AfterFluxWithPressureFD - c7AfterPredictorInputs)
              + (c7AfterPredictorInputs - c8AfterBoundaryEval)
              + (c8AfterBoundaryEval - c9AfterPredictorReverse)
              + (c9AfterPredictorReverse - c10OldState);
            const scalar telescopingTarget = c0Output - c10OldState;
            const scalar telescopingGap = telescopingSum - telescopingTarget;
            const scalar telescopingActualSum =
                (c0Output - c1AfterFinal)
              + (c1AfterFinal - c2AfterRelax)
              + (c2AfterRelax - c3AfterPressureFixed)
              + (c3AfterPressureFixed - c4ActualPressureReverse)
              + (c4ActualPressureReverse - c5AfterConsistentActual)
              + (c5AfterConsistentActual - c6AfterCoeffAlgebra)
              + (c6AfterCoeffAlgebra - c7AfterFlux)
              + (c7AfterFlux - c7AfterPredictorInputs)
              + (c7AfterPredictorInputs - c8AfterBoundaryEval)
              + (c8AfterBoundaryEval - c9AfterPredictorReverse)
              + (c9AfterPredictorReverse - c10OldState);
            const scalar telescopingActualGap =
                telescopingActualSum - telescopingTarget;
            const scalar telescopingL1 =
                mag(c0Output - c1AfterFinal)
              + mag(c1AfterFinal - c2AfterRelax)
              + mag(c2AfterRelax - c3AfterPressureFixed)
              + mag(c3AfterPressureFixed - c4AfterPressureCoeffFD)
              + mag(c4AfterPressureCoeffFD - c5AfterConsistentWithPressureFD)
              + mag
                (
                    c5AfterConsistentWithPressureFD
                  - c6AfterCoeffAlgebraWithPressureFD
                )
              + mag
                (
                    c6AfterCoeffAlgebraWithPressureFD
                  - c7AfterFluxWithPressureFD
                )
              + mag(c7AfterFluxWithPressureFD - c7AfterPredictorInputs)
              + mag(c7AfterPredictorInputs - c8AfterBoundaryEval)
              + mag(c8AfterBoundaryEval - c9AfterPredictorReverse)
              + mag(c9AfterPredictorReverse - c10OldState)
              + VSMALL;
            const scalar telescopingActualL1 =
                mag(c0Output - c1AfterFinal)
              + mag(c1AfterFinal - c2AfterRelax)
              + mag(c2AfterRelax - c3AfterPressureFixed)
              + mag(c3AfterPressureFixed - c4ActualPressureReverse)
              + mag(c4ActualPressureReverse - c5AfterConsistentActual)
              + mag(c5AfterConsistentActual - c6AfterCoeffAlgebra)
              + mag(c6AfterCoeffAlgebra - c7AfterFlux)
              + mag(c7AfterFlux - c7AfterPredictorInputs)
              + mag(c7AfterPredictorInputs - c8AfterBoundaryEval)
              + mag(c8AfterBoundaryEval - c9AfterPredictorReverse)
              + mag(c9AfterPredictorReverse - c10OldState)
              + VSMALL;

            printAdjacentStage("C0_output", "C1_afterFinalCorrection", c0Output, c1AfterFinal);
            printAdjacentStage("C1_afterFinalCorrection", "C2_afterPressureRelaxation", c1AfterFinal, c2AfterRelax);
            printAdjacentStage("C2_afterPressureRelaxation", "C3_afterPressureStageFixedCoefficients", c2AfterRelax, c3AfterPressureFixed);
            printAdjacentStage("C3_afterPressureStageFixedCoefficients", "C4_afterPressureStageIncludingRAtUCoefficientFD", c3AfterPressureFixed, c4AfterPressureCoeffFD);
            printAdjacentStage("C3_afterPressureStageFixedCoefficients", "C4_actualPressureReverse", c3AfterPressureFixed, c4ActualPressureReverse);
            printAdjacentStage("C4_afterPressureStageIncludingRAtUCoefficientFD", "C5_afterConsistentCorrectionFDVariant", c4AfterPressureCoeffFD, c5AfterConsistentWithPressureFD);
            printAdjacentStage("C4_actualPressureReverse", "C5_afterConsistentCorrectionActual", c4ActualPressureReverse, c5AfterConsistentActual);
            printAdjacentStage("C5_afterConsistentCorrectionCurrent", "C6_afterQ_RAtU_RAUAlgebraCurrent", c5AfterConsistent, c6AfterCoeffAlgebra);
            printAdjacentStage("C5_afterConsistentCorrectionActual", "C6_afterQ_RAtU_RAUAlgebraActual", c5AfterConsistentActual, c6AfterCoeffAlgebra);
            printAdjacentStage("C5_afterConsistentCorrectionFDVariant", "C6_afterQ_RAtU_RAUAlgebraFDVariant", c5AfterConsistentWithPressureFD, c6AfterCoeffAlgebraWithPressureFD);
            printAdjacentStage("C6_afterQ_RAtU_RAUAlgebraCurrent", "C7_afterFluxCurrent", c6AfterCoeffAlgebra, c7AfterFlux);
            printAdjacentStage("C6_afterQ_RAtU_RAUAlgebraFDVariant", "C7_afterFluxFDVariant", c6AfterCoeffAlgebraWithPressureFD, c7AfterFluxWithPressureFD);
            printAdjacentStage("C7_afterFluxFDVariant", "C8_afterConstrainHbyA", c7AfterFluxWithPressureFD, c7AfterPredictorInputs);
            printAdjacentStage("C8_afterConstrainHbyA", "C8b_afterVelocityBoundaryEvaluation", c7AfterPredictorInputs, c8AfterBoundaryEval);
            printAdjacentStage("C8b_afterVelocityBoundaryEvaluation", "C9_afterPredictorReverse", c8AfterBoundaryEval, c9AfterPredictorReverse);
            printAdjacentStage("C9_afterPredictorReverse", "C10_oldStateContraction", c9AfterPredictorReverse, c10OldState);

            Info<< "ATC-T smooth-phi coefficient trace: "
                << seedName << " " << directionName << " eps " << eps
                << " traceSeedMaxDiff " << traceSeedDiff
                << " nMomentumPhiSignCrossings "
                << countMomentumSignCrossings(dir, eps)
                << " dDscalarL1 " << dDNorm
                << " dDscalarMax " << dDMax
                << " dDdiagL1 " << dDdiagNorm
                << " dDdiagMax " << dDdiagMax
                << " dDboundaryL1 " << dDboundaryNorm
                << " dDboundaryMax " << dDboundaryMax
                << " dH1L1 " << dH1Norm
                << " dH1Max " << dH1Max
                << " dYInternalL1 " << dYInternalNorm
                << " dYBoundaryL1 " << dYBoundaryNorm
                << " barYBoundaryL1 " << barYBoundaryNorm
                << " boundaryDfinalActualTotal "
                << dFinalBoundaryActualTotal
                << " boundaryDfinalCurrentTotal "
                << dFinalBoundaryCurrentTotal
                << " boundaryDfinalCandidateTotal "
                << dFinalBoundaryCandidateTotal
                << " boundaryDfinalCandidateGapL1 "
                << dFinalBoundaryCandidateGapL1
                << " drAUL1 " << drAUNorm
                << " drAUMax " << drAUMax
                << " drAtUL1 " << drAtUNorm
                << " drAtUMax " << drAtUMax
                << " dqL1 " << dqNorm
                << " dqMax " << dqMax
                << " lhsDscalar " << lhsD
                << " lhsDdiag " << lhsDdiag
                << " lhsDboundary " << lhsDboundary
                << " rhsDInternal " << dRhsPhiI
                << " rhsDBoundary " << dRhsPhiB
                << " rhsDTotal " << rhsD
                << " gapD " << gapD
                << " gapDOverL1 " << mag(gapD)/scaleD
                << " rhsDiagCandidate " << rhsDiagCandidate
                << " gapDiagCandidateToDiag "
                << lhsDdiag - rhsDiagCandidate
                << " gapDiagCandidateToDscalar "
                << lhsD - rhsDiagCandidate
                << " lhsH1 " << lhsH1
                << " rhsH1Internal " << h1RhsPhiI
                << " rhsH1Boundary " << h1RhsPhiB
                << " rhsH1Total " << rhsH1
                << " gapH1 " << gapH1
                << " gapH1OverL1 " << mag(gapH1)/scaleH1
                << " lhsYInternal " << lhsYInternal
                << " lhsYBoundary " << lhsYBoundary
                << " rhsYU " << yRhsU
                << " rhsYInternalPhi " << yRhsPhiI
                << " rhsYBoundaryPhi " << yRhsPhiB
                << " rhsYTotal " << rhsY
                << " gapY " << gapY
                << " gapYOverL1 " << mag(gapY)/scaleY
                << " lhsPredictor " << lhsPredictor
                << " rhsPredictorU " << combinedRhsU
                << " rhsPredictorInternalPhi " << combinedRhsPhiI
                << " rhsPredictorBoundaryPhi " << combinedRhsPhiB
                << " rhsPredictorTotal " << rhsPredictor
                << " gapPredictorCombined " << gapPredictorCombined
                << " predictorLinearityGap " << predictorLinearityGap
                << " predictorLinearityRel " << predictorLinearityRel
                << " lhsAdjust " << lhsAdjust
                << " rhsAdjust " << rhsAdjust
                << " gapAdjust " << gapAdjust
                << " gapAdjustOverL1 " << traceGapOverL1(lhsAdjust, rhsAdjust)
                << " directH0Contraction " << barH0DirectContraction
                << " h0FromFluxContraction " << barH0FromFluxContraction
                << " h0TotalContraction " << barH0TotalContraction
                << " h0TotalSplitGap " << barH0SplitGap
                << " lhsFlux " << lhsFlux
                << " rhsFlux " << rhsFlux
                << " gapFlux " << gapFlux
                << " gapFluxOverL1 " << traceGapOverL1(lhsFlux, rhsFlux)
                << " lhsConstrain " << lhsConstrain
                << " rhsConstrain " << rhsConstrain
                << " gapConstrain " << gapConstrain
                << " gapConstrainOverL1 "
                << traceGapOverL1(lhsConstrain, rhsConstrain)
                << " lhsQH1 " << lhsQH1
                << " lhsQInternal " << lhsQInternal
                << " lhsQBoundary " << lhsQBoundary
                << " lhsQTotal " << lhsQ
                << " rhsQCurrent " << rhsQCurrent
                << " rhsQZeroBoundary " << rhsQZeroBoundary
                << " rhsQBoundaryCurrent " << rhsQBoundaryCurrent
                << " gapQBoundaryCurrent " << gapQBoundaryCurrent
                << " gapQBoundaryCurrentRel " << qCurrentRel
                << " gapQBoundaryZero " << gapQBoundaryZero
                << " gapQBoundaryZeroRel " << qZeroRel
                << " qSplitMaxDiff " << qSplitMaxDiff
                << " qBoundarySeedL1 " << qBoundarySeedL1
                << " qBoundaryForwardL1 " << qBoundaryForwardL1
                << " gapQBoundaryVsFullRel "
                << mag(fullSignedGap - gapQBoundaryCurrent)
                  /(mag(fullSignedGap) + mag(gapQBoundaryCurrent) + VSMALL)
                << " remainingAfterQ "
                << fullSignedGap - gapQBoundaryCurrent
                << " lhsConsistent " << lhsCons
                << " rhsConsistent " << rhsConsTotal
                << " gapConsistent " << gapCons
                << " gapConsistentRel "
                << mag(gapCons)/max(scaleCons, VSMALL)
                << " pressurePSolveCoeff " << pSolveCoeff
                << " pressurePhiInternalCoeff " << phiCoeffI
                << " pressurePhiBoundaryCoeff " << phiCoeffB
                << " lhsPressureCoeff " << lhsPressureCoeff
                << " rhsPressureCoeff " << rhsPressureCoeff
                << " gapPressureCoeff " << gapPressure
                << " gapPressureCoeffRel "
                << mag(gapPressure)/max(scalePressure, VSMALL)
                << " fixedPFluxInternalCoeff " << fixedPFluxInternal
                << " fixedPFluxBoundaryCoeff " << fixedPFluxBoundary
                << " fixedPFluxTotalCoeff " << fixedPFluxTotal
                << " pressureCoefficientOnly "
                << pressureCoefficientOnly
                << " pressureBoundaryStateEffect "
                << pressureBoundaryStateEffect
                << " pressureCoefficientMatrixFlux "
                << pressureCoefficientMatrixFlux
                << " pressureFaceFluxCorrection "
                << pressureFaceFluxCorrection
                << " pressureFluxDecompositionGap "
                << pressureFluxDecompositionGap
                << " fixedPIdentityGap " << fixedPIdentityGap
                << " fixedPReverseGap " << fixedPReverseGap
                << " pressureFaceSpaceLhsInternal " << faceLhsInternal
                << " pressureFaceSpaceLhsBoundary " << faceLhsBoundary
                << " pressureFaceSpaceRhsInternal " << faceRhsInternal
                << " pressureFaceSpaceRhsBoundary " << faceRhsBoundary
                << " pressureFaceSpaceGap "
                << faceSpaceLhs - faceSpaceRhs
                << " pressureFaceSpaceRel "
                << mag(faceSpaceLhs - faceSpaceRhs)/faceSpaceScale
                << " pressureInterpOrientARel "
                << internalOrientADiff/internalOrientScale
                << " pressureInterpOrientBRel "
                << internalOrientBDiff/internalOrientScale
                << " pressureInterpUseOrientationA "
                << useOrientationA
                << " pressureInternalFaceSeedContribution "
                << internalCellContribution
                << " pressureBoundaryFaceSeedContribution "
                << boundaryCellContribution
                << " pressureCellSeedTotal " << totalCellContribution
                << " pressureCellSeedGap "
                << faceSpaceLhs - totalCellContribution
                << " pressureCellSeedRel "
                << mag(faceSpaceLhs - totalCellContribution)/cellSpaceScale
                << " pressureDirectOwnerVariant "
                << directOwnerCellContribution
                << " pressureOmitBoundaryVariant "
                << omittedBoundaryContribution
                << " pressureBoundaryOwnerUnresolvedL1 "
                << boundaryOwnerMapUnresolved
                << " finalRAtUContribution " << finalRAtUContribution
                << " fullSignedGap " << fullSignedGap
                << " isolatedGapSum " << isolatedGapSum
                << " reconstructionRel " << reconstructionRel
                << " correctedReconstructedGap " << correctedGapSum
                << " correctedReconstructionRel "
                << correctedReconstructionRel
                << " telescopingSum " << telescopingSum
                << " telescopingTarget " << telescopingTarget
                << " telescopingGap " << telescopingGap
                << " telescopingGapOverL1 "
                << mag(telescopingGap)/telescopingL1
                << " telescopingActualSum " << telescopingActualSum
                << " telescopingActualGap " << telescopingActualGap
                << " telescopingActualGapOverL1 "
                << mag(telescopingActualGap)/telescopingActualL1
                << " C0 " << c0Output
                << " C1 " << c1AfterFinal
                << " C2 " << c2AfterRelax
                << " C3 " << c3AfterPressureFixed
                << " C4ActualPressureReverse "
                << c4ActualPressureReverse
                << " C4PressureCoeffFD " << c4AfterPressureCoeffFD
                << " C5Current " << c5AfterConsistent
                << " C5ActualPressureReverse "
                << c5AfterConsistentActual
                << " C5PressureCoeffFD " << c5AfterConsistentWithPressureFD
                << " C6Current " << c6AfterCoeffAlgebra
                << " C6PressureCoeffFD "
                << c6AfterCoeffAlgebraWithPressureFD
                << " C7Current " << c7AfterFlux
                << " C7PressureCoeffFD " << c7AfterFluxWithPressureFD
                << " C8 " << c7AfterPredictorInputs
                << " C8b " << c8AfterBoundaryEval
                << " C9 " << c9AfterPredictorReverse
                << " C10 " << c10OldState
                << endl;
        }
    };

    Info<< "ATC-T full-state map transpose columns: seed direction eps "
        << "lhsU lhsP lhsPhiInternal lhsPhiBoundary "
        << "rhsU rhsP rhsPhiInternal rhsPhiBoundary lhs rhs rel "
        << "absoluteGap symmetricRelativeGap contractionL1Scale "
        << "gapOverL1Scale"
        << endl;
    Info<< "ATC-T full-state base defect note: values are checked through "
        << "the existing tangent diagnostic; this transpose check uses the "
        << "stored state as the linearisation point."
        << endl;

    forAll(epsList, epsi)
    {
        const scalar eps = epsList[epsi];

        runSeedDirection("randomSeed", randomSeed, "combined", direction, eps);
        runSeedDirection("randomSeed", randomSeed, "UOnly", directionU, eps);
        runSeedDirection("randomSeed", randomSeed, "POnly", directionP, eps);
        runSeedDirection
        (
            "randomSeed",
            randomSeed,
            "PhiInternalOnly",
            directionPhiI,
            eps
        );
        runSeedDirection
        (
            "randomSeed",
            randomSeed,
            "PhiBoundaryOnly",
            directionPhiB,
            eps
        );
        runSeedDirection
        (
            "randomSeed",
            randomSeed,
            "PhiInternalSmoothOnly",
            directionPhiISmooth,
            eps
        );
        runSeedDirection
        (
            "randomSeed",
            randomSeed,
            "PhiBoundarySmoothOnly",
            directionPhiBSmooth,
            eps
        );
        runSeedDirection
        (
            "randomSeed",
            randomSeed,
            "combinedSmooth",
            directionSmooth,
            eps
        );

        runSeedDirection("thermalSeed", thermalSeed, "combined", direction, eps);
        runSeedDirection("thermalSeed", thermalSeed, "UOnly", directionU, eps);
        runSeedDirection("thermalSeed", thermalSeed, "POnly", directionP, eps);
        runSeedDirection
        (
            "thermalSeed",
            thermalSeed,
            "PhiInternalOnly",
            directionPhiI,
            eps
        );
        runSeedDirection
        (
            "thermalSeed",
            thermalSeed,
            "PhiBoundaryOnly",
            directionPhiB,
            eps
        );
        runSeedDirection
        (
            "thermalSeed",
            thermalSeed,
            "PhiInternalSmoothOnly",
            directionPhiISmooth,
            eps
        );
        runSeedDirection
        (
            "thermalSeed",
            thermalSeed,
            "PhiBoundarySmoothOnly",
            directionPhiBSmooth,
            eps
        );
        runSeedDirection
        (
            "thermalSeed",
            thermalSeed,
            "combinedSmooth",
            directionSmooth,
            eps
        );

        runSeedDirection
        (
            "seedAfterOneReverse",
            seedAfterOne,
            "combined",
            direction,
            eps
        );
        runSeedDirection
        (
            "seedAfterFourReverse",
            seedAfterFour,
            "combined",
            direction,
            eps
        );
    }

    restoreState();
}


void Foam::thermalAdjointSimple::checkPressureMapTranspose
(
    const vectorField& barUnewInt,
    const volScalarField& rAtU
)
{
    if (!checkPressureMapTranspose_)
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const scalarField& V = mesh_.V();
    volVectorField& UBase = primalVars_.U();
    const volScalarField& pBase = primalVars_.p();
    const surfaceScalarField& phi = primalVars_.phi();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        primalVars_.turbulence();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    auto removeMean = [](scalarField& fld)
    {
        const label n =
            returnReduce(fld.size(), sumOp<label>());

        if (n > 0)
        {
            fld -= gSum(fld)/scalar(n);
        }
    };

    const bool pressureNeedsReference = pBase.needReference();
    const label pRefCell = solverControl_().pRefCell();

    auto projectPressureDirection = [&](scalarField& fld)
    {
        if (!pressureNeedsReference)
        {
            return;
        }

        if (pRefCell >= 0 && pRefCell < mesh_.nCells())
        {
            fld[pRefCell] = Zero;
        }
        else
        {
            removeMean(fld);
        }
    };

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, UBase)
      + turbulence->divDevReff(UBase)
     ==
        fvOptions(UBase)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
    fvOptions.constrain(UEqn);

    volScalarField rAU(1.0/UEqn.A());

    volScalarField pWork
    (
        IOobject
        (
            "ATCTPressureMapBasep",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pBase
    );

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), UBase, pWork));
    surfaceScalarField F2("ATCTPressureMapF2Base", fvc::flux(HbyA));
    adjustPhi(F2, UBase, pWork);

    if (solverControl_().consistent())
    {
        tmp<volScalarField> trAtU =
            1.0/(1.0/rAU - UEqn.H1());
        F2 +=
            fvc::interpolate(trAtU() - rAU)
           *fvc::snGrad(pWork)
           *mesh_.magSf();
    }

    surfaceScalarField dF2
    (
        IOobject
        (
            "ATCTPressureMapdF2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F2
    );
    dF2 = dimensionedScalar(dF2.dimensions(), Zero);

    forAll(dF2.primitiveFieldRef(), facei)
    {
        dF2.primitiveFieldRef()[facei] =
            scalar((307*facei + 233) % 383)/scalar(383) - scalar(0.5);
    }

    forAll(dF2.boundaryFieldRef(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            dF2.boundaryFieldRef()[patchi] == Zero;
            continue;
        }

        fvsPatchScalarField& dF2p = dF2.boundaryFieldRef()[patchi];
        forAll(dF2p, facei)
        {
            dF2p[facei] =
                scalar((311*(facei + 1) + 239*(patchi + 1)) % 389)
               /scalar(389) - scalar(0.5);
        }
    }

    volVectorField dUold
    (
        IOobject
        (
            "ATCTPressureMapdUold",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        UBase
    );
    dUold = dimensionedVector(dUold.dimensions(), vector::zero);

    forAll(dUold.primitiveFieldRef(), celli)
    {
        dUold.primitiveFieldRef()[celli] =
            vector
            (
                scalar((313*celli + 241) % 397)/scalar(397) - scalar(0.5),
                scalar((317*celli + 251) % 401)/scalar(401) - scalar(0.5),
                scalar((331*celli + 257) % 409)/scalar(409) - scalar(0.5)
            );
    }
    forAll(dUold.boundaryFieldRef(), patchi)
    {
        if (!mesh_.boundary()[patchi].coupled())
        {
            dUold.boundaryFieldRef()[patchi] == vector::zero;
        }
    }

    volScalarField dpold
    (
        IOobject
        (
            "ATCTPressureMapdpold",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pBase
    );
    dpold = dimensionedScalar(dpold.dimensions(), Zero);

    forAll(dpold.primitiveFieldRef(), celli)
    {
        dpold.primitiveFieldRef()[celli] =
            scalar((337*celli + 263) % 419)/scalar(419) - scalar(0.5);
    }
    projectPressureDirection(dpold.primitiveFieldRef());
    forAll(dpold.boundaryFieldRef(), patchi)
    {
        if (!mesh_.boundary()[patchi].coupled())
        {
            dpold.boundaryFieldRef()[patchi] == Zero;
        }
    }

    auto reversePressureMap =
        [&]
        (
            const scalarField& rawBarpNew,
            surfaceScalarField& barF2,
            vectorField& barUold,
            scalarField& barpOld
        )
    {
        scalarField barpNew(rawBarpNew);
        projectPressureDirection(barpNew);

        volScalarField mu
        (
            IOobject
            (
                "paATCTPressureMapSeed",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pBase
        );
        mu = dimensionedScalar(mu.dimensions(), Zero);

        fvScalarMatrix muEqn(fvm::laplacian(rAtU, mu));
        muEqn.source() = barpNew;
        muEqn.setReference(solverControl_().pRefCell(), Zero);
        dictionary muSolver(muEqn.solverDict("p"));
        muSolver.set("relTol", scalar(0));
        muSolver.set("tolerance", scalar(1e-12));
        muEqn.solve(muSolver);

        barF2 = dimensionedScalar(barF2.dimensions(), Zero);
        barUold = vector::zero;
        barpOld = Zero;

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];
            barF2.primitiveFieldRef()[facei] = mu[P] - mu[N];
        }

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            fvsPatchScalarField& barF2p =
                barF2.boundaryFieldRef()[patchi];
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(barF2p, facei)
            {
                barF2p[facei] = mu[faceCells[facei]];
            }
        }
    };

    auto pressureDot = [&](const scalarField& a, const scalarField& b)
    {
        scalar value = Zero;
        forAll(a, celli)
        {
            value += a[celli]*b[celli];
        }
        return value;
    };

    auto faceDot =
        [&](const surfaceScalarField& a, const surfaceScalarField& b)
    {
        scalar value = Zero;
        forAll(a.primitiveField(), facei)
        {
            value += a[facei]*b[facei];
        }

        forAll(mesh_.boundary(), patchi)
        {
            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            const fvsPatchScalarField& ap = a.boundaryField()[patchi];
            const fvsPatchScalarField& bp = b.boundaryField()[patchi];
            forAll(ap, facei)
            {
                value += ap[facei]*bp[facei];
            }
        }

        return value;
    };

    auto vectorDot = [&](const vectorField& a, const volVectorField& b)
    {
        scalar value = Zero;
        forAll(a, celli)
        {
            value += a[celli] & b[celli];
        }
        return value;
    };

    auto runCase =
        [&]
        (
            const word& seedName,
            const scalarField& rawBarpNew,
            const word& caseName,
            const bool useInternalF,
            const bool useBoundaryF,
            const bool useU,
            const bool useP
        )
    {
        scalarField barpNew(rawBarpNew);
        projectPressureDirection(barpNew);

        surfaceScalarField dFCase
        (
            IOobject
            (
                "ATCTPressureMapdFCase",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dF2
        );
        dFCase = dimensionedScalar(dFCase.dimensions(), Zero);

        if (useInternalF)
        {
            dFCase.primitiveFieldRef() = dF2.primitiveField();
        }

        if (useBoundaryF)
        {
            forAll(dFCase.boundaryFieldRef(), patchi)
            {
                dFCase.boundaryFieldRef()[patchi] ==
                    dF2.boundaryField()[patchi];
            }
        }

        volVectorField dUCase
        (
            IOobject
            (
                "ATCTPressureMapdUCase",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dUold
        );
        dUCase = dimensionedVector(dUCase.dimensions(), vector::zero);
        if (useU)
        {
            dUCase.primitiveFieldRef() = dUold.primitiveField();
        }

        volScalarField dpCase
        (
            IOobject
            (
                "ATCTPressureMapdpCase",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dpold
        );
        dpCase = dimensionedScalar(dpCase.dimensions(), Zero);
        if (useP)
        {
            dpCase.primitiveFieldRef() = dpold.primitiveField();
            projectPressureDirection(dpCase.primitiveFieldRef());
        }

        surfaceScalarField FPlus("ATCTPressureMapFPlus", F2);
        surfaceScalarField FMinus("ATCTPressureMapFMinus", F2);
        FPlus.primitiveFieldRef() =
            F2.primitiveField() + scalar(1e-6)*dFCase.primitiveField();
        FMinus.primitiveFieldRef() =
            F2.primitiveField() - scalar(1e-6)*dFCase.primitiveField();

        forAll(FPlus.boundaryFieldRef(), patchi)
        {
            FPlus.boundaryFieldRef()[patchi] ==
                F2.boundaryField()[patchi]
              + scalar(1e-6)*dFCase.boundaryField()[patchi];
            FMinus.boundaryFieldRef()[patchi] ==
                F2.boundaryField()[patchi]
              - scalar(1e-6)*dFCase.boundaryField()[patchi];
        }

        volVectorField UPlus("ATCTPressureMapUPlus", UBase);
        volVectorField UMinus("ATCTPressureMapUMinus", UBase);
        UPlus.primitiveFieldRef() =
            UBase.primitiveField() + scalar(1e-6)*dUCase.primitiveField();
        UMinus.primitiveFieldRef() =
            UBase.primitiveField() - scalar(1e-6)*dUCase.primitiveField();
        forAll(UPlus.boundaryFieldRef(), patchi)
        {
            UPlus.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
            UMinus.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
        }

        volScalarField pPlusOld("ATCTPressureMappPlusOld", pBase);
        volScalarField pMinusOld("ATCTPressureMappMinusOld", pBase);
        pPlusOld.primitiveFieldRef() =
            pBase.primitiveField() + scalar(1e-6)*dpCase.primitiveField();
        pMinusOld.primitiveFieldRef() =
            pBase.primitiveField() - scalar(1e-6)*dpCase.primitiveField();
        forAll(pPlusOld.boundaryFieldRef(), patchi)
        {
            pPlusOld.boundaryFieldRef()[patchi] ==
                pBase.boundaryField()[patchi];
            pMinusOld.boundaryFieldRef()[patchi] ==
                pBase.boundaryField()[patchi];
        }

        tmp<volScalarField> tpPlus =
            computePressureMap(FPlus, UPlus, pPlusOld, rAtU, "ATCPMapPlus");
        tmp<volScalarField> tpMinus =
            computePressureMap(FMinus, UMinus, pMinusOld, rAtU, "ATCPMapMinus");

        scalarField dpnew
        (
            (tpPlus().primitiveField() - tpMinus().primitiveField())
           /(2e-6)
        );
        projectPressureDirection(dpnew);

        surfaceScalarField barF2("ATCTPressureMapbarF2", F2);
        vectorField barUold(mesh_.nCells(), vector::zero);
        scalarField barpOld(mesh_.nCells(), Zero);
        reversePressureMap(barpNew, barF2, barUold, barpOld);

        const scalar lhs = pressureDot(barpNew, dpnew);
        const scalar rhsF = faceDot(barF2, dFCase);
        const scalar rhsU = vectorDot(barUold, dUCase);
        const scalar rhsP = pressureDot(barpOld, dpCase.primitiveField());
        const scalar rhs = rhsF + rhsU + rhsP;
        const scalar scale = max(max(mag(lhs), mag(rhs)), VSMALL);

        Info<< "ATC-T pressure-map transpose check: "
            << seedName << "/" << caseName
            << " lhs = " << lhs
            << ", rhsF = " << rhsF
            << ", rhsU = " << rhsU
            << ", rhsP = " << rhsP
            << ", rhs = " << rhs
            << ", rel = " << mag(lhs - rhs)/scale
            << endl;
    };

    scalarField randomBarp(mesh_.nCells(), Zero);
    forAll(randomBarp, celli)
    {
        randomBarp[celli] =
            scalar((347*celli + 269) % 421)/scalar(421) - scalar(0.5);
    }
    projectPressureDirection(randomBarp);

    scalarField actualBarp(mesh_.nCells(), Zero);
    vectorField source(mesh_.nCells(), vector::zero);
    forAll(source, celli)
    {
        source[celli] = barUnewInt[celli]/V[celli];
    }

    const surfaceScalarField& w = mesh_.weights();
    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar SfSource =
            Sf[facei] & (rAtU[P]*source[P] - rAtU[N]*source[N]);

        actualBarp[P] -= w[facei]*SfSource;
        actualBarp[N] -= (scalar(1) - w[facei])*SfSource;
    }

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            mesh_.boundary()[patchi].type() == "empty"
         || mesh_.boundary()[patchi].coupled()
        )
        {
            continue;
        }

        const fvPatch& patch = mesh_.boundary()[patchi];
        const labelUList& faceCells = patch.faceCells();
        const fvsPatchVectorField& SfPatch = Sf.boundaryField()[patchi];

        tmp<scalarField> tValueInternalCoeffs =
            pBase.boundaryField()[patchi].valueInternalCoeffs
            (
                patch.weights()
            );
        const scalarField& valueInternalCoeffs = tValueInternalCoeffs();

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            actualBarp[celli] -=
                valueInternalCoeffs[facei]
               *(SfPatch[facei] & (rAtU[celli]*source[celli]));
        }
    }
    projectPressureDirection(actualBarp);

    Info<< "ATC-T pressure-map transpose perturbation: eps = 1e-06"
        << ", needReference = " << pressureNeedsReference
        << ", pRefCell = " << pRefCell
        << ", pressureGauge = "
        <<
        (
            !pressureNeedsReference
          ? "fixed"
          : (pRefCell >= 0 && pRefCell < mesh_.nCells() ? "reference" : "mean")
        )
        << ". Analytic reverse currently includes pressure solve D^T only."
        << endl;

    const wordList seedNames({"randomBarp", "actualFinalCorrectionBarp"});
    const List<const scalarField*> seeds({&randomBarp, &actualBarp});

    forAll(seedNames, seedi)
    {
        runCase(seedNames[seedi], *seeds[seedi], "internalF", true, false, false, false);
        runCase(seedNames[seedi], *seeds[seedi], "boundaryF", false, true, false, false);
        runCase(seedNames[seedi], *seeds[seedi], "Uold", false, false, true, false);
        runCase(seedNames[seedi], *seeds[seedi], "pold", false, false, false, true);
    }
}


void Foam::thermalAdjointSimple::checkFixedPointMapAdjointSensitivity
(
    const vectorField& barUnewInt,
    const volScalarField& rAtU,
    const word& designVariablesName,
    const scalar dt
)
{
    if (!checkFixedPointMapAdjoint_)
    {
        return;
    }

    if (!mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        Info<< "ATC-T fixed-point map adjoint check: no topOVars object"
            << endl;
        return;
    }

    const topOVariablesBase& vars =
        mesh_.lookupObject<topOVariablesBase>("topOVars");
    const topOZones& zones = vars.getTopOZones();
    const scalarField& V = mesh_.V();
    volScalarField& beta = const_cast<volScalarField&>(vars.beta());

    bitSet isDesign(mesh_.nCells(), false);
    if (zones.adjointPorousZoneIDs().empty())
    {
        isDesign.fill(true);
    }
    else
    {
        for (const label zoneID : zones.adjointPorousZoneIDs())
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.set(celli);
            }
        }
    }

    auto unsetZone = [&](const label zoneID)
    {
        if (zoneID != -1)
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.unset(celli);
            }
        }
    };

    for (const label zoneID : zones.fixedPorousZoneIDs())
    {
        unsetZone(zoneID);
    }
    for (const label zoneID : zones.fixedZeroPorousZoneIDs())
    {
        unsetZone(zoneID);
    }
    for (const label celli : zones.IOCells())
    {
        isDesign.unset(celli);
    }

    labelList reportCells(fixedPointTangentCells_.size(), -1);
    const labelList& requestedCells = fixedPointTangentCells_;

    label nReport = 0;
    forAll(requestedCells, i)
    {
        const label celli = requestedCells[i];
        if (celli >= 0 && celli < mesh_.nCells() && isDesign.test(celli))
        {
            reportCells[nReport++] = celli;
        }
    }

    if (nReport == 0)
    {
        Info<< "ATC-T fixed-point map adjoint check: no requested cells in "
            << "the active design set. Diagnostic only; betaMult is unchanged."
            << endl;
        return;
    }

    scalarField conductivitySens(mesh_.nCells(), Zero);
    {
        const volScalarField& T = TRef();
        const volScalarField& Ta = TaPtr_();
        const autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();
        const volScalarField Df(props_.DFluidField(T, "DFluidFieldFPMapCheck"));
        const volScalarField Ds(props_.DSolidField(T, "DSolidFieldFPMapCheck"));

        scalarField dJdDEff(mesh_.nCells(), Zero);
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const surfaceScalarField& w = mesh_.weights();
        const surfaceScalarField& deltaCoeffs = mesh_.nonOrthDeltaCoeffs();
        const surfaceScalarField& magSf = mesh_.magSf();

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];

            const scalar dJdDf =
                magSf[facei]*deltaCoeffs[facei]
               *(T[N] - T[P])*(Ta[N] - Ta[P]);

            dJdDEff[P] += w[facei]*dJdDf;
            dJdDEff[N] += (scalar(1) - w[facei])*dJdDf;
        }

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
                dJdDEff[fc[i]] -= Ta[fc[i]]*magSfb[i]*snGradT[i];
            }
        }

        conductivitySens =
            thermalSensScale_
           *dJdDEff/V
           *(
                Ds.primitiveField() - Df.primitiveField()
              - turbulence->nut()().primitiveField()/Prt_
            )
           *dt;

        vars.sourceTermSensitivities
        (
            conductivitySens,
            kInterpolation_(),
            scalar(1),
            designVariablesName,
            "beta"
        );
    }

    SimpleMapSeed seed(mesh_);
    seed.barU = barUnewInt;
    scalarField accumMapVol(mesh_.nCells(), Zero);

    scalar initialSeedNorm = VSMALL;
    forAll(seed.barU, celli)
    {
        initialSeedNorm += magSqr(seed.barU[celli]) + sqr(seed.barp[celli]);
    }
    initialSeedNorm = sqrt(initialSeedNorm);

    Info<< "ATC-T fixed-point map adjoint samples: "
        << "iter cell beta eps betaPartUVol/dt betaPartPVol/dt "
        << "betaPartTotalVol/dt accumMapVol/dt conductivityVol/dt "
        << "accumPlusCond/dt seedUNormRel seedPNormRel"
        << endl;

    for (label iter = 0; iter < fixedPointMapAdjointIters_; ++iter)
    {
        scalar seedUNorm = VSMALL;
        scalar seedPNorm = VSMALL;
        forAll(seed.barU, celli)
        {
            seedUNorm += magSqr(seed.barU[celli]);
            seedPNorm += sqr(seed.barp[celli]);
        }
        seedUNorm = sqrt(seedUNorm);
        seedPNorm = sqrt(seedPNorm);
        const scalar seedUNormRel = seedUNorm/initialSeedNorm;
        const scalar seedPNormRel = seedPNorm/initialSeedNorm;

        for (label reporti = 0; reporti < nReport; ++reporti)
        {
            const label cellj = reportCells[reporti];
            const scalar oldBeta = beta[cellj];
            const scalar room = min(oldBeta, scalar(1) - oldBeta);
            const scalar eps =
                min(fixedPointMapAdjointFDEps_, scalar(0.25)*max(SMALL, room));

            if (eps <= SMALL)
            {
                continue;
            }

            beta.primitiveFieldRef()[cellj] = oldBeta + eps;
            beta.correctBoundaryConditions();
            SimpleMapState statePlus =
                primalSimpleMapStateAtFrozenState("ATCFPMapPlus");

            beta.primitiveFieldRef()[cellj] = oldBeta - eps;
            beta.correctBoundaryConditions();
            SimpleMapState stateMinus =
                primalSimpleMapStateAtFrozenState("ATCFPMapMinus");

            beta.primitiveFieldRef()[cellj] = oldBeta;
            beta.correctBoundaryConditions();

            scalar betaPartU = Zero;
            scalar betaPartP = Zero;
            forAll(V, celli)
            {
                const vector dUnew =
                    (statePlus.U[celli] - stateMinus.U[celli])/(2*eps);
                const scalar dpnew =
                    (statePlus.p[celli] - stateMinus.p[celli])/(2*eps);
                betaPartU += seed.barU[celli] & dUnew;
                betaPartP += seed.barp[celli]*dpnew;
            }

            const scalar betaPart = betaPartU + betaPartP;
            const scalar betaPartUVol = betaPartU/V[cellj];
            const scalar betaPartPVol = betaPartP/V[cellj];
            const scalar betaPartVol = betaPart/V[cellj];
            accumMapVol[cellj] += betaPartVol;

            Info<< "ATC-T fixed-point map adjoint sample: "
                << iter << " "
                << cellj << " "
                << oldBeta << " "
                << eps << " "
                << betaPartUVol/dt << " "
                << betaPartPVol/dt << " "
                << betaPartVol/dt << " "
                << accumMapVol[cellj]/dt << " "
                << conductivitySens[cellj]/dt << " "
                << (accumMapVol[cellj] + conductivitySens[cellj])/dt << " "
                << seedUNormRel << " "
                << seedPNormRel
                << endl;
        }

        seed = reverseOneSimpleMapSeed(seed, rAtU);
    }

    beta.correctBoundaryConditions();

    Info<< "ATC-T fixed-point map adjoint check: iters = "
        << fixedPointMapAdjointIters_
        << ", eps = " << fixedPointMapAdjointFDEps_
        << ". Diagnostic only; betaMult is unchanged."
        << endl;
}


void Foam::thermalAdjointSimple::checkFixedPointTangentSensitivity
(
    const vectorField& barUnewInt,
    const word& designVariablesName,
    const scalar dt
)
{
    if (!checkFixedPointTangentSensitivity_)
    {
        return;
    }

    if (!mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        Info<< "ATC-T fixed-point tangent check: no topOVars object"
            << endl;
        return;
    }

    const topOVariablesBase& vars =
        mesh_.lookupObject<topOVariablesBase>("topOVars");
    const topOZones& zones = vars.getTopOZones();
    const scalarField& V = mesh_.V();
    volScalarField& beta = const_cast<volScalarField&>(vars.beta());
    volVectorField& U = primalVars_.U();
    volScalarField& p =
        const_cast<volScalarField&>(primalVars_.p());
    surfaceScalarField& phi = primalVars_.phi();

    static_cast<void>(barUnewInt);
    tmp<surfaceScalarField> tGPhi = thermalFluxSensitivity();
    const surfaceScalarField& gPhi = tGPhi();

    bitSet isDesign(mesh_.nCells(), false);
    if (zones.adjointPorousZoneIDs().empty())
    {
        isDesign.fill(true);
    }
    else
    {
        for (const label zoneID : zones.adjointPorousZoneIDs())
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.set(celli);
            }
        }
    }

    auto unsetZone = [&](const label zoneID)
    {
        if (zoneID != -1)
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.unset(celli);
            }
        }
    };

    for (const label zoneID : zones.fixedPorousZoneIDs())
    {
        unsetZone(zoneID);
    }
    for (const label zoneID : zones.fixedZeroPorousZoneIDs())
    {
        unsetZone(zoneID);
    }
    for (const label celli : zones.IOCells())
    {
        isDesign.unset(celli);
    }

    labelList reportCells(fixedPointTangentCells_.size(), -1);
    const labelList& requestedCells = fixedPointTangentCells_;

    label nReport = 0;
    forAll(requestedCells, i)
    {
        const label celli = requestedCells[i];
        if (celli >= 0 && celli < mesh_.nCells() && isDesign.test(celli))
        {
            reportCells[nReport++] = celli;
        }
    }

    if (nReport == 0)
    {
        Info<< "ATC-T fixed-point tangent check: no requested cells in "
            << "the active design set. Diagnostic only; betaMult is "
            << "unchanged." << endl;
        return;
    }

    scalarField conductivitySens(mesh_.nCells(), Zero);
    {
        const volScalarField& T = TRef();
        const volScalarField& Ta = TaPtr_();
        const autoPtr<incompressible::turbulenceModel>& turbulence =
            primalVars_.turbulence();
        const volScalarField Df
        (
            props_.DFluidField(T, "DFluidFieldFPTangentCheck")
        );
        const volScalarField Ds
        (
            props_.DSolidField(T, "DSolidFieldFPTangentCheck")
        );

        scalarField dJdDEff(mesh_.nCells(), Zero);
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const surfaceScalarField& w = mesh_.weights();
        const surfaceScalarField& deltaCoeffs = mesh_.nonOrthDeltaCoeffs();
        const surfaceScalarField& magSf = mesh_.magSf();

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];

            const scalar dJdDf =
                magSf[facei]*deltaCoeffs[facei]
               *(T[N] - T[P])*(Ta[N] - Ta[P]);

            dJdDEff[P] += w[facei]*dJdDf;
            dJdDEff[N] += (scalar(1) - w[facei])*dJdDf;
        }

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
                dJdDEff[fc[i]] -= Ta[fc[i]]*magSfb[i]*snGradT[i];
            }
        }

        conductivitySens =
            thermalSensScale_
           *dJdDEff/V
           *(
                Ds.primitiveField() - Df.primitiveField()
              - turbulence->nut()().primitiveField()/Prt_
            )
           *dt;

        vars.sourceTermSensitivities
        (
            conductivitySens,
            kInterpolation_(),
            scalar(1),
            designVariablesName,
            "beta"
        );
    }

    const volVectorField UBase
    (
        IOobject
        (
            "ATCTFPTangentUBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );
    const volScalarField pBase
    (
        IOobject
        (
            "ATCTFPTangentpBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p
    );
    const surfaceScalarField phiBase
    (
        IOobject
        (
            "ATCTFPTangentPhiBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );
    const volScalarField betaBase
    (
        IOobject
        (
            "ATCTFPTangentBetaBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        beta
    );
    const volScalarField TBase
    (
        IOobject
        (
            "ATCTFPTangentTBase",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        TRef()
    );

    const bool pressureNeedsReference = p.needReference();
    const label pRefCell = solverControl_().pRefCell();
    const scalar pRefValue = solverControl_().pRefValue();

    auto removeMean = [](scalarField& fld)
    {
        const label n =
            returnReduce(fld.size(), sumOp<label>());

        if (n > 0)
        {
            fld -= gSum(fld)/scalar(n);
        }
    };

    Info<< "ATC-T fixed-point tangent pressure reference: "
        << "needReference " << pressureNeedsReference
        << " pRefCell " << pRefCell
        << " pRefValue " << pRefValue
        << endl;

    forAll(p.boundaryField(), patchi)
    {
        if (mesh_.boundary()[patchi].type() == "empty")
        {
            continue;
        }

        Info<< "ATC-T fixed-point tangent pressure patch: "
            << mesh_.boundary()[patchi].name()
            << " type " << p.boundaryField()[patchi].type()
            << " fixesValue " << p.boundaryField()[patchi].fixesValue()
            << endl;
    }

    auto normalizePressureState = [&](SimpleMapState& state)
    {
        if (!pressureNeedsReference)
        {
            return;
        }

        if (pRefCell >= 0 && pRefCell < mesh_.nCells())
        {
            const scalar shift = state.p[pRefCell] - pRefValue;
            state.p -= shift;
        }
        else
        {
            removeMean(state.p);
        }
    };

    auto projectPressureDirection = [&](scalarField& dp)
    {
        if (!pressureNeedsReference)
        {
            return;
        }

        if (pRefCell >= 0 && pRefCell < mesh_.nCells())
        {
            dp[pRefCell] = Zero;
        }
        else
        {
            removeMean(dp);
        }
    };

    auto projectDirection = [&](SimpleMapState& state)
    {
        projectPressureDirection(state.p);
    };

    auto enforcePressureBoundaryState = [&]()
    {
        bool hasUpdateablePressureSnGrad = false;

        forAll(p.boundaryField(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            if
            (
                isA<updateablePatchTypes::updateableSnGrad>
                (
                    p.boundaryField()[patchi]
                )
            )
            {
                hasUpdateablePressureSnGrad = true;
            }
        }

        forAll(p.boundaryFieldRef(), patchi)
        {
            if
            (
                mesh_.boundary()[patchi].type() == "empty"
             || mesh_.boundary()[patchi].coupled()
            )
            {
                continue;
            }

            if (pBase.boundaryField()[patchi].fixesValue())
            {
                p.boundaryFieldRef()[patchi] ==
                    pBase.boundaryField()[patchi];
            }
        }

        if (!hasUpdateablePressureSnGrad)
        {
            p.correctBoundaryConditions();
        }
    };

    auto restoreState = [&]()
    {
        U.primitiveFieldRef() = UBase.primitiveField();
        forAll(U.boundaryFieldRef(), patchi)
        {
            U.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
        }
        U.correctBoundaryConditions();

        p.primitiveFieldRef() = pBase.primitiveField();
        forAll(p.boundaryFieldRef(), patchi)
        {
            p.boundaryFieldRef()[patchi] == pBase.boundaryField()[patchi];
        }
        enforcePressureBoundaryState();

        phi.primitiveFieldRef() = phiBase.primitiveField();
        forAll(phi.boundaryFieldRef(), patchi)
        {
            phi.boundaryFieldRef()[patchi] ==
                phiBase.boundaryField()[patchi];
        }
    };

    auto restoreBeta = [&]()
    {
        beta.primitiveFieldRef() = betaBase.primitiveField();
        forAll(beta.boundaryFieldRef(), patchi)
        {
            beta.boundaryFieldRef()[patchi] ==
                betaBase.boundaryField()[patchi];
        }
        beta.correctBoundaryConditions();
    };

    auto restoreTemperature = [&]()
    {
        volScalarField& T = const_cast<volScalarField&>(TRef());
        T.primitiveFieldRef() = TBase.primitiveField();
        forAll(T.boundaryFieldRef(), patchi)
        {
            T.boundaryFieldRef()[patchi] == TBase.boundaryField()[patchi];
        }
        T.correctBoundaryConditions();
    };

    auto setStatePerturbation =
        [&](const scalar coeff, const SimpleMapState& direction)
    {
        if (direction.phiInternal.size() != phi.primitiveField().size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiInternal size mismatch in "
                << "setStatePerturbation: "
                << direction.phiInternal.size() << " vs "
                << phi.primitiveField().size()
                << exit(FatalError);
        }
        if (direction.phiBoundary.size() != phi.boundaryField().size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiBoundary list size mismatch in "
                << "setStatePerturbation: "
                << direction.phiBoundary.size() << " vs "
                << phi.boundaryField().size()
                << exit(FatalError);
        }

        U.primitiveFieldRef() =
            UBase.primitiveField() + coeff*direction.U;
        p.primitiveFieldRef() =
            pBase.primitiveField() + coeff*direction.p;
        phi.primitiveFieldRef() = phiBase.primitiveField();
        phi.primitiveFieldRef() += coeff*direction.phiInternal;

        forAll(U.boundaryFieldRef(), patchi)
        {
            U.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
        }
        U.correctBoundaryConditions();

        forAll(p.boundaryFieldRef(), patchi)
        {
            p.boundaryFieldRef()[patchi] == pBase.boundaryField()[patchi];
        }
        forAll(phi.boundaryFieldRef(), patchi)
        {
            phi.boundaryFieldRef()[patchi] ==
                phiBase.boundaryField()[patchi];

            if
            (
                direction.phiBoundary[patchi].size()
             != phi.boundaryField()[patchi].size()
            )
            {
                FatalErrorInFunction
                    << "SimpleMapState phiBoundary patch size mismatch in "
                    << "setStatePerturbation on patch " << patchi << ": "
                    << direction.phiBoundary[patchi].size() << " vs "
                    << phi.boundaryField()[patchi].size()
                    << exit(FatalError);
            }

            phi.boundaryFieldRef()[patchi] +=
                coeff*direction.phiBoundary[patchi];
        }
        enforcePressureBoundaryState();
    };

    auto stateZero = [&](SimpleMapState& state)
    {
        state.U = vector::zero;
        state.p = Zero;
        state.phiInternal = Zero;
        forAll(state.phiBoundary, patchi)
        {
            state.phiBoundary[patchi] = Zero;
        }
    };

    auto checkPhiStateSizes =
        [&](const SimpleMapState& a, const SimpleMapState& b, const word& op)
    {
        if (a.phiInternal.size() != b.phiInternal.size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiInternal size mismatch in " << op
                << ": " << a.phiInternal.size() << " vs "
                << b.phiInternal.size()
                << exit(FatalError);
        }
        if (a.phiBoundary.size() != b.phiBoundary.size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiBoundary list size mismatch in " << op
                << ": " << a.phiBoundary.size() << " vs "
                << b.phiBoundary.size()
                << exit(FatalError);
        }
        forAll(a.phiBoundary, patchi)
        {
            if (a.phiBoundary[patchi].size() != b.phiBoundary[patchi].size())
            {
                FatalErrorInFunction
                    << "SimpleMapState phiBoundary patch size mismatch in "
                    << op << " on patch " << patchi << ": "
                    << a.phiBoundary[patchi].size() << " vs "
                    << b.phiBoundary[patchi].size()
                    << exit(FatalError);
            }
        }
    };

    auto stateAxpy =
        [&](SimpleMapState& y, const scalar a, const SimpleMapState& x)
    {
        checkPhiStateSizes(y, x, "stateAxpy");

        forAll(y.U, celli)
        {
            y.U[celli] += a*x.U[celli];
            y.p[celli] += a*x.p[celli];
        }
        y.phiInternal += a*x.phiInternal;
        forAll(y.phiBoundary, patchi)
        {
            y.phiBoundary[patchi] += a*x.phiBoundary[patchi];
        }
    };

    auto stateScale = [&](SimpleMapState& y, const scalar a)
    {
        forAll(y.U, celli)
        {
            y.U[celli] *= a;
            y.p[celli] *= a;
        }
        y.phiInternal *= a;
        forAll(y.phiBoundary, patchi)
        {
            y.phiBoundary[patchi] *= a;
        }
    };

    auto stateDot =
        [&]
        (
            const SimpleMapState& a,
            const SimpleMapState& b,
            const scalar UScale,
            const scalar pScale,
            const scalar phiScale
        )
    {
        checkPhiStateSizes(a, b, "stateDot");

        scalar val = Zero;
        const scalar invUSqr = scalar(1)/sqr(max(UScale, VSMALL));
        const scalar invpSqr = scalar(1)/sqr(max(pScale, VSMALL));
        const scalar invPhiSqr = scalar(1)/sqr(max(phiScale, VSMALL));

        forAll(a.U, celli)
        {
            val += (a.U[celli] & b.U[celli])*invUSqr;
            val += a.p[celli]*b.p[celli]*invpSqr;
        }
        forAll(a.phiInternal, facei)
        {
            val +=
                a.phiInternal[facei]
               *b.phiInternal[facei]
               *invPhiSqr;
        }
        forAll(a.phiBoundary, patchi)
        {
            forAll(a.phiBoundary[patchi], facei)
            {
                val +=
                    a.phiBoundary[patchi][facei]
                   *b.phiBoundary[patchi][facei]
                   *invPhiSqr;
            }
        }

        return returnReduce(val, sumOp<scalar>());
    };

    auto stateNorm =
        [&]
        (
            const SimpleMapState& a,
            const scalar UScale,
            const scalar pScale,
            const scalar phiScale
        )
    {
        return sqrt(max(stateDot(a, a, UScale, pScale, phiScale), Zero));
    };

    auto stateBlockNorms =
        [&]
        (
            const SimpleMapState& a,
            scalar& UNorm,
            scalar& pNorm,
            scalar& phiInternalNorm,
            scalar& phiBoundaryNorm
        )
    {
        scalar USqr = Zero;
        scalar pSqr = Zero;
        scalar phiInternalSqr = Zero;
        scalar phiBoundarySqr = Zero;

        forAll(a.U, celli)
        {
            USqr += magSqr(a.U[celli]);
            pSqr += sqr(a.p[celli]);
        }
        forAll(a.phiInternal, facei)
        {
            phiInternalSqr += sqr(a.phiInternal[facei]);
        }
        forAll(a.phiBoundary, patchi)
        {
            forAll(a.phiBoundary[patchi], facei)
            {
                phiBoundarySqr += sqr(a.phiBoundary[patchi][facei]);
            }
        }

        UNorm = sqrt(max(returnReduce(USqr, sumOp<scalar>()), Zero));
        pNorm = sqrt(max(returnReduce(pSqr, sumOp<scalar>()), Zero));
        phiInternalNorm =
            sqrt(max(returnReduce(phiInternalSqr, sumOp<scalar>()), Zero));
        phiBoundaryNorm =
            sqrt(max(returnReduce(phiBoundarySqr, sumOp<scalar>()), Zero));
    };

    auto physicalSeedContraction =
        [&](const SimpleMapState& state, scalar& internalFlow, scalar& boundaryFlow)
    {
        if (state.phiInternal.size() != gPhi.primitiveField().size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiInternal size mismatch in "
                << "physicalSeedContraction: "
                << state.phiInternal.size() << " vs "
                << gPhi.primitiveField().size()
                << exit(FatalError);
        }
        if (state.phiBoundary.size() != gPhi.boundaryField().size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiBoundary list size mismatch in "
                << "physicalSeedContraction: "
                << state.phiBoundary.size() << " vs "
                << gPhi.boundaryField().size()
                << exit(FatalError);
        }

        internalFlow = Zero;
        boundaryFlow = Zero;

        forAll(state.phiInternal, facei)
        {
            internalFlow +=
                couplingSign_
               *gPhi.primitiveField()[facei]
               *state.phiInternal[facei];
        }
        forAll(state.phiBoundary, patchi)
        {
            if
            (
                state.phiBoundary[patchi].size()
             != gPhi.boundaryField()[patchi].size()
            )
            {
                FatalErrorInFunction
                    << "SimpleMapState phiBoundary patch size mismatch in "
                    << "physicalSeedContraction on patch " << patchi << ": "
                    << state.phiBoundary[patchi].size() << " vs "
                    << gPhi.boundaryField()[patchi].size()
                    << exit(FatalError);
            }

            forAll(state.phiBoundary[patchi], facei)
            {
                boundaryFlow +=
                    couplingSign_
                   *gPhi.boundaryField()[patchi][facei]
                   *state.phiBoundary[patchi][facei];
            }
        }

        internalFlow = returnReduce(internalFlow, sumOp<scalar>());
        boundaryFlow = returnReduce(boundaryFlow, sumOp<scalar>());
        return internalFlow + boundaryFlow;
    };

    auto mapStateDifference =
        [&](const SimpleMapState& plus, const SimpleMapState& minus, const scalar eps)
    {
        checkPhiStateSizes(plus, minus, "mapStateDifference");

        SimpleMapState diff(mesh_);

        forAll(diff.U, celli)
        {
            diff.U[celli] = (plus.U[celli] - minus.U[celli])/(2*eps);
            diff.p[celli] = (plus.p[celli] - minus.p[celli])/(2*eps);
        }
        forAll(diff.phiInternal, facei)
        {
            diff.phiInternal[facei] =
                (plus.phiInternal[facei] - minus.phiInternal[facei])
               /(2*eps);
        }
        forAll(diff.phiBoundary, patchi)
        {
            forAll(diff.phiBoundary[patchi], facei)
            {
                diff.phiBoundary[patchi][facei] =
                    (
                        plus.phiBoundary[patchi][facei]
                      - minus.phiBoundary[patchi][facei]
                    )
                   /(2*eps);
            }
        }
        projectDirection(diff);

        return diff;
    };

    auto baseStateFromFields = [&]()
    {
        SimpleMapState state(mesh_);

        state.U = UBase.primitiveField();
        state.p = pBase.primitiveField();
        state.phiInternal = phiBase.primitiveField();
        forAll(state.phiBoundary, patchi)
        {
            state.phiBoundary[patchi] = phiBase.boundaryField()[patchi];
        }

        normalizePressureState(state);

        return state;
    };

    auto installState = [&](const SimpleMapState& state)
    {
        if (state.U.size() != U.primitiveField().size())
        {
            FatalErrorInFunction
                << "SimpleMapState U size mismatch in installState: "
                << state.U.size() << " vs " << U.primitiveField().size()
                << exit(FatalError);
        }
        if (state.p.size() != p.primitiveField().size())
        {
            FatalErrorInFunction
                << "SimpleMapState p size mismatch in installState: "
                << state.p.size() << " vs " << p.primitiveField().size()
                << exit(FatalError);
        }
        if (state.phiInternal.size() != phi.primitiveField().size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiInternal size mismatch in "
                << "installState: " << state.phiInternal.size()
                << " vs " << phi.primitiveField().size()
                << exit(FatalError);
        }
        if (state.phiBoundary.size() != phi.boundaryField().size())
        {
            FatalErrorInFunction
                << "SimpleMapState phiBoundary list size mismatch in "
                << "installState: " << state.phiBoundary.size()
                << " vs " << phi.boundaryField().size()
                << exit(FatalError);
        }

        U.primitiveFieldRef() = state.U;
        forAll(U.boundaryFieldRef(), patchi)
        {
            U.boundaryFieldRef()[patchi] == UBase.boundaryField()[patchi];
        }
        U.correctBoundaryConditions();

        p.primitiveFieldRef() = state.p;
        forAll(p.boundaryFieldRef(), patchi)
        {
            p.boundaryFieldRef()[patchi] == pBase.boundaryField()[patchi];
        }
        enforcePressureBoundaryState();

        phi.primitiveFieldRef() = state.phiInternal;
        forAll(phi.boundaryFieldRef(), patchi)
        {
            if (state.phiBoundary[patchi].size() != phi.boundaryField()[patchi].size())
            {
                FatalErrorInFunction
                    << "SimpleMapState phiBoundary patch size mismatch in "
                    << "installState on patch " << patchi << ": "
                    << state.phiBoundary[patchi].size() << " vs "
                    << phi.boundaryField()[patchi].size()
                    << exit(FatalError);
            }

            phi.boundaryFieldRef()[patchi] == state.phiBoundary[patchi];
        }
    };

    auto stateDifference =
        [&](const SimpleMapState& a, const SimpleMapState& b)
    {
        SimpleMapState diff(a);
        stateAxpy(diff, scalar(-1), b);
        projectDirection(diff);
        return diff;
    };

    auto rawBlockScale = [&](const SimpleMapState& state)
    {
        scalar UNorm = Zero;
        scalar pNorm = Zero;
        scalar phiInternalNorm = Zero;
        scalar phiBoundaryNorm = Zero;
        stateBlockNorms
        (
            state,
            UNorm,
            pNorm,
            phiInternalNorm,
            phiBoundaryNorm
        );

        return max
        (
            sqrt
            (
                sqr(UNorm)
              + sqr(pNorm)
              + sqr(phiInternalNorm)
              + sqr(phiBoundaryNorm)
            ),
            VSMALL
        );
    };

    auto iterateFullSimpleMap =
        [&]
        (
            const label cellj,
            const scalar betaValue,
            const SimpleMapState& initialState,
            const scalar UScale,
            const scalar pScale,
            const scalar phiScale,
            label& nIters,
            scalar& finalRel
        )
    {
        restoreState();
        restoreBeta();

        beta.primitiveFieldRef()[cellj] = betaValue;
        beta.correctBoundaryConditions();

        SimpleMapState current(initialState);
        normalizePressureState(current);
        nIters = 0;
        finalRel = VGREAT;

        const label maxIters =
            max(label(1), fixedPointValidationMaxIters_);

        for (label iter = 0; iter < maxIters; ++iter)
        {
            installState(current);

            SimpleMapState next =
                primalSimpleMapStateAtFrozenState("ATCFPValidationIter");
            normalizePressureState(next);

            SimpleMapState change = stateDifference(next, current);
            finalRel =
                stateNorm(change, UScale, pScale, phiScale)
               /max(stateNorm(next, UScale, pScale, phiScale), VSMALL);

            current = next;
            nIters = iter + 1;

            if (finalRel < fixedPointValidationRelTol_)
            {
                break;
            }
        }

        restoreState();
        restoreBeta();

        return current;
    };

    auto thermalStateRelChange =
        [&](const volScalarField& TNew, const scalarField& oldT)
    {
        scalar diff = Zero;
        scalar norm = VSMALL;

        forAll(oldT, celli)
        {
            diff += sqr(TNew[celli] - oldT[celli]);
            norm += sqr(TNew[celli]);
        }

        reduce(diff, sumOp<scalar>());
        reduce(norm, sumOp<scalar>());

        return sqrt(diff)/max(sqrt(norm), VSMALL);
    };

    auto solveDiagnosticTEqnOnce = [&]()
    {
        const surfaceScalarField& phiT = primalVars_.phi();
        volScalarField& T = const_cast<volScalarField&>(TRef());
        fv::options& fvOptions(fv::options::New(this->mesh_));

        const volScalarField DEff(this->DEff());

        tmp<fvScalarMatrix> tTEqn;
        if (props_.variableRhoCp())
        {
            const volScalarField C(props_.CField(T, "rhoCpNormFPTangentFD"));
            tTEqn =
            (
                C.internalField()*fvm::div(phiT, T)
              - fvm::laplacian(DEff, T)
             ==
                fvOptions(T)
            );
        }
        else
        {
            tTEqn =
            (
                fvm::div(phiT, T)
              - fvm::laplacian(DEff, T)
             ==
                fvOptions(T)
            );
        }

        fvScalarMatrix& TEqn = tTEqn.ref();
        TEqn.relax();
        fvOptions.constrain(TEqn);
        TEqn.solve().initialResidual();
        fvOptions.correct(T);
        T.correctBoundaryConditions();
    };

    auto convergeDiagnosticTemperature =
        [&](label& nIters, scalar& finalRel)
    {
        volScalarField& T = const_cast<volScalarField&>(TRef());
        nIters = 0;
        finalRel = VGREAT;

        const label maxIters =
            max(label(1), fixedPointValidationMaxIters_);

        for (label iter = 0; iter < maxIters; ++iter)
        {
            const scalarField oldT(T.primitiveField());
            solveDiagnosticTEqnOnce();

            finalRel = thermalStateRelChange(T, oldT);
            nIters = iter + 1;

            if (finalRel < fixedPointValidationRelTol_)
            {
                break;
            }
        }
    };

    auto weightedObjectiveValue = [&]()
    {
        scalar J = Zero;
        PtrList<objective>& functions =
            objectiveManager_.getObjectiveFunctions();

        for (objective& func : functions)
        {
            J += func.weight()*func.J();
        }

        return J;
    };

    auto objectiveAtConvergedFlow =
        [&]
        (
            const label cellj,
            const scalar betaValue,
            const SimpleMapState& state,
            label& thermalIters,
            scalar& thermalRel
        )
    {
        restoreState();
        restoreBeta();
        restoreTemperature();

        installState(state);
        beta.primitiveFieldRef()[cellj] = betaValue;
        beta.correctBoundaryConditions();

        convergeDiagnosticTemperature(thermalIters, thermalRel);
        const scalar J = weightedObjectiveValue();

        restoreState();
        restoreBeta();
        restoreTemperature();

        return J;
    };

    auto applyMxWithEps =
        [&](const SimpleMapState& direction, const scalar eps)
    {
        setStatePerturbation(eps, direction);
        SimpleMapState plus =
            primalSimpleMapStateAtFrozenState("ATCFPTangentMxPlus");

        setStatePerturbation(-eps, direction);
        SimpleMapState minus =
            primalSimpleMapStateAtFrozenState("ATCFPTangentMxMinus");

        restoreState();

        return mapStateDifference(plus, minus, eps);
    };

    auto applyMx =
        [&](const SimpleMapState& direction)
    {
        return applyMxWithEps(direction, fixedPointTangentFDEps_);
    };

    auto applyA =
        [&](const SimpleMapState& direction)
    {
        SimpleMapState result(direction);
        SimpleMapState Mx = applyMx(direction);

        stateAxpy(result, scalar(-1), Mx);
        projectDirection(result);

        return result;
    };

    auto applyAWithEps =
        [&](const SimpleMapState& direction, const scalar eps)
    {
        SimpleMapState result(direction);
        SimpleMapState Mx = applyMxWithEps(direction, eps);

        stateAxpy(result, scalar(-1), Mx);
        projectDirection(result);

        return result;
    };

    auto solveDense =
        [&](List<scalarField>& A, scalarField b, const label n)
    {
        scalarField x(n, Zero);

        for (label k = 0; k < n; ++k)
        {
            label pivot = k;
            scalar pivotMag = mag(A[k][k]);

            for (label i = k + 1; i < n; ++i)
            {
                if (mag(A[i][k]) > pivotMag)
                {
                    pivot = i;
                    pivotMag = mag(A[i][k]);
                }
            }

            if (pivot != k)
            {
                scalarField tmpRow(A[k]);
                A[k] = A[pivot];
                A[pivot] = tmpRow;

                const scalar tmpB = b[k];
                b[k] = b[pivot];
                b[pivot] = tmpB;
            }

            const scalar diag = A[k][k];
            if (mag(diag) <= VSMALL)
            {
                continue;
            }

            for (label i = k + 1; i < n; ++i)
            {
                const scalar factor = A[i][k]/diag;
                A[i][k] = Zero;

                for (label j = k + 1; j < n; ++j)
                {
                    A[i][j] -= factor*A[k][j];
                }
                b[i] -= factor*b[k];
            }
        }

        for (label i = n - 1; i >= 0; --i)
        {
            scalar sum = b[i];
            for (label j = i + 1; j < n; ++j)
            {
                sum -= A[i][j]*x[j];
            }

            if (mag(A[i][i]) > VSMALL)
            {
                x[i] = sum/A[i][i];
            }
        }

        return x;
    };

    auto gmresSolve =
        [&]
        (
            const label reportCell,
            const SimpleMapState& rhs,
            const scalar UScale,
            const scalar pScale,
            const scalar phiScale,
            label& nIters,
            scalar& estimatedRel,
            scalar& finalRel,
            scalarField& residualBlockRel,
            scalarField& rhsBlockRel
        )
    {
        const label maxIters = max(label(1), fixedPointTangentMaxIters_);

        SimpleMapState x(mesh_);
        stateZero(x);

        const scalar rhsNorm =
            max(stateNorm(rhs, UScale, pScale, phiScale), VSMALL);

        scalar rhsUNorm = Zero;
        scalar rhsPNorm = Zero;
        scalar rhsPhiInternalNorm = Zero;
        scalar rhsPhiBoundaryNorm = Zero;
        stateBlockNorms
        (
            rhs,
            rhsUNorm,
            rhsPNorm,
            rhsPhiInternalNorm,
            rhsPhiBoundaryNorm
        );
        const scalar rhsRawNorm =
            max
            (
                sqrt
                (
                    sqr(rhsUNorm)
                  + sqr(rhsPNorm)
                  + sqr(rhsPhiInternalNorm)
                  + sqr(rhsPhiBoundaryNorm)
                ),
                VSMALL
            );
        rhsBlockRel.setSize(4);
        rhsBlockRel[0] = rhsUNorm/rhsRawNorm;
        rhsBlockRel[1] = rhsPNorm/rhsRawNorm;
        rhsBlockRel[2] = rhsPhiInternalNorm/rhsRawNorm;
        rhsBlockRel[3] = rhsPhiBoundaryNorm/rhsRawNorm;

        List<SimpleMapState> basis
        (
            maxIters + 1,
            SimpleMapState(mesh_)
        );
        List<scalarField> H(maxIters + 1);
        forAll(H, row)
        {
            H[row].setSize(maxIters);
            H[row] = Zero;
        }

        basis[0] = rhs;
        stateScale(basis[0], scalar(1)/rhsNorm);

        {
            SimpleMapState A3 =
                applyAWithEps(basis[0], scalar(3e-6));
            SimpleMapState A1 =
                applyAWithEps(basis[0], scalar(1e-6));
            SimpleMapState A03 =
                applyAWithEps(basis[0], scalar(3e-7));

            SimpleMapState diff31(A3);
            stateAxpy(diff31, scalar(-1), A1);
            SimpleMapState diff103(A1);
            stateAxpy(diff103, scalar(-1), A03);

            const scalar A1Norm =
                max(stateNorm(A1, UScale, pScale, phiScale), VSMALL);

            Info<< "ATC-T fixed-point tangent Aeps check: cell "
                << reportCell
                << " rel(A3e-6-A1e-6) = "
                << stateNorm(diff31, UScale, pScale, phiScale)/A1Norm
                << ", rel(A1e-6-A3e-7) = "
                << stateNorm(diff103, UScale, pScale, phiScale)/A1Norm
                << endl;
        }

        scalarField bestY(0);
        label usedIters = 0;

        for (label k = 0; k < maxIters; ++k)
        {
            SimpleMapState w = applyA(basis[k]);

            for (label i = 0; i <= k; ++i)
            {
                H[i][k] =
                    stateDot(w, basis[i], UScale, pScale, phiScale);
                stateAxpy(w, -H[i][k], basis[i]);
            }

            H[k + 1][k] = stateNorm(w, UScale, pScale, phiScale);
            if (H[k + 1][k] > VSMALL && k + 1 < maxIters + 1)
            {
                basis[k + 1] = w;
                stateScale(basis[k + 1], scalar(1)/H[k + 1][k]);
            }

            const label m = k + 1;
            List<scalarField> normal(m);
            forAll(normal, i)
            {
                normal[i].setSize(m);
                normal[i] = Zero;
            }
            scalarField normalRhs(m, Zero);

            for (label i = 0; i < m; ++i)
            {
                normalRhs[i] = H[0][i]*rhsNorm;

                for (label j = 0; j < m; ++j)
                {
                    scalar accum = Zero;
                    for (label row = 0; row <= m; ++row)
                    {
                        accum += H[row][i]*H[row][j];
                    }
                    normal[i][j] = accum;
                }

                normal[i][i] += VSMALL;
            }

            bestY = solveDense(normal, normalRhs, m);
            usedIters = m;

            scalar lsNormSqr = Zero;
            for (label row = 0; row <= m; ++row)
            {
                scalar residual = (row == 0 ? rhsNorm : Zero);
                for (label col = 0; col < m; ++col)
                {
                    residual -= H[row][col]*bestY[col];
                }
                lsNormSqr += sqr(residual);
            }

            estimatedRel = sqrt(lsNormSqr)/rhsNorm;

            Info<< "ATC-T fixed-point tangent GMRES iter: cell "
                << reportCell
                << " iteration " << m
                << " estimatedRel " << estimatedRel
                << " hNext " << H[k + 1][k]
                << endl;

            const bool doTrueResidual =
                m == 5
             || m == 10
             || m == 15
             || m == maxIters
             || estimatedRel < fixedPointTangentRelTol_;

            if (doTrueResidual)
            {
                SimpleMapState xCheck(mesh_);
                stateZero(xCheck);
                for (label i = 0; i < m; ++i)
                {
                    stateAxpy(xCheck, bestY[i], basis[i]);
                }
                projectDirection(xCheck);

                SimpleMapState AxCheck = applyA(xCheck);
                SimpleMapState residualCheck(rhs);
                stateAxpy(residualCheck, scalar(-1), AxCheck);

                const scalar trueRel =
                    stateNorm
                    (
                        residualCheck,
                        UScale,
                        pScale,
                        phiScale
                    )/rhsNorm;

                scalar resUNorm = Zero;
                scalar resPNorm = Zero;
                scalar resPhiInternalNorm = Zero;
                scalar resPhiBoundaryNorm = Zero;
                stateBlockNorms
                (
                    residualCheck,
                    resUNorm,
                    resPNorm,
                    resPhiInternalNorm,
                    resPhiBoundaryNorm
                );

                Info<< "ATC-T fixed-point tangent GMRES true checkpoint: "
                    << "cell " << reportCell
                    << " iteration " << m
                    << " trueRel " << trueRel
                    << " residualUFraction " << resUNorm/rhsRawNorm
                    << " residualPFraction " << resPNorm/rhsRawNorm
                    << " residualPhiInternalFraction "
                    << resPhiInternalNorm/rhsRawNorm
                    << " residualPhiBoundaryFraction "
                    << resPhiBoundaryNorm/rhsRawNorm
                    << endl;
            }

            if (estimatedRel < fixedPointTangentRelTol_)
            {
                break;
            }
        }

        stateZero(x);
        for (label i = 0; i < usedIters; ++i)
        {
            stateAxpy(x, bestY[i], basis[i]);
        }
        projectDirection(x);

        SimpleMapState Ax = applyA(x);
        SimpleMapState residual(rhs);
        stateAxpy(residual, scalar(-1), Ax);
        finalRel = stateNorm(residual, UScale, pScale, phiScale)/rhsNorm;
        nIters = usedIters;

        scalar resUNorm = Zero;
        scalar resPNorm = Zero;
        scalar resPhiInternalNorm = Zero;
        scalar resPhiBoundaryNorm = Zero;
        stateBlockNorms
        (
            residual,
            resUNorm,
            resPNorm,
            resPhiInternalNorm,
            resPhiBoundaryNorm
        );
        residualBlockRel.setSize(4);
        residualBlockRel[0] = resUNorm/rhsRawNorm;
        residualBlockRel[1] = resPNorm/rhsRawNorm;
        residualBlockRel[2] = resPhiInternalNorm/rhsRawNorm;
        residualBlockRel[3] = resPhiBoundaryNorm/rhsRawNorm;

        Info<< "ATC-T fixed-point tangent GMRES solution sizes: U="
            << x.U.size()
            << ", p=" << x.p.size()
            << ", phiInternal=" << x.phiInternal.size()
            << ", phiPatches=" << x.phiBoundary.size()
            << endl;

        return x;
    };

    scalar UScale = max(gMax(mag(UBase.primitiveField())), SMALL);
    scalarField pBaseGauge(pBase.primitiveField());
    if (pressureNeedsReference)
    {
        if (pRefCell >= 0 && pRefCell < mesh_.nCells())
        {
            const scalar shift = pBaseGauge[pRefCell] - pRefValue;
            pBaseGauge -= shift;
        }
        else
        {
            removeMean(pBaseGauge);
        }
    }
    scalar pScale = max(gMax(mag(pBaseGauge)), SMALL);
    scalar phiScale = max(gMax(mag(phiBase.primitiveField())), SMALL);
    forAll(phiBase.boundaryField(), patchi)
    {
        if (phiBase.boundaryField()[patchi].size())
        {
            phiScale =
                max(phiScale, gMax(mag(phiBase.boundaryField()[patchi])));
        }
    }

    restoreState();
    restoreBeta();

    SimpleMapState baseState = baseStateFromFields();
    installState(baseState);
    SimpleMapState mapBase =
        primalSimpleMapStateAtFrozenState("ATCFPBaseDefect");
    normalizePressureState(mapBase);
    SimpleMapState baseDefect = stateDifference(mapBase, baseState);

    const scalar baseDefectRel =
        stateNorm(baseDefect, UScale, pScale, phiScale)
       /max(stateNorm(mapBase, UScale, pScale, phiScale), VSMALL);

    scalar baseDefectUNorm = Zero;
    scalar baseDefectPNorm = Zero;
    scalar baseDefectPhiInternalNorm = Zero;
    scalar baseDefectPhiBoundaryNorm = Zero;
    stateBlockNorms
    (
        baseDefect,
        baseDefectUNorm,
        baseDefectPNorm,
        baseDefectPhiInternalNorm,
        baseDefectPhiBoundaryNorm
    );
    const scalar baseRawScale = rawBlockScale(baseState);

    Info<< "ATC-T fixed-point tangent base defect: "
        << "baseDefectRel " << baseDefectRel
        << " baseDefectU " << baseDefectUNorm/baseRawScale
        << " baseDefectP " << baseDefectPNorm/baseRawScale
        << " baseDefectPhiInternal "
        << baseDefectPhiInternalNorm/baseRawScale
        << " baseDefectPhiBoundary "
        << baseDefectPhiBoundaryNorm/baseRawScale
        << endl;

    restoreState();
    restoreBeta();

    if (fixedPointTangentMaxIters_ <= 0)
    {
        Info<< "ATC-T fixed-point tangent check: base-defect-only mode. "
            << "Diagnostic only; betaMult is unchanged."
            << endl;
        return;
    }

    Info<< "ATC-T fixed-point tangent samples: "
        << "cell beta epsBeta gmresIters estimatedRel finalTrueRel "
        << "rhsU rhsP rhsPhiInternal rhsPhiBoundary "
        << "resU resP resPhiInternal resPhiBoundary "
        << "flowDerivativePhiInternal/dt flowDerivativePhiBoundary/dt "
        << "flowDerivative/dt "
        << "flowDerivativeVol/dt conductivityVol/dt "
        << "totalDerivativeVol/dt totalIntegrated/dt "
        << "physicalBetaObjectiveFD/dt signedRatio"
        << endl;

    restoreState();
    restoreBeta();

    for (label reporti = 0; reporti < nReport; ++reporti)
    {
        const label cellj = reportCells[reporti];
        const scalar oldBeta = beta[cellj];
        const scalar room = min(oldBeta, scalar(1) - oldBeta);
        const scalar epsBeta =
            min(fixedPointTangentFDEps_, scalar(0.25)*max(SMALL, room));

        if (epsBeta <= SMALL)
        {
            continue;
        }

        restoreState();
        beta.primitiveFieldRef()[cellj] = oldBeta + epsBeta;
        beta.correctBoundaryConditions();
        SimpleMapState statePlus =
            primalSimpleMapStateAtFrozenState("ATCFPTangentBetaPlus");

        restoreState();
        beta.primitiveFieldRef()[cellj] = oldBeta - epsBeta;
        beta.correctBoundaryConditions();
        SimpleMapState stateMinus =
            primalSimpleMapStateAtFrozenState("ATCFPTangentBetaMinus");

        restoreState();
        restoreBeta();

        SimpleMapState rhs =
            mapStateDifference(statePlus, stateMinus, epsBeta);

        label gmresIters = 0;
        scalar estimatedRel = VGREAT;
        scalar finalRel = VGREAT;
        scalarField residualBlockRel(4, Zero);
        scalarField rhsBlockRel(4, Zero);
        SimpleMapState delta =
            gmresSolve
            (
                cellj,
                rhs,
                UScale,
                pScale,
                phiScale,
                gmresIters,
                estimatedRel,
                finalRel,
                residualBlockRel,
                rhsBlockRel
            );

        scalar flowInternal = Zero;
        scalar flowBoundary = Zero;
        const scalar flowDerivative =
            physicalSeedContraction(delta, flowInternal, flowBoundary);
        const scalar flowDerivativeVol = flowDerivative/V[cellj];
        const scalar totalDerivativeVol =
            flowDerivativeVol + conductivitySens[cellj];
        const scalar totalIntegrated =
            flowDerivative + V[cellj]*conductivitySens[cellj];

        Info<< "ATC-T fixed-point tangent sample: "
            << cellj << " "
            << oldBeta << " "
            << epsBeta << " "
            << gmresIters << " "
            << estimatedRel << " "
            << finalRel << " "
            << rhsBlockRel[0] << " "
            << rhsBlockRel[1] << " "
            << rhsBlockRel[2] << " "
            << rhsBlockRel[3] << " "
            << residualBlockRel[0] << " "
            << residualBlockRel[1] << " "
            << residualBlockRel[2] << " "
            << residualBlockRel[3] << " "
            << flowInternal/dt << " "
            << flowBoundary/dt << " "
            << flowDerivative/dt << " "
            << flowDerivativeVol/dt << " "
            << conductivitySens[cellj]/dt << " "
            << totalDerivativeVol/dt << " "
            << totalIntegrated/dt << " "
            << "notComputed "
            << "notComputed"
            << endl;

        if (checkFixedPointTangentAgainstFD_)
        {
            scalarList validationEps(4, Zero);
            validationEps[0] = scalar(10)*fixedPointValidationBetaEps_;
            validationEps[1] =
                (scalar(10)/scalar(3))*fixedPointValidationBetaEps_;
            validationEps[2] = fixedPointValidationBetaEps_;
            validationEps[3] =
                (scalar(1)/scalar(3))*fixedPointValidationBetaEps_;

            Info<< "ATC-T fixed-point tangent FD validation samples: "
                << "cell epsBeta plusIterations plusFinalRel "
                << "minusIterations minusFinalRel "
                << "tangentVsFDRelTotal tangentVsFDRelU "
                << "tangentVsFDRelP tangentVsFDRelPhiInternal "
                << "tangentVsFDRelPhiBoundary "
                << "flowGMRES/dt flowFD/dt flowSignedRatio "
                << "flowRelativeDifference flowGMRESVol/dt flowFDVol/dt"
                << endl;

            Info<< "ATC-T fixed-point tangent objective FD validation "
                << "samples: cell epsBeta plusTIters plusTRel "
                << "minusTIters minusTRel Jplus Jminus objectiveFDVol/dt "
                << "flowDerivativeVol/dt conductivityVol/dt "
                << "tangentTotalVol/dt signedRatio relativeDifference"
                << endl;

            Info<< "ATC-T fixed-point tangent conductivity-only FD "
                << "validation samples: cell epsBeta plusTIters plusTRel "
                << "minusTIters minusTRel Jplus Jminus "
                << "conductivityFDVol/dt analyticConductivityVol/dt "
                << "signedRatio relativeDifference"
                << endl;

            Info<< "ATC-T fixed-point tangent flow-only objective FD "
                << "validation samples: cell epsBeta plusTIters plusTRel "
                << "minusTIters minusTRel Jplus Jminus flowObjectiveFDVol/dt "
                << "gPhiFlowFDVol/dt signedRatio relativeDifference"
                << endl;

            Info<< "ATC-T scalar flux residual direction checks: "
                << "cell epsBeta mode assembledTa/V/dt gPhi/V/dt "
                << "coupledGphi/V/dt flowObjectiveFDVol/dt "
                << "ratioGphiToAssembled ratioAssembledToObjectiveFD "
                << "ratioGphiToObjectiveFD"
                << endl;

            Info<< "ATC-T scalar flux residual patch checks: "
                << "cell epsBeta patch TpatchType split nFaces "
                << "nSignCross assembledTa/V/dt gPhi/V/dt"
                << endl;

            Info<< "ATC-T scalar tangent flow checks: "
                << "cell epsBeta objectiveLinear/V/dt assembledTa/V/dt "
                << "flowObjectiveFDVol/dt ratioLinearToAssembled "
                << "ratioLinearToObjectiveFD ratioAssembledToObjectiveFD"
                << endl;

            Info<< "ATC-T sign-preserving scalar flux checks: "
                << "cell epsBeta h nearZeroPhi mode nIncluded nExcluded "
                << "nSignCross assembledTa/V/dt gPhi/V/dt ratioGphiToAssembled"
                << endl;

            Info<< "ATC-T sign-preserving excluded flux contribution: "
                << "cell epsBeta nearZeroPhi nExcluded gPhiOriginal/V/dt"
                << endl;

            auto makePhiProbe =
                [&]
                (
                    const word& name,
                    const SimpleMapState& state,
                    const bool useInternal,
                    const bool useBoundary,
                    const label patchFilter,
                    const label signFilter
                )
            {
                surfaceScalarField phiProbe
                (
                    IOobject
                    (
                        name,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    phiBase
                );

                if (useInternal)
                {
                    phiProbe.primitiveFieldRef() = state.phiInternal;
                }

                surfaceScalarField::Boundary& phib =
                    phiProbe.boundaryFieldRef();

                forAll(phib, patchi)
                {
                    phib[patchi] == phiBase.boundaryField()[patchi];

                    if (!useBoundary)
                    {
                        continue;
                    }
                    if (patchFilter != -1 && patchi != patchFilter)
                    {
                        continue;
                    }

                    const scalarField& basePatch =
                        phiBase.boundaryField()[patchi];

                    forAll(phib[patchi], facei)
                    {
                        const bool includeFace =
                            signFilter == 0
                         || (signFilter < 0 && basePatch[facei] < 0)
                         || (signFilter > 0 && basePatch[facei] >= 0);

                        if (includeFace)
                        {
                            phib[patchi][facei] =
                                state.phiBoundary[patchi][facei];
                        }
                    }
                }

                return phiProbe;
            };

            auto thermalConvIntegratedResidual =
                [&](const surfaceScalarField& phiProbe)
            {
                volScalarField& TMutable = const_cast<volScalarField&>(TRef());

                tmp<fvScalarMatrix> tConvEqn;
                if (props_.variableRhoCp())
                {
                    const volScalarField C
                    (
                        props_.CField(TMutable, "rhoCpNormFluxResidualCheck")
                    );

                    tConvEqn =
                        C.internalField()
                       *fvm::div(phiProbe, TMutable, "div(phi,T)");
                }
                else
                {
                    tConvEqn = fvm::div(phiProbe, TMutable, "div(phi,T)");
                }

                fvScalarMatrix& convEqn = tConvEqn.ref();
                scalarField residual(-convEqn.residual());

                return residual;
            };

            auto thermalFluxResidualDifference =
                [&]
                (
                    const SimpleMapState& plusState,
                    const SimpleMapState& minusState,
                    const scalar eps,
                    const bool useInternal,
                    const bool useBoundary,
                    const label patchFilter,
                    const label signFilter
                )
            {
                surfaceScalarField phiPlus =
                    makePhiProbe
                    (
                        "ATCTScalarFluxResidualPhiPlus",
                        plusState,
                        useInternal,
                        useBoundary,
                        patchFilter,
                        signFilter
                    );
                surfaceScalarField phiMinus =
                    makePhiProbe
                    (
                        "ATCTScalarFluxResidualPhiMinus",
                        minusState,
                        useInternal,
                        useBoundary,
                        patchFilter,
                        signFilter
                    );

                const scalarField Rplus =
                    thermalConvIntegratedResidual(phiPlus);
                const scalarField Rminus =
                    thermalConvIntegratedResidual(phiMinus);

                scalarField dRphi(Rplus.size(), Zero);
                forAll(dRphi, celli)
                {
                    dRphi[celli] =
                        (Rplus[celli] - Rminus[celli])/(2*eps);
                }

                return dRphi;
            };

            auto assembledFluxResidualContraction =
                [&]
                (
                    const SimpleMapState& plusState,
                    const SimpleMapState& minusState,
                    const scalar eps,
                    const bool useInternal,
                    const bool useBoundary,
                    const label patchFilter,
                    const label signFilter
                )
            {
                surfaceScalarField phiPlus =
                    makePhiProbe
                    (
                        "ATCTScalarFluxResidualPhiPlus",
                        plusState,
                        useInternal,
                        useBoundary,
                        patchFilter,
                        signFilter
                    );
                surfaceScalarField phiMinus =
                    makePhiProbe
                    (
                        "ATCTScalarFluxResidualPhiMinus",
                        minusState,
                        useInternal,
                        useBoundary,
                        patchFilter,
                        signFilter
                    );

                const scalarField Rplus =
                    thermalConvIntegratedResidual(phiPlus);
                const scalarField Rminus =
                    thermalConvIntegratedResidual(phiMinus);

                const volScalarField& Ta = TaPtr_();
                scalar contraction = Zero;

                forAll(Rplus, celli)
                {
                    contraction +=
                        Ta[celli]*(Rplus[celli] - Rminus[celli])/(2*eps);
                }

                return returnReduce(contraction, sumOp<scalar>());
            };

            auto gPhiFluxContraction =
                [&]
                (
                    const SimpleMapState& plusState,
                    const SimpleMapState& minusState,
                    const scalar eps,
                    const bool useInternal,
                    const bool useBoundary,
                    const label patchFilter,
                    const label signFilter
                )
            {
                scalar contraction = Zero;

                if (useInternal)
                {
                    forAll(gPhi.primitiveField(), facei)
                    {
                        contraction +=
                            gPhi.primitiveField()[facei]
                           *(
                                plusState.phiInternal[facei]
                              - minusState.phiInternal[facei]
                            )
                           /(2*eps);
                    }
                }

                if (useBoundary)
                {
                    forAll(gPhi.boundaryField(), patchi)
                    {
                        if (patchFilter != -1 && patchi != patchFilter)
                        {
                            continue;
                        }

                        const fvsPatchScalarField& gPatch =
                            gPhi.boundaryField()[patchi];
                        const scalarField& basePatch =
                            phiBase.boundaryField()[patchi];

                        forAll(gPatch, facei)
                        {
                            const bool includeFace =
                                signFilter == 0
                             || (signFilter < 0 && basePatch[facei] < 0)
                             || (signFilter > 0 && basePatch[facei] >= 0);

                            if (includeFace)
                            {
                                contraction +=
                                    gPatch[facei]
                                   *(
                                        plusState.phiBoundary[patchi][facei]
                                      - minusState.phiBoundary[patchi][facei]
                                    )
                                   /(2*eps);
                            }
                        }
                    }
                }

                return returnReduce(contraction, sumOp<scalar>());
            };

            auto printScalarFluxResidualCheck =
                [&]
                (
                    const label cellj,
                    const scalar eps,
                    const word& mode,
                    const scalar assembled,
                    const scalar gPhiTerm,
                    const scalar objectiveFDVol
                )
            {
                const scalar assembledVol = assembled/V[cellj]/dt;
                const scalar gPhiVol = gPhiTerm/V[cellj]/dt;
                const scalar coupledGphiVol = couplingSign_*gPhiVol;
                const scalar ratioGphiToAssembled =
                    mag(assembledVol) > VSMALL ? gPhiVol/assembledVol : VGREAT;
                const scalar ratioAssembledToObjective =
                    mag(objectiveFDVol) > VSMALL
                  ? assembledVol/objectiveFDVol
                  : VGREAT;
                const scalar ratioGphiToObjective =
                    mag(objectiveFDVol) > VSMALL
                  ? gPhiVol/objectiveFDVol
                  : VGREAT;

                Info<< "ATC-T scalar flux residual direction check: "
                    << cellj << " "
                    << eps << " "
                    << mode << " "
                    << assembledVol << " "
                    << gPhiVol << " "
                    << coupledGphiVol << " "
                    << objectiveFDVol << " "
                    << ratioGphiToAssembled << " "
                    << ratioAssembledToObjective << " "
                    << ratioGphiToObjective
                    << endl;
            };

            auto objectiveTangentContraction =
                [&](const volScalarField& dT)
            {
                scalar objectiveVolume = Zero;
                scalar objectiveBoundary = Zero;
                PtrList<objective>& functions =
                    objectiveManager_.getObjectiveFunctions();

                const scalarField& V = mesh_.V();

                for (objective& func : functions)
                {
                    objectiveIncompressible& funcI =
                        refCast<objectiveIncompressible>(func);

                    if (funcI.hasdJdT())
                    {
                        const volScalarField& dJdT = funcI.dJdT();
                        forAll(dT.primitiveField(), celli)
                        {
                            objectiveVolume +=
                                func.weight()*V[celli]*dJdT[celli]*dT[celli];
                        }
                    }

                    if (!funcI.hasBoundarydJdT())
                    {
                        continue;
                    }

                    forAll(mesh_.boundary(), patchi)
                    {
                        if (mesh_.boundary()[patchi].type() == "empty")
                        {
                            continue;
                        }

                        const scalarField& bdJdT = funcI.boundarydJdT(patchi);
                        const fvPatchScalarField& dTp =
                            dT.boundaryField()[patchi];
                        const fvsPatchScalarField& magSfp =
                            mesh_.magSf().boundaryField()[patchi];

                        forAll(bdJdT, facei)
                        {
                            objectiveBoundary +=
                                func.weight()
                               *bdJdT[facei]
                               *dTp[facei]
                               *magSfp[facei];
                        }
                    }
                }

                return returnReduce
                (
                    objectiveVolume + objectiveBoundary,
                    sumOp<scalar>()
                );
            };

            auto solveScalarTangent =
                [&](const scalarField& dRphi)
            {
                const volScalarField& T = TRef();

                volScalarField dT
                (
                    IOobject
                    (
                        "ATCTScalarTangentdT",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    T
                );
                dT = dimensionedScalar(T.dimensions(), Zero);

                volScalarField::Boundary& dTb = dT.boundaryFieldRef();
                forAll(dTb, patchi)
                {
                    if (mesh_.boundary()[patchi].type() == "empty")
                    {
                        continue;
                    }

                    if (isA<fixedGradientFvPatchScalarField>(dTb[patchi]))
                    {
                        refCast<fixedGradientFvPatchScalarField>
                        (
                            dTb[patchi]
                        ).gradient() =
                            scalarField(dTb[patchi].size(), Zero);
                    }

                    if (dTb[patchi].fixesValue())
                    {
                        dTb[patchi] == scalar(0);
                    }
                }
                dT.correctBoundaryConditions();
                forAll(dTb, patchi)
                {
                    if
                    (
                        mesh_.boundary()[patchi].type() != "empty"
                     && dTb[patchi].fixesValue()
                    )
                    {
                        dTb[patchi] == scalar(0);
                    }
                }

                const volScalarField DEffBase(this->DEff());

                tmp<fvScalarMatrix> tdTEqn;
                if (props_.variableRhoCp())
                {
                    const volScalarField C
                    (
                        props_.CField(T, "rhoCpNormScalarTangentCheck")
                    );

                    tdTEqn =
                    (
                        C.internalField()*fvm::div(phiBase, dT, "div(phi,T)")
                      - fvm::laplacian(DEffBase, dT)
                    );
                }
                else
                {
                    tdTEqn =
                    (
                        fvm::div(phiBase, dT, "div(phi,T)")
                      - fvm::laplacian(DEffBase, dT)
                    );
                }

                fvScalarMatrix& dTEqn = tdTEqn.ref();
                dTEqn.source() = -dRphi;
                dTEqn.solve("T");
                dT.correctBoundaryConditions();

                return dT;
            };

            forAll(validationEps, epsi)
            {
                const scalar validationRoom =
                    scalar(0.25)*max(SMALL, room);
                const scalar validationBetaEps =
                    min(validationEps[epsi], validationRoom);

                if (validationBetaEps <= SMALL)
                {
                    continue;
                }

                label plusIters = 0;
                scalar plusRel = VGREAT;
                SimpleMapState xPlus =
                    iterateFullSimpleMap
                    (
                        cellj,
                        oldBeta + validationBetaEps,
                        baseState,
                        UScale,
                        pScale,
                        phiScale,
                        plusIters,
                        plusRel
                    );

                label minusIters = 0;
                scalar minusRel = VGREAT;
                SimpleMapState xMinus =
                    iterateFullSimpleMap
                    (
                        cellj,
                        oldBeta - validationBetaEps,
                        baseState,
                        UScale,
                        pScale,
                        phiScale,
                        minusIters,
                        minusRel
                    );

                SimpleMapState deltaFD =
                    mapStateDifference(xPlus, xMinus, validationBetaEps);
                SimpleMapState error = stateDifference(delta, deltaFD);

                const scalar fdNorm =
                    max(stateNorm(deltaFD, UScale, pScale, phiScale), VSMALL);
                const scalar errorRel =
                    stateNorm(error, UScale, pScale, phiScale)/fdNorm;

                scalar errorUNorm = Zero;
                scalar errorPNorm = Zero;
                scalar errorPhiInternalNorm = Zero;
                scalar errorPhiBoundaryNorm = Zero;
                stateBlockNorms
                (
                    error,
                    errorUNorm,
                    errorPNorm,
                    errorPhiInternalNorm,
                    errorPhiBoundaryNorm
                );

                scalar fdUNorm = Zero;
                scalar fdPNorm = Zero;
                scalar fdPhiInternalNorm = Zero;
                scalar fdPhiBoundaryNorm = Zero;
                stateBlockNorms
                (
                    deltaFD,
                    fdUNorm,
                    fdPNorm,
                    fdPhiInternalNorm,
                    fdPhiBoundaryNorm
                );

                scalar flowFDInternal = Zero;
                scalar flowFDBoundary = Zero;
                const scalar flowFD =
                    physicalSeedContraction
                    (
                        deltaFD,
                        flowFDInternal,
                        flowFDBoundary
                    );

                const scalar flowRatio =
                    mag(flowFD) > VSMALL ? flowDerivative/flowFD : VGREAT;
                const scalar flowRelDiff =
                    mag(flowDerivative - flowFD)/max(mag(flowFD), VSMALL);

                label plusTIters = 0;
                scalar plusTRel = VGREAT;
                const scalar Jplus =
                    objectiveAtConvergedFlow
                    (
                        cellj,
                        oldBeta + validationBetaEps,
                        xPlus,
                        plusTIters,
                        plusTRel
                    );

                label minusTIters = 0;
                scalar minusTRel = VGREAT;
                const scalar Jminus =
                    objectiveAtConvergedFlow
                    (
                        cellj,
                        oldBeta - validationBetaEps,
                        xMinus,
                        minusTIters,
                        minusTRel
                    );

                const scalar objectiveFDVol =
                    (Jplus - Jminus)/(2*validationBetaEps*V[cellj]);
                const scalar tangentTotalVol = totalDerivativeVol/dt;
                const scalar objectiveRatio =
                    mag(objectiveFDVol) > VSMALL
                  ? tangentTotalVol/objectiveFDVol
                  : VGREAT;
                const scalar objectiveRelDiff =
                    mag(tangentTotalVol - objectiveFDVol)
                   /max(mag(objectiveFDVol), VSMALL);

                Info<< "ATC-T fixed-point tangent FD validation sample: "
                    << cellj << " "
                    << validationBetaEps << " "
                    << plusIters << " "
                    << plusRel << " "
                    << minusIters << " "
                    << minusRel << " "
                    << errorRel << " "
                    << errorUNorm/max(fdUNorm, VSMALL) << " "
                    << errorPNorm/max(fdPNorm, VSMALL) << " "
                    << errorPhiInternalNorm
                       /max(fdPhiInternalNorm, VSMALL) << " "
                    << errorPhiBoundaryNorm
                       /max(fdPhiBoundaryNorm, VSMALL) << " "
                    << flowDerivative/dt << " "
                    << flowFD/dt << " "
                    << flowRatio << " "
                    << flowRelDiff << " "
                    << flowDerivativeVol/dt << " "
                    << flowFD/V[cellj]/dt
                    << endl;

                Info<< "ATC-T fixed-point tangent objective FD validation "
                    << "sample: "
                    << cellj << " "
                    << validationBetaEps << " "
                    << plusTIters << " "
                    << plusTRel << " "
                    << minusTIters << " "
                    << minusTRel << " "
                    << Jplus << " "
                    << Jminus << " "
                    << objectiveFDVol << " "
                    << flowDerivativeVol/dt << " "
                    << conductivitySens[cellj]/dt << " "
                    << tangentTotalVol << " "
                    << objectiveRatio << " "
                    << objectiveRelDiff
                    << endl;

                label plusCondTIters = 0;
                scalar plusCondTRel = VGREAT;
                const scalar JplusCondOnly =
                    objectiveAtConvergedFlow
                    (
                        cellj,
                        oldBeta + validationBetaEps,
                        baseState,
                        plusCondTIters,
                        plusCondTRel
                    );

                label minusCondTIters = 0;
                scalar minusCondTRel = VGREAT;
                const scalar JminusCondOnly =
                    objectiveAtConvergedFlow
                    (
                        cellj,
                        oldBeta - validationBetaEps,
                        baseState,
                        minusCondTIters,
                        minusCondTRel
                    );

                const scalar conductivityFDVol =
                    (JplusCondOnly - JminusCondOnly)
                   /(2*validationBetaEps*V[cellj]);
                const scalar analyticConductivityVol =
                    conductivitySens[cellj]/dt;
                const scalar conductivityRatio =
                    mag(conductivityFDVol) > VSMALL
                  ? analyticConductivityVol/conductivityFDVol
                  : VGREAT;
                const scalar conductivityRelDiff =
                    mag(analyticConductivityVol - conductivityFDVol)
                   /max(mag(conductivityFDVol), VSMALL);

                Info<< "ATC-T fixed-point tangent conductivity-only FD "
                    << "validation sample: "
                    << cellj << " "
                    << validationBetaEps << " "
                    << plusCondTIters << " "
                    << plusCondTRel << " "
                    << minusCondTIters << " "
                    << minusCondTRel << " "
                    << JplusCondOnly << " "
                    << JminusCondOnly << " "
                    << conductivityFDVol << " "
                    << analyticConductivityVol << " "
                    << conductivityRatio << " "
                    << conductivityRelDiff
                    << endl;

                label plusFlowOnlyTIters = 0;
                scalar plusFlowOnlyTRel = VGREAT;
                const scalar JplusFlowOnly =
                    objectiveAtConvergedFlow
                    (
                        cellj,
                        oldBeta,
                        xPlus,
                        plusFlowOnlyTIters,
                        plusFlowOnlyTRel
                    );

                label minusFlowOnlyTIters = 0;
                scalar minusFlowOnlyTRel = VGREAT;
                const scalar JminusFlowOnly =
                    objectiveAtConvergedFlow
                    (
                        cellj,
                        oldBeta,
                        xMinus,
                        minusFlowOnlyTIters,
                        minusFlowOnlyTRel
                    );

                const scalar flowObjectiveFDVol =
                    (JplusFlowOnly - JminusFlowOnly)
                   /(2*validationBetaEps*V[cellj]);
                const scalar gPhiFlowFDVol = flowFD/V[cellj]/dt;
                const scalar flowObjectiveRatio =
                    mag(flowObjectiveFDVol) > VSMALL
                  ? gPhiFlowFDVol/flowObjectiveFDVol
                  : VGREAT;
                const scalar flowObjectiveRelDiff =
                    mag(gPhiFlowFDVol - flowObjectiveFDVol)
                   /max(mag(flowObjectiveFDVol), VSMALL);

                Info<< "ATC-T fixed-point tangent flow-only objective FD "
                    << "validation sample: "
                    << cellj << " "
                    << validationBetaEps << " "
                    << plusFlowOnlyTIters << " "
                    << plusFlowOnlyTRel << " "
                    << minusFlowOnlyTIters << " "
                    << minusFlowOnlyTRel << " "
                    << JplusFlowOnly << " "
                    << JminusFlowOnly << " "
                    << flowObjectiveFDVol << " "
                    << gPhiFlowFDVol << " "
                    << flowObjectiveRatio << " "
                    << flowObjectiveRelDiff
                    << endl;

                restoreTemperature();
                restoreBeta();
                installState(baseState);

                label nInternalSignCrossings = 0;
                forAll(xPlus.phiInternal, facei)
                {
                    if
                    (
                        (
                            xPlus.phiInternal[facei] > 0
                         && xMinus.phiInternal[facei] < 0
                        )
                     || (
                            xPlus.phiInternal[facei] < 0
                         && xMinus.phiInternal[facei] > 0
                        )
                    )
                    {
                        ++nInternalSignCrossings;
                    }
                }
                reduce(nInternalSignCrossings, sumOp<label>());

                label nBoundarySignCrossings = 0;
                forAll(xPlus.phiBoundary, patchi)
                {
                    forAll(xPlus.phiBoundary[patchi], facei)
                    {
                        if
                        (
                            (
                                xPlus.phiBoundary[patchi][facei] > 0
                             && xMinus.phiBoundary[patchi][facei] < 0
                            )
                         || (
                                xPlus.phiBoundary[patchi][facei] < 0
                             && xMinus.phiBoundary[patchi][facei] > 0
                            )
                        )
                        {
                            ++nBoundarySignCrossings;
                        }
                    }
                }
                reduce(nBoundarySignCrossings, sumOp<label>());

                Info<< "ATC-T scalar flux residual direction setup: "
                    << "cell " << cellj
                    << " epsBeta " << validationBetaEps
                    << " divScheme " << mesh_.divScheme("div(phi,T)")
                    << " boundedConvection " << boundedConvection_
                    << " nInternalSignCrossings "
                    << nInternalSignCrossings
                    << " nBoundarySignCrossings "
                    << nBoundarySignCrossings
                    << endl;

                const scalar assembledFull =
                    assembledFluxResidualContraction
                    (
                        xPlus,
                        xMinus,
                        validationBetaEps,
                        true,
                        true,
                        -1,
                        0
                    );
                const scalar gPhiFull =
                    gPhiFluxContraction
                    (
                        xPlus,
                        xMinus,
                        validationBetaEps,
                        true,
                        true,
                        -1,
                        0
                    );
                printScalarFluxResidualCheck
                (
                    cellj,
                    validationBetaEps,
                    "full",
                    assembledFull,
                    gPhiFull,
                    flowObjectiveFDVol
                );

                const scalar assembledInternal =
                    assembledFluxResidualContraction
                    (
                        xPlus,
                        xMinus,
                        validationBetaEps,
                        true,
                        false,
                        -1,
                        0
                    );
                const scalar gPhiInternalOnly =
                    gPhiFluxContraction
                    (
                        xPlus,
                        xMinus,
                        validationBetaEps,
                        true,
                        false,
                        -1,
                        0
                    );
                printScalarFluxResidualCheck
                (
                    cellj,
                    validationBetaEps,
                    "internal",
                    assembledInternal,
                    gPhiInternalOnly,
                    flowObjectiveFDVol
                );

                const scalar assembledBoundary =
                    assembledFluxResidualContraction
                    (
                        xPlus,
                        xMinus,
                        validationBetaEps,
                        false,
                        true,
                        -1,
                        0
                    );
                const scalar gPhiBoundaryOnly =
                    gPhiFluxContraction
                    (
                        xPlus,
                        xMinus,
                        validationBetaEps,
                        false,
                        true,
                        -1,
                        0
                    );
                printScalarFluxResidualCheck
                (
                    cellj,
                    validationBetaEps,
                    "boundary",
                    assembledBoundary,
                    gPhiBoundaryOnly,
                    flowObjectiveFDVol
                );

                forAll(mesh_.boundary(), patchi)
                {
                    const fvPatch& patch = mesh_.boundary()[patchi];
                    if (patch.type() == "empty" || patch.size() == 0)
                    {
                        continue;
                    }

                    for (label signFilter = -1; signFilter <= 1; signFilter += 2)
                    {
                        label nFaces = 0;
                        label nSignCross = 0;
                        const scalarField& basePatch =
                            phiBase.boundaryField()[patchi];

                        forAll(basePatch, facei)
                        {
                            const bool includeFace =
                                (signFilter < 0 && basePatch[facei] < 0)
                             || (signFilter > 0 && basePatch[facei] >= 0);

                            if (!includeFace)
                            {
                                continue;
                            }

                            ++nFaces;
                            if
                            (
                                (
                                    xPlus.phiBoundary[patchi][facei] > 0
                                 && xMinus.phiBoundary[patchi][facei] < 0
                                )
                             || (
                                    xPlus.phiBoundary[patchi][facei] < 0
                                 && xMinus.phiBoundary[patchi][facei] > 0
                                )
                            )
                            {
                                ++nSignCross;
                            }
                        }

                        reduce(nFaces, sumOp<label>());
                        reduce(nSignCross, sumOp<label>());

                        if (nFaces == 0)
                        {
                            continue;
                        }

                        const scalar assembledPatch =
                            assembledFluxResidualContraction
                            (
                                xPlus,
                                xMinus,
                                validationBetaEps,
                                false,
                                true,
                                patchi,
                                signFilter
                            );
                        const scalar gPhiPatch =
                            gPhiFluxContraction
                            (
                                xPlus,
                                xMinus,
                                validationBetaEps,
                                false,
                                true,
                                patchi,
                                signFilter
                            );

                        Info<< "ATC-T scalar flux residual patch check: "
                            << cellj << " "
                            << validationBetaEps << " "
                            << patch.name() << " "
                            << TRef().boundaryField()[patchi].type() << " "
                            << (signFilter < 0 ? "inflow" : "outflow") << " "
                            << nFaces << " "
                            << nSignCross << " "
                            << assembledPatch/V[cellj]/dt << " "
                            << gPhiPatch/V[cellj]/dt
                            << endl;
                    }
                }

                scalar maxPhi = gMax(mag(phiBase.primitiveField()));
                forAll(phiBase.boundaryField(), patchi)
                {
                    if (phiBase.boundaryField()[patchi].size())
                    {
                        maxPhi =
                            max
                            (
                                maxPhi,
                                gMax(mag(phiBase.boundaryField()[patchi]))
                            );
                    }
                }
                reduce(maxPhi, maxOp<scalar>());

                const scalar nearZeroPhi =
                    max(scalar(1e-12)*max(maxPhi, scalar(1)), VSMALL);

                scalar signRoom = VGREAT;
                label nIncludedInternal = 0;
                label nExcludedInternal = 0;
                label nIncludedBoundary = 0;
                label nExcludedBoundary = 0;

                forAll(deltaFD.phiInternal, facei)
                {
                    const scalar dphif = deltaFD.phiInternal[facei];
                    if (mag(dphif) <= SMALL)
                    {
                        continue;
                    }

                    if (mag(phiBase.primitiveField()[facei]) <= nearZeroPhi)
                    {
                        ++nExcludedInternal;
                    }
                    else
                    {
                        ++nIncludedInternal;
                        signRoom =
                            min
                            (
                                signRoom,
                                mag(phiBase.primitiveField()[facei])
                               /mag(dphif)
                            );
                    }
                }

                forAll(deltaFD.phiBoundary, patchi)
                {
                    const scalarField& basePatch =
                        phiBase.boundaryField()[patchi];

                    forAll(deltaFD.phiBoundary[patchi], facei)
                    {
                        const scalar dphif =
                            deltaFD.phiBoundary[patchi][facei];

                        if (mag(dphif) <= SMALL)
                        {
                            continue;
                        }

                        if (mag(basePatch[facei]) <= nearZeroPhi)
                        {
                            ++nExcludedBoundary;
                        }
                        else
                        {
                            ++nIncludedBoundary;
                            signRoom =
                                min
                                (
                                    signRoom,
                                    mag(basePatch[facei])/mag(dphif)
                                );
                        }
                    }
                }

                reduce(signRoom, minOp<scalar>());
                reduce(nIncludedInternal, sumOp<label>());
                reduce(nExcludedInternal, sumOp<label>());
                reduce(nIncludedBoundary, sumOp<label>());
                reduce(nExcludedBoundary, sumOp<label>());

                const label nIncludedTotal =
                    nIncludedInternal + nIncludedBoundary;
                const scalar signH =
                    nIncludedTotal
                  ? min(scalar(1e-5), scalar(0.1)*signRoom)
                  : scalar(1e-5);

                SimpleMapState signPlus(baseState);
                SimpleMapState signMinus(baseState);
                SimpleMapState excludedPlus(baseState);
                SimpleMapState excludedMinus(baseState);

                label nSignCrossInternal = 0;
                label nSignCrossBoundary = 0;
                scalar gPhiExcluded = Zero;

                forAll(deltaFD.phiInternal, facei)
                {
                    const scalar dphif = deltaFD.phiInternal[facei];
                    if (mag(dphif) <= SMALL)
                    {
                        continue;
                    }

                    const scalar basePhi = phiBase.primitiveField()[facei];

                    if (mag(basePhi) <= nearZeroPhi)
                    {
                        excludedPlus.phiInternal[facei] =
                            basePhi + signH*dphif;
                        excludedMinus.phiInternal[facei] =
                            basePhi - signH*dphif;
                        gPhiExcluded +=
                            gPhi.primitiveField()[facei]*dphif;
                    }
                    else
                    {
                        signPlus.phiInternal[facei] = basePhi + signH*dphif;
                        signMinus.phiInternal[facei] = basePhi - signH*dphif;

                        if
                        (
                            signPlus.phiInternal[facei]
                           *signMinus.phiInternal[facei] < 0
                        )
                        {
                            ++nSignCrossInternal;
                        }
                    }
                }

                forAll(deltaFD.phiBoundary, patchi)
                {
                    const scalarField& basePatch =
                        phiBase.boundaryField()[patchi];
                    const fvsPatchScalarField& gPatch =
                        gPhi.boundaryField()[patchi];

                    forAll(deltaFD.phiBoundary[patchi], facei)
                    {
                        const scalar dphif =
                            deltaFD.phiBoundary[patchi][facei];

                        if (mag(dphif) <= SMALL)
                        {
                            continue;
                        }

                        const scalar basePhi = basePatch[facei];

                        if (mag(basePhi) <= nearZeroPhi)
                        {
                            excludedPlus.phiBoundary[patchi][facei] =
                                basePhi + signH*dphif;
                            excludedMinus.phiBoundary[patchi][facei] =
                                basePhi - signH*dphif;
                            gPhiExcluded += gPatch[facei]*dphif;
                        }
                        else
                        {
                            signPlus.phiBoundary[patchi][facei] =
                                basePhi + signH*dphif;
                            signMinus.phiBoundary[patchi][facei] =
                                basePhi - signH*dphif;

                            if
                            (
                                signPlus.phiBoundary[patchi][facei]
                               *signMinus.phiBoundary[patchi][facei] < 0
                            )
                            {
                                ++nSignCrossBoundary;
                            }
                        }
                    }
                }

                reduce(nSignCrossInternal, sumOp<label>());
                reduce(nSignCrossBoundary, sumOp<label>());
                gPhiExcluded = returnReduce(gPhiExcluded, sumOp<scalar>());

                auto printSignPreservingCheck =
                    [&]
                    (
                        const word& mode,
                        const scalar assembled,
                        const scalar gPhiTerm,
                        const label nIncluded,
                        const label nExcluded,
                        const label nSignCross
                    )
                {
                    const scalar assembledVol = assembled/V[cellj]/dt;
                    const scalar gPhiVol = gPhiTerm/V[cellj]/dt;
                    const scalar ratio =
                        mag(assembledVol) > VSMALL
                      ? gPhiVol/assembledVol
                      : VGREAT;

                    Info<< "ATC-T sign-preserving scalar flux check: "
                        << cellj << " "
                        << validationBetaEps << " "
                        << signH << " "
                        << nearZeroPhi << " "
                        << mode << " "
                        << nIncluded << " "
                        << nExcluded << " "
                        << nSignCross << " "
                        << assembledVol << " "
                        << gPhiVol << " "
                        << ratio
                        << endl;
                };

                const scalar signAssembledFull =
                    assembledFluxResidualContraction
                    (
                        signPlus,
                        signMinus,
                        signH,
                        true,
                        true,
                        -1,
                        0
                    );
                const scalar signGPhiFull =
                    gPhiFluxContraction
                    (
                        signPlus,
                        signMinus,
                        signH,
                        true,
                        true,
                        -1,
                        0
                    );
                printSignPreservingCheck
                (
                    "full",
                    signAssembledFull,
                    signGPhiFull,
                    nIncludedTotal,
                    nExcludedInternal + nExcludedBoundary,
                    nSignCrossInternal + nSignCrossBoundary
                );

                const scalar signAssembledInternal =
                    assembledFluxResidualContraction
                    (
                        signPlus,
                        signMinus,
                        signH,
                        true,
                        false,
                        -1,
                        0
                    );
                const scalar signGPhiInternal =
                    gPhiFluxContraction
                    (
                        signPlus,
                        signMinus,
                        signH,
                        true,
                        false,
                        -1,
                        0
                    );
                printSignPreservingCheck
                (
                    "internal",
                    signAssembledInternal,
                    signGPhiInternal,
                    nIncludedInternal,
                    nExcludedInternal,
                    nSignCrossInternal
                );

                const scalar signAssembledBoundary =
                    assembledFluxResidualContraction
                    (
                        signPlus,
                        signMinus,
                        signH,
                        false,
                        true,
                        -1,
                        0
                    );
                const scalar signGPhiBoundary =
                    gPhiFluxContraction
                    (
                        signPlus,
                        signMinus,
                        signH,
                        false,
                        true,
                        -1,
                        0
                    );
                printSignPreservingCheck
                (
                    "boundary",
                    signAssembledBoundary,
                    signGPhiBoundary,
                    nIncludedBoundary,
                    nExcludedBoundary,
                    nSignCrossBoundary
                );

                Info<< "ATC-T sign-preserving excluded flux contribution: "
                    << cellj << " "
                    << validationBetaEps << " "
                    << nearZeroPhi << " "
                    << nExcludedInternal + nExcludedBoundary << " "
                    << gPhiExcluded/V[cellj]/dt
                    << endl;

                const scalarField dRphiFull =
                    thermalFluxResidualDifference
                    (
                        xPlus,
                        xMinus,
                        validationBetaEps,
                        true,
                        true,
                        -1,
                        0
                    );
                volScalarField dTFlow = solveScalarTangent(dRphiFull);
                const scalar objectiveLinear =
                    objectiveTangentContraction(dTFlow);

                const scalar objectiveLinearVol =
                    objectiveLinear/V[cellj]/dt;
                const scalar assembledFullVol =
                    assembledFull/V[cellj]/dt;
                const scalar ratioLinearToAssembled =
                    mag(assembledFullVol) > VSMALL
                  ? objectiveLinearVol/assembledFullVol
                  : VGREAT;
                const scalar ratioLinearToObjective =
                    mag(flowObjectiveFDVol) > VSMALL
                  ? objectiveLinearVol/flowObjectiveFDVol
                  : VGREAT;
                const scalar ratioAssembledToObjective =
                    mag(flowObjectiveFDVol) > VSMALL
                  ? assembledFullVol/flowObjectiveFDVol
                  : VGREAT;

                Info<< "ATC-T scalar tangent flow check: "
                    << cellj << " "
                    << validationBetaEps << " "
                    << objectiveLinearVol << " "
                    << assembledFullVol << " "
                    << flowObjectiveFDVol << " "
                    << ratioLinearToAssembled << " "
                    << ratioLinearToObjective << " "
                    << ratioAssembledToObjective
                    << endl;
            }
        }
    }

    restoreState();
    restoreBeta();
    restoreTemperature();

    Info<< "ATC-T fixed-point tangent check: maxIters = "
        << fixedPointTangentMaxIters_
        << ", relTol = " << fixedPointTangentRelTol_
        << ", eps = " << fixedPointTangentFDEps_
        << ". Diagnostic only; betaMult is unchanged."
        << endl;
}


void Foam::thermalAdjointSimple::addProjectionCoeffFDSensitivity
(
    scalarField& betaMult,
    const word& designVariablesName,
    const scalar dt
)
{
    if (!projectionCoeffFDTopOSens_)
    {
        return;
    }

    if (!mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        Info<< "ATC-T projection-coeff FD check: no topOVars object" << endl;
        return;
    }

    if (designVariablesName != "beta")
    {
        Info<< "ATC-T projection-coeff FD check: perturbing physical beta "
            << "while designVariablesName is " << designVariablesName << endl;
    }

    const scalarField& V = mesh_.V();
    const topOVariablesBase& vars =
        mesh_.lookupObject<topOVariablesBase>("topOVars");
    const topOZones& zones = vars.getTopOZones();

    bitSet isDesign(mesh_.nCells(), false);
    if (!zones.adjointPorousZoneIDs().empty())
    {
        for (const label zoneID : zones.adjointPorousZoneIDs())
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.set(celli);
            }
        }
    }
    else
    {
        isDesign = bitSet(mesh_.nCells(), true);
    }

    for (const label zoneID : zones.fixedPorousZoneIDs())
    {
        for (const label celli : mesh_.cellZones()[zoneID])
        {
            isDesign.unset(celli);
        }
    }
    for (const label zoneID : zones.fixedZeroPorousZoneIDs())
    {
        for (const label celli : mesh_.cellZones()[zoneID])
        {
            isDesign.unset(celli);
        }
    }
    for (const label celli : zones.IOCells())
    {
        isDesign.unset(celli);
    }

    volScalarField& beta =
        const_cast<volScalarField&>(vars.beta());
    tmp<surfaceScalarField> tGPhi = thermalFluxSensitivity();
    const surfaceScalarField& gPhi = tGPhi();

    scalarField coeffSens(mesh_.nCells(), Zero);
    label nChecked = 0;
    label nSkipped = 0;

    FixedList<label, 12> topCells(-1);
    FixedList<scalar, 12> topMag(-1);

    auto addReportCell =
        [&](const label celli, const scalar magValue)
        {
            forAll(topCells, slot)
            {
                if (magValue > topMag[slot])
                {
                    for (label j = topCells.size() - 1; j > slot; --j)
                    {
                        topCells[j] = topCells[j - 1];
                        topMag[j] = topMag[j - 1];
                    }
                    topCells[slot] = celli;
                    topMag[slot] = magValue;
                    break;
                }
            }
        };

    forAll(beta.primitiveField(), cellj)
    {
        if (!isDesign.test(cellj))
        {
            continue;
        }

        const scalar oldBeta = beta[cellj];
        const scalar eps =
            min
            (
                projectionCoeffFDEps_,
                scalar(0.25)*max(SMALL, min(oldBeta, scalar(1) - oldBeta))
            );

        if (eps <= SMALL)
        {
            ++nSkipped;
            continue;
        }

        beta.primitiveFieldRef()[cellj] = oldBeta + eps;
        beta.correctBoundaryConditions();
        tmp<surfaceScalarField> tPhiPlus =
            primalProjectedFluxAtFrozenState("ATCCoeffPlus");
        surfaceScalarField phiPlus(tPhiPlus());

        beta.primitiveFieldRef()[cellj] = oldBeta - eps;
        beta.correctBoundaryConditions();
        tmp<surfaceScalarField> tPhiMinus =
            primalProjectedFluxAtFrozenState("ATCCoeffMinus");
        surfaceScalarField phiMinus(tPhiMinus());

        beta.primitiveFieldRef()[cellj] = oldBeta;
        beta.correctBoundaryConditions();

        scalar contraction = Zero;

        forAll(gPhi.primitiveField(), facei)
        {
            contraction +=
                gPhi[facei]*(phiPlus[facei] - phiMinus[facei])/(2*eps);
        }

        forAll(mesh_.boundary(), patchi)
        {
            if (mesh_.boundary()[patchi].type() == "empty")
            {
                continue;
            }

            const fvsPatchScalarField& gPatch =
                gPhi.boundaryField()[patchi];
            const fvsPatchScalarField& plusPatch =
                phiPlus.boundaryField()[patchi];
            const fvsPatchScalarField& minusPatch =
                phiMinus.boundaryField()[patchi];

            forAll(gPatch, facei)
            {
                contraction +=
                    gPatch[facei]
                   *(plusPatch[facei] - minusPatch[facei])/(2*eps);
            }
        }

        coeffSens[cellj] = couplingSign_*dt*contraction/V[cellj];
        addReportCell(cellj, mag(coeffSens[cellj]));
        ++nChecked;
    }

    beta.correctBoundaryConditions();

    Info<< "ATC-T projection-coeff FD check: nChecked = " << nChecked
        << ", nSkipped = " << nSkipped
        << ", eps = " << projectionCoeffFDEps_
        << ", add = " << addProjectionCoeffFDSens_
        << endl;

    Info<< "ATC-T projection-coeff FD samples: "
        << "cell beta coeffSens/dt"
        << endl;

    forAll(topCells, i)
    {
        const label celli = topCells[i];
        if (celli == -1)
        {
            continue;
        }

        Info<< "ATC-T projection-coeff FD sample: "
            << celli << " "
            << beta[celli] << " "
            << coeffSens[celli]/dt
            << endl;
    }

    if (addProjectionCoeffFDSens_)
    {
        betaMult += coeffSens;
    }
}


void Foam::thermalAdjointSimple::addResidualFDMomentumSensitivity
(
    scalarField& betaMult,
    const scalarField& upstreamContribution,
    const word& designVariablesName,
    const scalar dt
)
{
    if (!residualFDTopOSens_)
    {
        return;
    }

    if (!mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        Info<< "ATC-T residual-FD topO check: no topOVars object" << endl;
        return;
    }

    const incompressibleAdjointVars& adjointVars = getAdjointVars();
    const volVectorField& Ua = adjointVars.Ua();
    const scalarField& V = mesh_.V();
    const topOVariablesBase& vars =
        mesh_.lookupObject<topOVariablesBase>("topOVars");
    const topOZones& zones = vars.getTopOZones();

    if (designVariablesName != "beta")
    {
        Info<< "ATC-T residual-FD topO check: perturbing physical beta while "
            << "designVariablesName is " << designVariablesName << endl;
    }

    volScalarField& beta = const_cast<volScalarField&>(vars.beta());

    bitSet isDesign(mesh_.nCells(), false);
    if (zones.adjointPorousZoneIDs().empty())
    {
        isDesign.fill(true);
    }
    else
    {
        for (const label zoneID : zones.adjointPorousZoneIDs())
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.set(celli);
            }
        }
    }

    auto unsetZone = [&](const label zoneID)
    {
        if (zoneID != -1)
        {
            for (const label celli : mesh_.cellZones()[zoneID])
            {
                isDesign.unset(celli);
            }
        }
    };

    for (const label zoneID : zones.fixedPorousZoneIDs())
    {
        unsetZone(zoneID);
    }
    for (const label zoneID : zones.fixedZeroPorousZoneIDs())
    {
        unsetZone(zoneID);
    }
    for (const label celli : zones.IOCells())
    {
        isDesign.unset(celli);
    }

    scalarField residualMomSens(mesh_.nCells(), Zero);

    scalar sumSqrUpstream = Zero;
    scalar sumSqrDiff = Zero;
    scalar sumSqrDiffFlip = Zero;
    label nChecked = 0;
    label nSkipped = 0;

    const label nReport = 12;
    labelList topCells(nReport, -1);
    scalarField topMetric(nReport, -GREAT);

    auto addReportCell = [&](const label celli, const scalar metric)
    {
        forAll(topCells, i)
        {
            if (metric > topMetric[i])
            {
                for (label j = topCells.size() - 1; j > i; --j)
                {
                    topMetric[j] = topMetric[j - 1];
                    topCells[j] = topCells[j - 1];
                }
                topMetric[i] = metric;
                topCells[i] = celli;
                break;
            }
        }
    };

    forAll(beta.primitiveField(), cellj)
    {
        if (!isDesign.test(cellj))
        {
            continue;
        }

        const scalar oldBeta = beta[cellj];
        const scalar room = min(oldBeta, scalar(1) - oldBeta);

        if (room <= SMALL)
        {
            ++nSkipped;
            continue;
        }

        const scalar eps = min(residualFDEps_, scalar(0.25)*room);

        if (eps <= SMALL)
        {
            ++nSkipped;
            continue;
        }

        beta.primitiveFieldRef()[cellj] = oldBeta + eps;
        beta.correctBoundaryConditions();
        tmp<vectorField> tRPlus = primalMomentumResidualLHS();
        vectorField RPlus(tRPlus());

        beta.primitiveFieldRef()[cellj] = oldBeta - eps;
        beta.correctBoundaryConditions();
        tmp<vectorField> tRMinus = primalMomentumResidualLHS();
        vectorField RMinus(tRMinus());

        beta.primitiveFieldRef()[cellj] = oldBeta;
        beta.correctBoundaryConditions();

        const vectorField dRdbeta((RPlus - RMinus)/(2*eps));

        scalar contraction = Zero;
        forAll(dRdbeta, celli)
        {
            contraction += Ua[celli] & dRdbeta[celli];
        }

        residualMomSens[cellj] = dt*contraction/V[cellj];

        const scalar upstream = upstreamContribution[cellj];
        const scalar residual = residualMomSens[cellj];

        sumSqrUpstream += sqr(upstream);
        sumSqrDiff += sqr(residual - upstream);
        sumSqrDiffFlip += sqr(-residual - upstream);
        addReportCell(cellj, max(mag(upstream), mag(residual)));
        ++nChecked;
    }

    beta.correctBoundaryConditions();

    const scalar l2Scale = sqrt(max(sumSqrUpstream, VSMALL));

    Info<< "ATC-T residual-FD topO check: nChecked = " << nChecked
        << ", nSkipped = " << nSkipped
        << ", relL2(residual-upstream) = " << sqrt(sumSqrDiff)/l2Scale
        << ", relL2(-residual-upstream) = "
        << sqrt(sumSqrDiffFlip)/l2Scale
        << ", eps = " << residualFDEps_
        << ", replace = " << replaceUpstreamMomentumSens_
        << endl;

    Info<< "ATC-T residual-FD topO samples: "
        << "cell beta upstream/dt residualFD/dt ratio ratioFlipped"
        << endl;

    forAll(topCells, i)
    {
        const label celli = topCells[i];
        if (celli == -1)
        {
            continue;
        }

        const scalar upstream = upstreamContribution[celli]/dt;
        const scalar residual = residualMomSens[celli]/dt;
        const scalar scale = (mag(upstream) > VSMALL) ? upstream : VSMALL;

        Info<< "ATC-T residual-FD topO sample: "
            << celli << " "
            << beta[celli] << " "
            << upstream << " "
            << residual << " "
            << residual/scale << " "
            << -residual/scale
            << endl;
    }

    if (replaceUpstreamMomentumSens_)
    {
        betaMult += residualMomSens;
    }
}


void Foam::thermalAdjointSimple::addMomentumSource(fvVectorMatrix& matrix)
{
    adjointSimple::addMomentumSource(matrix);

    if (exactMomentumATC_)
    {
        tmp<volVectorField> tExactATC = exactMomentumATCSource();
        matrix += tExactATC();
    }

    if (exactConvectionTranspose_)
    {
        tmp<volVectorField> tExactConvT = exactConvectionTransposeSource();
        matrix += tExactConvT();
    }

    const volScalarField& T = TRef();
    const volScalarField& Ta = TaPtr_();

    // (ATC-T) the thermal coupling into adjoint momentum: d/du of the primal
    // convective term C(T) u.grad(T), contracted with Ta.
    //
    // Two continuous forms are available, related by
    //     Ta grad(T) = grad(Ta T) - T grad(Ta)
    // so they differ by a pure gradient, which an incompressible adjoint can in
    // principle absorb into the adjoint pressure. That cancellation is exact in
    // the continuum but NOT guaranteed discretely, especially across sharp
    // Brinkman/conductivity interfaces. They are offered as a switch because
    // this is exactly the term under investigation (see
    // examples/coOptimiseChannelAndBody): the two forms are the cheap proxy for
    // "is our continuous source the discrete transpose of fvm::div(phi,T)?".
    if (couplingForm_ == "none")
    {
        // Null test: no thermal coupling into adjoint momentum at all. Used to
        // prove that Ua (and hence the Brinkman sensitivity) really is driven
        // only by this term when there is no flow objective.
        return;
    }

    if (couplingForm_ == "exactFluxTranspose")
    {
        tmp<surfaceScalarField> tGPhi = thermalFluxSensitivity();
        tmp<volVectorField> tSource = projectedFluxMomentumSource(tGPhi());

        matrix += couplingSign_*tSource();

        return;
    }

    tmp<volVectorField> tSource;

    if (couplingForm_ == "negTGradTa")
    {
        tSource = -T*fvc::grad(Ta);
    }
    else    // "TaGradT", the default
    {
        tSource = Ta*fvc::grad(T);
    }

    if (props_.variableRhoCp())
    {
        const volScalarField C(props_.CField(T, "rhoCpNorm"));

        matrix += couplingSign_*(C*tSource())();
    }
    else
    {
        matrix += couplingSign_*tSource();
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
    couplingForm_("TaGradT"),
    boundedConvection_(true),
    checkProjectionTranspose_(false),
    checkScalarThermalTranspose_(false),
    exactScalarThermalTranspose_(false),
    exactMomentumATC_(false),
    exactConvectionTranspose_(false),
    residualFDTopOSens_(false),
    replaceUpstreamMomentumSens_(false),
    residualFDEps_(1e-6),
    projectionCoeffFDTopOSens_(false),
    addProjectionCoeffFDSens_(false),
    projectionCoeffFDEps_(1e-6),
    usePredictorReverseMomentumSens_(false),
    checkOneStepBetaMapSensitivity_(false),
    oneStepBetaMapFDEps_(1e-6),
    checkFixedPointMapAdjoint_(false),
    checkStateMapTranspose_(false),
    checkFullStateMapTranspose_(false),
    checkPressureMapTranspose_(false),
    checkFixedPointTangentSensitivity_(false),
    fixedPointMapAdjointIters_(5),
    fixedPointMapAdjointFDEps_(1e-6),
    fixedPointTangentMaxIters_(30),
    fixedPointTangentRelTol_(1e-7),
    fixedPointTangentFDEps_(1e-6),
    fixedPointTangentCells_({707, 708, 709}),
    checkFixedPointTangentAgainstFD_(false),
    fixedPointValidationMaxIters_(200),
    fixedPointValidationRelTol_(1e-9),
    fixedPointValidationBetaEps_(1e-6),
    kInterpolation_(nullptr)
{
    const dictionary& thermalDict = dict.subDict("thermal");

    props_.read(thermalDict, mesh_);
    Prt_ = thermalDict.getOrDefault<scalar>("Prt", 0.85);
    couplingSign_ = thermalDict.getOrDefault<scalar>("couplingSign", 1);
    thermalSensScale_ =
        thermalDict.getOrDefault<scalar>("thermalSensScale", 1);
    couplingForm_ =
        thermalDict.getOrDefault<word>("couplingForm", "TaGradT");
    boundedConvection_ =
        thermalDict.getOrDefault<bool>("boundedConvection", true);
    checkProjectionTranspose_ =
        thermalDict.getOrDefault<bool>("checkProjectionTranspose", false);
    checkScalarThermalTranspose_ =
        thermalDict.getOrDefault<bool>("checkScalarThermalTranspose", false);
    exactScalarThermalTranspose_ =
        thermalDict.getOrDefault<bool>("exactScalarThermalTranspose", false);
    exactMomentumATC_ =
        thermalDict.getOrDefault<bool>("exactMomentumATC", false);
    exactConvectionTranspose_ =
        thermalDict.getOrDefault<bool>("exactConvectionTranspose", false);
    residualFDTopOSens_ =
        thermalDict.getOrDefault<bool>("residualFDTopOSens", false);
    replaceUpstreamMomentumSens_ =
        thermalDict.getOrDefault<bool>("replaceUpstreamMomentumSens", false);
    residualFDEps_ =
        thermalDict.getOrDefault<scalar>("residualFDEps", 1e-6);
    projectionCoeffFDTopOSens_ =
        thermalDict.getOrDefault<bool>("projectionCoeffFDTopOSens", false);
    addProjectionCoeffFDSens_ =
        thermalDict.getOrDefault<bool>("addProjectionCoeffFDSens", false);
    projectionCoeffFDEps_ =
        thermalDict.getOrDefault<scalar>("projectionCoeffFDEps", 1e-6);
    usePredictorReverseMomentumSens_ =
        thermalDict.getOrDefault<bool>
        (
            "usePredictorReverseMomentumSens",
            false
        );
    checkOneStepBetaMapSensitivity_ =
        thermalDict.getOrDefault<bool>
        (
            "checkOneStepBetaMapSensitivity",
            false
        );
    oneStepBetaMapFDEps_ =
        thermalDict.getOrDefault<scalar>("oneStepBetaMapFDEps", 1e-6);
    checkFixedPointMapAdjoint_ =
        thermalDict.getOrDefault<bool>
        (
            "checkFixedPointMapAdjoint",
            false
        );
    checkStateMapTranspose_ =
        thermalDict.getOrDefault<bool>
        (
            "checkStateMapTranspose",
            false
        );
    checkFullStateMapTranspose_ =
        thermalDict.getOrDefault<bool>
        (
            "checkFullStateMapTranspose",
            false
        );
    checkPressureMapTranspose_ =
        thermalDict.getOrDefault<bool>
        (
            "checkPressureMapTranspose",
            false
        );
    checkFixedPointTangentSensitivity_ =
        thermalDict.getOrDefault<bool>
        (
            "checkFixedPointTangentSensitivity",
            false
        );
    fixedPointMapAdjointIters_ =
        thermalDict.getOrDefault<label>("fixedPointMapAdjointIters", 5);
    fixedPointMapAdjointFDEps_ =
        thermalDict.getOrDefault<scalar>("fixedPointMapAdjointFDEps", 1e-6);
    fixedPointTangentMaxIters_ =
        thermalDict.getOrDefault<label>("fixedPointTangentMaxIters", 30);
    fixedPointTangentRelTol_ =
        thermalDict.getOrDefault<scalar>("fixedPointTangentRelTol", 1e-7);
    fixedPointTangentFDEps_ =
        thermalDict.getOrDefault<scalar>("fixedPointTangentFDEps", 1e-6);
    fixedPointTangentCells_ =
        thermalDict.getOrDefault<labelList>
        (
            "fixedPointTangentCells",
            labelList({707, 708, 709})
        );
    checkFixedPointTangentAgainstFD_ =
        thermalDict.getOrDefault<bool>
        (
            "checkFixedPointTangentAgainstFD",
            false
        );
    fixedPointValidationMaxIters_ =
        thermalDict.getOrDefault<label>
        (
            "fixedPointValidationMaxIters",
            200
        );
    fixedPointValidationRelTol_ =
        thermalDict.getOrDefault<scalar>
        (
            "fixedPointValidationRelTol",
            1e-9
        );
    fixedPointValidationBetaEps_ =
        thermalDict.getOrDefault<scalar>
        (
            "fixedPointValidationBetaEps",
            1e-6
        );
    kInterpolation_ =
        topOInterpolationFunction::New
        (
            mesh_,
            thermalDict.subDict("kInterpolation")
        );

    validateExactThermalCouplingOptions();

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
        couplingForm_ =
            thermalDict.getOrDefault<word>("couplingForm", "TaGradT");
        boundedConvection_ =
            thermalDict.getOrDefault<bool>("boundedConvection", true);
        checkProjectionTranspose_ =
            thermalDict.getOrDefault<bool>("checkProjectionTranspose", false);
        checkScalarThermalTranspose_ =
            thermalDict.getOrDefault<bool>
            (
                "checkScalarThermalTranspose",
                false
            );
        exactScalarThermalTranspose_ =
            thermalDict.getOrDefault<bool>
            (
                "exactScalarThermalTranspose",
                false
            );
        exactMomentumATC_ =
            thermalDict.getOrDefault<bool>("exactMomentumATC", false);
        exactConvectionTranspose_ =
            thermalDict.getOrDefault<bool>("exactConvectionTranspose", false);
        residualFDTopOSens_ =
            thermalDict.getOrDefault<bool>("residualFDTopOSens", false);
        replaceUpstreamMomentumSens_ =
            thermalDict.getOrDefault<bool>
            (
                "replaceUpstreamMomentumSens",
                false
            );
        residualFDEps_ =
            thermalDict.getOrDefault<scalar>("residualFDEps", 1e-6);
        projectionCoeffFDTopOSens_ =
            thermalDict.getOrDefault<bool>
            (
                "projectionCoeffFDTopOSens",
                false
            );
        addProjectionCoeffFDSens_ =
            thermalDict.getOrDefault<bool>("addProjectionCoeffFDSens", false);
        projectionCoeffFDEps_ =
            thermalDict.getOrDefault<scalar>("projectionCoeffFDEps", 1e-6);
        usePredictorReverseMomentumSens_ =
            thermalDict.getOrDefault<bool>
            (
                "usePredictorReverseMomentumSens",
                false
            );
        checkOneStepBetaMapSensitivity_ =
            thermalDict.getOrDefault<bool>
            (
                "checkOneStepBetaMapSensitivity",
                false
            );
        oneStepBetaMapFDEps_ =
            thermalDict.getOrDefault<scalar>("oneStepBetaMapFDEps", 1e-6);
        checkFixedPointMapAdjoint_ =
            thermalDict.getOrDefault<bool>
            (
                "checkFixedPointMapAdjoint",
                false
            );
        checkStateMapTranspose_ =
            thermalDict.getOrDefault<bool>
            (
                "checkStateMapTranspose",
                false
            );
        checkFullStateMapTranspose_ =
            thermalDict.getOrDefault<bool>
            (
                "checkFullStateMapTranspose",
                false
            );
        checkPressureMapTranspose_ =
            thermalDict.getOrDefault<bool>
            (
                "checkPressureMapTranspose",
                false
            );
        checkFixedPointTangentSensitivity_ =
            thermalDict.getOrDefault<bool>
            (
                "checkFixedPointTangentSensitivity",
                false
            );
        fixedPointMapAdjointIters_ =
            thermalDict.getOrDefault<label>("fixedPointMapAdjointIters", 5);
        fixedPointMapAdjointFDEps_ =
            thermalDict.getOrDefault<scalar>
            (
                "fixedPointMapAdjointFDEps",
                1e-6
            );
        fixedPointTangentMaxIters_ =
            thermalDict.getOrDefault<label>
            (
                "fixedPointTangentMaxIters",
                30
            );
        fixedPointTangentRelTol_ =
            thermalDict.getOrDefault<scalar>
            (
                "fixedPointTangentRelTol",
                1e-7
            );
        fixedPointTangentFDEps_ =
            thermalDict.getOrDefault<scalar>
            (
                "fixedPointTangentFDEps",
                1e-6
            );
        fixedPointTangentCells_ =
            thermalDict.getOrDefault<labelList>
            (
                "fixedPointTangentCells",
                labelList({707, 708, 709})
            );
        checkFixedPointTangentAgainstFD_ =
            thermalDict.getOrDefault<bool>
            (
                "checkFixedPointTangentAgainstFD",
                false
            );
        fixedPointValidationMaxIters_ =
            thermalDict.getOrDefault<label>
            (
                "fixedPointValidationMaxIters",
                200
            );
        fixedPointValidationRelTol_ =
            thermalDict.getOrDefault<scalar>
            (
                "fixedPointValidationRelTol",
                1e-9
            );
        fixedPointValidationBetaEps_ =
            thermalDict.getOrDefault<scalar>
            (
                "fixedPointValidationBetaEps",
                1e-6
            );
        // rebuild: otherwise RAMP continuation through a dictionary re-read
        // would be silently ignored
        kInterpolation_ =
            topOInterpolationFunction::New
            (
                mesh_,
                thermalDict.subDict("kInterpolation")
            );

        validateExactThermalCouplingOptions();

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
    const scalarField betaBefore(betaMult);

    // Brinkman momentum contribution (upstream machinery)
    adjointSimple::topOSensMultiplier(betaMult, designVariablesName, dt);

    const scalarField upstreamContribution(betaMult - betaBefore);

    if (residualFDTopOSens_ && replaceUpstreamMomentumSens_)
    {
        betaMult = betaBefore;
    }

    addResidualFDMomentumSensitivity
    (
        betaMult,
        upstreamContribution,
        designVariablesName,
        dt
    );

    addProjectionCoeffFDSensitivity(betaMult, designVariablesName, dt);

    if (checkFixedPointTangentSensitivity_)
    {
        const vectorField dummySeed(mesh_.nCells(), vector::zero);
        checkFixedPointTangentSensitivity(dummySeed, designVariablesName, dt);
    }

    if (checkFullStateMapTranspose_)
    {
        checkFullStateMapTranspose();
    }

    const bool needExactFluxDiagnostic =
        couplingForm_ == "exactFluxTranspose"
     && (
            checkProjectionTranspose_
         || checkPressureMapTranspose_
         || checkStateMapTranspose_
         || checkFullStateMapTranspose_
         || usePredictorReverseMomentumSens_
         || checkOneStepBetaMapSensitivity_
         || checkFixedPointMapAdjoint_
        );

    if (needExactFluxDiagnostic)
    {
        tmp<surfaceScalarField> tGPhi = thermalFluxSensitivity();
        tmp<volVectorField> tSource = projectedFluxMomentumSource(tGPhi());

        if (checkProjectionTranspose_)
        {
            checkFlowBlockTranspose(tSource());
        }
    }

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
