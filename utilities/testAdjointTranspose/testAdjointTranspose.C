/*---------------------------------------------------------------------------*\
    thermalTopO utility: is fvm::div(-phi,Ta) the algebraic transpose of
    fvm::div(phi,T)?

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    Pure operator test. No optimisation, no primal or adjoint solve. Reads a
    frozen (phi, T, Ta) state, assembles the primal convection matrix A with the
    case's own scheme, and compares

        A^T Ta          (exact transpose, from A's ldu coefficients)
        div(-phi,Ta)    (what the adjoint solver actually assembles)

    The continuous adjoint of u.grad(T) is the reversed-flux operator, but the
    ALGEBRAIC transpose of an upwind matrix is a downwind operator, so the two
    need not agree. Where they disagree, Ta is wrong before the momentum coupling
    is even considered.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"
    runTime.setTime(timeDirs.last(), timeDirs.size() - 1);

    Info<< "Reading frozen state (phi, T, Ta) at time " << runTime.timeName()
        << nl << endl;

    volScalarField T(IOobject("T", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    volScalarField Ta(IOobject("Ta", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    surfaceScalarField phi(IOobject("phi", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::NO_WRITE), mesh);

    // Primal convection matrix, assembled exactly as thermalSimple does
    fvScalarMatrix A(fvm::div(phi, T));

    // Adjoint convection operator, assembled exactly as thermalAdjointSimple does
    fvScalarMatrix B(fvm::div(-phi, Ta));

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const label nCells = mesh.nCells();

    // y = A^T * Ta   (exact algebraic transpose from the ldu coefficients)
    scalarField yT(nCells, Zero);
    {
        const scalarField& d = A.diag();
        const scalarField& up = A.upper();
        const scalarField& lo = A.lower();
        forAll(d, c) { yT[c] = d[c]*Ta[c]; }
        forAll(own, f)
        {
            // A: (A x)_own += upper[f]*x[nei];  (A x)_nei += lower[f]*x[own]
            // A^T: swap the roles of upper and lower
            yT[own[f]] += lo[f]*Ta[nei[f]];
            yT[nei[f]] += up[f]*Ta[own[f]];
        }
    }

    // y = B * Ta   (what the adjoint actually applies)
    scalarField yB(nCells, Zero);
    {
        const scalarField& d = B.diag();
        const scalarField& up = B.upper();
        const scalarField& lo = B.lower();
        forAll(d, c) { yB[c] = d[c]*Ta[c]; }
        forAll(own, f)
        {
            yB[own[f]] += up[f]*Ta[nei[f]];
            yB[nei[f]] += lo[f]*Ta[own[f]];
        }
    }

    // ---------------------------------------------------------------------
    // (1) OFF-DIAGONAL test: the interior convection stencil, free of any BC.
    //     A^T has upper^T = lower(A), lower^T = upper(A). So the adjoint
    //     operator B is the algebraic transpose iff
    //         B.upper == A.lower   AND   B.lower == A.upper
    // ---------------------------------------------------------------------
    {
        const scalarField dU(B.upper() - A.lower());
        const scalarField dL(B.lower() - A.upper());
        const scalar sU = max(gMax(mag(A.lower())), SMALL);
        const scalar sL = max(gMax(mag(A.upper())), SMALL);

        Info<< "=== (1) OFF-DIAGONAL (interior stencil, no BC involved) ===" << nl
            << "  max|B.upper - A.lower| / scale = " << gMax(mag(dU))/sU << nl
            << "  max|B.lower - A.upper| / scale = " << gMax(mag(dL))/sL << nl
            << (max(gMax(mag(dU))/sU, gMax(mag(dL))/sL) < 1e-10
                ? "  -> off-diagonals MATCH: the interior stencil IS the transpose."
                : "  -> off-diagonals DIFFER: the interior stencil is NOT the transpose.")
            << nl << endl;
    }

    // ---------------------------------------------------------------------
    // (2) FULL operator, restricted to cells with NO boundary face. The
    //     diagonal legitimately differs on boundary cells (T has fixedValue at
    //     the inlet, Ta has fixedValue 0 -- that is the correct adjoint BC, not
    //     an error), so those cells must be excluded to test the operator.
    // ---------------------------------------------------------------------
    boolList nearBoundary(nCells, false);
    forAll(mesh.boundary(), p)
    {
        if (mesh.boundary()[p].size() && !mesh.boundary()[p].coupled())
        {
            for (const label c : mesh.boundary()[p].faceCells())
            {
                nearBoundary[c] = true;
            }
        }
    }

    scalar maxDiffInt = 0, maxRefInt = 0;
    label nInt = 0, worst = -1;
    forAll(yT, c)
    {
        if (nearBoundary[c]) continue;
        ++nInt;
        maxRefInt = max(maxRefInt, mag(yT[c]));
        if (mag(yT[c] - yB[c]) > maxDiffInt)
        {
            maxDiffInt = mag(yT[c] - yB[c]);
            worst = c;
        }
    }

    Info<< "=== (2) FULL A^T Ta vs div(-phi,Ta), INTERIOR cells only ("
        << nInt << " of " << nCells << ") ===" << nl
        << "  max|A^T Ta|        = " << maxRefInt << nl
        << "  max|A^T Ta - B Ta| = " << maxDiffInt << nl
        << "  relative           = "
        << (maxRefInt > SMALL ? maxDiffInt/maxRefInt : 0) << nl;
    if (worst >= 0)
    {
        Info<< "  worst interior cell " << worst << " at " << mesh.C()[worst]
            << ": A^T Ta = " << yT[worst] << ", B Ta = " << yB[worst] << endl;
    }
    Info<< nl;

    // continuity: the bounded scheme adds Sp(div(faceFlux)) with OPPOSITE sign
    // for A and B, so a non-solenoidal phi breaks the transpose by 2*div(phi).
    {
        volScalarField divPhi(fvc::div(phi));
        Info<< "=== continuity (bounded schemes add Sp(div(phi)) with opposite"
            << " signs) ===" << nl
            << "  max|div(phi)| = " << gMax(mag(divPhi.primitiveField())) << nl
            << endl;
    }

    // ---------------------------------------------------------------------
    // (3) MOMENTUM SOURCE. The exact derivative of the primal convective term
    //     contracted with Ta, with respect to the face flux, is
    //
    //         d/dphi_f [ Ta^T R_conv ] = T_f * (Ta_P - Ta_N)
    //
    //     with T_f the SAME frozen face value the primal upwind scheme uses.
    //     Mapping phi_f = Sf_f . U_f and U_f = w_f U_P + (1-w_f) U_N back to the
    //     cells gives the exact face-based adjoint momentum source:
    //
    //         S_P += w_f     * T_f * (Ta_P - Ta_N) * Sf_f / V_P
    //         S_N += (1-w_f) * T_f * (Ta_P - Ta_N) * Sf_f / V_N
    //
    //     thermalAdjointSimple currently inserts the CONTINUOUS form
    //     Ta*grad(T) instead. Compare them.
    // ---------------------------------------------------------------------
    {
        const surfaceScalarField& w = mesh.weights();
        const surfaceVectorField& Sf = mesh.Sf();
        const scalarField& V = mesh.V();

        vectorField Sexact(nCells, Zero);

        forAll(own, f)
        {
            const label P = own[f];
            const label N = nei[f];

            // frozen upwind face value, exactly as "Gauss upwind" picks it
            const scalar Tf = (phi[f] > 0) ? T[P] : T[N];
            const scalar dTa = Ta[P] - Ta[N];
            const vector base = Tf*dTa*Sf[f];

            Sexact[P] += w[f]*base/V[P];
            Sexact[N] += (1 - w[f])*base/V[N];
        }

        // what the solver actually adds
        const volVectorField Scont(Ta*fvc::grad(T));

        scalar mE = 0, mC = 0, mD = 0;
        label worstC = -1;
        forAll(Sexact, c)
        {
            if (nearBoundary[c]) continue;
            mE = max(mE, mag(Sexact[c]));
            mC = max(mC, mag(Scont[c]));
            if (mag(Sexact[c] - Scont[c]) > mD)
            {
                mD = mag(Sexact[c] - Scont[c]);
                worstC = c;
            }
        }

        Info<< "=== (3) ADJOINT MOMENTUM SOURCE (interior cells) ===" << nl
            << "  max|S exact face-based|      = " << mE << nl
            << "  max|S continuous Ta*grad(T)| = " << mC << nl
            << "  max|difference|              = " << mD << nl
            << "  relative                     = "
            << (mE > SMALL ? mD/mE : 0) << nl;
        if (worstC >= 0)
        {
            Info<< "  worst cell " << worstC << " at " << mesh.C()[worstC] << nl
                << "    exact      = " << Sexact[worstC] << nl
                << "    continuous = " << Scont[worstC] << endl;
        }
        Info<< nl;
    }

    Info<< "End" << nl << endl;
    return 0;
}
