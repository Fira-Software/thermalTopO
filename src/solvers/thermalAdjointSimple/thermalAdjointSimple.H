/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.

Class
    Foam::thermalAdjointSimple

Description
    Adjoint SIMPLE solver with an adjoint temperature (Ta) equation for
    conjugate-heat-transfer topology optimisation. Companion of the
    thermalSimple primal solver; see docs/derivation.md for the formulation.

    Solves, before each adjoint momentum iteration,

        div(-C(T) phi, Ta) - laplacian(DEff(beta, T), Ta) = - sum_k w_k dJdT_k

    and injects the thermal coupling source couplingSign*(C(T) Ta grad(T))
    into the adjoint momentum equation through the addMomentumSource hook.

    C(T) = rho(T) cp(T)/(rho cp)_ref is the normalised volumetric heat
    capacity of the primal equation; it is unity for constant properties, and
    both expressions then reduce to the constant-property adjoint.

    Boundary conditions on Ta: the solver drives every fixedGradient Ta
    patch with

        gradient = ( - sum_k w_k bdJdT_k ) / DEff

    realising the objective flux BC on objective patches and the adiabatic
    adjoint wall (zero gradient) elsewhere. Inlets and outlets take
    fixedValue 0 in 0/Ta (the outlet value is the high-Peclet limit of the
    exact Robin condition).

    Contributes the conductivity-interpolation sensitivity term through the
    topOSensMultiplier hook using its own kInterpolation RAMP.

    Supports temperature-dependent fluid and solid properties (tables rhoTable,
    cpTable, kFluidTable, DSolidTable), evaluated at the converged primal T.

    Frozen-property linearisation (v1, documented and FD-gated): the property
    derivatives dC/dT, dD/dT and dnu/dT are omitted from the adjoint operator,
    i.e. C, DEff and nu are treated as frozen coefficient fields at the primal
    solution rather than as functions of the state. Also frozen: the turbulent
    thermal diffusivity (constant Prt) and the turbulence model itself. The
    resulting gradient error is measured, not assumed: see docs/derivation.md
    section 6 and cases/varprops.

    Dictionary: same "thermal" sub-dictionary as thermalSimple, plus
    optional "couplingSign" (default 1) and "thermalSensScale" (default 1).

SourceFiles
    thermalAdjointSimple.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_thermalAdjointSimple_H
#define Foam_thermalAdjointSimple_H

#include "adjointSimple.H"
#include "topOInterpolationFunction.H"
#include "thermalPropertyTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class thermalAdjointSimple Declaration
\*---------------------------------------------------------------------------*/

class thermalAdjointSimple
:
    public adjointSimple
{
private:

    // Private Member Functions

        //- No copy construct
        thermalAdjointSimple(const thermalAdjointSimple&) = delete;

        //- No copy assignment
        void operator=(const thermalAdjointSimple&) = delete;


protected:

    // Protected data

        //- Adjoint temperature field
        autoPtr<volScalarField> TaPtr_;

        //- Fluid and solid properties, evaluated at the frozen primal T
        thermalPropertyTables props_;

        //- Turbulent Prandtl number (constant)
        scalar Prt_;

        //- Sign of the Ta grad(T) coupling into adjoint momentum
        scalar couplingSign_;

        //- Scale (incl. sign) of the conductivity sensitivity term
        scalar thermalSensScale_;

        //- Form of the (ATC-T) coupling into adjoint momentum:
        //- "TaGradT" (default, continuous), "negTGradTa" (continuous),
        //- "exactFluxTranspose" (discrete), or "none". See addMomentumSource().
        word couplingForm_;

        //- Does the primal use a "bounded" convection scheme (chi = 1)?
        //- Only affects couplingForm exactFluxTranspose.
        bool boundedConvection_;

        //- Run local projection-transpose inner-product checks
        bool checkProjectionTranspose_;

        //- Run the scalar thermal adjoint transpose check
        bool checkScalarThermalTranspose_;

        //- Solve Ta with the exact algebraic transpose of the primal TEqn
        bool exactScalarThermalTranspose_;

        //- Add exact discrete momentum ATC from the primal flux dependence.
        //- Use with ATCModel cancel when comparing against upstream ATC.
        bool exactMomentumATC_;

        //- Correct fvm::div(-phi,Ua) to the exact transpose of
        //- fvm::div(phi,U) for the frozen transported-velocity operator.
        bool exactConvectionTranspose_;

        //- Expensive diagnostic: replace or compare the upstream Brinkman
        //- multiplier with a finite difference of the assembled residual
        bool residualFDTopOSens_;

        //- If true, replace the upstream contribution with residual FD
        bool replaceUpstreamMomentumSens_;

        //- Physical-beta perturbation used by residualFDTopOSens
        scalar residualFDEps_;

        //- Expensive diagnostic for beta dependence of the primal SIMPLE
        //- projection coefficients at frozen U,p
        bool projectionCoeffFDTopOSens_;

        //- If true, add the coefficient-FD term to betaMult
        bool addProjectionCoeffFDSens_;

        //- Physical-beta perturbation for projectionCoeffFDTopOSens
        scalar projectionCoeffFDEps_;

        //- Diagnostic only: write the exact predictor-reverse momentum seed
        //- for comparison with Ua. Default false; does not replace betaMult.
        bool usePredictorReverseMomentumSens_;

        //- Diagnostic only: finite-difference one frozen SIMPLE map with
        //- respect to physical beta and contract with the thermal map seed.
        bool checkOneStepBetaMapSensitivity_;

        //- Physical-beta perturbation for checkOneStepBetaMapSensitivity
        scalar oneStepBetaMapFDEps_;

        //- Diagnostic only: truncated fixed-point SIMPLE-map adjoint series
        bool checkFixedPointMapAdjoint_;

        //- Diagnostic only: finite-difference transpose check for the full
        //- (U,p) SIMPLE state map before attempting a fixed-point solve
        bool checkStateMapTranspose_;

        //- Diagnostic only: finite-difference transpose check for the full
        //- (U,p,phi) SIMPLE state map
        bool checkFullStateMapTranspose_;

        //- Diagnostic only: finite-difference transpose check for the
        //- pressure sub-map including constrainPressure and boundary rows
        bool checkPressureMapTranspose_;

        //- Diagnostic only: forward tangent fixed-point SIMPLE-map solve
        bool checkFixedPointTangentSensitivity_;

        //- Number of reverse-map steps for checkFixedPointMapAdjoint
        label fixedPointMapAdjointIters_;

        //- Physical-beta perturbation for checkFixedPointMapAdjoint
        scalar fixedPointMapAdjointFDEps_;

        //- Maximum matrix-free GMRES iterations for tangent sensitivity
        label fixedPointTangentMaxIters_;

        //- Relative tolerance for matrix-free tangent GMRES
        scalar fixedPointTangentRelTol_;

        //- State and beta perturbation for tangent sensitivity diagnostics
        scalar fixedPointTangentFDEps_;

        //- Physical beta cells sampled by the tangent sensitivity diagnostic
        labelList fixedPointTangentCells_;

        //- Diagnostic only: compare tangent GMRES against converged map FD
        bool checkFixedPointTangentAgainstFD_;

        //- Maximum iterations for converged map FD validation
        label fixedPointValidationMaxIters_;

        //- Relative state-change tolerance for map FD validation
        scalar fixedPointValidationRelTol_;

        //- Base physical beta perturbation for map FD validation
        scalar fixedPointValidationBetaEps_;

        //- Interpolation function for the conductivity RAMP
        autoPtr<topOInterpolationFunction> kInterpolation_;

        //- Integrated reverse seed for one frozen SIMPLE map
        struct SimpleMapSeed
        {
            vectorField barU;
            scalarField barp;
            scalarField barPhiInternal;
            List<scalarField> barPhiBoundary;

            SimpleMapSeed(const label nCells)
            :
                barU(nCells, vector::zero),
                barp(nCells, Zero),
                barPhiInternal(),
                barPhiBoundary()
            {}

            SimpleMapSeed(const fvMesh& mesh)
            :
                barU(mesh.nCells(), vector::zero),
                barp(mesh.nCells(), Zero),
                barPhiInternal(mesh.nInternalFaces(), Zero),
                barPhiBoundary(mesh.boundary().size())
            {
                forAll(barPhiBoundary, patchi)
                {
                    barPhiBoundary[patchi].setSize
                    (
                        mesh.boundary()[patchi].size()
                    );
                    barPhiBoundary[patchi] = Zero;
                }
            }
        };

        //- Reverse seed returned by the SIMPLE predictor map
        struct PredictorReverseSeed
        {
            vectorField barUold;
            scalarField barPhiInternal;
            List<scalarField> barPhiBoundary;

            explicit PredictorReverseSeed(const fvMesh& mesh)
            :
                barUold(mesh.nCells(), vector::zero),
                barPhiInternal(mesh.nInternalFaces(), Zero),
                barPhiBoundary(mesh.boundary().size())
            {
                forAll(barPhiBoundary, patchi)
                {
                    barPhiBoundary[patchi].setSize
                    (
                        mesh.boundary()[patchi].size()
                    );
                    barPhiBoundary[patchi] = Zero;
                }
            }
        };

        //- Reverse seed for the final relaxed momentum diagonal's phi path
        struct FinalDiagonalPhiSeed
        {
            scalarField internal;
            List<scalarField> boundary;

            explicit FinalDiagonalPhiSeed(const fvMesh& mesh)
            :
                internal(mesh.nInternalFaces(), Zero),
                boundary(mesh.boundary().size())
            {
                forAll(boundary, patchi)
                {
                    boundary[patchi].setSize(mesh.boundary()[patchi].size());
                    boundary[patchi] = Zero;
                }
            }
        };

        //- Integrated seed on a volVectorField's explicit value space
        struct VolVectorValueSeed
        {
            vectorField internal;
            List<vectorField> boundary;

            explicit VolVectorValueSeed(const fvMesh& mesh)
            :
                internal(mesh.nCells(), vector::zero),
                boundary(mesh.boundary().size())
            {
                forAll(boundary, patchi)
                {
                    boundary[patchi].setSize
                    (
                        mesh.boundary()[patchi].size()
                    );
                    boundary[patchi] = vector::zero;
                }
            }
        };

        //- Reverse split for constrainHbyA's explicit boundary map
        struct ConstrainHbyAReverseSeed
        {
            vectorField barYInternal;
            List<vectorField> barYBoundary;
            List<vectorField> barUBoundary;

            explicit ConstrainHbyAReverseSeed(const fvMesh& mesh)
            :
                barYInternal(mesh.nCells(), vector::zero),
                barYBoundary(mesh.boundary().size()),
                barUBoundary(mesh.boundary().size())
            {
                forAll(barYBoundary, patchi)
                {
                    const label nFaces = mesh.boundary()[patchi].size();
                    barYBoundary[patchi].setSize(nFaces);
                    barUBoundary[patchi].setSize(nFaces);
                    barYBoundary[patchi] = vector::zero;
                    barUBoundary[patchi] = vector::zero;
                }
            }
        };

        //- Diagnostic coefficient seeds in the SIMPLE-map reverse
        struct SmoothPhiCoefficientTrace
        {
            scalarField barDfinal;
            scalarField barDFromY;
            scalarField barDDownstream;
            scalarField barDCombined;
            scalarField barUEqnH1Coeff;

            scalarField barrAUFromQ;
            scalarField barrAtUFromQ;
            scalarField barrAtUFromFinalCorrection;
            scalarField barrAtUFromPressureStage;

            scalarField barq;
            scalarField barqFromH1;
            scalarField barqFromF2InternalFaces;
            scalarField barqFromF2BoundaryCurrent;
            scalarField barQFaceInternal;
            List<scalarField> barQFaceBoundary;

            scalarField barrAtUBeforeQ;
            scalarField barrAtUAfterQ;
            scalarField barrAUBeforeQ;
            scalarField barrAUAfterQ;
            scalarField barrAUBeforeRAtUReverse;
            scalarField barrAUAfterRAtUReverse;
            scalarField barH1Expected;
            scalarField barDExpected;
            scalarField rAUUsed;
            scalarField rAtUUsed;

            vectorField barH0Internal;
            List<vectorField> barH0Boundary;
            vectorField barH0DirectInternal;
            List<vectorField> barH0DirectBoundary;
            vectorField barH0FromFluxInternal;
            List<vectorField> barH0FromFluxBoundary;
            vectorField barH0TotalInternal;
            List<vectorField> barH0TotalBoundary;

            vectorField barYInternal;
            List<vectorField> barYBoundary;
            List<vectorField> barUBoundary;

            vectorField barUFromBoundary;

            scalarField barpNew;
            scalarField barpSolveFromRelaxOnly;
            scalarField barpSolveFromPhiFlux;
            scalarField barpSolveTotal;
            scalarField barpSolve;

            scalarField barF2Internal;
            List<scalarField> barF2Boundary;

            scalarField barF1Internal;
            List<scalarField> barF1Boundary;

            scalarField barF0Internal;
            List<scalarField> barF0Boundary;

            explicit SmoothPhiCoefficientTrace(const fvMesh& mesh)
            :
                barDfinal(mesh.nCells(), Zero),
                barDFromY(mesh.nCells(), Zero),
                barDDownstream(mesh.nCells(), Zero),
                barDCombined(mesh.nCells(), Zero),
                barUEqnH1Coeff(mesh.nCells(), Zero),
                barrAUFromQ(mesh.nCells(), Zero),
                barrAtUFromQ(mesh.nCells(), Zero),
                barrAtUFromFinalCorrection(mesh.nCells(), Zero),
                barrAtUFromPressureStage(mesh.nCells(), Zero),
                barq(mesh.nCells(), Zero),
                barqFromH1(mesh.nCells(), Zero),
                barqFromF2InternalFaces(mesh.nCells(), Zero),
                barqFromF2BoundaryCurrent(mesh.nCells(), Zero),
                barQFaceInternal(mesh.nInternalFaces(), Zero),
                barQFaceBoundary(mesh.boundary().size()),
                barrAtUBeforeQ(mesh.nCells(), Zero),
                barrAtUAfterQ(mesh.nCells(), Zero),
                barrAUBeforeQ(mesh.nCells(), Zero),
                barrAUAfterQ(mesh.nCells(), Zero),
                barrAUBeforeRAtUReverse(mesh.nCells(), Zero),
                barrAUAfterRAtUReverse(mesh.nCells(), Zero),
                barH1Expected(mesh.nCells(), Zero),
                barDExpected(mesh.nCells(), Zero),
                rAUUsed(mesh.nCells(), Zero),
                rAtUUsed(mesh.nCells(), Zero),
                barH0Internal(mesh.nCells(), vector::zero),
                barH0Boundary(mesh.boundary().size()),
                barH0DirectInternal(mesh.nCells(), vector::zero),
                barH0DirectBoundary(mesh.boundary().size()),
                barH0FromFluxInternal(mesh.nCells(), vector::zero),
                barH0FromFluxBoundary(mesh.boundary().size()),
                barH0TotalInternal(mesh.nCells(), vector::zero),
                barH0TotalBoundary(mesh.boundary().size()),
                barYInternal(mesh.nCells(), vector::zero),
                barYBoundary(mesh.boundary().size()),
                barUBoundary(mesh.boundary().size()),
                barUFromBoundary(mesh.nCells(), vector::zero),
                barpNew(mesh.nCells(), Zero),
                barpSolveFromRelaxOnly(mesh.nCells(), Zero),
                barpSolveFromPhiFlux(mesh.nCells(), Zero),
                barpSolveTotal(mesh.nCells(), Zero),
                barpSolve(mesh.nCells(), Zero),
                barF2Internal(mesh.nInternalFaces(), Zero),
                barF2Boundary(mesh.boundary().size()),
                barF1Internal(mesh.nInternalFaces(), Zero),
                barF1Boundary(mesh.boundary().size()),
                barF0Internal(mesh.nInternalFaces(), Zero),
                barF0Boundary(mesh.boundary().size())
            {
                forAll(barH0Boundary, patchi)
                {
                    const label nFaces = mesh.boundary()[patchi].size();
                    barH0Boundary[patchi].setSize(nFaces, vector::zero);
                    barH0DirectBoundary[patchi].setSize(nFaces, vector::zero);
                    barH0FromFluxBoundary[patchi].setSize(nFaces, vector::zero);
                    barH0TotalBoundary[patchi].setSize(nFaces, vector::zero);
                    barYBoundary[patchi].setSize(nFaces, vector::zero);
                    barUBoundary[patchi].setSize(nFaces, vector::zero);
                    barQFaceBoundary[patchi].setSize(nFaces, Zero);
                    barF2Boundary[patchi].setSize(nFaces, Zero);
                    barF1Boundary[patchi].setSize(nFaces, Zero);
                    barF0Boundary[patchi].setSize(nFaces, Zero);
                }
            }
        };

        //- Cell fields returned by one frozen SIMPLE map
        struct SimpleMapState
        {
            vectorField U;
            scalarField p;
            scalarField phiInternal;
            List<scalarField> phiBoundary;

            SimpleMapState()
            :
                U(),
                p(),
                phiInternal(),
                phiBoundary()
            {}

            SimpleMapState(const label nCells)
            :
                U(nCells, vector::zero),
                p(nCells, Zero),
                phiInternal(),
                phiBoundary()
            {}

            SimpleMapState(const fvMesh& mesh)
            :
                U(mesh.nCells(), vector::zero),
                p(mesh.nCells(), Zero),
                phiInternal(mesh.nInternalFaces(), Zero),
                phiBoundary(mesh.boundary().size())
            {
                forAll(phiBoundary, patchi)
                {
                    phiBoundary[patchi].setSize
                    (
                        mesh.boundary()[patchi].size()
                    );
                    phiBoundary[patchi] = Zero;
                }
            }
        };


    // Protected Member Functions

        //- Return the primal temperature field (registered by thermalSimple)
        const volScalarField& TRef() const;

        //- Conductivity indicator field Ik(beta) in [0,1]
        tmp<volScalarField> kIndicator() const;

        //- Effective thermal diffusivity DEff(beta)
        tmp<volScalarField> DEff() const;

        //- Drive fixedGradient Ta boundary conditions from the objectives
        //- and the Robin outlet relation
        void updateTaBCs(const volScalarField& DEff);

        //- Assemble and solve the adjoint temperature equation
        virtual void solveTaEqn();

        //- Optional check that Ta transposes the frozen primal thermal block
        void checkScalarThermalTranspose
        (
            const volScalarField& DEffField,
            const volScalarField& source
        );

        //- Guard exact discrete diagnostics against the legacy scalar adjoint
        void validateExactThermalCouplingOptions() const;

        //- Exact thermal sensitivity with respect to the pressure-corrected
        //- primal face flux
        tmp<surfaceScalarField> thermalFluxSensitivity() const;

        //- Apply the transpose of the SIMPLE pressure projection to gPhi and
        //- return the corresponding corrected-velocity cell source
        tmp<volVectorField> projectedFluxMomentumSource
        (
            const surfaceScalarField& gPhi,
            const word& namePrefix = "ATCT"
        );

        //- Exact sensitivity of primal momentum convection to phi
        tmp<surfaceScalarField> momentumFluxSensitivity() const;

        //- Exact discrete ATC source from momentum flux dependence
        tmp<volVectorField> exactMomentumATCSource();

        //- Exact discrete transpose correction for frozen momentum convection
        tmp<volVectorField> exactConvectionTransposeSource();

        //- Optional local algebra checks for the projection and velocity maps
        void checkProjectedFluxTranspose
        (
            const surfaceScalarField& gPhi,
            const surfaceScalarField& gF,
            const volVectorField& source,
            const volScalarField& rAU,
            const volScalarField& rAtU,
            const fvScalarMatrix& pEqn
        );

        //- Rebuild the exact SIMPLE predictor map Y = rAU*UEqn.H()
        tmp<volVectorField> computePredictorY
        (
            const volVectorField& Uprobe,
            const word& namePrefix,
            const bool doRelax = true,
            const bool doFvOptionsConstrain = true
        );

        //- Rebuild predictor pieces: total, ldu-only, source-only
        void computePredictorYPieces
        (
            const volVectorField& Uprobe,
            const word& namePrefix,
            volVectorField& Ytotal,
            volVectorField& Yldu,
            volVectorField& Ysource,
            const bool doRelax = true,
            const bool doFvOptionsConstrain = true,
            scalarField* rAUI = nullptr,
            volVectorField* YsourceRaw = nullptr,
            volVectorField* YsourceRelax = nullptr,
            volVectorField* YsourceConstrain = nullptr
        );

        //- Finite-difference truth-model check for the predictor-map transpose
        void checkPredictorMapTranspose
        (
            const vectorField& barYInt,
            const word& checkName,
            const bool doRelax = true,
            const bool doFvOptionsConstrain = true
        );

        //- Analytic reverse of Y = rAU*UEqn.H() for the predictor-map check
        tmp<vectorField> reversePredictorY
        (
            const vectorField& barYInt,
            const bool doRelax = true,
            const bool doFvOptionsConstrain = true
        );

        //- Full-state analytic reverse of Y = rAU*UEqn.H()
        PredictorReverseSeed reversePredictorYFullState
        (
            const vectorField& barYInt,
            const bool doRelax = true,
            const bool doFvOptionsConstrain = true
        );

        //- Full-state analytic reverse of Y with extra UEqn coefficient seeds
        PredictorReverseSeed reversePredictorYFullState
        (
            const vectorField& barYInt,
            const scalarField& extraBarDfinal,
            const scalarField& extraBarUEqnH1Coeff,
            const bool doRelax = true,
            const bool doFvOptionsConstrain = true
        );

        //- Full-state analytic reverse of Y with explicit boundary value seeds
        PredictorReverseSeed reversePredictorYFullState
        (
            const vectorField& barYInt,
            const List<vectorField>& barYBoundary,
            const scalarField& extraBarDfinal,
            const scalarField& extraBarUEqnH1Coeff,
            const bool doRelax = true,
            const bool doFvOptionsConstrain = true
        );

        //- Reverse extra downstream seeds on the final relaxed UEqn diagonal
        FinalDiagonalPhiSeed reverseExtraFinalDiagonalToPhi
        (
            const scalarField& extraBarDfinal,
            const scalar alphaU
        ) const;

        //- Reverse the fixed-pressure pEqn.flux coefficient dependence on rAtU
        scalarField reversePressureRAtUExact
        (
            const volScalarField& rAtUBase,
            const volScalarField& pSolveBase,
            const surfaceScalarField& barF2
        );

        //- Reverse explicit Stokes source-like term in the SIMPLE predictor
        vectorField reverseExplicitStokesViscousSource
        (
            const vectorField& barSourceIntegrated
        );

        //- One state reverse pass of the frozen SIMPLE map
        tmp<surfaceScalarField> reverseAdjustPhiExact
        (
            const surfaceScalarField& phiBeforeAdjust,
            const surfaceScalarField& barPhiAfterAdjust,
            const volVectorField& U
        );

        //- Internal implementation optionally recording coefficient seeds
        SimpleMapSeed reverseOneSimpleMapSeedImpl
        (
            const SimpleMapSeed& seedNew,
            const volScalarField& rAtU,
            SmoothPhiCoefficientTrace* trace
        );

        //- One state reverse pass of the frozen SIMPLE map
        SimpleMapSeed reverseOneSimpleMapSeed
        (
            const SimpleMapSeed& seedNew,
            const volScalarField& rAtU
        );

        //- One velocity-only reverse pass of the frozen SIMPLE map
        tmp<vectorField> reverseSimpleMapVelocitySeed
        (
            const vectorField& barUnewInt,
            const volScalarField& rAtU
        );

        //- Deprecated hybrid residual/SIMPLE-map diagnostic, not a pass/fail gate
        void checkFlowBlockTranspose(const volVectorField& source);

        //- Integrated left-hand-side primal momentum residual at frozen U,p
        tmp<vectorField> primalMomentumResidualLHS();

        //- Rebuild the pressure-projected primal flux at frozen U,p
        tmp<surfaceScalarField> primalProjectedFluxAtFrozenState
        (
            const word& namePrefix
        );

        //- Rebuild one frozen SIMPLE map and return the corrected U update
        tmp<volVectorField> primalSimpleMapUAtFrozenState
        (
            const word& namePrefix
        );

        //- Rebuild one frozen SIMPLE map and return corrected U and p
        SimpleMapState primalSimpleMapStateAtFrozenState
        (
            const word& namePrefix
        );

        //- Rebuild the pressure sub-map pnew = P(F2,Uold,pold)
        tmp<volScalarField> computePressureMap
        (
            const surfaceScalarField& F2,
            const volVectorField& Uold,
            const volScalarField& pold,
            const volScalarField& rAtU,
            const word& namePrefix
        );

        //- Diagnostic physical-beta derivative of one frozen SIMPLE map
        void checkOneStepBetaMapSensitivity
        (
            const vectorField& barUnewInt,
            const vectorField& UbarExactInt,
            const word& designVariablesName,
            const scalar dt
        );

        //- Diagnostic truncated fixed-point map-adjoint beta sensitivity
        void checkFixedPointMapAdjointSensitivity
        (
            const vectorField& barUnewInt,
            const volScalarField& rAtU,
            const word& designVariablesName,
            const scalar dt
        );

        //- Diagnostic forward tangent fixed-point beta sensitivity
        void checkFixedPointTangentSensitivity
        (
            const vectorField& barUnewInt,
            const word& designVariablesName,
            const scalar dt
        );

        //- Diagnostic transpose identity for one frozen (U,p) SIMPLE map
        void checkStateMapTranspose
        (
            const vectorField& barUnewInt,
            const volScalarField& rAtU
        );

        //- Diagnostic transpose identity for one frozen (U,p,phi) SIMPLE map
        void checkFullStateMapTranspose();

        //- Diagnostic transpose identity for the pressure sub-map
        void checkPressureMapTranspose
        (
            const vectorField& barUnewInt,
            const volScalarField& rAtU
        );

        //- Diagnostic residual-FD replacement for the Brinkman multiplier
        void addResidualFDMomentumSensitivity
        (
            scalarField& betaMult,
            const scalarField& upstreamContribution,
            const word& designVariablesName,
            const scalar dt
        );

        //- Diagnostic direct beta sensitivity of SIMPLE flux coefficients
        void addProjectionCoeffFDSensitivity
        (
            scalarField& betaMult,
            const word& designVariablesName,
            const scalar dt
        );

        //- Add the thermal coupling source to the adjoint momentum equation
        virtual void addMomentumSource(fvVectorMatrix& matrix);


public:

    // Static Data Members

        //- Run-time type information
        TypeName("thermalAdjointSimple");


    // Constructors

        //- Construct from mesh and dictionary
        thermalAdjointSimple
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName,
            const word& solverName
        );


    //- Destructor
    virtual ~thermalAdjointSimple() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void mainIter();


        // Sensitivity contributions

            //- Add the conductivity-interpolation term to the topology
            //- sensitivities
            virtual void topOSensMultiplier
            (
                scalarField& betaMult,
                const word& designVariablesName,
                const scalar dt
            );


        // Access

            //- Const access to the adjoint temperature
            const volScalarField& Ta() const { return TaPtr_(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
