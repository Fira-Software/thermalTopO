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

        //- Interpolation function for the conductivity RAMP
        autoPtr<topOInterpolationFunction> kInterpolation_;


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
