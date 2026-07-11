/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details, <https://www.gnu.org/licenses/>.

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.

Class
    Foam::thermalSimple

Description
    Primal SIMPLE solver for topology optimisation with conjugate heat
    transfer in a single-domain (Brinkman-penalised) formulation.

    Extends the adjointOptimisationFoam 'simple' primal solver with a
    temperature equation

        div(phi, T) - laplacian(DEff(beta), T) = 0

    where the effective diffusivity interpolates between fluid and solid
    conductivity over the topology-optimisation design field beta, looked up
    from the mesh registry ("topOVars", as populated by topODesignVariables),
    plus a turbulent contribution nut/Prt masked out of solidified cells:

        DEff(beta) = DFluid + (DSolid - DFluid)*Ik(beta)
                   + (1 - Ik(beta))*nut/Prt

    All quantities are kinematic: DFluid = k_f/(rho_f cp_f),
    DSolid = k_s/(rho_f cp_f). See docs/derivation.md section 1 for why the
    single (rho cp)_f scaling is exact at steady state.

    Dictionary entries (in the primal solver dict, subDict "thermal"):
    \verbatim
    thermal
    {
        DFluid          1.72e-7;   // k_f/(rho_f cp_f)
        DSolid          3.55e-5;   // k_s/(rho_f cp_f)
        // optional temperature dependence, overrides DSolid:
        // DSolidTable  table ((373 3.9e-5) (673 3.3e-5));
        Prt             0.85;
        kInterpolation             // RAMP for conductivity interpolation
        {
            function    BorrvallPetersson;
            b           20;
        }
    }
    \endverbatim

SourceFiles
    thermalSimple.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_thermalSimple_H
#define Foam_thermalSimple_H

#include "simple.H"
#include "topOInterpolationFunction.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class thermalSimple Declaration
\*---------------------------------------------------------------------------*/

class thermalSimple
:
    public simple
{
private:

    // Private Member Functions

        //- No copy construct
        thermalSimple(const thermalSimple&) = delete;

        //- No copy assignment
        void operator=(const thermalSimple&) = delete;


protected:

    // Protected data

        //- Temperature field
        autoPtr<volScalarField> TPtr_;

        //- Molecular diffusivities, fluid and solid [m2/s]
        scalar DFluid_;
        scalar DSolid_;

        //- Optional temperature-dependent solid diffusivity D_s(T) [m2/s];
        //- overrides DSolid when present (entry "DSolidTable")
        autoPtr<Function1<scalar>> DSolidTablePtr_;

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Interpolation function for the conductivity RAMP
        autoPtr<topOInterpolationFunction> kInterpolation_;


    // Protected Member Functions

        //- Allocate the temperature field
        void allocateThermalFields();

        //- Conductivity indicator field Ik(beta) in [0,1].
        //  Zero everywhere if no topOVars in the registry (plain CHT run).
        tmp<volScalarField> kIndicator() const;

        //- Effective thermal diffusivity DEff(beta)
        tmp<volScalarField> DEff() const;

        //- Assemble and solve the temperature equation
        virtual void solveTEqn();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("thermalSimple");


    // Constructors

        //- Construct from mesh and dictionary
        thermalSimple
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& solverName
        );


    //- Destructor
    virtual ~thermalSimple() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void mainIter();


        // Access

            //- Const access to the temperature field
            const volScalarField& T() const { return TPtr_(); }

            //- Non-const access to the temperature field
            volScalarField& T() { return TPtr_.ref(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
