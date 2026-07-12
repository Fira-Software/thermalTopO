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

        C(T) div(phi, T) - laplacian(DEff(beta, T), T) = 0

    where the effective diffusivity interpolates between fluid and solid
    conductivity over the topology-optimisation design field beta, looked up
    from the mesh registry ("topOVars", as populated by topODesignVariables),
    plus a turbulent contribution nut/Prt masked out of solidified cells:

        DEff(beta, T) = DFluid(T) + (DSolid(T) - DFluid(T))*Ik(beta)
                      + (1 - Ik(beta))*nut/Prt

    All quantities are kinematic, scaled by a *constant* reference volumetric
    heat capacity (rho cp)_ref: DFluid = k_f/(rho cp)_ref,
    DSolid = k_s/(rho cp)_ref, and

        C(T) = rho(T) cp(T)/(rho cp)_ref

    is the normalised volumetric heat capacity multiplying convection. With
    constant properties C == 1 and the equation collapses to the form of
    docs/derivation.md section 1; with tabulated rho, cp, k it is the variable
    -property form of section 6. The reference must be constant for continuity
    of DEff grad(T) to remain equivalent to continuity of the physical flux
    k grad(T) at the fluid/solid interface.

    Temperature-dependent properties are handled by thermalPropertyTables;
    mu(T) enters the momentum equation separately, through the
    viscosityModels::temperatureTable transport model.

    Dictionary entries (in the primal solver dict, subDict "thermal"):
    \verbatim
    thermal
    {
        DFluid          1.72e-7;   // k_f/(rho cp)_ref
        DSolid          3.55e-5;   // k_s/(rho cp)_ref
        // optional temperature dependence; a table overrides its constant:
        // TRef         300;
        // rhoTable     table ((300 1000)    (340 850));      // [kg/m3]
        // cpTable      table ((300 4000)    (340 3600));     // [J/kg/K]
        // kFluidTable  table ((300 40)      (340 44));       // [W/m/K]
        // DSolidTable  table ((373 3.9e-5)  (673 3.3e-5));   // [m2/s]
        Prt             0.85;      // constant turbulent Prandtl number
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
#include "thermalPropertyTables.H"

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

        //- Fluid and solid properties, constant or tabulated in T
        thermalPropertyTables props_;

        //- Turbulent Prandtl number (constant; see docs/derivation.md 6.4)
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
