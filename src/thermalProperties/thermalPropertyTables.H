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
    Foam::thermalPropertyTables

Description
    Temperature-dependent material properties for the single-domain
    (Brinkman-penalised) conjugate heat transfer formulation, shared by the
    thermalSimple primal and thermalAdjointSimple adjoint solvers so that the
    two can never disagree about the property state.

    The steady energy equation with variable properties is

        rho(T) cp(T) u.grad(T) = div(k(T) grad(T))

    Dividing through by a *constant* reference volumetric heat capacity
    (rho cp)_ref gives the form actually solved:

        C(T) u.grad(T) = div(D(T) grad(T))

        C(T) = rho(T) cp(T) / (rho cp)_ref      [-]
        D(T) = k(T) / (rho cp)_ref              [m2/s]

    The reference must be a constant, not the local field: only then is
    continuity of D grad(T) equivalent to continuity of the physical heat
    flux k grad(T) across the fluid/solid interface. Dividing by the local
    rho cp would break conjugate flux continuity. See docs/derivation.md
    section 6.

    With no tables present, C == 1 and D == DFluid/DSolid, recovering the
    constant-property formulation of section 1 exactly.

    Dictionary entries (subDict "thermal" of the primal/adjoint solver dict):
    \verbatim
    thermal
    {
        DFluid          1.0e-5;     // k_f/(rho cp)_ref; ignored if kFluidTable
        DSolid          1.0e-3;     // k_s/(rho cp)_ref; ignored if DSolidTable

        // optional temperature dependence. Units: rho [kg/m3], cp [J/kg/K],
        // k [W/m/K], DSolid [m2/s]. rhoTable and cpTable come as a pair.
        TRef            300;        // reference T for (rho cp)_ref
        rhoTable        table ((300 1000) (340 850));
        cpTable         table ((300 4000) (340 3600));
        kFluidTable     table ((300 40)   (340 44));
        DSolidTable     table ((300 1.2e-3) (340 0.8e-3));
        // rhoCpRef      4.0e6;     // optional explicit override [J/m3/K]
    }
    \endverbatim

    Note: rho(T) and mu(T) are *also* needed by the momentum equation, where
    they enter as nu(T) = mu(T)/rho(T) through the viscosityModels::
    temperatureTable transport model selected in constant/transportProperties.
    The rhoTable given there must match the one given here; they are read from
    different dictionaries because OpenFOAM constructs the transport model and
    the primal solver from different files.

SourceFiles
    thermalPropertyTables.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_thermalPropertyTables_H
#define Foam_thermalPropertyTables_H

#include "volFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class thermalPropertyTables Declaration
\*---------------------------------------------------------------------------*/

class thermalPropertyTables
{
    // Private Data

        //- Constant diffusivities [m2/s], used where no table is given
        scalar DFluid_;
        scalar DSolid_;

        //- Reference volumetric heat capacity (rho cp)_ref [J/m3/K]
        scalar rhoCpRef_;

        //- Reference temperature at which (rho cp)_ref is evaluated [K]
        scalar TRef_;

        //- Optional property tables
        autoPtr<Function1<scalar>> rhoPtr_;
        autoPtr<Function1<scalar>> cpPtr_;
        autoPtr<Function1<scalar>> kFluidPtr_;
        autoPtr<Function1<scalar>> DSolidPtr_;

        // Optional multiple solid materials, labelled by cellZone.
        // The design variable still chooses fluid vs SOLID; it does not choose
        // WHICH solid. Material labels are fixed for the whole optimisation, and
        // design cells that solidify take defaultGeneratedSolidMaterial.

        //- True if a "solidMaterials" sub-dictionary was supplied
        bool multiMaterial_;

        //- Material names, in dictionary order
        wordList materialNames_;

        //- Constant diffusivity per material [m2/s] (used where no table)
        scalarList materialConst_;

        //- Optional D_s(T) table per material (entry may be unset)
        PtrList<Function1<scalar>> materialTable_;

        //- Material index per cell
        labelList materialID_;

        //- Material assigned to design cells that become solid
        label defaultMaterial_;


    // Private Member Functions

        //- Tabulated field: scale*f(T), evaluated on cells and boundary faces
        tmp<volScalarField> evaluate
        (
            const Function1<scalar>& f,
            const volScalarField& T,
            const word& name,
            const scalar scale,
            const dimensionSet& dims
        ) const;

        //- Uniform field of the given value
        tmp<volScalarField> uniform
        (
            const volScalarField& T,
            const word& name,
            const scalar value,
            const dimensionSet& dims
        ) const;

        //- Solid diffusivity of one material at one temperature [m2/s]
        scalar materialDSolid(const label matID, const scalar T) const;


public:

    // Constructors

        //- Construct null; call read() before use
        thermalPropertyTables();


    // Member Functions

        //- Read the constants and any optional tables
        void read(const dictionary& thermalDict, const fvMesh& mesh);

        //- True if k_f is tabulated, i.e. D_f varies with T
        bool variableDFluid() const noexcept { return bool(kFluidPtr_); }

        //- True if rho cp is tabulated, i.e. C(T) is not unity
        bool variableRhoCp() const noexcept { return bool(rhoPtr_); }

        //- True if the solid diffusivity is tabulated
        bool variableDSolid() const noexcept;

        //- True if any thermal coefficient is temperature-dependent
        bool hasTemperatureDependentThermalCoefficients() const noexcept;

        //- True if several solid materials are labelled in the domain
        bool multiMaterial() const noexcept { return multiMaterial_; }

        //- Material names, in dictionary order
        const wordList& materialNames() const noexcept
        {
            return materialNames_;
        }

        //- Material index per cell
        const labelList& materialID() const noexcept { return materialID_; }

        //- Material index per cell, as a writable field (diagnostic/reporting)
        tmp<volScalarField> materialIDField
        (
            const volScalarField& T,
            const word& name
        ) const;

        //- Reference volumetric heat capacity [J/m3/K]
        scalar rhoCpRef() const noexcept { return rhoCpRef_; }

        //- Fluid diffusivity D_f(T) = k_f(T)/(rho cp)_ref [m2/s]
        tmp<volScalarField> DFluidField
        (
            const volScalarField& T,
            const word& name
        ) const;

        //- Solid diffusivity D_s(T) [m2/s]
        tmp<volScalarField> DSolidField
        (
            const volScalarField& T,
            const word& name
        ) const;

        //- Normalised volumetric heat capacity C(T) = rho cp/(rho cp)_ref [-]
        tmp<volScalarField> CField
        (
            const volScalarField& T,
            const word& name
        ) const;

        //- Local thermal diffusivity alpha(T) = k/(rho cp) = D_f/C [m2/s].
        //  Diagnostic only: the equation is assembled from C and D, not alpha.
        tmp<volScalarField> alphaField
        (
            const volScalarField& T,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
