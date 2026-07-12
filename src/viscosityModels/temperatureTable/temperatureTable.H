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
    Foam::viscosityModels::temperatureTable

Description
    Incompressible transport model with a tabulated temperature-dependent
    kinematic viscosity

        nu(T) = mu(T)/rho(T)

    from Function1 tables for the dynamic viscosity mu [Pa s] and the density
    rho [kg/m3], or directly from a table of nu [m2/s]. The temperature field
    is looked up from the object registry by name (default "T"), so the model
    works with any solver that registers a temperature: here, thermalSimple.

    nu is recomputed from the current T by correct(), which the SIMPLE loop of
    adjointOptimisationFoam calls once per primal iteration.

    This is the same registry-lookup pattern as the upstream
    viscosityModels::Arrhenius model, but nu is rebuilt from the tables on
    every call rather than accumulated onto the previous value, so repeated
    correct() calls are idempotent.

    Selection in constant/transportProperties:
    \verbatim
    transportModel  temperatureTable;

    temperatureTableCoeffs
    {
        TRef        300;        // nu used before T is available
        muTable     table ((300 1.0e-3) (340 5.0e-4));   // [Pa s]
        rhoTable    table ((300 1000)   (340 850));      // [kg/m3]
        // nuTable   table ((300 1e-6) (340 5.88e-7));   // [m2/s], alternative
        // field     T;         // name of the temperature field
    }
    \endverbatim

    rhoTable must match the one given to the primal/adjoint solver's "thermal"
    subdictionary: the two dictionaries are read separately by OpenFOAM and are
    not cross-checked.

SourceFiles
    temperatureTable.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_viscosityModels_temperatureTable_H
#define Foam_viscosityModels_temperatureTable_H

#include "viscosityModel.H"
#include "volFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                      Class temperatureTable Declaration
\*---------------------------------------------------------------------------*/

class temperatureTable
:
    public viscosityModel
{
protected:

    // Protected Data

        //- Model coefficients
        dictionary coeffs_;

        //- Dynamic viscosity mu(T) [Pa s]
        autoPtr<Function1<scalar>> muPtr_;

        //- Density rho(T) [kg/m3]
        autoPtr<Function1<scalar>> rhoPtr_;

        //- Kinematic viscosity nu(T) [m2/s], as an alternative to mu and rho
        autoPtr<Function1<scalar>> nuPtr_;

        //- Name of the temperature field in the registry
        word TName_;

        //- Temperature used before the field is available [K]
        scalar TRef_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Kinematic viscosity field
        volScalarField nu_;


    // Protected Member Functions

        //- nu at a single temperature
        scalar nuValue(const scalar T) const;

        //- Rebuild nu_ from the current temperature field
        void updateNu();


public:

    //- Run-time type information
    TypeName("temperatureTable");


    // Constructors

        //- Construct from components
        temperatureTable
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~temperatureTable() = default;


    // Member Functions

        //- Return the kinematic viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the kinematic viscosity on a patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Update nu from the current temperature
        virtual void correct()
        {
            updateNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
