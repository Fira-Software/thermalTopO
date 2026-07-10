/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.

Class
    Foam::objectives::objectiveMeanTemperature

Description
    Volume-averaged temperature over a set of cellZones,

        J = sum_z int_z T dV / V_total

    The simplest thermally-coupled objective; primarily used for
    finite-difference verification of the adjoint gradients (V1 gate).

    Requires a primal solver providing a registered temperature field "T"
    (thermalSimple).

    Dictionary example (inside objectives sub-dict):
    \verbatim
    meanT
    {
        weight  1;
        type    meanTemperature;
        zones   (heatedZone);
    }
    \endverbatim

SourceFiles
    objectiveMeanTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_objectiveMeanTemperature_H
#define Foam_objectiveMeanTemperature_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace objectives
{

/*---------------------------------------------------------------------------*\
                   Class objectiveMeanTemperature Declaration
\*---------------------------------------------------------------------------*/

class objectiveMeanTemperature
:
    public objectiveIncompressible
{
    // Private data

        //- Cell zones defining the integration region
        labelList zones_;

        //- Total volume of the integration region
        scalar VTotal_;


    // Private Member Functions

        //- Return the registered temperature field
        const volScalarField& TRef() const;


public:

    //- Runtime type information
    TypeName("meanTemperature");


    // Constructors

        //- Construct from components
        objectiveMeanTemperature
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveMeanTemperature() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Volume source term for the adjoint temperature equation
        virtual void update_dJdT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
