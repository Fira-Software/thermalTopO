/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.

Class
    Foam::objectives::objectivePatchTemperaturePNorm

Description
    Differentiable surrogate for the peak temperature on a set of patches:

        J = ( int_Gamma T^P dGamma / |Gamma| )^{1/P}

    J -> max(T) as P -> inf; typical P = 8..20 with continuation. Internally
    scaled by the current patch maximum so arbitrary P remains
    well-conditioned.

    The adjoint temperature equation receives the boundary contribution

        bdJdT = J^{1-P} T^{P-1} / |Gamma|     (per unit area)

    used by the adjoint solver as a flux-type boundary condition on Gamma
    (see docs/derivation.md section 3.1).

    Dictionary example:
    \verbatim
    peakT
    {
        weight  1;
        type    patchTemperaturePNorm;
        patches (plasmaFacing);
        P       12;
    }
    \endverbatim

SourceFiles
    objectivePatchTemperaturePNorm.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_objectivePatchTemperaturePNorm_H
#define Foam_objectivePatchTemperaturePNorm_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace objectives
{

/*---------------------------------------------------------------------------*\
              Class objectivePatchTemperaturePNorm Declaration
\*---------------------------------------------------------------------------*/

class objectivePatchTemperaturePNorm
:
    public objectiveIncompressible
{
    // Private data

        //- Patches defining Gamma
        labelList patches_;

        //- p-norm exponent
        scalar P_;

        //- Total patch area (global)
        scalar area_;


    // Private Member Functions

        //- Return the registered temperature field
        const volScalarField& TRef() const;


public:

    //- Runtime type information
    TypeName("patchTemperaturePNorm");


    // Constructors

        //- Construct from components
        objectivePatchTemperaturePNorm
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectivePatchTemperaturePNorm() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Boundary contribution to the adjoint temperature BC on Gamma
        virtual void update_boundarydJdT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
