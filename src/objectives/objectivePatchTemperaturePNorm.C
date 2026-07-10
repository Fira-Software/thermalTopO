/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "objectivePatchTemperaturePNorm.H"
#include "createZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace objectives
{

defineTypeNameAndDebug(objectivePatchTemperaturePNorm, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectivePatchTemperaturePNorm,
    dictionary
);


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField& objectivePatchTemperaturePNorm::TRef() const
{
    if (!mesh_.foundObject<volScalarField>("T"))
    {
        FatalErrorInFunction
            << "No registered temperature field 'T'. "
            << "objectivePatchTemperaturePNorm requires a thermal primal "
            << "solver (thermalSimple)."
            << exit(FatalError);
    }

    return mesh_.lookupObject<volScalarField>("T");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

objectivePatchTemperaturePNorm::objectivePatchTemperaturePNorm
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    patches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        ).sortedToc()
    ),
    P_(dict.getOrDefault<scalar>("P", 12)),
    area_(Zero)
{
    if (patches_.empty())
    {
        FatalErrorInFunction
            << "No valid patches for " << dict.get<wordRes>("patches")
            << exit(FatalError);
    }

    for (const label patchI : patches_)
    {
        area_ += gSum(mesh_.boundary()[patchI].magSf());
    }

    // Allocate the boundary contribution w.r.t. T
    bdJdTPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

scalar objectivePatchTemperaturePNorm::J()
{
    const volScalarField& T = TRef();

    // conditioning scale: current global patch maximum
    scalar TScale = SMALL;
    for (const label patchI : patches_)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchI];
        if (!Tp.empty())
        {
            TScale = max(TScale, max(Tp));
        }
    }
    reduce(TScale, maxOp<scalar>());

    scalar S = Zero;
    for (const label patchI : patches_)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchI];
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
        S += sum(pow(Tp/TScale, P_)*magSf);
    }
    reduce(S, sumOp<scalar>());
    S /= area_;

    J_ = TScale*pow(S, scalar(1)/P_);

    return J_;
}


void objectivePatchTemperaturePNorm::update_boundarydJdT()
{
    const volScalarField& T = TRef();

    // Ensure J_ (and its internal scale) reflect the current T
    const scalar Jcur = J();
    const scalar TScale = Jcur > SMALL ? Jcur : SMALL;

    for (const label patchI : patches_)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchI];

        // bdJdT = J^{1-P} T^{P-1} / area, computed in scaled form
        // (T/J)^{P-1} / area  ==  J^{1-P} T^{P-1} / area
        bdJdTPtr_()[patchI] = pow(Tp/TScale, P_ - scalar(1))/area_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
