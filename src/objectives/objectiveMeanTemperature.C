/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "objectiveMeanTemperature.H"
#include "createZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace objectives
{

defineTypeNameAndDebug(objectiveMeanTemperature, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveMeanTemperature,
    dictionary
);


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField& objectiveMeanTemperature::TRef() const
{
    if (!mesh_.foundObject<volScalarField>("T"))
    {
        FatalErrorInFunction
            << "No registered temperature field 'T'. "
            << "objectiveMeanTemperature requires a thermal primal solver "
            << "(thermalSimple)."
            << exit(FatalError);
    }

    return mesh_.lookupObject<volScalarField>("T");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

objectiveMeanTemperature::objectiveMeanTemperature
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    zones_(mesh_.cellZones().indices(dict.get<wordRes>("zones"))),
    VTotal_(Zero)
{
    // Total zone volume (global)
    const scalarField& V = mesh_.V();
    for (const label zI : zones_)
    {
        const cellZone& zone = mesh_.cellZones()[zI];
        for (const label cellI : zone)
        {
            VTotal_ += V[cellI];
        }
    }
    reduce(VTotal_, sumOp<scalar>());

    if (VTotal_ < SMALL)
    {
        FatalErrorInFunction
            << "Zero total volume for zones " << dict.get<wordRes>("zones")
            << exit(FatalError);
    }

    // Allocate the volume source w.r.t. T
    dJdTPtr_.reset(createZeroFieldPtr<scalar>(mesh_, "dJdT", dimless));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

scalar objectiveMeanTemperature::J()
{
    J_ = Zero;

    const volScalarField& T = TRef();
    const scalarField& V = mesh_.V();

    for (const label zI : zones_)
    {
        const cellZone& zone = mesh_.cellZones()[zI];
        for (const label cellI : zone)
        {
            J_ += T[cellI]*V[cellI];
        }
    }
    reduce(J_, sumOp<scalar>());
    J_ /= VTotal_;

    return J_;
}


void objectiveMeanTemperature::update_dJdT()
{
    volScalarField& dJdT = dJdTPtr_();
    dJdT.primitiveFieldRef() = Zero;

    for (const label zI : zones_)
    {
        const cellZone& zone = mesh_.cellZones()[zI];
        for (const label cellI : zone)
        {
            dJdT[cellI] = 1.0/VTotal_;
        }
    }

    dJdT.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
