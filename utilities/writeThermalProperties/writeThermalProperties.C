/*---------------------------------------------------------------------------*\
    thermalTopO utility: write the thermal property fields the solver assembles.

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    Reads the "thermal" sub-dictionary of the primal solver from
    system/optimisationDict, plus a temperature field, and writes the property
    fields that thermalSimple/thermalAdjointSimple would assemble from them:

        DSolid      solid diffusivity  [m2/s]  (per-material, per-cell, D_s(T))
        DFluid      fluid diffusivity  [m2/s]
        materialID  solid material index per cell (multi-material runs)

    Reporting and verification only: it solves nothing.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "IOdictionary.H"
#include "thermalPropertyTables.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"
    runTime.setTime(timeDirs.last(), timeDirs.size() - 1);

    IOdictionary optDict
    (
        IOobject
        (
            "optimisationDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& primalSolvers = optDict.subDict("primalSolvers");
    const dictionary& op1 = primalSolvers.subDict(primalSolvers.toc()[0]);
    const dictionary& thermalDict = op1.subDict("thermal");

    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );

    thermalPropertyTables props;
    props.read(thermalDict, mesh);

    volScalarField DSolid(props.DSolidField(T, "DSolid"));
    volScalarField DFluid(props.DFluidField(T, "DFluid"));

    DSolid.write();
    DFluid.write();

    Info<< nl << "Wrote DSolid, DFluid at time " << runTime.timeName() << endl;

    if (props.multiMaterial())
    {
        volScalarField matID(props.materialIDField(T, "materialID"));
        matID.write();

        Info<< "Solid materials: " << props.materialNames() << nl
            << "Wrote materialID" << endl;
    }

    Info<< nl << "End" << nl << endl;
    return 0;
}
