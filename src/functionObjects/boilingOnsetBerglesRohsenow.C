/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "boilingOnsetBerglesRohsenow.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(boilingOnsetBerglesRohsenow, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        boilingOnsetBerglesRohsenow,
        dictionary
    );
}
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::functionObjects::boilingOnsetBerglesRohsenow::qWall
(
    const label patchi
) const
{
    const auto& turb =
        mesh_.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    const volScalarField& he = turb.transport().he();

    return
        turb.alphaEff(patchi)
       *he.boundaryField()[patchi].snGrad();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::functionObjects::boilingOnsetBerglesRohsenow::
boilingOnsetBerglesRohsenow
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    patchSet_(),
    pAbs_(Zero),
    Tsat_(Zero),
    writeFields_(false)
{
    read(dict);

    writeHeader(file(), "Bergles-Rohsenow onset-of-nucleate-boiling check");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "minMargin[W/m2]");
    writeTabbed(file(), "nFacesONB");
    writeTabbed(file(), "maxTw[K]");
    file() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

bool Foam::functionObjects::boilingOnsetBerglesRohsenow::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    patchSet_ =
        mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches"));
    pAbs_ = dict.get<scalar>("pAbs");
    Tsat_ = dict.get<scalar>("Tsat");
    writeFields_ = dict.getOrDefault<bool>("writeFields", false);

    return true;
}


bool Foam::functionObjects::boilingOnsetBerglesRohsenow::execute()
{
    return true;
}


bool Foam::functionObjects::boilingOnsetBerglesRohsenow::write()
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const scalar pBar = pAbs_/1.0e5;
    const scalar expn = 2.16/pow(pBar, 0.0234);
    const scalar coeff = 1082.0*pow(pBar, 1.156);

    for (const label patchi : patchSet_.sortedToc())
    {
        const word& patchName = mesh_.boundary()[patchi].name();
        const scalarField Tw(T.boundaryField()[patchi]);
        const scalarField qw(qWall(patchi));

        scalar minMargin = GREAT;
        label nONB = 0;
        forAll(Tw, facei)
        {
            if (Tw[facei] > Tsat_)
            {
                const scalar qONB =
                    coeff*pow(1.8*(Tw[facei] - Tsat_), expn);
                const scalar margin = qONB - qw[facei];
                minMargin = min(minMargin, margin);
                if (margin < 0)
                {
                    ++nONB;
                }
            }
        }
        reduce(minMargin, minOp<scalar>());
        reduce(nONB, sumOp<label>());
        const scalar maxTw = gMax(Tw);

        Log << type() << " " << name() << " patch " << patchName
            << ": ONB margin(min) = "
            << (minMargin == GREAT ? word("inf (Tw < Tsat everywhere)")
                                   : word(Foam::name(minMargin)))
            << " W/m2, faces in nucleate regime = " << nONB
            << ", max Tw = " << maxTw << " K" << endl;

        writeCurrentTime(file());
        file() << tab << patchName << tab
               << (minMargin == GREAT ? scalar(-1) : minMargin) << tab
               << nONB << tab << maxTw << endl;
    }

    return true;
}


// ************************************************************************* //
