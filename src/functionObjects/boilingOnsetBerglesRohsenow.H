/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.

Class
    Foam::functionObjects::boilingOnsetBerglesRohsenow

Description
    Checks wall heat flux against the Bergles-Rohsenow onset-of-nucleate-
    boiling (ONB) correlation on selected wall patches:

        q_ONB = 1082 P^1.156 [1.8 (T_w - T_sat)]^(2.16 / P^0.0234)   [W/m2]

    with P the absolute pressure in bar and temperatures in Kelvin (the
    correlation's temperature difference in K equals the original degC
    form). Valid only where T_w > T_sat; walls below T_sat are reported as
    forced-convection (infinite margin).

    Writes, per patch: min ONB margin (q_ONB - q_w) and the number of faces
    in the nucleate-boiling regime (q_w > q_ONB and T_w > T_sat), plus
    optional patch fields qWallBR and qONB.

    Intended for conjugate-heat-transfer verification runs (compressible
    thermo) - e.g. UKAEA ISSF Case Study 9, Annex 1 section 2.4.5.

    Example:
    \verbatim
    boilingCheck
    {
        type            boilingOnsetBerglesRohsenow;
        libs            (thermalTopO);
        region          fluid;
        patches         (fluid_to_solid);
        pAbs            9.5e6;      // absolute pressure [Pa]
        Tsat            580.4;      // saturation temperature [K]
        writeFields     false;
    }
    \endverbatim

SourceFiles
    boilingOnsetBerglesRohsenow.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_boilingOnsetBerglesRohsenow_H
#define Foam_functionObjects_boilingOnsetBerglesRohsenow_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
              Class boilingOnsetBerglesRohsenow Declaration
\*---------------------------------------------------------------------------*/

class boilingOnsetBerglesRohsenow
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- Patches to check
        labelHashSet patchSet_;

        //- Absolute operating pressure [Pa]
        scalar pAbs_;

        //- Saturation temperature at pAbs [K]
        scalar Tsat_;

        //- Write patch fields
        bool writeFields_;


    // Private Member Functions

        //- Wall heat flux [W/m2] on a patch (into the fluid domain boundary)
        tmp<scalarField> qWall(const label patchi) const;


public:

    //- Runtime type information
    TypeName("boilingOnsetBerglesRohsenow");


    // Constructors

        boilingOnsetBerglesRohsenow
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~boilingOnsetBerglesRohsenow() = default;


    // Member Functions

        virtual bool read(const dictionary& dict);

        virtual bool execute();

        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
