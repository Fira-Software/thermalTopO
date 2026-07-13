/*---------------------------------------------------------------------------*\
    thermalTopO: thermal extensions to OpenFOAM's adjointOptimisationFoam

    Copyright (C) 2026 Fira Software Ltd
    SPDX-License-Identifier: GPL-3.0-or-later

    This offering is not approved or endorsed by OpenCFD Limited, producer
    and distributor of the OpenFOAM software via www.openfoam.com, and owner
    of the OPENFOAM and OpenCFD trade marks.
\*---------------------------------------------------------------------------*/

#include "thermalPropertyTables.H"
#include "fvMesh.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::thermalPropertyTables::thermalPropertyTables()
:
    DFluid_(Zero),
    DSolid_(Zero),
    rhoCpRef_(Zero),
    TRef_(Zero),
    rhoPtr_(nullptr),
    cpPtr_(nullptr),
    kFluidPtr_(nullptr),
    DSolidPtr_(nullptr),
    multiMaterial_(false),
    materialNames_(),
    materialConst_(),
    materialTable_(),
    materialID_(),
    defaultMaterial_(0)
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::uniform
(
    const volScalarField& T,
    const word& name,
    const scalar value,
    const dimensionSet& dims
) const
{
    const fvMesh& mesh = T.mesh();

    return tmp<volScalarField>::New
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(name, dims, value)
    );
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::evaluate
(
    const Function1<scalar>& f,
    const volScalarField& T,
    const word& name,
    const scalar scale,
    const dimensionSet& dims
) const
{
    tmp<volScalarField> tfld(uniform(T, name, Zero, dims));
    volScalarField& fld = tfld.ref();

    scalarField& fi = fld.primitiveFieldRef();
    const scalarField& Ti = T.primitiveField();
    forAll(fi, celli)
    {
        fi[celli] = scale*f.value(Ti[celli]);
    }

    // Boundary faces are evaluated at the boundary temperature, not left at
    // the internal or a constant value: on a patch where solid meets the
    // domain boundary (a heated wall, say) the wrong face diffusivity would
    // otherwise be used in the Laplacian flux.
    volScalarField::Boundary& bfld = fld.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchi];
        fvPatchScalarField& fp = bfld[patchi];

        forAll(fp, facei)
        {
            fp[facei] = scale*f.value(Tp[facei]);
        }
    }

    return tfld;
}


Foam::scalar Foam::thermalPropertyTables::materialDSolid
(
    const label matID,
    const scalar T
) const
{
    if (materialTable_.set(matID))
    {
        return materialTable_[matID].value(T);
    }

    return materialConst_[matID];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::thermalPropertyTables::read
(
    const dictionary& thermalDict,
    const fvMesh& mesh
)
{
    rhoPtr_.reset(nullptr);
    cpPtr_.reset(nullptr);
    kFluidPtr_.reset(nullptr);
    DSolidPtr_.reset(nullptr);

    if (thermalDict.found("rhoTable") != thermalDict.found("cpTable"))
    {
        FatalIOErrorInFunction(thermalDict)
            << "rhoTable and cpTable must be given together: the equation is "
            << "scaled by the volumetric heat capacity rho*cp, so one without "
            << "the other is not meaningful."
            << exit(FatalIOError);
    }

    if (thermalDict.found("rhoTable"))
    {
        rhoPtr_ = Function1<scalar>::New("rhoTable", thermalDict, &mesh);
        cpPtr_ = Function1<scalar>::New("cpTable", thermalDict, &mesh);
    }
    if (thermalDict.found("kFluidTable"))
    {
        kFluidPtr_ = Function1<scalar>::New("kFluidTable", thermalDict, &mesh);
    }
    if (thermalDict.found("DSolidTable"))
    {
        DSolidPtr_ = Function1<scalar>::New("DSolidTable", thermalDict, &mesh);
    }

    // Constants: required only where the corresponding table is absent
    DFluid_ =
        kFluidPtr_
      ? thermalDict.getOrDefault<scalar>("DFluid", Zero)
      : thermalDict.get<scalar>("DFluid");

    DSolid_ =
        DSolidPtr_
      ? thermalDict.getOrDefault<scalar>("DSolid", Zero)
      : thermalDict.get<scalar>("DSolid");

    // Reference volumetric heat capacity
    if (thermalDict.found("rhoCpRef"))
    {
        rhoCpRef_ = thermalDict.get<scalar>("rhoCpRef");
        TRef_ = thermalDict.getOrDefault<scalar>("TRef", Zero);
    }
    else if (rhoPtr_)
    {
        TRef_ = thermalDict.get<scalar>("TRef");
        rhoCpRef_ = rhoPtr_->value(TRef_)*cpPtr_->value(TRef_);
    }
    else if (kFluidPtr_)
    {
        FatalIOErrorInFunction(thermalDict)
            << "kFluidTable gives k [W/m/K], which is converted to a "
            << "diffusivity by (rho cp)_ref. Supply either rhoCpRef, or "
            << "rhoTable + cpTable + TRef."
            << exit(FatalIOError);
    }

    if (rhoCpRef_ < SMALL && (rhoPtr_ || kFluidPtr_))
    {
        FatalIOErrorInFunction(thermalDict)
            << "Non-positive reference volumetric heat capacity "
            << rhoCpRef_ << exit(FatalIOError);
    }

    // ---- optional: several solid materials, labelled by cellZone -----------
    // The design variable chooses fluid vs SOLID. It does NOT choose WHICH
    // solid: material labels are fixed for the whole optimisation, and design
    // cells that solidify take defaultGeneratedSolidMaterial.
    multiMaterial_ = thermalDict.found("solidMaterials");

    if (multiMaterial_)
    {
        const dictionary& matDict = thermalDict.subDict("solidMaterials");

        materialNames_ = matDict.toc();
        const label nMat = materialNames_.size();

        if (!nMat)
        {
            FatalIOErrorInFunction(thermalDict)
                << "solidMaterials is empty." << exit(FatalIOError);
        }

        materialConst_.setSize(nMat, Zero);
        materialTable_.clear();
        materialTable_.setSize(nMat);

        forAll(materialNames_, i)
        {
            const dictionary& md = matDict.subDict(materialNames_[i]);

            if (md.found("DSolidTable"))
            {
                materialTable_.set
                (
                    i,
                    Function1<scalar>::New("DSolidTable", md, &mesh).ptr()
                );
            }
            else
            {
                materialConst_[i] = md.get<scalar>("DSolid");
            }
        }

        // material assigned to generated solid
        const word defName
        (
            thermalDict.get<word>("defaultGeneratedSolidMaterial")
        );
        defaultMaterial_ = materialNames_.find(defName);

        if (defaultMaterial_ == -1)
        {
            FatalIOErrorInFunction(thermalDict)
                << "defaultGeneratedSolidMaterial " << defName
                << " is not one of solidMaterials " << materialNames_
                << exit(FatalIOError);
        }

        // label every cell: default first, then override by zone
        materialID_.setSize(mesh.nCells(), defaultMaterial_);

        if (thermalDict.found("solidMaterialZones"))
        {
            const dictionary& zd = thermalDict.subDict("solidMaterialZones");

            for (const word& matName : zd.toc())
            {
                const label matID = materialNames_.find(matName);

                if (matID == -1)
                {
                    FatalIOErrorInFunction(thermalDict)
                        << "solidMaterialZones names material " << matName
                        << ", which is not in solidMaterials "
                        << materialNames_ << exit(FatalIOError);
                }

                for (const word& zoneName : zd.get<wordList>(matName))
                {
                    const label zoneID =
                        mesh.cellZones().findZoneID(zoneName);

                    if (zoneID == -1)
                    {
                        FatalIOErrorInFunction(thermalDict)
                            << "No cellZone " << zoneName
                            << " for solid material " << matName
                            << exit(FatalIOError);
                    }

                    for (const label celli : mesh.cellZones()[zoneID])
                    {
                        materialID_[celli] = matID;
                    }
                }
            }
        }

        Info<< "thermalTopO: " << nMat << " solid materials "
            << materialNames_ << "; generated solid = " << defName << endl;
    }

    if (rhoPtr_ || kFluidPtr_ || DSolidPtr_)
    {
        Info<< "thermalTopO: temperature-dependent properties active"
            << " [fluid k: " << Switch(bool(kFluidPtr_))
            << ", fluid rho*cp: " << Switch(bool(rhoPtr_))
            << ", solid D: " << Switch(bool(DSolidPtr_)) << "]" << nl;

        if (rhoPtr_ || kFluidPtr_)
        {
            Info<< "    (rho cp)_ref = " << rhoCpRef_ << " J/m3/K";
            if (rhoPtr_)
            {
                Info<< " at TRef = " << TRef_ << " K";
            }
            Info<< endl;
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::DFluidField
(
    const volScalarField& T,
    const word& name
) const
{
    if (kFluidPtr_)
    {
        return evaluate(*kFluidPtr_, T, name, scalar(1)/rhoCpRef_, dimViscosity);
    }

    return uniform(T, name, DFluid_, dimViscosity);
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::DSolidField
(
    const volScalarField& T,
    const word& name
) const
{
    if (multiMaterial_)
    {
        // Spatially varying solid diffusivity: each cell uses its own labelled
        // material, evaluated at the local temperature. The topology
        // sensitivity term needs no change, because it already consumes D_s as
        // a FIELD (Ds - Df - nut/Prt).
        tmp<volScalarField> tDs(uniform(T, name, Zero, dimViscosity));
        volScalarField& Ds = tDs.ref();

        scalarField& Dsi = Ds.primitiveFieldRef();
        const scalarField& Ti = T.primitiveField();
        forAll(Dsi, celli)
        {
            Dsi[celli] = materialDSolid(materialID_[celli], Ti[celli]);
        }

        volScalarField::Boundary& bDs = Ds.boundaryFieldRef();
        forAll(bDs, patchi)
        {
            const labelUList& fc = T.mesh().boundary()[patchi].faceCells();
            const fvPatchScalarField& Tp = T.boundaryField()[patchi];
            fvPatchScalarField& Dsp = bDs[patchi];

            forAll(Dsp, facei)
            {
                Dsp[facei] =
                    materialDSolid(materialID_[fc[facei]], Tp[facei]);
            }
        }

        return tDs;
    }

    if (DSolidPtr_)
    {
        return evaluate(*DSolidPtr_, T, name, scalar(1), dimViscosity);
    }

    return uniform(T, name, DSolid_, dimViscosity);
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::materialIDField
(
    const volScalarField& T,
    const word& name
) const
{
    tmp<volScalarField> tId(uniform(T, name, Zero, dimless));

    if (multiMaterial_)
    {
        scalarField& f = tId.ref().primitiveFieldRef();
        forAll(f, celli)
        {
            f[celli] = scalar(materialID_[celli]);
        }
    }

    return tId;
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::CField
(
    const volScalarField& T,
    const word& name
) const
{
    if (!rhoPtr_)
    {
        // Constant properties: C == 1, and the formulation collapses to the
        // constant-property one exactly.
        return uniform(T, name, scalar(1), dimless);
    }

    tmp<volScalarField> tC(uniform(T, name, scalar(1), dimless));
    volScalarField& C = tC.ref();

    const scalar rrhoCpRef = scalar(1)/rhoCpRef_;

    scalarField& Ci = C.primitiveFieldRef();
    const scalarField& Ti = T.primitiveField();
    forAll(Ci, celli)
    {
        const scalar Tc = Ti[celli];
        Ci[celli] = rrhoCpRef*rhoPtr_->value(Tc)*cpPtr_->value(Tc);
    }

    volScalarField::Boundary& bC = C.boundaryFieldRef();
    forAll(bC, patchi)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchi];
        fvPatchScalarField& Cp = bC[patchi];

        forAll(Cp, facei)
        {
            const scalar Tf = Tp[facei];
            Cp[facei] = rrhoCpRef*rhoPtr_->value(Tf)*cpPtr_->value(Tf);
        }
    }

    return tC;
}


Foam::tmp<Foam::volScalarField> Foam::thermalPropertyTables::alphaField
(
    const volScalarField& T,
    const word& name
) const
{
    tmp<volScalarField> tAlpha(DFluidField(T, name));

    if (rhoPtr_)
    {
        tAlpha.ref() /= CField(T, "alphaRhoCpNorm")();
    }

    return tAlpha;
}


// ************************************************************************* //
