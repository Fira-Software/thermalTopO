# multiSolidMaterials

**Different solid materials in the domain, throughout the topology optimisation.**

**Verified example:** bundled checks assert that each labelled zone and the
generated solid use the intended material property fields.

The design variable chooses **fluid vs solid**. It does **not** choose *which*
solid. Material labels are **fixed** for the whole optimisation, assigned by
cellZone, and every material carries its own diffusivity, constant or tabulated
in temperature. Design cells that solidify take
`defaultGeneratedSolidMaterial`.

This is deliberately *not* multi-material topology optimisation (the optimiser
selecting between tungsten and copper). It is the far more common industrial
case: a component made of several fixed materials, inside which a coolant path
is being optimised.

## Setup

```
solidMaterials
{
    tungsten { DSolidTable table ((300 1.2e-3) (340 0.8e-3)); }   // D_s(T)
    CuCrZr   { DSolid 4.0e-4; }                                    // constant
    copper   { DSolid 2.0e-3; }
}
solidMaterialZones
{
    CuCrZr  (CuCrZrZone);
    copper  (copperZone);
}
defaultGeneratedSolidMaterial tungsten;
```

Unlabelled cells take the default material. The solver reports:

```
thermalTopO: 3 solid materials 3(tungsten CuCrZr copper); generated solid = tungsten
```

## Why the sensitivities need no change

`thermalPropertyTables::DSolidField()` already returns D_s as a **field**, and
the conductivity sensitivity already consumes it as one:

    dJ/dbeta  ~  I_k'(beta) * (D_s - D_f - nu_t/Pr_t) * ...

Once `D_s` is locally material-dependent, the topology sensitivity is
material-aware for free. No change to the adjoint.

## Run and check

    ./Allrun

`writeThermalProperties` writes the assembled `DSolid`, `DFluid` and
`materialID` fields, and `checkMaterials.py` asserts that each labelled zone
really uses its own material:

```
                      zone  cells     mean Ds    expected  result
                    CuCrZr     60  4.0000e-04  4.0000e-04  OK
                    copper     60  2.0000e-03  2.0000e-03  OK
   designSpace(->tungsten)    720  1.1070e-03  1.1070e-03  OK (tabulated D_s(T))
         tungsten(default)    600  1.1135e-03  1.1135e-03  OK (tabulated D_s(T))

PASS
```

The generated solid is tungsten and stays temperature-dependent; the fixed
zones keep their own properties, throughout the optimisation.

## Scope

Verified: the property field is assembled correctly per material, per cell, at
the local temperature, and is used by the primal, the adjoint and the
sensitivity assembly. For body-fitted re-verification of a finished design, use
OpenFOAM's pre-existing multi-region conjugate-heat-transfer capability with one
solid region per material.
