# coOptimiseChannelAndBody — co-optimisation of a coolant path and its conductive body

**STATUS: WORK IN PROGRESS. NOT a verified demonstrator. Do not cite its
results as evidence.** The gradient has the correct *sign* here but does not
yet match finite differences in *magnitude* (see "Verification status"). It is
published because the case, the harness and the findings below are useful, not
because the optimisation is trusted.

## What it is

A generic 2-D cooled slab. Nothing in it comes from any specific component
geometry.

```
 y=0.024  +-----------------------------------------+  heatedSurface (fixed flux in)
          |             crust  (FIXED SOLID)         |
 y=0.020  +-----------------------------------------+
          |          upper body  (DESIGN)            |
 y=0.016  +======+---------------------------+======+
   inlet ->      |      channel (DESIGN)      |      -> outlet
 y=0.008  +======+---------------------------+======+
          |          lower body  (DESIGN)            |
 y=0.004  +-----------------------------------------+
          |              base  (FIXED SOLID)         |
 y=0      +-----------------------------------------+  bottomWall (adiabatic)
        x=0    x=0.012                    x=0.108   x=0.12
               |<-- plug -->|         |<-- plug -->|      (FIXED FLUID)
```

The point: **one design field spans the coolant path and the surrounding
conductive body**, so both co-evolve, while a solid crust under the heated
surface is pinned and the coolant can never be sealed off from its inlet.

| Requirement | How |
|---|---|
| designable fluid **+** solid body region | `adjointPorousZones (designSpace)` — 1536 cells spanning channel and body |
| fixed solid crust | `fixedPorousZones (crustZone baseZone)` + `fixedPorousValues (1 1)` |
| fixed inlet/outlet zones | `fixedZeroPorousZones (inletPlug outletPlug)` (pinned fluid) |
| p-norm hot-wall objective | `patchTemperaturePNorm` (P=8) on `heatedSurface` |
| pressure loss | `PtLosses` on (inlet outlet) — **monitored**, see below |
| design history every N cycles | `writeInterval` in `controlDict`; β iso-surface written **every** cycle |
| iso-surface / STL export | `optimisation/topOIsoSurfaces/*.stl` and `*.vtp`, written by the framework each cycle |

Run: `./Allrun`. Gradient check: `./fd_check.py`.

## Verification status — why this is not yet a demonstrator

`fd_check.py` / `fd_sweep.py` compare the adjoint sensitivity against central
finite differences on single design cells, against the **converged** primal.

First: the FD itself must be trusted. A step-size sweep shows a clean plateau
at eps <= 0.002 (FD stable to 3 significant figures); at eps >= 0.01 the tanh
projection (b = 20) is so nonlinear that the FD can even change sign. Anything
outside the plateau is an artefact, not a finding.

Inside the plateau:

- **Sign: correct.**
- **Magnitude: WRONG, and not yet explained.** With `weight 1e-6` and no
  normalisation, `topOSens` must equal `weight * dJ/dalpha`, so the
  adjoint/FD ratio must be **1e-6 in every cell** — which is what
  `cases/fdcheck` shows (fitted scale 8.93e-7, cells agreeing to 0.4–1.6 %).
  Here the ratio is **1e-5 to 4e-5**: both ~10–40x too large *and* varying by
  3.3x between neighbouring cells.

**What has been ruled out:**

- *Not* FD noise or step size — the plateau is flat and reproducible.
- *Not* an unconverged primal — `residualControl` is met and the FD harness
  refuses any sample whose objective is still moving.
- *Not* cell-volume/geometry bookkeeping — the mesh is uniform (every cell
  1.5 x 1.0 x 1.0 mm), so a `V`/`magSf` convention error cannot produce a
  cell-to-cell spread.
- *Not* the cell-centred vs face-based conductivity sensitivity. The
  sensitivity is now computed **face-based** (differentiating the discrete
  `fvm::laplacian(DEff, T)` operator itself, see the code comment in
  `topOSensMultiplier`), which is the correct discrete derivative and is
  kept. It moved the ratios from 1.244e-5 / 4.132e-5 to 1.074e-5 / 3.678e-5:
  the 3.3x spread is **unchanged**. So the material-interface interpolation is
  not the cause.

**Leading remaining suspect: the regularisation/projection chain.** A pure
interface error would not inflate the overall scale by an order of magnitude.
With `tanh` projection at b = 20, dbeta/dalpha_filtered peaks near 10 and
swings hard wherever the filtered field crosses the 0.5 threshold — which is
exactly the body/channel interface in this case, and never happens in
`cases/fdcheck`, whose design field sits at alpha = 0.05–0.35, far from the
threshold. The next step is to instrument the chain: write `alphaTilda`,
`beta`, and `dbeta/dalphaTilda`, and FD the chain itself
(`dbeta_j/dalpha_i`) rather than the whole objective.

Until that is closed, **nothing from this case is evidence of anything.**

## Findings that came out of building it (these ARE solid)

1. **`kIndicator()` boundary bug (fixed).** The conductivity indicator was
   built with *calculated* patch fields initialised to zero and then
   `correctBoundaryConditions()`, which is a no-op for calculated patches. So
   `I_k = 0` on every boundary face regardless of the design, and `D_eff` at
   any wall was always the *fluid* value. Where solid meets a boundary — a
   heated surface on a pinned crust, exactly this case — the wall heat flux
   was understated by the full solid/fluid ratio: this case rose 0.6 K instead
   of 61 K. Now built with `extrapolatedCalculated` so the indicator
   extrapolates from the adjacent cell.

2. **Adjoint objective BC used the wrong diffusivity (fixed).** The
   `fixedGradient` Ta BC divided by `DEff.boundaryField()[p].patchInternalField()`
   (the adjacent *cell*) instead of the patch *face* value, which is what the
   primal Laplacian actually uses. Now uses the face value.

3. **Thermal sensitivities are explicitly masked** in pinned and non-design
   zones rather than relying on upstream to ignore them. (Upstream does
   restrict the *update* to active design variables, so the written
   `topOSens` still carries upstream's own Brinkman term outside the design
   space — that part is written but never applied.)

4. **Primal convergence gate — the trap that cost the most time.** With
   `limitedLinear` on `div(phi,T)`, T's outer residual sat on a *flat plateau*
   (1.696e-3, unchanged from iteration 100 to 6000) while the linear solve
   converged happily to 1e-13. The primal hit the `nIters` cap with
   `residualControl` never satisfied, and the objective was still drifting:
   J = 368 K at 6000 iterations, 425 K at 40 000, converged ~426 K. Across
   optimisation cycles this *looked exactly like the optimiser going uphill*
   (J climbing 368 → 396 → 410 → 424) and was completely insensitive to the
   step size — because the design step was irrelevant next to the primal
   drift. This is the limiter flip-flop at sharp Brinkman interfaces that the
   top-level README warns about. `bounded Gauss upwind` converges in **171**
   iterations. **Always confirm `residualControl` was met, not just that the
   run finished.**

5. **Adjoint flow relaxation.** The stock upstream `PtLosses` adjoint
   (`adjointSimple`, no thermal coupling at all) stalls at residual ~0.1 in
   this stiff porous geometry with primal-like relaxation, and its
   sensitivities blow up to 1e29. Dropping adjoint velocity relaxation from
   0.6 to **0.3** fixes it (`max|topOSensdp|` → 1.9e-2, sane) with **full
   `standard` ATC retained** — no ATC approximation needed. Adjoint relaxation
   wants the reverse pattern from the primal: low on velocity, high on
   pressure.

## Known-open: the constrained update

With the pressure loss as an active inequality constraint, `nullSpace` SIGFPEs
in its dual Newton solve and `ISQP` stalls, in this geometry. Objective and
constraint gradients differ by ~5 orders (O(1e-7) vs O(1e-2)); normalising the
objective (`normalise`/`normFactor` rather than a tiny `weight`) helps but does
not fix it. The pressure loss is therefore **monitored**, not enforced, here.
The *active* pressure-drop constraint is demonstrated, and works, in
`cases/demo2d`.

Left unconstrained, the optimiser does what you would expect: it drives the
peak surface temperature down (426.3 → 397.4 K over four cycles) by pushing
conductive material into the coolant, and the pressure loss climbs ~700x. That
is a correct response to the objective it was given, and it is exactly why the
constraint matters.
