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

### What has been ruled out (each by measurement, not argument)

Reproduce with `fd_sweep.py`, `fd_sweep_noreg.py`, `fd_chain.py`.

1. **FD noise / step size.** `fd_sweep.py` shows a flat plateau at
   eps <= 0.002 (stable to 3 s.f.). The eps >= 0.01 artefacts — including a
   *sign flip* — were the step size, not the code.
2. **Unconverged primal.** `residualControl` is met and the harness rejects
   any sample whose objective is still moving.
3. **Volume / geometry bookkeeping.** The mesh is uniform (every cell
   1.5 x 1.0 x 1.0 mm), so a `V`/`magSf` convention error cannot produce a
   *cell-to-cell* spread.
4. **Cell-centred vs face-based conductivity sensitivity.** The sensitivity is
   now computed **face-based**, differentiating the discrete
   `fvm::laplacian(DEff, T)` operator itself (see the derivation in the code
   comment in `topOSensMultiplier`). It is the correct discrete derivative and
   is kept — but it moved the ratios only from 1.244e-5 / 4.132e-5 to
   1.074e-5 / 3.678e-5. **The 3.3x spread is unchanged**, so the material-
   interface interpolation is NOT the cause.
5. **The thermal sensitivity itself (pre-chain dJ/dbeta) is CORRECT.**
   `fd_sweep_noreg.py` reruns the gate with the Helmholtz filter off. Then:

       cell 628: ratio = 1.1286e-6      (FD 1.82170e-2 / 1.82165e-2 / 1.82173e-2)
       cell 626: ratio = 1.1700e-6

   Both equal the objective **weight (1e-6)** and agree with each other to
   3.7 % — the same quality as `cases/fdcheck`. Whatever is broken, it is not
   the thermal adjoint, not the objective, and not the conductivity term.
6. **The tanh projection derivative.** In the filter-off ablation the
   projection is still active, and `dbeta/dalphaTilda` cancels exactly between
   `topOSens` and the FD (it appears once in the forward map and once in the
   transpose). That ablation passes, so the projection derivative is fine.
7. **`postProcessSens` IS the exact transpose of the true beta(alpha) map.**
   `fd_chain.py` is a chain-only test with no PDE physics: freeze the
   pre-chain field `gbeta = topologySens`, measure the true Jacobian
   `dbeta_i/dalpha_k` by perturbing raw alpha and reading the solver's own
   `beta`, and compare `sum_i gbeta_i * dbeta_i/dalpha_k` against `topOSens_k`.
   Result: the two agree to 5 significant figures for both probe cells, up to
   a single global factor of 6.6667e8 = **1/V** (the volume weight
   `postProcessSens` applies at the end). A single raw-alpha perturbation
   correctly spreads to ~420 beta cells, so the filter really is active.

### The remaining suspect: the filter path

Turning the Helmholtz filter on is what breaks it:

| cell | ratio, filter OFF | ratio, filter ON | inflation |
|---|---|---|---|
| 628 | 1.13e-6 | 1.07e-5 | **9.5x** |
| 626 | 1.17e-6 | 3.68e-5 | **31.5x** |

The relevant upstream code is `fieldRegularisation::postProcessSens`
(`.../topODesignVariables/regularisation/fieldRegularisation.C:154`):

    if (project_)    sens *= sharpenFunction_->derivative(betaArg_);  // dbeta/dalphaTilda
    if (regularise_) regularise(sens, sens, false);                   // filter transpose
    sens *= mesh_.V();                                                // volume

against the forward map in `updateBeta()`:

    if (regularise_) regularise(alpha_, alphaTilda_(), true);         // filter
    if (project_)    sharpenFunction_->interpolate(betaArg_, beta_);  // projection

Note `project_` and `regularise_` are **separate flags**: the filter-off
ablation above still projects.

**Open logical tension, for whoever picks this up.** Test 7 says the transpose
is exact; test 5 says the pre-chain sensitivity is correct; yet composing them
fails. One of the two must be state-dependent — they were measured on
*different* beta fields (near-binary with the filter off, smooth with it on).
That is the thread to pull.

**Next experiment (ablation C):** filter ON, projection OFF
(`regularise true`, `project false`). That isolates the filter transpose with
nothing else in the path. If the ratio stays at 1e-6 the bug is in the
interaction; if it inflates, it is the Helmholtz transpose or its volume
weighting.

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
