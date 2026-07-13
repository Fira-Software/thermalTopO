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
   `fd_sweep_noreg.py` reruns the gate with filter AND projection off
   (ablation A -- see the table below; `project` defaults to `regularise`
   upstream, so turning one off turns both off). Then:

       cell 628: ratio = 1.1286e-6      (FD 1.82170e-2 / 1.82165e-2 / 1.82173e-2)
       cell 626: ratio = 1.1700e-6

   Both equal the objective **weight (1e-6)** and agree with each other to
   3.7 % — the same quality as `cases/fdcheck`. Whatever is broken, it is not
   the thermal adjoint, not the objective, and not the conductivity term.
6. **The tanh projection derivative, on its own.** Ablation B (projection on,
   filter off) passes at 1.297e-6 / 1.311e-6. So the projection derivative is
   not wrong by itself.
7. **`postProcessSens` IS the exact transpose of the true beta(alpha) map.**
   `fd_chain.py` is a chain-only test with no PDE physics: freeze the
   pre-chain field `gbeta = topologySens`, measure the true Jacobian
   `dbeta_i/dalpha_k` by perturbing raw alpha and reading the solver's own
   `beta`, and compare `sum_i gbeta_i * dbeta_i/dalpha_k` against `topOSens_k`.
   Result: the two agree to 5 significant figures for both probe cells, up to
   a single global factor of 6.6667e8 = **1/V** (the volume weight
   `postProcessSens` applies at the end). A single raw-alpha perturbation
   correctly spreads to ~420 beta cells, so the filter really is active.

### The full ablation: only the COMBINATION fails

Filter and projection are separate flags upstream
(`project_(dict.getOrDefault<bool>("project", regularise_))`, so `project`
*defaults to* `regularise` — a trap: switching `regularise` off silently
switches the projection off too). Running all four combinations
(`fd_ablation_B.py`, `fd_ablation_C.py`; eps = 1e-3):

| | filter | projection | cell 628 | cell 626 | verdict |
|---|---|---|---|---|---|
| **A** | off | off | 1.13e-6 | 1.17e-6 | OK — equals the weight |
| **B** | off | **on** | 1.297e-6 | 1.311e-6 | OK (cells agree to 1.1 %) |
| **C** | **on** | off | 8.92e-7 | 9.14e-7 | OK (cells agree to 2.5 %) |
| **D** | **on** | **on** | 1.07e-5 | 3.68e-5 | **FAILS: 10–40x off, 3.3x spread** |

The filter alone is fine. The projection alone is fine. **Only the two together
fail.** So neither operator is individually wrong, and the earlier readings
("it's the projection", then "it's the filter") were both wrong.

The one thing unique to D: beta is the *filtered-then-sharpened* field, so it
is the only configuration carrying genuinely **intermediate (grey) beta at the
interfaces**. In A and B beta is near-binary; in C beta is smooth but never
sharpened. Note also that C can *mask* a pointwise error, because the filter
transpose smooths it away, whereas in D the projection derivative (tanh, b=20,
peaking near 10) amplifies exactly the cells where beta is intermediate. So the
error most likely lives in the sensitivity **at intermediate beta**, and D is
simply the only case that both produces it and amplifies it.

The relevant upstream code, `fieldRegularisation::postProcessSens`
(`.../topODesignVariables/regularisation/fieldRegularisation.C:154`):

    if (project_)    sens *= sharpenFunction_->derivative(betaArg_);  // dbeta/dalphaTilda
    if (regularise_) regularise(sens, sens, false);                   // filter transpose
    sens *= mesh_.V();                                                // volume

against the forward map in `updateBeta()`:

    if (regularise_) regularise(alpha_, alphaTilda_(), true);         // alphaTilda = H alpha
    if (project_)    sharpenFunction_->interpolate(betaArg_, beta_);  // beta = P(alphaTilda)

The ordering is correct (forward H then P; transpose P' then H^T), and
`fd_chain.py` confirms `postProcessSens` reproduces the *measured* Jacobian
transpose to 5 s.f. So the composition is right, which is what makes the D
failure interesting rather than trivial.

### Two more hypotheses killed (13 Jul, later)

**FD secant vs tangent through the steep tanh — DEAD.** If the FD step were
crossing the sharp projection, the adjoint (a local tangent) would exceed the
FD (a secant), and the ratio would fall toward the weight as eps -> 0. It does
not. Deep sweep in case D, eps from 3e-5 to 3e-3:

    cell 628:  9.87e-6, 1.071e-5, 1.073e-5, 1.077e-5, 1.101e-5
    cell 626:  3.58e-5, 3.60e-5,  3.61e-5,  3.61e-5,  3.85e-5

Flat across two decades, no trend toward 1e-6. It is a real chain error, not
FD resolution.

**Chain order / double application — NOT VISIBLE IN THE CODE.** The chain is
fully determined and is correct:

- `topOVariablesBase::sourceTerm` interpolates with `this->beta()`, the
  PROJECTED field, so `Ik = K(P(H alpha))` -- filter, then projection, then
  conductivity interpolation. That is the physically right order (not
  `H(K(...))`).
- `sourceTermSensitivities` is literally `sens *= betaMax * K'(beta)`
  (betaMax = 1 for us). Nothing else.
- `postProcessSens` then applies `P'(alphaTilda)`, then `H^T`, then `* V`.

so `topOSens = V * H^T * diag(P'(alphaTilda)) * diag(K'(beta)) * dJ/dIk`,
which IS `dIk/dalpha = K'(beta) * P'(alphaTilda) * H` transposed. No double
application, correct order. (The mesh is uniform, so `*V` commuting past `H^T`
is harmless.)

**The smoking-gun table, which half-fires.** Both probe cells sit essentially
ON the projection threshold, where P' is at its maximum of ~10:

| cell | alphaTilda | beta | P'(alphaTilda) | K'(beta) | ratio/weight |
|---|---|---|---|---|---|
| 628 | 0.5019 | 0.5193 | 9.99 | 0.162 | 10.7 |
| 626 | 0.5158 | 0.6528 | 9.07 | 0.106 | 36.1 |

Cell 628's ratio/weight matches its own P' almost exactly. Cell 626's does not.

CORRECTION (an earlier version of this file argued that 36.1 > max(P') = 10 was
impossible and therefore paradoxical -- that was WRONG). topOSens is a
filter-weighted sum over NEIGHBOURS,

    topOSens_k = V * sum_i  H_ik * P'_i * K'_i * g_i

so the amplification relative to the probe cell's own P' can exceed max(P')
freely. There is no paradox.

### Localisation: it is the near-zero-beta cells

`fd_localise.py` ranks the actual contributions, using the numerically measured
Jacobian dbeta_i/dalpha_k. For probe cell 626:

     cell    zone     beta   P'(aT)  dbeta/da       gbeta     contrib     %
      786  design   0.0002    0.007    0.0001  -1.618e+05  -2.17e+01  31.8
      785  design   0.0002    0.008    0.0001  -1.348e+05  -1.42e+01  20.7
      787  design   0.0001    0.005    0.0001  -1.901e+05  -1.35e+01  19.8
      626  design   0.6535    9.074    0.6596   1.704e+01   1.12e+01 -16.4  <- probe cell
      708  design   0.0036    0.141    0.0014  -8.197e+03  -1.12e+01  16.4

The probe cell contributes only -16 %, with the OPPOSITE sign to the total. The
sum is dominated by **near-zero-beta fluid cells** (beta ~ 1e-4) whose pre-chain
sensitivity gbeta is about **10^4 times larger** than the probe cell's
(-1.6e5 vs +17).

Those are exactly the cells where the BorrvallPetersson derivative
`K'(beta) = q(1+q)/(beta+q)^2` and the Brinkman `a'(beta) = betaMax * I'(beta)`
are maximal (q = 1/b = 0.05, betaMax = 2500).

**And no ablation has ever tested gbeta there.** A and B only probe the
*perturbed* cell's own gbeta (the chain is local, so neighbours never enter),
and in C the projection is off, so beta never gets driven to ~1e-4 -- the
filtered field alphaTilda stays moderate in the channel. Case D is the only one
that both creates near-zero beta (projection sharpens alphaTilda ~ 0.1 down to
~1e-4) AND pulls those cells into the probe's stencil (filter). That is why only
D fails.

**Fixed-zone ordering is NOT the cause.** beta does move in 48-52 fixed-zone
cells under a raw-alpha perturbation, but their combined contribution is
+0.01 to +0.02 % of the total. Ruled out.

### Next diagnostic

Verify the pre-chain sensitivity `gbeta` **in a near-zero-beta cell**. That
region is untested by every ablation so far, and it dominates the failing sum.
It cannot be FD'd directly by perturbing its own raw alpha (those channel cells
sit at alpha = 0, so the minus side clips), so use one of:

- a design whose channel cells start at alpha = 0.05 rather than 0, so a
  central difference is admissible there;
- or split the written sensitivity into its two parts --

      momPostChain   (from adjointSimple::topOSensMultiplier, the Brinkman term)
      thermPostChain (ours, after sourceTermSensitivities)

  and check which of the two carries the 1e5 magnitude at beta ~ 1e-4. The
  Brinkman term is upstream code and scales as betaMax * I'(beta) = 2500 * 21 at
  beta -> 0, so it is the more likely carrier -- in which case this is NOT a
  thermalTopO bug at all.

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
