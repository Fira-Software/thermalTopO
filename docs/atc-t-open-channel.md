# Verifying the thermal-to-momentum adjoint coupling in open-channel flow

**Working note. Fira Software Ltd, 13 July 2026.**

This note records, in full, the verification campaign behind `examples/
coOptimiseChannelAndBody`. It is published because the method is reusable and
because the finding is exactly the kind of capability gap that open verification
is supposed to surface: a continuous-adjoint coupling term that is correct in the
smooth limit, passes every gradient test in a porous-medium regime, and is only
exposed by an open channel.

**Summary.** The (ATC-T) source that couples a thermal objective into the adjoint
momentum equation is the continuous form `Ta grad(T)`. That is correct in the
smooth, low-velocity limit -- which is why it verifies to 0.2-1.6% against finite
differences in every case in `cases/`, all of which place the design region in a
low-velocity Brinkman sponge. In an OPEN CHANNEL it is not the discrete transpose
of the primal convection operator, because OpenFOAM's `bounded Gauss <scheme>`
assembles `div(phi,T) - Sp(div(phi),T)`: at convergence `div(phi) ~ 0`, so the
primal RESIDUAL is unaffected, but the DERIVATIVE with respect to the face flux is
not. The exact flux sensitivity has been derived, and PROVEN against a
finite-difference of the face Lagrangian to 10-11 significant figures:

    gPhi_f = T_f (A_P - A_N) - chi (A_P T_P - A_N T_N),    A_P = C_P Ta_P

with chi = 1 for a bounded scheme. `utilities/testAdjointTranspose` is the
no-solve gate that checks it. The remaining work is to inject `gPhi_f` through
the framework's Rhie-Chow flux linearisation rather than as a cell momentum
source. Until that is complete, `examples/coOptimiseChannelAndBody` is an
experimental case and its gradients are not used for any quantitative claim.

The verified capability in `cases/` is unaffected: those campaigns are unchanged
and independently reproducible.

---

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

### It is NOT upstream OpenFOAM, and NOT the coupling sign

**Stock upstream passes the same gate.** `fd_flowonly.py` runs the identical
mesh, zones, filter and projection with the STOCK solvers only -- `simple` +
`adjointSimple` + `PtLosses` -- with no thermalTopO in the sensitivity path at
all:

    cell 628:  adj/FD = 1.359e-6, 1.359e-6, 1.359e-6   (eps 5e-4 .. 2e-3)
    cell 626:  adj/FD = 1.218e-6, 1.218e-6, 1.218e-6

Both equal the objective weight. So the upstream Brinkman sensitivity, the
filter, the projection and their composition are all CORRECT for a flow
objective in this geometry. The bug is not upstream.

**The coupling sign is right.** Flipping `couplingSign` to -1 in case D gives
NEGATIVE ratios (-1.86e-5, -6.01e-5), i.e. worse. (The FD is identical for both
signs, confirming `couplingSign` only affects the adjoint.)

### What that leaves: the thermal -> momentum adjoint coupling

The Brinkman sensitivity is `a'(beta) * (u . u_a)`. Its formula is upstream and
now proven correct -- but the `u_a` it consumes is produced by OUR coupling term
in `addMomentumSource()`:

    matrix += couplingSign_ * (C(T) * Ta * grad(T))      // the ATC-T term

`a'(beta) = betaMax * I'(beta)` is maximal as beta -> 0, so open-fluid cells
amplify `u_a` enormously (a' = 2500 * 21 at beta = 0). Those are exactly the
cells shown above to dominate the failing sum, and they carry ~100 % Brinkman /
~0 % thermal sensitivity.

**And this term has never been exercised where it matters.** `cases/fdcheck`,
`cases/varprops` and `cases/demo2d` all put the design region in a low-velocity
Brinkman *sponge* (alpha 0.05-0.35), so `u` is small, `u . u_a` is small, and
the Brinkman contribution is negligible -- the thermal term dominates and the FD
gate passes. This example is the first with an **open channel** carrying real
velocity through the design space. If the ATC-T coupling is wrong in magnitude
or form (rather than sign), this is the only case in the repo that would ever
reveal it.

Decomposing the sensitivity by coupling sign (the Brinkman term is linear in
`u_a`, which is linear in the coupling source):

    cell 628:  S_nocoupling/FD = -4.1e-6    S_coupling/FD = +1.45e-5
    cell 626:  S_nocoupling/FD = -1.2e-5    S_coupling/FD = +4.8e-5

Neither part equals the weight, so both the coupled and uncoupled halves are off
at these cells.

### ATC-T: sign AND continuous form both tested. Neither is the fix.

**Label correction.** Setting `thermalSensScale 0` removes only our **direct
conductivity sensitivity**. It does NOT remove thermalTopO from the picture: the
thermal objective still drives the adjoint momentum through
`addMomentumSource()`, and the upstream Brinkman formula then consumes that
`u_a`. So the split above is

    direct conductivity sensitivity        (thermalTopO, the term we add)
    thermal-INDUCED Brinkman sensitivity   (upstream formula, OUR u_a)

and NOT "ours vs upstream". The upstream formula is proven correct by the
flow-only test; the suspect was always the `u_a` our coupling produces.

**Both continuous forms fail.** The two candidate ATC-T sources differ by a pure
gradient (`Ta grad(T) = grad(Ta T) - T grad(Ta)`), which an incompressible
adjoint can absorb into adjoint pressure in the continuum but not necessarily
discretely. Both are now implemented behind `couplingForm` and both were run in
case D (eps = 1e-3; target ratio = 1e-6):

| couplingForm | cell 628 | cell 626 |
|---|---|---|
| `TaGradT` (default) | 1.04e-5 | 3.55e-5 |
| `negTGradTa` | 1.24e-5 | 4.52e-5 |

`negTGradTa` is slightly WORSE. Combined with the earlier `couplingSign = -1`
test (which gave negative ratios), the discrepancy is largely **insensitive to
the ATC-T form**: sign flip, form swap, and default all land in the same 10-45x
band. So "the continuous coupling source is simply the wrong expression" does
not explain it either.

The `couplingForm` switch is retained -- it is cheap, documented, and the right
place to plug in a discrete face-based transpose when one is written.

### Null test PASSES (no contamination)

With `couplingForm none` + `thermalSensScale 0` and only the thermal objective
active (no flow objective, no PtLosses adjoint in the sensitivity field):

    1/Ua          = uniform (0 0 0)     EXACTLY
    1/topOSensas1 = uniform 0           EXACTLY

So the adjoint velocity, and therefore the entire Brinkman contribution to the
thermal objective, is driven **solely** by our ATC-T coupling. The comparison is
not contaminated by a stale Ua, a pressure-loss adjoint or a default objective
source. `couplingForm none` is retained as a permanent guard.

### Convection scheme: error shrinks with diffusivity, but the scheme test is INCONCLUSIVE

The case uses `bounded Gauss upwind` for `div(phi,T)` and `div(-phi,Ta)` (forced
by the limiter stall). The discrete transpose of an UPWIND convection matrix is a
downwind operator, not the reversed-flux upwind operator, so
`fvm::div(-phi,Ta)` is not guaranteed to be the algebraic transpose of
`fvm::div(phi,T)`. To test, the case was smoothed (D_f 1e-5 -> 1e-4,
D_s 1e-3 -> 1e-2, contrast unchanged; cell Peclet 3 -> 0.3) so that a
near-self-adjoint scheme is admissible:

| convection scheme | cell 628 | cell 626 |
|---|---|---|
| upwind (smoothed) | 5.73e-6 | 2.79e-6 |
| linear (smoothed) | 5.28e-6 | 2.79e-6 |
| *(original, Pe ~3, upwind)* | *1.04e-5* | *3.55e-5* |

Two readings:

1. **The scheme barely matters here -- but the test is WEAK.** At Pe ~ 0.3 upwind
   and linear nearly coincide, so this cannot discriminate between them. A real
   scheme test needs HIGH Pe with a self-adjoint scheme, which is exactly where
   central convection oscillates. The scalar-transpose hypothesis is therefore
   **not killed; it is untested.**

2. **Raising the diffusivity 10x collapsed the error from 10-36x to 2.8-5.7x.**
   That confirms the PATH: more diffusion -> weaker grad(T) -> smaller Ta ->
   weaker ATC-T source -> smaller u_a -> smaller thermal-induced Brinkman term.
   The error scales with the strength of exactly the path localised above.

But it remains **3-6x off even in a smooth, low-Peclet, nearly self-adjoint
configuration**, so the fault is not merely an upwind-transpose artefact.
Something in the ATC-T -> Brinkman path is wrong independently of the convection
scheme.

### ROOT CAUSE FOUND: the adjoint momentum source is not the discrete transpose

`utilities/testAdjointTranspose` is a pure operator test: no optimisation, no
primal or adjoint solve. It reads a frozen (phi, T, Ta) state, assembles the
primal convection matrix `A = fvm::div(phi,T)` with the case's own scheme, and
compares against what the adjoint actually applies.

**(1) The SCALAR convection transpose is CORRECT. Hypothesis dead.**

    max|B.upper - A.lower| / scale = 0        (exactly)
    max|B.lower - A.upper| / scale = 0        (exactly)
    full A^T Ta vs div(-phi,Ta), interior cells: relative 1.6e-8

So `fvm::div(-phi,Ta)` IS the exact algebraic transpose of `fvm::div(phi,T)` for
upwind -- the reversed-flux upwind operator has exactly the transposed
off-diagonals. The residual 1.6e-8 is the `bounded` Sp(div(phi)) term, which
enters A and B with opposite signs (max|div(phi)| = 6e-8, i.e. continuity error).
`Ta` is therefore the correct discrete adjoint temperature.

(An earlier run of this test reported "relative 4.27" and looked damning. That
was WRONG: it included the diagonal, which differs *legitimately* on boundary
cells because T has `fixedValue 300` at the inlet while Ta has `fixedValue 0` --
the correct adjoint BC, not an error.)

**(2) The MOMENTUM SOURCE is NOT the transpose. This is the defect.**

The exact discrete derivative of the primal convective term contracted with Ta is

    d/dphi_f [ Ta^T R_conv ] = T_f * (Ta_P - Ta_N)

with `T_f` the same frozen upwind face value the primal uses. Against the
continuous `Ta*grad(T)` that `addMomentumSource()` actually inserts:

    max|S exact face-based|      = 819363
    max|S continuous Ta*grad(T)| = 265729
    max|difference|              = 835500
    relative                     = 1.02          <-- 102 % wrong

    worst cell 1278:
      exact      = (-237586,  784162, 0)
      continuous = (   -645,  -17037, 0)

Wrong by a factor of ~3 in magnitude and, in the worst cells, **wrong in
direction** -- the y-component flips sign.

**This explains everything.** In a low-velocity Brinkman sponge (`cases/fdcheck`,
`cases/varprops`, `cases/demo2d`) `u . u_a` is tiny, so this error is invisible
and every FD gate passes. In an **open channel** with upwind convection it is
O(1), and `a'(beta) = betaMax * I'(beta)` at beta -> 0 (2500 * 21) then amplifies
it into the 10-40x sensitivity error. It also explains why raising the
diffusivity 10x shrank the error to 2.8-5.7x: weaker grad(T) -> smaller Ta ->
weaker source -> smaller u_a.

### The fix is NOT yet correct

A first `couplingForm exactFaceTranspose` is implemented (see
`addMomentumSource()`), mapping `phi_f = Sf_f . U_f` with
`U_f = w_f U_P + (1-w_f) U_N` back onto the cells. **It makes the FD gate
worse**, so it is not right:

| couplingForm | cell 628 | cell 626 |
|---|---|---|
| `TaGradT` (current default) | 1.04e-5 | 3.55e-5 |
| `exactFaceTranspose` | 1.79e-5 | 8.02e-5 |
| `exactFaceTranspose`, sign -1 | -2.61e-5 | -1.05e-4 |

Target is 1e-6. The likely reasons, in order:

1. **Rhie-Chow.** In SIMPLE the face flux is NOT `interpolate(U).Sf` -- it
   carries the Rhie-Chow pressure-smoothing correction. So `dphi_f/dU_P` is not
   `w_f Sf_f`, and the chain from face flux back to cell velocity is wrong. The
   framework's own adjoint machinery must already contain the correct
   `dphi/dU` linearisation (that is what makes the stock PtLosses adjoint work);
   the fix should reuse it rather than re-derive it.
2. Boundary faces are omitted from the source.
3. The `bounded` scheme's `Sp(div(phi))` term is omitted.

So: the DEFECT is proven and located, but the REPLACEMENT is still open. The
`couplingForm` switch (`TaGradT` | `negTGradTa` | `exactFaceTranspose` | `none`)
is the place to land it, and `utilities/testAdjointTranspose` is the gate it must
pass BEFORE any optimisation cycle is run.

### CONFIRMED: the exact flux sensitivity, including the bounded-convection term

`utilities/testAdjointTranspose` test (4) finite-differences the face Lagrangian
`L(phi) = Ta . (A(phi) T)` directly, with NO solve, and compares against the
closed form. Writing `A_P = C_P Ta_P` (the primal ROW-scales by C_P, so the jump
is `C_P Ta_P - C_N Ta_N`, NOT `C_f (Ta_P - Ta_N)`), and with `chi = 1` for a
`bounded` scheme:

    gPhi_f = T_f (A_P - A_N) - chi (A_P T_P - A_N T_N)

Measured:

       face      FD dL/dphi        gPhi(chi=1)        gPhi(chi=0)
       2388   -304.074666184   -304.074666183       6.34
       2390   -259.703889269   -259.703889270      90.73
       2392   -230.509321407   -230.509321409     112.55
       2547   -222.420199711   -222.420199709     -11.43

The FD matches `chi = 1` to **10-11 significant figures**; `chi = 0` is
completely wrong. So:

- The **bounded-convection correction is real and essential.** OpenFOAM's
  `bounded Gauss upwind` assembles `div(phi,T) - Sp(div(phi),T)`. At convergence
  `div(phi) ~ 0`, so the primal RESIDUAL is unchanged -- but the DERIVATIVE with
  respect to the flux is not: `d/dphi[-T div(phi)] = -T div(dphi) != 0`. Every
  continuous form (`Ta grad(T)`, `-T grad(Ta)`) misses this term entirely, which
  is why they all failed by the same O(1) factor.
- For bounded upwind, `gPhi` collapses to the DOWNWIND adjoint value times the
  temperature jump: `Ta_N (T_N - T_P)` for phi > 0, `Ta_P (T_N - T_P)` for
  phi < 0. A completely different object from `Ta grad(T)`.

`couplingForm exactFluxTranspose` implements exactly this and is committed.

### The one remaining gap: dphi/dU is NOT w_f * Sf (Rhie-Chow)

Even with the flux sensitivity proven correct, the FD gate still fails:

| couplingForm | cell 628 | cell 626 |
|---|---|---|
| `TaGradT` | 1.04e-5 | 3.55e-5 |
| `exactFluxTranspose` | 1.20e-5 | 3.83e-5 |
| target | 1e-6 | 1e-6 |

The flux sensitivity `gPhi_f` is right, so the error must now be in the MAPPING
from face flux back to cell velocity. `exactFluxTranspose` currently assumes

    phi_f = Sf_f . U_f,   U_f = w_f U_P + (1 - w_f) U_N
    =>  dphi_f/dU_P = w_f Sf_f

But in SIMPLE the face flux is **Rhie-Chow corrected**:

    phi = interpolate(HbyA) . Sf - rAUf * snGrad(p) * |Sf|

so `phi` is not a plain interpolation of U, and `dphi/dU != w_f Sf_f`. A flux
sensitivity therefore cannot be expressed correctly as a pure CELL MOMENTUM
source at all: in a SIMPLE/continuous-adjoint loop it has to enter through the
adjoint flux / adjoint pressure coupling, which is precisely the machinery that
makes the stock `PtLosses` adjoint work (and which passes the same FD gate at
1.22-1.36e-6).

**Next:** find where upstream injects a flux-level sensitivity for flow
objectives (the adjoint pressure/flux equation, not `addMomentumSource`) and add
`gPhi_f` there, rather than converting it to a cell source. Do NOT re-derive the
Rhie-Chow linearisation -- reuse the framework's.

### Next diagnostic

The remaining rigorous test is an **operator-only FD**, with no optimisation and
no primal/adjoint solve. Do it in TWO parts, scalar first:

**(1) Is `fvm::div(-phi,Ta)` the algebraic transpose of `fvm::div(phi,T)`?**
Freeze phi, T, Ta, C, mesh and schemes. For an internal face f (owner P,
neighbour N) with frozen primal face value `T_f = cP_f T_P + cN_f T_N` (freeze the
limiter; for upwind cP/cN are 1/0 by flux sign), the primal contributes
`R_P += phi_f T_f`, `R_N -= phi_f T_f`. The algebraic transpose is therefore

    (A^T Ta)_P += phi_f * cP_f * (Ta_P - Ta_N)
    (A^T Ta)_N += phi_f * cN_f * (Ta_P - Ta_N)

which is NOT generally `fvm::div(-phi,Ta)` under the same named scheme. If this
fails, u_a is wrong no matter how good addMomentumSource() is.

**(2) Then the momentum source.** Freeze T, Ta, U, phi, C and the mesh.
Define the functional corresponding to the primal thermal convection residual,

    M(U) = sum_cells Ta_P * R_conv,P(U, T)      // R_conv from fvm::div(phi(U), T)

and finite-difference M with respect to a velocity component in the suspect
open-fluid cells. Compare against the vector source that `addMomentumSource()`
actually inserts. The exact face form should be, conceptually,

    d/dphi_f [Ta^T R_T] = T_f * (Ta_P - Ta_N)
    S_U,P += w_f       * T_f * (Ta_P - Ta_N) * Sf_f / V_P
    S_U,N += (1 - w_f) * T_f * (Ta_P - Ta_N) * Sf_f / V_N

with the sign settled by the FD, not by inspection (the fvMatrix source-side
convention can flip what looks obvious on paper). For variable rho*cp, replace
`(Ta_P - Ta_N)` with `(C_P Ta_P - C_N Ta_N)`, because the primal row-scales by
`C.internalField()`, not by a face-interpolated `C_f`.

The production candidate is `couplingForm exactFaceTranspose`. `TaGradT` and
`negTGradTa` are DIAGNOSTIC ONLY: neither uses the same face interpolation,
upwind selection, limiter, flux orientation, boundary treatment or row-scaling
as the primal `fvm::div(phi,T)`, so neither can be expected to be its transpose. That settles unambiguously whether the continuous source is
the discrete transpose of `fvm::div(phi,T)`, and it fixes the sign convention
without guesswork. It is the exact analogue of `fd_chain.py`, but for

    U -> phi(U) -> fvm::div(phi,T) -> thermal residual -> adjoint momentum source

rather than the design-variable chain.

Companion test case: an **open-channel variant of `cases/fdcheck`** (initial
alpha 0 through the design box, so u.u_a is large and the Brinkman term carries
real weight). Every existing case puts the design in a low-velocity Brinkman
sponge, which is why none of them exercises this path.

### Technical debt found on the way (variable rho cp)

Not the current bug -- this case is constant-property -- but it must be on the
Q7 list before variable fluid properties are trusted in OPEN-flow optimisation.
The primal row-scales the convection matrix, `C.internalField()*fvm::div(phi,T)`,
while the adjoint uses a face-interpolated flux `C_f phi_f`. The exact transpose
of a cell-row-scaled convection matrix is not generally the face-interpolated
form: it should involve owner/neighbour cell factors, closer to

    C_P Ta_P - C_N Ta_N     rather than     C_f (Ta_P - Ta_N)

In `cases/varprops` the design sits in a sponge, so this is not exercised there
either.

### Status

Audit the **ATC-T coupling term** `couplingSign * C(T) * Ta * grad(T)` in
`thermalAdjointSimple::addMomentumSource()` against the derivation
(docs/derivation.md section 3.2), in a case where the Brinkman sensitivity is
NOT negligible. Concretely:

1. Build a variant of `cases/fdcheck` with an **open channel** through the
   design box (initial alpha 0 there, not 0.05-0.35) so that `u . u_a` is
   large and the Brinkman term carries real weight. Run the FD gate. If that
   fails too, the ATC-T term is confirmed and the co-optimisation geometry is
   incidental.
2. Check the residual convention: the framework's adjoint momentum equation may
   expect the source on the opposite side, or scaled by the volume, relative to
   what `addMomentumSource()` adds. A constant factor would show up as a
   constant ratio error; a `C(T)`- or `grad(T)`-dependent one would not.
3. Only then revisit `postProcessSens`.

Note the derivation note (section 9, "Open questions for reviewers") already
flags exactly this: *"Sign/normalisation of (ATC-T) against the NTUA residual
convention -- checked in code against FD, but an independent eye is wanted."*
The FD check referenced there could not have caught a magnitude error, because
every case in the repo at that time had a negligible Brinkman contribution.

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

---

## Solution roadmap (13 Jul 2026) — the missing object is the adjoint of the SIMPLE pressure projection

`gPhi_f` is proven (10-11 significant figures against an FD of the face
Lagrangian). What is missing is not a coefficient but an entire elimination step.

The primal flux that the energy equation actually convects with is the
**pressure-corrected** flux:

    phi = phiHbyA - pEqn.flux()

Writing `F = phiHbyA`, `C` for the pressure-flux correction operator
(`p -> pEqn.flux()`), `D` for the cell divergence and `L = D C` for the pressure
equation operator:

    dp   = L^-1 D dF
    dphi = dF - C dp

so, for a thermal flux sensitivity `g = gPhi = dL_T/dphi`,

    dJ = g^T dphi = [ g - D^T lambda ]^T dF,      L^T lambda = C^T g

**The correct predictor-flux sensitivity is therefore**

    gF = gPhi - D^T lambda

not `gPhi`. `exactFluxTranspose` currently maps `gPhi` straight onto the cells
with `dphi_f/dU_P ~ w_f Sf_f`, which is equivalent to assuming
`phi = interpolate(U) . Sf`. That injects the right face sensitivity into the
**wrong variable**, and the omitted term is the whole adjoint of the pressure
projection -- which is why `TaGradT`, `-T grad(Ta)` and `exactFluxTranspose` all
sit in the same 10-40x band rather than differing.

### Route A is NOT available (checked at source, 13 Jul)

The natural fix would be to hand `gPhi` to upstream and let the existing adjoint
flow machinery project it. **There is no such hook.** The complete set of
derivative hooks an `objectiveIncompressible` may supply is

    dJdv  dJdp  dJdT  dJdb  dJdbField  dJdnut  dJdGradU
    + the boundary variants (boundarydJdv, boundarydJdp, ...)

There is **no `dJdphi`, no flux-sensitivity concept anywhere in
`optimisation/adjointOptimisation/adjoint`.** `objectivePtLosses` never needs
one, because it is a *boundary* objective: it enters the adjoint through `dJdp`
and `dJdv` on the inlet/outlet patches and never touches the interior flux
operator. That is precisely why the stock `PtLosses` adjoint passes the same FD
gate (1.22-1.36e-6) in this geometry under a continuous adjoint, while a thermal
objective does not: **a thermal objective depends on velocity through the
interior convection operator, evaluated on a SIMPLE-projected flux.**

**This is a genuine capability gap in the framework, not a defect in it.** The
continuous-adjoint architecture has no mechanism to express an objective
sensitivity with respect to the discrete face flux -- which is exactly what any
objective convected by `phi` requires once the discrete and continuous adjoints
are asked to agree at O(1) accuracy in convection-dominated flow.

### Route B — eliminated projection (the implementation target)

    1. gPhi_f      from the proven chi=1 formula (done)
    2. rhs         = C^T gPhi          (face -> cell)
    3. solve       L^T lambda = rhs    (same operator as the pressure equation)
    4. gF_f        = gPhi_f - (lambda_P - lambda_N)
    5. Su_P       += w_f     gF_f Sf_f / V_P
       Su_N       += (1-w_f) gF_f Sf_f / V_N
    6. matrix     += couplingSign_ * Su

with `c_f` in step 2 taken from the coefficient `pEqn.flux()` actually uses
(`rAUf * magSf * nonOrthDeltaCoeffs`), **not** a guessed one. Do NOT combine this
with a Route-A-style uneliminated coupling: that would double-count the
projection.

### Acceptance criteria — both are no-solve identities, and both must pass BEFORE any topology FD run

**(i) Projection-adjoint identity.** Freeze `rAUf`, `F`, the mesh and the
pressure-equation coefficients. For random `dF`:

    L dp = D dF ;  dphi = dF - C dp
    check:  sum_f gPhi_f dphi_f  ==  sum_f gF_f dF_f        (to roundoff)

**(ii) HbyA mapping identity.** For random cell-vector `dH`:

    dF_f = Sf_f . (w_f dH_P + (1-w_f) dH_N)
    check:  sum_f gF_f dF_f  ==  sum_P V_P Su_P . dH_P      (to roundoff)

Only then is the open-channel FD gate meaningful again. `utilities/
testAdjointTranspose` is the natural home for both.

### Scope note

The scalar-adjoint transpose test stays in the suite regardless: for
constant-property bounded upwind `fvm::div(-phi,Ta)` **is** the exact transpose
(verified: off-diagonals match to 0), but that is a property of that scheme pair,
not a general guarantee.
