# Continuous adjoint for conjugate-heat-transfer topology optimisation

**thermalTopO derivation note — draft v0.1, 10 July 2026**
**For review: adjoint/optimisation formulation (O.D.), CHT physics (M.G.)**

This note derives the adjoint system implemented by `thermalTopO`, an
extension of OpenFOAM's `adjointOptimisationFoam` (v2512) adding an energy
equation, its adjoint, and thermal objectives to the existing porosity-based
topology-optimisation framework. Sign conventions follow the residual forms
used by the NTUA/PCOpt implementation so that the new terms drop into the
existing sensitivity assembly.

## 1. Primal problem

Steady, incompressible, turbulent flow with Brinkman penalisation over a
design field β ∈ [0, 1] (β = 0 fluid, β = 1 solid), single-domain
conjugate heat transfer. With kinematic pressure p and I(β) a RAMP
interpolation (I(0) = 0, I(1) = 1):

    R_p ≡ ∇·u = 0
    R_u ≡ (u·∇)u + ∇p − ∇·[2 ν_eff D(u)] + a(β) u = 0
    R_T ≡ u·∇T − ∇·[D_eff(β) ∇T] = 0

where

    D(u)      = ½(∇u + ∇uᵀ)
    ν_eff     = ν + ν_t                      (ν_t from primal k-ω SST)
    a(β)      = a_max I(β)                   (Brinkman inverse permeability)
    D_lam(β)  = D_f + (D_s − D_f) I_k(β),    D_f = k_f/(ρ_f c_p,f),
                                             D_s = k_s/(ρ_f c_p,f)
    D_eff(β)  = D_lam(β) + (1 − I_k(β)) ν_t/Pr_t

Notes:

1. **Single-field CHT scaling.** The energy equation is written per unit
   (ρ_f c_p,f). In solidified cells u → 0 (Brinkman), so the equation
   degenerates to ∇·(D_s ∇T) = 0, which at steady state is exactly the solid
   conduction equation — the ρc_p factor is immaterial without a transient
   or source term. Interface flux continuity holds because the *same*
   constant 1/(ρ_f c_p,f) multiplies the conductivity on both sides of any
   interior face; face diffusivities use harmonic interpolation so the
   fluid/solid conductivity jump (k_s/k_f ≈ 206 for tungsten/water) is
   respected.
2. **Interpolations.** a(β) uses the framework's existing Brinkman RAMP
   (`betaMax`, interpolation function objects). I_k(β) for conductivity is a
   separate RAMP with its own convexity parameter, because the optimal
   penalisation of conduction and of momentum differ (large k-contrast is the
   known hard case: continuation on the RAMP parameter is required; see §7).
3. **Turbulence masking.** ν_t/Pr_t is multiplied by (1 − I_k(β)) so
   solidified cells carry no turbulent diffusivity even where the turbulence
   model leaves residual ν_t.
4. **Volumetric heating / boundary flux.** The target problem class has a fixed
   surface heat flux q_pf on the heat-incident patch Γ_pf and adiabatic
   remaining boundaries; no volumetric source. BCs: T = T_in at inlet, ∇T·n = 0 at outlet and
   adiabatic walls, −k ∂T/∂n = −q_pf on Γ_pf (flux into domain).

## 2. Objective

Peak temperature is non-differentiable, so we minimise a p-norm surrogate on
the plasma-facing patch:

    J = ( 1/|Γ_pf| ∫_{Γ_pf} T^P dΓ )^{1/P},    P ≈ 8–20 (continuation on P)

For derivation we use the monotone transform J_P = ∫ T^P dΓ (identical
minimiser). A cellZone-mean variant `objectiveMeanTemperature`
(J = ∫_Z T dΩ / |Z|) is also implemented — it is the simplest correct
objective and is used for finite-difference verification.

## 3. Adjoint system

Lagrangian L = J_P + ∫_Ω (u_a·R_u + p_a R_p + T_a R_T) dΩ.

### 3.1 Variation in T → adjoint energy equation

Using ∇·u = 0 and integration by parts:

    δ_T L = ∫_Ω δT [ −u·∇T_a − ∇·(D_eff ∇T_a) ] dΩ
          + ∮_Γ δT [ T_a (u·n) + D_eff ∂T_a/∂n ] dΓ
          − ∮_Γ D_eff T_a ∂(δT)/∂n dΓ
          + ∫_{Γ_pf} P T^{P−1} δT dΓ

Setting each independent variation to zero:

**Field equation**

    −u·∇T_a − ∇·(D_eff ∇T_a) = 0                                   (AT)

(reversed convection; no volumetric source for boundary objectives).

**Adjoint temperature BCs**

| Primal patch | Primal T BC | Adjoint T_a BC |
|---|---|---|
| inlet | T = T_in (δT = 0) | T_a = 0 |
| outlet | ∇T·n = 0 | T_a = 0 (Robin relation degenerates to this for Pe >> 1; error O(D/uL)) |
| adiabatic walls (u = 0) | ∂T/∂n = 0 | ∂T_a/∂n = 0 |
| Γ_pf (fixed flux, u = 0) | −k ∂T/∂n = −q_pf | D_eff ∂T_a/∂n = −P T^{P−1} |

The objective enters solely through the flux-type BC on Γ_pf. For
`objectiveMeanTemperature` the Γ_pf BC becomes homogeneous and a volume
source −1/|Z| appears in (AT) instead.

### 3.2 Variation in u → thermal coupling into adjoint momentum

δ_u of ∫ T_a R_T dΩ contributes ∫ T_a (δu·∇T) dΩ, i.e. the adjoint momentum
equation of the existing framework gains the source term

    + T_a ∇T                                                        (ATC-T)

on its right-hand side (same side as the existing adjoint transpose
convection term, in the framework's residual convention). All other adjoint
momentum/pressure terms, BCs, and the differentiated k-ω SST adjoint are
untouched — they are inherited from `adjointOptimisationFoam`.

### 3.3 Frozen terms in v1 (documented approximations)

1. **∂ν_t/∂(u, β) in the energy equation**: the dependence of the turbulent
   diffusivity in R_T on the flow (through ν_t) is neglected in the adjoint
   ("frozen α_t"). The momentum-side SST adjoint remains fully differentiated
   (upstream capability). Consequence: gradients are inexact in strongly
   turbulence-dominated heat transfer; the FD check quantifies this (§5).
2. **Temperature-dependent k_s(T)** (implemented 11 Jul 2026): both solvers
   accept an optional `DSolidTable` (Function1); the primal evaluates
   D_s(T) per cell each iteration, the adjoint and the conductivity
   sensitivity term use the same field at the frozen primal T. The
   ∂D_s/∂T linearisation is NOT differentiated (frozen). Quantified with
   the table active in primal, adjoint AND sensitivities (an earlier run
   had it in the primal only and is superseded): with an exaggerated -33%
   D_s variation across the case's temperature range, per-cell FD errors
   equal the constant-property production campaign (solid 0.5% median /
   1.0% max; edge 0.6%/1.5%; sponge 1.2%/1.7%) - the frozen term
   contributes no measurable gradient error at these conditions.
   **Extended 12 Jul 2026 to the fluid side** — ρ(T), c_p(T), k_f(T), μ(T)
   as well as k_s(T), under the same frozen-property linearisation. The
   formulation, its interface-continuity constraint and the verification
   are in §6, which supersedes this note.

## 4. Sensitivities

With the adjoint fields available, the gradient of J w.r.t. the design
variable in each cell is the partial derivative of L through the explicit β
dependence:

    dJ/dβ_c = [ a'(β) (u·u_a)
              − I_k'(β) (D_s − D_f − ν_t/Pr_t) (∇T·∇T_a) ] V_c   (signs in
    the framework's convention; the diffusive term enters via integration by
    parts of −∇·(D' ∇T) against T_a)

The first term is the framework's existing Brinkman sensitivity; the second
is new and is added through the `designVariables`/sensitivity hooks
introduced by the v2506 re-architecture (no core-code modification).

## 5. Verification protocol (V1 gate)

Coarse 2D heated channel (~2k cells), laminar first, then SST:

1. Solve primal + adjoint; assemble dJ/dβ for all design cells.
2. For each of ~50 randomly selected design cells c: perturb β_c by ±ε
   (central differences, ε sweep 1e-3..1e-5), re-solve primal, evaluate J.
3. Report max relative deviation |FD − adjoint|/|FD| over the sample; accept
   < 1% laminar (exact adjoint), < 5–10% SST (frozen-α_t effect,
   objective-dependent). Plot published in the repository README.
   Compliance against this pre-registered gate, stated precisely: the
   p-norm campaign passes it outright (max 0.9%); the zone-mean campaign
   passes on median (0.5%) with max 2.1%, marginally outside the strict
   max-1% reading; the production configuration's max is 1.0% (solid) to
   1.7% (sponge). We report the gate as met in median terms and near-met
   in max terms, rather than claiming strict max-compliance everywhere.

### 5.1 Results (11 July 2026) — gate PASSED

Verification case: 2D channel 0.1 x 0.02 m, laminar (Re = 400), Dirichlet
hot wall, Brinkman sponge (alpha 0.05 background / 0.35 blob,
betaMax 2500), conductivity ratio 100, `objectiveMeanTemperature` on a
downstream cellZone. All runs converged to machine-level residuals with
automated guards (no unconverged samples admitted).

1. **Directional derivatives** (7 directions, 720 design cells): adjoint
   vs central FD agree to 1-3% on the dominant directions after fitting
   the single expected bookkeeping constant; near-gradient-orthogonal
   random directions are noise-dominated, as expected.
2. **Single-cell central FD** (20 cells across blob/edge/sponge regimes,
   first-order advection): 100% sign agreement; median |error| 12.3%,
   uniform +11% bias in solid-interior cells.
3. **Scheme attribution**: with second-order (limitedLinear) advection on
   T and T_a, the solid-interior per-cell errors collapse to
   **median 0.5%, max 2.1%** - the 1% laminar target met in median terms
   (see the compliance statement under the protocol above; the p-norm
   campaign of item 4 meets it outright, max 0.9%). The first-order bias
   is therefore advection-scheme consistency, not a formulation error.
4. **Patch p-norm objective (11 Jul 2026)**: the same single-cell protocol
   applied to `objectivePatchTemperaturePNorm` (P = 8, adiabatic wall
   patch), which exercises the boundary-driven adjoint flux path
   (update_boundarydJdT -> fixedGradient drive): solid interior
   **median 0.2%, max 0.9%**; high-signal sponge cells 2-3%; blob-edge
   cells 6-10% with consistent sign (unregularised design steps, same
   class as item 3's localisation note). Both objective classes are
   FD-verified.
5. **Production configuration (11 Jul 2026)**: regularisation ON (Helmholtz
   filter + tanh projection) with the p-norm objective, FD compared against
   the projection-chained sensitivities: **solid interior median 0.5% (max
   1.0%), edge median 0.6% (max 1.5%), sponge median 1.2% (max 1.7%)** -
   the design-step localisation of items 3-4 vanishes under regularisation,
   as predicted. No weak regime remains.
   **Harness caveat for reproducers**: with regularisation active, the
   per-cycle `topologySens<solver>` field is written BEFORE the
   filter/projection chain rule (upstream behaviour, documented in
   topODesignVariables.C); FD comparisons in design-variable space must
   enable `writeAllFields true` and use `topOSens<solver>` instead. The
   optimiser itself always receives the chained gradients.
4. Confirmed conventions: couplingSign = +1, thermalSensScale = +1
   (as derived); design-variable-to-indicator map verified as identity
   with regularisation off; sensitivity scale constant uniform across
   cells and attributable to the framework's volume/step bookkeeping.

Practical guidance derived from the campaign: solve T and T_a with tight
inner tolerances (they are linear given the flow); do not under-relax T_a
(fvMatrix::relax at factor 1.0 still clips and biases); prefer bounded
upwind or limitedLinear over cellLimited-gradient linearUpwind near sharp
Brinkman interfaces (limiter flip-flop stalls SIMPLE); primal outlet
T_a = 0 (high-Peclet limit of the Robin condition).

## 6. Temperature-dependent material properties

Implemented 12 July 2026, in **both** the fluid and the solid, active
throughout the optimisation (primal, adjoint and sensitivity assembly).
Supersedes the solid-only treatment of §3.3.2, which it subsumes.

### 6.1 Primal formulation

With variable properties the steady energy equation is

    ρ(T) c_p(T) u·∇T = ∇·[k(T) ∇T]

Divide through by a **constant** reference volumetric heat capacity
(ρ c_p)_ref = ρ(T_ref) c_p(T_ref):

    R_T ≡ C(T) u·∇T − ∇·[D_eff(β, T) ∇T] = 0

    C(T)        = ρ(T) c_p(T) / (ρ c_p)_ref                       [−]
    D_f(T)      = k_f(T) / (ρ c_p)_ref                            [m²/s]
    D_s(T)      = k_s(T) / (ρ c_p)_ref     (tabulated directly)   [m²/s]
    D_eff(β,T)  = D_f + (D_s − D_f) I_k(β) + (1 − I_k(β)) ν_t/Pr_t

With constant properties C ≡ 1 and this collapses **exactly** to §1; the
constant-property code path is unchanged and reproduces its earlier results
bitwise.

**Why the reference must be a constant, not the local ρc_p.** §1 note 1
guarantees interface flux continuity because the *same constant* 1/(ρ_f c_p,f)
multiplies the conductivity on both sides of every interior face: continuity
of D∇T is then equivalent to continuity of the physical flux k∇T. That
argument survives variable properties **only if the divisor stays constant**.
Dividing instead by the local field ρ(T)c_p(T) would put the local thermal
diffusivity α(T) = k/(ρ c_p) in the Laplacian, and

    ∇·(α ∇T) ≠ (1/ρc_p) ∇·(k ∇T)    whenever ρc_p varies in space,

the difference being a spurious flux α ∇ln(ρc_p)·∇T. It would also be
meaningless across the fluid/solid interface, where the *fluid's* ρc_p has no
definition on the solid side. So the ρc_p variation is carried on the
convective term, where it belongs, and never inside the divergence. α(T) =
D_f/C is exposed as a diagnostic output only; it is never used to assemble
the equation.

**Momentum.** ν(T) = μ(T)/ρ(T) enters through a new incompressible transport
model, `viscosityModels::temperatureTable`, which looks T up from the object
registry and rebuilds ν from the μ and ρ tables. `adjointOptimisationFoam`'s
SIMPLE loop already calls `laminarTransport().correct()` once per iteration,
so ν is refreshed every primal iteration with no change to the upstream loop.
ν_eff = ν(T) + ν_t as before. (This is the registry-lookup pattern of the
upstream `viscosityModels::Arrhenius` model, but ν is rebuilt from the tables
each call rather than accumulated onto its previous value.)

**Discretisation.** C multiplies the convection matrix row-wise (diagonal,
off-diagonals, source and boundary coefficients), which is exactly
C_P (u·∇T)_P V_P because the discrete face fluxes sum to zero. All properties
are re-evaluated from the current T at every primal iteration, on cells *and*
on boundary faces — the latter matters wherever solid meets a domain boundary
(a heat-incident surface on a fixed solid crust), since a face diffusivity left at a constant
would corrupt the wall heat flux.

### 6.2 Adjoint: frozen-property linearisation

C, D_eff and ν are treated as **frozen coefficient fields** evaluated at the
converged primal T: the property derivatives ∂C/∂T, ∂D/∂T and ∂ν/∂T are not
differentiated. Under that linearisation, δ_T of ∫ T_a R_T gives

    −∇·(C u T_a) − ∇·(D_eff ∇T_a) = − Σ_k w_k dJ/dT_k

and δ_u gives the modified thermal coupling into adjoint momentum

    + C(T) T_a ∇T                                              (ATC-T)

Both reduce to §3.1 and §3.2 when C ≡ 1.

Note the C-weighted adjoint flux is **not** solenoidal: ∇·(C u) ≠ 0 where C
varies. The adjoint convection must therefore use a *conservative*
(unbounded) scheme. A "bounded" scheme subtracts Sp(∇·(C u), T_a), which
would silently replace the conservative operator by the non-conservative
form C u·∇T_a — not the transpose of the primal term, and a wrong gradient.
The solver gives this flux its own scheme key, `div(-phiC,Ta)`, so a
variable-property case must declare it explicitly rather than inherit a
`bounded` entry by accident.

Sensitivities (§4) are unchanged in form, with the properties evaluated at
the frozen T:

    dJ/dβ_c ⊃ − I_k'(β) [D_s(T) − D_f(T) − ν_t/Pr_t] (∇T·∇T_a) V_c

C carries no β dependence: it multiplies convection, which the Brinkman term
drives to zero in solidified cells.

Omitted, and declared: ∂ν/∂T (the route by which adjoint momentum would feed
back into T_a through viscosity), ∂C/∂T, ∂D/∂T, and — as in §3.3.1 — ν_t
itself. The cost of these omissions is measured, not assumed (§6.3).

### 6.3 Verification (12 July 2026)

`cases/varprops`: the fdcheck geometry, laminar, production configuration
(regularisation on, patch p-norm objective), with **all five tables active
simultaneously** — ρ(T), c_p(T), k_f(T), μ(T) in the fluid and D_s(T) in the
solid. A manufactured verification fluid, anchored so that at 300 K it
reduces exactly to the constant-property case (ν = 1e-6, D_f = 1e-5) and the
Péclet/Reynolds regime is unchanged. Variation across the tabulated 300–340 K
span (the case's temperature field reaches 338.8 K, so it realises very
nearly all of it): ρ −15 %, c_p −10 %, ρc_p −23.5 % (so C: 1.000 → 0.765),
k_f +10 %, μ −50 % (so ν: 1.0e-6 → 5.9e-7, Re 400 → 680), D_s −33 %.

That μ(T) genuinely reaches the momentum equation was checked directly
against an otherwise identical constant-ν run: 22 % change in peak pressure,
11.5 % in peak velocity.

Single-cell central differences, ε = 0.008, compared against `topOSensas1`:

| Regime | median \|err\| | max \|err\| | constant-property campaign |
|---|---|---|---|
| solid interior | 0.4 % | 0.8 % | 0.5 % / 1.0 % |
| interface edge | 0.7 % | 0.9 % | 0.6 % / 1.5 % |
| sponge | 1.4 % | 1.6 % | 1.2 % / 1.7 % |

Sign agreement: **18/18** cells. Per-cell accuracy with every property
tabulated is statistically indistinguishable from the constant-property
production campaign — the frozen-property terms contribute no measurable
per-cell gradient error at these conditions.

Directional derivatives over the 720 design cells (blob, band, three random)
agree in sign in all five directions, with magnitude errors of 2.0–5.8 %
against 1–3 % for the same directions with constant properties. This modest
growth is the honest cost of the frozen-property linearisation and is
reported as such; it is well inside what a descent direction needs, and the
per-cell gradients that drive the filtered update are unaffected.

Reproduce: `cases/varprops/Allrun` then `./fd_varprops.py`.

### 6.4 Turbulent thermal diffusivity

Pr_t is a documented constant (default 0.85), with ν_t/Pr_t masked by
(1 − I_k(β)) as in §1 note 3, and frozen in the adjoint as in §3.3.1. The
verification case is laminar (ν_t = 0), so this campaign does not exercise
it. A tabulated or locally-varying Pr_t is a drop-in extension of the same
Function1 machinery; quantification on turbulent cases remains roadmap, as
for the frozen-ν_t term itself.

### 6.5 Production data

The tables are the mechanism; the values above are chosen for verification
conditioning. Production runs load real coolant and structural data through
the same entries: for a water-cooled tungsten component, IAPWS water at the
operating pressure and a tabulated tungsten conductivity. Note ρ and μ appear in *two* dictionaries
(`constant/transportProperties` for the viscosity model,
`optimisationDict`'s `thermal` subdict for the energy equation) because
OpenFOAM constructs the transport model and the primal solver from different
files; they are not cross-checked, and must be kept consistent by hand.

## 7. Constraints and optimiser

- **Pressure-drop cap**: existing `objectivePtLosses` as an inequality
  constraint under the framework's ISQP/MMA update methods (upstream since
  v2312). No new code.
- **No boiling (Bergles–Rohsenow ONB)**: v1 enforces a proxy constraint — a
  p-norm of T over near-interface fluid cells kept below T_ONB estimated at
  the local wall flux — plus an exact post-check via the
  `boilingOnsetBerglesRohsenow` functionObject on the body-fitted verified
  design. The exact ONB aggregate as a differentiable constraint is roadmap.
- **No enclosed voids / small-scale porosity**: enforced by the framework's
  regularisation (Helmholtz filter radius) + projection continuation, and by
  the extraction pipeline: β iso-surface → STL → snappyHexMesh →
  `chtMultiRegionSimpleFoam` (k-ω SST, wall-resolved) recompute of all
  reported quantities. Only body-fitted numbers are reported as final.
- **Continuation**: RAMP convexity for I_k(β) relaxed → sharpened over
  optimisation cycles (the 206:1 conductivity contrast makes early
  intermediate densities essential for a usable design space).

## 8. Mapping to a high-heat-flux component problem

The formulation targets the following generic problem class: a coolant channel
inside a conducting block, with a fixed heat flux on one surface, a protected
solid layer beneath it, and caps on pumping cost and on boiling margin.

| Problem element | Formulation element |
|---|---|
| Design domain: channel interior only | topOZones: design cells = fluid region; frozen fluid elsewhere |
| Design domain: channel + surrounding body, minus a protected surface crust | design cells extended into the solid body; crust cells fixed solid (`fixedPorousZones`) |
| Objective: peak temperature on the heat-incident surface | `objectivePatchTemperaturePNorm` on Γ_pf |
| Pressure-drop cap | `objectivePtLosses` as an ISQP/nullSpace constraint |
| No nucleate boiling | proxy constraint + ONB post-check (§7) |
| Iterative geometry history | framework's per-cycle β / iso-surface writes (upstream) |
| Temperature-dependent solid AND fluid properties | tabulated in primal, adjoint and sensitivities (§6) |

## 9. Open questions for reviewers

1. Sign/normalisation of (ATC-T) against the NTUA residual convention —
   checked in code against FD, but an independent eye is wanted.
2. Robin outlet BC for T_a: implementation as mixed BC with f = u·n/D_eff —
   any stability concerns at recirculating outlets? (Outlet extension makes
   recirculation unlikely here.)
3. Is the (1 − I_k(β)) masking of ν_t/Pr_t the right choice vs masking ν_t
   inside the turbulence model itself? (v1 masks diffusivity only; the
   Brinkman term already kills u in solid cells.)
4. P-continuation schedule for the p-norm (start P = 4, double to 16?) vs
   fixed P with normalisation.
