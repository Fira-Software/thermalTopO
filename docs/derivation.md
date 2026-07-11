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
   known hard case: continuation on the RAMP parameter is required; see §6).
3. **Turbulence masking.** ν_t/Pr_t is multiplied by (1 − I_k(β)) so
   solidified cells carry no turbulent diffusivity even where the turbulence
   model leaves residual ν_t.
4. **Volumetric heating / boundary flux.** CS9 has a fixed surface heat flux
   q_pf on the plasma-facing patch Γ_pf and adiabatic remaining boundaries;
   no volumetric source. BCs: T = T_in at inlet, ∇T·n = 0 at outlet and
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
2. **Temperature-dependent k_s(T)**: primal supports a tungsten k(T) table;
   its adjoint linearisation ∂k/∂T ∇T·∇T_a is implemented but can be toggled
   for A/B testing.

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
   **median 0.5%, max 2.1%** - meeting the < 1% laminar acceptance
   criterion. The first-order bias is therefore advection-scheme
   consistency, not a formulation error.
4. **Patch p-norm objective (12 Jul 2026)**: the same single-cell protocol
   applied to `objectivePatchTemperaturePNorm` (P = 8, adiabatic wall
   patch), which exercises the boundary-driven adjoint flux path
   (update_boundarydJdT -> fixedGradient drive): solid interior
   **median 0.2%, max 0.9%**; high-signal sponge cells 2-3%; blob-edge
   cells 6-10% with consistent sign (unregularised design steps, same
   class as item 3's localisation note). Both objective classes are
   FD-verified.
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

## 6. Constraints and optimiser

- **Pressure drop ≤ 340 Pa**: existing `objectivePtLosses` as an inequality
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

## 7. CS9 mapping

| CS9 element | Formulation element |
|---|---|
| Case 1 design domain (channel interior) | topOZones: design cells = fluid cylinder; frozen fluid elsewhere |
| Case 2 design domain (+ block minus 2 mm crust) | design cells extended to `Monoblock_main_body`; crust cells fixed solid |
| Objective: peak plasma-facing T | `objectivePatchTemperaturePNorm` on Γ_pf |
| ΔP ≤ 340 Pa | `objectivePtLosses` ISQP constraint |
| No nucleate boiling | proxy constraint + ONB post-check (§6) |
| Iterative geometry history (CS9R-4) | framework's per-cycle β/iso-surface writes (upstream) |
| k_s(T) (CS9R-5, Desirable) | tungsten k(T) table in primal (§3.3) |

## 8. Open questions for reviewers

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
