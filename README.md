# thermalTopO

**Conjugate-heat-transfer topology optimisation for OpenFOAM's adjoint framework.**

`thermalTopO` extends the topology-optimisation capability of
`adjointOptimisationFoam` (OpenFOAM v2512) with the thermal physics it
currently lacks: a temperature equation in the primal, the corresponding
adjoint temperature equation with its coupling into adjoint momentum,
thermal objective functions, and the conductivity-interpolation term in the
topology sensitivities. With it, the existing porosity-based topology
optimisation — including the upstream constrained update methods
(ISQP / nullSpace / MMA) and the fully differentiated k-ω SST adjoint —
can optimise cooling geometries for thermal objectives such as *minimise
peak wall temperature subject to a pressure-drop cap*.

Developed and released by [Fira Software Ltd](https://firasoftware.com)
(author: Dr S. Kalogerakos), July 2026. Motivated by heat-exchanger and
high-heat-flux cooling design problems, including the topology optimisation
of coolant channels in fusion plasma-facing components.

> This offering is not approved or endorsed by OpenCFD Limited, producer
> and distributor of the OpenFOAM software via www.openfoam.com, and owner
> of the OPENFOAM and OpenCFD trade marks.

## What it adds

| Component | Class | Status |
|---|---|---|
| Primal SIMPLE + energy equation over the design field | `thermalSimple` | working, verified |
| Adjoint SIMPLE + adjoint energy equation and Ta∇T momentum coupling | `thermalAdjointSimple` | working, verified |
| Temperature-dependent properties ρ(T), c_p(T), k_f(T), k_s(T) in primal, adjoint and sensitivities | `thermalPropertyTables` | working, FD-verified |
| Temperature-dependent viscosity μ(T) in the momentum equation | `viscosityModels::temperatureTable` | working, verified |
| Zone-mean temperature objective | `objectiveMeanTemperature` | working, FD-verified |
| Patch p-norm (peak-surrogate) temperature objective | `objectivePatchTemperaturePNorm` | working |
| Conductivity term in topology sensitivities | via `topOSensMultiplier` hook | working, FD-verified |
| Bergles–Rohsenow onset-of-nucleate-boiling monitor | `boilingOnsetBerglesRohsenow` functionObject | working |

Everything plugs into unmodified OpenFOAM v2512 through its runtime-selection
and sensitivity extension points; no core sources are patched.

## Verification (gradient correctness)

Continuous-adjoint gradients are verified against central finite
differences on a 2D laminar conjugate-heat-transfer case (Brinkman sponge,
conductivity ratio 100, hot wall, mean-temperature objective on a
downstream zone); every constituent solve converged to machine-level
residuals with automated guards:

![FD verification](docs/figures/fd_verification.png)

- **Single-cell central FD, solid interior, consistent 2nd-order schemes:
  zone-mean objective median |error| 0.5 % (max 2.1 %); patch p-norm
  objective median |error| 0.2 % (max 0.9 %).** Against the derivation
  note's pre-registered < 1 % laminar gate: the p-norm campaign passes it
  outright; the zone-mean campaign passes on median with max 2.1 %. The
  p-norm campaign exercises the boundary-driven adjoint flux path used by
  peak-temperature objectives.
- **Production configuration (regularisation on, p-norm objective,
  chained sensitivities): all regimes verify to max 1.7 % per cell**
  (solid 0.5 % median, edge 0.6 %, sponge 1.2 %) - the design-step error
  localisation observed without regularisation vanishes, as predicted.
  Note for reproducers: with regularisation active compare FD against
  `topOSens<solver>` (enable `writeAllFields`), not `topologySens<solver>`,
  which upstream writes before the filter/projection chain rule.
- Directional derivatives across 720 design cells: 1–3 % on
  gradient-dominant directions.
- With 1st-order upwind advection the solid-interior error grows to a
  uniform ≈ +11 % — an advection-scheme-consistency effect, not a
  formulation error. Cells adjacent to *unfiltered* step discontinuities
  in the design field carry larger local errors; regularisation (on in
  any production run) removes that configuration by construction.
- Known v1 approximation: frozen turbulent thermal diffusivity in the
  adjoint energy equation (the momentum-side k-ω SST adjoint is the fully
  differentiated upstream one). Quantification on turbulent cases is on
  the roadmap.

Full derivation and verification protocol: [docs/derivation.md](docs/derivation.md).

## Temperature-dependent fluid and solid properties

Properties can be tabulated in temperature **in both the fluid and the solid**,
and they stay active **throughout** the optimisation — primal, adjoint and
sensitivity assembly — not just in a final verification run.

| Property | Entry | Where it acts |
|---|---|---|
| ρ(T) [kg/m³] | `rhoTable` | volumetric heat capacity, and ν = μ/ρ |
| c_p(T) [J/kg/K] | `cpTable` | volumetric heat capacity |
| k_f(T) [W/m/K] | `kFluidTable` | fluid diffusivity D_f = k_f/(ρc_p)_ref |
| μ(T) [Pa s] | `muTable` | momentum, via ν(T) = μ(T)/ρ(T) |
| k_s(T), as D_s [m²/s] | `DSolidTable` | solid diffusivity |
| Pr_t | `Prt` | constant (default 0.85), documented, frozen in the adjoint |

The energy equation is solved as

    C(T) u·∇T = ∇·[D_eff(β, T) ∇T],    C(T) = ρ(T)c_p(T)/(ρc_p)_ref,
                                        D_f(T) = k_f(T)/(ρc_p)_ref

i.e. the ρc_p variation is carried on the **convective** term, scaled by a
*constant* reference (ρc_p)_ref. That is not cosmetic: dividing by the *local*
ρc_p instead — putting α(T) = k/(ρc_p) straight into the Laplacian — breaks
continuity of the physical heat flux k∇T at the fluid/solid interface, which
is the one thing a conjugate formulation may not get wrong. The local thermal
diffusivity α(T) = k/(ρc_p) = D_f/C is computed by the property class
(`thermalPropertyTables::alphaField`) for diagnostics, and is deliberately
never used to assemble the equation. All properties are re-evaluated from the
current T every primal iteration, on cells *and* boundary faces. Derivation:
[§6](docs/derivation.md).

μ(T) reaches the momentum equation through a new incompressible transport
model, `temperatureTable`, selected in `constant/transportProperties`:

```
transportModel  temperatureTable;
temperatureTableCoeffs
{
    TRef        300;    // nu used before T exists
    muTable     table ((300 1.0e-3) (340 5.0e-4));   // [Pa s]
    rhoTable    table ((300 1000)   (340 850));      // [kg/m3]
    // nuTable  table ((300 1e-6) (340 5.88e-7));    // [m2/s], alternative
}
```

and the thermal tables go in the `thermal` subdict of **both** solvers:

```
thermal
{
    TRef            300;    // (rho cp)_ref = rho(TRef)*cp(TRef)
    rhoTable        table ((300 1000) (340 850));
    cpTable         table ((300 4000) (340 3600));
    kFluidTable     table ((300 40)   (340 44));
    DSolidTable     table ((300 1.2e-3) (340 0.8e-3));
    Prt             0.85;
    kInterpolation  { function BorrvallPetersson; b 20; }
}
```

A table overrides its constant (`DFluid`, `DSolid`), so constant-property
cases are unaffected and reproduce their earlier results bitwise. `rhoTable`
appears in both dictionaries because OpenFOAM builds the transport model and
the primal solver from different files; **keep them consistent — they are not
cross-checked.** A variable-property case must also declare the adjoint's
C-weighted convective flux in `fvSchemes`, which must *not* be `bounded`:

```
div(-phiC,Ta)   Gauss limitedLinear 1;
```

**Adjoint treatment: frozen-property linearisation.** ∂C/∂T, ∂D/∂T and ∂ν/∂T
are not differentiated; C, D_eff and ν are frozen coefficient fields at the
primal solution. The cost of that is measured, not assumed.

**Verification** (`cases/varprops`, laminar, production configuration, all
five tables active at once: ρ, c_p, k_f, μ *and* D_s. Across the tabulated
300–340 K span, which the case's temperature field realises to 338.8 K:
ρc_p −23.5 %, ν −41 %, k_f +10 %, D_s −33 %):

| Regime | median \|err\| | max \|err\| | same case, constant properties |
|---|---|---|---|
| solid interior | 0.4 % | 0.8 % | 0.5 % / 1.0 % |
| interface edge | 0.7 % | 0.9 % | 0.6 % / 1.5 % |
| sponge | 1.4 % | 1.5 % | 1.2 % / 1.7 % |

Sign agreement 18/18 cells. Per-cell gradient accuracy with every property
tabulated is indistinguishable from the constant-property campaign. Across
the five directional derivatives the sign is likewise always correct, with
magnitude errors of 2.0–5.8 % against 1–3 % constant-property — a modest,
honestly-reported cost of the frozen-property terms. That μ(T) genuinely
reaches momentum was confirmed against an otherwise identical constant-ν run
(22 % change in peak pressure, 11.5 % in peak velocity).

The tables are the mechanism; the verification values are chosen for FD
conditioning. Production runs load real coolant/structural data through the
same entries (for a fusion monoblock: IAPWS water at 95 bar, tungsten k(T)).

## Build

Requires an OpenFOAM v2512 installation (openfoam.com line) with
development headers (`openfoam2512-default` Debian/Ubuntu package or
equivalent).

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
cd src && wmake
```

Produces `$FOAM_USER_LIBBIN/libthermalTopO.so`. Load it per case via
`libs ("libthermalTopO.so");` in `system/controlDict`.

## Usage sketch

In `system/optimisationDict`, select the thermal solver pair and add a
`thermal` sub-dictionary to both:

```
primalSolvers
{
    op1
    {
        type    incompressible;
        solver  thermalSimple;
        thermal
        {
            DFluid  1.72e-7;   // k_f / (rho_f cp_f)  [m2/s]
            DSolid  3.55e-5;   // k_s / (rho_f cp_f)  [m2/s]
            Prt     0.85;
            kInterpolation { function BorrvallPetersson; b 20; }
        }
        ...
    }
}
adjointManagers
{
    am1
    {
        primalSolver op1;
        adjointSolvers
        {
            as1
            {
                type    incompressible;
                solver  thermalAdjointSimple;
                thermal { /* same entries */ }
                objectives
                {
                    type incompressible;
                    objectiveNames
                    {
                        peakT
                        {
                            weight  -1e-6;   // small weight keeps adjoint
                                             // magnitudes solver-friendly
                            type    patchTemperaturePNorm;
                            patches (heatedWall);
                            P       12;
                        }
                    }
                }
                ...
            }
            // pressure-drop cap: standard adjointSimple + PtLosses with
            // isConstraint true; volume cap: null solver + topOVolume
        }
    }
}
```

Fields: `0/T` (primal temperature) and `0/Ta` (adjoint temperature,
dimensions `[0 2 -2 -1 0 0 0]`; fixedValue 0 at inlets *and* outlets,
fixedGradient elsewhere — the solver drives objective-patch gradients).
The Brinkman source is the upstream `topOSource` fvOption on `(U Ua)`;
initial designs are supplied through the `alpha` field.

Working examples under `cases/`: `fdcheck` (verification),
`demo2d` (heat-extraction maximisation under pressure-drop and volume
constraints). In the demo, the optimiser grows a sharp conductive fin on
the hot wall, improving the objective by +5.3 K while the nullSpace update
brings total-pressure losses to the 2x-baseline cap, strictly feasible to
within 0.03 % at cycle 80:

![fin growth demo](docs/figures/demo_fin_growth.png)

## Practical notes (learned the hard way)

- Solve `T` and `Ta` with tight inner tolerances (`relTol 0;`) — they are
  linear given the flow; loose inner solves turn convergence into a crawl.
- Do **not** under-relax `Ta`: even `relax(1.0)` clips for diagonal
  dominance against the current field and biases the fixed point. Leave
  `Ta` out of `relaxationFactors`.
- Near sharp Brinkman interfaces, prefer `bounded Gauss upwind` /
  `limitedLinear` over `cellLimited`-gradient `linearUpwind`: limiter
  flip-flop can stall SIMPLE at a constant residual plateau.
- Keep the thermal objective weight small (≈1e-6 for SI Kelvin-scale
  objectives): the adjoint fields of mean/peak-temperature objectives are
  physically enormous, and segregated adjoint loops are happier at small
  magnitudes. Scales cancel in the optimiser.
- The `nullSpace` update method has proven more robust than ISQP on
  badly scale-separated constraint sets in our tests.

## Licence

GPL-3.0-or-later (as required for OpenFOAM-linking code). Copyright (C)
2026 Fira Software Ltd.
