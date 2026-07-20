# thermalTopO

**Conjugate-heat-transfer topology optimisation for OpenFOAM's adjoint framework.**

`thermalTopO` extends the topology-optimisation capability of
`adjointOptimisationFoam` (OpenFOAM v2512) with the thermal physics it
currently lacks: a temperature equation in the primal, the corresponding
adjoint temperature equation with its coupling into adjoint momentum,
thermal objective functions, and the conductivity-interpolation term in the
topology sensitivities. With it, the existing porosity-based topology
optimisation, including the upstream constrained update methods
(ISQP / nullSpace / MMA) and the momentum-side adjoint k-ω SST implementation,
can optimise cooling geometries for thermal objectives such as *minimise peak
wall temperature subject to a pressure-drop cap*.

Developed and released by [Fira Software Ltd](https://firasoftware.com)
(author: Dr S. Kalogerakos), July 2026. Motivated by heat-exchanger and
high-heat-flux cooling design problems, including the topology optimisation
of coolant channels in fusion plasma-facing components.

> This offering is not approved or endorsed by OpenCFD Limited, producer
> and distributor of the OpenFOAM software via www.openfoam.com, and owner
> of the OPENFOAM and OpenCFD trade marks.

## What it adds

<table>
  <thead>
    <tr>
      <th>Component</th>
      <th>Class</th>
      <th>Status</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Primal SIMPLE + energy equation over the design field</td>
      <td><code>thermalSimple</code></td>
      <td>working, verified</td>
    </tr>
    <tr>
      <td>Adjoint SIMPLE + adjoint energy equation and thermal-to-momentum coupling</td>
      <td><code>thermalAdjointSimple</code></td>
      <td>working and FD-verified in published porous verification cases</td>
    </tr>
    <tr>
      <td>Exact matrix-free fixed-point beta sensitivity for the guarded Stokes branch</td>
      <td><code>useAnalyticFixedPointAdjointBeta</code></td>
      <td>production-integrated and runtime-guarded in v0.4.2</td>
    </tr>
    <tr>
      <td>Temperature-dependent properties ρ(T), c_p(T), k_f(T), k_s(T) in primal, adjoint and sensitivities</td>
      <td><code>thermalPropertyTables</code></td>
      <td>working, FD-verified</td>
    </tr>
    <tr>
      <td>Several fixed solid materials, labelled by cellZone, each with its own D_s or D_s(T)</td>
      <td><code>thermalPropertyTables</code></td>
      <td>working, verified</td>
    </tr>
    <tr>
      <td>Temperature-dependent viscosity μ(T) in the momentum equation</td>
      <td><code>viscosityModels::temperatureTable</code></td>
      <td>working, verified</td>
    </tr>
    <tr>
      <td>Zone-mean temperature objective</td>
      <td><code>objectiveMeanTemperature</code></td>
      <td>working, FD-verified</td>
    </tr>
    <tr>
      <td>Patch p-norm (peak-surrogate) temperature objective</td>
      <td><code>objectivePatchTemperaturePNorm</code></td>
      <td>working</td>
    </tr>
    <tr>
      <td>Conductivity term in topology sensitivities</td>
      <td>via <code>topOSensMultiplier</code> hook</td>
      <td>working, FD-verified</td>
    </tr>
    <tr>
      <td>Bergles-Rohsenow onset-of-nucleate-boiling monitor</td>
      <td><code>boilingOnsetBerglesRohsenow</code> functionObject</td>
      <td>working</td>
    </tr>
    <tr>
      <td>Discrete adjoint-transpose gate (no-solve operator test)</td>
      <td><code>utilities/testAdjointTranspose</code></td>
      <td>working</td>
    </tr>
    <tr>
      <td>Property-field reporting (<code>DSolid</code>, <code>DFluid</code>, <code>materialID</code>)</td>
      <td><code>utilities/writeThermalProperties</code></td>
      <td>working</td>
    </tr>
  </tbody>
</table>

Everything plugs into unmodified OpenFOAM v2512 through its runtime-selection
and sensitivity extension points; no core sources are patched.

## Status

**Verified** (reproducible from a clean checkout, in `cases/`):

- Conjugate-heat-transfer topology optimisation over a Brinkman design field.
- Thermal objectives: zone mean and patch p-norm (peak-temperature surrogate).
- Direct thermal conductivity topology sensitivities.
- Temperature-dependent property tables, fluid **and** solid, active in primal,
  adjoint and sensitivity assembly.
- Fixed-zone design controls (pinned solid, pinned fluid, explicit design space).
- Several fixed solid materials in one domain, each with its own constant or
  temperature-dependent diffusivity, active throughout the optimisation
  (`examples/multiSolidMaterials`). The design variable chooses fluid vs solid,
  not which solid.
- Geometry and design-history output, every cycle.
- Constrained demonstration: fin growth under a pressure-drop cap.
- Onset-of-nucleate-boiling monitor.
- Production fixed-point beta sensitivity for the guarded constant-property,
  one-way-coupled, serial Stokes branch. It solves
  `(I - M_x)^T psi = J_x`, evaluates `J_beta + psi^T M_beta`, replaces the
  stock Brinkman sensitivity without double-counting, and retains the direct
  conductivity contribution.

**Experimental** (in `examples/`, published for transparency and reuse):

- Open-channel channel-and-body co-optimisation. This exercises the
  thermal-induced Brinkman path in open-flow design cells. The exact
  fixed-point beta mode above is a guarded constant-property Stokes production
  path and a verification reference; it is not presented as an exact
  variable-property SST adjoint for this example. See
  [`docs/atc-t-open-channel.md`](docs/atc-t-open-channel.md).

## Exact fixed-point beta mode, v0.4.2

`thermalAdjointSimple` includes an opt-in production mode:

```
solveFixedPointAdjoint              true;
computeAnalyticFixedPointAdjointBeta true;
useAnalyticFixedPointAdjointBeta    true;
```

In that mode the solver uses the validated one-step SIMPLE transpose as a
matrix-free action, solves the fixed-point adjoint equation

    (I - M_x)^T psi = J_x

and assembles the topology sensitivity as

    dJ/dbeta = J_beta + psi^T M_beta.

The production implementation replaces the stock Brinkman momentum sensitivity
with `M_beta^T psi`, then adds the existing direct conductivity contribution
unchanged. It is deliberately guarded to the scope that has been validated:
OpenFOAM v2512, serial execution, non-coupled patches, Stokes flow, consistent
SIMPLE, zero pressure non-orthogonal correctors, constant thermal properties,
one-way flow-to-thermal coupling, one active `topOSource` momentum option, no
active temperature `fvOptions`, and zero explicit thermal non-orthogonal
diffusion correction. A structural lower-bound relaxation tie at `beta=0` is
handled as the feasible right derivative, not as a two-sided classical
derivative.

<table>
  <thead>
    <tr>
      <th>Directory</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>src/</code></td>
      <td>the module</td>
    </tr>
    <tr>
      <td><code>cases/</code></td>
      <td>verified cases; every gradient claim in this README is reproduced here</td>
    </tr>
    <tr>
      <td><code>utilities/</code></td>
      <td>verification tooling, incl. a no-solve operator-transpose gate</td>
    </tr>
    <tr>
      <td><code>examples/</code></td>
      <td>experimental work in progress; not evidence</td>
    </tr>
  </tbody>
</table>

Verification is treated as a gate, not a report: an operator is checked
algebraically before it is trusted, and a case is not cited until its gradients
pass finite differences.

## Verification (gradient correctness)

Continuous-adjoint gradients are verified against central finite
differences on a laminar conjugate-heat-transfer case with a Brinkman sponge,
high conductivity contrast, a hot wall and a downstream mean-temperature
objective. Every constituent solve converged to machine-level
residuals with automated guards:

![FD verification](docs/figures/fd_verification.png)

One global scale factor is applied consistently to account for the objective
weight/sign convention; agreement is assessed cell-by-cell after that single
factor.

- **Single-cell central FD, solid interior, consistent 2nd-order schemes:
  zone-mean and patch p-norm objectives verify at sub-percent median
  accuracy.** The patch p-norm campaign exercises the boundary-driven adjoint
  flux path used by peak-temperature objectives.
- **Production configuration, regularisation on, p-norm objective and chained
  sensitivities: all tested regimes verify at low single-digit percent or
  better.** The design-step error localisation observed without regularisation
  vanishes, as predicted.
  Note for reproducers: with regularisation active compare FD against
  `topOSens<solver>` (enable `writeAllFields`), not `topologySens<solver>`,
  which upstream writes before the filter/projection chain rule.
- Directional derivatives across the full design field retain the correct sign
  and low single-digit magnitude error on gradient-dominant directions.
- With first-order upwind advection the verification shows a uniform bias. This
  is an advection-scheme-consistency effect, not a formulation error. Cells
  adjacent to *unfiltered* step discontinuities in the design field carry larger
  local errors; regularisation, on in production runs, removes that
  configuration by construction.
- Known variable-property/SST approximation: turbulent thermal diffusivity,
  temperature-dependent coefficient fields and property derivatives are frozen
  at the converged primal state in the thermal adjoint. The exact fixed-point
  beta mode is a separate guarded constant-property Stokes path.

Full derivation and verification protocol: [docs/derivation.md](docs/derivation.md).

## Temperature-dependent fluid and solid properties

Properties can be tabulated in temperature **in both the fluid and the solid**,
and they stay active **throughout** the optimisation: primal, adjoint and
sensitivity assembly, not only a final verification run.

<table>
  <thead>
    <tr>
      <th>Property</th>
      <th>Entry</th>
      <th>Where it acts</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>ρ(T) [kg/m³]</td>
      <td><code>rhoTable</code></td>
      <td>volumetric heat capacity, and ν = μ/ρ</td>
    </tr>
    <tr>
      <td>c_p(T) [J/kg/K]</td>
      <td><code>cpTable</code></td>
      <td>volumetric heat capacity</td>
    </tr>
    <tr>
      <td>k_f(T) [W/m/K]</td>
      <td><code>kFluidTable</code></td>
      <td>fluid diffusivity D_f = k_f/(ρc_p)_ref</td>
    </tr>
    <tr>
      <td>μ(T) [Pa s]</td>
      <td><code>muTable</code></td>
      <td>momentum, via ν(T) = μ(T)/ρ(T)</td>
    </tr>
    <tr>
      <td>k_s(T), as D_s [m²/s]</td>
      <td><code>DSolidTable</code></td>
      <td>solid diffusivity</td>
    </tr>
    <tr>
      <td>Pr_t</td>
      <td><code>Prt</code></td>
      <td>constant (default 0.85), documented, frozen in the adjoint</td>
    </tr>
  </tbody>
</table>

The energy equation is solved as

    C(T) u·∇T = ∇·[D_eff(β, T) ∇T],    C(T) = ρ(T)c_p(T)/(ρc_p)_ref,
                                        D_f(T) = k_f(T)/(ρc_p)_ref

i.e. the ρc_p variation is carried on the **convective** term, scaled by a
*constant* reference (ρc_p)_ref. That is not cosmetic: dividing by the *local*
ρc_p instead, putting α(T) = k/(ρc_p) straight into the Laplacian, breaks
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
cases take an unchanged code path. `rhoTable`
appears in both dictionaries because OpenFOAM builds the transport model and
the primal solver from different files; **keep them consistent. They are not
cross-checked.** A variable-property case must also declare the adjoint's
C-weighted convective flux in `fvSchemes`, which must *not* be `bounded`:

```
div(-phiC,Ta)   Gauss limitedLinear 1;
```

**Adjoint treatment: frozen-property linearisation.** ∂C/∂T, ∂D/∂T and ∂ν/∂T
are not differentiated; C, D_eff and ν are frozen coefficient fields at the
primal solution. The cost of that is measured, not assumed.

**Verification** (`cases/varprops`, laminar production configuration, all five
tables active at once: ρ, c_p, k_f, μ *and* D_s): single-cell and directional
finite-difference checks retain the correct sign and remain close to the
constant-property campaign. The public scripts and logs in `cases/varprops`
contain the exact numerical rows. This README keeps only the generic conclusion
so it is not confused with any project-specific validation report. The
manufactured property span is deliberately strong enough to exercise the
variable-property paths, including μ(T) in the momentum equation through
`viscosityModels::temperatureTable`.

The tables are the mechanism; the verification values are chosen for FD
conditioning. Production runs load real coolant and structural data through the
same entries, for example representative high-pressure water and tabulated
solid conductivity.

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
fixedGradient elsewhere. The solver drives objective-patch gradients).
The Brinkman source is the upstream `topOSource` fvOption on `(U Ua)`;
initial designs are supplied through the `alpha` field.

Working examples under `cases/`: `fdcheck` (verification),
`demo2d` (heat-extraction maximisation under pressure-drop and volume
constraints). In the demo, the optimiser grows a conductive fin on the hot wall
while the nullSpace update brings total-pressure losses to the configured cap:

![fin growth demo](docs/figures/demo_fin_growth.png)

## Practical notes (learned the hard way)

- Solve `T` and `Ta` with tight inner tolerances (`relTol 0;`). They are
  linear given the flow; loose inner solves turn convergence into a crawl.
- Do **not** under-relax `Ta`: even `relax(1.0)` clips for diagonal
  dominance against the current field and biases the fixed point. Leave
  `Ta` out of `relaxationFactors`.
- Near sharp Brinkman interfaces, prefer `bounded Gauss upwind` /
  `limitedLinear` over `cellLimited`-gradient `linearUpwind`: limiter
  flip-flop can stall SIMPLE at a constant residual plateau.
- Small thermal objective weights can improve adjoint solver conditioning.
  In constrained optimisation, use `normalise` or `normFactor` so objective
  and constraint gradients are comparably scaled.
- The `nullSpace` update method has proven more robust than ISQP on
  badly scale-separated constraint sets in our tests.

## Licence

GPL-3.0-or-later (as required for OpenFOAM-linking code). Copyright (C)
2026 Fira Software Ltd.
