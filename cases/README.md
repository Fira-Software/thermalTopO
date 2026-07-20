# Cases

## fdcheck (60x24, laminar 2D CHT verification)
The gradient-verification harness. Dictionaries:
- `optimisationDict.full` / `.primalonly`: zone-mean objective, no regularisation
- `optimisationDict.pnorm` / `.pnorm_primalonly`: patch p-norm objective
- `optimisationDict.pnorm_reg*`: production configuration, regularisation on
- `optimisationDict.pnorm_regT*`: as above plus tabulated D_s(T)

Drivers (run after `./Allrun` has set up mesh/zones/fields):
- `fd_driver.py`: directional derivatives, zone-mean objective
- `fd_cells.py`: single-cell central differences, three regimes
- `fd_pnorm.py`, `fd_pnorm_reg.py`, `fd_kT.py`: p-norm campaigns
  (unregularised / production / D_s(T))

Published results are quoted in the top-level README and derivation note.
With regularisation active, compare FD against `topOSens<solver>`
(`writeAllFields true`), not `topologySens<solver>`.

## fdcheck-fine (120x48)
4x refinement of fdcheck for the scheme-consistency study.

## varprops (60x24, laminar 2D CHT, temperature-dependent properties)
Gradient verification with rho(T), cp(T), k_f(T), mu(T) in the fluid AND
D_s(T) in the solid, all active in primal, adjoint and sensitivities.
Same geometry, design field and production configuration (regularisation +
patch p-norm) as `fdcheck`, so the numbers are directly comparable.

- `constant/transportProperties` selects the `temperatureTable` viscosity
  model (mu(T)/rho(T)); switch it back to `Newtonian` for a constant-nu
  comparison run.
- `system/fvSchemes` declares `div(-phiC,Ta)` for the C(T)-weighted adjoint
  flux. It must NOT be `bounded`; see derivation note section 6.2.
- Dictionaries: `optimisationDict.varprops` / `.varprops_primalonly`
- Driver: `./fd_varprops.py` (after `./Allrun`)

Manufactured verification fluid: at the reference temperature it reduces to
the constant-property fdcheck case, while the tabulated ranges exercise
volumetric heat capacity, viscosity, fluid conductivity and solid conductivity.
The exact numerical rows are kept in the reproducible driver outputs rather
than in this public summary.

## demo2d
Constrained heat-transfer enhancement demonstration: maximise downstream
zone temperature subject to a total-pressure-loss cap and a solid-volume cap;
nullSpace update. Reproduces the fin-growth figure via `./Allrun`.
