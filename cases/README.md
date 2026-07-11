# Cases

## fdcheck (60x24, laminar 2D CHT verification)
The gradient-verification harness. Dictionaries:
- `optimisationDict.full` / `.primalonly` — zone-mean objective, no regularisation
- `optimisationDict.pnorm` / `.pnorm_primalonly` — patch p-norm objective
- `optimisationDict.pnorm_reg*` — production configuration (regularisation on)
- `optimisationDict.pnorm_regT*` — as above plus tabulated D_s(T)

Drivers (run after `./Allrun` has set up mesh/zones/fields):
- `fd_driver.py` — directional derivatives, zone-mean objective
- `fd_cells.py` — single-cell central differences, three regimes
- `fd_pnorm.py`, `fd_pnorm_reg.py`, `fd_kT.py` — p-norm campaigns
  (unregularised / production / D_s(T))

Published results are quoted in the top-level README and derivation note.
With regularisation active, compare FD against `topOSens<solver>`
(`writeAllFields true`), not `topologySens<solver>`.

## fdcheck-fine (120x48)
4x refinement of fdcheck for the scheme-consistency study.

## demo2d
Constrained heat-transfer enhancement demonstration: maximise downstream
zone temperature subject to total-pressure losses <= 2x baseline (active
constraint) and a solid-volume cap; nullSpace update. Reproduces the
fin-growth figure: `./Allrun` (about 2-3 hours serial, 80 cycles).
