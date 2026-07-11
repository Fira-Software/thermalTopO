# Contributing to thermalTopO

Issues and pull requests are welcome: bug reports with a minimal case,
verification results on other configurations, new thermal objectives, and
extensions toward the roadmap items in the README (turbulent-diffusivity
adjoint chaining, refined boiling-onset constraints).

Ground rules:

- The code targets unmodified OpenFOAM v2512 (openfoam.com line) and follows
  the upstream adjointOptimisationFoam conventions; please match the
  surrounding style.
- Gradient-affecting changes must come with a finite-difference check
  (see cases/fdcheck and docs/derivation.md section 5 for the protocol and
  the convergence-guard requirements; note the topOSens vs topologySens
  caveat when regularisation is active).
- Licence: GPL-3.0-or-later. By submitting a contribution you agree to
  license it under the same terms. Keep copyright lines accurate.
- This project is not approved or endorsed by OpenCFD Limited.
