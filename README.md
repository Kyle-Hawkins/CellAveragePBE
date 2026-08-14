# CellAveragePBE

A cell-average **population-balance solver** for yeast fermentation, written in Julia. The package
couples a size-structured cell population — a mesoscopic **integro-PDE** for the number density
`N(x, t)` over cell mass `x` — to a well-mixed **bioreactor chemistry** model (glucose, ethanol,
oxygen, CO₂). The fully coupled system reproduces the sustained **metabolic oscillations**
characteristic of the Henson-group yeast model (Zhu et al. 2000; Zhang et al. 2002).

An independent Julia implementation of the published model described under
[References](#references) — written from the literature, not ported from an existing
codebase.

## What it demonstrates

- **Cell Average Technique** (Kumar et al. 2006) for the division/breakage integro-PDE — non-integer
  daughter masses are redistributed onto neighboring grid nodes with consistent number and mass placement.
- **Upwind finite-volume advection** for cell growth up the mass axis.
- **Coupled stiff ODE integration** — `n` population cells + 8 chemistry ODEs (gas–liquid mass transfer,
  substrate uptake integrated over the population), solved with a BDF method.
- **Idiomatic, documented, tested Julia** — parameters in one immutable struct, small composable functions,
  a 26-test suite.

## The model

```
∂N(x,t)/∂t + ∂[g·N(x,t)]/∂x  =  B_div(x,t)  −  [ S(x) + D ]·N(x,t)
```

- `g = K_net(...)` — mass-specific growth rate (Monod kinetics on the local chemistry)
- `B_div` — cell-average division birth from the breakage kernel `b(x, y)` and selection rate `S(x)`
- `D` — dilution rate

coupled to 8 chemistry ODEs for glucose, ethanol, their effective (sensed) concentrations, and
dissolved/gas-phase O₂ and CO₂.

## Layout

```
src/
  parameters.jl   all physical/numerical constants in one keyword struct
  grid.jl         mass grid + cell-average interpolation weights
  kernels.jl      division-mass thresholds, breakage, selection, ethanol-formation
  kinetics.jl     Monod growth kinetics (concentration-clamped for stability)
  cellaverage.jl  Cell Average Technique birth term
  chemistry.jl    the 8 coupled bioreactor ODEs
  model.jl        RHS assembly, initial conditions, solve drivers
test/runtests.jl  26 tests (grid, kernels, kinetics, quadrature, RHS)
notebook/         population_balance.ipynb + rendered PDF
```

## Running it

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()

using CellAveragePBE

# Coupled bioreactor — sustained metabolic oscillations
out = solve_pbe(n = 60, tspan = (0.0, 25.0), p = PBEParameters(D = 0.3), N_scale = 1/500)
N, chem = split_state(out.sol.u[end], out.grid.n)

# Population balance alone, in a fixed chemical environment
outp = solve_population(n = 80, tspan = (0.0, 4.0))
```

Run the tests with `julia --project=. -e 'using Pkg; Pkg.test()'` (or `include("test/runtests.jl")`).

The `notebook/` directory contains a narrated walk-through (`population_balance.ipynb`) and its rendered
PDF, covering the kernels, the growth-plus-division mechanism, and the coupled oscillations.

## References

The numerics and the yeast model are both taken from the literature; this package is an
independent Julia implementation of them. BibTeX entries for everything below are in
[`docs/references.bib`](docs/references.bib).

**Population balance discretization**

- J. Kumar, M. Peglow, G. Warnecke, S. Heinrich and L. Mörl, "Improved accuracy and
  convergence of discretized population balance for aggregation: the cell average technique,"
  *Chemical Engineering Science* **61**(10), 3327–3342, 2006. — **the cell average technique
  itself**, i.e. the scheme in `src/cellaverage.jl`.
- J. Kumar, M. Peglow, G. Warnecke and S. Heinrich, "An efficient numerical technique for
  solving population balance equation involving aggregation, breakage, growth and nucleation,"
  *Powder Technology*, 2008. — cell average extended to combined growth and breakage, which is
  the combination solved here.
- S. Kumar and D. Ramkrishna, "On the solution of population balance equations by
  discretization — I. A fixed pivot technique," *Chemical Engineering Science* **51**(8),
  1311–1332, 1996. — the *fixed pivot* technique: the predecessor the cell average method
  improves on, **not** the method implemented here. Note this is a different Kumar (S., with
  Ramkrishna) from the cell-average papers (J. Kumar) — the two are easy to conflate.
- D. Ramkrishna, *Population Balances: Theory and Applications to Particulate Systems in
  Engineering*, Academic Press, 2000. — general reference for the population balance itself.

**Yeast cell population balance model**

- G.-Y. Zhu, A. Zamamiri, M. A. Henson and M. A. Hjortsø, "Model predictive control of
  continuous yeast bioreactors using cell population balance models," *Chemical Engineering
  Science* **55**, 6155–6167, 2000.
- Y. Zhang, A. M. Zamamiri, M. A. Henson and M. A. Hjortsø, "Cell population models for
  bifurcation analysis and nonlinear control of continuous yeast bioreactors," *Journal of
  Process Control* **12**, 721–734, 2002.
- M. A. Henson, "Dynamic modeling of microbial cell populations," *Current Opinion in
  Biotechnology*, 2003. — review covering this model class.

The division kernels, the piecewise-linear transition/division masses, the three-mode Monod
kinetics and the effective (sensed) concentrations all follow the Henson-group model; the
kinetic and yield parameters in `src/parameters.jl` are that model's published values rather
than fitted ones.
