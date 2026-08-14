"""
    CellAveragePBE

A cell-average population-balance solver for yeast fermentation, coupling a
size-structured cell population (mesoscopic integro-PDE, discretized by the Cell
Average Technique of Kumar et al., 2006) to a well-mixed bioreactor chemistry
model (differential-algebraic-flavored ODEs for glucose, ethanol, oxygen, and
CO2). The yeast kinetics and division kernels follow the Henson-group cell
population balance model (Zhu et al. 2000; Zhang et al. 2002). See the README for
full citations.

Quick start:

```julia
using CellAveragePBE
out = solve_pbe(n = 60)                 # build + integrate the coupled system
N, chem = split_state(out.sol.u[end], out.grid.n)   # final population + chemistry
```
"""
module CellAveragePBE

using QuadGK: quadgk
using OrdinaryDiffEq
using ADTypes: AutoFiniteDiff

include("parameters.jl")
include("grid.jl")
include("kernels.jl")
include("kinetics.jl")
include("cellaverage.jl")
include("chemistry.jl")
include("model.jl")

export PBEParameters, MassGrid, ChemState
export m_ts, m_ds, breakage, selection, ethanol_formation
export K_gf, K_go, K_eo, K_net
export breakage_births, birth_CA, chemistry_derivs, trapz
export pbe_rhs!, pbe_ode_term, initial_state, solve_pbe, split_state
export pbe_population_rhs!, solve_population

end # module
