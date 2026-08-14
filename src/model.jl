# Assembly of the full coupled system: n population-balance ODEs (one per mass
# cell) coupled to the 8 chemistry ODEs, plus initial conditions and a solve
# driver. State layout: u[1:n] = number density N, u[n+1:n+8] = chemistry.

"""
    pbe_ode_term(i, N, chem, Bb, Vb, grid, p)

Net rate of change of the number density in interior cell `i`:

  * **division birth**  — cell-average redistribution of breakage products,
  * **growth advection** — first-order upwind flux `-∂(g·N)/∂x` that moves cells
    up the mass axis at the mass-specific growth rate `g = K_net` (uniform in
    size, so the flux is `-(g·N_i - g·N_{i-1})/Δx`),
  * **death** — loss to dilution `D` and to division `selection(x_i)`.

Growth is now a genuine advection term rather than the self-canceling birth/death
of the original formulation, so cells grow toward the division mass and sustain
the population.
"""
function pbe_ode_term(i::Int, N, chem::ChemState, Bb, Vb, grid::MassGrid, p::PBEParameters)
    S = chem.G_eff + chem.E_eff
    g = K_net(chem.G_eff, chem.E_eff, chem.O, p)
    dx = grid.bounds[i + 1] - grid.bounds[i]
    mt = m_ts(S, p); md = m_ds(S, p)

    birth = birth_CA(i, Bb, Vb, grid, p)
    advection = -g * (N[i] - N[i - 1]) / dx            # upwind (g ≥ 0)
    death = (p.D + selection(grid.centers[i], mt, md, p)) * N[i]
    return birth + advection - death
end

"""
    pbe_rhs!(du, u, params, t)

In-place right-hand side of the coupled PBE + chemistry system. `params` is the
tuple `(grid, p)`. Boundary cells `1` and `n` are held fixed (no-flux end caps),
matching the original formulation.
"""
function pbe_rhs!(du, u, params, t)
    grid, p = params
    n = grid.n
    N = @view u[1:n]
    chem = ChemState(u[n + 1], u[n + 2], u[n + 3], u[n + 4],
                     u[n + 5], u[n + 6], u[n + 7], u[n + 8])

    S = chem.G_eff + chem.E_eff
    Bb, Vb = breakage_births(N, S, grid, p)   # once per evaluation, shared by all cells

    du[1] = 0.0
    du[n] = 0.0
    @inbounds for i in 2:n - 1
        du[i] = pbe_ode_term(i, N, chem, Bb, Vb, grid, p)
    end

    dchem = chemistry_derivs(N, chem, grid, p)
    @inbounds for j in 1:8
        du[n + j] = dchem[j]
    end
    return nothing
end

"""
    initial_state(grid, p; kwargs...)

Bimodal initial cell distribution (two Gaussian sub-populations) stacked with the
initial chemistry vector `(G, G_eff, E, E_eff, O, O_out, C, C_out)`.
"""
function initial_state(grid::MassGrid, p::PBEParameters;
                       center1 = 4.0, center2 = 7.5, sigma = 1.0, N1 = 350.0, N2 = 450.0,
                       chem0 = (4.0, 3.0, 2.0, 1.0, 0.1, 0.2, 0.05, 0.1))
    x = grid.centers
    N0 = @. N1 * exp(-(center1 - x)^2 / sigma^2) + N2 * exp(-(center2 - x)^2 / sigma^2)
    return vcat(N0, collect(float.(chem0)))
end

"""
    solve_pbe(; n, tspan, saveat, p, solver, kwargs...)

Build and solve the coupled model. Returns a named tuple `(sol, grid, p)` where
`sol` is the `OrdinaryDiffEq` solution (state length `n + 8`). A stiff solver is
used by default because the gas-liquid mass-transfer terms are stiff.
"""
function solve_pbe(; n::Int = 60, tspan = (0.0, 20.0),
                   saveat = range(tspan[1], tspan[2]; length = 200),
                   p::PBEParameters = PBEParameters(), N_scale = 1.0,
                   solver = FBDF(autodiff = AutoFiniteDiff()), reltol = 1e-6, abstol = 1e-8, kwargs...)
    grid = MassGrid(n)
    u0 = initial_state(grid, p)
    u0[1:n] .*= N_scale   # scale population to the chemistry (cell-number-to-concentration factor)
    prob = ODEProblem(pbe_rhs!, u0, tspan, (grid, p))
    sol = solve(prob, solver; saveat = saveat, reltol = reltol, abstol = abstol, kwargs...)
    return (sol = sol, grid = grid, p = p)
end

"""
    pbe_population_rhs!(du, N, params, t)

Right-hand side of the population balance alone, evolving the size distribution
`N` in a *fixed* chemical environment `env::ChemState`. `params = (grid, p, env)`.
Isolating the population from the (stiff, tightly coupled) chemistry gives a clean,
well-posed demonstration of the Cell Average Technique — division redistributes
newborns while growth advects mass up the grid.
"""
function pbe_population_rhs!(du, N, params, t)
    grid, p, env = params
    n = grid.n
    S = env.G_eff + env.E_eff
    Bb, Vb = breakage_births(N, S, grid, p)
    du[1] = 0.0
    du[n] = 0.0
    @inbounds for i in 2:n - 1
        du[i] = pbe_ode_term(i, N, env, Bb, Vb, grid, p)
    end
    return nothing
end

"""
    solve_population(; n, tspan, saveat, env, p, solver, kwargs...)

Integrate the population balance in a fixed chemical environment `env`. Returns
`(sol, grid, p, env)`; each `sol.u` is the length-`n` number density `N(t)`.
"""
function solve_population(; n::Int = 60, tspan = (0.0, 4.0),
                          saveat = range(tspan[1], tspan[2]; length = 100),
                          env::ChemState = ChemState(5.0, 5.0, 0.5, 0.5, 0.05, 0.0, 0.0, 0.0),
                          p::PBEParameters = PBEParameters(),
                          N0::Union{Nothing,AbstractVector} = nothing,
                          solver = FBDF(autodiff = AutoFiniteDiff()),
                          reltol = 1e-6, abstol = 1e-8, kwargs...)
    grid = MassGrid(n)
    u0 = N0 === nothing ? initial_state(grid, p)[1:n] : collect(float.(N0))
    length(u0) == n || throw(ArgumentError("N0 must have length n = $n"))
    prob = ODEProblem(pbe_population_rhs!, u0, tspan, (grid, p, env))
    sol = solve(prob, solver; saveat = saveat, reltol = reltol, abstol = abstol, kwargs...)
    return (sol = sol, grid = grid, p = p, env = env)
end

"Split a solution state vector into `(N, ChemState)` for cell count `n`."
function split_state(u, n::Int)
    N = u[1:n]
    chem = ChemState(u[n + 1], u[n + 2], u[n + 3], u[n + 4],
                     u[n + 5], u[n + 6], u[n + 7], u[n + 8])
    return N, chem
end
