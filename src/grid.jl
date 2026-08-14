"""
    MassGrid(n; xmax = 12.0)

Uniform mass grid of `n` cells on `[0, xmax]`. Stores the `n+1` cell edges
(`bounds`) and the `n` cell midpoints (`centers`). Cell `i` occupies
`[bounds[i], bounds[i+1]]` with representative mass `centers[i]`.

The grid is uniform in mass. Storing edges and centers separately keeps the
finite-volume bookkeeping explicit: advective fluxes are evaluated on `bounds`,
while populations and kernels are evaluated at `centers`.
"""
struct MassGrid
    n::Int
    bounds::Vector{Float64}   # length n+1
    centers::Vector{Float64}  # length n
end

function MassGrid(n::Int; xmax::Float64 = 12.0)
    n ≥ 3 || throw(ArgumentError("need at least 3 cells (got $n)"))
    bounds = collect(range(0.0, xmax; length = n + 1))
    centers = [(bounds[i] + bounds[i + 1]) / 2 for i in 1:n]
    return MassGrid(n, bounds, centers)
end

"Upper integration limit for daughters of mother-cell `k` that land in cell `i` (`i ≤ k`)."
p_upper(grid::MassGrid, i::Int, k::Int) = k == i ? grid.centers[i] : grid.bounds[i + 1]

"""
    lambda_weight(grid, vbar, i, forward)

Linear cell-average interpolation weight used by the Cell Average Technique.
`forward = true` interpolates toward the upper neighbor `i+1`, `false` toward
the lower neighbor `i-1`. Valid for interior cells `2 ≤ i ≤ n-1`.
"""
function lambda_weight(grid::MassGrid, vbar::Float64, i::Int, forward::Bool)
    xref = forward ? grid.centers[i + 1] : grid.centers[i - 1]
    return (vbar - xref) / (grid.centers[i] - xref)
end
