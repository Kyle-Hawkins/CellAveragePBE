# Cell Average Technique (Kumar, Peglow, Warnecke, Heinrich & Moerl, Chem. Eng.
# Sci. 61(10):3327-3342, 2006; extended to growth + breakage in Kumar, Peglow,
# Warnecke & Heinrich, Powder Technol., 2008) for the breakage + growth birth
# terms. The key numerical idea: all newborn cells produced anywhere in a bin are
# first averaged into a single representative mass (vbar = Vb/Bb below), and only
# then redistributed onto the two neighboring nodes, so that both number and mass
# are conserved on the coarse grid.
#
# This is what distinguishes it from the FIXED PIVOT technique of S. Kumar &
# Ramkrishna (Chem. Eng. Sci. 51(8):1311-1332, 1996), which splits each birth
# across its bracketing nodes individually. Averaging first retains information
# about where inside the cell the births actually concentrated; the fixed pivot
# discards it and over-predicts the number density in the large-mass tail on
# coarse grids. That matters here, where n is 60-80 cells.

"""
    breakage_births(N, S, grid, p) -> (Bb, Vb)

Number-birth `Bb[i]` and mass-birth `Vb[i]` due to division, for every cell `i`,
computed once per right-hand-side evaluation. Each is a sum over mother cells
`k ≥ i` of `N[k] · selection(x_k) · ∫ breakage` over the daughter range for cell
`i`. Cells with zero selection rate or zero population are skipped, which prunes
most of the nominal O(n²) work since only heavy cells divide.
"""
function breakage_births(N, S, grid::MassGrid, p::PBEParameters)
    mt = m_ts(S, p)
    md = m_ds(S, p)
    n = grid.n
    Bb = zeros(eltype(N), n)
    Vb = zeros(eltype(N), n)
    for k in 1:n
        Nk = N[k]
        Nk == 0 && continue
        yk = grid.centers[k]
        Sf = selection(yk, mt, md, p)
        Sf == 0 && continue
        w = Nk * Sf
        for i in 1:k
            lo = grid.bounds[i]
            hi = p_upper(grid, i, k)
            hi ≤ lo && continue
            bint, _ = quadgk(x -> breakage(x, yk, mt, p), lo, hi; rtol = 1e-6)
            vint, _ = quadgk(x -> x * breakage(x, yk, mt, p), lo, hi; rtol = 1e-6)
            Bb[i] += w * bint
            Vb[i] += w * vint
        end
    end
    return Bb, Vb
end

"""
    birth_CA(i, Bb, Vb, grid, p)

Cell-average division birth into interior cell `i` from the precomputed breakage
births (`Bb`, `Vb`). The average newborn mass in each of cells `i-1, i, i+1` is
`vbar = Vb / Bb`; those newborns are redistributed onto node `i` with linear
weights, gated by the location of each average mass relative to the cell nodes,
so that both number and mass are placed consistently on the coarse grid. (Growth
is handled separately as an advection flux; see `pbe_ode_term`.)
"""
function birth_CA(i::Int, Bb, Vb, grid::MassGrid, p::PBEParameters)
    c = grid.centers
    Bm, Vm = Bb[i - 1], Vb[i - 1]
    Bi, Vi = Bb[i],     Vb[i]
    Bp, Vp = Bb[i + 1], Vb[i + 1]

    vbar_m = iszero(Bm) ? zero(Vm) : Vm / Bm  # average newborn mass from cell i-1
    vbar_i = iszero(Bi) ? zero(Vi) : Vi / Bi
    vbar_p = iszero(Bp) ? zero(Vp) : Vp / Bp

    T1 = Bm * lambda_weight(grid, vbar_m, i, false) * ((vbar_m - c[i - 1]) > 0)
    T2 = Bi * lambda_weight(grid, vbar_i, i, false) * ((c[i] - vbar_i) > 0)
    T3 = Bi * lambda_weight(grid, vbar_i, i, true)  * ((vbar_i - c[i]) > 0)
    T4 = Bp * lambda_weight(grid, vbar_p, i, true)  * ((c[i + 1] - vbar_p) > 0)
    return T1 + T2 + T3 + T4
end
