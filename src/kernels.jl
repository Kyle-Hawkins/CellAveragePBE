# Population-balance kernels: division-mass thresholds, breakage (daughter
# distribution), selection (division rate), and ethanol formation. Each is written
# as plain branchy Julia rather than as a product of step functions: the branches
# state the gating conditions directly, and they are faster.

"""
    m_ts(S, p)

Transition (minimum-division) mass as a piecewise-linear function of substrate
concentration `S`: constant below `S_l`, linear between `S_l` and `S_h`, constant
above `S_h`. Continuous at both breakpoints.
"""
function m_ts(S, p::PBEParameters)
    if S < p.S_l
        return p.m_t0 + p.K_t * (p.S_l - p.S_h)
    elseif S ≤ p.S_h
        return p.m_t0 + p.K_t * (S - p.S_h)
    else
        return p.m_t0
    end
end

"Mean daughter (division) mass, piecewise-linear in substrate `S`."
function m_ds(S, p::PBEParameters)
    if S < p.S_l
        return p.m_d0 + p.K_d * (p.S_l - p.S_h)
    elseif S ≤ p.S_h
        return p.m_d0 + p.K_d * (S - p.S_h)
    else
        return p.m_d0
    end
end

"""
    breakage(x, y, mt, p)

Daughter-mass distribution: probability density that a dividing mother of mass
`y` produces a daughter of mass `x`. Two Gaussian lobes (at `mt` and at the
complementary mass `y - mt`) enforce mass-conserving asymmetric division. Returns
0 unless the daughter is lighter than the mother and the mother exceeds the
minimum divisible mass `mt + m_0`.
"""
function breakage(x, y, mt, p::PBEParameters)
    (y - x) < 0 && return zero(x)             # daughter cannot exceed mother
    (y - mt - p.m_0) ≤ 0 && return zero(x)    # mother below minimum divisible mass
    return p.A * exp(-p.beta * (x - mt)^2) + p.A * exp(-p.beta * (x - y + mt)^2)
end

"""
    selection(x, mt, md, p)

Division-rate kernel: rate at which a cell of mass `x` divides. A Gaussian
centered on the division mass `md` (active between `mt - m_0` and `md`) plus a
constant plateau `gamma` for cells heavier than `md`.
"""
function selection(x, mt, md, p::PBEParameters)
    term1 = ((x - mt + p.m_0) ≥ 0 && (md - x) ≥ 0) ?
            p.gamma * exp(-p.eps * (x - md)^2) : zero(x)
    term2 = (x - md) > 0 ? p.gamma : zero(x)
    return term1 + term2
end

"""
    ethanol_formation(x, mt, p)

Mass-specific ethanol-formation weight for a cell of mass `x`; a Gaussian offset
above the transition mass `mt`, active only for cells heavier than `mt`.
"""
function ethanol_formation(x, mt, p::PBEParameters)
    (x - mt) > 0 || return zero(x)
    return p.gamma_e * exp(-p.eps_e * (x - mt - p.me)^2)
end
