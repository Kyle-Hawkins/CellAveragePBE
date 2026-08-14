# Well-mixed bioreactor chemistry: 8 ODEs for liquid- and gas-phase species
# coupled to the cell population through mass-specific uptake/production rates
# integrated over the population (∫ N(x) · rate(x) dx via the trapezoid rule).

"Composite trapezoidal integral of `y` over abscissae `x` (both length `n`)."
function trapz(x::AbstractVector, y::AbstractVector)
    s = zero(eltype(y))
    @inbounds for i in 1:length(x) - 1
        s += (x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2
    end
    return s
end

"State of the 8 bioreactor chemistry species (liquid G/E, effective G/E, O/CO2 liquid & gas)."
struct ChemState{T}
    G::T; G_eff::T; E::T; E_eff::T; O::T; O_out::T; C::T; C_out::T
end

"""
    chemistry_derivs(N, chem, grid, p) -> NTuple{8}

Time derivatives of the 8 chemistry species. Liquid glucose/ethanol/oxygen/CO2
exchange with feed, gas phase, and the population's metabolic demand; the
"effective" concentrations relax toward the bulk values (a first-order sensing
lag). Returns derivatives in the order
`(G, G_eff, E, E_eff, O, O_out, C, C_out)`.
"""
function chemistry_derivs(N, chem::ChemState, grid::MassGrid, p::PBEParameters)
    G, G_eff, E, E_eff = chem.G, chem.G_eff, chem.E, chem.E_eff
    O, O_out, C, C_out = chem.O, chem.O_out, chem.C, chem.C_out
    x = grid.centers
    S = G_eff + E_eff
    mt = m_ts(S, p)

    kgf = K_gf(G_eff, p)
    kgo = K_go(G_eff, O, p)
    keo = K_eo(G_eff, E_eff, O, p)
    fE = ethanol_formation.(x, mt, Ref(p))   # mass-specific ethanol weight per cell

    # Glucose: feed/dilution minus fermentative + oxidative uptake
    dG = p.D * (p.G_f - G) - trapz(x, N .* (kgf / p.Y_gf + kgo / p.Y_go))
    # Ethanol: dilution + fermentative production minus oxidative consumption
    dE = p.D * (p.E_f - E) + trapz(x, N .* (92 / 180 .* fE .* kgf / p.Y_gf .- keo / p.Y_eo))
    # First-order effective-concentration relaxation
    dG_eff = p.alpha_g * (G - G_eff)
    dE_eff = p.alpha_e * (E - E_eff)
    # Dissolved O2: gas-liquid transfer minus respiratory demand
    dO = p.Kloa * (p.H_O * O_out - O) - trapz(x, N .* (96 / 46 * keo / p.Y_eo + 192 / 180 * kgo / p.Y_go))
    dO_out = (p.F * (p.O_in - O_out) - p.Kloa * (p.H_O * O_out - O) * p.Vl) / p.Vg
    # Dissolved CO2: gas-liquid transfer plus respiratory + fermentative production
    dC = p.Klca * (p.H_C * C_out - C) +
         trapz(x, N .* (264 / 180 * kgo / p.Y_go + 88 / 46 * keo / p.Y_eo .+ 88 / 180 .* fE .* kgf / p.Y_gf))
    dC_out = (p.F * (p.C_in - C_out) - p.Klca * (p.H_C * C_out - C) * p.Vl) / p.Vg

    return (dG, dG_eff, dE, dE_eff, dO, dO_out, dC, dC_out)
end
