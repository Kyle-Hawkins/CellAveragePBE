"""
    PBEParameters

Physical and numerical constants for the cell-average population-balance model of
yeast fermentation with Henson-type kinetics. The kinetic, kernel and yield values
are the published parameters of the Henson-group yeast cell population balance model
(Zhu et al., Chem. Eng. Sci. 55:6155-6167, 2000; Zhang et al., J. Process Control
12:721-734, 2002) -- they were NOT fitted to data here. Collecting every constant in one
immutable, keyword-constructed struct keeps them out of the function bodies, which
makes the model reproducible and easy to perturb for sensitivity studies: any subset
is overridden by keyword, e.g. `PBEParameters(D = 0.3)`.

Units: masses in 1e-13 g, concentrations in g/L, rates in 1/h.
"""
Base.@kwdef struct PBEParameters
    # --- Monod growth kinetics ---
    mu_mgf::Float64 = 30.0    # max glucose->biomass rate (fermentative)
    mu_mgo::Float64 = 3.25    # max glucose->biomass rate (oxidative)
    mu_meo::Float64 = 7.0     # max ethanol->biomass rate (oxidative)
    K_mgf::Float64  = 40.0    # glucose affinity, fermentative
    K_mgo::Float64  = 2.0     # glucose affinity, oxidative
    K_meo::Float64  = 1.3     # ethanol affinity, oxidative
    K_mgd::Float64  = 0.001   # oxygen affinity, glucose-oxidative
    K_med::Float64  = 0.001   # oxygen affinity, ethanol-oxidative
    K_inhib::Float64 = 0.4    # glucose inhibition of ethanol uptake

    # --- Division mass thresholds (piecewise-linear in substrate S) ---
    S_l::Float64  = 0.1       # lower substrate breakpoint
    S_h::Float64  = 2.0       # upper substrate breakpoint
    K_t::Float64  = 0.01      # slope of transition mass m_t(S)
    m_t0::Float64 = 4.55      # baseline transition mass
    K_d::Float64  = 3.83      # slope of division mass m_d(S)
    m_d0::Float64 = 10.25     # baseline division mass
    m_0::Float64  = 1.0       # minimum mass increment required to divide

    # --- Breakage (daughter-mass) kernel: sum of two Gaussians ---
    A::Float64    = sqrt(10 / π)  # amplitude
    beta::Float64 = 4.0           # width

    # --- Selection (division-rate) kernel ---
    gamma::Float64 = 400.0    # amplitude
    eps::Float64   = 7.0      # width

    # --- Ethanol-formation kernel ---
    gamma_e::Float64 = 8.0
    eps_e::Float64   = 20.0
    me::Float64      = 1.54

    # --- Growth-birth regularization ---
    x0::Float64  = 1e-5
    nt0::Float64 = 1.0

    # --- Bioreactor chemistry ---
    G_f::Float64     = 30.0   # feed glucose
    D::Float64       = 0.15   # dilution rate
    Y_gf::Float64    = 0.15   # yield glucose->biomass (fermentative)
    Y_go::Float64    = 0.65   # yield glucose->biomass (oxidative)
    Y_eo::Float64    = 0.5    # yield ethanol->biomass (oxidative)
    E_f::Float64     = 0.0    # feed ethanol
    alpha_g::Float64 = 20.0   # glucose effective-concentration relaxation
    alpha_e::Float64 = 20.0   # ethanol effective-concentration relaxation
    Kloa::Float64    = 1500.0 # O2 mass-transfer coefficient
    H_O::Float64     = 0.032  # O2 Henry constant
    F::Float64       = 90.0   # gas flow rate
    O_in::Float64    = 0.21   # inlet O2 partial pressure
    Vl::Float64      = 0.1    # liquid volume
    Vg::Float64      = 0.9    # gas volume
    C_in::Float64    = 0.0003 # inlet CO2 partial pressure
    Klca::Float64    = 1500.0 # CO2 mass-transfer coefficient
    H_C::Float64     = 0.83   # CO2 Henry constant
end
