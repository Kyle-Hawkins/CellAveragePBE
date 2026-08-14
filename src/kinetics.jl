# Monod-type specific growth kinetics. Arguments are *effective* concentrations
# (glucose `G`, ethanol `E`, dissolved oxygen `O`). Each rate is defined once, at
# top level, so the several call sites that need them cannot drift apart.
#
# Concentrations are clamped at zero (`pos`) before entering the saturation terms.
# The oxygen affinities are O(1e-3), so a transient overshoot to a slightly
# negative O would otherwise make `O/(K+O)` blow up — clamping keeps the rates
# physical (they vanish as a substrate is depleted) and the coupled ODE stable.

@inline pos(x) = max(x, zero(x))

"Specific fermentative glucose->biomass rate."
K_gf(G, p::PBEParameters) = p.mu_mgf * pos(G) / (p.K_mgf + pos(G))

"Specific oxidative glucose->biomass rate (oxygen-limited)."
function K_go(G, O, p::PBEParameters)
    g, o = pos(G), pos(O)
    return p.mu_mgo * g / (p.K_mgo + g) * o / (p.K_mgd + o)
end

"Specific oxidative ethanol->biomass rate (oxygen-limited, glucose-inhibited)."
function K_eo(G, E, O, p::PBEParameters)
    g, e, o = pos(G), pos(E), pos(O)
    return p.mu_meo * e / (p.K_meo + e) * o / (p.K_med + o) * p.K_inhib / (p.K_inhib + g)
end

"Net specific growth rate, summed over the three metabolic modes."
K_net(G, E, O, p::PBEParameters) = K_gf(G, p) + K_go(G, O, p) + K_eo(G, E, O, p)
