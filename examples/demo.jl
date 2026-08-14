# Minimal runnable demo of CellAveragePBE (no plotting dependencies).
# Run:  julia --project=. examples/demo.jl

using CellAveragePBE

println("Coupled bioreactor -- metabolic oscillations")
out = solve_pbe(n = 60, tspan = (0.0, 25.0), p = PBEParameters(D = 0.3), N_scale = 1/500)
sol = out.sol; n = out.grid.n
G = [split_state(u, n)[2].G for u in sol.u]
E = [split_state(u, n)[2].E for u in sol.u]
println("  retcode         = ", sol.retcode)
println("  glucose range   = ", round(minimum(G), sigdigits = 3), " .. ",
        round(maximum(G), sigdigits = 3), " g/L")
println("  ethanol maximum = ", round(maximum(E), sigdigits = 3), " g/L")

println("\nPopulation balance alone (fixed environment)")
outp = solve_population(n = 80, tspan = (0.0, 4.0))
println("  retcode         = ", outp.sol.retcode)
println("  final cell total = ", round(sum(outp.sol.u[end]), sigdigits = 4))
