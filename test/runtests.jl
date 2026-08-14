using CellAveragePBE
using Test

const p = PBEParameters()

@testset "MassGrid" begin
    g = MassGrid(60; xmax = 12.0)
    @test length(g.bounds) == g.n + 1
    @test length(g.centers) == g.n
    @test g.bounds[1] == 0.0 && g.bounds[end] ≈ 12.0
    # centers strictly inside their cells, monotone increasing
    @test all(g.bounds[i] < g.centers[i] < g.bounds[i + 1] for i in 1:g.n)
    @test issorted(g.centers)
    @test_throws ArgumentError MassGrid(2)
end

@testset "Kernel continuity & gating" begin
    # m_ts / m_ds continuous at the two substrate breakpoints
    for m in (m_ts, m_ds)
        @test m(p.S_l - 1e-8, p) ≈ m(p.S_l + 1e-8, p) atol = 1e-6
        @test m(p.S_h - 1e-8, p) ≈ m(p.S_h + 1e-8, p) atol = 1e-6
    end
    mt = m_ts(1.0, p)
    # breakage vanishes when daughter exceeds mother or mother too small to divide
    @test breakage(9.0, 6.0, mt, p) == 0.0            # daughter > mother
    @test breakage(1.0, mt + 0.5 * p.m_0, mt, p) == 0.0  # mother below mt + m_0
    @test breakage(mt, 12.0, mt, p) > 0.0             # peak near transition mass
    # selection is nonnegative and plateaus above the division mass
    md = m_ds(1.0, p)
    @test selection(md + 1.0, mt, md, p) ≈ p.gamma
    @test selection(0.0, mt, md, p) == 0.0
    # ethanol formation gated above transition mass
    @test ethanol_formation(mt - 1.0, mt, p) == 0.0
    @test ethanol_formation(mt + p.me, mt, p) > 0.0
end

@testset "Monod kinetics" begin
    @test K_gf(0.0, p) == 0.0
    @test K_gf(1e6, p) ≈ p.mu_mgf rtol = 1e-3      # saturates at mu_mgf
    @test K_net(4.0, 2.0, 0.1, p) ≈
          K_gf(4.0, p) + K_go(4.0, 0.1, p) + K_eo(4.0, 2.0, 0.1, p)
    @test K_net(4.0, 2.0, 0.1, p) > 0.0
end

@testset "Trapezoid rule" begin
    x = collect(range(0, 1; length = 101))
    @test trapz(x, x)     ≈ 0.5   rtol = 1e-6       # ∫₀¹ x dx
    @test trapz(x, x .^ 2) ≈ 1/3  rtol = 1e-3       # ∫₀¹ x² dx
end

@testset "RHS well-formed" begin
    g = MassGrid(40)
    u0 = initial_state(g, p)
    du = similar(u0)
    pbe_rhs!(du, u0, (g, p), 0.0)
    @test all(isfinite, du)
    @test du[1] == 0.0 && du[g.n] == 0.0            # boundary cells held fixed
    @test length(du) == g.n + 8
end
