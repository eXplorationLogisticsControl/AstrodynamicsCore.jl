"""Tests for mean / eccentric / true anomaly conversions."""

using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end

_wrap_angle(Δ) = acos(cos(Δ))

@testset "anomaly conversions" begin
    for (M, e) in [
        (deg2rad(30), 0.0),
        (deg2rad(120), 0.3),
        (deg2rad(200), 0.7),
        (1.25, 0.9),
    ]
        M_mod = mod(M, 2π)
        E = AstrodynamicsCore.ma2ea(M, e)
        @test E - e * sin(E) ≈ M_mod atol = 1e-10

        ta = AstrodynamicsCore.ma2ta(M, e)
        E_ta = AstrodynamicsCore.ta2ea(ta, e)
        @test isfinite(E_ta)
        @test _wrap_angle((E_ta - e * sin(E_ta)) - M_mod) ≈ 0.0 atol = 1e-10

        M_back = AstrodynamicsCore.ta2ma(ta, e)
        @test _wrap_angle(M_back - M_mod) ≈ 0.0 atol = 1e-10
    end
end
