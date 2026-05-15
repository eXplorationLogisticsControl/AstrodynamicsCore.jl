"""Test Kepler problem function"""

using LinearAlgebra
using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end

@testset "kepler" begin
    MU = 398600.435507
    R_CHECK = [-4022.024182950954, -6221.706935396668, -4255.655925778433]
    V_CHECK = [-1.519176635355084, -3.564858730446071, 5.33536331860247]
    RV_CHECK = [R_CHECK; V_CHECK]

    Δt = 1000.0
    RV_final = AstrodynamicsCore.propagate_lagrangian(MU, RV_CHECK, 0.0, Δt)

    RV0_back = AstrodynamicsCore.propagate_lagrangian(MU, RV_final, Δt, 0.0)
    Δr = RV0_back[1:3] - RV_CHECK[1:3]
    Δv = RV0_back[4:6] - RV_CHECK[4:6]
    @test norm(Δr) < 1e-11
    @test norm(Δv) < 1e-11
end

