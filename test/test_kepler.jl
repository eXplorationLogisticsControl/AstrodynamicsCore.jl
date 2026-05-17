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

@testset "hypertrig" begin
    @test isfinite(AstrodynamicsCore.hypertrig_s(1.0))
    @test isfinite(AstrodynamicsCore.hypertrig_s(-1.0))
    @test isfinite(AstrodynamicsCore.hypertrig_s(1e-10))
    @test isfinite(AstrodynamicsCore.hypertrig_c(1.0))
    @test isfinite(AstrodynamicsCore.hypertrig_c(-1.0))
    @test isfinite(AstrodynamicsCore.hypertrig_c(1e-10))
end

@testset "propagate_lagrangian stm" begin
    MU = 398600.435507
    R_CHECK = [-4022.024182950954, -6221.706935396668, -4255.655925778433]
    V_CHECK = [-1.519176635355084, -3.564858730446071, 5.33536331860247]
    RV_CHECK = [R_CHECK; V_CHECK]
    Δt = 1000.0
    state1, Φ = AstrodynamicsCore.propagate_lagrangian(
        MU,
        RV_CHECK,
        0.0,
        Δt,
        1e-12,
        20,
        false,
        true,
    )
    @test length(state1) == 6
    @test size(Φ) == (6, 6)
    rv_prop = AstrodynamicsCore.propagate_lagrangian(MU, RV_CHECK, 0.0, Δt)
    @test state1 ≈ rv_prop rtol = 1e-10
end

@testset "propagate_lagrangian batch" begin
    MU = 398600.435507
    RV_CHECK = [
        -4022.024182950954
        -6221.706935396668
        -4255.655925778433
        -1.519176635355084
        -3.564858730446071
        5.33536331860247
    ]
    ts = LinRange(0.0, 500.0, 4)
    rvs = AstrodynamicsCore.propagate_lagrangian(MU, RV_CHECK, 0.0, ts)
    @test size(rvs) == (6, length(ts))
    @test rvs[:, 1] ≈ RV_CHECK rtol = 1e-10
end

@testset "propagate_lagrangian verbose" begin
    MU = 398600.435507
    RV_CHECK = [
        -4022.024182950954
        -6221.706935396668
        -4255.655925778433
        -1.519176635355084
        -3.564858730446071
        5.33536331860247
    ]
    redirect_stdout(devnull) do
        rv = AstrodynamicsCore.propagate_lagrangian(
            MU,
            RV_CHECK,
            0.0,
            1000.0,
            1e-12,
            20,
            true,
            false,
        )
        @test length(rv) == 6
    end
end

@testset "propagate_lagrangian hyperbolic" begin
    MU = 398600.435507
    kep_hyp = [-25000.0, 1.4, deg2rad(25), deg2rad(40), deg2rad(15), deg2rad(5)]
    rv0 = AstrodynamicsCore.kep2rv(kep_hyp, MU)
    rv1 = AstrodynamicsCore.propagate_lagrangian(MU, rv0, 0.0, 200.0)
    @test length(rv1) == 6
    @test norm(rv1[1:3]) > 0
end

