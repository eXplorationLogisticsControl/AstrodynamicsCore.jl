"""Test ordinary equinoctial element conversions."""

using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end


angle_difference_eq(x, y) = atan(sin(x - y), cos(x - y))


@testset "equinoctial elements conversions" begin
    cases = [
        ([7000.0, 0.0, 0.0, 0.0, 0.0, 0.3], 398600.435507),
        ([7000.0, 1e-12, 0.7, 0.4, 0.0, 1.2], 398600.435507),
        ([9000.0, 0.2, 1e-6, 0.2, 0.7, 2.1], 398600.435507),
        ([12000.0, 0.4, 0.8, 1.0, 0.5, 2.4], 398600.435507),
        ([30000.0, 0.85, 0.3, 2.0, 1.1, 0.2], 398600.435507),
        ([2000.0, 0.05, deg2rad(30.0), deg2rad(80.0), deg2rad(20.0), deg2rad(45.0)],
            4902.800066),
        ([10000.0, 0.1, deg2rad(170.0), deg2rad(40.0), deg2rad(15.0), deg2rad(5.0)],
            398600.435507),
    ]

    for (kep, mu) in cases
        rv = AstrodynamicsCore.kep2rv(kep, mu)
        mee = AstrodynamicsCore.kep2mee(kep)
        eq = AstrodynamicsCore.kep2eq(kep)

        kep_back = AstrodynamicsCore.eq2kep(eq)
        @test isapprox(kep_back, kep; atol = 1e-12)

        rv_from_eq = AstrodynamicsCore.eq2rv(eq, mu)
        eq_from_rv = AstrodynamicsCore.rv2eq(rv, mu)
        @test isapprox(rv_from_eq, rv; rtol = 2e-11, atol = 2e-11)
        @test isapprox(AstrodynamicsCore.eq2rv(eq_from_rv, mu), rv; rtol = 2e-11, atol = 2e-11)

        mee_from_eq = AstrodynamicsCore.eq2mee(eq)
        eq_from_mee = AstrodynamicsCore.mee2eq(mee)
        @test isapprox(AstrodynamicsCore.mee2rv(mee_from_eq, mu), rv; rtol = 2e-11, atol = 2e-11)
        @test isapprox(AstrodynamicsCore.eq2rv(eq_from_mee, mu), rv; rtol = 2e-11, atol = 2e-11)
        @test abs(angle_difference_eq(AstrodynamicsCore.mee2eq(mee_from_eq)[6], eq[6])) < 2e-11
    end
end


@testset "equinoctial elements domain" begin
    @test_throws DomainError AstrodynamicsCore.kep2eq(
        [-25000.0, 1.4, 0.2, 0.3, 0.4, 0.5])
    @test_throws DomainError AstrodynamicsCore.kep2eq(
        [7000.0, 0.1, pi, 0.3, 0.4, 0.5])
    @test_throws DomainError AstrodynamicsCore.eq2rv(
        [7000.0, 1.0, 0.0, 0.0, 0.0, 0.0], 398600.435507)
    @test_throws DomainError AstrodynamicsCore.eq2rv(
        [-7000.0, 0.1, 0.0, 0.0, 0.0, 0.0], 398600.435507)
end


@testset "equinoctial elements types" begin
    eq = BigFloat[7000, 0.1, -0.05, 0.2, -0.1, 0.7]
    mu = BigFloat("398600.435507")
    mee = AstrodynamicsCore.eq2mee(eq)
    rv = AstrodynamicsCore.eq2rv(eq, mu)
    @test eltype(mee) == BigFloat
    @test eltype(rv) == BigFloat
    @test isapprox(AstrodynamicsCore.eq2rv(AstrodynamicsCore.rv2eq(rv, mu), mu), rv; rtol = 1e-25)
end
