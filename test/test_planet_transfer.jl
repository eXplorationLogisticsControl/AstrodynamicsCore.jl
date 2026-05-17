"""Tests for Planet, transfer, and time-to-anomaly helpers."""

using LinearAlgebra
using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end

const MU_pt = 398600.435507
const RV_pt = [
    -4022.024182950954
    -6221.706935396668
    -4255.655925778433
    -1.519176635355084
    -3.564858730446071
    5.33536331860247
]

@testset "Planet" begin
    planet = AstrodynamicsCore.Planet(MU_pt, 0.0, RV_pt, "test")
    @test planet.name == "test"
    @test planet.period > 0
    io = IOBuffer()
    show(io, planet)
    @test occursin("test", String(take!(io)))
    rv = AstrodynamicsCore.eph(planet, 500.0)
    @test length(rv) == 6
    @test norm(rv[1:3]) > 0
end

@testset "time_until_θ" begin
    kep = AstrodynamicsCore.rv2kep(RV_pt, MU_pt)
    @test AstrodynamicsCore.time_until_θ(kep, kep[6], MU_pt) ≈ 0.0 atol = 1e-12
    tf = AstrodynamicsCore.time_until_θ(kep, kep[6] + deg2rad(45), MU_pt)
    @test tf > 0
end

@testset "tangential_circularize" begin
    kep = AstrodynamicsCore.rv2kep(RV_pt, MU_pt)
    kep_apo = copy(kep)
    kep_apo[6] = π
    ΔV = AstrodynamicsCore.tangential_circularize(kep_apo, MU_pt)
    @test ΔV > 0
    kep_circ = copy(kep)
    kep_circ[2] = 0.0
    kep_circ[6] = 0.0
    @test AstrodynamicsCore.tangential_circularize(kep_circ, MU_pt) ≈ 0.0 atol = 1e-10
end
