"""Tests for miscellaneous helpers."""

using GeometryBasics
using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end

@testset "get_sphere_mesh" begin
    mesh = AstrodynamicsCore.get_sphere_mesh(1.0, 12)
    @test mesh isa Mesh
    @test length(mesh.position) > 0
end

@testset "moving_average" begin
    vs = collect(1.0:10)
    avg = AstrodynamicsCore.moving_average(vs, 3)
    @test length(avg) == length(vs) - 2
    @test avg[1] ≈ 2.0
    @test AstrodynamicsCore.moving_average([4.0], 3) ≈ [4.0]
    @test_throws ArgumentError AstrodynamicsCore.moving_average(vs, 0)
end

@testset "_moving_average_padded" begin
    vs = collect(1.0:5)
    padded = AstrodynamicsCore._moving_average_padded(vs, 3)
    @test length(padded) == length(vs)
end
