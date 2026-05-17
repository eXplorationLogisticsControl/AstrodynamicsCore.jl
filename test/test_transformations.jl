"""Tests for rotation helpers."""

using LinearAlgebra
using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end

@testset "rotation matrices" begin
    φ = deg2rad(25)
    @test AstrodynamicsCore._rotmat_ax1(φ) ≈ [
        1 0 0
        0 cos(φ) sin(φ)
        0 -sin(φ) cos(φ)
    ]
    @test AstrodynamicsCore._rotmat_ax2(φ) ≈ [
        cos(φ) 0 -sin(φ)
        0 1 0
        sin(φ) 0 cos(φ)
    ]
    @test AstrodynamicsCore._rotmat_ax3(φ) ≈ [
        cos(φ) sin(φ) 0
        -sin(φ) cos(φ) 0
        0 0 1
    ]
end

@testset "perifocal2geocentric" begin
    vec_pf = [1.0, 2.0, 3.0]
    ω, i, Ω = deg2rad(10), deg2rad(20), deg2rad(30)
    v_gec = AstrodynamicsCore.perifocal2geocentric(vec_pf, ω, i, Ω)
    v_back =
        AstrodynamicsCore._rotmat_ax3(ω) *
        AstrodynamicsCore._rotmat_ax1(i) *
        AstrodynamicsCore._rotmat_ax3(Ω) *
        v_gec
    @test v_back ≈ vec_pf rtol = 1e-14
end
