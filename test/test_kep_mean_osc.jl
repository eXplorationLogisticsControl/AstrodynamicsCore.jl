"""Tests for Keplerian mean / osculating conversions under J2."""

using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end

const Re_test = 6378.1363
const J2_test = 1082.63e-6

@testset "kep_osc2mean" begin
    osc_coe = [
        42164.0
        0.0003
        deg2rad(0.14)
        deg2rad(30)
        0.0
        0.0
    ]
    mean_coe = AstrodynamicsCore.kep_osc2mean(osc_coe, Re_test, J2_test)
    mean_coe_check = [
        4.216399905013438e4
        2.628234529136558e-4
        0.002443415535344
        0.523598775598787
        3.644373691713554e-11
        -3.693134686955091e-11
    ]
    @test mean_coe ≈ mean_coe_check rtol = 1e-12
end

@testset "kep_osc2mean" begin
    osc_coe = [
        42165.727783078
        0.0003
        0.1389338295990
        deg2rad(30)
        0.0
        0.0
    ]
    mean_coe = AstrodynamicsCore.kep_osc2mean(osc_coe, Re_test, J2_test)
    osc_coe_back = AstrodynamicsCore.kep_mean2osc(mean_coe, Re_test, J2_test)
    @test osc_coe_back[1] ≈ osc_coe[1] rtol = 1e-7
    @test osc_coe_back[2:end] ≈ osc_coe[2:end] atol = 1e-7
end
