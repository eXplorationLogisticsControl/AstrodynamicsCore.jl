"""Test Lambert problem function"""

using LinearAlgebra
using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end

function _concat_lambert_v1v2(r1, r2, tof, m, μ)
    res = AstrodynamicsCore.lambert(r1, r2, tof, m, μ)
    res.exitflag == 1 || error("Lambert failed")
    return vcat(res.v1, res.v2)
end

function fd_step(x)
    cbrt(eps(Float64)) * max(1.0, abs(x))
end

function fd_dv_dt(r1, r2, tof, m, μ)
    h = fd_step(tof)
    dv_dt1 = (_concat_lambert_v1v2(r1, r2, tof - h, m, μ) - _concat_lambert_v1v2(r1, r2, tof + h, m, μ)) / (2h)
    dv_dt2 = (_concat_lambert_v1v2(r1, r2, tof + h, m, μ) - _concat_lambert_v1v2(r1, r2, tof - h, m, μ)) / (2h)
    return hcat(dv_dt1, dv_dt2)
end

function fd_dv_dr(r1, r2, tof, m, μ)
    r = vcat(r1, r2)
    dv_dr = zeros(6, 6)
    for j in 1:6
        h = fd_step(r[j])
        r_plus = copy(r)
        r_plus[j] += h
        r_minus = copy(r)
        r_minus[j] -= h
        dv_dr[:, j] = (_concat_lambert_v1v2(r_plus[1:3], r_plus[4:6], tof, m, μ) -
                       _concat_lambert_v1v2(r_minus[1:3], r_minus[4:6], tof, m, μ)) / (2h)
    end
    return dv_dr
end

@testset "lambert" begin
    # initial and final condition
    r1vec = [0.79, 0.0, 0.0]
    r2vec = [-0.6, -0.17, 0.015]
    mu =  1.0
    m = 0
    tofs = LinRange(1.2, 10.0, 20)

    for tof in tofs
        res = AstrodynamicsCore.lambert(r1vec, r2vec, tof, m, mu)

        x0 = [r1vec; res.v1]
        RV_final = AstrodynamicsCore.propagate_lagrangian(mu, x0, 0.0, tof)
        @test res.exitflag == 1
        @test norm(RV_final - [r2vec; res.v2]) < 1e-12
    end
end


@testset "lambert_cw" begin
    r1vec = [0.79, 0.0, 0.0]
    r2vec = [-0.6, -0.17, 0.015]
    mu = 1.0
    m = 0
    tofs = LinRange(1.2, 10.0, 10)

    for tof in tofs
        res_ccw = AstrodynamicsCore.lambert(r1vec, r2vec, tof, m, mu, false)
        res_cw = AstrodynamicsCore.lambert(r1vec, r2vec, tof, m, mu, true)

        @test res_ccw.exitflag == 1
        @test res_cw.exitflag == 1
        @test res_ccw.v1 != res_cw.v1

        for res in (res_ccw, res_cw)
            x0 = [r1vec; res.v1]
            RV_final = AstrodynamicsCore.propagate_lagrangian(mu, x0, 0.0, tof)
            @test norm(RV_final - [r2vec; res.v2]) < 1e-12
        end
    end
end


@testset "lambert_multi_rev" begin
    r1vec = [1.0, 0.0, 0.0]
    r2vec = [0.0, 1.0, 0.0]
    mu = 1.0
    tof = 8.0

    res = AstrodynamicsCore.lambert(r1vec, r2vec, tof, 1, mu)
    @test res isa AstrodynamicsCore.LambertMultiResults
    @test res.exitflag == 1
    @test length(res.v1) == 3
    @test length(res.v2) == 3
    @test res.revs == [0, 1, 1]
    @test res.branch == [:zero, :left, :right]

    for k in eachindex(res.v1)
        x0 = [r1vec; res.v1[k]]
        RV_final = AstrodynamicsCore.propagate_lagrangian(mu, x0, 0.0, tof)
        @test norm(RV_final[1:3] - r2vec) < 1e-10
        @test norm(RV_final[4:6] - res.v2[k]) < 1e-10
    end
end


@testset "lambert_jacobians" begin
    # initial and final condition
    r1vec = [0.79, 0.0, 0.0]
    r2vec = [-0.6, -0.17, 0.015]
    mu =  1.0
    m = 0
    tofs = LinRange(1.2, 10.0, 20)

    for tof in tofs
        res, dv_dt, dv_dr = AstrodynamicsCore.lambert_jac(r1vec, r2vec, tof, m, mu)
        @test res.exitflag == 1
        @test dv_dt[:, 2] ≈ -dv_dt[:, 1] atol=1e-12

        dv_dt_fd = fd_dv_dt(r1vec, r2vec, tof, m, mu)
        dv_dr_fd = fd_dv_dr(r1vec, r2vec, tof, m, mu)

        @test dv_dt ≈ dv_dt_fd atol=1e-8 rtol=1e-8
        @test dv_dr ≈ dv_dr_fd atol=1e-8 rtol=1e-8
    end
end
