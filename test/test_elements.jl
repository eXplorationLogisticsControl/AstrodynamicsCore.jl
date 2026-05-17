"""Test for elements functions"""

using Test

if !@isdefined AstrodynamicsCore
    include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))
end

@testset "kep2rv" begin
    # from https://naif.jpl.nasa.gov/pub/naif/PSYCHE/kernels/spk/de440.bsp.lbl
    # test case 1
    MU = 398600.435507
    COE = [8000, 0.12, deg2rad(85), deg2rad(240), deg2rad(200), deg2rad(130)]
    R_CHECK = [-4022.024182950954, -6221.706935396668, -4255.655925778433]
    V_CHECK = [-1.519176635355084, -3.564858730446071, 5.33536331860247]
    RV_CHECK = [R_CHECK; V_CHECK]
    RV = AstrodynamicsCore.kep2rv(COE, MU)
    @test RV ≈ RV_CHECK rtol=1e-8

    COE_BACK = AstrodynamicsCore.rv2kep(RV, MU)
    ΔCOE = [COE_BACK[1] - COE[1]; 
            COE_BACK[2] - COE[2]; 
            acos(cos(COE_BACK[3] - COE[3]));
            acos(cos(COE_BACK[4] - COE[4]));
            acos(cos(COE_BACK[5] - COE[5]));
            acos(cos(COE_BACK[6] - COE[6]))];
    @test ΔCOE ≈ zeros(6) atol=1e-10

    # test case 2
    COE = [26590, 0.06, deg2rad(15), deg2rad(109), deg2rad(27), deg2rad(58)];
    R_CHECK = [-24090.96849278776, -5928.260586177573, 6620.625587995996];
    V_CHECK = [0.7942641884144943, -3.924902041686212, 0.1411640937441778];
    RV_CHECK = [R_CHECK; V_CHECK]
    RV = AstrodynamicsCore.kep2rv(COE, MU)
    @test RV ≈ RV_CHECK rtol=1e-8

    COE_BACK = AstrodynamicsCore.rv2kep(RV, MU)
    ΔCOE = [COE_BACK[1] - COE[1]; 
            COE_BACK[2] - COE[2]; 
            acos(cos(COE_BACK[3] - COE[3]));
            acos(cos(COE_BACK[4] - COE[4]));
            acos(cos(COE_BACK[5] - COE[5]));
            acos(cos(COE_BACK[6] - COE[6]))];
    @test ΔCOE ≈ zeros(6) atol=1e-10
end

@testset "mee2rv" begin
    # test case 1
    MU = 398600.435507
    COE = [8000, 0.12, deg2rad(85), deg2rad(240), deg2rad(200), deg2rad(130)]
    R_CHECK = [-4022.024182950954, -6221.706935396668, -4255.655925778433]
    V_CHECK = [-1.519176635355084, -3.564858730446071, 5.33536331860247]
    RV_CHECK = [R_CHECK; V_CHECK]

    MEE = AstrodynamicsCore.kep2mee(COE)
    RV = AstrodynamicsCore.mee2rv(MEE, MU)
    @test RV ≈ RV_CHECK rtol=1e-8

    MEE_BACK = AstrodynamicsCore.rv2mee(RV, MU)
    ΔMEE = [MEE_BACK[1] - MEE[1];
            MEE_BACK[2] - MEE[2];
            MEE_BACK[3] - MEE[3];
            MEE_BACK[4] - MEE[4];
            MEE_BACK[5] - MEE[5];
            acos(cos(MEE_BACK[6] - MEE[6]))];
    @test ΔMEE ≈ zeros(6) atol=1e-10

    # test case 2
    COE = [26590, 0.06, deg2rad(15), deg2rad(109), deg2rad(27), deg2rad(58)];
    R_CHECK = [-24090.96849278776, -5928.260586177573, 6620.625587995996];
    V_CHECK = [0.7942641884144943, -3.924902041686212, 0.1411640937441778];
    RV_CHECK = [R_CHECK; V_CHECK]

    MEE = AstrodynamicsCore.kep2mee(COE)
    RV = AstrodynamicsCore.mee2rv(MEE, MU)
    @test RV ≈ RV_CHECK rtol=1e-8

    MEE_BACK = AstrodynamicsCore.rv2mee(RV, MU)
    ΔMEE = [MEE_BACK[1] - MEE[1];
            MEE_BACK[2] - MEE[2];
            MEE_BACK[3] - MEE[3];
            MEE_BACK[4] - MEE[4];
            MEE_BACK[5] - MEE[5];
            acos(cos(MEE_BACK[6] - MEE[6]))];
    @test ΔMEE ≈ zeros(6) atol=1e-10
end

@testset "mee2kep" begin
    MU = 398600.435507
    COE = [8000, 0.12, deg2rad(85), deg2rad(240), deg2rad(200), deg2rad(130)]
    MEE = AstrodynamicsCore.kep2mee(COE)
    COE_back = AstrodynamicsCore.mee2kep(MEE)
    RV = AstrodynamicsCore.kep2rv(COE, MU)
    RV_back = AstrodynamicsCore.kep2rv(COE_back, MU)
    @test RV ≈ RV_back rtol = 1e-10
end

@testset "kep2rv_perifocal" begin
    MU = 398600.435507
    θ = deg2rad(130)
    kep = [8000.0, 0.12, 0.0, 0.0, 0.0, θ]
    rv_pf = AstrodynamicsCore.kep2rv_perifocal(kep, MU)
    rv = AstrodynamicsCore.kep2rv(kep, MU)
    @test rv ≈ rv_pf rtol = 1e-12
end

@testset "rv2kep equatorial" begin
    MU = 398600.435507
    r = 7000.0
    v = sqrt(MU / r)
    rv = [r, 0.0, 0.0, 0.0, v, 0.0]
    kep = AstrodynamicsCore.rv2kep(rv, MU)
    @test kep[3] ≈ 0.0 atol = 1e-12
    RV_back = AstrodynamicsCore.kep2rv(kep, MU)
    @test RV_back ≈ rv rtol = 1e-10
end

@testset "rv2mee retrograde" begin
    MU = 398600.435507
    COE = [8000, 0.12, deg2rad(85), deg2rad(240), deg2rad(200), deg2rad(130)]
    RV = AstrodynamicsCore.kep2rv(COE, MU)
    MEE_retro = AstrodynamicsCore.rv2mee(RV, MU, true)
    @test length(MEE_retro) == 6
    @test all(isfinite, MEE_retro)
end

@testset "hyperbolic kep2rv" begin
    MU = 398600.435507
    kep_hyp = [-25000.0, 1.4, deg2rad(25), deg2rad(40), deg2rad(15), deg2rad(5)]
    rv = AstrodynamicsCore.kep2rv(kep_hyp, MU)
    kep_back = AstrodynamicsCore.rv2kep(rv, MU)
    @test kep_back[2] ≈ kep_hyp[2] rtol = 1e-8
    rv_pf = AstrodynamicsCore.kep2rv_perifocal(kep_hyp, MU)
    @test length(rv_pf) == 6
end
