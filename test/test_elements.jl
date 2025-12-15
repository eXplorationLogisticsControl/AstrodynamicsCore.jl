"""Test for elements functions"""

using Test

include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))

@testset "Elements" begin
    # from https://naif.jpl.nasa.gov/pub/naif/PSYCHE/kernels/spk/de440.bsp.lbl
    # test case 1
    MU = 398600.435507
    COE = [8000, 0.12, deg2rad(85), deg2rad(240), deg2rad(200), deg2rad(130)]
    R_CHECK = [-4022.024182950954, -6221.706935396668, -4255.655925778433]
    V_CHECK = [-1.519176635355084, -3.564858730446071, 5.33536331860247]
    RV_CHECK = [R_CHECK; V_CHECK]
    RV = AstrodynamicsCore.kep2rv(COE, MU)
    @test RV ≈ RV_CHECK rtol=1e-8

    # test case 2
    COE = [26590, 0.06, deg2rad(15), deg2rad(109), deg2rad(27), deg2rad(58)];
    R_CHECK = [-24090.96849278776, -5928.260586177573, 6620.625587995996];
    V_CHECK = [0.7942641884144943, -3.924902041686212, 0.1411640937441778];
    RV_CHECK = [R_CHECK; V_CHECK]
    RV = AstrodynamicsCore.kep2rv(COE, MU)
    @test RV ≈ RV_CHECK rtol=1e-8
end