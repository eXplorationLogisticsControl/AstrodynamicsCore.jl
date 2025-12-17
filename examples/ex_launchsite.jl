"""Example for launch site assessment"""

using GLMakie
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))

GM = 398600.44
DU = 6378.137
VU = sqrt(GM / DU)

MU = 1.0

LAT = deg2rad(28.6)
LON = deg2rad(0.0)
T_Inr2site = AstrodynamicsCore._rotmat_ax2(LAT) * AstrodynamicsCore._rotmat_ax3(LON)
R0 = T_Inr2site * [1.0, 0.0, 0.0]

function get_launch_trajectory(R0, RA::Real, launch_bearing::Real)
    # get transfer rv's until apogee 
    SMA = (norm(R0) + RA) / 2
    V0mag = sqrt(MU*(2/norm(R0) - 1/SMA))
    V0 = V0mag * T_Inr2site * [0.0, sin(launch_bearing), cos(launch_bearing)]
    RV0 = [R0; V0]
    kep0 = AstrodynamicsCore.rv2kep(RV0, MU)

    tf = AstrodynamicsCore.time_until_θ(kep0, π, MU)
    transfer = AstrodynamicsCore.Planet(MU, 0.0, RV0, "transfer")
    times = LinRange(0.0, tf, 100)
    transfer_rvs = hcat([AstrodynamicsCore.eph(transfer, t) for t in times]...)
    
    # circularize at apogee
    ΔV = AstrodynamicsCore.tangential_circularize(AstrodynamicsCore.rv2kep(transfer_rvs[:,end], MU), MU)
    RV0_orbit = transfer_rvs[:,end]
    RV0_orbit[4:6] = (norm(RV0_orbit[4:6]) + ΔV) * RV0_orbit[4:6]/norm(RV0_orbit[4:6])
    orbit = AstrodynamicsCore.Planet(MU, 0.0, RV0_orbit, "final")
    orbit_rvs = hcat([AstrodynamicsCore.eph(orbit, t) for t in LinRange(0.0, orbit.period, 200)]...)
    kep_final = AstrodynamicsCore.rv2kep(orbit_rvs[:,1], MU)
    return V0mag, ΔV, transfer_rvs, orbit_rvs, kep_final
end

# plot planet
fig = Figure(size=(800,400))
ax = Axis3(fig[1,1]; aspect = :data)

# plot planet mesh
gb_mesh = AstrodynamicsCore.get_sphere_mesh(1.0, 30)
wireframe!(ax, gb_mesh, color=(:grey, 0.2), linewidth=0.3, transparency=true)

# launch site
lines!(ax, [0.0, R0[1]], [0.0, R0[2]], [0.0, R0[3]], color=:blue)

# get transfer and orbit rv's
V0mags = Float64[]
keps_final = []
linewidth = 0.75
bearings = deg2rad.(LinRange(45,135,15))
for launch_bearing in bearings
    _V0mag, _ΔV, transfer_rvs, orbit_rvs, kep_final = get_launch_trajectory(R0, RA, launch_bearing)
    lines!(ax, transfer_rvs[1,:], transfer_rvs[2,:], transfer_rvs[3,:], color=:green, linestyle=:dash, linewidth=linewidth)
    lines!(ax, orbit_rvs[1,:], orbit_rvs[2,:], orbit_rvs[3,:], color=:green, linewidth=linewidth)
    push!(V0mags, _V0mag)
    push!(keps_final, kep_final)
end
keps_final = hcat(keps_final...)
ax_DV = Axis(fig[1,2]; xlabel="RAAN, deg", ylabel="Inclination, deg")
scatterlines!(ax_DV, rad2deg.(keps_final[4,:]), rad2deg.(keps_final[3,:]), color=:black, linewidth=linewidth)

display(fig)