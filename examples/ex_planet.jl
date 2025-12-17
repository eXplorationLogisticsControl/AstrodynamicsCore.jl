"""Example for planets"""

using GLMakie
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))

GM = 398600.44
DU = 6378.137
VU = sqrt(GM / DU)

MU = 1.0

LAT = deg2rad(23)
LON = deg2rad(0.0)
R0 = DU * [cos(LON)*cos(LAT), sin(LON)*cos(LAT), sin(LAT)] / DU

# define launch
launch_bearing = deg2rad(97)
V0 = 8.7 * [0.0, sin(launch_bearing), cos(launch_bearing)] / VU
RV0 = [R0; V0]
kep0 = AstrodynamicsCore.rv2kep(RV0, MU)

# get transfer rv's until apogee 
tf = AstrodynamicsCore.time_until_θ(kep0, π, MU)
transfer = AstrodynamicsCore.Planet(MU, t0, RV0, "transfer")
times = LinRange(0.0, tf, 100)
transfer_rvs = hcat([AstrodynamicsCore.eph(transfer, t) for t in times]...)

# circularize at apogee
ΔV = AstrodynamicsCore.tangential_circularize(AstrodynamicsCore.rv2kep(transfer_rvs[:,end], MU), MU)
RV0_orbit = transfer_rvs[:,end]
RV0_orbit[4:6] = (norm(RV0_orbit[4:6]) + ΔV) * RV0_orbit[4:6]/norm(RV0_orbit[4:6])
orbit = AstrodynamicsCore.Planet(MU, t0, RV0_orbit, "final")
orbit_rvs = hcat([AstrodynamicsCore.eph(orbit, t) for t in LinRange(0.0, orbit.period, 200)]...)

# plot planet
fig = Figure(size=(400,400))
ax = Axis3(fig[1,1]; aspect = :data)

# plot planet mesh
gb_mesh = AstrodynamicsCore.get_sphere_mesh(1.0, 30)
wireframe!(ax, gb_mesh, color=(:grey, 0.2), linewidth=0.3, transparency=true)

# plot orbit
lines!(ax, [0.0, R0[1]], [0.0, R0[2]], [0.0, R0[3]], color=:blue)
lines!(ax, transfer_rvs[1,:], transfer_rvs[2,:], transfer_rvs[3,:], color=:red)
lines!(ax, orbit_rvs[1,:], orbit_rvs[2,:], orbit_rvs[3,:], color=:green)

display(fig)