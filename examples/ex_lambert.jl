"""Example solving Lambert's problem."""

using GLMakie
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))


r1 = [0.79, 0.0, 0.0]       # initial position vector
r2 = [-0.6, 0.27, 0.15]     # final position vector
tof = 12.21                 # time of fight
mu =  1.0                   # gravitational parameter
m = 1                       # max number of revolutions
cw = false                  # whether to take clockwise path (default: false)

res = AstrodynamicsCore.lambert(r1, r2, tof, m, mu, cw)

if m == 0
    x1 = [r1; res.v1]
else
    x1 = [r1; res.v1[end]]
    x2 = [r2; res.v2[end]]
end
times = LinRange(0.0, tof, 100)

rvs = AstrodynamicsCore.propagate_lagrangian(mu, x1, 0.0, times)

# plot
fig = Figure(size=(600,500))
ax = Axis3(fig[1,1]; aspect = :data, azimuth=deg2rad(260), elevation=deg2rad(15))
lines!(ax, rvs[1,:], rvs[2,:], rvs[3,:], color=:black)
lines!(ax, [0.0, r1[1]], [0.0, r1[2]], [0.0, r1[3]], color=:blue)
lines!(ax, [0.0, r2[1]], [0.0, r2[2]], [0.0, r2[3]], color=:green)

scatter!(ax, [r1[1]], [r1[2]], [r1[3]], color=:blue, markersize=10)
scatter!(ax, [r2[1]], [r2[2]], [r2[3]], color=:green, markersize=10)

#save(joinpath(@__DIR__, "plots/ex_lambert.png"), fig; px_per_unit=5)
display(fig)