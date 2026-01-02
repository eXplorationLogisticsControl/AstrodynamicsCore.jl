"""Orbital elements diagram"""

using GLMakie
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))

GM = 398600.44
DU = 6378.137
VU = sqrt(GM / DU)

MU = 1.0

KEP0 = [1.0, 0.6, deg2rad(60), deg2rad(300), deg2rad(50), deg2rad(72)]
RV0 = AstrodynamicsCore.kep2rv(KEP0, MU)
RV_perigee = AstrodynamicsCore.kep2rv([KEP0[1:5]; 0.0], MU)
RV_rising = AstrodynamicsCore.kep2rv([KEP0[1:5]; deg2rad(-50)], MU)
RV_falling = AstrodynamicsCore.kep2rv([KEP0[1:5]; deg2rad(130)], MU)
orbit = AstrodynamicsCore.Planet(MU, 0.0, RV0, "final")
orbit_rvs = hcat([AstrodynamicsCore.eph(orbit, t) for t in LinRange(0.0, orbit.period, 200)]...)

h_dir = cross(RV0[1:3], RV0[4:6])
h_dir = h_dir / norm(h_dir)

# plot planet
fig = Figure(size=(800,800))
ax = Axis3(fig[1,1]; aspect = :data, azimuth=deg2rad(260), elevation=deg2rad(15))

# plot planet mesh
gb_mesh = AstrodynamicsCore.get_sphere_mesh(0.1, 30)
wireframe!(ax, gb_mesh, color=(:grey, 0.8), linewidth=0.3, transparency=true)

# # plot inertial axes
# arrow_length = 0.4
# arrows3d!(ax, [0.0, 0.0, 0.0], arrow_length*[1.0, 0.0, 0.0], color=:red)
# arrows3d!(ax, [0.0, 0.0, 0.0], arrow_length*[0.0, 1.0, 0.0], color=:green)
# arrows3d!(ax, [0.0, 0.0, 0.0], arrow_length*[0.0, 0.0, 1.0], color=:blue)

# plot orbit
lines!(ax, orbit_rvs[1,:], orbit_rvs[2,:], orbit_rvs[3,:], color=:black)

# plot orbital plane
# lower = [Point3f(-1,-1,-0.1), Point3f(-1,1,-0.1), Point3f(1,1,-0.1), Point3f(1,-1,-0.1)]
# upper = [Point3f(-1,-1,0), Point3f(-1,1,0), Point3f(1,1,0), Point3f(1,-1,0)]
# band!(ax, lower, upper, color=:blue, alpha=0.1)
surface!(ax, LinRange(-1.5,0.8,100), LinRange(-1.5,1.5,100), (x,y) -> -1e-2, alpha=0.6)

# plot satellite position
lines!(ax, [0.0, RV_rising[1]], [0.0, RV_rising[2]], [0.0, RV_rising[3]], color=:blue)
lines!(ax, [0.0, RV_falling[1]], [0.0, RV_falling[2]], [0.0, RV_falling[3]], color=:blue, linestyle=:dash)
lines!(ax, [0.0, RV_perigee[1]], [0.0, RV_perigee[2]], [0.0, RV_perigee[3]], color=:red)
lines!(ax, [0.0, RV0[1]], [0.0, RV0[2]], [0.0, RV0[3]], color=:black)
scatter!(ax, [RV_rising[1]], [RV_rising[2]], [RV_rising[3]], color=:blue, marker=:circle)
scatter!(ax, [RV0[1]], [RV0[2]], [RV0[3]], color=:black, marker=:utriangle, markersize=15)
hidedecorations!(ax)

save(joinpath(@__DIR__, "plots/ex_elements.png"), fig; px_per_unit=5)
display(fig)