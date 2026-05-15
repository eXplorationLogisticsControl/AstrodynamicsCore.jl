"""Make plots of Hohmann transfers"""

using GLMakie

include(joinpath(@__DIR__, "..", "src", "AstrodynamicsCore.jl"))

# --------------------------------------------------------------------------------------- #
# circular to circular
MU = 1.0
R0 = 1.0
rv0 = AstrodynamicsCore.kep2rv([R0, 0.0, 0.0, 0.0, 0.0, 0.0], MU)
orbit0 = AstrodynamicsCore.Planet(MU, 0.0, rv0, "initial")
rvs_initial = hcat([AstrodynamicsCore.eph(orbit0, t) for t in LinRange(0.0, orbit0.period, 100)]...)

RF = 4.0
rvf = AstrodynamicsCore.kep2rv([RF, 0.0, 0.0, 0.0, 0.0, 0.0], MU)
orbitf = AstrodynamicsCore.Planet(MU, 0.0, rvf, "final")
rvs_final = hcat([AstrodynamicsCore.eph(orbitf, t) for t in LinRange(0.0, orbitf.period, 100)]...)

SMA_hohmann = (R0 + RF) / 2.0
ECC_hohmann = (RF - R0) / (RF + R0)
rv0_hohmann = AstrodynamicsCore.kep2rv([SMA_hohmann, ECC_hohmann, 0.0, 0.0, 0.0, 0.0], MU)
orbit0_hohmann = AstrodynamicsCore.Planet(MU, 0.0, rv0_hohmann, "hohmann")
rvs_hohmann = hcat([AstrodynamicsCore.eph(orbit0_hohmann, t) for t in LinRange(0.0, orbit0_hohmann.period/2, 100)]...)
rvs_hohmann_2 = hcat([AstrodynamicsCore.eph(orbit0_hohmann, t) for t in LinRange(orbit0_hohmann.period/2, orbit0_hohmann.period, 100)]...)

fontsize = 18
fig = Figure(size=(400,400))
ax = Axis(fig[1,1]; xlabel="x", ylabel="y", aspect=DataAspect(), xlabelsize=fontsize, ylabelsize=fontsize, xticklabelsize=fontsize, yticklabelsize=fontsize)
poly!(Circle(Point2f(0, 0), 0.5), color = :dodgerblue, label=nothing)
hidedecorations!(ax)
lines!(ax, rvs_initial[1,:], rvs_initial[2,:], color=:blue, label="Initial")
lines!(ax, rvs_final[1,:], rvs_final[2,:], color=:green, label="Final", linewidth=2.0)
lines!(ax, rvs_hohmann[1,:], rvs_hohmann[2,:], color=:red, label="HT", linewidth=2.0)
lines!(ax, rvs_hohmann_2[1,:], rvs_hohmann_2[2,:], color=:red, label=nothing, linewidth=2.0, linestyle=:dash)
scatter!(ax, rvs_hohmann[1,1], rvs_hohmann[2,1], color=:red, markersize=10)
scatter!(ax, rvs_hohmann[1,end], rvs_hohmann[2,end], color=:red, markersize=10)
Legend(fig[2,1], ax, orientation=:horizontal, labelsize=fontsize)
save(joinpath(@__DIR__, "plots/hohmann_circ2circ.png"), fig; px_per_unit=5)


# --------------------------------------------------------------------------------------- #
# elliptical to elliptical (opposite apsis)
MU = 1.0
R0 = 1.9
ECC0 = 0.45
rv0 = AstrodynamicsCore.kep2rv([R0, ECC0, 0.0, 0.0, 0.0, 0.0], MU)
orbit0 = AstrodynamicsCore.Planet(MU, 0.0, rv0, "initial")
rvs_initial = hcat([AstrodynamicsCore.eph(orbit0, t) for t in LinRange(0.0, orbit0.period, 100)]...)

RF = 7.0
ECCf = 0.33
rvf = AstrodynamicsCore.kep2rv([RF, ECCf, 0.0, 0.0, π, 0.0], MU)
orbitf = AstrodynamicsCore.Planet(MU, 0.0, rvf, "final")
rvs_final = hcat([AstrodynamicsCore.eph(orbitf, t) for t in LinRange(0.0, orbitf.period, 100)]...)

# hohmann option 1
RP0 = R0 * (1 - ECC0)
RA0 = RF * (1 - ECCf)
SMA_hohmann = (RP0 + RA0) / 2.0
ECC_hohmann = (RA0 - RP0) / (RA0 + RP0)
rv0_hohmann = AstrodynamicsCore.kep2rv([SMA_hohmann, ECC_hohmann, 0.0, 0.0, 0.0, 0.0], MU)
orbit0_hohmann = AstrodynamicsCore.Planet(MU, 0.0, rv0_hohmann, "hohmann")
rvs_hohmann_1 = hcat([AstrodynamicsCore.eph(orbit0_hohmann, t) for t in LinRange(0.0, orbit0_hohmann.period/2, 100)]...)

# hohmann option 2
RP0 = R0 * (1 + ECC0)
RA0 = RF * (1 + ECCf)
SMA_hohmann = (RP0 + RA0) / 2.0
ECC_hohmann = (RA0 - RP0) / (RA0 + RP0)
rv0_hohmann = AstrodynamicsCore.kep2rv([SMA_hohmann, ECC_hohmann, 0.0, 0.0, π, 0.0], MU)
orbit0_hohmann = AstrodynamicsCore.Planet(MU, 0.0, rv0_hohmann, "hohmann")
rvs_hohmann_2 = hcat([AstrodynamicsCore.eph(orbit0_hohmann, t) for t in LinRange(0.0, orbit0_hohmann.period/2, 100)]...)

# plot
fontsize = 18
fig = Figure(size=(400,400))
ax = Axis(fig[1,1]; xlabel="x", ylabel="y", aspect=DataAspect(), xlabelsize=fontsize, ylabelsize=fontsize, xticklabelsize=fontsize, yticklabelsize=fontsize)
poly!(Circle(Point2f(0, 0), 0.5), color = :dodgerblue, label=nothing)
hidedecorations!(ax)
lines!(ax, rvs_initial[1,:], rvs_initial[2,:], color=:blue, label="Initial")
lines!(ax, rvs_final[1,:], rvs_final[2,:], color=:green, label="Final", linewidth=2.0)

lines!(ax, rvs_hohmann_1[1,:], rvs_hohmann_1[2,:], color=:red, label="HT 1", linewidth=2.0)
scatter!(ax, rvs_hohmann_1[1,1], rvs_hohmann_1[2,1], color=:red, markersize=10)
scatter!(ax, rvs_hohmann_1[1,end], rvs_hohmann_1[2,end], color=:red, markersize=10)

lines!(ax, rvs_hohmann_2[1,:], rvs_hohmann_2[2,:], color=:orange, label="HT 2", linewidth=2.0)
scatter!(ax, rvs_hohmann_2[1,1], rvs_hohmann_2[2,1], color=:orange, markersize=10)
scatter!(ax, rvs_hohmann_2[1,end], rvs_hohmann_2[2,end], color=:orange, markersize=10)

Legend(fig[2,1], ax, orientation=:horizontal, labelsize=fontsize)
save(joinpath(@__DIR__, "plots/hohmann_ell2ell_opposite.png"), fig; px_per_unit=5)
display(fig)


# --------------------------------------------------------------------------------------- #
# elliptical to elliptical (same side apsis)
MU = 1.0
R0 = 1.9
ECC0 = 0.45
rv0 = AstrodynamicsCore.kep2rv([R0, ECC0, 0.0, 0.0, 0.0, 0.0], MU)
orbit0 = AstrodynamicsCore.Planet(MU, 0.0, rv0, "initial")
rvs_initial = hcat([AstrodynamicsCore.eph(orbit0, t) for t in LinRange(0.0, orbit0.period, 100)]...)

RF = 7.0
ECCf = 0.33
rvf = AstrodynamicsCore.kep2rv([RF, ECCf, 0.0, 0.0, 0.0, 0.0], MU)
orbitf = AstrodynamicsCore.Planet(MU, 0.0, rvf, "final")
rvs_final = hcat([AstrodynamicsCore.eph(orbitf, t) for t in LinRange(0.0, orbitf.period, 100)]...)

# hohmann option 1
RP0 = R0 * (1 - ECC0)
RA0 = RF * (1 + ECCf)
SMA_hohmann = (RP0 + RA0) / 2.0
ECC_hohmann = (RA0 - RP0) / (RA0 + RP0)
rv0_hohmann = AstrodynamicsCore.kep2rv([SMA_hohmann, ECC_hohmann, 0.0, 0.0, 0.0, 0.0], MU)
orbit0_hohmann = AstrodynamicsCore.Planet(MU, 0.0, rv0_hohmann, "hohmann")
rvs_hohmann_1 = hcat([AstrodynamicsCore.eph(orbit0_hohmann, t) for t in LinRange(0.0, orbit0_hohmann.period/2, 100)]...)

# hohmann option 2
RP0 = R0 * (1 + ECC0)
RA0 = RF * (1 - ECCf)
SMA_hohmann = (RP0 + RA0) / 2.0
ECC_hohmann = (RA0 - RP0) / (RA0 + RP0)
rv0_hohmann = AstrodynamicsCore.kep2rv([SMA_hohmann, ECC_hohmann, 0.0, 0.0, π, 0.0], MU)
orbit0_hohmann = AstrodynamicsCore.Planet(MU, 0.0, rv0_hohmann, "hohmann")
rvs_hohmann_2 = hcat([AstrodynamicsCore.eph(orbit0_hohmann, t) for t in LinRange(0.0, orbit0_hohmann.period/2, 100)]...)

# plot
fontsize = 18
fig = Figure(size=(400,400))
ax = Axis(fig[1,1]; xlabel="x", ylabel="y", aspect=DataAspect(), xlabelsize=fontsize, ylabelsize=fontsize, xticklabelsize=fontsize, yticklabelsize=fontsize)
poly!(Circle(Point2f(0, 0), 0.5), color = :dodgerblue, label=nothing)
hidedecorations!(ax)
lines!(ax, rvs_initial[1,:], rvs_initial[2,:], color=:blue, label="Initial")
lines!(ax, rvs_final[1,:], rvs_final[2,:], color=:green, label="Final", linewidth=2.0)

lines!(ax, rvs_hohmann_1[1,:], rvs_hohmann_1[2,:], color=:red, label="HT 1", linewidth=2.0)
scatter!(ax, rvs_hohmann_1[1,1], rvs_hohmann_1[2,1], color=:red, markersize=10)
scatter!(ax, rvs_hohmann_1[1,end], rvs_hohmann_1[2,end], color=:red, markersize=10)

lines!(ax, rvs_hohmann_2[1,:], rvs_hohmann_2[2,:], color=:orange, label="HT 2", linewidth=2.0)
scatter!(ax, rvs_hohmann_2[1,1], rvs_hohmann_2[2,1], color=:orange, markersize=10)
scatter!(ax, rvs_hohmann_2[1,end], rvs_hohmann_2[2,end], color=:orange, markersize=10)

Legend(fig[2,1], ax, orientation=:horizontal, labelsize=fontsize)
save(joinpath(@__DIR__, "plots/hohmann_ell2ell_same.png"), fig; px_per_unit=5)
display(fig)