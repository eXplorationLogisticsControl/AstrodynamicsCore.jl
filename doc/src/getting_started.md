# Getting Started

## Element conventions

Keplerian elements are stored as six-vectors

```julia
kep = [a, e, i, Ω, ω, θ]
```

| Index | Symbol | Description |
|------:|--------|-------------|
| 1 | ``a`` | Semimajor axis (km); negative for hyperbolic orbits |
| 2 | ``e`` | Eccentricity |
| 3 | ``i`` | Inclination (rad) |
| 4 | ``Ω`` | Right ascension of ascending node (rad) |
| 5 | ``ω`` | Argument of perigee (rad) |
| 6 | ``θ`` | True anomaly (rad) |

States are Cartesian vectors `[r_x, r_y, r_z, v_x, v_y, v_z]` in consistent length and time units.

## Keplerian ↔ Cartesian

```julia
using AstrodynamicsCore

μ = 398600.435507
kep = [7000.0, 0.05, deg2rad(28), deg2rad(40), deg2rad(10), deg2rad(75)]

rv = kep2rv(kep, μ)
kep ≈ rv2kep(rv, μ)  # round-trip (compare angles with care near 0/π)
```

## Propagation

Propagate an initial state with the universal-variable Lagrange formulation:

```julia
rv0 = kep2rv(kep, μ)
rv1 = propagate_lagrangian(μ, rv0, 0.0, 600.0)

# optional state transition matrix
rv1, Φ = propagate_lagrangian(μ, rv0, 0.0, 600.0, 1e-12, 20, false, true)

# multiple epochs
rvs = propagate_lagrangian(μ, rv0, 0.0, LinRange(0.0, 600.0, 50))
```

## J2 mean–osculating mapping

```julia
Re = 6378.1363
J2 = 1082.63e-6
osc = [42164.0, 3e-4, deg2rad(0.14), deg2rad(30), 0.0, 0.0]

mean = kep_osc2mean(osc, Re, J2)
osc_back = kep_mean2osc(mean, Re, J2)
```

## Examples

Runnable scripts live under `examples/` in the repository (`ex_elements.jl`, `ex_planet.jl`, `plot_hohmann.jl`, etc.).
