# Orbital Elements

```@meta
CurrentModule = AstrodynamicsCore
```

AstrodynamicsCore converts between Keplerian elements, modified equinoctial elements (MEE), and
Cartesian states `[r_x, r_y, r_z, v_x, v_y, v_z]`.

## Conventions

**Keplerian elements** are stored as six-elements vector with **true anomaly** in the sixth component:

```julia
kep = [a, e, i, Ω, ω, θ]
```

**Modified equinoctial elements** are stored as six-elements vector:

```julia
mee = [p,f,g,h,k,L]
```

where `L` is the osculating true longitude.

See [Getting Started](@ref "Getting Started") for the full symbol table. Length and time units must be
consistent with the gravitational parameter `μ` (km and seconds in the examples below).

## Keplerian ↔ Cartesian

The NAIF reference case below is taken from the package test suite (`test/test_elements.jl`).

```julia
using AstrodynamicsCore

μ = 398600.435507
kep = [8000, 0.12, deg2rad(85), deg2rad(240), deg2rad(200), deg2rad(130)]

rv = kep2rv(kep, μ)
# rv ≈ [-4022.02, -6221.71, -4255.66, -1.51918, -3.56486, 5.33536]

kep_back = rv2kep(rv, μ)
```

Angles returned by `rv2kep` may differ from the input by multiples of 2π. Compare with care:

```julia
Δθ = acos(cos(kep_back[6] - kep[6]))  # ≈ 0
```

A second test case with different geometry:

```julia
kep2 = [26590, 0.06, deg2rad(15), deg2rad(109), deg2rad(27), deg2rad(58)]
rv2 = kep2rv(kep2, μ)
kep2 ≈ rv2kep(rv2, μ)  # round-trip (check angles individually)
```

## Modified equinoctial elements (MEE) ↔ Cartesian

To convert between MEE and Cartesian state vector, use

```julia
rv  = mee2rv(mee, μ)
```

and

```julia
mee = rv2mee(rv, μ)
```

For retrograde orbits, `rv2mee` accepts an optional `retrograde` flag (defaults to `false`):

```julia
mee_retro = rv2mee(rv, μ, true)
```

## Modified equinoctial elements (MEE) ↔ Keplerian elements

To convert between MEE and Keplerian elements, use

```julia
mee = kep2mee(kep)
```

and 

```julia
kep = mee2kep(mee)
```


## Perifocal frame

[`kep2rv_perifocal`](@ref) returns the state expressed in the perifocal frame. For an equatorial
orbit (`i = Ω = ω = 0`) it agrees with the inertial-frame result from `kep2rv`:

```julia
θ = deg2rad(130)
kep_eq = [8000.0, 0.12, 0.0, 0.0, 0.0, θ]

rv_pf = kep2rv_perifocal(kep_eq, μ)
rv    = kep2rv(kep_eq, μ)
rv ≈ rv_pf
```

## Edge cases

**Equatorial circular orbit** — inclination and node are zero:

```julia
r = 7000.0
v = sqrt(μ / r)
rv_circ = [r, 0.0, 0.0, 0.0, v, 0.0]

kep_circ = rv2kep(rv_circ, μ)
kep_circ[3] ≈ 0.0  # i ≈ 0
kep2rv(kep_circ, μ) ≈ rv_circ
```

**Hyperbolic orbit** — negative semimajor axis:

```julia
kep_hyp = [-25000.0, 1.4, deg2rad(25), deg2rad(40), deg2rad(15), deg2rad(5)]
rv_hyp = kep2rv(kep_hyp, μ)
rv2kep(rv_hyp, μ)[2] ≈ kep_hyp[2]  # eccentricity preserved
```

## See also

- [Orbital Elements API](@ref "Orbital Elements API") — full docstrings for all conversion functions
- [Kepler Propagation](@ref "Kepler Propagation") — propagate Cartesian states in time
- [Mean Elements (J2)](@ref "Mean Elements (J2)") — mean–osculating mappings (uses mean anomaly, not true anomaly)
