# AstrodynamicsCore.jl

```@meta
CurrentModule = AstrodynamicsCore
```

```@docs
AstrodynamicsCore
```

AstrodynamicsCore provides foundational orbit mechanics routines in Julia:

- Keplerian and modified equinoctial element conversions
- Mean / true / eccentric anomaly maps
- Universal-variable propagation (Déprit formulation) with optional STM
- First-order J2 mean–osculating element mappings
- Lightweight `Planet` ephemeris helpers and tangential transfer ΔV estimates

## Installation

The package is developed in-tree. From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using AstrodynamicsCore
```

## Documentation

Build the manual locally from the `doc` directory:

```julia
using Pkg
Pkg.activate("doc")
Pkg.instantiate()
include("make.jl")
```

The HTML site is written to `doc/build/`.

## Quick example

```julia
using AstrodynamicsCore

μ = 398600.435507
kep = [8000.0, 0.12, deg2rad(85), deg2rad(240), deg2rad(200), deg2rad(130)]

rv = kep2rv(kep, μ)
kep_back = rv2kep(rv, μ)
```

See [Getting Started](@ref) for additional workflows.
