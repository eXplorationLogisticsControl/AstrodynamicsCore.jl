# AstrodynamicsCore.jl

Core astrodynamics routines for orbit mechanics in Julia: element conversions, two-body propagation, J2 mean–osculating mappings, and simple transfer helpers.

<p align="center">
  <a href="https://github.com/eXplorationLogisticsControl/AstrodynamicsCore.jl/actions/workflows/test.yml">
    <img src="https://github.com/eXplorationLogisticsControl/AstrodynamicsCore.jl/actions/workflows/test.yml/badge.svg" alt="CI tests"/>
  </a>
  <a href="https://github.com/eXplorationLogisticsControl/AstrodynamicsCore.jl/actions/workflows/docs.yml">
    <img src="https://github.com/eXplorationLogisticsControl/AstrodynamicsCore.jl/actions/workflows/docs.yml/badge.svg" alt="Documentation build"/>
  </a>
</p>

## Overview

`AstrodynamicsCore` is a lightweight library for foundational orbital mechanics. It targets mission analysis workflows that need reliable Keplerian and modified equinoctial element conversions, universal-variable propagation (Déprit formulation), and first-order J2 mean–osculating element mappings (Schaub, Appendix G).

The package also includes utilities for sampling trajectories on a `Planet` model, estimating tangential circularization ΔV, and building simple visualization meshes.

## Features

- **Keplerian & MEE conversions** — `kep2rv`, `rv2kep`, `kep2mee`, `mee2rv`, and related helpers
- **Anomaly maps** — mean, eccentric, and true anomaly (`ma2ea`, `ma2ta`, `ta2ea`, `ta2ma`)
- **Propagation** — `propagate_lagrangian` with optional state transition matrix; batch propagation over time grids
- **J2 mean elements** — `kep_osc2mean` and `kep_mean2osc`
- **Transfers & ephemeris** — `Planet`, `eph`, `time_until_θ`, `tangential_circularize`

## Installation

Clone the repository and activate the project environment:

```julia
using Pkg
Pkg.activate("/path/to/AstrodynamicsCore.jl")
Pkg.instantiate()
```

## Quick example

```julia
using AstrodynamicsCore

μ = 398600.435507
kep = [8000.0, 0.12, deg2rad(85), deg2rad(240), deg2rad(200), deg2rad(130)]

rv = kep2rv(kep, μ)
rv1 = propagate_lagrangian(μ, rv, 0.0, 1000.0)
kep_back = rv2kep(rv1, μ)
```

Runnable scripts are in the `examples/` directory (`ex_elements.jl`, `ex_planet.jl`, `plot_hohmann.jl`, and others).

## Documentation

The API manual is built with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). The **Documentation build** badge above reflects the [`docs`](.github/workflows/docs.yml) GitHub Actions workflow.

Build locally from the repository root:

```julia
using Pkg
Pkg.activate("doc")
Pkg.instantiate()
include("doc/make.jl")
```

Open `doc/build/index.html` in a browser. On CI, `deploydocs` in `doc/make.jl` can publish to GitHub Pages when `GITHUB_ACTIONS` is set.

## Testing

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

## License

This package is released under the [MIT License](LICENSE).
