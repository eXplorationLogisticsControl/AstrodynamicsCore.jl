# `AstrodynamicsCore.jl`: core astrodynamics routines

<p align="center">
  <img src="https://github.com/eXplorationLogisticsControl/AstrodynamicsCore.jl/actions/workflows/test.yml/badge.svg" alt="test workflow"/>
</p>

## Documentation

Build the [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) manual from the `doc` folder:

```julia
using Pkg
Pkg.activate("doc")
Pkg.instantiate()
include("doc/make.jl")
```

Open `doc/build/index.html` in a browser. On CI, the site can be deployed to GitHub Pages via `deploydocs` in `doc/make.jl`.