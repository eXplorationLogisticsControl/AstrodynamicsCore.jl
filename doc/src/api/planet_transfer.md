# Planets & Transfers

```@meta
CurrentModule = AstrodynamicsCore
```

```@docs
Planet
eph
tangential_circularize
```

`Planet` stores `μ`, epoch `t0`, initial state `rv0`, inferred `period`, and a `name` label. Use [`eph`](@ref) to sample trajectories along the Keplerian arc.
