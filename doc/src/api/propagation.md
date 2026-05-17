# Propagation

```@meta
CurrentModule = AstrodynamicsCore
```

Two-body propagation uses the universal-variable Lagrange coefficients from Déprit (1996).

```@docs
propagate_lagrangian
time_until_θ
hypertrig_s
hypertrig_c
```

The vectorized call `propagate_lagrangian(μ, rv0, t0, ts)` accepts a `Vector` or `LinRange` of times and returns a `6 × n` matrix of states.
