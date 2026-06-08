# Kepler Propagation

```@meta
CurrentModule = AstrodynamicsCore
```

Two-body propagation uses the universal-variable Lagrange formulation. The main
entry point is [`propagate_lagrangian`](@ref), which accepts a Cartesian state `[r; v]` and
propagates it from time `t0` to time `t1`.

## Single-epoch propagation (no STM)

The simplest call returns only the final state:

```julia
using AstrodynamicsCore

μ = 398600.435507
rv0 = [
    -4022.024182950954, -6221.706935396668, -4255.655925778433,
    -1.519176635355084, -3.564858730446071, 5.33536331860247,
]

Δt = 1000.0
rv1 = propagate_lagrangian(μ, rv0, 0.0, Δt)
```

A round-trip check confirms the solver is reversible: propagate forward, then backward to recover
the initial state.

```julia
rv0_back = propagate_lagrangian(μ, rv1, Δt, 0.0)
norm(rv0_back[1:3] - rv0[1:3])  # < 1e-11
norm(rv0_back[4:6] - rv0[4:6])  # < 1e-11
```

## With state transition matrix

Set `stm=true` to also obtain the 6×6 state transition matrix Φ, which maps perturbations in the
initial state to the final state: `δx₁ = Φ δx₀`.

```julia
rv1, Φ = propagate_lagrangian(μ, rv0, 0.0, Δt, 1e-12, 20, false, true)

size(Φ)  # (6, 6)

# Final state matches the no-STM call
rv1_no_stm = propagate_lagrangian(μ, rv0, 0.0, Δt)
rv1 ≈ rv1_no_stm
```

The full signature is:

```julia
propagate_lagrangian(μ, rv0, t0, t1, tol, maxiter, verbose, stm)
```

| Argument | Default | Description |
|----------|---------|-------------|
| `tol` | `1e-14` | Laguerre-correction tolerance |
| `maxiter` | `20` | Maximum iterations |
| `verbose` | `false` | Print convergence messages |
| `stm` | `false` | Return state transition matrix |

## Batch propagation

Pass a `Vector` or `LinRange` of target times to obtain a `6 × n` matrix of states at each epoch.
The first column corresponds to `t0`.

```julia
ts = LinRange(0.0, 500.0, 4)
rvs = propagate_lagrangian(μ, rv0, 0.0, ts)

size(rvs)        # (6, 4)
rvs[:, 1] ≈ rv0  # initial state at t = 0
```

## Hyperbolic orbits

`propagate_lagrangian` handles hyperbolic trajectories when the initial state comes from a negative
semimajor axis:

```julia
kep_hyp = [-25000.0, 1.4, deg2rad(25), deg2rad(40), deg2rad(15), deg2rad(5)]
rv_hyp = kep2rv(kep_hyp, μ)
rv1 = propagate_lagrangian(μ, rv_hyp, 0.0, 200.0)
norm(rv1[1:3])  # > 0
```

## See also

- [Propagation API](@ref "Propagation API") — full docstrings for `propagate_lagrangian` and related utilities
- [Lambert Problem](@ref "Lambert Problem") — solve boundary-value problems and verify with propagation
- [Orbital Elements](@ref "Orbital Elements") — build initial states from Keplerian elements
