# Solving Lambert Problems

```@meta
CurrentModule = AstrodynamicsCore
```

The Lambert boundary-value problem finds the departure and arrival velocities that connect two position
vectors in a given time of flight (TOF). AstrodynamicsCore implements Izzo's fast Lambert solver and
analytical Jacobians from Arora et al. (2015).

## Problem setup

States are Cartesian position vectors in consistent length and time units, together with the maximum number of revolutions `m`:

```julia
using AstrodynamicsCore

r1 = [0.79, 0.0, 0.0]
r2 = [-0.6, -0.17, 0.015]
μ = 1.0
m = 0
tof = 2.5
```

We then call [`lambert`](@ref) with initial and final positions, TOF, revolution count, and gravitational parameter `μ`. 
A successful solve returns `exitflag == 1` with velocity vectors `v1`, `v2`.

```julia
res = lambert(r1, r2, tof, m, μ)
res.exitflag  # 1 on success
res.v1        # departure velocity
res.v2        # arrival velocity
```


## Jacobians

[`lambert_jac`](@ref) returns the Lambert solution together with partial derivatives of the stacked
velocity vector `v = [v1; v2]`:

| Matrix | Size | Meaning |
|--------|------|---------|
| `dv_dt` | 6×2 | ∂v/∂t1, ∂v/∂t2 with TOF = t2 − t1 |
| `dv_dr` | 6×6 | ∂v/∂r1 (cols 1–3), ∂v/∂r2 (cols 4–6) |

```julia
res, dv_dt, dv_dr = lambert_jac(r1, r2, tof, m, μ)

size(dv_dt)  # (6, 2)
size(dv_dr)  # (6, 6)

# Only TOF couples to time, so the two time columns are opposite:
dv_dt[:, 2] ≈ -dv_dt[:, 1]
```

## See also

- [Kepler Propagation](@ref "Kepler Propagation") — propagate Lambert departure states with `propagate_lagrangian`
- [Getting Started](@ref "Getting Started") — element conventions and quick examples
