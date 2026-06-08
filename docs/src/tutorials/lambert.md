# Solving Lambert Problems

```@meta
CurrentModule = AstrodynamicsCore
```

The Lambert boundary-value problem finds the departure and arrival velocities that connect two position
vectors in a given time of flight (TOF). AstrodynamicsCore implements Izzo's fast Lambert solver (kep3
port) and analytical Jacobians from Arora et al. (2015).

## Problem setup

States are Cartesian position vectors in consistent length and time units. The argument `m` is the
**maximum number of full revolutions** to search for (kep3 `multi_revs` semantics). Set `m = 0` for the
standard single-arc transfer.

```julia
using AstrodynamicsCore

r1 = [0.79, 0.0, 0.0]
r2 = [-0.6, -0.17, 0.015]
μ = 1.0
m = 0
tof = 2.5
```

Call [`lambert`](@ref) with initial and final positions, TOF, `m`, and gravitational parameter `μ`.
When `m == 0`, a successful solve returns a [`LambertResults`](@ref) with `exitflag == 1` and velocity
vectors `v1`, `v2`.

```julia
res = lambert(r1, r2, tof, m, μ)
res.exitflag  # 1 on success
res.v1        # departure velocity
res.v2        # arrival velocity
```

Verify the transfer by propagating the departure state with [`propagate_lagrangian`](@ref):

```julia
x0 = [r1; res.v1]
rv_final = propagate_lagrangian(μ, x0, 0.0, tof)
rv_final ≈ [r2; res.v2]  # position and velocity at arrival
```

## Clockwise (retrograde) transfers

Pass `cw = true` to select the retrograde transfer arc. The short way is inferred automatically from
the geometry (the angular momentum vector must have a non-zero `z` component).

```julia
res_ccw = lambert(r1, r2, tof, 0, μ, false)  # default: counter-clockwise
res_cw  = lambert(r1, r2, tof, 0, μ, true)   # retrograde
```

Both solutions satisfy the same boundary conditions but follow different arcs.

## Multi-revolution solutions

When `m > 0`, [`lambert`](@ref) returns a [`LambertMultiResults`](@ref) containing up to `2 * Nmax + 1`
solution branches, where `Nmax ≤ m` is the largest revolution count feasible for the given TOF:

| Index pattern | `revs[k]` | `branch[k]` | Meaning |
|---------------|-----------|-------------|---------|
| 1 | 0 | `:zero` | Direct (0-rev) transfer |
| 2, 3 | 1 | `:left`, `:right` | 1-rev left and right branches |
| 4, 5 | 2 | `:left`, `:right` | 2-rev left and right branches |
| … | … | … | … |

```julia
r1 = [1.0, 0.0, 0.0]
r2 = [0.0, 1.0, 0.0]
tof = 8.0

res = lambert(r1, r2, tof, 1, μ)   # search up to 1 full revolution
res isa LambertMultiResults        # true when m > 0
length(res.v1)                     # 3 when Nmax == 1

res.revs    # [0, 1, 1]
res.branch  # [:zero, :left, :right]
res.x       # Izzo x variable per branch

# Pick a branch and propagate
k = 2   # e.g. 1-rev left branch
x0 = [r1; res.v1[k]]
rv_final = propagate_lagrangian(μ, x0, 0.0, tof)
```

Not every branch may exist for a given TOF: if the non-dimensional time is too short, `Nmax` is reduced
automatically and fewer than `2 * m + 1` solutions are returned.

## Jacobians

[`lambert_jac`](@ref) returns the Lambert solution together with partial derivatives of the stacked
velocity vector `v = [v1; v2]`. Jacobian support is limited to the single 0-revolution solution
(`m == 0`); the `cw` keyword is supported.

| Matrix | Size | Meaning |
|--------|------|---------|
| `dv_dt` | 6×2 | ∂v/∂t1, ∂v/∂t2 with TOF = t2 − t1 |
| `dv_dr` | 6×6 | ∂v/∂r1 (cols 1–3), ∂v/∂r2 (cols 4–6) |

```julia
res, dv_dt, dv_dr = lambert_jac(r1, r2, tof, 0, μ)

size(dv_dt)  # (6, 2)
size(dv_dr)  # (6, 6)

# Only TOF couples to time, so the two time columns are opposite:
dv_dt[:, 2] ≈ -dv_dt[:, 1]
```

## See also

- [Kepler Propagation](@ref "Kepler Propagation") — propagate Lambert departure states with `propagate_lagrangian`
- [Getting Started](@ref "Getting Started") — element conventions and quick examples
