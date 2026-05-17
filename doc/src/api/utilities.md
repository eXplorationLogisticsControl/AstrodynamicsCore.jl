# Utilities

```@meta
CurrentModule = AstrodynamicsCore
```

```@docs
get_sphere_mesh
moving_average
```

`get_sphere_mesh` returns a `GeometryBasics.Mesh` suitable for Makie wireframe plots. `moving_average` computes a trailing mean with optional padding via `_moving_average_padded` (internal).
