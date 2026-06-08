# Mean Elements (J2)

```@meta
CurrentModule = AstrodynamicsCore
```

Under J2 perturbation, osculating Keplerian elements oscillate about slowly varying mean values.
AstrodynamicsCore provides first-order mean–osculating mappings following Schaub, *Appendix G*.

## Convention

Mean-element vectors use **mean anomaly** `M` in the sixth component:

```julia
mean_kep = [a, e, i, Ω, ω, M]
```

This differs from the osculating convention `[a, e, i, Ω, ω, θ]` used in
[Orbital Elements](@ref "Orbital Elements") and [Getting Started](@ref "Getting Started"), where the
sixth entry is true anomaly `θ`.

## Physical constants

The mappings require equatorial radius `Re` and J2 coefficient. The test suite uses Earth values:

```julia
Re = 6378.1363    # km
J2 = 1082.63e-6
```

## Osculating → mean

Convert osculating elements to their J2-averaged counterparts with [`kep_osc2mean`](@ref):

```julia
using AstrodynamicsCore

osc_coe = [
    42164.0
    0.0003
    deg2rad(0.14)
    deg2rad(30)
    0.0
    0.0
]

mean_coe = kep_osc2mean(osc_coe, Re, J2)
# mean_coe ≈ [4.2164e4, 2.628e-4, 0.00244, 0.5236, 3.64e-11, -3.69e-11]
```

The GEO-like example above has small inclination and near-zero `ω` and `M`, so the mean elements
are close to the osculating values but not identical.

## Mean → osculating round-trip

[`kep_mean2osc`](@ref) maps mean elements back to osculating elements. A round-trip should recover
the original osculating set:

```julia
osc_coe = [
    42165.727783078
    0.0003
    0.1389338295990
    deg2rad(30)
    0.0
    0.0
]

mean_coe = kep_osc2mean(osc_coe, Re, J2)
osc_back = kep_mean2osc(mean_coe, Re, J2)

osc_back[1] ≈ osc_coe[1]        # semimajor axis
osc_back[2:end] ≈ osc_coe[2:end]
```

## See also

- [Mean Elements (J2) API](@ref "Mean Elements (J2) API") — full docstrings for `kep_osc2mean` and `kep_mean2osc`
- [Orbital Elements](@ref "Orbital Elements") — osculating Keplerian and MEE conversions
