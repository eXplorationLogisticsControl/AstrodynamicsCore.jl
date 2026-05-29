# Mean Elements (J2)

```@meta
CurrentModule = AstrodynamicsCore
```

First-order mean–osculating mappings follow Schaub, *Appendix G*.

```@docs
kep_osc2mean
kep_mean2osc
```

Element vectors use mean anomaly in the sixth component: `[a, e, i, Ω, ω, M]`.
