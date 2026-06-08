"""
    AstrodynamicsCore

Core astrodynamics routines for orbital element conversions, Keplerian propagation,
J2 mean/osculating mappings, and simple transfer helpers.
"""
module AstrodynamicsCore

using FileIO
using GeometryBasics
using LinearAlgebra
using Printf

include("misc.jl")
include("transformations.jl")
include("elements.jl")
include("kepler.jl")
include("lambert.jl")
include("elements_mean.jl")

include("planet.jl")
include("transfer.jl")

export perifocal2geocentric, kep2rv, rv2kep, kep2mee, mee2kep, mee2rv, rv2mee
export ma2ea, ma2ta, ta2ea, ta2ma
export kep_osc2mean, kep_mean2osc
export propagate_lagrangian
export LambertResults, LambertMultiResults, lambert, lambert_jac
export Planet, eph

end # module AstrodynamicsCore
