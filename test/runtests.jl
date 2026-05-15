"""Run tests"""

using Test

include(joinpath(@__DIR__, "test_elements.jl"))
include(joinpath(@__DIR__, "test_kepler.jl"))
include(joinpath(@__DIR__, "test_kep_mean_osc.jl"))