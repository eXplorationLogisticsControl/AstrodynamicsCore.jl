"""Analytical orbit transfer expressions"""


"""
    tangential_circularize(kep, μ)

ΔV required for tangential circularization at true anomaly `θ` in `kep`.

Returns ``v_{circ} - v_{vis}`` at the current radius (apogee/perigee burns).
"""
function tangential_circularize(kep::Array{<:Real,1}, μ::Real)
    a,e,_,_,_,θ = kep
    #@assert mod(θ, π) < 1e-8 "θ must be 0 or π"
    ra = a * (1 + e)
    v_vis  = sqrt(μ*(2/ra - 1/a))
    v_circ = sqrt(μ/ra)
    return v_circ - v_vis
end

