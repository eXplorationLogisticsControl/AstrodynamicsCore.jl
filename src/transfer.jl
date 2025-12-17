"""Analytical orbit transfer expressions"""


function tangential_circularize(kep::Array{<:Real,1}, μ::Real)
    a,e,_,_,_,θ = kep
    #@assert mod(θ, π) < 1e-8 "θ must be 0 or π"
    ra = a * (1 + e)
    v_vis  = sqrt(μ*(2/ra - 1/a))
    v_circ = sqrt(μ/ra)
    return v_circ - v_vis
end

