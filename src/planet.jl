"""Planet and ephemeris helpers."""

"""
    Planet(μ, t0, rv0, name="Planet")

Simple Keplerian planet model with constant `μ`, epoch `t0`, and initial state `rv0`.

The orbital period is inferred from `rv0` and stored in `planet.period`.
"""
mutable struct Planet
    μ::Real
    t0::Real
    rv0::Array{<:Real,1}
    period::Real
    name::String

    function Planet(μ::Real, t0::Real, rv0::Array{<:Real,1}, name::String="Planet")
        # compute period
        r, v = norm(rv0[1:3]), norm(rv0[4:6])
        a = (r * μ) / (2*μ - r*v^2)
        period = 2π * sqrt(a^3 / μ)
        new(μ, t0, rv0, period, name)
    end
end


function Base.show(io::IO, planet::Planet)
    @printf(io, "Planet %s\n", planet.name)
    @printf(io, "  μ      : %1.4e\n", planet.μ)
    @printf(io, "  t0     : %1.4e\n", planet.t0)
    @printf(io, "  period : %1.4e\n", planet.period)
end


"""
    eph(planet, t)

Propagate `planet` to time `t` and return the inertial state `[r; v]`.
"""
function eph(planet::Planet, t::Real)
    Δt = t - planet.t0
    return propagate_lagrangian(planet.μ, planet.rv0, planet.t0, Δt)
end