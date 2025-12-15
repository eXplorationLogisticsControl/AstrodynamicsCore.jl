"""Elements functions"""


"""
Convert 3-element vector from perifocal to geocentric frame
"""
function perifocal2geocentric(vec_pf::Vector, ω::Real, i::Real, Ω::Real)
    # rotate by ω
    v1 = _rotmat_ax3(-ω) * vec_pf
    # rotate by inclination
    v2 = _rotmat_ax1(-i) * v1
    # rotate by Ω
    v_gec = _rotmat_ax3(-Ω) * v2
    return v_gec
end


function kep2rv(kep::Array{<:Real,1}, μ::Real)
	# unpack
	a, e, i, Ω, ω, θ = kep
    # angular momentum
    if e <= 1.0
        h = sqrt(a * μ * (1 - e^2))
    else
        h = sqrt(abs(a) * μ * (e^2 - 1))   # hyperbolic case
    end
    # perifocal vector
    x = (h^2 / μ) * (1 / (1 + e * cos(θ))) * cos(θ)
    y = (h^2 / μ) * (1 / (1 + e * cos(θ))) * sin(θ)
    vx = (μ / h) * (-sin(θ))
    vy = (μ / h) * (e + cos(θ))
    rpf = [x, y, 0.0]
    vpf = [vx, vy, 0.0]

    r_rot3 = perifocal2geocentric(rpf, ω, i, Ω)
    v_rot3 = perifocal2geocentric(vpf, ω, i, Ω)

    # save inertial state
    state_inr = vcat(r_rot3, v_rot3)[:]
    return state_inr
end


function kep2meesma(kep::Array{<:Real,1})
    # unpack
    a, e, i, Ω, ω, θ = kep
    # compute MEEs
    f = e*cos(Ω+ω)
    g = e*sin(Ω+ω)
    h = tan(i/2)*cos(Ω)
    k = tan(i/2)*sin(Ω)
    l = Ω + ω + θ
    return [a,f,g,h,k,l]
end


function meesma2kep(meesma::Array{<:Real,1})
    # unpack
    a, f, g, h, k, l = meesma
    # compute Keplerian elements
    e = sqrt(f^2 + g^2)
    i = atan(2*sqrt(h^2+k^2), 1-h^2-k^2)
    Ω = atan(k,h)
    ω = atan(g*h-f*k, f*h+g*k)
    θ = l - Ω - ω
    return [a,e,i,Ω,ω,θ]
end


function meesma2rv(meesma::Array{<:Real,1}, μ::Real)
    return kep2rv(meesma2kep(meesma), μ)
end
